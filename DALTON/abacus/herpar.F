
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
C File: herpar.F
C
C Some history:
C Original : PVM
C June 1995: ported to MPI by Kenneth Ruud
C
C  /* Deck her_pardrv */
      SUBROUTINE HER_PARDRV(FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                      WORK,LWORK,ITYPE,MAXDER,IATOM,NODV,
     &                      NOPV,NOCONT,TKTIME,RETUR,IPRINT,
     &                      ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,
     &                      DINTSKP,RELCALC,GENCTR)

C  Hermit parallel driver.  Only called from PARDRV (checked Nov. 2016).

C Note Nov. 2016 hjaaj: GENCTR is apparently related to some BUILD_OPENRSP
C work, but no comment explains the meaning. Was introduced in Jan. 2014
C according to "git blame".  Currently it is **UNDEFINED** in the calls of PARDRV !!!

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
Cef begin
#include "incore.h"
Cef end
C
      LOGICAL   NODV,NOPV,NOCONT,RETUR,TKTIME,RELCALC, rma_model_save
      LOGICAL   GENCTR
      DIMENSION FMAT(*), DMAT(N2BASX,NDMAT), ISYMDM(*), IFCTYP(*),
     &          NSTAT(NODTOT), WORK(LWORK), HESSEE(*)
      DIMENSION GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)

C
C Used from common blocks:
C
C AOVEC : MXAOVC (for BLOCKS)
C MXCENT: MXCOOR
C MXORB : MXSHEL (for INFPAR)
C INFPAR: NDEGDI,NTASK
C BLOCKS: MAXSHL
C NUCLEI: NUCDEG
C INFORB: NNBASX,N2BASX
C
#include "infpar.h"
#include "blocks.h"
#include "nuclei.h"
#include "inforb.h"
#include "memint.h"
C
      CALL QENTER('HER_PARDRV')
      IF (IPRINT .GT. 2) CALL TITLER('Output from HER_PARDRV','*',103)
C
C     Determine total number of tasks - MTOTTK.
C
      MTOTTK = MAXSHL*(MAXSHL+1)/2
C
      IF (MTOTTK .GT. MAXTSK) THEN
         WRITE(LUPRI,'(//5X,A,/5X,A,/5X,2(A,I5))')
     &    'ERROR : Number of tasks greater than maximum.',
     &    'Recompilation will be necessary',
     &    'MTOTTK =',MTOTTK, 'MAXTSK =',MAXTSK
         CALL RECOMPI
         CALL QUIT('ERROR: MTOTTK greater than MAXTSK in HER_PARDRV')
      ENDIF
C
C     Determine dimension of array TMPMAT
C
C     Expectation values of differentiated integrals (all atoms)
      IF (ITYPE .EQ. 2) THEN
         L_GRA = MXCOOR
         L_HES = L_GRA*L_GRA
         L_TMPMT = MAX(L_GRA,L_HES,MXPRIM) ! for gradient, Hessian, expgrd
C     Direct calculation of Fock matrices
C        using skeleton matrix approach (AO basis)
      ELSE IF (ITYPE .EQ. 3) THEN
         L_FCK = NDMAT*N2BASX
         NUMSKP = 8
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
         L_TMPMT = NUMSKP
#else
         L_TMPMT = L_FCK + NUMSKP
#endif
C     Derivatives with respect to magnetic field(DDFOCK=T)
      ELSE IF (ITYPE .EQ. -5) THEN
         L_SUS = 9
C         L_TMPMT = 3*N2BASX + L_SUS
C modified by Bin Gao, June 3, 2009, from linsca
         NFMAT = 3*NDMAT
#ifdef BUILD_OPENRSP
         IF (USE_OPENRSP) THEN
            IF (MAXDER .EQ. 2 .AND. .NOT.GENCTR) NFMAT = 6*NDMAT
         END IF
#endif
         L_FCK = NFMAT*N2BASX
         L_TMPMT = L_FCK + L_SUS
      ELSE IF (ITYPE .EQ. 8) THEN
C     Derivative Fock matrix for specified atom
         L_TMPMT = 3*NUCDEG(IATOM)*N2BASX
      END IF
C
C     Determine number of tasks pr. batch - NTASK.
Cef begin
C     The case of incore calculations does not support NTASK > 1 yet.
C     The routines DOPAR and PARLOP must be modified to support NTASK > 1.
      IF (AOSAVE) THEN
         NTASK = 1
      ELSE
         NTASK = MAX(1,MIN(MTOTTK,INT(NDEGDI*MTOTTK/(100*NODTOT))))
      END IF
Cef end
C
C     Determine dimension of array INDEX
C
      IBLOCK = NODTOT*NTASK
      NINDEX = (IBLOCK + 1)*INT(MTOTTK/IBLOCK + 1)
C
C     Be sure to allocate enough workspace for array SORTED
C
      IF (L_TMPMT .LT. (MTOTTK + 2)) L_TMPMT = MTOTTK + 2
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK

!     rma model parallelization only for itype == 3, so:
      rma_model_save = rma_model
      if(itype /= 3) rma_model = .false.

      if(.not.rma_model)then
        CALL MEMGET('REAL',KTMPMT,L_TMPMT,WORK,KFREE,LFREE)
      else
        if(itype == 2)then ! expectation values
          CALL MEMGET('REAL',KTMPMT,L_TMPMT,WORK,KFREE,LFREE)
        else ! fock matrices in some way
          CALL MEMGET('REAL',KTMPMT,     0,WORK,KFREE,LFREE)
        end if
      end if

      CALL MEMGET('INTE',KINDEX,3*NINDEX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTIMES,MTOTTK,  WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTKCPU,MTOTTK,  WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTSORT,MTOTTK+2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KWHICH,MTOTTK,  WORK,KFREE,LFREE)
C
      CALL HER_PARDRV_1(FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                WORK(KTMPMT),WORK(KINDEX),WORK(KTIMES),
     &                WORK(KTKCPU),WORK(KWHICH),WORK(KFREE),LFREE,
     &                NINDEX,L_TMPMT,L_GRA,L_HES,L_SUS,L_FCK,
     &                ITYPE,MAXDER,IATOM,NODV,NOPV,
     &                NOCONT,TKTIME,RETUR,IPRINT,
     &                ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,
     &                WORK(KTSORT),GENCTR)
!     restore rma_model in case it had been reset
      rma_model = rma_model_save

      CALL MEMREL('HER_PARDRV',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('HER_PARDRV')
      RETURN
      END
C  /* Deck her_pardr1 */
      SUBROUTINE HER_PARDRV_1(
     &   FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &   TMPMAT,INDEX,TIMES,TSKCPU,IWHICH,WORK,LWORK,
     &   NINDEX,L_TMPMT,L_GRA,L_HES,L_SUS,L_FCK,
     &   ITYPE,MAXDIF,IATOM,NODV,NOPV,
     &   NOCONT,TKTIME,RETUR,IPRINT,
     &   ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,SORTING,
     &   GENCTR)
C
C
C     **************************************************************
C     *     Driver routine for the parallel HERMIT calculation.    *
C     **************************************************************
C
C
C      A short survey of the different messagetags (MTAGX) :
C
C
C      NEW LABELS:
C
C      10 - Tell nodes this is a HERMIT run, send NCODE and IPRINT
C      30 - Send initialization to nodes.
C      40   >Label used in ERI only<
C      50 - Receive (IJ|**) request from node
C      60 - Send new (IJ|**) labels to node
C      70 - Receive results from nodes
C      80 - Receive overall timings
C
C
C      OLD LABELS:
C
C      10 - Send initialization to nodes.
C      20 - Receive request from node.
C      30 - Send batch to node.
C      40 - Receive final results from nodes.
C      50 - Get timing results from nodes.
C
C           NCODE  - PVM-message encoding option.
C           ITYPE  - Calc. type: 2 = expectation values.
C                                3 = fock matrices.
C                               -5 = Integrals derivated with respect
C                                    to magnetic field.
C                                8 = Fock matrices derivated with respect
C                                    to atomic coordinates.
C           NTASK  - Number of tasks
C           IPRPAR - Print level in parallelization routines on master.
C           MAXDIF -
C           MAXREP -
C           ISYMDM - Symmetries of density matrices
C           IFCTYP - Fock matrix types (see twoint).
C           IATOM  -
C           NODV   - 1-el. Density matrix neglected in TWOEXP
C           NOPV   - 2-el.           -   "   -
C           NOCONT -           -   "   -
C           RETUR  - Program will exit after spec. shells
C           TKTIME - Take time in TWOINT
C           TIMING - Take time for each integralbatch IJ
C           DOREPS -
C           DOCOOR -
C           DOSYM  -
C           MULD2H - Start address for common block INFORB.
C           NMLINE - Number of lines in MOLECULE.INP
C           MLINE  - MOLECULE.INP as internal file
C           NDMAT  - Number of density-matrices
C           DMAT   - Density-matrices.
C           HFXFAC - HF exchange factor
C
C
C
C NOTE: If the program get short of memory in calculation of
C       right hand sides (itype = 8) one can reduce the temporary
C       matrix TMPMAT used in HER_RVRES from 3*N2BASX to NNBASX by
C       looping over Fock matrices in HER_SDRES and HER_RVRES
C       and then performing a send each time rather than transferring
C       all matrices in one batch.
C
C
!     module dependencies
      use memory_parallel
      use parallel_communication_models_mpi
      use rma_windows

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
      use mpi
#include "implicit.h"
#else
#include "implicit.h"
#include "mpif.h"
#endif
#else
#include "implicit.h"
#endif

C
#include "maxorb.h"
C
C Used from common blocks
C
C MXORB : MXSHEL (used for INFPAR)
C GNRINF: NEWBAS
C SIRIPC: NEWGEO
C INFPAR: NODTOT, NTASK
C INFORB: N2BASX
C
#include "gnrinf.h"
#include "siripc.h"
#include "inforb.h"
#include "infpar.h"
Cef begin
#include "incore.h"
Cef end
C
Cef Dimension of array INDEX changed from INDEX(NINDEX) to INDEX(3,NINDEX)
      LOGICAL FIRST,NODV,NOPV,NOCONT,RETUR,TKTIME,GENCTR
      DIMENSION FMAT(*), DMAT(N2BASX,NDMAT), ISYMDM(*), IFCTYP(*),
     &          NSTAT(NODTOT), TMPMAT(L_TMPMT), INDEX(3,NINDEX),
     &          TIMES(MTOTTK), TSKCPU(MTOTTK), IWHICH(MTOTTK),
     &          WORK(LWORK), NPOS(0:MAXTSK+1), HESSEE(*), SORTING(*)
CTROND
      DIMENSION GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)
CTROND
#ifdef VAR_MPI
      integer, parameter             :: dummy_buff = 2
      integer(kind=MPI_ADDRESS_KIND) :: nelement
      integer(kind=MPI_OFFSET_KIND)  :: displacement
      character(len=16)              :: shared_dmat_file
      integer                        :: istat(mpi_status_size)
#endif
C
C     NB! NPOS cannot be allocated from WORK due to the SAVEing below.
C
      SAVE FIRST,NPOS
      DATA FIRST/.TRUE./
C
      CALL QENTER('HER_PARDRV_1')
      IF (IPRINT .GT. 4) CALL TITLER('Output from HER_PARDRV_1','*',103)
C
C     Take calculation time of each IJ pair?
C
      TIMING = FIRST .OR. NEWBAS .OR. NEWGEO
C
C--------------------------
C     Make the index array.
C--------------------------
C
      CALL HER_INDEKS(TIMING,INDEX,NINDEX,MTOTTK,NODTOT,NTASK,NBATCH,
     &            NPOS,IPRINT)
C
Cef begin
C 10.mai     INITX = .TRUE.
Cef end
C---------------------------------------
C     Send initialization data to nodes.
C---------------------------------------
C
      CALL HER_SDINIT(DMAT,NDMAT,ISYMDM,IFCTYP,ITYPE,MAXDIF,IATOM,
     &                NODV,NOPV,NOCONT,TKTIME,RETUR,FIRST,IPRINT,
     &                ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,GENCTR)
C
Cef begin
C---------------------------------------------------------------
C   INDEX may need to be reinitialized after CALL HER_SDINIT.
C---------------------------------------------------------------
C
      IF (AOSAVE .AND. (.NOT. INITX)) THEN
         CALL HER_INDEKS(TIMING,INDEX,NINDEX,MTOTTK,NODTOT,NTASK,NBATCH,
     &        NPOS,IPRINT)
C 10.mai         INITX = .TRUE.
      END IF
Cef end
C
!--------------------------------------------------------------------------------------
!     initialize memory windows for rma access within a given shared-memory / NUMA node
!--------------------------------------------------------------------------------------
      if(rma_model)then
         IF (IPRINT .GT. 2) then
            write(lupri,*) 'Using rma model'
            flush(lupri)
         end if

        rma_win_info%dmat_win = -1
        rma_win_info%fmat_win = -1

        nelement            = n2basx*ndmat ! for itype == 3: ndmat == nfmat


!       sknecht: idea for RMA window on master: assign shmem-MASTER to some other member of the master group (via key-value in comm_split)
!                                               and introduce an integer list which keeps track of which FMAT(i)
!                                               was assigned to which process thus introducing a memory split within the numa-node(s) which could be advantageous.
!                                               caution: not making MASTER to shmem-master may hold some pitfalls...


#ifdef VAR_MPI
!       test dmat on shared file
        write(shared_dmat_file,'(a12,i4)')
     &  'dmat_shared.',communication_info_mpi%intra_node_group_id
        call mpi_file_open(
     &                     communication_info_mpi%
     &                     communication_intranode,
     &                     shared_dmat_file,
     &                     mpi_mode_CREATE +
     &                     mpi_mode_RDWR   +
     &                     mpi_mode_DELETE_ON_CLOSE,
     &                     mpi_info_null,
     &                     rma_win_info%dmat_fh,
     &                     ierr
     &                    )
        displacement = 0
        call mpi_file_set_view(
     &                         rma_win_info%dmat_fh,
     &                         displacement,
     &                         mpi_real8,
     &                         mpi_real8,
     &                         'native',
     &                         mpi_info_null,
     &                         ierr
     &                        )
        if(communication_info_mpi%my_intra_node_id == 0)then
          nelement4 = nelement
          call mpi_file_write(
     &                        rma_win_info%dmat_fh,
     &                        dmat,
     &                        nelement4,
     &                        mpi_real8,
     &                        istat,
     &                        ierr
     &                       )
        end if
        call mpi_barrier(communication_info_mpi%
     &                   communication_intranode,ierr)
!       end test dmat on shared file


        call set_rma_window(
     &                      fmat,
     &                      nelement,
     &                      master,
     &                      communication_info_mpi%
     &                      communication_shmemnode,
     &                      rma_win_info%fmat_win,
!    &                      .true.,
     &                      .false.,
     &                      'no_locks',
     &                      'true')

#endif
!       print *, 'fmat_win, dmat_win',rma_win_info%fmat_win,
!    &            rma_win_info%dmat_win,mynum
      end if

C-----------------------
C     Start calculation.
C-----------------------
C
      CALL DOPAR(INDEX,MTOTTK,NBATCH,NSTAT,IPRINT)
C
C--------------------------------------
C     Receive final results from nodes.
C--------------------------------------
C
      CALL HER_RVRES(FMAT,TMPMAT,HESSEE,TIMES,TSKCPU,IWHICH,L_TMPMT,
     &               L_GRA,L_HES,L_SUS,L_FCK,ITYPE,IPRINT,DINTSKP,
     &               GENCTR)
C
C-----------------------------------------
C     Sort the integral-calculation times.
C-----------------------------------------
C
C       NB! Array TMPMAT is only used as
C           workspace for array SORTED.
      IF (TIMING) CALL PARSRT(MTOTTK,NPOS,SORTING,TIMES,IPRINT)

!------------------------
!     free memory windows
!------------------------
      if(rma_model)then
#ifdef VAR_MPI
!       call free_rma_window(rma_win_info%dmat_win,master)
        call free_rma_window(rma_win_info%fmat_win,master)
!       test dmat on shared file
        call mpi_file_close(rma_win_info%dmat_fh,ierr)
!       end test dmat on shared file
#endif
      end if

      IF (FIRST)  FIRST  = .FALSE.
      IF (NEWGEO) NEWGEO = .FALSE.
      IF (NEWBAS) NEWBAS = .FALSE.
C
      CALL QEXIT('HER_PARDRV_1')
      RETURN
      END
C  /* Deck her_indeks */
      SUBROUTINE HER_INDEKS(TIMING,INDEX,NINDEX,MTOTTK,NODTOT,
     &                  NTASK,NBATCH,NPOS,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
Cef begin
#include "maxorb.h"
#include "incore.h"
Cef end
C
Cef Dimension of array INDEX changed from INDEX(NINDEX) to INDEX(3,NINDEX)
C
      DIMENSION INDEX(3,NINDEX), NPOS(0:MTOTTK+1)
      LOGICAL TIMING
C
      IF (IPRINT .GT. 4) CALL TITLER('Output from HER_INDEKS','*',103)
C
C     This subroutine determines in which order integrals will be calculated.
C
      IF (AOSAVE .AND. NINDEX .GT. MXTSK) THEN
         CALL QUIT('INDX_C in incore.h must be enlarged.')
      END IF
C     Avoid reinitialization of INDEX when it has been initialized in
C     a previous iteration
      IF (.NOT.AOSAVE .OR. .NOT.INITX ) THEN
         CALL IZERO(INDEX(1,1),(3*NINDEX))
      ELSE IF (AOSAVE .AND. INITX) THEN
         CALL ICOPY((3*NINDEX),INDX_C(1,1),1,INDEX(1,1),1)
Cef 10.mai:
      END IF
C
      ICOUNT = 0
      IBLOCK = NODTOT*NTASK
C
      DO 100 I = MTOTTK, 1, -1
         MORE1  = (ICOUNT/IBLOCK)*IBLOCK
         MORE2  = MOD(ICOUNT,NODTOT)*NTASK
         MORE3  = ((ICOUNT-MORE1)/NODTOT) + 1
         IPLACE = MORE1+MORE2+MORE3
C
         IF (.NOT.AOSAVE) THEN
            IF (TIMING) THEN
               INDEX(1,IPLACE) = I
            ELSE
               INDEX(1,IPLACE) = NPOS(I)
            END IF
         ELSE IF (AOSAVE .AND. .NOT.INITX) THEN
            INDEX(1,IPLACE) = I
         END IF
C
         ICOUNT=ICOUNT+1
 100  CONTINUE
C
C The code piece below is a trick to make the 2-dim INDEX-array work with old code for
C the routines PARLOP()/DOPAR() which assume a 1-dim INDEX array
      IF (.NOT.AOSAVE) THEN
         DO 111 I = 1, NINDEX
            INDEX(1+MOD(I-1,3),1+INT((I-1)/3)) = INDEX(1,I)
 111     CONTINUE
      END IF
C
C     Determine number of batches (NBATCH).
C
      NOTALL = (MTOTTK/IBLOCK)*IBLOCK
      NREST  = MTOTTK - NOTALL
      IF (NREST .GT. NODTOT) NREST = NODTOT
      NBATCH = NOTALL/NTASK + NREST
C
C     Print current settings.
C
      NEACH = NBATCH/NODTOT
      NREST = NBATCH - NEACH*NODTOT
C
      IF (IPRINT .GT. 4) THEN
         WRITE(LUPRI,'(5(/5X,A,I5),2X,A7,I3,A1/)')
     &        'Number of nodes   :',NODTOT,
     &        'Number of tasks   :',MTOTTK,
     &        'Tasks pr. batch   :',NTASK,
     &        'Number of batches :',NBATCH,
     &        'Batches pr. node  :',NEACH,
     &        '(rest =',NREST,')'
C
C        Print the index array.
C
         WRITE(LUPRI,'(/5X,A/)') 'The order of task-distribution :'
         DO 400 I=1, NBATCH*NTASK
Cef begin
            WRITE(LUPRI,'(5X,A3,I5,A15,I5)')
     &           'I =',I,'index(i) =',INDEX(1,I)
Cef end
 400      CONTINUE
      END IF
C
C      do 314 jj = 1,nbatch*ntask
C         do 315 ii = 1,3
C            write(lupri,*),'i,j = ',ii,jj,'INDEX(i,j)=',INDEX(ii,jj)
C         write(lupri,*),'j = ',jj,'INDEX(1,j)=',INDEX(1,jj)
C 315     continue
C 314  continue
C
Cef 10.mai begin
      INITX = .TRUE.
Cef 10.mai end
      RETURN
      END
C  /* Deck her_sdinit */
      SUBROUTINE HER_SDINIT(DMAT,NDMAT,ISYMDM,IFCTYP,ITYPE,MAXDIF,
     &                      IATOM,NODV,NOPV,NOCONT,TKTIME,RETUR,
     &                      FIRST,IPRINT,
     &                      ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,GENCTR)

!     module dependencies
      use parallel_communication_models_mpi
      use rma_windows

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
      use mpi
#include "implicit.h"
#else
#include "implicit.h"
#include "mpif.h"
#endif
#else
#include "implicit.h"
#endif
C
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mtags.h"
#include "expopt.h"
C
      LOGICAL NODV,NOPV,NOCONT,TKTIME,RETUR,FIRST,GENCTR
      DIMENSION DMAT(N2BASX,NDMAT), ISYMDM(*), IFCTYP(*)
CTROND
CKR      DIMENSION GABRAO(*), DMRAO(*), DMRSO(*), DINTSKP(*)
      DIMENSION GABRAO(*), DMRAO(*), DMRSO(*)
CTROND
#include "infpar.h"
C
C Used from common blocks:
C
C MXCENT : MXCENT (for DORPS, etc.)
C MXORB  : MXSHEL (for INFPAR)
C MAXAQN : *      (for SYMMET)
C SIRIPC : NEWGEO (PVM: BASDIR, LENBAS)
C INFORB : N2BASX
C MOLINP : NMLINE, MLINE()
C CBIREA : UNCONT
C DORPS  : DOREPS(), DOCOOR()
C ABAINF : DOSYM()
C SYMMET : MAXREP
C NUCLEI : NUCDEG()
C INFPAR : *
C HRUNIT : *
C PARINT : *
C
#include "gnrinf.h"
#include "inftap.h"
#include "siripc.h"
#include "inforb.h"
#include "molinp.h"
#include "cbirea.h"
#include "dorps.h"
#include "abainf.h"
#include "symmet.h"
#include "nuclei.h"
#include "dftcom.h"
#include "parint.h"
#include "blocks.h"
#include "r12int.h"
C
#include "cbisol.h"
C
Cef begin
#include "incore.h"
Cef end
      character(len=6) :: nmat_wo_win_env

      CALL QENTER('HER_SDINIT')
      IF (IPRINT .GT. 3) THEN
         CALL TITLER('Output from HER_SDINIT','*',103)
         WRITE(LUPRI,'(/A,I5)') 'Parallel type (ITYPE) is',ITYPE
      END IF
C
C     Calculate the number of Fock matrices that will be
C     returned from slaves to master.
C
      IF (ITYPE .EQ. 2) THEN
         NFMAT = 0
      ELSE IF (ITYPE .EQ. 3) THEN
         NFMAT = NDMAT
      ELSE IF (ITYPE .EQ. -5) THEN
         NFMAT = 3*NDMAT
#ifdef BUILD_OPENRSP
!radovan: this breaks dft_camb3lyp_magsus test on master
!         after the problem is fixed, remove the ifdef!
        IF (USE_OPENRSP) THEN
C added by Bin Gao, June 3, 2009
         IF (MAXDIF .EQ. 2 .AND. .NOT. GENCTR)
     &      NFMAT = 6*NDMAT
        END IF
#endif
      ELSE IF (ITYPE .EQ. 8) THEN
         NFMAT = 3*NDMAT*NUCDEG(IATOM)
      ELSE
         WRITE (LUPRI,'(//1X,A,/1X,A,I2)')
     &    'ERROR: specified ITYPE for TWOINT not defined in HER_SDINIT',
     &    'ITYPE =',ITYPE
         CALL QUIT('Specified ITYPE not defined in HER_SDINIT.')
      END IF
C
C----------------------------
C     Set common-block PARINT
C----------------------------
C
      JATOM  = IATOM
      JLUDAS = LUDASP
      JLUINT = LUINTA
      JLUONE = LUONEL
      JLUSOL = LUSOL
      JLUSUP = LUSUPM
      JMXDIF = MAXDIF
      JMXREP = MAXREP
      JNFMAT = NFMAT
      JTASK  = NTASK
      JTYPE  = ITYPE
CTROND: Modify when merging with DIRAC...
CTROND      J2TYP  = I2TYP
      J2TYP  = 0
      JCEDIF = ICEDIF
      JFTHRS = IFTHRS
C     logicals
      JNODV  = NODV
      JNOPV  = NOPV
      JNOCNT = NOCONT
      JRETUR = RETUR
      JTKTIM = TKTIME
      JSOLVN = SOLVNT
      JRELCL = RELCAL
C
#if defined (VAR_MPI)
      CALL MPIXBCAST(FIRST ,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWBAS,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWGEO,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(HFXFAC,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(HFXMU, 1,'DOUBLE',MASTER)
      CALL MPIXBCAST(DOSRIN,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(SRINTS,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(ERFEXP(0),3,'LOGICAL',MASTER)
      CALL MPIXBCAST(CHIVAL,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(EXPGRA,1,'LOGICAL',MASTER)
#ifdef BUILD_OPENRSP
      CALL MPIXBCAST(USE_OPENRSP,1,'LOGICAL',MASTER)
#endif
C
Cef begin
      IF (AOSAVE .AND. NEWGEO .AND. (.NOT. FIRST)) THEN
C     IF (AOSAVE .AND. NEWGEO) THEN
C         CALL IZERO(INDX_SHL,MXTSK)
         LINTSV = .FALSE.
         LINTMP = .FALSE.
         INITX = .FALSE.
         MSAVE = .TRUE.
         MMCORE = MMWORK
         LMCORE = MMCORE
         ISCORE = 1
         JSCORE = ISCORE
         N_SHL = 1
         I_SHL = 1
         INDX_SHL1 = 0
         INDX_SHL2 = 0
         INDX_SHL3 = 0
         INDX_SHL4 = 0
C         CALL CLEAR_INCOREMEM()
      END IF
Cef end
      IF (FIRST .OR. NEWBAS .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/A/A,3L10)')
     &        ' Sending new orbital information',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C  /molinp/
         CALL MPIXBCAST(NMLINE,1,'INTEGER',MASTER)
C Arnfinn, nov -08: Sending NONTYP_QM, numbers of atomtypes except MM
         CALL MPIXBCAST(NONTYP_QM,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NMLAU,1,'INTEGER',MASTER)
C  /molinc/
         CALL MPIXBCAST(MLINE,len_MLINE*NMLINE,'STRING',MASTER)
C  /cbirea/
         CALL MPIXBCAST(LCMMAX,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMSTR,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMEND,1,'INTEGER',MASTER)
         CALL MPIXBCAST(LCNTNUUM,1,'LOGICAL',MASTER)
         CALL MPIXBCAST(UNCONT,1,'LOGICAL',MASTER)
C  /inforb/
         CALL MPIXBCAST(MULD2H,NINFI,'INTEGER',MASTER)
C  /infpar/
         CALL MPIXBCAST(NODTOT,NPARI,'INTEGER',MASTER)
C  /dorps/ DOREPS, DOCOOR (NDORL covers both)
         CALL MPIXBCAST(DOREPS(0),NDORL,'LOGICAL',MASTER)
C  dosym()
         CALL MPIXBCAST(DOSYM,NSYML,'LOGICAL',MASTER)
C  comr12
         CALL MPIXBCAST(R12INT,NR12L,'LOGICAL',MASTER)
         CALL MPIXBCAST(LMULBS,1,'LOGICAL',MASTER)
         CALL MPIXBCAST(INTGAC,NR12I,'INTEGER',MASTER)
      END IF
C
C  /parint/
      CALL MPIXBCAST(JATOM,NINTI,'INTEGER',MASTER)
      CALL MPIXBCAST(JNODV,NINTL,'LOGICAL',MASTER)

!     rma_model
      CALL MPIXBCAST(rma_model,1,'LOGICAL',MASTER)
C
!     update info about type of fock-matrix construction
      CALL MPIXBCAST(ISYMDM,NFMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(IFCTYP,NFMAT,'INTEGER',MASTER)
!     distribute density matrices - that can be huge chunks... # first
      CALL MPIXBCAST(NDMAT,1,'INTEGER',MASTER)

      if(rma_model)then

        ! set # of fmat/dmat outside the RMA window ==> min = 1 otherwise we introduce some stupid communication
        !                                                       overhead in regular SCF calculations...
        nmat_wo_win_env = ' '
        call getenv('NO_RMA_MTX',nmat_wo_win_env)
        read(nmat_wo_win_env,'(i6)') rma_win_info%nmat_max_wo_win

!       sanity checks: minimum == 1; maximum == ndmat
        if(rma_win_info%nmat_max_wo_win <= 0 )
     &     rma_win_info%nmat_max_wo_win  = 1
        if(rma_win_info%nmat_max_wo_win  > ndmat)
     &     rma_win_info%nmat_max_wo_win  = ndmat

        CALL MPIXBCAST(rma_win_info%nmat_max_wo_win,1,'INTEGER',MASTER)

#ifdef VAR_MPI
!       distribute density matrix(ces) among node masters
        call mpi_bcast(dmat,n2basx*ndmat,mpi_real8,master,
     &                 communication_info_mpi%communication_internode,
     &                 ierr)
!       distribute density matrix(ces) among node/NUMA node slaves
        call mpi_bcast(dmat,
     &                 n2basx*rma_win_info%nmat_max_wo_win,
     &                 mpi_real8,
     &                 0,
     &                 communication_info_mpi%communication_shmemnode,
     &                 ierr)
#endif

      else

!       write(lupri,*) 'master dmat',nbast
!       call output(dMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
!       call flshfo(lupri)
        CALL MPIXBCAST(DMAT,N2BASX*NDMAT,'DOUBLE',MASTER)
      end if
C
C added by Bin Gao, June 3, 2009
      CALL MPIXBCAST(GENCTR,1,'LOGICAL',MASTER)
C
      IF(ITYPE.EQ.3) THEN
        N2RED = NSYMBL*NSYMBL
        CALL MPIXBCAST(N2RED,1,'INTEGER',MASTER)
        CALL MPIXBCAST(DMRAO,N2RED*NDMAT,'DOUBLE',MASTER)
        CALL MPIXBCAST(GABRAO,N2RED,'DOUBLE',MASTER)
      ENDIF
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(FIRST ,1,'LOGICAL')
      CALL PVMXPACK(NEWBAS,1,'LOGICAL')
      CALL PVMXPACK(NEWGEO,1,'LOGICAL')
      CALL PVMXPACK(BASDIR,60,'STRING')
      CALL PVMXPACK(LENBAS,1,'INTEGER')
C
      IF (FIRST .OR. NEWBAS .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Sending new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C  /molinc/
         CALL PVMXPACK(NMLINE,1,'INTEGER')
C  /molinp/
         CALL PVMXPACK(MLINE,len_MLINE*NMLINE,'STRING')
C  /inforb/
         CALL PVMXPACK(MULD2H,NINFI,'INTEGER')
C  /infpar/
         CALL PVMXPACK(NODTOT,NPARI,'INTEGER')
C  /dorps/
         CALL PVMXPACK(DOREPS(0),NDORL,'LOGICAL')
C  dosym()
         CALL PVMXPACK(DOSYM,NSYML,'LOGICAL')
C  comr12
         CALL PVMXBCAST(R12INT,NR12L,'LOGICAL')
         CALL PVMXBCAST(INTGAC,NR12I,'INTEGER')
      END IF
C  /parint/
      CALL PVMXPACK(JATOM,NINTI,'INTEGER')
      CALL PVMXPACK(JNODV,NINTL,'LOGICAL')
C
C     Pack the density matrices, their symmetries and the Fock matrix types.
C
      call quit('pvm version of Dalton SCF not working - '//
     &          'ask s. knecht at stefan.knecht@phys.chem.ethz.ch '//
     &          'for help.')
      CALL PVMXPACK(NDMAT,1,'INTEGER')
      CALL PVMXPACK(DMAT,N2BASX*NDMAT,'DOUBLE')
      CALL PVMXPACK(ISYMDM,NFMAT,'INTEGER')
      CALL PVMXPACK(IFCTYP,NFMAT,'INTEGER')
C
      IF(ITYPE.EQ.3) THEN
        N2RED = NSYMBL*NSYMBL
        CALL PVMXPACK(N2RED,1,'INTEGER')
        CALL PVMXPACK(DMRAO,N2RED*NDMAT,'DOUBLE')
        CALL PVMXPACK(GABRAO,N2RED,'DOUBLE')
      ENDIF
C
C     -----------------------------------
C     Multicast initialization to nodes.
C     -----------------------------------
C
      CALL PVMXMCAST(NODTOT,NODEID(1),MTAG3)
C
#endif
      CALL QEXIT('HER_SDINIT')
      RETURN
      END
C  /* Deck her_sendreceive_molinfo */
      SUBROUTINE HER_sendreceive_molinfo(IPRINT,WORK,KFREE,LFREE)
!
!     Hans Joergen Aa. Jensen, Nov. 2016
!
!     New routine (2016) called by both master and slaves to transfer
!     molecule information from master to slaves.
!

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
      use mpi
#include "implicit.h"
#else
#include "implicit.h"
#include "mpif.h"
#endif
#else
#include "implicit.h"
#endif
C
      INTEGER IPRINT, KFREE, LFREE
      REAL*8  WORK(*)

#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mtags.h"
C
C Used from common blocks:
C
C MXCENT : MXCENT (for DORPS, etc.)
C MXORB  : MXSHEL (for INFPAR)
C MAXAQN : *      (for SYMMET)
C
C INFPAR : MASTER, MYNUM, ?
C
C gnrinf : a lot
C SIRIPC : NEWGEO
C INFORB : MULD2H, ...
C MOLINP : NMLINE, MLINE()
C CBIREA : UNCONT, LMULBS
C DORPS  : DOREPS(), DOCOOR()
C ABAINF : DOSYM()
C SYMMET : MAXREP
C DFTCOM : HFXFAC, HFXMU, ...
C EXPOPT : EXPGRA
C R12INT : R12INT, INTGAC
C
#include "infpar.h"
#include "gnrinf.h"
#include "siripc.h"
#include "inforb.h"
#include "molinp.h"
#include "cbirea.h"
#include "dorps.h"
#include "abainf.h"
#include "symmet.h"
#include "dftcom.h"
#include "expopt.h"
#include "r12int.h"

      CALL QENTER('HER_sendreceive_molinfo')
      IF (IPRINT .GT. 3) THEN
         CALL TITLER('Output from HER_sendreceive_molinfo','*',103)
         WRITE(LUPRI,'(A,I0)') ' From node number ',MYNUM
      END IF

#if defined (VAR_MPI)
      CALL MPIXBCAST(NEWBAS,1,'LOGICAL',MASTER)! gnrinf.h
      CALL MPIXBCAST(NEWGEO,1,'LOGICAL',MASTER)! siripc.h
      CALL MPIXBCAST(HFXFAC,1,'DOUBLE',MASTER) ! dftcom.h
      CALL MPIXBCAST(HFXMU, 1,'DOUBLE',MASTER) ! dftcom.h
      CALL MPIXBCAST(DOSRIN,1,'LOGICAL',MASTER)! gnrinf.h
      CALL MPIXBCAST(SRINTS,1,'LOGICAL',MASTER)! gnrinf.h
      CALL MPIXBCAST(ERFEXP,1,'LOGICAL',MASTER)! gnrinf.h
      CALL MPIXBCAST(CHIVAL,1,'DOUBLE',MASTER) ! gnrinf.h
      CALL MPIXBCAST(EXPGRA,1,'LOGICAL',MASTER)! expopt.h
C
! hjaaj Nov 2016: always send basis set information
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(A,2L10)')
     &        ' Sending/receiving orbital information.'//
     &        ' NEWBAS,NEWGEO :',NEWBAS,NEWGEO
C  /molinp/
         CALL MPIXBCAST(NMLINE,1,'INTEGER',MASTER)
C Arnfinn, nov -08: Sending NONTYP_QM, numbers of atomtypes except MM
         CALL MPIXBCAST(NONTYP_QM,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NMLAU,1,'INTEGER',MASTER)
C  /molinc/
         CALL MPIXBCAST(MLINE,len_MLINE*NMLINE,'STRING',MASTER)
C  /cbirea/
         CALL MPIXBCAST(LCMMAX,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMSTR,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMEND,1,'INTEGER',MASTER)
         CALL MPIXBCAST(LCNTNUUM,1,'LOGICAL',MASTER)
         CALL MPIXBCAST(UNCONT,1,'LOGICAL',MASTER)
         CALL MPIXBCAST(LMULBS,1,'LOGICAL',MASTER)     ! cbirea.h
C  /inforb/
         CALL MPIXBCAST(MULD2H,NINFI,'INTEGER',MASTER)
C  /infpar/
         CALL MPIXBCAST(NODTOT,NPARI,'INTEGER',MASTER)
C  /dorps/ DOREPS, DOCOOR (NDORL covers both)
         CALL MPIXBCAST(DOREPS(0),NDORL,'LOGICAL',MASTER)
C  dosym()
         CALL MPIXBCAST(DOSYM,NSYML,'LOGICAL',MASTER)
C  comr12
         CALL MPIXBCAST(R12INT,NR12L,'LOGICAL',MASTER) ! r12int.h
         CALL MPIXBCAST(INTGAC,NR12I,'INTEGER',MASTER) ! r12int.h
!     END IF

      IF (MYNUM .NE. MASTER) THEN ! I am a slave
C
C        Since REAINI is not initialized outside of SETHER,
C        make sure slaves initialize it properly
C        (REAINI will return immediately if already called).
C
         RELCAL = .FALSE.
         TSTINP = .FALSE.
         CALL REAINI(0,RELCAL,TSTINP)
C
C        Other defaults for the slaves
C
         USE_OPENRSP = .FALSE.
C
C        Set basis set information always (with NEWGEO .true.)
C
         CALL SETHER(0,.true.,WORK(KFREE),LFREE)
!        CALL SETHER(0,NEWGEO,WORK(KFREE),LFREE)
C
      END IF
C
#else
      CALL QUIT('ERROR: '//
     &   'HER_sendreceive_molinfo called for sequential Dalton')
#endif
      CALL QEXIT('HER_sendreceive_molinfo')
      RETURN
      END
C  /* Deck her_rvres */
      SUBROUTINE HER_RVRES(FMAT,TMPMAT,HESSEE,TIMES,TSKCPU,IWHICH,
     &                     L_TMPMT,L_GRA,L_HES,L_SUS,L_FCK,ITYPE,
     &                     IPRINT,DINTSKP,GENCTR)
!     module dependencies
      use rma_windows
#ifdef VAR_MPI
      use parallel_communication_models_mpi
      use one_sided_communication_wrappers
#ifdef USE_MPI_MOD_F90
      use mpi
#endif
#endif
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "mtags.h"
      PARAMETER (D1 = 1.0D0)
C added by Bin Gao, June 3, 2009
      logical GENCTR
      DIMENSION FMAT(*), TMPMAT(L_TMPMT), TIMES(MTOTTK),
     &          TSKCPU(MTOTTK), IWHICH(MTOTTK),DSKPBF(8),
     &          DINTSKP(8), HESSEE(*)
C
C Used from common blocks:
C
C MXCENT : MXCENT (for ENERGY)
C MXORB  : MXSHEL (for INFPAR)
C ENERGY : GRADEE,HESSEE
C SUSCPT : SUS2EL
C INFORB : NNBASX,N2BASX
C INFPAR : TIMING, ???
C
#include "energy.h"
#include "suscpt.h"
#include "inforb.h"
#include "infpar.h"
#include "parint.h"
#include "expopt.h"
#include "gnrinf.h"

#ifdef VAR_MPI
#ifndef USE_MPI_MOD_F90
#include "mpif.h"
#endif
       integer (kind=mpi_offset_kind) :: offset_mat
#endif
C
      CALL QENTER('HER_RVRES')
      IF (IPRINT .GT. 4) THEN
         CALL TITLER('Output from HER_RVRES','*',103)
         WRITE(LUPRI,'(/A,I5)') 'Parallel type (ITYPE) is',ITYPE
      END IF
C
C     Clean matrices.
C
      IF (TIMING) CALL DZERO(TIMES,MTOTTK)
C
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,L_GRA)
C modified by Bin Gao, December 24, 2009
         CALL DZERO(HESSEE,L_HES)
      ELSE IF (ITYPE .EQ. 3) THEN
         if(.not.rma_model)then
           CALL DZERO(FMAT,L_FCK)
           CALL DZERO(DINTSKP,8)
         end if
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
         CALL BIN_SUM(L_FCK, FMAT)
         CALL BIN_SUM(8, DINTSKP)
#endif
      ELSE IF (ITYPE .EQ. -5) THEN
#ifdef BUILD_OPENRSP
!radovan: this breaks dft_camb3lyp_magsus test on master
!         after the problem is fixed, remove the ifdef!

C modified by Bin Gao, December 21, 2009
        IF (USE_OPENRSP) THEN
         CALL DZERO(FMAT,L_FCK)
         CALL DZERO(SUS2EL,L_SUS)
        ELSE
#endif
         CALL DZERO(FMAT,3*N2BASX)
#ifdef BUILD_OPENRSP
        END IF
#endif
COBL/AMT
C for CAMB3LYP, we must not zero the susceptibilities the second run in TWOLOP
         IF (.NOT. DONETWOLOP)THEN
             DONETWOLOP = .TRUE.
             CALL DZERO(SUS2EL,L_SUS)
         END IF
      ELSE IF (ITYPE .EQ. 8) THEN
         CALL DZERO(FMAT,L_TMPMT)
      END IF
C
C---------------------------
C     Start loop over nodes.
C---------------------------
C
      ICOUNT = 0
      DO 100 I = 1,NODTOT
C
#if defined (VAR_MPI)
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MTAG7)
         IF (ITYPE .EQ. 2) THEN
            CALL MPIXRECV(TMPMAT(1),L_GRA,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(L_GRA,D1,TMPMAT,1,GRADEE,1)
            CALL MPIXRECV(TMPMAT(1),L_HES,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(L_HES,D1,TMPMAT(1),1,HESSEE,1)
            CALL MPIXRECV(TMPMAT(1),MXPRIM,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(MXPRIM,D1,TMPMAT(1),1,ALPGRD,1)
#if !defined(_CRAYT3E) && !defined(NO_BINSUM)
         ELSE IF (ITYPE .EQ. 3) THEN
            if(.not.rma_model)then
              CALL MPIXRECV(TMPMAT,L_FCK,'DOUBLE',NWHO,MTAG7)
              CALL MPIXRECV(TMPMAT(1+L_FCK),8,'DOUBLE',NWHO,MTAG7)
              CALL DAXPY(L_FCK,D1,TMPMAT,1,FMAT,1)
              CALL DAXPY(8,D1,TMPMAT(L_FCK+1),1,DINTSKP,1)
            end if
#endif
         ELSE IF (ITYPE .EQ. 8) THEN
            CALL MPIXRECV(TMPMAT(1),L_TMPMT,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(L_TMPMT,D1,TMPMAT,1,FMAT,1)
         ELSE IF (ITYPE .EQ. -5) THEN
#ifdef BUILD_OPENRSP
!radovan: this breaks dft_camb3lyp_magsus test on master
!         after the problem is fixed, remove the ifdef!
           IF (USE_OPENRSP) THEN
            CALL MPIXRECV(TMPMAT(1),L_FCK,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(L_FCK,D1,TMPMAT,1,FMAT,1)
           ELSE
#endif
            CALL MPIXRECV(TMPMAT(1),3*N2BASX,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(3*N2BASX,D1,TMPMAT,1,FMAT,1)
#ifdef BUILD_OPENRSP
           END IF
#endif
            CALL MPIXRECV(TMPMAT(1),L_SUS,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(L_SUS,D1,TMPMAT,1,SUS2EL,1)
         END IF
C
         CALL MPIXRECV(ITASKS,1,'INTEGER',NWHO,MTAG7)
C
         IF (TIMING) THEN
            CALL MPIXRECV(IWHICH(ICOUNT+1),ITASKS,'INTEGER',NWHO,MTAG7)
            CALL MPIXRECV(TSKCPU(ICOUNT+1),ITASKS,'DOUBLE' ,NWHO,MTAG7)
            DO 200 J=1,ITASKS
               K = IWHICH(J)
               IF (K .GT. 0) TIMES(K) = TSKCPU(J)
 200        CONTINUE
         END IF
C
         ICOUNT = ICOUNT + ITASKS
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG7)
         CALL PVMXUNPACK(TMPMAT,L_TMPMT,'DOUBLE')
C
         IF (ITYPE .EQ. 2) THEN
            CALL DAXPY(L_GRA,D1,TMPMAT(1),1,GRADEE,1)
            CALL DAXPY(L_HES,D1,TMPMAT(L_GRA+1),1,HESSEE,1)
C
         ELSE IF (ITYPE .EQ. 3) THEN
            CALL DAXPY(L_FCK,D1,TMPMAT,1,FMAT,1)
            CALL DAXPY(8,D1,TMPMAT(L_FCK+1),1,DINTSKP,1)
C
         ELSE IF (ITYPE .EQ. 8) THEN
            CALL DAXPY(L_TMPMT,D1,TMPMAT(1),1,FMAT,1)
         ELSE IF (ITYPE .EQ. -5) THEN
            CALL DAXPY(3*N2BASX,D1,TMPMAT(1),1,FMAT,1)
            CALL DAXPY(L_SUS,D1,TMPMAT(3*N2BASX+1),1,SUS2EL,1)
         END IF
C
C        Unpack the number of tasks received by the node.
C
         CALL PVMXUNPACK(ITASKS,1,'INTEGER')
C
         IF (TIMING) THEN
            CALL PVMXUNPACK(IWHICH(ICOUNT+1),ITASKS,'INTEGER')
            CALL PVMXUNPACK(TSKCPU(ICOUNT+1),ITASKS,'DOUBLE')
C
            DO 200 J=1,ITASKS
               K = IWHICH(J)
               IF (K .GT. 0) TIMES(K) = TSKCPU(J)
 200        CONTINUE
         END IF
C
         ICOUNT = ICOUNT + ITASKS
#endif
  100 CONTINUE

!     collect data for rma_model (only itask == 3 at present)
      if(rma_model)then

#ifdef VAR_MPI
        call mpi_reduce(mpi_in_place,
     &                  dintskp,
     &                  8,
     &                  mpi_real8,
     &                  mpi_sum,
     &                  master,
     &                  communication_info_mpi%
     &                  communication_glb_world,
     &                  ierr)

!       all fmatS outside the RMA window...
        call mpi_reduce(
     &                  mpi_in_place,
     &                  fmat,
     &                  rma_win_info%nmat_max_wo_win * n2basx,
     &                  mpi_real8,
     &                  mpi_sum,
     &                  master,
     &                  communication_info_mpi%
     &                  communication_glb_world,
     &                  ierr)

!       now all fmatS inside the RMA window... if any
        if(rma_win_info%nmat_max_wo_win < jnfmat )then
          call mpi_barrier(communication_info_mpi%
     &                     communication_shmemnode,
     &                     ierr)
          call mpi_reduce(
     &                    mpi_in_place,
     &                    fmat(1 + rma_win_info%
     &                    nmat_max_wo_win * n2basx),
     &                    (jnfmat - rma_win_info%
     &                    nmat_max_wo_win) * n2basx,
     &                    mpi_real8,
     &                    mpi_sum,
     &                    master,
     &                    communication_info_mpi%
     &                    communication_internode,
     &                    ierr)
        end if
#endif
      end if
C
      IF (ICOUNT .NE. MTOTTK) THEN
        WRITE(LUPRI,'(/5X,A)') 'Error in parallel calculation!'
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks to calculate  : ',
     &       MTOTTK
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks    calculated : ',
     &       ICOUNT
c        CALL QUIT('Inconsistence in HER_RVRES -> MTOTTK .NE. ICOUNT')
      ENDIF
C
      CALL QEXIT('HER_RVRES')
      RETURN
      END

C  /* Deck parsrt */
      SUBROUTINE PARSRT(MTOTTK,NPOS,SORTED,TIMES,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION  TIMES(MTOTTK), NPOS(0:MTOTTK+1), SORTED(0:MTOTTK+1)
C
      IF (IPRINT .GT. 4) CALL TITLER('Output from PARSRT','*',103)
C
C------------------------------------------------------
C     The array TIMES containing CPU-times is sorted by
C     increasing CPUs and stored in the array SORTED.
C------------------------------------------------------
C
C     SORTED(0) must be less than any possible CPU-time,
C     and TIMMAX must be greater than any possible CPU-time
C
      ITEMS     = MTOTTK
      TIMMAX    = 99999D00
      SORTED(0) =  -1.0D00
      SORTED(1) = TIMES(1)
      NPOS  (0) = 0
      NPOS  (1) = 1
C
      DO 100,ITEM = 2,ITEMS
C
         LASTLO       = 0
         LASTHI       = ITEM
         NUMBER       = INT(ITEM/2)
         SORTED(ITEM) = TIMMAX
C
 200     CONTINUE
C
         IF (TIMES(ITEM) .EQ. SORTED(NUMBER) .OR.
     &      (TIMES(ITEM) .GT. SORTED(NUMBER) .AND.
     &       TIMES(ITEM) .LT. SORTED(NUMBER+1))) THEN
C
            DO 300, NCOUNT = ITEM,NUMBER,-1
               NPOS  (NCOUNT+1) = NPOS  (NCOUNT)
               SORTED(NCOUNT+1) = SORTED(NCOUNT)
 300        CONTINUE
C
            NPOS  (NUMBER+1) = ITEM
            SORTED(NUMBER+1) = TIMES(ITEM)
            GO TO 100
C
         ELSE IF (TIMES(ITEM) .LT. SORTED(NUMBER)) THEN
C
            LASTHI = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         ELSE IF (TIMES(ITEM) .GT. SORTED(NUMBER)) THEN
C
            LASTLO = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         END IF
 100  CONTINUE
C
      IF (IPRINT .GT. 4) THEN
         DO 400 I = 1, ITEMS
            WRITE(LUPRI,'(15X,A4,I5,A9,F8.3)')
     &           'IJ =', NPOS(I), '   time =', SORTED(I)
 400      CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck her_nodstr */
      SUBROUTINE HER_NODSTR(LWORK_save,wrkdlm,IPRINT)
C
C    *****************************************************************
C    *    This is the node program for the construction of fock-     *
C    *   matrices, derivated fock-matrices and expectation values.   *
C    *****************************************************************
C
!     module dependencies
      use memory_parallel
      use parallel_communication_models_mpi
      use rma_windows
#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
      use mpi
#endif
#endif

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"
C
      PARAMETER (MXDMAT = 2*MXCOOR)
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,GENCTR
      DIMENSION ISYMDM(MXDMAT), IFCTYP(MXDMAT),
     &          DINTSKP(2,4)
      real(8), allocatable :: work(:)
      integer, intent(in)  :: lwork_save
      integer              :: lwork
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C MXCENT : MXCOOR (for ENERGY)
C ENERGY : GRADEE, HESSEE
C SUSCPT : SUS2EL
C INFORB : N2BASX
C INFPAR : SLAVE
C
#include "energy.h"
#include "suscpt.h"
#include "inforb.h"
#include "infpar.h"
#include "gnrinf.h"
C
#ifdef VAR_MPI
#ifndef USE_MPI_MOD_F90
#include "mpif.h"
#endif
      integer, parameter             :: dummy_buff = 2
      integer(kind=MPI_ADDRESS_KIND) :: nelement
      integer(kind=MPI_OFFSET_KIND)  :: displacement
      integer(kind=MPI_ADDRESS_KIND) :: nelement_d, nelement_f
      integer(kind=MPI_ADDRESS_KIND) :: dmat_p, fmat_p
      real(8)                        :: dmat_buff(dummy_buff)
      real(8)                        :: fmat_buff(dummy_buff)
      pointer                        (dmat_p, dmat_buff)
      pointer                        (fmat_p, fmat_buff)
      character(len=16)              :: shared_dmat_file
      integer                        :: istat(mpi_status_size)
#endif
      integer, allocatable           :: JSTRS(:)
      integer, allocatable           :: NPRIM(:)
      integer, allocatable           :: NCONT(:)
      integer, allocatable           :: IORBS(:)
      integer, allocatable           :: JORBS(:)
      integer, allocatable           :: KORBS(:)

C     This is a slave in a parallel run; check for programming error
      IF (.NOT.SLAVE) THEN
         CALL QUIT('ERROR: in HER_NODSTR but not SLAVE')
      END IF

      CALL QENTER('HER_NODSTR')

      IF (IPRINT .GE. 2) THEN
         CALL TITLER('Output from HER_NODSTR','*',103)
      END IF

!     Start timing
      CALL GETTIM(CPU1,WALL1)

!     allocate local variables
      NPAO   = MXSHEL*MXAOVC
      allocate(JSTRS(NPAO*2))
      allocate(NPRIM(NPAO*2))
      allocate(NCONT(NPAO*2))
      allocate(IORBS(NPAO  ))
      allocate(JORBS(NPAO  ))
      allocate(KORBS(NPAO  ))

!     set local work memory size
      lwork = lwork_save - ( 9 * npao / (8/sizeof(int)))

      allocate(work(0:lwork+1),stat=i)
!     Set memory traps
      work(0)       = wrkdlm
      work(1+lwork) = wrkdlm
!     set WORK pointers
      kwork = 1
      kfree = kwork
      lfree = lwork

      if(i /= 0)then
        write (luerr,*) mynum,
     &      ': ALLOCATE command failed to allocate'//
     &      ' the requested memory. Error code:',i
        write (luerr,*)
     &      'Reduce the memory demands and be welcome back'
        call quit('Failed to allocate memory')
      end if

!     Receive initialization from master - 1 - common block information
!     =================================================================
!
      CALL HER_RVINIT(WORK(1),KFREE,LFREE,
     &                NDMAT,ISYMDM,IFCTYP,ITYPE,
     &                IATOM,MAXDIF,NODV,NOPV,NOCONT,
     &                RETUR,TKTIME,JSTRS,NPRIM,NCONT,
     &                IORBS,JORBS,KORBS,IPRINT,
     &                I2TYP,ICEDIF,IFTHRS,
     &                N2RED,dmat_buff,fmat_buff,KHESEE,KDMRAO,
     &                KGABAO,KDMRSO,GENCTR,.true.)
      deallocate(work)

!     Receive initialization from master - 2 - data on allocated matrices
!     ===================================================================
!
      nelement_d     = n2basx*ndmat
      nelement_f     = n2basx*nfmat
      if(rma_model .and.
     &   communication_info_mpi%my_shmem_node_id /= 0)then
          nelement_d = n2basx*rma_win_info%nmat_max_wo_win
          nelement_f = n2basx*rma_win_info%nmat_max_wo_win
      end if
#ifdef VAR_MPI
      call memallocmpi(nelement_d*8,dmat_p)
      call memallocmpi(nelement_f*8,fmat_p)
#endif
      call dzero(fmat_buff,nelement_f)

!     set local work memory size
      lwork = lwork_save - ( 9 * npao / (8/sizeof(int))) -
     &        nelement_d - nelement_f

      allocate(work(0:lwork+1),stat=i)
!     Set memory traps
      work(0)       = wrkdlm
      work(1+lwork) = wrkdlm
!     set WORK pointers
      kwork = 1
      kfree = kwork
      lfree = lwork

      CALL HER_RVINIT(WORK(1),KFREE,LFREE,
     &                NDMAT,ISYMDM,IFCTYP,ITYPE,IATOM,
     &                MAXDIF,NODV,NOPV,NOCONT,RETUR,TKTIME,
     &                JSTRS,NPRIM,NCONT,
     &                IORBS,JORBS,KORBS,IPRINT,
     &                I2TYP,ICEDIF,IFTHRS,
     &                N2RED,dmat_buff,fmat_buff,KHESEE,KDMRAO,
     &                KGABAO,KDMRSO,GENCTR,.false.)
      IF (IPRINT .GE. 1) THEN
         WRITE(LUPRI,'(/A,4I5)')
     &   'After HER_RVINIT; ITYPE, IATOM, NDMAT, MAXDIF =',
     &   ITYPE,IATOM,NDMAT,MAXDIF
      END IF
C
C     Initialize matrices
C     ===================
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,MXCOOR)
         CALL DZERO(WORK(KHESEE),MXCOOR*MXCOOR)
      ELSE IF (ITYPE .EQ. -5) THEN
         CALL DZERO(SUS2EL,9)
      END IF

!--------------------------------------------------------------------------------------
!     initialize memory windows for rma access within a given shared-memory / NUMA node
!--------------------------------------------------------------------------------------
      if(rma_model)then

        rma_win_info%dmat_win = -1
        rma_win_info%fmat_win = -1


#ifdef VAR_MPI
        nelement            = n2basx*ndmat
        if(communication_info_mpi%my_shmem_node_id /= 0) nelement = 0

!       call set_rma_window(
!    &                      dmat_buff,
!    &                      nelement,
!    &                      communication_info_mpi%my_shmem_node_id,
!    &                      communication_info_mpi%
!    &                      communication_shmemnode,
!    &                      rma_win_info%dmat_win,
!!   &                      .true.,
!    &                      .false.,
!    &                      'no_locks',
!    &                      'true')

!       test dmat on shared file
        write(shared_dmat_file,'(a12,i4)')
     &  'dmat_shared.',communication_info_mpi%intra_node_group_id
        call mpi_file_open(
     &                     communication_info_mpi%
     &                     communication_intranode,
     &                     shared_dmat_file,
     &                     mpi_mode_CREATE +
     &                     mpi_mode_RDWR   +
     &                     mpi_mode_DELETE_ON_CLOSE,
     &                     mpi_info_null,
     &                     rma_win_info%dmat_fh,
     &                     ierr )
        displacement = 0
        call mpi_file_set_view(
     &                         rma_win_info%dmat_fh,
     &                         displacement,
     &                         mpi_real8,
     &                         mpi_real8,
     &                         'native',
     &                         mpi_info_null,
     &                         ierr )
        if(communication_info_mpi%my_intra_node_id == 0)then
          nelement4 = nelement
          call mpi_file_write(
     &                        rma_win_info%dmat_fh,
     &                        dmat_buff,
     &                        nelement4,
     &                        mpi_real8,
     &                        istat,
     &                        ierr )
        end if
        call mpi_barrier(communication_info_mpi%
     &                   communication_intranode,ierr)
!       end test dmat on shared file

        nelement            = n2basx*nfmat
        if(communication_info_mpi%my_shmem_node_id /= 0) nelement = 0

        call set_rma_window(
     &                      fmat_buff,
     &                      nelement,
     &                      communication_info_mpi%my_shmem_node_id,
     &                      communication_info_mpi%
     &                      communication_shmemnode,
     &                      rma_win_info%fmat_win,
!    &                      .true.,
     &                      .false.,
     &                      'no_locks',
     &                      'true')
#endif
!       print *, 'fmat_win, dmat_win',rma_win_info%fmat_win,
!    &            rma_win_info%dmat_win,mynum

      end if

C
C     ***********************************
C     ***** Calculate the integrals *****
C     ***********************************
C
      IF (IPRINT .GE. 1) THEN
         write(lupri,*) 'HER_NODSTR: Calling TWOINT, NDMAT=',NDMAT
         flush(lupri)
      END IF
      GMAT   = 0
      INDXAB = 0
      NUMDIS = 0
      MAXDIS = 1
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      ISHLA  = 0
C
      CALL TWOINT(WORK(KFREE),LFREE,
     &     WORK(KHESEE),FMAT_buff,DMAT_buff,NDMAT,
     &     ISYMDM,IFCTYP,GMAT,INDXAB,NUMDIS,MAXDIS,ITYPE,
     &     MAXDIF,IATOM,NODV,NOPV,NOCONT,TKTIME,IPRINT,
     &     IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,ISHLA,I2TYP,
     &     JSTRS,NPRIM,NCONT,IORBS,
     &     ICEDIF,IFTHRS,WORK(KGABAO),WORK(KDMRAO),
     &     WORK(KDMRSO),DINTSKP,RELCALC,GENCTR)

!------------------------
!     free memory windows
!------------------------
      if(rma_model)then
#ifdef VAR_MPI
!       call free_rma_window(rma_win_info%dmat_win,master)
        call free_rma_window(rma_win_info%fmat_win,master)
!       test dmat on shared file
        call mpi_file_close(rma_win_info%dmat_fh,ierr)
!       end test dmat on shared file
#endif
      end if

      CALL MEMREL('HER_NODSTR',WORK(1),KWORK,KWORK,KFREE,LFREE)
C
C     Send overall timings.
C
      CALL SDTIM(CPU1,WALL1,IPRINT)
C
      deallocate(work)

#ifdef VAR_MPI
      call memfreempi(fmat_buff)
      call memfreempi(dmat_buff)
#endif

      deallocate(JSTRS)
      deallocate(NPRIM)
      deallocate(NCONT)
      deallocate(IORBS)
      deallocate(JORBS)
      deallocate(KORBS)
C
      IF (IPRINT .GE. 1) THEN
         write(lupri,*) 'Exit from HER_NODSTR'
         flush(lupri)
      END IF
      CALL QEXIT('HER_NODSTR')
      RETURN
      END
C  /* Deck her_rvinit */
      SUBROUTINE HER_RVINIT(WORK,KFREE,LFREE,
     &     NDMAT,ISYMDM,IFCTYP,ITYPE,
     &     IATOM,MAXDIF,NODV,NOPV,NOCONT,
     &     RETUR,TKTIME,JSTRSH,NPRIMS,NCONTS,IORBSH,
     &     JORBSH,KORBSH,IPRINT,I2TYP,ICEDIF,IFTHRS,
     &     N2RED,dmat_buff,fmat_buff,KHESEE,KDMRAO,
     &     KGABAO,KDMRSO,GENCTR,set_cb_blocks)
C

!     module dependencies
      use parallel_communication_models_mpi
      use rma_windows

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
      use mpi
#include "implicit.h"
#else
#include "implicit.h"
#include "mpif.h"
#endif
#else
#include "implicit.h"
#endif

#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mtags.h"
      PARAMETER (D0 = 0.0D0)
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,FIRST,NOSYM,GENCTR
      LOGICAL AUTOSY,DOCART,DOOWN
      CHARACTER*1 KASYM(3,3),ID3
C
      DIMENSION WORK(*), ISYMDM(*), IFCTYP(*), JSTRSH(*), NPRIMS(*),
     &          NCONTS(*), IORBSH(*), JORBSH(*), KORBSH(*)

      real(8), intent(inout)    :: fmat_buff(*), dmat_buff(*)
      logical, intent(in)       :: set_cb_blocks

C
#include "infpar.h"
C
C Used from common blocks:
C
C MXCENT : MXCENT (for DORPS, etc.)
C MXORB  : MXSHEL (for INFPAR)
C MAXAQN : *      (for SYMMET)
C GNRINF : NEWBAS, (PVM: BASDIR, LENBAS)
C SIRIPC : NEWGEO
C INFORB : N2BASX
C MOLINP : NMLINE, MLINE()
C CBIREA : UNCONT
C DORPS  : DOREPS(), DOCOOR()
C ABAINF : DOSYM()
C SYMMET : MAXREP
C NUCLEI : NUCDEG()
C INFPAR : *
C HRUNIT : *
C PARINT : *
C
#include "gnrinf.h"
#include "siripc.h"
#include "inforb.h"
#include "inftap.h"
#include "molinp.h"
#include "cbirea.h"
#include "dorps.h"
#include "abainf.h"
#include "symmet.h"
#include "nuclei.h"
#include "dftcom.h"
#include "parint.h"
#include "r12int.h"
#include "expopt.h"
C
#include "cbisol.h"
#include "incore.h"
C
      IF (IPRINT .GT. 2) THEN
         CALL TITLER('Output from HER_RVINIT','*',103)
         WRITE(LUPRI,*) 'Set Common blocks is ',set_cb_blocks
         WRITE(LUPRI,*) 'IPRINT            is ',IPRINT
      END IF

!     part - 1 - receive data on common blocks
!     ----------------------------------------
      if(set_cb_blocks)then

        call qenter('HER_RVINIT - 1 -')
C
C       Since REAINI is not initialized outside of SETHER, make sure slaves
C       initialize it properly (REAINI will return immediately if
C       already called).
C
        RELCAL = .FALSE.
        TSTINP = .FALSE.
        CALL REAINI(0,RELCAL,TSTINP)
C
#if defined (VAR_MPI)
        CALL MPIXBCAST(FIRST ,   1,'LOGICAL',MASTER)
        CALL MPIXBCAST(NEWBAS,   1,'LOGICAL',MASTER)
        CALL MPIXBCAST(NEWGEO,   1,'LOGICAL',MASTER)
        CALL MPIXBCAST(HFXFAC,   1,'DOUBLE' ,MASTER)
        CALL MPIXBCAST(HFXMU,    1,'DOUBLE' ,MASTER)
        CALL MPIXBCAST(DOSRIN,   1,'LOGICAL',MASTER)
        CALL MPIXBCAST(SRINTS,   1,'LOGICAL',MASTER)
        CALL MPIXBCAST(ERFEXP(0),3,'LOGICAL',MASTER)
        CALL MPIXBCAST(CHIVAL,   1,'DOUBLE' ,MASTER)
        CALL MPIXBCAST(EXPGRA,   1,'LOGICAL',MASTER)
#ifdef BUILD_OPENRSP
        CALL MPIXBCAST(USE_OPENRSP,1,'LOGICAL',MASTER)
#endif
C
C       IF (CHIVAL .EQ. -1.0D0) THEN calculate 1/R12
C                               ELSE calculate g_lr(R12)
C       IF (SRINTS) THEN  calculate 1/R12 - g_lr(R12)
C
C
Cef begin
        IF (AOSAVE .AND. NEWGEO .AND. (.NOT. FIRST)) THEN
C         IF (AOSAVE .AND. NEWGEO) THEN
C         CALL IZERO(INDX_SHL,MXTSK)
          LINTSV = .FALSE.
          LINTMP = .FALSE.
          INITX = .FALSE.
          MSAVE = .TRUE.
          MMCORE = MMWORK
          LMCORE = MMCORE
          ISCORE = 1
          JSCORE = ISCORE
          N_SHL = 1
          I_SHL = 1
          INDX_SHL1 = 0
          INDX_SHL2 = 0
          INDX_SHL3 = 0
          INDX_SHL4 = 0
C         CALL CLEAR_INCOREMEM()
        END IF
Cef end
        TIMING = FIRST .OR. NEWBAS .OR. NEWGEO
        IF (IPRINT .GE. 2) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &       ' Receive new orbital and geometry information ?',
     &       ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C
        IF (TIMING) THEN
C  /molinp/
              CALL MPIXBCAST(NMLINE,1,'INTEGER',MASTER)
C Arnfinn nov -08: Receiving number of QM atom types:
              CALL MPIXBCAST(NONTYP_QM,1,'INTEGER',MASTER)
              CALL MPIXBCAST(NMLAU,1,'INTEGER',MASTER)
C  /molinc/
              CALL MPIXBCAST(MLINE,len_MLINE*NMLINE,'STRING',MASTER)
C  /cbirea/
              CALL MPIXBCAST(LCMMAX,1,'INTEGER',MASTER)
              CALL MPIXBCAST(NCMSTR,1,'INTEGER',MASTER)
              CALL MPIXBCAST(NCMEND,1,'INTEGER',MASTER)
              CALL MPIXBCAST(LCNTNUUM,1,'LOGICAL',MASTER)
              CALL MPIXBCAST(UNCONT,1,'LOGICAL',MASTER)
C  /inforb/
              CALL MPIXBCAST(MULD2H,NINFI,'INTEGER',MASTER)
C  /infpar/
              CALL MPIXBCAST(NODTOT,NPARI,'INTEGER',MASTER)
C  /dorps/ DOREPS, DOCOOR (NDORL covers both)
              CALL MPIXBCAST(DOREPS(0),NDORL,'LOGICAL',MASTER)
C  dosym()
              CALL MPIXBCAST(DOSYM,NSYML,'LOGICAL',MASTER)
C  /comr12/
              CALL MPIXBCAST(R12INT,NR12L,'LOGICAL',MASTER)
              CALL MPIXBCAST(LMULBS,1,'LOGICAL',MASTER)
              CALL MPIXBCAST(INTGAC,NR12I,'INTEGER',MASTER)
        END IF
C  /parint/
        CALL MPIXBCAST(JATOM,NINTI,'INTEGER',MASTER)
        CALL MPIXBCAST(JNODV,NINTL,'LOGICAL',MASTER)

!       rma_model
        CALL MPIXBCAST(rma_model,1,'LOGICAL',MASTER)
#endif
#if defined (VAR_PVM)
        CALL PVMXRECV(MASTER,MTAG3)
        CALL PVMXUNPACK(FIRST ,1,'LOGICAL')
        CALL PVMXUNPACK(NEWBAS,1,'LOGICAL')
        CALL PVMXUNPACK(NEWGEO,1,'LOGICAL')
        CALL PVMXUNPACK(BASDIR,60,'STRING')
        CALL PVMXUNPACK(LENBAS,1,'INTEGER')
C
C       --------------------------------------------------
C       Unpack common block MOLINP and INFORB (only once).
C       --------------------------------------------------
C
        TIMING = FIRST .OR. NEWBAS .OR. NEWGEO
C
        IF (TIMING) THEN
           IF (IPRINT .GE. 2) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &          ' Receiving new orbital inf.',
     &          ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C  /molinp/
           CALL PVMXUNPACK(NMLINE,1,'INTEGER')
C  /molinc/
           CALL PVMXUNPACK(MLINE,len_MLINE*NMLINE,'STRING')
C  /inforb/
           CALL PVMXUNPACK(MULD2H,NINFI,'INTEGER')
C  /infpar/
           CALL PVMXUNPACK(NODTOT,NPARI,'INTEGER')
C
           DO 100 I = 1, NODTOT
              IF (NODEID(I) .EQ. MYTID) MYNUM = I
 100       CONTINUE
C  /dorps/ DOREPS, DOCOOR (NDORL covers both)
           CALL PVMXUNPACK(DOREPS(0),NDORL,'LOGICAL')
C  dosym()
           CALL PVMXUNPACK(DOSYM,NSYML,'LOGICAL')
C  comr12
           CALL PVMXBCAST(R12INT,NR12L,'LOGICAL')
           CALL PVMXBCAST(INTGAC,NR12I,'INTEGER')
        END IF
C  /parint/
        CALL PVMXUNPACK(JATOM,NINTI,'INTEGER')
        CALL PVMXUNPACK(JNODV,NINTL,'LOGICAL')
#endif
C
C       In case of QM3 calculation, update number of atom types
C       in "line 4" in the molecule file
C       (only # of QM atom types, not # of MM atom types)
        CALL LINE4(MLINE(NMLAU),NONTYP,NSYMOP,CRT,KCHARG,THRS,AUTOSY,
     &             KASYM,ID3,DOCART,DOOWN)
        AUTOSY = .FALSE.
        NOSYM  = NSYMOP.EQ.0
!
        CALL LINE4W(MLINE(NMLAU),NONTYP_QM,NSYMOP,KCHARG,THRS,AUTOSY,
     &              NOSYM,KASYM,ID3,DOCART,DOOWN)

C       set local/common block information from PARINT.
        IATOM  = JATOM
        ITYPE  = JTYPE
        LUDASP = JLUDAS
        LUINTA = JLUINT
        LUONEL = JLUONE
        LUSOL  = JLUSOL
        LUSUPM = JLUSUP
        MAXDIF = JMXDIF
        MAXREP = JMXREP
        NFMAT  = JNFMAT
        NTASK  = JTASK
        I2TYP  = J2TYP
        ICEDIF = JCEDIF
        IFTHRS = JFTHRS
        NOCONT = JNOCNT
        NODV   = JNODV
        NOPV   = JNOPV
        RETUR  = JRETUR
        TKTIME = JTKTIM
        SOLVNT = JSOLVN
        RELCAL = JRELCL
C
C       Set hermit.
C
        IF (NEWGEO .AND. SOLVNT) THEN
           NEWGEO = .TRUE.
           NUCIND = NUCIND + 1
           NUCDEP = NUCDEP + 1
           NATOMS = NATOMS + 1
           NCNTCV = NUCIND
           NCLINE(NUCIND) = 0
           NAMN(NUCIND)       = 'cav '
           NAMEX(3*NUCIND-2)  = 'cav  x'
           NAMEX(3*NUCIND-1)  = 'cav  y'
           NAMEX(3*NUCIND)    = 'cav  z'
           NAMDEP(NUCDEP)     = 'cavity'
           NAMDPX(3*NUCDEP-2) = 'cavity x'
           NAMDPX(3*NUCDEP-1) = 'cavity y'
           NAMDPX(3*NUCDEP  ) = 'cavity z'
           IF (NUCDEP .GT. MXCENT) THEN
              WRITE (LUPRI,'(//2A,/A,I5)')
     &         ' Too many atomic centers: MXCENT exceed in READIN for',
     &         ' solvent cavity,',' Current limit:',MXCENT
              CALL QUIT('*** ERROR *** MXCENT exceeded in READIN')
           END IF
           CORD(1,NUCIND) = D0
           CORD(2,NUCIND) = D0
           CORD(3,NUCIND) = D0
           ISTBNU(NUCIND) = 7
           CHARGE(NUCIND) = D0
           CALL NUCPRO(WORK(KFREE),LFREE)
           NEWGEO = .FALSE.
        END IF
C
        IF (NEWGEO) RDINPC = .FALSE.
C
C
        CALL SETHER(0,NEWGEO,WORK(KFREE),LFREE)
C
C       Set common block BLOCKS
C
        CALL PAOVEC(JSTRSH,NPRIMS,NCONTS,IORBSH,JORBSH,KORBSH,
     &              0,.FALSE.,0)
C
#ifdef VAR_MPI

!       update info about type of fock-matrix construction
        CALL MPIXBCAST(ISYMDM,NFMAT,'INTEGER',MASTER)
        CALL MPIXBCAST(IFCTYP,NFMAT,'INTEGER',MASTER)

!       receive density matrices - that can be huge chunks... # first
        CALL MPIXBCAST(NDMAT,1,'INTEGER',MASTER)

!       set rma_model according to jtype
        if(jtype /= 3) rma_model = .false.

        if(rma_model)then
          ! set # of fmat/dmat outside the RMA window: min = 1
          CALL MPIXBCAST(rma_win_info%nmat_max_wo_win,1,'INTEGER',
     &                   MASTER)
        end if
#endif
        if (jtype .eq. 2 .and. iprint .gt. 5) then
           write(lupri,*) ' HER_RVINIT -1 - type 2, NDMAT =', NDMAT
           write(lupri,*) ' HER_RVINIT -1 - type 2, DOCOOR:'
           write(lupri,'(10(3L2,4X))') ((DOCOOR(i,j),i=1,3), j=1,mxcent)
        end if
        call qexit('HER_RVINIT - 1 -')
        return

      end if ! set_cb_blocks

!     part - 2 - receive data on allocated matrices
!     ---------------------------------------------
      call qenter('HER_RVINIT - 2 -')

!     set local information from PARINT
      itype  = jtype
!     write(lupri,*) 'again ... jtype rma_model',jtype,rma_model
!     call flshfo(lupri)

#ifdef VAR_MPI
      if(rma_model)then

        if(communication_info_mpi%my_shmem_node_id == 0)then
!         distribute density matrix(ces) among NUMA/node masters
          call mpi_bcast(dmat_buff,n2basx*ndmat,mpi_real8,master,
     &                   communication_info_mpi%communication_internode,
     &                   ierr)
        end if

!       distribute density matrix(ces) among node/NUMA node slaves
        call mpi_bcast(dmat_buff,
     &                 n2basx*rma_win_info%nmat_max_wo_win,
     &                 mpi_real8,
     &                 0,
     &                 communication_info_mpi%communication_shmemnode,
     &                 ierr)


      else
!       write(lupri,*) 'slave ndmat, nbast',ndmat,nbast
!       call flshfo(lupri)
        call mpixbcast(dmat_buff,n2basx*ndmat,'DOUBLE',master)
!       ncols = MIN(5,nbast)
!       call output(dmat_buff,1,nbast,1,ncols,nbast,nbast,-1,lupri)
!       call flshfo(lupri)
      end if

      CALL MEMGET('REAL',KHESEE,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
C added by Bin Gao, June 3, 2009
      CALL MPIXBCAST(GENCTR,1,'LOGICAL',MASTER)
C

      IF(ITYPE.EQ.3) THEN
        CALL MPIXBCAST(N2RED,1,'INTEGER',MASTER)
        CALL MEMGET('REAL',KDMRAO,NDMAT*N2RED,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KGABAO,      N2RED,WORK,KFREE,LFREE)
        CALL MPIXBCAST(WORK(KDMRAO),NDMAT*N2RED,'DOUBLE',MASTER)
        CALL MPIXBCAST(WORK(KGABAO),N2RED,'DOUBLE',MASTER)
      ELSE
        KDMRAO = KFREE
        KGABAO = KFREE
      ENDIF
      KDMRSO = KFREE
#endif
#if defined (VAR_PVM)
      call quit('pvm version of Dalton SCF not working - '//
     &          'ask s. knecht at stefan.knecht@phys.chem.ethz.ch '//
     &          'for help.')
C     Unpack density matrices.
C
      CALL PVMXUNPACK(NDMAT,1,'INTEGER')
      CALL MEMGET('REAL',KDMAT,NDMAT*N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFMAT,NFMAT*N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHESEE,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL PVMXUNPACK(WORK(KDMAT),N2BASX*NDMAT,'DOUBLE')
      CALL PVMXUNPACK(ISYMDM,NFMAT,'INTEGER')
      CALL PVMXUNPACK(IFCTYP,NFMAT,'INTEGER')
C
      IF(ITYPE.EQ.3) THEN
        CALL PVMXUNPACK(N2RED,1,'INTEGER')
        CALL MEMGET('REAL',KDMRAO,NDMAT*N2RED,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KGABAO,N2RED,WORK,KFREE,LFREE)
        CALL PVMXUNPACK(WORK(KDMRAO),N2RED*NDMAT,'DOUBLE')
        CALL PVMXUNPACK(WORK(KGABAO),N2RED,'DOUBLE')
        KLAST  = KGABAO + N2RED
      ELSE
        KDMRAO = KFREE
        KGABAO = KFREE
      ENDIF
      KDMRSO = KFREE
#endif
C
      call qexit('HER_RVINIT - 2 -')
      END

C     /* Deck parlop */
Cef begin
      SUBROUTINE PARLOP(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,GMAT,
     &     MAXDER,EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &     SOFOCK,DISTRI,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &     PERTUR,IATOM,MULE,MULTE,NODV,NOPV,NOCONT,
     &     THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &     SQ12EL,INDHER,INDHSQ,IODDHR,TSKCPU,IWHICH,IJS,
     &     ISYMDM,IFCTYP,ADISTR,ITYPE,JSTRSH,NPRIMS,NCONTS,
     &     IORBSH,I2TYP,ICEDIF,IFTHRS,
     &     GABRAO,DMRAO,DMRSO,DINTSKP,RELCALC,GENCTR)
C
C     Copied from TWOLOP and rewritten by Paal Dahle Nov.1994

!     module dependencies
      use parallel_communication_models_mpi
      use pelib_interface, only: use_pelib
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "dummy.h"
#include "mtags.h"
#include "pcmlog.h"
#include "gnrinf.h"
      LOGICAL PRINTA, PRINTB, PRINTC, PRINTD, NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, DIA2SO, ZFS2EL,
     &        DISTRI, NOCONT, SPNORB,
     &        RETUR, FIRST, SQ12EL, LONDON, SUSCEP, ADISTR, DONE,
     &        SOFOCK,RELCALC,GENCTR
      DIMENSION DMAT(*), FMAT(*), GMAT(*),INDHSQ(*),
     & IODDHR(*), INDHER(*), WORK(LWORK), ISYMDM(*), IFCTYP(*),
     &          JSTRSH(*), NPRIMS(*), NCONTS(*), IORBSH(MXSHEL,MXAOVC),
     &          TSKCPU(MTOTTK), IWHICH(MTOTTK), IJS(NTASK),
     &          GABRAO(*), DMRAO(*),DMRSO(*),DINTSKP(2,4), HESSEE(*)
#include "cbisol.h"
#include "twocom.h"
#include "nuclei.h"
#include "energy.h"
#include "taysol.h"
#include "suscpt.h"
#include "blocks.h"
#include "symmet.h"
#include "infpar.h"
#include "inftap.h"
#include "inforb.h"
#include "incore.h"

#include "trkoor.h"
      REAL*8 ERGMOL, GRDMOL(NCOOR), HESMOL(NCOOR,NCOOR) ! automatic arrays
C
C For MPE event logging
C      evIDb2 = MPE_Log_get_event_number()
C      evIDe2 = MPE_Log_get_event_number()
C      call MPE_Describe_state(evIDb2,evIDe2,
C     &     "PARLOP","red")
C
C Begin logging the event
C      call MPE_Log_event(evIDb2,0,'')
C
      CALL QENTER('PARLOP')
      IF (JPRINT .GT. 3) THEN
         CALL TITLER('Output from PARLOP','*',103)
         WRITE(LUPRI,'(/A,I5)') 'Parallel type (ITYPE) is',ITYPE
         WRITE(LUPRI,*) 'MAXDER, EXPECT',MAXDER,EXPECT
      END IF
C
      FIRST = .TRUE.
      DIRAC = RELCALC
      IF (EXPECT .AND. .NOT.NOPV) THEN
         CALL REWSPL(LUPAO)
      END IF
      IF (SUSCEP) THEN
         IF (.NOT.NOPV) CALL REWSPL(LUPAO)
         IF (.NOT.NOPV) CALL REWSPL(LUPAS)
COBL/AMT
C for CAMB3LYP, we must not zero the susceptibilities the second run in TWOLOP
         IF (.NOT. DONETWOLOP)THEN
             DONETWOLOP = .TRUE.
             CALL DZERO(SUS2EL,9)
         END IF
      END IF
C
C
C     For direct contributions: some parameters
C       FCKTHR is threshold for screening
C       ICEFLG gives information about separate screening of
C              Coulomb/exchange for each DMAT
C       NCM    is the number of DMAT requiring Coulomb-contributions
C       NEM    is the number og DMAT requiring exchange-contributions
C       Screening proceeds in three steps as documented by DINTSKP:
C         Step 1: Screening on integral batches
C           DINTSKP(1,1) - total number of integrals
C           DINTSKP(2,1) - number of integrals skipped (batchwise)
C         Step 2: Screening on individual integrals
C                 while unpacking indices
C           DINTSKP(1,2) - number of integrals remaining after step 1
C           DINTSKP(2,2) - number of integrals skipped
C         Step 3a: Screening on Coulomb contributions
C           DINTSKP(1,3) - NCM times number of integrals remaining
C                         after step 2
C           DINTSKP(2,3) - NCM times number of integrals skipped
C         Step 3b: Screening on exchange contributions
C           DINTSKP(1,4) - NEM times number of integrals remaining
C                         after step 2
C           DINTSKP(2,4) - NEM times number of integrals skipped
C
      DOSCRN = .FALSE.
      CALL DZERO(DINTSKP,8)
      IF(DIRFCK.OR.SOFOCK)  THEN
        IF(IFTHRS.LT.16) THEN
          DOSCRN = .TRUE.
          FCKTHR = -IFTHRS
          FCKTHR = 10.0D0**FCKTHR
          ICEFLG = ICEDIF
          NCM = 0
          NEM = 0
          DO I = 1,NDMAT
            IY  = MOD(IFCTYP(I),10)
            IC  = MOD(IY,2)
            NCM = NCM + IC
            IE  = (IY - IC)/2
            NEM = NEM + IE
          ENDDO
        ENDIF
      ENDIF
C
      IF(I2TYP.EQ.0) THEN
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'PARLOP: Unknown I2TYP =' ,I2TYP
        CALL QUIT('Unknown I2TYP !!!')
      ENDIF
C
C     ***********************************
C     ***** Send request for a task *****
C     ***********************************
C
      ICOUNT = 0
      CALL IZERO(IWHICH,MTOTTK)
      CALL DZERO(TSKCPU,MTOTTK)
C
      IPLACE = 0
      ITOTNT = 1
      LINTMP = .FALSE.
      I_SHL = 1
C
C
 100  CONTINUE
C
      IF (AOSAVE) THEN
C
         IF (NTASK .GT. 1) THEN
            CALL QUIT('NTASK > 1 is not ported to AOSAVE = TRUE.')
         END IF
C
#if defined (VAR_MPI)
         CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MTAG5)
         CALL MPIXSEND(IPLACE,1,'INTEGER',MASTER,MTAG51)
C
         IF (IPLACE .NE. 0) THEN
            INTSAV = 0
            IF (LINTSV) THEN
               INTSAV = 1
            END IF
            CALL MPIXSEND(INTSAV,1,'INTEGER',MASTER,MTAG52)
         END IF
         LINTSV = .FALSE.
C
#endif
#if defined (VAR_PVM)
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(MYNUM,1,'INTEGER')
         CALL PVMXSEND(MASTER,MTAG5)
#endif
C
C     Receive IJ labels for NTASK new (IJ|**) integrals
C
         CALL IZERO(IJS,NTASK)
C
#if defined (VAR_MPI)
         CALL MPIXRECV(DONE,1,'LOGICAL',MASTER,MTAG6)
         IF (DONE) GOTO 300
         CALL MPIXRECV(IJS,NTASK,'INTEGER',MASTER,MTAG61)
         IF (IJS(1) .EQ. 0) THEN
            IPLACE = 0
            GOTO 100
         END IF
         CALL MPIXRECV(IINTSV,NTASK,'INTEGER',MASTER,MTAG64)
C         write(lupri,*) 'MYNUM,IJS(1),IINTSV',MYNUM,IJS,IINTSV
C
         IF (IINTSV .EQ. 1) THEN
            LINTSV = .TRUE.
         ELSE
            LINTSV = .FALSE.
         END IF
         CALL MPIXRECV(IPLACE,1,'INTEGER',MASTER,MTAG62)
C
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(MASTER,MTAG6)
         CALL PVMXUNPACK(DONE,1,'LOGICAL')
         IF (DONE) GOTO 300
         CALL PVMXUNPACK(IJS,NTASK,'INTEGER')
#endif
C
C ELSE below is for the compability to old code: AOSAVE = FALSE
C
      ELSE
C
#if defined (VAR_MPI)
         CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MTAG5)
#endif
#if defined (VAR_PVM)
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(MYNUM,1,'INTEGER')
         CALL PVMXSEND(MASTER,MTAG5)
#endif
C
C     Receive IJ labels for NTASK new (IJ|**) integrals
C
         CALL IZERO(IJS,NTASK)
C
#if defined (VAR_MPI)
         CALL MPIXRECV(DONE,1,'LOGICAL',MASTER,MTAG6)
         IF (DONE) GOTO 300
         CALL MPIXRECV(IJS,NTASK,'INTEGER',MASTER,MTAG6)
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(MASTER,MTAG6)
         CALL PVMXUNPACK(DONE,1,'LOGICAL')
         IF (DONE) GOTO 300
         CALL PVMXUNPACK(IJS,NTASK,'INTEGER')
#endif
C
C     END IF below corresponds to: IF (AOSAVE)... ELSE
      END IF
C
      IF (JPRINT .GT. 4) THEN
         WRITE(LUPRI,'(5X,A,I4)')
     &        'Receiving IJ-batches: ',(IJS(KK),KK=1,NTASK)
      END IF
C
      DO 200, I = 1,NTASK
C
         IJSHEL = IJS(I)
C
C        print *,'MYNUM,NTASK',MYNUM,NTASK
C	print *,'IJS(1:NTASK)'
C         write(lupri,*), 'NTASK',NTASK
C   DO 299 LL=1,NTASK
C       DO 299 LL=1,55
C        write(lupri,*),'MYNUM, LL =  ',MYNUM,LL,'IJS(LL) = ',IJS(LL)
C 299     CONTINUE
         IF (IJSHEL .EQ. 0) THEN
            ICOUNT = ICOUNT + I - 1
            GO TO 100
         END IF
C
C        The compressed index IJSHEL is split to indices ISHELA and ISHELB.
C
         CALL UNPKIJ(IJSHEL,ISHELA,ISHELB)
C
C        Check if this distribution is requested, otherwise
C        go for next task [ FOR EFFICIENCY THIS OUGHT NOT TO
C        HAPPEN! ]  950428-K.Ruud+H.J.Aa.Jensen
C
         ICA = LCLASH(ISHELA)
         ICB = LCLASH(ISHELB)
C
         IF (TIMING) CPUTK1 = SECOND()
C
C        *****************************
C        ***** First Shell Index *****
C        *****************************
C
         NHKTA  = NHKTSH(ISHELA)
         KHKTA  = KHKTSH(ISHELA)
         KCKTA  = KCKTSH(ISHELA)
         SPHRA  = SPHRSH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
         MULA   = ISTBSH(ISHELA)
         MULTA  = MULT(MULA)
         NSTRA  = IORBSB(IORBSH(ISHELA,1))
         NUCA   = NUCOSH(ISHELA)
         NORBA  = NORBSH(ISHELA)
         IF (.NOT.BIGVEC) THEN
            CORAX0 = CENTSH(ISHELA,1)
            CORAY0 = CENTSH(ISHELA,2)
            CORAZ0 = CENTSH(ISHELA,3)
         END IF
         PRINTA = .TRUE.
         IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINTA = .FALSE.
C
C        ******************************
C        ***** Second Shell Index *****
C        ******************************
C
         NHKTB  = NHKTSH(ISHELB)
         KHKTB  = KHKTSH(ISHELB)
         KCKTB  = KCKTSH(ISHELB)
         SPHRB  = SPHRSH(ISHELB)
         NCENTB = NCNTSH(ISHELB)
         MULB   = ISTBSH(ISHELB)
         MULTB  = MULT(MULB)
         NSTRB  = IORBSB(IORBSH(ISHELB,1))
         NUCB   = NUCOSH(ISHELB)
         NORBB  = NORBSH(ISHELB)
         IF (.NOT.BIGVEC) THEN
            CORBX0 = CENTSH(ISHELB,1)
            CORBY0 = CENTSH(ISHELB,2)
            CORBZ0 = CENTSH(ISHELB,3)
         END IF
         GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
         IGENAB = 1
         IF (.NOT.GENAB) IGENAB = 2
         NSETA  = NSETSH(ISHELA,IGENAB)
         NSETB  = NSETSH(ISHELB,IGENAB)
         PRINTB = PRINTA
         IF ((ISHELB.NE.IPRNTB).AND.(IPRNTB.NE.0)) PRINTB = .FALSE.
C
C        *****************************
C        ***** Third Shell Index *****
C        *****************************
C
         ICMAX = ISHELA
         IF (SPNORB) ICMAX = MAXSHL
         IF(I2TYP.EQ.2) ICMAX = NLRGBL
         DO 400 ISHELC = ICSTRT, ICMAX
C
            ICC    = LCLASH(ISHELC)
C
            NHKTC  = NHKTSH(ISHELC)
            KHKTC  = KHKTSH(ISHELC)
            KCKTC  = KCKTSH(ISHELC)
            SPHRC  = SPHRSH(ISHELC)
            NCENTC = NCNTSH(ISHELC)
            MULC   = ISTBSH(ISHELC)
            MULTC  = MULT(MULC)
            NSTRC  = IORBSB(IORBSH(ISHELC,1))
            NUCC   = NUCOSH(ISHELC)
            NORBC  = NORBSH(ISHELC)
            IF (.NOT.BIGVEC) THEN
               CORCX0 = CENTSH(ISHELC,1)
               CORCY0 = CENTSH(ISHELC,2)
               CORCZ0 = CENTSH(ISHELC,3)
            END IF
            PRINTC = PRINTB
            IF ((ISHELC.NE.IPRNTC).AND.(IPRNTC.NE.0)) PRINTC=.FALSE.
C
C           ******************************
C           ***** Fourth Shell Index *****
C           ******************************
C
            IDMAX = ISHELC
            IF (.NOT.SPNORB.AND.(ISHELA.EQ.ISHELC)) IDMAX = ISHELB
            DO 500 ISHELD = IDSTRT,IDMAX
C
               ICD    = LCLASH(ISHELD)
C
               NHKTD  = NHKTSH(ISHELD)
               KHKTD  = KHKTSH(ISHELD)
               KCKTD  = KCKTSH(ISHELD)
               SPHRD  = SPHRSH(ISHELD)
               NCENTD = NCNTSH(ISHELD)
               MULD   = ISTBSH(ISHELD)
               MULTD  = MULT(MULD)
               NSTRD  = IORBSB(IORBSH(ISHELD,1))
               NUCD   = NUCOSH(ISHELD)
               NORBD  = NORBSH(ISHELD)
               IF (.NOT.BIGVEC) THEN
                  CORDX0 = CENTSH(ISHELD,1)
                  CORDY0 = CENTSH(ISHELD,2)
                  CORDZ0 = CENTSH(ISHELD,3)
               END IF
               GENCD = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
               IGENCD = 1
               IF (.NOT.GENCD) IGENCD = 2
               NSETC = NSETSH(ISHELC,IGENCD)
               NSETD = NSETSH(ISHELD,IGENCD)
               PRINTD = PRINTC
               IF ((ISHELD .NE. IPRNTD).AND.(IPRNTD .NE. 0))
     &              PRINTD = .FALSE.
C
               SHAEQB = ISHELA .EQ. ISHELB
               SHCEQD = ISHELC .EQ. ISHELD
               SHABAB = (ISHELA.EQ.ISHELC) .AND. (ISHELB.EQ.ISHELD)
C
C              *******************************
C              ***** Calculate integrals *****
C              *******************************
C
c              IPRINT = 0
               IPRINT = JPRINT
C
C     Note DUMMY should be HESSEE, K.Ruud-96
C
               CALL TWOODS(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK,LWORK,
     &                 UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                 EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,
     &                 IATOM,MULE,MULTE,
     &                 MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,
     &                 FIRST,SQ12EL,INDHSQ,IODDHR,INDHER,IFCTYP,
     &                 ADISTR,JSTRSH,NPRIMS,NCONTS,IORBSH,DUMMY,
     &                 ICEDIF,IFTHRS,DINTSKP,
     &                 GABRAO,DMRAO,DMRSO,IREPDM,GENCTR)
C
               IF (RETUR) THEN
                  IF (ISHELA .EQ. IPRNTA .AND.
     &                ISHELB .EQ. IPRNTB .AND.
     &                ISHELC .EQ. IPRNTC .AND.
     &                ISHELD .EQ. IPRNTD) GOTO 9999
               END IF
C
 500        CONTINUE
C
 400  CONTINUE
C
         IF (TIMING) THEN
            CPUTK2 = SECOND()
            CPUTSK = CPUTK2 - CPUTK1
            IF (.NOT.AOSAVE) IPLACE = ICOUNT + I
C            IPLACE = ICOUNT + I
            IWHICH(IPLACE) = IJSHEL
            TSKCPU(IPLACE) = CPUTSK
         END IF
 200  CONTINUE
C
      ICOUNT = ICOUNT + NTASK
      LINTSV = LINTMP
      LINTMP = .FALSE.
      GOTO 100
C
 300  CONTINUE
C     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C     Note that DUMMY here should be HESSEE,K.Ruud, Jan.-97
C
CTROND      IF (DIRFCK .OR. DDFOCK .OR. EXPECT) THEN
CTROND         CALL SKLFCK(FMAT,DUMMY,WORK,LWORK,JPRINT,DIRFCK,DDFOCK,EXPECT,
CTROND     &               PERTUR,NODV,MAXDER,LONDON,NDMAT,ISYMDM,IFCTYP,
CTROND     &               IATOM,.FALSE.)
CTROND      END IF
C
C     ----- Print Section - Gradient and Hessian Elements -----
C
      IF (EXPECT) THEN
         IF (JPRINT .GT. 5) THEN
            KCSTRA = 1
            KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
            KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
            IF (KLAST .GT. LWORK) CALL STOPIT('PARLOP',' ',KLAST,LWORK)
            CALL HEADER('Two-electron integral gradient',-1)
            CALL PRIGRD(GRADEE,WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('Potential energy (NN + NE + EE) gradient',-1)
            CALL ZERGRD
            CALL ADDGRD(GRADNN)
            CALL ADDGRD(GRADNA)
            CALL ADDGRD(GRADEE)
            CALL ABAREAD_TAYMOL(ERGMOL,GRDMOL,HESMOL,NCOOR)
            CALL PRIGRD(GRDMOL,WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('Molecular gradient',-1)
            CALL ADDGRD(GRADFS)
            CALL ADDGRD(GRADKE)
            IF (SOLVNT .OR. PCM) THEN
               CALL ADDGRD(GSOLTT)
               CALL ADDGRD(GSOLNN)
            END IF
            IF (PCM) CALL ADDGRD(GSOLCV)
            IF (USE_PELIB()) CALL ADDGRD(PEGRAD)
            CALL ABAREAD_TAYMOL(ERGMOL,GRDMOL,HESMOL,NCOOR)
            CALL PRIGRD(GRDMOL,WORK(KCSTRA),WORK(KSCTRA))
            NCDEP3 = 3*NUCDEP
            GRDNRM = DDOT(NCDEP3,GRDMOL,1,GRDMOL,1)
            GRDNRM = SQRT(GRDNRM)
            WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *         'Molecular gradient norm:', GRDNRM
            CALL ZERGRD
         END IF
      END IF
      IF (LONDON .AND. MAXDER.EQ.2) THEN
         SUS2EL(2,1) = SUS2EL(1,2)
         SUS2EL(3,1) = SUS2EL(1,3)
         SUS2EL(3,2) = SUS2EL(2,3)
         IF (JPRINT .GT. 5) THEN
            CALL HEADER('Two-electron integral susceptibilities',-1)
            CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
C
C     Print Fock matrices (only node masters)
C
      IF(DIRFCK .AND. JPRINT.GT.5
     &   .and. communication_info_mpi%my_shmem_node_id == 0) THEN
         CALL HEADER('Fock matrix in PARLOP',-1)
         DO 600 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
!           WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
!           CALL OUTPUT(DMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
  600    CONTINUE
      END IF
C
C     Send final results to master
C
      CALL HER_SDRES(FMAT,HESSEE,NDMAT,DINTSKP,TSKCPU,IWHICH,ICOUNT,
     &     ITYPE,JPRINT,GENCTR)
Cef 10.mai begin
C This is a fix for an unresolved error in the AOSAVE-related code.
C The price of the fix is a bit lower memory efficiency.
        MSAVE = .FALSE.
Cef 10.mai end
C
C End logging the event
C      call MPE_Log_event(evIDe2,0,'')
C
 9999 CALL QEXIT('PARLOP')
      RETURN
      END
Cef end
C  /* Deck her_sdres */
      SUBROUTINE HER_SDRES(FMAT,HESSEE,NDMAT,DINTSKP,TSKCPU,IWHICH,
     &                     ICOUNT,ITYPE,IPRINT,GENCTR)
!     module dependencies
      use rma_windows

#ifdef VAR_MPI
      use parallel_communication_models_mpi
      use one_sided_communication_wrappers
#ifdef USE_MPI_MOD_F90
      use mpi
#endif
#endif
!
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "mtags.h"
C added by Bin Gao, June 3, 2009
      logical GENCTR
      DIMENSION FMAT(*), TSKCPU(ICOUNT), IWHICH(ICOUNT),
     &          DINTSKP(8), HESSEE(*)
C
C Used from common blocks:
C
C MXCENT: MXCOOR (for ENERGY)
C MXORB : MXSHEL (for INFPAR)
C ENERGY: GRADEE,HESSEE
C SUSCPT: SUS2EL
C INFORB: N2BASX
C INFPAR: MTOTTK, MASTER, TIMING
C
#include "energy.h"
#include "suscpt.h"
#include "inforb.h"
#include "infpar.h"
#include "expopt.h"
#ifdef VAR_MPI
#ifndef USE_MPI_MOD_F90
#include "mpif.h"
#endif
       integer (kind=mpi_offset_kind) :: offset_mat
#endif
C
      IF (IPRINT .GT. 3) THEN
         CALL TITLER('Output from HER_SDRES','*',103)
         WRITE(LUPRI,'(/A,I5)') 'Parallel type (ITYPE) is',ITYPE
      END IF
C
      L_GRA = MXCOOR
      L_HES = MXCOOR*MXCOOR
      L_SUS = 9
C
#if defined (VAR_MPI)
      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MTAG7)
C
      IF (ITYPE .EQ. 2) THEN
         CALL MPIXSEND(GRADEE,L_GRA,'DOUBLE',MASTER,MTAG7)
         CALL MPIXSEND(HESSEE,L_HES,'DOUBLE',MASTER,MTAG7)
         CALL MPIXSEND(ALPGRD,MXPRIM,'DOUBLE',MASTER,MTAG7)
      ELSE IF (ITYPE .EQ. 3) THEN
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
         CALL BIN_SUM(NDMAT*N2BASX, FMAT)
         CALL BIN_SUM(8, DINTSKP)
#else
         if(.not.rma_model)then
           CALL MPIXSEND(FMAT,NDMAT*N2BASX,'DOUBLE',MASTER,MTAG7)
           CALL MPIXSEND(DINTSKP,8,'DOUBLE',MASTER,MTAG7)
         end if
#endif
      ELSE IF (ITYPE .EQ. -5 .OR. ITYPE .EQ. 8) THEN
         CALL MPIXSEND(FMAT,NFMAT*N2BASX,'DOUBLE',MASTER,MTAG7)
         IF (ITYPE .EQ. -5) THEN
            CALL MPIXSEND(SUS2EL,L_SUS ,'DOUBLE',MASTER,MTAG7)
         END IF
      END IF
C
      CALL MPIXSEND(ICOUNT,1,'INTEGER',MASTER,MTAG7)
C
      IF (TIMING) THEN
         CALL MPIXSEND(IWHICH,ICOUNT,'INTEGER',MASTER,MTAG7)
         CALL MPIXSEND(TSKCPU,ICOUNT,'DOUBLE' ,MASTER,MTAG7)
      END IF

!     collect data for rma_model
      if(rma_model)then
        call mpi_reduce(dintskp,
     &                  mpi_in_place,
     &                  8,
     &                  mpi_real8,
     &                  mpi_sum,
     &                  master,
     &                  communication_info_mpi%
     &                  communication_glb_world,
     &                  ierr)
!       all fmatS outside the RMA window...
        call mpi_reduce(fmat,
     &                  mpi_in_place,
     &                  rma_win_info%nmat_max_wo_win * n2basx,
     &                  mpi_real8,
     &                  mpi_sum,
     &                  master,
     &                  communication_info_mpi%
     &                  communication_glb_world,
     &                  ierr)

!       now all fmatS inside the RMA window... if any
        if(rma_win_info%nmat_max_wo_win < ndmat)then

          call mpi_barrier(communication_info_mpi%
     &                     communication_shmemnode,
     &                     ierr)
          if(communication_info_mpi%my_shmem_node_id == 0)then

            call mpi_reduce(fmat(1 + rma_win_info%
     &                      nmat_max_wo_win * n2basx),
     &                      mpi_in_place,
     &                      (ndmat - rma_win_info%
     &                      nmat_max_wo_win) * n2basx,
     &                      mpi_real8,
     &                      mpi_sum,
     &                      master,
     &                      communication_info_mpi%
     &                      communication_internode,
     &                      ierr)
          end if
        end if
      end if

#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      IF (ITYPE .EQ. 2) THEN
         CALL PVMXPACK(GRADEE,L_GRA,'DOUBLE')
         CALL PVMXPACK(HESSEE,L_HES,'DOUBLE')
      ELSE IF (ITYPE .EQ. 3) THEN
         CALL PVMXPACK(FMAT,NDMAT*N2BASX,'DOUBLE')
         CALL PVMXPACK(DINTSKP,8,'DOUBLE')
      ELSE IF (ITYPE .EQ. -5 .OR. ITYPE .EQ. 8) THEN
         CALL PVMXPACK(FMAT,NFMAT*N2BASX,'DOUBLE')
         IF (ITYPE .EQ. -5) THEN
            CALL PVMXPACK(SUS2EL,L_SUS,'DOUBLE')
         END IF
      END IF
C
C     Pack the number of tasks received.
C
      CALL PVMXPACK(ICOUNT,1,'INTEGER')
C
C     Pack the integral-calculation times.
C
      IF (TIMING) THEN
         CALL PVMXPACK(IWHICH,ICOUNT,'INTEGER')
         CALL PVMXPACK(TSKCPU,ICOUNT,'DOUBLE')
      END IF
C
      CALL PVMXSEND(MASTER,MTAG7)
#endif
C
      RETURN
      END
C  /* Deck recompi */
      SUBROUTINE RECOMPI
C
#include "implicit.h"
#include "priunit.h"
C
      WRITE(LUPRI,'(//2(/5X,A))')
     &  'FOR RECOMPILATION:',
     &  'Change parameter in infpar.h and rebuild with make'
C
      RETURN
      END
C  /* Deck prlinp */
      SUBROUTINE PRLINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "mtags.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      PARAMETER (NTABLE = 6)
      LOGICAL FIRST, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "infpar.h"
      SAVE FIRST
      DATA TABLE /'.NODES ','.DEGREE','.PRINT ','.ENCODE','.DEBUG ',
     &            '.RMA-MD'/
      DATA FIRST /.TRUE./
C
C-------------------------
C     Initialize /PRLINP/.
C-------------------------
C
      IF (.NOT. FIRST) THEN
         IF (WORD .NE. '*END OF' .OR. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      FIRST = .FALSE.
C
C     For infpar.h
C
      SLAVE  = .FALSE.
C     ... this is the master in a parallel run
      INFPAR_DEBUG  = .FALSE.
      NTASK  = 1
      IPRPAR = 0
      NDEGDI = 5
#if defined (VAR_PVM)
      NODTOT = 0
      NCODE  = PVMDATADEFAULT
#endif
      rma_model  = .false.
C
      NEWDEF = (WORD .EQ. '*PARALL')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in PRLINP')
    1          CONTINUE
                  READ(LUCMD,*) NODS
#if defined (VAR_MPI)
                  IF (NODS .NE. NODTOT) THEN
                     WRITE (LUPRI,'(/A,I3,A,I3)')
     &                 'Number of nodes in input ', NODS,
     &                 ' .ne. nodes in MPI calculation ',NODTOT
                     CALL QUIT('Wrong number of nodes in input')
                  END IF
#endif
#if defined (VAR_PVM)
                  NODTOT = NODS
#endif
               GO TO 100
    2          CONTINUE
                  READ(LUCMD,*) NDEGDI
               GO TO 100
    3          CONTINUE
                  READ(LUCMD,*) IPRPAR
               GO TO 100
    4          CONTINUE
#if defined (VAR_PVM)
                  READ (LUCMD,'(A7)') WORD
                  CALL UPCASE(WORD)
                  IF (WORD.EQ.' DATADE') THEN
                     NCODE = PVMDATADEFAULT
                  ELSE IF (WORD.EQ.' DATARA') THEN
                     NCODE = PVMDATARAW
                  ELSE IF (WORD.EQ.' DATAIN') THEN
                     NCODE = PVMDATAINPLACE
                  ELSE
                     WRITE(LUPRI,'(/1X,A,A,A//1X,A/15X,A/15X,A/)')
     &              '.ENCODE type ',word,' unknown!!!',
     &              'Options are : DATADEFAULT','DATARAW','DATAINPLACE'
                     CALL QUIT('Nonexistent encoding type in input.')
                  END IF
#else
                  WRITE(LUPRI,'(/1X,A)')'.ENCODE keyword is not in use.'
#endif
               GO TO 100
    5          CONTINUE   ! .DEBUG
                  INFPAR_DEBUG = .TRUE.
               GO TO 100
    6          CONTINUE   ! .RMA-MD
                  rma_model = .true.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in PRLINP')
            END IF
      END IF
  300 CONTINUE
C
      IF (NODTOT .GT. MAXNOD) THEN
         WRITE(LUPRI,'(/5X,A,/12X,A,2(/5X,A20,I8))')
     &    'ERROR: Number of nodes allocated has exceeded hardcoded',
     &    'limit. Reduce number of nodes or recompile program. ',
     &    'Allocated (NODTOT) :',NODTOT,'Limit (MAXNOD)     :',MAXNOD
         CALL RECOMPI
         CALL QUIT('ERROR: Number of nodes requested is too large.')
      ENDIF
C
      RETURN
      END
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
      SUBROUTINE BIN_SUM(N, SUM)
C     ******************************************************************
C     *
C     *  Performs a binary tree sum of SUM(0:N-1) across all PE's and
C     *  leaves the result in SUM at node 0. This is equivalent to
C     *
C     *     MPI_Reduce(SUM, SUM, N, MPI_DOUBLE_PRECISION, MPI_SUM, 0,
C     *                MPI_COMM_WORLD, IERR)
C     *
C     * Input:
C     *  N      - dimension of input array SUM.
C     *  SUM    - REAL*8 array to be summed across nodes.
C     *
C     * Output:
C     *  SUM   - sum of SUM() across all nodes is returned on PE 0,
C     *          the arrays on other PEs are not modified.
C     *
C     * NOTES:
C     *        - stream safe SHMEM implementation.
C     *        - use only 1072 words of buffer space.
C     *        - calls BARRIER() which implies this is a globally
C     *          synchronizing call.
C     *        - tested with CRAY T3E only. At least some udcflushes or
C     *          some other magics will be necessary on the T3D.
C     *
C     * -BEGIN-TERMS-
C     *  Redistribution is disallowed without notifying the author. In
C     *  no event shall the lines between and including -BEGIN-TERMS-
C     *  and --END-TERMS-- be modified.
C     *
C     *  (C) 1997, Jorn Amundsen, University Administration, NTNU.
C     * --END-TERMS--
C     *
C     ******************************************************************

      INTEGER N
      REAL*8 SUM(0:N)

C---  Cray SHMEM routines
      INTEGER SHMEM_MY_PE, SHMEM_N_PES, SHMEM_INT8_CSWAP

      INTEGER NS, NR, NPROC, ME
      LOGICAL DOSUM

C---  Communication control data (NPAD is for streams padding)
      INTEGER NSEG, NPAD, M_RD, M_WRT
      PARAMETER (NSEG = 512, NPAD = 24, M_RD = 0, M_WRT = -1)
      INTEGER ISWP, ISEG, I, J
CDIR $SUPPRESS STAT
      INTEGER STAT(2)
      REAL TMP(0:2*(NSEG+NPAD))
CDIR$ CACHE_ALIGN /BINSUM__INTERNAL/
      COMMON /BINSUM__INTERNAL/ TMP,STAT



C---  Initialization
      ME = SHMEM_MY_PE()
      NPROC = SHMEM_N_PES()
      IF (NPROC .LE. 1) RETURN
      DOSUM = .TRUE.
      NS = 1


      DO WHILE (NS .LT. NPROC)

         NR = NS * 2
         CALL BARRIER() ! if DOSUM is false do barriers only

         IF (DOSUM) THEN
            IF (MOD(ME,NR) .EQ. 0 .AND. ME+NS .LT. NPROC) THEN
C---           Receiver: add up SUM and TMP in size NSEG chunks
               ISWP = 0
               DO ISEG = 0, N-1, NSEG
C                 Busy wait on STAT(ISWP) == M_WRT
 100              CONTINUE
                  IF (SHMEM_INT8_CSWAP(STAT(ISWP), M_WRT, M_WRT, ME)
     $                 .NE. M_WRT) GOTO 100

                  J = ISWP*(NSEG+NPAD)
CDIR$             UNROLL 8
                  DO I = ISEG, MIN(ISEG+NSEG-1,N-1)
                     SUM(I) = SUM(I) + TMP(J)
                     J = J + 1
                  ENDDO

C                 Update the status word
                  CALL SHMEM_INT8_SWAP(STAT(ISWP), M_RD, ME)
                  ISWP = XOR(ISWP,1)
               ENDDO
            ELSE IF (MOD(ME,NR) .NE. 0) THEN
C---           Sender: write SUM in segments of size NSEG to node ME-NS
               ISWP = 0
               DO ISEG = 0, N-1, NSEG
C                 Busy wait on STAT(ISWP) == M_RD
 200              CONTINUE
                  IF (SHMEM_INT8_CSWAP(STAT(ISWP), M_RD, M_RD, ME-NS)
     $                 .NE. M_RD) GOTO 200

                  J = ISWP*(NSEG+NPAD)
                  CALL SHMEM_PUT(TMP(J), SUM(ISEG), MIN(NSEG,N-ISEG),
     $                 ME-NS)

C                 Update the status word
                  CALL SHMEM_INT8_SWAP(STAT(ISWP), M_WRT, ME-NS)
                  ISWP = XOR(ISWP,1)
               ENDDO
               DOSUM = .FALSE.
            ENDIF
         ENDIF

         NS = NS * 2
      ENDDO

      CALL BARRIER()
      END
#endif
#if defined (VAR_MPI)
C /* Deck J1INTP */
      SUBROUTINE J1INTP(NBAST,NOSIM,KSYMP,DMAT,EXP1VL,EXPVAL,TOFILE,
     &                  IPRTYP,MATDIM,WORK,LWORK)
C
C     Master routine for distributing tesseraes to slaves
C     K.Ruud, June 10 2005, Pisa
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "pcmdef.h"
C
      LOGICAL EXP1VL, TOFILE
      DIMENSION DMAT(MATDIM,NOSIM), WORK(LWORK)
#include "priunit.h"
#include "pcm.h"
      DIMENSION EXPVAL(NTS,NOSIM)
#include "mtags.h"
#include "infpar.h"
#include "mpif.h"
C
      IF (TOFILE) CALL QUIT('Parallel calculations do not allow for '//
     &     'storing PCM-integrals on disk')
C
C     Wake up the slaves
C
      NNBASX = NBAST*(NBAST + 1)/2
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRPCM,1,'INTEGER',MASTER)
C
      CALL MPIXBCAST(KSYMP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(EXP1VL,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NTSIRR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      IF (EXP1VL) THEN
         CALL MPIXBCAST(DMAT,NNBASX,'DOUBLE',MASTER)
      ELSE
         CALL MPIXBCAST(EXPVAL,NTS*NOSIM,'DOUBLE',MASTER)
         CALL DZERO(DMAT,MATDIM*NOSIM)
      END IF
      CALL MPIXBCAST(XTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(YTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
C
C     Loop over all tesserae
C
      DO ITS = 1, NTSIRR
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(ITS,1,'INTEGER',NWHO,MPTAG2)
      END DO
C
C     Send end message to all slaves
C
      ITS = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(ITS ,1,'INTEGER',NWHO,MPTAG2)
      END DO
C
C     Collect data from all slaves
C
      IF (EXP1VL) THEN
         CALL DZERO(WORK,NTS)
         CALL MPI_REDUCE(WORK,EXPVAL,NTS,MPI_DOUBLE_PRECISION,
     &                   MPI_SUM,0,MPI_COMM_WORLD,IERR)
      ELSE IF (.NOT. TOFILE) THEN
C
C        WORK(1) = DMAT(1)
C
         CALL DZERO(WORK(NOSIM*MATDIM+1),MATDIM*NOSIM)
         CALL MPI_REDUCE(WORK(NOSIM*MATDIM+1),DMAT,MATDIM*NOSIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      END IF
      RETURN
      END
C
C /* Deck JINTS */
      SUBROUTINE J1INTS(WORK,LWORK,IPRTYP,IPRTMP)
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "pcmdef.h"
#include "iratdef.h"
      LOGICAL TOFILE, EXP1VL, TRIMAT
      CHARACTER*7 INTLAB
      CHARACTER*8 LABINT(3*MXCOOR)
      DIMENSION WORK(LWORK)
C
#include "priunit.h"
#include "orgcom.h"
#include "symmet.h"
#include "infpar.h"
#include "pcm.h"
#include "mtags.h"
#include "mpif.h"
C
      IPRPCM = IPRTMP
      CALL MPIXBCAST(KSYMP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(EXP1VL,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NTSIRR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      NNBASX = NBAST*(NBAST + 1)/2
      N2BASX = NBAST*NBAST
      IF (IPRTYP .EQ. 11) THEN
         MATDIM = NNBASX
         TRIMAT = .TRUE.
         NCOMP  = 1
         INTLAB = 'NPETES '
      ELSE
         MATDIM = N2BASX
         TRIMAT = .FALSE.
         NCOMP = 3
         INTLAB = 'PCMBSOL'
      END IF
C
      KDEN   = 1
      KINTRP = KDEN + MATDIM*NOSIM
      KINTAD = KINTRP + (3*MXCOOR + 1)/IRAT
      KEXPVL = KINTAD + (3*MXCOOR + 1)/IRAT
      KLST1  = KEXPVL + NTSIRR*(MAXREP + 1)*NOSIM
C
      IF (EXP1VL) THEN
         CALL MPIXBCAST(WORK(KDEN),NNBASX,'DOUBLE',MASTER)
         CALL DZERO(WORK(KEXPVL),NTSIRR*(MAXREP + 1))
      ELSE
         CALL MPIXBCAST(WORK(KEXPVL),NTSIRR*(MAXREP + 1)*NOSIM,'DOUBLE',
     &                  MASTER)
         CALL DZERO(WORK(KDEN),MATDIM*NOSIM)
      END IF
      CALL MPIXBCAST(XTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(YTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
C
C     Loop over tesseraes to be calculated on this node
C
 10   CONTINUE
      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(ITS,1,'INTEGER',MASTER,MPTAG2)
C
      IF (ITS .GT. 0) THEN
         DIPORG(1) = XTSCOR(ITS)
         DIPORG(2) = YTSCOR(ITS)
         DIPORG(3) = ZTSCOR(ITS)
         NTESP = 1
         KPATOM = 0
C
C        Calculates nuclear potential energy integrals (in AO basis) for
C        the given tessera
C
         L=1
         KTMP = KLST1
         IF (.NOT. TOFILE .AND. .NOT. EXP1VL) THEN
            KMAT = KTMP + 8
            IF (IPRTYP .EQ. 11) THEN
               KLAST = KMAT + (MAXREP + 1)*MATDIM
            ELSE
               KLAST = KMAT + (MAXREP + 1)*MATDIM*NCOMP
            END IF
            NCOMP = (MAXREP + 1)
         ELSE
            KMAT = KTMP + 8
            KLAST = KMAT
            NCOMP = 0
         END IF
         LWRK = LWORK - KLAST + 1
         CALL GET1IN(WORK(KMAT),INTLAB,NCOMP,WORK(KLAST),LWRK,LABINT,
     &               WORK(KINTRP),WORK(KINTAD),L,TOFILE,KPATOM,TRIMAT,
     &               WORK(KTMP),EXP1VL,WORK(KDEN),IPRPCM)
         IF (IPRTYP .EQ. 13) THEN
            DO IOSIM = 1, NOSIM
               IADR = KEXPVL + ITS - 1 + NTS*(IOSIM - 1)
               JADR = KMAT + (IOSIM - 1)*MATDIM
               MADR = KDEN + (IOSIM - 1)*MATDIM
               CALL DAXPY(MATDIM,WORK(IADR),WORK(JADR),1,
     &                 WORK(MADR),1)
            END DO
         ELSE IF (EXP1VL) THEN
            DO I = 1, NCOMP
               IADR = KEXPVL + ITS - 1 + (I-1)*NTSIRR
               WORK(IADR) = -WORK(KTMP+I-1)
            END DO
         ELSE
            DO IOSIM = 1, NOSIM
               IADR = KEXPVL + ITS - 1 + NTSIRR*(MAXREP + 1)*(IOSIM - 1)
               KADR = KMAT + (KSYMP - 1)*MATDIM
               MADR = KDEN + (IOSIM - 1)*MATDIM
               CALL DAXPY(MATDIM,-WORK(IADR),WORK(KADR),1,
     &                 WORK(MADR),1)
            END DO
         END IF
         GO TO 10
      END IF
C
C     No more tesseraes to calculate
C
      IF (EXP1VL) THEN
         CALL MPI_REDUCE(WORK(KEXPVL),MPI_IN_PLACE,NTSIRR*(MAXREP + 1),
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      ELSE IF (.NOT. TOFILE) THEN
         CALL MPI_REDUCE(WORK(KDEN),MPI_IN_PLACE,MATDIM*NOSIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      END IF
      RETURN
      END
C /* Deck J1XP */
      SUBROUTINE J1XP(NOSIM,VTEX,UCMO,UBO,UDV,UDVTR,XJ1SQ,XJ1SQT,
     &                TOFILE,TRPLET,JWOPSY,WORK,LWORK)
C
C     Master routine for distributing tesseraes to slaves for
C     calculating one-index transformed potentials
C     K.Ruud, June 10 2005, Pisa
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxash.h"
#include "pcmdef.h"
C
      LOGICAL TOFILE, TRPLET
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "pcm.h"
#include "inforb.h"
#include "infind.h"
      DIMENSION VTEX(NTS,NOSIM), UCMO(NORBT*NBAST), UBO(NOSIM*N2ORBX),
     &          UDV(N2ASHX), UDVTR(N2ASHX), XJ1SQ(NOSIM*N2ORBX),
     &          XJ1SQT(NOSIM*N2ORBX)
#include "mtags.h"
#include "infpar.h"
#include "wrkrsp.h"
#include "mpif.h"
C defined parallel calculation types
#include "iprtyp.h"
C
      IF (TOFILE) CALL QUIT('Parallel calculations do not allow for '//
     &     'storing PCM-integrals on disk')
C
C     Wake up the slaves
C
      IPRTYP = J1XP_WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRPCM,1,'INTEGER',MASTER)
C
      CALL MPIXBCAST(NTSIRR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
C
      CALL MPIXBCAST(QSE,NTS,'DOUBLE',MASTER)
      CALL MPIXBCAST(QSN,NTS,'DOUBLE',MASTER)
      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(UBO,N2ORBX*NOSIM,'DOUBLE',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(XTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(YTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZTSCOR,NTSIRR,'DOUBLE',MASTER)
      IF (TRPLET) THEN
         CALL MPIXBCAST(UDVTR,N2ASHX,'DOUBLE',MASTER)
      ELSE
         CALL MPIXBCAST(UDV,N2ASHX,'DOUBLE',MASTER)
      END IF
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
C
C     Loop over all tesserae
C
      DO ITS = 1, NTSIRR
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(ITS,1,'INTEGER',NWHO,MPTAG2)
      END DO
C
C     Send end message to all slaves
C
      ITS = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(ITS ,1,'INTEGER',NWHO,MPTAG2)
      END DO
C
C     Collect data from all slaves
C
      NDIMZ = MAX(N2ORBX,NTS)
      CALL DZERO(WORK,NOSIM*NDIMZ)
      CALL DZERO(VTEX,NTS*NOSIM)
      CALL MPI_REDUCE(WORK,VTEX,NOSIM*NTS,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
      IF (TRPLET) THEN
         CALL MPI_REDUCE(WORK,XJ1SQT,NOSIM*N2ORBX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      ELSE
         CALL MPI_REDUCE(WORK,XJ1SQ,NOSIM*N2ORBX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      END IF
      RETURN
      END
C /* Deck J1XS */
      SUBROUTINE J1XS(WORK,LWORK,IPRTMP)
C
C     Parallell construction of one-index transformed potentials
C     K.Ruud, Pisa June 2005
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "maxash.h"
#include "pcmdef.h"
#include "iratdef.h"
      PARAMETER (D0 = 0.0D0)
      LOGICAL TOFILE, EXP1VL, TRIMAT, TRPLET
      CHARACTER*8 LABINT(3*MXCOOR)
      DIMENSION WORK(LWORK)
C
#include "priunit.h"
#include "orgcom.h"
#include "symmet.h"
#include "infpar.h"
#include "inforb.h"
#include "infind.h"
#include "pcm.h"
#include "mtags.h"
#include "mpif.h"
C
      IPRPCM = IPRTMP
      CALL MPIXBCAST(NTSIRR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
C
      KINTRP = 1
      KINTAD = KINTRP + (3*MXCOOR + 1)/IRAT
      KEXPVL = KINTAD + (3*MXCOOR + 1)/IRAT
      KJ1AO  = KEXPVL + NTSIRR*(MAXREP + 1)*NOSIM
      KJ1    = KJ1AO  + NNBASX*(MAXREP + 1)
      KJ1SQ  = KJ1    + NNORBX
      KJ1XSQ = KJ1SQ  + N2ORBX
      KJ1XAC = KJ1XSQ + N2ORBX*NOSIM
      KUCMO  = KJ1XAC + N2ASHX*NOSIM
      KUDV   = KUCMO  + NORBT*NBAST
      KUBO   = KUDV   + N2ASHX
      KPCMX  = KUBO   + N2ORBX*NOSIM
      KLAST  = KPCMX  + N2ORBX*NOSIM
      LWRK   = LWORK - KLAST + 1
C
      CALL MPIXBCAST(QSE,NTSIRR*(MAXREP + 1),'DOUBLE',MASTER)
      CALL MPIXBCAST(QSN,NTSIRR*(MAXREP + 1),'DOUBLE',MASTER)
      CALL MPIXBCAST(WORK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WORK(KUBO),N2ORBX*NOSIM,'DOUBLE',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(XTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(YTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZTSCOR,NTSIRR,'DOUBLE',MASTER)
      CALL MPIXBCAST(WORK(KUDV),N2ASHX,'DOUBLE',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
C
C     Loop over tesseraes to be calculated on this node
C
      CALL DZERO(WORK(KPCMX),NOSIM*N2ORBX)
      CALL DZERO(WORK(KEXPVL),NOSIM*NTSIRR*(MAXREP + 1))
 10   CONTINUE
      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(ITS,1,'INTEGER',MASTER,MPTAG2)
C
      IF (ITS .GT. 0) THEN
         L = 1
         NCOMP     = MAXREP + 1
         DIPORG(1) = XTSCOR(ITS)
         DIPORG(2) = YTSCOR(ITS)
         DIPORG(3) = ZTSCOR(ITS)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WORK(KJ1AO),'NPETES ',NCOMP,WORK(KLAST),LWRK,
     &               LABINT,WORK(KINTRP),WORK(KINTAD),L,TOFILE,KPATOM,
     &               TRIMAT,DUMMY,EXP1VL,DUMMY,IPRPCM)
         JJ1AO = KJ1AO
C     Transform AO pot. int. into MO basis  V(AO) --> V(MO)
         CALL UTHU(WORK(JJ1AO),WORK(KJ1),WORK(KUCMO),WORK(KLAST),
     *             NBAST,NORBT)
C Transform V(MO) from triangular to square format
         CALL DSPTSI(NORBT,WORK(KJ1),WORK(KJ1SQ))
         CALL DZERO(WORK(KJ1XSQ),N2ORBX*NOSIM)
         DO IOSIM = 1, NOSIM
            JUBO = KUBO + (IOSIM - 1) * N2ORBX
            JJ1X = KJ1XSQ + (IOSIM - 1) * N2ORBX
            JJ1XAC = KJ1XAC + (IOSIM - 1) * N2ASHX
            JPCMX  = KPCMX + (IOSIM - 1)*N2ORBX
            CALL ONEXH1(WORK(JUBO),WORK(KJ1SQ),WORK(JJ1X))
C
            IF (NASHT .GT. 0) CALL GETACQ(WORK(JJ1X),WORK(JJ1XAC))
            IF (IPRPCM .GE. 15) THEN
               WRITE (LUPRI,'(/A)') ' J1X_mo matrix:'
               CALL OUTPUT(WORK(JJ1X),1,NORBT,1,NORBT,
     &              NORBT,NORBT,1,LUPRI)
               IF (NASHT .GT. 0) THEN
                  WRITE (LUPRI,'(/A)') ' J1X_ac matrix:'
                  CALL OUTPUT(WORK(JJ1XAC),1,NASHT,1,NASHT,
     &                 NASHT,NASHT,1,LUPRI)
               END IF
            END IF
C
C     Expectation value of transformed potential on tesserae:
C                     <0|\tilde{V}|0>
C
            IF (KSYMOP .EQ. 1) THEN
               IF (TRPLET) THEN
                  FACTOR = SLVTLM(WORK(KUDV),WORK(JJ1XAC),WORK(JJ1X),
     &                            TJ1XAC)
               ELSE
                  FACTOR = SLVQLM(WORK(KUDV),WORK(JJ1XAC),WORK(JJ1X),
     &                            TJ1XAC)
               END IF
               IADR = KEXPVL + NTSIRR*(MAXREP + 1)*(IOSIM - 1)+ITS-1
               WORK(IADR) = FACTOR
               IF (IPRPCM .GE. 6) THEN
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- J1X expectation value :',WORK(IADR)
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- active part of J1X    :',TJ1XAC
               END IF
            END IF
C     KPCMX: \tilde{J} + \tilde{X}
            CALL DAXPY(N2ORBX,QSN(ITS)+QSE(ITS),WORK(JJ1X),
     &                 1,WORK(JPCMX),1)
         END DO
C
C     Non-totally symmetric perturbation operators
C
         IF (KSYMOP .GT. 1) THEN
            JJ1AO = KJ1AO + (KSYMOP - 1)*NNBASX
            JTS = (KSYMOP - 1)*NTSIRR + ITS
C     Transform AO pot. int. into MO basis  V(AO) --> V(MO)
            CALL UTHU(WORK(JJ1AO),WORK(KJ1),WORK(KUCMO),WORK(KLAST),
     *                NBAST,NORBT)
C Transform V(MO) from triangular to square format
            CALL DSPTSI(NORBT,WORK(KJ1),WORK(KJ1SQ))
            CALL DZERO(WORK(KJ1XSQ),N2ORBX*NOSIM)
            DO IOSIM = 1, NOSIM
               JUBO = KUBO + (IOSIM - 1) * N2ORBX
               JJ1X = KJ1XSQ + (IOSIM - 1) * N2ORBX
               JJ1XAC = KJ1XAC + (IOSIM - 1) * N2ASHX
               CALL ONEXH1(WORK(JUBO),WORK(KJ1SQ),WORK(JJ1X))
C
               IF (NASHT .GT. 0) CALL GETACQ(WORK(JJ1X),WORK(JJ1XAC))
               IF (IPRPCM .GE. 15) THEN
                  WRITE (LUPRI,'(/A)') ' J1X_mo matrix (KSYMOP):'
                  CALL OUTPUT(WORK(JJ1X),1,NORBT,1,NORBT,
     &                 NORBT,NORBT,1,LUPRI)
                  IF (NASHT .GT. 0) THEN
                     WRITE (LUPRI,'(/A)') ' J1X_ac matrix (KSYMOP):'
                     CALL OUTPUT(WORK(JJ1XAC),1,NASHT,1,NASHT,
     &                    NASHT,NASHT,1,LUPRI)
                  END IF
               END IF
C
C     Expectation value of transformed potential on tesserae:
C                     <0|\tilde{V}|0>
C
               IF (TRPLET) THEN
                  FACTOR = SLVTLM(WORK(KUDV),WORK(JJ1XAC),WORK(JJ1X),
     &                            TJ1XAC)
               ELSE
                  FACTOR = SLVQLM(WORK(KUDV),WORK(JJ1XAC),WORK(JJ1X),
     &                            TJ1XAC)
               ENDIF
               IADR = KEXPVL + NTSIRR*(MAXREP + 1)*(IOSIM - 1)+JTS-1
               WORK(IADR) = FACTOR
               IF (IPRPCM .GE. 6) THEN
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- J1X expectation value :',WORK(IADR)
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- active part of J1X    :',TJ1XAC
               END IF
            ENDDO
         END IF
         GO TO 10
      END IF
C
C     Return expectation values and one-index transformed potentials
C
      CALL MPI_REDUCE(WORK(KEXPVL),MPI_IN_PLACE,
     &                NOSIM*NTSIRR*(MAXREP + 1),MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(WORK(KPCMX),MPI_IN_PLACE,NOSIM*N2ORBX,
     &                MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                IERR)
      RETURN
      END
C /* Deck JX2P */
      SUBROUTINE J2XP(IKLVL,ISYMV1,ISYMV2,ISYMV3,ISYMDN,ZYM1,ZYM2,ZYM3,
     &                TCHG,OVLAP,CMO,DEN1,WORK,LWORK)
C
C     Master routine for calculating multiple one-index transformations
C     of PCM potentials, K.Ruud, Tromso June 2005
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "pcmdef.h"
C
#include "inforb.h"
#include "infdim.h"
      DIMENSION ZYM1(N2ORBX), ZYM2(N2ORBX), ZYM3(N2ORBX), TCHG(NTS),
     &          CMO(NORBT*NBAST), DEN1(NASHDI,NASHDI)
C
#include "priunit.h"
#include "infpar.h"
#include "mtags.h"
#include "pcm.h"
#include "mpif.h"
C defined parallel calculation types
#include "iprtyp.h"
C
C     Wake up the slaves
C
      IPRTYP = J2XP_WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRPCM,1,'INTEGER',MASTER)
C
C     Send necessary start-up data
C
      CALL MPIXBCAST(IKLVL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(CMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(DEN1,NASHDI*NASHDI,'DOUBLE',MASTER)
      IF (IKLVL .GE. 1) CALL MPIXBCAST(ZYM1,N2ORBX,'DOUBLE',MASTER)
      IF (IKLVL .GE. 2) CALL MPIXBCAST(ZYM2,N2ORBX,'DOUBLE',MASTER)
      IF (IKLVL .GE. 3) CALL MPIXBCAST(ZYM3,N2ORBX,'DOUBLE',MASTER)
      IF (IKLVL .GE. 1) CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      IF (IKLVL .GE. 2) CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      IF (IKLVL .GE. 3) CALL MPIXBCAST(ISYMV3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMDN,1,'INTEGER',MASTER)
      CALL MPIXBCAST(OVLAP,1,'DOUBLE',MASTER)
C
C     Loop over all tesserae
C
      DO ITS = 1, NTSIRR
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(ITS,1,'INTEGER',NWHO,MPTAG2)
      END DO
C
C     Send end message to all slaves
C
      ITS = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(ITS ,1,'INTEGER',NWHO,MPTAG2)
      END DO
C
C     Collect relevant results from the slaves
C
      CALL DZERO(WORK,NTS)
      CALL MPI_REDUCE(WORK,TCHG,NTS,MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                MPI_COMM_WORLD,IERR)
C
      RETURN
      END
C /* Deck JX2S */
      SUBROUTINE J2XS(WRK,LWORK,IPRRSP)
C
C     Slave routine for calculating multiple one-index transformations
C     of PCM potentials, K.Ruud, Tromso June 2005
C
#include "implicit.h"
#include "dummy.h"
#include "iratdef.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "pcmdef.h"
C
      CHARACTER*8 LABINT(3*MXCOOR)
      LOGICAL CLCCHG, TOFILE, TRIMAT, EXP1VL
      DIMENSION WRK(LWORK)
#include "priunit.h"
#include "inforb.h"
#include "symmet.h"
#include "infdim.h"
#include "orgcom.h"
#include "mtags.h"
#include "pcm.h"
#include "infpar.h"
#include "mpif.h"
C
      CALL MPIXBCAST(IKLVL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C
      KDEN1  = 1
      KINTRP = KDEN1  + NASHDI*NASHDI
      KINTAD = KINTRP + (3*MXCOOR + 1)/IRAT
      KTCHG  = KINTAD + (3*MXCOOR + 1)/IRAT
      KCHGAO = KTCHG  + NTSIRR*(MAXREP + 1)
      KUCMO  = KCHGAO + NNBASX*(MAXREP + 1)
      KTELM  = KUCMO  + NORBT*NBAST
      KTLMA  = KTELM  + N2ORBX
      KTLMB  = KTLMA  + N2ORBX
      KZYM1  = KTLMB  + N2ORBX
      IF (IKLVL .GE. 1) THEN
         KZYM2 = KZYM1  + N2ORBX
      ELSE
         KZYM2 = KZYM1
      END IF
      IF (IKLVL .GE. 2) THEN
         KZYM3 = KZYM2  + N2ORBX
      ELSE
         KZYM3 = KZYM2
      END IF
      IF (IKLVL .GE. 3) THEN
         KLAST = KZYM3 + N2ORBX
      ELSE
         KLAST = KZYM3
      END IF
      LWRK = LWORK - KLAST + 1
C
      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KDEN1),NASHDI*NASHDI,'DOUBLE',MASTER)
      IF (IKLVL .GE. 1) CALL MPIXBCAST(WRK(KZYM1),N2ORBX,
     &                                 'DOUBLE',MASTER)
      IF (IKLVL .GE. 2) CALL MPIXBCAST(WRK(KZYM2),N2ORBX,
     &                                 'DOUBLE',MASTER)
      IF (IKLVL .GE. 3) CALL MPIXBCAST(WRK(KZYM3),N2ORBX,
     &                                 'DOUBLE',MASTER)
      IF (IKLVL .GE. 1) CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      IF (IKLVL .GE. 2) CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      IF (IKLVL .GE. 3) CALL MPIXBCAST(ISYMV3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMDN,1,'INTEGER',MASTER)
      CALL MPIXBCAST(OVLAP,1,'DOUBLE',MASTER)
C
      NSIM = 1
      IF(IKLVL.LT.0) THEN
         CLCCHG = .FALSE.
      ELSE
         CLCCHG = .TRUE.
      END IF
      IF (IKLVL.LE.0) ISYM = 1
      IF (IKLVL.EQ.1) ISYM = ISYMV1
      IF (IKLVL.EQ.2) ISYM = MULD2H(ISYMV1,ISYMV2)
      IF (IKLVL.EQ.3) ISYM = MULD2H(MULD2H(ISYMV1,ISYMV2),ISYMV3)
      ISYMT = MULD2H(ISYM,ISYMDN)
      CALL DZERO(WRK(KTCHG),NTSIRR*(MAXREP + 1))
C
 10   CONTINUE
      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(ITS,1,'INTEGER',MASTER,MPTAG2)
C
      IF (ITS .GT. 0) THEN
         CALL DZERO(WRK(KCHGAO),NNBASX*(MAXREP + 1))
         L = 1
         NCOMP = (MAXREP + 1)
         DIPORG(1) = XTSCOR(ITS)
         DIPORG(2) = YTSCOR(ITS)
         DIPORG(3) = ZTSCOR(ITS)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WRK(KCHGAO),'NPETES ',NCOMP,WRK(KLAST),LWRK,
     &               LABINT,WRK(KINTRP),WRK(KINTAD),L,TOFILE,KPATOM,
     &               TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
         JCHGAO = KCHGAO
         DO ILOP = 1, (MAXREP + 1)
            ISYM = ILOP
            JTS = (ILOP - 1)*NTSIRR + ITS
            CALL UTHU(WRK(JCHGAO),WRK(KTLMA),WRK(KUCMO),WRK(KLAST),
     $                NBAST,NORBT)
            CALL DSPTSI(NORBT,WRK(KTLMA),WRK(KTELM))
C with no transformation we copy the array from KTELM to KTLMA
            IF (IKLVL.EQ.0) THEN
               CALL DZERO(WRK(KTLMA),N2ORBX)
               CALL DCOPY(N2ORBX,WRK(KTELM),1,WRK(KTLMA),1)
            END IF
C First transformation of charges: q^e_{ab} --> q^e_{ab}({}^1\kappa)
            IF (IKLVL.GE.1) THEN
               CALL DZERO(WRK(KTLMA),N2ORBX)
               CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KTELM),WRK(KTLMA),ISYM)
               ISYM = MULD2H(ISYM,ISYMV1)
            END IF
C Second transformation of charges: q^e_{ab}({}^1\kappa) --> q^e_{ab}({}^1\kappa {}^2\kappa)
            IF (IKLVL.GE.2) THEN
               CALL DZERO(WRK(KTLMB),N2ORBX)
               CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYM)
               CALL DCOPY(N2ORBX,WRK(KTLMB),1,WRK(KTLMA),1)
               ISYM = MULD2H(ISYM,ISYMV2)
            END IF
C Third transformation of charges: hope you can figure out the formula.....
            IF (IKLVL.GE.3) THEN
               CALL DZERO(WRK(KTLMA),N2ORBX)
               CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KTLMB),WRK(KTLMA),ISYM)
               ISYM = MULD2H(ISYM,ISYMV3)
            END IF
C Contract transformed charges with the density
            CALL MELONE(WRK(KTLMA),ISYM,WRK(KDEN1),OVLAP,FACT,
     &                  200,'POPGET')
            KTS = KTCHG + JTS - 1
            WRK(KTS) = FACT
            JCHGAO = JCHGAO + NNBASX
         END DO
         GOTO 10
      END IF
C
C     Return relevant data to the master process
C
      CALL MPI_REDUCE(WRK(KTCHG),MPI_IN_PLACE,NTSIRR*(MAXREP + 1),
     &                MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                IERR)
      RETURN
      END
#endif
