!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
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

!**********************************************************************

      SUBROUTINE TRACI_CTL(vec1,vec2,c2,mo2mo_mat,par_dist_block_list,
     &                     block_list,rcctos,grouplist,proclist)
*
* Master routine for transforming a CI vector to new orbital basis
*
* Input vector on LUC
* output vector delivered on LUHC
* (LUSC1,LUSC2) used as scratch files
*
* Jeppe Olsen, January 98
! re-factored and parallelized by Stefan Knecht, November 2011
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLH1SAVE
!               for addressing of WORK
      real(8), intent(inout) :: vec1(*)
      real(8), intent(inout) :: vec2(*)
      real(8), intent(inout) :: c2(*)
      real(8), intent(inout) :: mo2mo_mat(*)
      integer, intent(in)    :: par_dist_block_list(*)
      integer, intent(in)    :: block_list(*)
      integer, intent(in)    :: proclist(*)
      integer, intent(in)    :: grouplist(*)
      integer, intent(in)    :: rcctos(*)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "clunit.inc"
#include "crun.inc"
#include "cstate.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "parluci.h"
      integer, parameter  :: isigden = 2
      real(8)             :: cdummy(2), sdummy(2)
      integer             :: my_in_fh, my_out_fh, my_sc1_fh, my_sc2_fh
      integer             :: len_ilu1
      integer             :: len_ilu2
      integer             :: len_iluc
      integer(kind=8)     :: my_in_off
      integer(kind=8)     :: my_out_off
      integer(kind=8)     :: my_scr_off
*
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
!     NTEST = 1000 ! debug
      NTEST = 0    ! default
      IDUM  = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACIC')
*
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(luwrt,*) ' Welcome to TRACI_CTL'
      WRITE(luwrt,*) ' ===================='
      call flshfo(luwrt)
#endif

*. Copy of one-electron integrals
      CALL MEMMAN(KLH1SAVE,NTOOB**2,'ADDL  ',2,'H1SAVE')
*. We are going to mess with the one-electron integrals, take a copy
      call dcopy(NTOOB*NTOOB,WORK(KINT1),1,WORK(KLH1SAVE),1)


      len_ilu1   = 0
      len_ilu2   = 0
      len_iluc   = 0
      my_in_off  = 0
      my_out_off = 0
      my_scr_off = 0
#ifdef VAR_MPI
      my_in_fh   = ilu1
      my_out_fh  = ilu1
      my_sc1_fh  = ilu2
      my_sc2_fh  = ilu2
      my_in_off  = 0
      my_out_off = 0
      my_scr_off = 0

      len_ilu1   = iall_lu1
      len_ilu2   = iall_lu2
      len_iluc   = iall_luc
#else
      my_in_fh   = lusc1
      my_out_fh  = luhc
      my_sc1_fh  = lusc2
      my_sc2_fh  = lusc3

      call copvcd(luc,my_in_fh,vec1,1,-1)
      call rewine(my_in_fh,-1)
      CALL rewine(my_out_fh,-1)
#endif
      lu_ref     = luc
      lu_refout  = luhc

!     transform CI vector
      call tracid(mo2mo_mat,work(kint1),lu_ref,lu_refout,
     &            my_in_fh,my_out_fh,
     &            my_sc1_fh,my_sc2_fh,
     &            vec1,vec2,c2,par_dist_block_list,block_list,
     &            rcctos,grouplist,proclist,
     &            len_ilu1,len_ilu2,len_iluc,
     &            my_in_off,my_out_off,my_scr_off)

#ifdef LUCI_DEBUG
      if(luci_myproc .eq. luci_master)then
        WRITE(luwrt,*)
        WRITE(luwrt,*) ' Analysis of rotated state'
        WRITE(luwrt,*) ' ========================='
        WRITE(luwrt,*)
        CALL WRTVCD(VEC1,luhc,1,-1)
      end if
#undef LUCI_DEBUG
#endif

!     clean up time : restore 1-e integrals
      CALL dcopy(NTOOB*NTOOB,WORK(KLH1SAVE),1,WORK(KINT1),1)

!     release local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACIC')

      END
!**********************************************************************

      SUBROUTINE TRACID(T,H,lu_ref,lu_refout,LUCIN,LUCOUT,LUSC1,LUSC2,
     &                  VEC1,VEC2,C2,par_dist_block_list,block_list,
     &                  rcctos,grouplist,proclist,
     &                  len_ilu_inout,len_ilu_scr,len_iluc,
     &                  my_in_off,my_out_off,my_scr_off)
*
* Transform CI vector on LUCIN with T matrix after
* Docent Malmquist's recipe. Place result as next vector on LUOUT
*
* The transformation is done as a sequence of one-electron transformations
*
* with each orbital transformation being
*
* Sum(k=0,2) ( 1/k! sum(n'.ne.n) S(n'n) E_{n'n} ) Tnn^N_n
*
* with Sn'n = T(n'n)/Tnn
*
      use lucita_file_list_pointer
#ifdef VAR_MPI
      use par_mcci_io
#ifndef VAR_USE_MPIF
      use mpi
      IMPLICIT REAL*8(A-H,O-Z)
#else
      IMPLICIT REAL*8(A-H,O-Z)
#include "mpif.h"
#endif
#endif
      INTEGER*8 KLCBLTPx
!               for addressing of WORK
#include "parluci.h"
#include "mxpdim.inc"
#include "oper.inc"
#include "intform.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cicisp.inc"
#include "cstate.inc"
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      REAL*8 INPRDD
*. Input
      real(8), intent(inout)      :: T(*)
      real(8), intent(inout)      :: H(*)
      real(8), intent(inout)      :: VEC1(*)
      real(8), intent(inout)      :: VEC2(*)
      real(8), intent(inout)      :: C2(*)
      integer, intent(in)         :: par_dist_block_list(*)
      integer, intent(in)         :: block_list(*)
      integer, intent(in)         :: proclist(*)
      integer, intent(in)         :: grouplist(*)
      integer, intent(in)         :: rcctos(*)
      integer, intent(in)         :: len_ilu_inout
      integer, intent(in)         :: len_ilu_scr
      integer, intent(in)         :: len_iluc
      integer(kind=8)             :: my_in_off
      integer(kind=8)             :: my_out_off
      integer(kind=8)             :: my_scr_off
      integer, parameter          :: isigden = 1
      real(8)                     :: cdummy(2), sdummy(2)
      integer, allocatable        :: lu_inoutlst(:)
      integer, allocatable        :: lu_scrlst(:)
      integer, allocatable        :: lu_clst(:)
      integer, allocatable        :: blocks_per_batch(:)
      integer, allocatable        :: batch_length(:)
      integer, allocatable        :: block_offset_batch(:)
      integer, allocatable        :: block_info_batch(:,:)
      integer, allocatable        :: blocktype(:)

!     NTEST = 1000 ! debug
      NTEST = 0000 ! default

!     a bit of info for the sigma routine
      I_RES_AB = 0

!     1-electron integrals in complete block form (NTOOB,NTOOB)
      IH1FORM  = 2

!     set up block and batch structure of vector

      allocate(lu_inoutlst(len_ilu_inout))
      allocate(lu_scrlst(len_ilu_scr))
      allocate(lu_clst(len_iluc))
      allocate(blocks_per_batch(mxntts))
      allocate(batch_length(mxntts))
      allocate(block_offset_batch(mxntts))
      allocate(block_info_batch(8,mxntts))
      allocate(blocktype(nsmst))

      lu_inoutlst        = 0
      lu_scrlst          = 0
      lu_clst            = 0
      blocks_per_batch   = 0
      batch_length       = 0
      block_offset_batch = 0
      block_info_batch   = 0
      blocktype          = 0

      IATP     = 1
      IBTP     = 2
      CALL Z_BLKFO_partitioning_parallel(ICSPC,ICSM,iatp,ibtp,
     &                                   blocks_per_batch,
     &                                   batch_length,  
     &                                   block_offset_batch,
     &                                   block_info_batch,
     &                                   NBATCH,NBLOCK,
     &                                   par_dist_block_list)

      CALL ZBLTP_IDC(ISMOST(1,ICSM),NSMST,IDC,blocktype)

#ifdef LUCI_DEBUG
#ifndef VAR_MPI
      WRITE(luwrt,*) ' input vector to ntoob loop'
      CALL WRTVCD(VEC1,lucin,1,-1)
      call rewine(lucin,-1)
      call flshfo(luwrt)
#endif
#endif

#ifdef VAR_MPI
      if(luci_nmproc > 1)then

!       sequential --> MPI I/O
!       ----------------------
        call rewino(lu_ref)
!       step 1: the rhs vector
        call mcci_cp_vcd_mpi_2_seq_io_interface(vec1,lu_ref,lucin,
     &                                          my_in_off,
     &                                          LU_inoutlst,
     &                                          par_dist_block_list,
     &                                          block_list,
     &                                          MPI_COMM_WORLD,
     &                                          NUM_BLOCKS,1,1,1)
      end if
#endif

!     transform each orbital separately
      DO K = 1, NTOOB

!       place (T(P,K)/S(K,K) in one-electron integral list work(kint1)
        CALL T_ROW_TO_H(T,H,K,TKK)

!       scale elements with T_{kk}^Nk depending on occupation type Nk (Nk=0,1,2)
        CALL T_TO_NK_VEC(TKK,K,ICSM,ICSPC,lucin,LUSC1,VEC1,
     &                   block_info_batch,blocktype,nblock,
     &                   my_in_off,LU_inoutlst)

!       for each orbital calculate (1+T+1/2 T^2)|0>
!       step 1: + T
!       -------------------------------------------
#ifdef VAR_MPI
        if(luci_nmproc .gt. 1)then
          CALL MPI_BARRIER(MYNEW_COMM,IERR)
!
!         reset LUCLIST
          call izero(lu_clst,len_iluc)

          call mcci_cp_vcd_batch(lucin,iluc,vec1,
     &                           nbatch,blocks_per_batch,
     &                           batch_length,  
     &                           block_offset_batch,
     &                           block_info_batch,
     &                           block_list,my_in_off,my_luc_off,
     &                           LU_inoutlst,lu_clst,
     &                           my_vec2_ioff,my_act_blk2,0)
!         set offset for sigma-file
          JVEC_SF = 0

          lusc_vector_file = iluc
          luhc_vector_file = lusc1
        end if
#else
        call copvcd(lusc1,lucin,vec1,1,-1)
        call rewine(lucin,-1)
        call rewine(lusc1,-1)

        lusc_vector_file = lucin
        luhc_vector_file = lusc1
#endif

        CALL SIGDEN_CI(VEC1,VEC2,C2,lusc_vector_file,lusc1,
     &                 cdummy,sdummy,ISIGDEN,-1
#ifdef VAR_MPI
     &                 ,lu_clst,LU_scrlst,block_list,
     &                 par_dist_block_list,rcctos,grouplist,proclist
#endif
     &                 )
        
#ifdef VAR_MPI
        if(luci_nmproc .gt. 1)then

          CALL VECSUM_PP_B_RL(VEC1,VEC2,
     &                        LU_scrlst,  
     &                        LU_inoutlst,
     &                        NBATCH,blocks_per_batch,
     &                        batch_length,
     &                        block_offset_batch,
     &                        block_info_batch,
     &                        my_scr_off,
     &                        my_in_off,
     &                        1,
     &                        lusc1,
     &                        lucin,
     &                        1.0d0,
     &                        1.0d0)
        end if
#else

!       y (lusc2) := 1.0d0*lucin + 1.0d0*lusc1
        CALL VECSMD(VEC1,VEC2,1.0d0,1.0d0,lucin,LUSC1,lusc2,1,-1)
        CALL COPVCD(lusc2,lucin,VEC1,1,-1)
        CALL rewine(lusc2,-1)
        CALL rewine(lusc1,-1)

#ifdef LUCI_DEBUG
        CALL WRTVCD(VEC1,lucin,1,-1)
        call rewine(lucin,-1)
#endif

#endif

!       step 2: + 1/2 T^2
!       -----------------

#ifdef VAR_MPI
        if(luci_nmproc .gt. 1)then
          CALL MPI_BARRIER(MYNEW_COMM,IERR)
!
!         reset LUCLIST
          call izero(lu_clst,len_iluc)

          call mcci_cp_vcd_batch(lusc1,iluc,vec1,
     &                           NBATCH,blocks_per_batch,
     &                           batch_length,
     &                           block_offset_batch,
     &                           block_info_batch,
     &                           block_list,my_scr_off,my_luc_off,
     &                           LU_scrlst,lu_clst,
     &                           my_vec2_ioff,my_act_blk2,0)
!         set offset for sigma-file
          JVEC_SF = 0

          lusc_vector_file = iluc
          luhc_vector_file = lusc1
        end if
#else
        lusc_vector_file = lusc1
        luhc_vector_file = lusc2
#endif

        CALL SIGDEN_CI(VEC1,VEC2,C2,lusc_vector_file,luhc_vector_file,
     &                 cdummy,sdummy,ISIGDEN,-1
#ifdef VAR_MPI
     &                 ,lu_clst,LU_scrlst,block_list,
     &                 par_dist_block_list,rcctos,grouplist,proclist
#endif
     &           )

#ifdef VAR_MPI
        if(luci_nmproc .gt. 1)then
          CALL VECSUM_PP_B_RL(VEC1,VEC2,
     &                        LU_scrlst,
     &                        LU_inoutlst,
     &                        NBATCH,blocks_per_batch,
     &                        batch_length,
     &                        block_offset_batch,
     &                        block_info_batch,
     &                        my_scr_off,
     &                        my_out_off,
     &                        1,
     &                        lusc1,
     &                        lucout,
     &                        0.5d0,
     &                        1.0d0)
        end if
#else

!       y (lusc1) := 1.0d0*lucin + 0.5d0*lusc2
        CALL VECSMD(VEC1,VEC2,1.0d0,0.5d0,lucin,lusc2,LUSC1,1,-1)

#ifdef LUCI_DEBUG
        WRITE(luwrt,*) ' 2nd vecsmd: in 1'
        CALL WRTVCD(VEC1,lucin,1,-1)
        WRITE(luwrt,*) ' 2nd vecsmd: in 2'
        CALL WRTVCD(VEC1,lusc2,1,-1)
        WRITE(luwrt,*) ' 2nd vecsmd: out'
        CALL WRTVCD(VEC1,lusc1,1,-1)
        call rewine(lusc1,-1)
        call flshfo(luwrt)
#endif

!       transfer back to lucin
        CALL COPVCD(LUSC1,lucin,VEC1,1,-1)
        CALL rewine(lucin,-1)
#endif

      END DO

#ifdef VAR_MPI
      if(luci_nmproc > 1)then

!       collect the rotated vector
!       --------------------------
        call rewine(lu_refout,-1)

!       the lhs vector
        call mcci_cp_vcd_mpi_2_seq_io_interface(vec2,lu_refout,
     &                                          lucout,
     &                                          my_out_off,
     &                                          LU_inoutlst,
     &                                          par_dist_block_list,
     &                                          block_list,
     &                                          MPI_COMM_WORLD,
     &                                          num_blocks,1,1,2)

      end if ! luci_nmproc > 1
#else

!     transfer rotated vector to lu_refout
      CALL COPVCD(lucin,lu_refout,VEC1,1,-1)
#endif

#ifdef LUCI_DEBUG
      if(luci_myproc .eq. luci_master)then
        CNORM = INPRDD(VEC1,VEC2,lu_refout,lu_refout,1,-1)
        WRITE(luwrt,*) ' Norm of transformed vector', CNORM
      end if
#endif

      deallocate(lu_inoutlst)
      deallocate(lu_scrlst)
      deallocate(lu_clst)
      deallocate(blocks_per_batch)
      deallocate(batch_length)
      deallocate(block_offset_batch)
      deallocate(block_info_batch)
      deallocate(blocktype)

      END
!**********************************************************************

      SUBROUTINE T_ROW_TO_H(T,H,K,TKK)
*
* Set H integrals
*
*    Column K : H(P,K) = T(P,K)/T(K,K), P.NE.K
*    Other Columns     = 0
* - and return T_{kk} in TKK
*
*
* Jeppe Olsen, Jan 98
* For rotation of CI vectors
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "glbbas.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "parluci.h"
*. Input ( in blocked form)
      DIMENSION T(*)
*. Output ( also in blocked form)
      DIMENSION H(*)
*
      CALL DZERO(H,NTOOB**2)

      KSM  = ISMFSO(K)
      KOFF = IBSO(KSM)
      KREL = K - KOFF + 1
      NK   = NTOOBS(KSM)
!     symmetry offset in full 1e-matrix block
      IOFF = IFRMR(WORK(KPGINT1A(1)),1,KSM)

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(luwrt,*) 'ioff in T and H ==> ',IOFF+(KREL-1)*NK
      write(luwrt,*) 'NK elements     ==> ',NK
      write(luwrt,*) 'absolute ioff   ==> ',ioff
      write(luwrt,*) 'krel, ksm       ==> ',krel,ksm
      call flshfo(luwrt)
#endif
!     copy orbital block with symmetry ksm
      CALL DCOPY(NK,T(IOFF+(KREL-1)*NK),1,H(IOFF+(KREL-1)*NK),1)

!     set T_kk
      TKK = H(IOFF-1+(KREL-1)*NK+KREL)
      IF(TKK .NE. 0.0D0) THEN
        FAC = 1.0D0/TKK
        CALL DSCAL(NK,FAC,H(IOFF+(KREL-1)*NK),1)
        H(IOFF-1+(KREL-1)*NK+KREL) = 0.0D0
      END IF

#ifdef LUCI_DEBUG
      WRITE(luwrt,*) ' output from T_ROW_TO_H, TKK ==> ',TKK
      WRITE(luwrt,*) ' Row to be transferred ', KREL
      WRITE(luwrt,*) ' Updated H matrix '
      WRITE(luwrt,*) ' NTOOB, NTOOBS(1)', NTOOB, NTOOBS(1)
      CALL APRBLM2(H,NTOOBS,NTOOBS,NSMOB,0)
#undef LUCI_DEBUG
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C,
     &                       block_info,block_type_info,nblock,
     &                       my_in_off,LU_inoutlst)
*
* Evaluate T**(NK_operator) times vector on file LUIN
* to yield vector on file LUOUT
* (NK_operator is number operator for orbital K )
*
* Note LUCIN and LUCOUT are both rewinded before read/write
* Input
* =====
*  T : Input constant
*  KORB : Orbital in symmetry order
*
*  ISM,ISPC : Symmetry and space of state on LUIN
*  C : Scratch block
*
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLASTR, KLBSTR, KLKAOC, KLKBOC
!               for addressing of WORK
#include "priunit.h"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strinp.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "crun.inc"
#include "csm.inc"
      integer, intent(in)         :: nblock
      integer, intent(in)         :: block_type_info(*)
      integer, intent(in)         :: block_info(8,nblock)
      integer(kind=8), intent(in) :: my_in_off
      integer, intent(inout)      :: lu_inoutlst(*)

*. Scratch block, must hold a batch of blocks
      DIMENSION C(*)
*
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(lupri,*) ' T_TO_NK_VEC speaking '
      WRITE(lupri,*) ' ISM, ISPC = ', ISM,ISPC
      call flshfo(lupri)
#endif

      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'T_TO_N')

      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
      CALL MEMMAN(KLKAOC,MXNSTR,     'ADDL  ',1,'KLKAOC')
      CALL MEMMAN(KLKBOC,MXNSTR,     'ADDL  ',1,'KLKBOC')
      
!     orbital K in type ordering
      KKORB = IREOST(KORB)
#ifdef LUCI_DEBUG
      write(lupri,*) 'korb, kkorb ==> ',korb, kkorb
      call flshfo(lupri)
#endif

      CALL T_TO_NK_VECS(T,KKORB,C,LUCIN,LUCOUT,
     &                  WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &                  NBLOCK,block_info,
     &                  NAEL,NBEL,WORK(KLASTR),WORK(KLBSTR),
     &                  block_type_info,NSMST,
     &                  ICISTR,NTOOB,WORK(KLKAOC),WORK(KLKBOC),
     &                  my_in_off,lu_inoutlst)

      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'T_TO_N')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE T_TO_NK_VECS(T,KORB,C,LUCIN,LUCOUT,NSSOA,NSSOB,
     &                        NBLOCK,IBLOCK,
     &                        NAEL,NBEL,IASTR,IBSTR,IBLTP,NSMST,
     &                        ICISTR,NORB,IKAOCC,IKBOCC,
     &                        my_in_off,lu_inoutlst)
*
* Multiply Vector in LUCIN with t **NK_op to yield vector on LUCOUT
*
* Both files are initially rewinded
*
*
* Jeppe Olsen, Feb. 1998
*

#ifdef VAR_MPI
#ifndef VAR_USE_MPIF
      use mpi
      IMPLICIT REAL*8(A-H,O-Z)
#else
      IMPLICIT REAL*8(A-H,O-Z)
#include "mpif.h"
#endif
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#else
      IMPLICIT REAL*8(A-H,O-Z)
#endif
#include "maxorb.h"
#include "infpar.h"
#include "parluci.h"
#include "priunit.h"
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
*. Scratch
      DIMENSION C(*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IKAOCC(*),IKBOCC(*)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
      DIMENSION IBLTP(*)
      integer, intent(in)         :: nblock
      integer, intent(inout)      :: lu_inoutlst(*)
      integer(kind=8), intent(in) :: my_in_off

      integer(kind=8)             :: my_internal_in_off
      integer                     :: my_active_block
*
#ifdef VAR_MPI
      my_internal_in_off = 0
      my_internal_in_off = my_internal_in_off + my_in_off
      my_active_block    = 0
#else
      CALL REWINE(LUCIN,-1)
      CALL REWINE(LUCOUT,-1)
#endif
*
!
      T2 = T**2
      DO JBLOCK = 1, NBLOCK

        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
        WRITE(lupri,*) ' IATP IBTP IASM IBSM ', IATP,IBTP,IASM,IBSM
        WRITE(lupri,*) ' IBLTP(IASM)         ', IBLTP(IASM) 
        call flshfo(lupri)
#undef LUCI_DEBUG
#endif

        if(iatp.eq.0)then
          if(luci_nmproc .gt. 1)then
            cycle
          else
            IF(ICISTR .GE. 2)THEN
!             Read in a Type-Type-symmetry block
              CALL IFRMDS(LDET,1,-1,LUCIN)
              CALL FRMDSC_LUCI(C,LDET,-1,LUCIN,IMZERO,IAMPACK)
            END IF
            goto 999
          end if
        end if


!       obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
!       occupation of orb KORB
        DO JSTR = 1, NASTR1
          KOCC = 0
          DO JAEL = 1, NAEL
            IF(IASTR(JAEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKAOCC(JSTR) = KOCC
        END DO

#ifdef LUCI_DEBUG
        WRITE(lupri,*) ' IKAOCC array '
        CALL IWRTMAMN(IKAOCC,1,NASTR1,1,NASTR1,lupri)
        call flshfo(lupri)
#endif

!       Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                           NORB,0,IDUM,IDUM)
        DO JSTR = 1, NBSTR1
          KOCC = 0
          DO JBEL = 1, NBEL
            IF(IBSTR(JBEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKBOCC(JSTR) = KOCC
        END DO

#ifdef LUCI_DEBUG
        WRITE(lupri,*) ' IKBOCC array '
        CALL IWRTMAMN(IKBOCC,1,NBSTR1,1,NBSTR1,lupri)
        call flshfo(lupri)
#endif

        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF

        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)

!       WRITE(lupri,*) ' NIA NIB ', NIA,NIB

        IMZERO = 0
        if(luci_nmproc .gt. 1)then
#ifdef VAR_MPI
          my_active_block = my_active_block + 1

          if(lu_inoutlst(my_active_block) .eq. 0)then
            imzero = 1
            ldet   = 0
          else
            imzero = 0
            ldet   = iblock(7,jblock)
            CALL MPI_FILE_READ_AT(LUCIN,my_internal_in_off,C,LDET,
     &                            MPI_REAL8,my_STATUS,IERR)
!           new offset (not yet...)
            my_internal_in_off = my_internal_in_off
          end if
#endif
        else
          IF(ICISTR .GE. 2)THEN
*. Read in a Type-Type-symmetry block
            CALL IFRMDS(LDET,1,-1,LUCIN)
            CALL FRMDSC_LUCI(C,LDET,-1,LUCIN,IMZERO,IAMPACK)
          END IF
        end if

        IF(IMZERO.NE.1) THEN
*
          IDET = 0
          DO  IB = 1,NIB
            IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
              MINIA = IB
            ELSE
              MINIA = 1
            END IF
            DO  IA = MINIA,NIA

              IDET = IDET + 1
              KABOCC = IKAOCC(IA)+IKBOCC(IB)

#ifdef LUCI_DEBUG
              WRITE(lupri,*) ' KABOCC, C(IDET), T, T**2',
     &                         KABOCC, C(IDET), T, T**2
              WRITE(lupri,*) ' IA IB IDET',IA,IB,IDET
              call flshfo(lupri)
#endif

              IF(KABOCC.EQ.1) THEN
                C(IDET) = T*C(IDET)
              ELSE IF(KABOCC.EQ.2) THEN
                C(IDET) = T2 *C(IDET)
              END IF
            END DO
*           ^ End of loop over alpha strings
          END DO
*         ^ End of loop over beta strings
*
        END IF
*       ^ End of if statement for nonvanishing blocks

 999    if(luci_nmproc .gt. 1)then

#ifdef VAR_MPI
          CALL MPI_FILE_WRITE_AT(LUCIN,my_internal_in_off,C,LDET,
     &                           MPI_REAL8,my_STATUS,IERR)
#endif
!         new offset (now...)
          my_internal_in_off = my_internal_in_off + ldet
        else
          CALL ITODS(LDET,1,-1,LUCOUT)
          CALL TODSC_LUCI(C,LDET,-1,LUCOUT)
        end if

      END DO
*     ^ End of loop over blocks

#ifndef VAR_MPI
*. This is the end, the end of every file my friend, the end
      CALL ITODS(-1,1,-1,LUCOUT)
#endif

      END
***********************************************************************

      SUBROUTINE TRAH1(NBAS,NORB,NSM,HAO,C,HMO,IHSM,SCR)
*
*. Transform one-electron integrals from ao's to mo's.
*
*. Symmetry of integrals is IHSM, all integrals blocks assumed complete,
* i.e not packed to lower half
*
* Jeppe Olsen
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION HAO(*),C(*)
      DIMENSION NORB(*),NBAS(*)
#include "multd2h.inc"

*. Output
      DIMENSION  HMO(*)
*. Scratch
      DIMENSION SCR(*)
*. Loop over integral blocks
      IBHAO = 1
      IBHMO = 1
      DO IRSM = 1, NSM
        ICSM = MULTD2H(IRSM,IHSM)
*. Pointers to offsets in transformation matrices
        IBR = 1
        DO ISM = 1, IRSM-1
          IBR = IBR + NORB(ISM)*NBAS(ISM)
        END DO
        IBC = 1
        DO ISM = 1, ICSM-1
          IBC = IBC + NORB(ISM)*NBAS(ISM)
        END DO
*.
        LRMO = NORB(IRSM)
        LRAO = NBAS(IRSM)
*
        LCMO = NORB(ICSM)
        LCAO = NBAS(ICSM)
C       write(6,*) ' TRAH1 : IRSM ICSM ',IRSM,ICSM
C       WRITE(6,*) ' LRAO LRMO LCAO LCMO ',LRAO,LRMO,LCAO,LCMO

*
C            MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &             NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
        ZERO = 0.0D0
        ONE= 1.0D0
*.C(row)T*Hao
        CALL MATML7(SCR,C(IBR),HAO(IBHAO),
     &       LRMO,LCAO,LRAO,LRMO,LRAO,LCAO,ZERO,ONE,1)
*. (C(row)T*Hao)*C(column)
        CALL MATML7(HMO(IBHMO),SCR,C(IBC),
     &       LRMO,LCMO,LRMO,LCAO,LCAO,LCMO,ZERO,ONE,0)
*
        IBHAO = IBHAO + LRAO*LCAO
        IBHMO = IBHMO + LRMO*LCMO
*.
      END DO
*
      RETURN
      END
!**********************************************************************

      SUBROUTINE TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
*
* Transform a blocked symmetric matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LBLOCK(NBLOCK)
*. Output
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block
      DIMENSION SCR(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP = 1
         IOFFC = 1
       ELSE
         IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
         IOFFC = IOFFC + LBLOCK(IBLOCK-1)** 2
       END IF
       L = LBLOCK(IBLOCK)
       K1 = 1
       K2 = 1 + L **2
*. Unpack block of A
C      TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       SIGN = 1.0D0
       CALL TRIPAK_LUCI(SCR(K1),AIN(IOFFP),2,L,L,SIGN)
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL MATML7(SCR(K2),X(IOFFC),SCR(K1),L,L,L,L,L,L,
     &             ZERO,ONE,1)
*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFC),L,L,L,L,L,L,
     &             ZERO,ONE,0)
*. Pack and transfer
       CALL TRIPAK_LUCI(SCR(K1),AOUT(IOFFP),1,L,L,SIGN)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        CALL APRBLM2(AOUT,LBLOCK,LBLOCK,NBLOCK,1)
      END IF
*
      RETURN
      END
!**********************************************************************

      SUBROUTINE TRAN_SYM_BLOC_MAT2(AIN,X,NBLOCK,LBLOCK,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LBLOCK(NBLOCK)
*. Output
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block
      DIMENSION SCR(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP = 1
         IOFFC = 1
       ELSE
         IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
         IOFFC = IOFFC + LBLOCK(IBLOCK-1)** 2
       END IF
       L = LBLOCK(IBLOCK)
       K1 = 1
       K2 = 1 + L **2
*. Unpack block of A
C      TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AIN(IOFFP),2,L,L,SIGN)
       ELSE
         CALL COPVEC(AIN(IOFFC),SCR(K1),L*L)
       END IF
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL MATML7(SCR(K2),X(IOFFC),SCR(K1),L,L,L,L,L,L,
     &             ZERO,ONE,1)
*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFC),L,L,L,L,L,L,
     &             ZERO,ONE,0)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AOUT(IOFFP),1,L,L,SIGN)
       ELSE
         CALL COPVEC(SCR(K1),AOUT(IOFFC),L*L)
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        CALL APRBLM2(AOUT,LBLOCK,LBLOCK,NBLOCK,ISYM)
      END IF
*
      RETURN
      END
!**********************************************************************

      SUBROUTINE TRAN_SYM_BLOC_MAT3
     &   (AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LX_ROW(NBLOCK),LX_COL(NBLOCK)
*. Output
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block
      DIMENSION SCR(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP_IN = 1
         IOFFC_IN = 1
         IOFFP_OUT= 1
         IOFFC_OUT= 1
         IOFFX = 1
       ELSE
         IOFFP_IN = IOFFP_IN + LX_ROW(IBLOCK-1)*(LX_ROW(IBLOCK-1)+1)/2
         IOFFC_IN = IOFFC_IN + LX_ROW(IBLOCK-1) ** 2
         IOFFP_OUT= IOFFP_OUT+ LX_COL(IBLOCK-1)*(LX_COL(IBLOCK-1)+1)/2
         IOFFC_OUT= IOFFC_OUT+ LX_COL(IBLOCK-1) ** 2
         IOFFX = IOFFX + LX_ROW(IBLOCK-1)*LX_COL(IBLOCK-1)
       END IF
       LXR = LX_ROW(IBLOCK)
       LXC = LX_COL(IBLOCK)
       K1 = 1
       K2 = 1 + MAX(LXR,LXC) ** 2
*. Unpack block of A
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AIN(IOFFP_IN),2,LXR,LXR,SIGN)
C             TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       ELSE
         CALL COPVEC(AIN(IOFFC_IN),SCR(K1),LXR*LXR)
       END IF
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL SETVEC(SCR(K2),ZERO,LXR*LXC)
       CALL MATML7(SCR(K2),X(IOFFX),SCR(K1),LXC,LXR,LXR,LXC,LXR,LXR,
     &             ZERO,ONE,1)
C?     WRITE(6,*) ' Half transformed matrix '
C?     CALL WRTMT_LU(SCR(K2),LXC,LXR,LXC,LXR)

*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL SETVEC(SCR(K1),ZERO,LXC*LXC)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFX),LXC,LXC,LXC,LXR,LXR,LXC,
     &             ZERO,ONE,0)
C?     WRITE(6,*) ' Transformed matrix '
C?     CALL WRTMT_LU(SCR(K1),LXC,LXC,LXC,LXC)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AOUT(IOFFP_OUT),1,LXC,LXC,SIGN)
       ELSE
         CALL COPVEC(SCR(K1),AOUT(IOFFC_OUT),LXC*LXC)
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL APRBLM2(AOUT,LX_COL,LX_COL,NBLOCK,ISYM)
      END IF
*
      END
***********************************************************************

      SUBROUTINE SCLH2(XLAMBDA)
*
*. Scale two electron integrals
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      call dscal(nint2,xlambda,work(kint2),1)
*
      END

***********************************************************************
      SUBROUTINE TRA1_SIMPLE(C)
*
* Transform one-electron integrals
*
* Integrals from KINT1O used : Integrals without particle-hole corrections
*
* Input integrals in KINT1
* Output integrals also in KINT1
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLH1, KLSCR
!               for addressing of WORK
*. General Input
#include "mxpdim.inc"
#include "orbinp.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Specific Input
      DIMENSION C(*)
*
      NTEST = 100
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA1_S')
*. Largest symmetry block of orbitals
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Largest number of orbitals in sym block ', MXSOB
      END IF
*. Space for 1e- integrals and two sym blocks
      LENH1 = NTOOB ** 2
      CALL MEMMAN(KLH1,LENH1,'ADDL  ',2,'LH1   ')
      
      LENSCR = 2 * MXSOB ** 2
      CALL MEMMAN(KLSCR,LENSCR,'ADDL   ',2,'H1SCR ')
      WRITE(6,*) ' LENH1 LENSCR ', LENH1, LENSCR
*. and do it
C     TRAN_SYM_BLOC_MAT3 (AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      WRITE(6,*) ' Before TRAN_SYM '
      CALL COPVEC(WORK(KINT1O),WORK(KINT1),NINT1)
      CALL TRAN_SYM_BLOC_MAT3(WORK(KINT1),C,NSMOB,NTOOBS,NTOOBS,
     &                        WORK(KLH1),WORK(KLSCR),1)
      WRITE(6,*) ' After TRAN_SYM '
*. Copy to original array
      CALL COPVEC(WORK(KLH1),WORK(KINT1),NINT1)
*. Flush memory
      WRITE(6,*) ' Returning from TRA1 '
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'TRA1_S')
*
      RETURN
      END
!**********************************************************************

      SUBROUTINE TRA2_SIMPLE(C)
*
* Trivial 2-electron integral transformation
*
*. Input integrals in KINT2
*. Output integrals in KINT2
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLBL1, KLBL2
!               for addressing of WORK
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "multd2h.inc"
*. Specific input
      DIMENSION C(*)
*
      NTEST  = 0
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA2_S')
*. Largest symmetry block
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Two symmetry blocks
      LENBL = MXSOB ** 4
      IF(NTEST.GE.100)
     &WRITE(6,*) ' Size of symmetry block ', LENBL
      CALL  MEMMAN(KLBL1,LENBL,'ADDL ',2,'KLBL1 ')
      CALL  MEMMAN(KLBL2,LENBL,'ADDL ',2,'KLBL2 ')
C     WRITE(6,*) ' Two integral blocks allocated '
*. Loop over symmetry blocks
      DO ISM = 1, NSMOB
        DO JSM = 1, ISM
          DO KSM = 1, ISM
            IF(KSM.EQ.ISM) THEN
              LSM_MX = JSM
            ELSE
              LSM_MX = KSM
            END IF
            DO LSM = 1, LSM_MX
*. Ensure that integrals have correct symmetry
              INTSYM = 1
              IJSM = MULTD2H(ISM,JSM)
              IJKSM = MULTD2H(IJSM,KSM)
              IF(INTSYM.EQ.MULTD2H(IJKSM,LSM)) THEN
*. Fetch 2-e integral block (IJ!KL)
C             LGETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &               IXCHNG,IKSM,JLSM,ICOUL)
C             WRITE(6,*)  'ISM, JSM, KSM, LSM ',ISM,KSM,KSM,LSM
              CALL LGETINT(WORK(KLBL1),0,ISM,0,JSM,0,KSM,0,LSM,
     &                     0,0,0,1)
*. (Type = 0  => Complete symmetryblock)
*. Offsets and dimensions for symmetryblocks in C
              IOFF = 1
              DO IISM = 1, ISM-1
                IOFF = IOFF + NTOOBS(IISM)**2
              END DO
              NI = NTOOBS(ISM)
*
              JOFF = 1
              DO JJSM = 1, JSM-1
                JOFF = JOFF + NTOOBS(JJSM)**2
              END DO
              NJ = NTOOBS(JSM)
*
              KOFF = 1
              DO KKSM = 1, KSM-1
                KOFF = KOFF + NTOOBS(KKSM)**2
              END DO
              NK = NTOOBS(KSM)
*
              LOFF = 1
              DO LLSM = 1, LSM-1
                LOFF = LOFF + NTOOBS(LLSM)**2
              END DO
              NL = NTOOBS(LSM)
*. Transform 2-electron integral block
              CALL TRA_2EL_BLK_SIMPLE(WORK(KLBL1),
     &        C(IOFF),NI,C(JOFF),NJ,C(KOFF),NK,C(LOFF),NL,WORK(KLBL2))
*. Transfer symmetry block to integral list
              CALL PUTINT(WORK(KLBL1),0,ISM,0,JSM,0,KSM,0,LSM)
            END IF
*           ^ Check if integrals have correct symmetry
            END DO
          END DO
        END DO
      END DO
*     ^ End of loop over symmetries
*. Dump integrals to file LU90
      LU90 = 90
      CALL DMPINT(LU90)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRA2_S')
*
      RETURN
      END
!**********************************************************************

      SUBROUTINE TRA_2EL_BLK_SIMPLE(XINT,CI,NI,CJ,NJ,CK,NK,CL,NL,SCR)
*
* Transform 2-electron integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CI(NI,NI),CJ(NJ,NJ),CK(NK,NK),CL(NL,NL)
*. Input and output
      DIMENSION XINT(*)
* Matrix given in complete form XI(NI,NJ,NK,NL)
*. Scratch
      DIMENSION SCR(*)
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input integral block as X(IJ,KL)'
        CALL WRTMT_LU(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*. Transform first two indeces
      DO K = 1, NK
       DO L = 1, NL
         KLOF = 1 + ((L-1)*NK+K-1)*NI*NJ
         FACTORC = 0.0D0
         FACTORAB = 1.0D0
         CALL MATML7(SCR(KLOF),XINT(KLOF),CJ,NI,NJ,NI,NJ,NJ,NJ,
     &                FACTORC,FACTORAB,0)
         CALL MATML7(XINT(KLOF),CI,SCR(KLOF),NI,NJ,NI,NI,NI,NJ,
     &                FACTORC,FACTORAB,1)
       END DO
      END DO
*. Transpose Block X(IJ,KL)
      NIJ = NI*NJ
      NKL = NK*NL
      CALL TRPMT3(XINT,NIJ,NKL,SCR)
*. Matrix is now SCR(KL,IJ)
*. Transform last two indeces
      DO I = 1, NI
        DO J = 1, NJ
         IJOF = ((J-1)*NI+I-1)*NK*NL+1
         CALL MATML7(XINT(IJOF),SCR(IJOF),CL,NK,NL,NK,NL,NL,NL,
     &               FACTORC,FACTORAB,0)
         CALL MATML7(SCR(IJOF),CK,XINT(IJOF),NK,NL,NK,NK,NK,NL,
     &               FACTORC,FACTORAB,1)
        END DO
      END DO
*. Transpose to obtain XINT(IJ,KL)
      CALL TRPMT3(SCR,NKL,NIJ,XINT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output integral block as X(IJ,KL)'
        CALL WRTMT_LU(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END
