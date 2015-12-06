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

!**********************************************************************

      SUBROUTINE TRACID(T,H,lu_ref,lu_refout,LUCIN,LUCOUT,LUSC1,LUSC2,
     &                  BVC_fh,VEC1,VEC2,C2,
     &                  NBATCH,NBLOCK,blocks_per_batch,
     &                  batch_length,block_offset_batch,
     &                  block_info_batch,blocktype,
     &                  par_dist_block_list,block_list,
     &                  rcctos,grouplist,proclist,
     &                  lu_inoutlst,lu_scrlst,lu_BVClst,
     &                  len_ilu_inout,len_ilu_scr,len_iBVC,
     &                  my_in_off,my_out_off,my_scr_off,my_BVC_off)
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
#ifdef VAR_MPI
      use par_mcci_io
#ifdef USE_MPI_MOD_F90
      use mpi
      IMPLICIT REAL*8(A-H,O-Z)
#else
      IMPLICIT REAL*8(A-H,O-Z)
#include "mpif.h"
#endif
#endif
#include "parluci.h"
#include "mxpdim.inc"
#include "oper.inc"
#include "intform.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "cstate.inc"
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      REAL*8 INPRDD
*. Input
      real(8), intent(inout) :: T(*)
      real(8), intent(inout) :: H(*)
      real(8), intent(inout) :: VEC1(*)
      real(8), intent(inout) :: VEC2(*)
      real(8), intent(inout) :: C2(*)
      integer, intent(in)    :: nblock
      integer, intent(in)    :: nbatch
      integer, intent(in)    :: blocks_per_batch(*)
      integer, intent(in)    :: batch_length(*)
      integer, intent(in)    :: block_offset_batch(*)
      integer, intent(in)    :: block_info_batch(8,*)
      integer, intent(in)    :: blocktype(*)
      integer, intent(in)    :: par_dist_block_list(*)
      integer, intent(in)    :: block_list(*)
      integer, intent(in)    :: proclist(*)
      integer, intent(in)    :: grouplist(*)
      integer, intent(in)    :: rcctos(*)
      integer, intent(in)    :: BVC_fh
      integer, intent(in)    :: len_ilu_inout
      integer, intent(in)    :: len_ilu_scr
      integer, intent(in)    :: len_iBVC
      integer, intent(inout) :: lu_inoutlst(*)
      integer, intent(inout) :: lu_scrlst(*)
      integer, intent(inout) :: lu_BVClst(*)
      integer(kind=8)        :: my_in_off
      integer(kind=8)        :: my_out_off
      integer(kind=8)        :: my_scr_off
      integer(kind=8)        :: my_BVC_off
      integer, parameter     :: isigden = 1
      real(8)                :: cdummy(2), sdummy(2)
      logical                :: skip_orbital

!     NTEST = 1000 ! debug
      NTEST = 0000 ! default

!     a bit of info for the sigma routine
      I_RES_AB = 0

!     1-electron integrals in complete block form (NTOOB,NTOOB)
      IH1FORM  = 2

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
#ifndef VAR_MPI
      WRITE(luwrt,*) ' input vector to ntoob loop'
      CALL WRTVCD(VEC1,lucin,1,-1)
      call rewine(lucin,-1)
      call flshfo(luwrt)
#endif
#endif

#ifdef BLABLA_DOIT
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
#endif

      skip_orbital = .false.

!     transform each orbital separately
      DO K = 1, NTOOB

!       place (T(P,K)/S(K,K) in one-electron integral list work(kint1)
        CALL T_ROW_TO_H(T,H,K,TKK,skip_orbital)

!       proceed to next orbital if ddot(H(row)) equals 0.0d0
        if(skip_orbital) cycle

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
          call izero(lu_BVClst,len_iBVC)

          call mcci_cp_vcd_batch(lucin,BVC_fh,vec1,
     &                           nbatch,blocks_per_batch,
     &                           batch_length,  
     &                           block_offset_batch,
     &                           block_info_batch,
     &                           block_list,
     &                           my_in_off,
     &                           my_BVC_off,
     &                           LU_inoutlst,
     &                           lu_BVClst,
     &                           my_vec2_ioff,my_act_blk2,0)
!         set offset for sigma-file
          JVEC_SF = 0

          lusc_vector_file = BVC_fh
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
     &                 cdummy,sdummy,ISIGDEN,-1,xxs2dm
#ifdef VAR_MPI
     &                 ,lu_BVClst,
     &                 LU_scrlst,
     &                 block_list,
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

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
        WRITE(luwrt,*) ' 1st vecsmd: in 1'
        CALL WRTVCD(VEC1,lucin,1,-1)
        WRITE(luwrt,*) ' 1st vecsmd: in 2'
        CALL WRTVCD(VEC1,lusc1,1,-1)
        WRITE(luwrt,*) ' 1st vecsmd: out'
        CALL WRTVCD(VEC1,lusc2,1,-1)
        call rewine(lusc2,-1)
        call rewine(lusc1,-1)
        call rewine(lucin,-1)
        call flshfo(luwrt)
#endif

#endif

!       step 2: + 1/2 T^2
!       -----------------

#ifdef VAR_MPI
        if(luci_nmproc .gt. 1)then
          CALL MPI_BARRIER(MYNEW_COMM,IERR)
!
!         reset LUCLIST
          call izero(lu_BVClst,len_iBVC)

          call mcci_cp_vcd_batch(lusc1,BVC_fh,vec1,
     &                           NBATCH,blocks_per_batch,
     &                           batch_length,
     &                           block_offset_batch,
     &                           block_info_batch,
     &                           block_list,
     &                           my_scr_off,
     &                           my_BVC_off,
     &                           LU_scrlst,
     &                           lu_BVClst,
     &                           my_vec2_ioff,my_act_blk2,0)
!         set offset for sigma-file
          JVEC_SF = 0

          lusc_vector_file = BVC_fh
          luhc_vector_file = lusc1
        end if
#else
        lusc_vector_file = lusc1
        luhc_vector_file = lusc2
#endif

        CALL SIGDEN_CI(VEC1,VEC2,C2,lusc_vector_file,luhc_vector_file,
     &                 cdummy,sdummy,ISIGDEN,-1,xxs2dm
#ifdef VAR_MPI
     &                 ,lu_BVClst,
     &                  LU_scrlst,
     &                  block_list,
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

!#define LUCI_DEBUG
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
#ifdef BLABLA_DOIT
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
#endif
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


      END
!**********************************************************************

      SUBROUTINE T_ROW_TO_H(T,H,K,TKK,skip_orbital)
*
!     purpose: Set H integrals and return logical skip_orbital advising the calling
!              routine to proceed to the next orbital
*
*    Column K : H(P,K) = T(P,K)/T(K,K), P.NE.K
*    Other Columns     = 0
* - and return T_{kk} in TKK
*
*
* Jeppe Olsen, Jan 98
! Stefan Knecht, Jan 2012 revised for lucita-mcscf in dalton/dirac
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
      logical, intent(inout) :: skip_orbital
*
      CALL DZERO(H,NTOOB**2)

!     KSM  = ISMFSO(k)        
      KSM  = ISMFSO(ireots(k))
      KOFF = IBSO(KSM)
!     KREL = K - KOFF + 1
      KREL = ireots(k) - KOFF + 1
      NK   = NTOOBS(KSM)
!     symmetry offset in full 1e-matrix block
      IOFF = IFRMR(WORK(KPGINT1A(1)),1,KSM)

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(luwrt,*)        '*****************'
      write(luwrt,'(a,i3)') 'ORBITAL # == ', K
      write(luwrt,'(a,i3)') 'ORBITAL # reo', ireots(k)
      write(luwrt,*)        '*****************'
!     write(luwrt,*) 'ioff in T        ==> ',KOFF+(K-1)*NTOOB
      write(luwrt,*) 'ioff in T        ==> ',KOFF+(ireots(k)-1)*NTOOB
      write(luwrt,*) 'ioff in H        ==> ',IOFF+(KREL-1)*NK
      write(luwrt,*) 'NK elements      ==> ',NK
      write(luwrt,*) 'relative koff    ==> ',koff
      write(luwrt,*) 'relative ioff    ==> ',ioff
      write(luwrt,*) 'krel             ==> ',krel
      write(luwrt,*) 'orbital symmetry ==> ',ksm
      call flshfo(luwrt)
#undef LUCI_DEBUG
#endif
!     copy orbital block with symmetry ksm
!!    CALL DCOPY(NK,T(IOFF+(KREL-1)*NK),1,H(IOFF+(KREL-1)*NK),1)
!!!   CALL DCOPY(NK,T((K-1)*NTOOB+KOFF),1,H(IOFF+(KREL-1)*NK),1)
      CALL DCOPY(NK,T((ireots(k)-1)*NTOOB+KOFF),1,H(IOFF+(KREL-1)*NK),1)

!     set T_kk
      TKK = H(IOFF-1+(KREL-1)*NK+KREL)
      IF(TKK .NE. 0.0D0) THEN
        FAC = 1.0D0/TKK
        CALL DSCAL(NK,FAC,H(IOFF+(KREL-1)*NK),1)
        H(IOFF-1+(KREL-1)*NK+KREL) = 0.0D0
      END IF

!     check for zero row (in this case we can skip the whole loop
!     business in the calling routine) - stefan jan 2012
      skip_orbital = .false.
      zero_row = ddot(NK,H(IOFF+(KREL-1)*NK),1,H(IOFF+(KREL-1)*NK),1)
      if(zero_row == 0.0d0) skip_orbital = .true.

#ifdef LUCI_DEBUG
      WRITE(luwrt,*) ' output from T_ROW_TO_H, TKK ==> ',TKK
      WRITE(luwrt,*) ' Row to be transferred ', KREL
      WRITE(luwrt,*) ' Updated H matrix '
      WRITE(luwrt,*) ' NTOOB', NTOOB
      CALL APRBLM2(H,NTOOBS,NTOOBS,NSMOB,0)
#undef LUCI_DEBUG
#endif

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
      
!     orbital K in type ordering - stefan jan 2012 (in paris on the wy
!     to chile): not true any longer. SOT-mat entries are resorted to
!     sirius order (type order) therefore no reordering here needed.

!     KKORB = IREOST(KORB)
!     KKORB = IREOTS(KORB)
      KKORB = KORB
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(lupri,*) 'korb, kkorb (ST), (TS) ==> ',korb, kkorb,
     & IREOST(KORB), IREOTS(KORB)
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
#ifdef USE_MPI_MOD_F90
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

!#define LUCI_DEBUG
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

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
!             WRITE(lupri,*) ' KABOCC, C(IDET), T, T**2',
!    &                         KABOCC, C(IDET), T, T**2
              WRITE(lupri,*) ' KABOCC',KABOCC
              WRITE(lupri,*) ' IA IB IDET, C, T',IA,IB,IDET,C(IDET),T
#undef LUCI_DEBUG
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
