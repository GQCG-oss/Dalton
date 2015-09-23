!Simple tools common for cc routines
module cc_tools_module
   use,intrinsic :: iso_c_binding, only:c_f_pointer, c_loc, c_size_t

!`DIL backend (depends on Fortran-2003/2008, MPI-3):
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
#ifdef VAR_PTR_RESHAPE
#ifdef VAR_MPI
!#define DIL_ACTIVE
!#define DIL_DEBUG_ON
#endif
#endif
#endif

   use precision
   use memory_handling
#ifdef VAR_MPI
   use lsmpi_type
   use lsmpi_module
#endif
   use lstiming
   use typedeftype
   use integralinterfaceDEC, only: II_GET_ERI_INTEGRALBLOCK_INQUIRE
   use IntegralInterfaceMOD
   use dec_workarounds_module
   use dec_typedef_module
   use reorder_frontend_module
   use tensor_interface_module
   use gpu_interfaces
!#ifdef VAR_OPENACC
!   use openacc, only: acc_is_present
!#endif

   interface get_tpl_and_tmi
      module procedure get_tpl_and_tmi_fort, get_tpl_and_tmi_tensors
   end interface get_tpl_and_tmi
   
   abstract interface
      function ab_eq_c(a,b) result(c)
         use precision
         import
         implicit none
         real(realk), intent(in) :: a, b
         real(realk) :: c
      end function ab_eq_c
   end interface

  
   private
   public :: mo_work_dist,get_tpl_and_tmi,lspdm_get_tpl_and_tmi,&
        & get_a22_and_prepb22_terms_ex,&
        & get_a22_and_prepb22_terms_exd, combine_and_transform_sigma,&
        & get_I_cged, get_I_plusminus_le, add_int_to_sio4, get_B22_contrib_mo,&
        & print_tensor_unfolding_with_labels, &
        & squareup_block_triangular_squarematrix,&
        & calc_i_leq_j, calc_i_geq_j, simulate_intloop_and_get_worksize, &
        & successive_4ao_mo_trafo, solver_energy_full, solver_decnp_full,&
        & ccsdpt_energy_e4_full, ccsdpt_energy_e5_full, ccsdpt_decnp_e4_full,&
        & ccsdpt_decnp_e5_full, get_t1_matrices, get_nbuffs_scheme_0, &
        & get_split_scheme_0, get_tpl_and_tmi_fort, solver_energy_full_mod

   
   contains

   subroutine mo_work_dist(m,fai,tl,have_part,nod)
      implicit none
      integer,intent(in) :: m
      integer,intent(inout)::fai
      integer,intent(inout)::tl
      logical,intent(out)  :: have_part
      integer(kind=ls_mpik),optional,intent(in)::nod
      integer(kind=ls_mpik) :: nnod, me
      integer :: l,ml

      me   = 0
      nnod = 1
#ifdef VAR_MPI
      nnod = infpar%lg_nodtot
      me   = infpar%lg_mynum
#endif

      if(present(nod))me=nod

      !Setting transformation variables for each rank
      !**********************************************
      l   = (m) / nnod
      ml  = mod(m,nnod)
      fai = me * l + 1
      tl  = l

      if(ml>0)then
         if(me<ml)then
            fai = fai + me
            tl  = l + 1
         else
            fai = fai + ml
            tl  = l
         endif
      endif

      !If too many nodes are used set the first element to an invalid counter
      !**********************************************************************
      if(fai > m)then
         fai = -1
         have_part = .false.
      else
         have_part = .true.
      endif

   end subroutine mo_work_dist



   !> \brief Reorder t to use symmetry in both occupied and virtual indices,
   !> thereby restricting the first virtual to be less equal to the second and
   !> make symmetric and antisymmetric combinations of these 
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine get_tpl_and_tmi_tensors(t2,tpl,tmi)
      implicit none
      !>input the full doubles ampitudes in the order nv nv no no
      type(tensor),intent(inout) :: t2
      !> output amplitudes with indices restricted c<d,i<j
      type(tensor), intent(inout) :: tpl
      !> output amplitudes with indices restricted c>d,i<j
      type(tensor), intent(inout) :: tmi
      integer ::crvd,no,nv
      nv   = t2%dims(1)
      no   = t2%dims(3)

      if ( (t2%itype  == TT_DENSE .or. t2%itype  == TT_REPLICATED) .and. &
          &(tpl%itype == TT_DENSE .or. tpl%itype == TT_REPLICATED) .and. &
          &(tmi%itype == TT_DENSE .or. tmi%itype == TT_REPLICATED) ) then

         call get_tpl_and_tmi_fort(t2%elm1,nv,no,tpl%elm1,tmi%elm1)

         !maybe a sync replicated is needed here
 
      elseif( t2%itype == TT_TILED_DIST .and. tpl%itype == TT_TILED_DIST .and. tmi%itype == TT_TILED_DIST)then

         call lspdm_get_tpl_and_tmi(t2,tpl,tmi)

      else
         call lsquit("ERROR(get_tpl_and_tmi_tensors): combination of itypes not implemented",-1)
      endif

   end subroutine get_tpl_and_tmi_tensors

   subroutine lspdm_get_tpl_and_tmi(t2,tpl,tmi)
      implicit none
      type(tensor),intent(inout) :: t2 !intent(in)
      type(tensor),intent(inout) :: tpl, tmi !tpl[nor,nvr],tmi[nor,nvr]

      integer :: i,j,a,b,nocc,nvirt,da,db,di,dj,gtnr,lt,nelt
      integer :: otmi(2),otpl(2),ot2(4)
#ifdef VAR_PTR_RESHAPE
      real(realk), pointer, contiguous :: tt1(:,:,:,:),tt2(:,:,:,:),tpm(:,:)
      real(realk), pointer, contiguous :: buf1(:),buf2(:)
#else
      real(realk), pointer :: tt1(:,:,:,:),tt2(:,:,:,:),tpm(:,:)
      real(realk), pointer :: buf1(:),buf2(:)
#endif
      integer :: dcged,dilej,ccged,cilej,gcged,gilej
      real(realk) :: sol
      integer :: nor,no,nvr,nv,k,c,d
      integer :: maxtile(4), mintile(4),tdim(4)
      integer :: tctr1,tctr2,tctr3,tctr4
      integer :: nelms,nvrs,nors
      integer :: gc,gd,gi,gj
      integer :: mingilej,maxgilej,mingcged,maxgcged
      integer :: check, nfound, founds(t2%ntiles), mtile(t2%mode), ctile
      logical :: Already_in,contributed

#ifdef VAR_MPI
      if( t2%access_type /= AT_ALL_ACCESS .or. tpl%access_type /= AT_ALL_ACCESS .or. tmi%access_type /= AT_ALL_ACCESS  )then
         call lsquit("ERROR(lspdm_get_tpl_and_tmi): access types of the arrays not (yet) possible",-1)
      endif

      nv = t2%dims(1)
      no = t2%dims(3)

      nvr = nv * (nv + 1) / 2
      nor = no * (no + 1) / 2

      nvrs = tpl%tdim(2)
      nors = tpl%tdim(1)

      call mem_alloc(buf1,t2%tsize)
      call mem_alloc(buf2,t2%tsize)

      do lt = 1, tpl%nlti

         gtnr = tpl%ti(lt)%gt

         call get_midx(gtnr,otpl,tpl%ntpm,tpl%mode)

         !Facilitate access
!        call c_f_pointer(c_loc(tpl%ti(lt)%t),tpm,shape=tpl%ti(lt)%d(1:2)) !`DIL
#ifdef VAR_PTR_RESHAPE
         tpm(1:tpl%ti(lt)%d(1),1:tpl%ti(lt)%d(2)) => tpl%ti(lt)%t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
         call c_f_pointer(c_loc(tpl%ti(lt)%t(1)),tpm,tpl%ti(lt)%d)
#else
         call lsquit('ERROR(lspdm_get_tpl_and_tmi): unable to reshape pointers!',-1)
#endif

         !build list of tiles to get for the current tpl tile
         !get offset for tile counting
         do j=1,tpl%mode
            otpl(j)=(otpl(j)-1)*tpl%tdim(j)
         enddo

         dcged = tpl%ti(lt)%d(2)
         dilej = tpl%ti(lt)%d(1)

         nfound = 0
         do cilej=1,dilej

            gilej = otpl(1) + cilej

            call calc_i_leq_j(gilej,no,i,j)

            do ccged=1,dcged

               gcged = otpl(2) + ccged

               call calc_i_geq_j(gcged,nv,c,d)

               mtile = [ (c-1)/t2%tdim(1)+1 , (d-1)/t2%tdim(2)+1, (j-1)/t2%tdim(3)+1, (i-1)/t2%tdim(4)+1 ]

               ctile = get_cidx(mtile,t2%ntpm,t2%mode)

               Already_in = .false.

               do check = 1, nfound
                  if( founds(check) == ctile ) then
                     Already_in = .true.
                     exit
                  endif
               enddo

               if(.not.Already_in)then
                  founds(nfound+1) = ctile
                  nfound = nfound + 1
               endif

            enddo
         enddo

         mingilej = otpl(1) + 1
         maxgilej = otpl(1) + dilej
         mingcged = otpl(2) + 1
         maxgcged = otpl(2) + dcged

         do check = 1,nfound

            ctile = founds(check)

            call get_midx(ctile,mtile,t2%ntpm,t2%mode)

            call get_tile_dim(tdim,t2,mtile)

            nelms = 1
            do i=1,t2%mode
               nelms = nelms * tdim(i)
            enddo

            call tensor_get_tile(t2,mtile,buf1,nelms)

            if(mtile(2)/=mtile(1)) call tensor_get_tile(t2,[mtile(2),mtile(1),mtile(3),mtile(4)],buf2,nelms)

!           call c_f_pointer(c_loc(buf1),tt1,shape=tdim(1:4)) !`DIL
#ifdef VAR_PTR_RESHAPE
            tt1(1:tdim(1),1:tdim(2),1:tdim(3),1:tdim(4)) => buf1
            if(mtile(2)==mtile(1))then
               tt2(1:tdim(2),1:tdim(1),1:tdim(3),1:tdim(4)) => buf1
            else
               tt2(1:tdim(2),1:tdim(1),1:tdim(3),1:tdim(4)) => buf2
            endif
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
            call c_f_pointer(c_loc(buf1(1)),tt1,tdim)
            if(mtile(2)==mtile(1))then
!             call c_f_pointer(c_loc(buf1),tt2,shape=(/tdim(2),tdim(1),tdim(3),tdim(4)/)) !`DIL
              call c_f_pointer( c_loc(buf1(1)), tt2, [tdim(2),tdim(1),tdim(3),tdim(4)] )
            else
!             call c_f_pointer(c_loc(buf2),tt2,shape=(/tdim(2),tdim(1),tdim(3),tdim(4)/)) !`DIL
              call c_f_pointer( c_loc(buf2(1)), tt2, [tdim(2),tdim(1),tdim(3),tdim(4)] )
            endif
#else
            call lsquit('ERROR(lspdm_get_tpl_and_tmi): unable to reshape pointers!',-1)
#endif

            !get offset for tile counting
            ot2(1)=(mtile(1)-1)*t2%tdim(1)
            ot2(2)=(mtile(2)-1)*t2%tdim(2)
            ot2(3)=(mtile(3)-1)*t2%tdim(3)
            ot2(4)=(mtile(4)-1)*t2%tdim(4)

            da = tdim(1)
            db = tdim(2)
            di = tdim(3)
            dj = tdim(4)

            contributed  = .false.

            do j=1,dj
               gj = ot2(4)+j
               do i=1,di
                  gi = ot2(3)+i

                  if( gi <= gj )then
                     gilej = gi + ((gj-1)*gj)/2

                     if( mingilej <= gilej .and. gilej<= maxgilej )then

                        do b=1,db
                           gd = ot2(2)+b
                           do a=1,da
                              gc = ot2(1)+a

                              gcged = gc + (gd - 1) * nv - ((gd-1)*gd)/2

                              if( gc > gd .and. mingcged <= gcged .and. gcged<= maxgcged )then
                                 !print *,infpar%lg_mynum,"ONE"

                                 ccged = mod(gcged-1,nvrs)+1
                                 cilej = mod(gilej-1,nors)+1

                                 tpm(cilej,ccged) = tt1(a,b,i,j)+tt2(b,a,i,j)

                                 contributed  = .true.

                              else if( gc == gd .and. mingcged <= gcged .and. gcged<= maxgcged)then
                                 !print *,infpar%lg_mynum,"TWO"

                                 ccged = mod(gcged-1,nvrs)+1
                                 cilej = mod(gilej-1,nors)+1

                                 tpm(cilej,ccged) = 0.5E0_realk*(tt1(a,b,i,j)+tt2(b,a,i,j))

                                 contributed  = .true.

                              endif
                           end do
                        end do

                     endif

                  endif

               end do
            end do

         enddo

      enddo



      do lt = 1, tmi%nlti

         gtnr = tmi%ti(lt)%gt

         call get_midx(gtnr,otmi,tmi%ntpm,tmi%mode)

         !Facilitate access
!        call c_f_pointer(c_loc(tmi%ti(lt)%t),tpm,shape=tmi%ti(lt)%d(1:2)) !`DIL
#ifdef VAR_PTR_RESHAPE
         tpm(1:tmi%ti(lt)%d(1),1:tmi%ti(lt)%d(2)) => tmi%ti(lt)%t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
         call c_f_pointer( c_loc(tmi%ti(lt)%t(1)), tpm, tmi%ti(lt)%d )
#else
         call lsquit('ERROR(lspdm_get_tpl_and_tmi): unable to reshape pointers!',-1)
#endif

         !build list of tiles to get for the current tmi tile
         !get offset for tile counting
         do j=1,tmi%mode
            otmi(j)=(otmi(j)-1)*tmi%tdim(j)
         enddo

         dcged = tmi%ti(lt)%d(2)
         dilej = tmi%ti(lt)%d(1)

         nfound = 0
         do cilej=1,dilej

            gilej = otmi(1) + cilej

            call calc_i_leq_j(gilej,no,i,j)

            do ccged=1,dcged

               gcged = otmi(2) + ccged

               call calc_i_geq_j(gcged,nv,c,d)

               mtile = [ (c-1)/t2%tdim(1)+1 , (d-1)/t2%tdim(2)+1, (j-1)/t2%tdim(3)+1, (i-1)/t2%tdim(4)+1 ]

               ctile = get_cidx(mtile,t2%ntpm,t2%mode)

               Already_in = .false.

               do check = 1, nfound
                  if( founds(check) == ctile ) then
                     Already_in = .true.
                     exit
                  endif
               enddo

               if(.not.Already_in)then
                  founds(nfound+1) = ctile
                  nfound = nfound + 1
               endif

            enddo
         enddo

         mingilej = otmi(1) + 1
         maxgilej = otmi(1) + dilej
         mingcged = otmi(2) + 1
         maxgcged = otmi(2) + dcged

         do check = 1,nfound

            ctile = founds(check)

            call get_midx(ctile,mtile,t2%ntpm,t2%mode)

            call get_tile_dim(tdim,t2,mtile)

            nelms = 1
            do i=1,t2%mode
               nelms = nelms * tdim(i)
            enddo

            call tensor_get_tile(t2,mtile,buf1,nelms)

            if(mtile(2)/=mtile(1)) call tensor_get_tile(t2,[mtile(2),mtile(1),mtile(3),mtile(4)],buf2,nelms)

!           call c_f_pointer(c_loc(buf1),tt1,shape=tdim(1:4)) !`DIL
#ifdef VAR_PTR_RESHAPE
            tt1(1:tdim(1),1:tdim(2),1:tdim(3),1:tdim(4)) => buf1
            if(mtile(2)==mtile(1))then
               tt2(1:tdim(2),1:tdim(1),1:tdim(3),1:tdim(4)) => buf1
            else
               tt2(1:tdim(2),1:tdim(1),1:tdim(3),1:tdim(4)) => buf2
            endif
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
            call c_f_pointer( c_loc(buf1(1)), tt1, tdim )
            if(mtile(2)==mtile(1))then
!             call c_f_pointer(c_loc(buf1),tt2,shape=(/tdim(2),tdim(1),tdim(3),tdim(4)/)) !`DIL
              call c_f_pointer( c_loc( buf1(1) ) , tt2, [tdim(2),tdim(1),tdim(3),tdim(4)] )
            else
!             call c_f_pointer(c_loc(buf2),tt2,shape=(/tdim(2),tdim(1),tdim(3),tdim(4)/)) !`DIL
              call c_f_pointer( c_loc( buf2(1) ) , tt2, [tdim(2),tdim(1),tdim(3),tdim(4)] )
            endif
#else
            call lsquit('ERROR(lspdm_get_tpl_and_tmi): unable to reshape pointers!',-1)
#endif

            !get offset for tile counting
            ot2(1)=(mtile(1)-1)*t2%tdim(1)
            ot2(2)=(mtile(2)-1)*t2%tdim(2)
            ot2(3)=(mtile(3)-1)*t2%tdim(3)
            ot2(4)=(mtile(4)-1)*t2%tdim(4)

            da = tdim(1)
            db = tdim(2)
            di = tdim(3)
            dj = tdim(4)

            do j=1,dj
               gj = ot2(4)+j
               do i=1,di
                  gi = ot2(3)+i

                  if( gi <= gj )then
                     gilej = gi + ((gj-1)*gj)/2

                     if( mingilej <= gilej .and. gilej<= maxgilej )then

                        do b=1,db
                           gd = ot2(2)+b
                           do a=1,da
                              gc = ot2(1)+a

                              gcged = gc + (gd - 1) * nv - ((gd-1)*gd)/2

                              if( gc >= gd .and. mingcged <= gcged .and. gcged<= maxgcged )then

                                 ccged = mod(gcged-1,nvrs)+1
                                 cilej = mod(gilej-1,nors)+1

                                 tpm(cilej,ccged) = tt1(a,b,i,j) - tt2(b,a,i,j)

                              endif
                           end do
                        end do

                     endif

                  endif

               end do
            end do

         enddo

      enddo

      call mem_dealloc(buf1)
      call mem_dealloc(buf2)

      call lsmpi_barrier(infpar%lg_comm)
#endif
   end subroutine lspdm_get_tpl_and_tmi

   !> \brief Reorder t to use symmetry in both occupied and virtual indices,
   !> thereby restricting the first virtual to be less equal to the second and
   !> make symmetric and antisymmetric combinations of these 
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine get_tpl_and_tmi_fort(t2,nv,no,tpl,tmi)
      implicit none
      !> number of virtual and occupied orbitals
      integer, intent(in) :: nv, no
      !>input the full doubles ampitudes in the order nv nv no no
      real(realk),intent(in) :: t2(nv,nv,no,no)
      !> output amplitudes with indices restricted c>d,i<j
      real(realk), intent(inout) :: tpl(:)
      !> output amplitudes with indices restricted c>d,i<j
      real(realk), intent(inout) :: tmi(:)
      integer :: dims(4)
      integer ::i,j,d, ilej,cged,crvd
      integer :: counter,counter2
      dims(1)=nv
      dims(2)=nv
      dims(3)=no
      dims(4)=no

      !get the combined reduced virtual dimension
      crvd = dims(1) * (dims(1)+1) / 2
      !$OMP PARALLEL PRIVATE(i,j,d,cged,ilej) SHARED(crvd,t2,tmi,tpl,dims) DEFAULT(NONE)
      !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
      do j=1,dims(3)
         do i=1,dims(3)

            if(i>j)cycle

            cged=1
            ilej = i + (j-1) * j / 2

            do d=1,dims(2)
               !copy the elements c>d into tpl and tmi
               call dcopy(dims(2)-d+1,t2(d,d,i,j),1,tpl(cged+(ilej-1)*crvd),1)
               call dcopy(dims(2)-d+1,t2(d,d,i,j),1,tmi(cged+(ilej-1)*crvd),1)

               !add and subtract the counterparts for tpl and tmi respectively
               call daxpy(dims(1)-d+1,-1.0E0_realk,t2(d,d,i,j),dims(1),tmi(cged+(ilej-1)*crvd),1)
               call daxpy(dims(1)-d+1,1.0E0_realk,t2(d,d,i,j),dims(1),tpl(cged+(ilej-1)*crvd),1)

               !take care of plus combination diagonal elements to be only half
               tpl(cged+(ilej-1)*crvd)=0.5E0_realk*tpl(cged+(ilej-1)*crvd)

               !count in the triangular part of the matrix
               cged = cged+dims(2)-d+1
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
   end subroutine get_tpl_and_tmi_fort

   !> \brief calculate a and b terms in a kobayashi fashion
   !> \author Patrick Ettenhuber
   !> \date December 2012
   subroutine get_a22_and_prepb22_terms_ex(w0,w1,w2,w3,tpl,tmi,no,nv,nb,fa,fg,la,lg,&
         &xo,yo,xv,yv,om2,sio4,s,wszes,lo,twork,tcomm,order,rest_occ_om2,scal)
      implicit none
      !> W0 SIZE
      integer(kind=8),intent(in)  :: wszes(4)
      !> workspace with exchange integrals
      real(realk),intent(inout) :: w1(wszes(2))
      !> empty workspace of correct sizes
      real(realk),intent(inout),target :: w0(wszes(1)),w2(wszes(3))
      real(realk),intent(inout) :: w3(wszes(4))
      !> the t+ and t- combinations with a value of the amplitudes with the
      !diagonal elements divided by two
      type(tensor),intent(inout) :: tpl,tmi
      !> number of occupied, virutal and ao indices
      integer, intent(in) :: no,nv,nb
      !> first alpha and first gamma indices of the current loop
      integer, intent(in) :: fa,fg
      !> lengths of the alpha ang gamma batches in the currnet loop
      integer, intent(in) :: la,lg
      !> the doubles residual to update
      !real(realk), intent(inout) :: om2(:)
      type(tensor), intent(inout) :: om2
      !> the lambda transformation matrices
      real(realk), intent(in)    :: xo(:),yo(:),xv(:),yv(:)
      !> the sio4 matrix to calculate the b2.2 contribution
      type(tensor),intent(inout) ::sio4
      !> scheme
      integer,intent(in) :: s
      logical,intent(in) :: lo
      !timing information
      real(realk) :: twork,tcomm
      integer,optional,intent(in) :: order(4)
      logical,optional,intent(in) :: rest_occ_om2
      real(realk),optional :: scal
      integer :: goffs,aoffs,tlen,tred,nor,nvr
      integer(kind=8) :: s0, s2, s3
      character(80) :: msg
      real(realk), pointer :: buf(:)
      integer(kind=8) :: nbuf
      integer :: faleg,laleg,laleg_req,i,nerrors
      integer, parameter :: nids = 15
      integer :: nids_use
#ifdef VAR_OPENACC
      integer(kind=acc_handle_kind) :: acc_h(nids),transp
      integer,parameter             :: lsacc_sync    = acc_async_sync
      integer,parameter             :: lsacc_handlek = acc_handle_kind
      integer(c_size_t)             :: total_gpu,free_gpu
#else
      integer                       :: acc_h(nids),transp
      integer,parameter             :: lsacc_sync    = -4
      integer,parameter             :: lsacc_handlek = 4
#endif
      type(c_ptr)                   :: cub_h(nids),dummy47
      real(realk) :: p10, nul
      integer(4)  :: stat,curr_id
      logical :: one

      acc_h = 0
      cub_h = c_null_ptr

      p10 = 1.0E0_realk
      nul = 0.0E0_realk

      stat = 0

      call time_start_phase(PHASE_WORK)

      if (DECinfo%acc_sync) then
         acc_h  = lsacc_sync
         transp = lsacc_sync
      else
         do curr_id = 1,nids
            acc_h(curr_id) = int(curr_id-1,kind=lsacc_handlek)
         enddo
         transp = nids
      endif

#ifdef VAR_CUBLAS
      do curr_id = 1,nids
         ! initialize the CUBLAS context
         stat = cublasCreate_v2(cub_h(curr_id))
         if(stat/=0)print*,"warning: cublas create failed 1",stat
         dummy47 = acc_get_cuda_stream(acc_h(curr_id))
         stat = cublasSetStream_v2(cub_h(curr_id), dummy47)
         if(stat/=0)print*,"warning: cublas set stream failed 1",stat
      enddo
#endif

      s0 = wszes(1)
      s2 = wszes(3)
      s3 = wszes(4)

      nor=(no*(no+1))/2
      nvr=(nv*(nv+1))/2

      !Determine the offsets in the alpha and gamma indices, which arise from
      !non uniform batch distributions, e.g. consider:
      !Case 1:   gamma _____
      !               |\    |
      !               |_\___|
      !  alpha        |  \  |
      !               |___\_| here an offset in gamma for the second gamma batch is
      !               needed, since the elements to be considered do not
      !               begin with the first element in that batch
      goffs=0
      if(fa-fg>0)goffs=fa-fg
      !Case 2:   gamma ______
      !               |\ |   |
      !               | \|   |
      !  alpha        |  |\  |
      !               |__|_\_| here an offset in alpha for the second alpha batch is
      !               needed, since the elements to be considered do not
      !               begin with the first element in that batch
      aoffs=0
      if(fg-fa>0)aoffs=fg-fa

      !Determine the dimension of the triangular part in the current batch
      !and the total number of elements (tred) in the triangular plus
      !rectangular parts
      tred=0
      tlen=min(min(min(la,lg),fg+lg-fa),fa+la-fg)

      !calculate amount of triangular elements in the current batch
      if(fa+la-1>=fg)tred= tlen*(tlen+1)/2
      !add the rectangular contribution if lg is larger than tlen 
      if(fa>=fg.and.fg+lg-fa-tlen>0)tred=tred+(fg+lg-fa-tlen)*la
      if(fa<fg)tred=tred+lg*aoffs
      if(fa<fg.and.fa+la>fg) tred=tred+(lg-tlen)*(la-aoffs)
      !if only rectangular
      if(tlen<=0)then
         tlen=0
         aoffs=0
         goffs=0
         tred=la*lg
      endif

      !Request a fraction of tred resulting in intermediates that can be stored
      !in GPU memory if a GPU is available, if not just use the full tred
#ifdef VAR_OPENACC
      total_gpu = 0
      free_gpu  = 0
      call get_dev_mem(total_gpu,free_gpu)
      select case(s)
      case(4,3)

         !subtract the memory required to store yv and tpl/tmi
         free_gpu = free_gpu - 8 * (nb*nv + nor*nvr)

         !decrease the number of async ids with minimal batch size until all fits into GPU memory
         laleg_req = 1
         do nids_use = nids, 0
            !              #of async ids       w0 size          w2 size            w3 size
            if((free_gpu - nids_use * 8 * (nb*laleg_req*nv + nb*laleg_req*nb + laleg_req*nor*nvr) > 0) exit
         enddo

         if(nids_use <= 0)call lsquit("ERROR(get_a22_and_prepb22_terms_ex): not enough memory on the GPU available",-1)

         !increase the batch size until the intermediates to not fit into the GPU memory anymore
         do laleg_req = 2, tred + 1
            !            #of async ids      w0 size          w2 size            w3 size
            if((free_gpu - nids_use * 8 * (nb*laleg_req*nv + nb*laleg_req*nb + laleg_req*nor*nvr) < 0 ) exit
         enddo

      !reduce the found requested size by one which was the last size that fulfilled the memory requirements
         laleg_req = laleg_req - 1

      case default
         !default use the full batch
         laleg_req = tred
      end select
#else
      nids_use  = 1
      laleg_req = tred
#endif


      select case(s)
      case(4,3)


         !$acc enter data copyin(yv(1:nb*nv),tpl%elm1(1:nor*nvr)) create(w0(1:nb*laleg_req*nv)) async(transp)

         !!SYMMETRIC COMBINATION
         ! (w2): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
         curr_id = 0
         call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
         do faleg=1,tred,laleg_req
            

            laleg = laleg_req
            if(tred-faleg+1<laleg_req) laleg = tred-faleg+1

            curr_id = mod(curr_id,nids_use)+1

            !(w2):I+ [beta delta alpha<=gamma] <= (w2):I [beta delta alpha gamma ] + (w2):I[delta beta alpha gamma]
            call get_I_plusminus_le(w2,'+',fa,fg,la,lg,nb,tlen,tred,goffs,s2,faleg,laleg)

            !$acc data copyin(w2(1:nb*laleg*nb)) copyout(w3(1+(faleg-1)*nor:nor+(faleg+laleg-2)*nor)) &
            !$acc& wait(transp) async(acc_h(curr_id)) 

            !(w0):I+ [delta alpha<=gamma c] = (w2):I+ [beta, delta alpha<=gamma] * Lambda^h[beta c]
            call ls_dgemm_acc('t','n',nb*laleg,nv,nb,p10,w2,nb,yv,nb,nul,w0,nb*laleg,&
            &i8*nb*laleg*nb,i8*nv*nb,i8*nb*laleg*nv,acc_h(curr_id),cub_h(curr_id))

            !(w2):I+ [alpha<=gamma c d] = (w0):I+ [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
            call ls_dgemm_acc('t','n',laleg*nv,nv,nb,p10,w0,nb,yv,nb,nul,w2,nv*laleg,&
            &i8*laleg*nb*nv,i8*nb*nv,i8*laleg*nv*nv,acc_h(curr_id),cub_h(curr_id))

            !(w0):I+ [alpha<=gamma c>=d] <= (w2):I+ [alpha<=gamma c d] 
            call get_I_cged(w0,w2,laleg,nv,acc_h=acc_h(curr_id))

#ifdef VAR_OPENACC
            !(w3.1):sigma+ [i>= j alpha<=gamma] = t+ [c>=d i>=j]^T * (w2):I+ [alpha<=gamma c>=d]^T
            call ls_dgemm_acc('t','t',nor,laleg,nvr,0.5E0_realk,tpl%elm1,nvr,w0,laleg,nul,w3(1+(faleg-1)*nor),nor,&
            &i8*laleg*nvr,i8*nor*nvr,i8*laleg*nor,acc_h(curr_id),cub_h(curr_id))
#else
            !(w3.1):sigma+ [alpha<=gamma i>=j] = (w2):I+ [alpha<=gamma c>=d] * t+ [c>=d i>=j]
            call dgemm('n','n',laleg,nor,nvr,0.5E0_realk,w0,laleg,tpl%elm1,nvr,nul,w3(faleg),tred)
#endif
            !$acc end data 
         enddo

         !$acc wait async(transp)
         !$acc exit data delete(tpl%elm1(1:nor*nvr)) async(transp)
         !$acc enter data copyin(tmi%elm1(1:nor*nvr)) async(transp)

#ifdef VAR_OPENACC
         call array_reorder_2d(p10,w3,nor,tred,[2,1],nul,w2)
         call array_reorder_2d(p10,w2,tred,nor,[1,2],nul,w3)
#endif

         !!ANTI-SYMMETRIC COMBINATION
         ! (w2): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
         curr_id = 0
         call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
         do faleg=1,tred,laleg_req

            curr_id = mod(curr_id,nids_use)+1

            laleg = laleg_req
            if(tred-faleg+1<laleg_req) laleg = tred-faleg+1

            !(w2):I+ [beta delta alpha<=gamma] <= (w2):I [beta delta alpha gamma ] + (w2):I[delta beta alpha gamma]
            call get_I_plusminus_le(w2,'-',fa,fg,la,lg,nb,tlen,tred,goffs,s2,faleg,laleg)

            !$acc data copyin(w2(1:nb*laleg*nb)) copyout(w3(tred*nor+1+(faleg-1)*nor:tred*nor+nor+(faleg+laleg-2)*nor))&
            !$acc& wait(transp) async(acc_h(curr_id))

            !(w0):I+ [delta alpha<=gamma c] = (w2):I+ [beta, delta alpha<=gamma] * Lambda^h[beta c]
            call ls_dgemm_acc('t','n',nb*laleg,nv,nb,p10,w2,nb,yv,nb,nul,w0,nb*laleg,&
            &i8*nb*laleg*nb,i8*nv*nb,i8*nb*laleg*nv,acc_h(curr_id),cub_h(curr_id))

            !(w2):I+ [alpha<=gamma c d] = (w0):I+ [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
            call ls_dgemm_acc('t','n',laleg*nv,nv,nb,p10,w0,nb,yv,nb,nul,w2,nv*laleg,&
            &i8*laleg*nb*nv,i8*nb*nv,i8*laleg*nv*nv,acc_h(curr_id),cub_h(curr_id))

            !(w0):I+ [alpha<=gamma c>=d] <= (w2):I+ [alpha<=gamma c d] 
            call get_I_cged(w0,w2,laleg,nv,acc_h=acc_h(curr_id))

            !(w3.1):sigma+ [alpha<=gamma i>=j] = (w2):I+ [alpha<=gamma c>=d] * t+ [c>=d i>=j]
#ifdef VAR_OPENACC
            call ls_dgemm_acc('t','t',nor,laleg,nvr,0.5E0_realk,tmi%elm1,nvr,w0,laleg,nul,w3(tred*nor+1+(faleg-1)*nor),nor,&
            &i8*laleg*nvr,i8*nor*nvr,i8*laleg*nor,acc_h(curr_id),cub_h(curr_id))
#else
            call dgemm('n','n',laleg,nor,nvr,0.5E0_realk,w0,laleg,tmi%elm1,nvr,nul,w3(tred*nor+faleg),tred)
#endif
            !$acc end data 

         enddo
         !$acc wait
         !$acc exit data delete(yv(1:nb*nv),tmi%elm1(1:nor*nvr),w0(1:nb*laleg_req*nv))

#ifdef VAR_OPENACC
         call array_reorder_2d(p10,w3(tred*nor+1:tred*nor+tred*nor),nor,tred,[2,1],nul,w2)
         call array_reorder_2d(p10,w2,tred,nor,[1,2],nul,w3(tred*nor+1:tred*nor+tred*nor))

#ifdef VAR_CUBLAS
         ! Destroy the CUBLAS context
         do curr_id = 1,nids
            stat = cublasDestroy_v2(cub_h(curr_id))
            if(stat/=0)print*,"warning: cublas destroy failed 1",stat
         enddo
#endif
#endif

      case(2)

         if( s0 - tred*nvr > s2 - nor*nvr)then
            buf  => w0(tred*nvr+1:)
            nbuf = s0 - tred*nvr
         else
            buf  => w2(nor*nvr+1:)
            nbuf = s2 - nor*nvr
         endif

         !!SYMMETRIC COMBINATION
         ! (w2): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
         call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
         !(w2):I+ [delta alpha<=gamma beta] <= (w2):I [beta delta alpha gamma ] + (w2):I[beta delta alpha gamma]
         call get_I_plusminus_le(w2,'+',fa,fg,la,lg,nb,tlen,tred,goffs,s2,1,tred)
         !(w0):I+ [delta alpha<=gamma c] = (w2):I+ [beta, delta alpha<=gamma] * Lambda^h[beta c]
         call dgemm('t','n',nb*tred,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nb*tred)
         !(w0):I+ [alpha<=gamma c d] = (w2):I+ [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
         call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w0,nb,yv,nb,0.0E0_realk,w2,nv*tred)
         !(w0):I+ [alpha<=gamma c>=d] <= (w2):I+ [alpha<=gamma c d] 
         call get_I_cged(w0,w2,tred,nv)
         !(w2.2)tpl
         call tensor_gather(1.0E0_realk,tpl,0.0E0_realk,w2,(i8*nor)*nvr,oo=[2,1],wrk=buf,iwrk=nbuf)
         !(w3.1):sigma+ [alpha<=gamma i>=j] = (w2):I+ [alpha<=gamma c>=d] * (w2):t+ [c>=d i>=j]
         call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w0,tred,w2,nvr,0.0E0_realk,w3,tred)

         !!ANTI-SYMMETRIC COMBINATION
         ! (w2): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
         call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
         !(w2):I- [delta alpha<=gamma beta] <= (w2):I [beta delta alpha gamma ] - (w2):I[beta delta alpha gamma]
         call get_I_plusminus_le(w2,'-',fa,fg,la,lg,nb,tlen,tred,goffs,s2,1,tred)
         !(w0):I- [delta alpha<=gamma c] = (w2):I- [delta alpha<=gamma beta] * Lambda^h[beta c]
         call dgemm('t','n',nb*tred,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nb*tred)
         !(w2):I- [alpha<=gamma c d] = (w0):I- [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
         call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w0,nb,yv,nb,0.0E0_realk,w2,nv*tred)
         !(w0):I- [alpha<=gamma c<=d] <= (w2):I- [alpha<=gamma c d]
         call get_I_cged(w0,w2,tred,nv)
         !(w2)tmi
         call tensor_gather(1.0E0_realk,tmi,0.0E0_realk,w2,(i8*nor)*nvr,oo=[2,1],wrk=buf,iwrk=nbuf)
         !(w3.2):sigma- [alpha<=gamma i<=j] = (w2):I- [alpha<=gamma c>=d] * (w2):t- [c>=d i>=j]
         call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w0,tred,w2,nvr,0.0E0_realk,w3(tred*nor+1),tred)
      case default
         call lsquit("ERROR(get_a22_and_prepb22_terms_ex): wrong scheme on input",-1)
      end select

      !COMBINE THE TWO SIGMAS OF W3 IN W2
      !(w2):sigma[alpha<=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] + 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
      !(w2):sigma[alpha>=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] - 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
      call combine_and_transform_sigma(om2,w0,w2,w3,s0,s2,s3,xv,xo,sio4,nor,tlen,tred,fa,fg,la,lg,&
         &no,nv,nb,goffs,aoffs,s,lo,twork,tcomm,order=order,rest_occ_om2=rest_occ_om2,scal=scal,sio4_ilej = (s/=2))

      call time_start_phase(PHASE_WORK, at=twork)
   end subroutine get_a22_and_prepb22_terms_ex

#ifdef DIL_ACTIVE
   !> \brief calculate a and b terms in a Kobayashi fashion
   !> \author Patrick Ettenhuber, Dmitry I. Lyakh (`DIL)
   !> \date December 2012
   subroutine get_a22_and_prepb22_terms_exd(w0,w1,w2,w3,tpl,tmi,no,nv,nb,fa,fg,la,lg,&
         &xo,yo,xv,yv,om2,sio4,s,wszes,o2ilej,lo,dil_lock_out,twork,tcomm,order,rest_occ_om2,scal)
      implicit none
      !> workspace with exchange integrals
      real(realk),intent(inout) :: w1(:)
      !> empty workspace of correct sizes
      real(realk),intent(inout) :: w0(:),w2(:),w3(:)
      !> the t+ and t- combinations with a value of the amplitudes with the
      !diagonal elements divided by two
      type(tensor),intent(inout) :: tpl,tmi,o2ilej !tpl[nor,nvr],tmi[nor,nvr],o2ilej[nv,nv,nor]
      !> number of occupied, virutal and ao indices
      integer, intent(in) :: no,nv,nb
      !> first alpha and first gamma indices of the current loop
      integer, intent(in) :: fa,fg
      !> lengths of the alpha ang gamma batches in the currnet loop
      integer, intent(in) :: la,lg
      !> the doubles residual to update
      !real(realk), intent(inout) :: om2(:)
      type(tensor), intent(inout) :: om2
      !> the lambda transformation matrices
      real(realk), intent(in)    :: xo(:),yo(:),xv(:),yv(:)
      !> the sio4 matrix to calculate the b2.2 contribution
      type(tensor),intent(inout) :: sio4
      !> scheme
      integer,intent(in) :: s
      logical,intent(in) :: lo,dil_lock_out
      !timing information
      real(realk) :: twork,tcomm
      !> W0 SIZE
      integer(kind=8),intent(in)  :: wszes(4)
      integer,optional,intent(in) :: order(4)
      logical,optional,intent(in) :: rest_occ_om2
      real(realk),optional :: scal
      integer :: goffs,aoffs,tlen,tred,nor,nvr
      integer(kind=8) :: s0, s2, s3
!{`DIL:
#ifdef DIL_ACTIVE
     integer:: nors,nvrs,scheme
     character(256):: tcs
     type(dil_tens_contr_t):: tch
     integer(INTL):: dil_mem,l0
     integer(INTD):: i0,i1,i2,i3,errc,tens_rank,tens_dims(MAX_TENSOR_RANK),tens_bases(MAX_TENSOR_RANK)
     integer(INTD):: ddims(MAX_TENSOR_RANK),ldims(MAX_TENSOR_RANK),rdims(MAX_TENSOR_RANK)
     integer(INTD):: dbase(MAX_TENSOR_RANK),lbase(MAX_TENSOR_RANK),rbase(MAX_TENSOR_RANK)
     real(realk):: r0
     integer(INTD):: sch_sym=1
#endif
!}

      call time_start_phase(PHASE_WORK)

      s0 = wszes(1)
      s2 = wszes(3)
      s3 = wszes(4)

      nor=(no*(no+1))/2
      nvr=(nv*(nv+1))/2

      !Determine the offsets in the alpha and gamma indices, which arise from
      !non uniform batch distributions, e.g. consider:
      !Case 1:   gamma _____
      !               |\    |
      !               |_\___|
      !  alpha        |  \  |
      !               |___\_| here an offset in gamma for the second gamma batch is
      !               needed, since the elements to be considered do not
      !               begin with the first element in that batch
      goffs=0
      if(fa-fg>0)goffs=fa-fg
      !Case 2:   gamma ______
      !               |\ |   |
      !               | \|   |
      !  alpha        |  |\  |
      !               |__|_\_| here an offset in alpha for the second alpha batch is
      !               needed, since the elements to be considered do not
      !               begin with the first element in that batch
      aoffs=0
      if(fg-fa>0)aoffs=fg-fa

      !Determine the dimension of the triangular part in the current batch
      !and the total number of elements (tred) in the triangular plus
      !rectangular parts
      tred=0
      tlen=min(min(min(la,lg),fg+lg-fa),fa+la-fg)

      !calculate amount of triangular elements in the current batch
      if(fa+la-1>=fg)tred= tlen*(tlen+1)/2
      !add the rectangular contribution if lg is larger than tlen 
      if(fa>=fg.and.fg+lg-fa-tlen>0)tred=tred+(fg+lg-fa-tlen)*la
      if(fa<fg)tred=tred+lg*aoffs
      if(fa<fg.and.fa+la>fg) tred=tred+(lg-tlen)*(la-aoffs)
      !if only rectangular
      if(tlen<=0)then
         tlen=0
         aoffs=0
         goffs=0
         tred=la*lg
      endif

      select case(s)
#ifdef DIL_ACTIVE
      case(1) !`DIL: Scheme 1 only
         !!SYMMETRIC COMBINATION:
         ! (w2): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
         call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
         !(w0):I+ [delta alpha<=gamma beta] <= (w1):I [alpha beta gamma delta] + (w1):I[alpha delta gamma beta]
         call get_I_plusminus_le(w2,'+',fa,fg,la,lg,nb,tlen,tred,goffs,s2,1,tred)
         !(w0):I+ [delta alpha<=gamma c] = (w0):I+ [delta alpha<=gamma beta] * Lambda^h[beta c]
         call dgemm('t','n',nb*tred,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nb*tred)
         !(w2):I+ [alpha<=gamma c d] = (w0):I+ [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
         call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w0,nb,yv,nb,0.0E0_realk,w2,nv*tred)
         !(w0):I+ [alpha<=gamma c>=d] <= (w0):I+ [alpha<=gamma c d] 
         call get_I_cged(w0,w2,tred,nv)
         !(w3.1):sigma+ [alpha<=gamma i>=j] = (w0):I+ [alpha<=gamma c>=d] * (w0):t+ [c>=d i>=j]
!        call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w0,tred,tpl,nvr,0.0E0_realk,w3,tred)
         if(DIL_DEBUG) then !`DIL: Tensor contraction 6
          write(DIL_CONS_OUT,'("#DEBUG(DIL): Process ",i6,"[",i6,"] starting tensor contraction 6:",3(1x,i7))')&
          &infpar%lg_mynum,infpar%mynum,nor,nvr,tred
         endif
         call dil_array_init(w3,i8*tred*nor)
         tcs='D(z,y)+=L(z,x)*R(y,x)'
         call dil_clean_tens_contr(tch)
         tens_rank=2; tens_dims(1:tens_rank)=(/int(tred,INTD),int(nor,INTD)/)
         call dil_set_tens_contr_args(tch,'d',errc,tens_rank,tens_dims,w3)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC6: DA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC6: Destination arg set failed!',-1)
         tens_rank=2; tens_dims(1:tens_rank)=(/int(tred,INTD),int(nvr,INTD)/)
         call dil_set_tens_contr_args(tch,'l',errc,tens_rank,tens_dims,w0)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC6: LA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC6: Left arg set failed!',-1)
         call dil_set_tens_contr_args(tch,'r',errc,tens_distr=tpl)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC6: RA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC6: Right arg set failed!',-1)
         call dil_set_tens_contr_spec(tch,tcs,errc,alpha=0.5E0_realk)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC6: CC: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC6: Contr spec set failed!',-1)
         dil_mem=dil_get_min_buf_size(tch,errc)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC6: BS: ',infpar%lg_mynum,errc,dil_mem
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC6: Buf size set failed!',-1)
         dil_mem=dil_buf_size*8_INTL; call dil_prepare_buffer(tch,dil_mem,errc,dil_buffer)
         if(errc.ne.0) call lsquit('ERROR(ccsd_residual_integral_driven): TC6: Buf alloc failed!',-1)
         call dil_tensor_contract(tch,DIL_TC_EACH,dil_mem,errc,locked=dil_lock_out)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC6: TC: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC6: Tens contr failed!',-1)

         !!ANTI-SYMMETRIC COMBINATION:
         ! (w2): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
         call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
         !(w0):I- [delta alpha<=gamma beta] <= (w1):I [alpha beta gamma delta] + (w1):I[alpha delta gamma beta]
         call get_I_plusminus_le(w2,'-',fa,fg,la,lg,nb,tlen,tred,goffs,s2,1,tred)
         !(w0):I- [delta alpha<=gamma c] = (w0):I- [delta alpha<=gamma beta] * Lambda^h[beta c]
         call dgemm('t','n',nb*tred,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nb*tred)
         !(w2):I- [alpha<=gamma c d] = (w0):I- [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
         call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w0,nb,yv,nb,0.0E0_realk,w2,nv*tred)
         !(w0):I- [alpha<=gamma c<=d] <= (w2):I- [alpha<=gamma c d]
         call get_I_cged(w0,w2,tred,nv)
         !(w3.2):sigma- [alpha<=gamma i<=j] = (w0):I- [alpha<=gamma c>=d] * (w0):t- [c>=d i>=j]
!        call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w0,tred,tmi,nvr,0.0E0_realk,w3(tred*nor+1),tred)
         if(DIL_DEBUG) then !`DIL: Tensor contraction 7
          write(DIL_CONS_OUT,'("#DEBUG(DIL): Process ",i6,"[",i6,"] starting tensor contraction 7:",3(1x,i7))')&
          &infpar%lg_mynum,infpar%mynum,nor,nvr,tred
         endif
         call dil_array_init(w3(tred*nor+1:),i8*tred*nor)
         tcs='D(z,y)+=L(z,x)*R(y,x)'
         call dil_clean_tens_contr(tch)
         tens_rank=2; tens_dims(1:tens_rank)=(/int(tred,INTD),int(nor,INTD)/)
         call dil_set_tens_contr_args(tch,'d',errc,tens_rank,tens_dims,w3(tred*nor+1:))
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC7: DA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC7: Destination arg set failed!',-1)
         tens_rank=2; tens_dims(1:tens_rank)=(/int(tred,INTD),int(nvr,INTD)/)
         call dil_set_tens_contr_args(tch,'l',errc,tens_rank,tens_dims,w0)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC7: LA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC7: Left arg set failed!',-1)
         call dil_set_tens_contr_args(tch,'r',errc,tens_distr=tmi)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC7: RA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC7: Right arg set failed!',-1)
         call dil_set_tens_contr_spec(tch,tcs,errc,alpha=0.5E0_realk)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC7: CC: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC7: Contr spec set failed!',-1)
         dil_mem=dil_get_min_buf_size(tch,errc)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC7: BS: ',infpar%lg_mynum,errc,dil_mem
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC7: Buf size set failed!',-1)
         dil_mem=dil_buf_size*8_INTL; call dil_prepare_buffer(tch,dil_mem,errc,dil_buffer)
         if(errc.ne.0) call lsquit('ERROR(ccsd_residual_integral_driven): TC7: Buf alloc failed!',-1)
         call dil_tensor_contract(tch,DIL_TC_EACH,dil_mem,errc,locked=dil_lock_out)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC7: TC: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_a22_and_prepb22_terms_exd): TC7: Tens contr failed!',-1)
#endif
      case default
         call lsquit("ERROR(get_a22_and_prepb22_terms_exd): wrong scheme on input",-1)
      end select

      !COMBINE THE TWO SIGMAS OF W3 IN W2
      !(w2):sigma[alpha<=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] + 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
      !(w2):sigma[alpha>=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] - 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
      s0 = wszes(1)
      s2 = wszes(3)
      s3 = wszes(4)
      if(s==1) then; scheme=sch_sym; else; scheme=s; endif  !`DIL: remove
      call combine_and_transform_sigma(om2,w0,w2,w3,s0,s2,s3,xv,xo,sio4,nor,tlen,tred,fa,fg,la,lg,&
         &no,nv,nb,goffs,aoffs,scheme,lo,twork,tcomm,order=order,rest_occ_om2=rest_occ_om2,scal=scal,& !`DIL: scheme --> s
         &sio4_ilej=(s/=2.and.s/=1),o2tens=o2ilej)

      call time_start_phase(PHASE_WORK, at=twork)
   end subroutine get_a22_and_prepb22_terms_exd
#endif

   !> \brief Combine sigma matrixes in symmetric and antisymmetric combinations 
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine combine_and_transform_sigma(omega,w0,w2,w3,s0,s2,s3,xvirt,xocc,sio4,nor,tlen,tred,fa,fg,la,lg,&
         &no,nv,nb,goffs,aoffs,s,lock_outside,twork,tcomm,order,rest_occ_om2,scal,act_no,sio4_ilej,query,o2tens)
      implicit none
      !> size of w0
      integer(kind=8),intent(inout)   :: s0,s2,s3
      !\> omega should be the residual matrix which contains the second parts
      !of the A2 and B2 term
      !real(realk),intent(inout) :: omega(nv*nv*no*no)
      type(tensor),intent(inout) :: omega
      !> w0 is just some workspace on input
      real(realk),intent(inout) :: w0(s0)
      !> w2 is just some workspace on input
      real(realk),intent(inout),target :: w2(s2)
      !> w3 contains the symmetric and antisymmetric combinations 
      real(realk),intent(inout) :: w3(s3)
      !> sio4 are the reduced o4 integrals which are used to calculate the B2.2
      !contribution after the loop, update them in the loops
      type(tensor),intent(inout) :: sio4
      !> Lambda p virutal part
      real(realk),intent(in) :: xvirt(:)
      !> Lambda p occupied part
      real(realk),intent(in) :: xocc(:)
      !> number of reduced occupied indices 
      integer,intent(in) :: nor
      !> total number of upper triangular elements in the batch
      integer,intent(in) :: tlen
      !> first element of alpha and gamma in the current batch
      integer,intent(in) :: fa,fg
      !> length of alpha and gamma in the current batch
      integer,intent(in) :: la, lg
      !> number of triangular elements in the batch
      integer,intent(in) :: tred
      !>number of occupied, virtual and ao indices
      integer,intent(in) :: no, nv, nb
      !>offsets in the currnt alpha and gamma batches to get to the first upper
      !triangular element
      integer,intent(in) :: goffs, aoffs
      !> scheme
      integer,intent(in) :: s
      logical,intent(in) :: lock_outside
      integer,optional,intent(in)  :: order(4),act_no
      !restricted i<=j in the omega2 and or sio4
      logical,optional, intent(in) :: rest_occ_om2, sio4_ilej, query
      real(realk), optional, intent(in) :: scal
      type(tensor), optional, intent(inout) :: o2tens
      !timing information
      real(realk) :: twork,tcomm
      !> the doubles amplitudes
      !real(realk),intent(in) :: amps(nv*nv*no*no)
      !type(tensor),intent(in) :: amps
      real(realk) :: scaleitby
      integer(kind=8)       :: pos,pos2,pos21,i,j,dim_big,dim_small,ttri,tsq,nel2cp,ncph,pos1
      integer ::occ,gamm,alpha,case_sel,full1,full2,offset1,offset2
      integer :: l1,l2,lsa,lsg,gamm_i_b,a,b,full1T,full2T,jump,ft1,ft2
      logical               :: second_trafo_step
      real(realk),pointer   :: dumm(:)
      integer               :: mv((nv*nv)/2),st,dims(2),no2
      real(realk),pointer   :: source(:,:),drain(:,:)
      integer(kind=ls_mpik) :: mode
      integer(kind=long)    :: o2v2
      logical               :: rest_o2_occ, rest_sio4,qu
#ifdef VAR_PTR_RESHAPE
      real(realk), pointer, contiguous  :: h1(:,:,:,:), t1(:,:,:), h(:)
#else
      real(realk), pointer  :: h1(:,:,:,:), t1(:,:,:), h(:)
#endif
      !$ integer, external  :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads
#ifdef DIL_ACTIVE
!{`DIL:
     character(256):: tcs
     type(dil_tens_contr_t):: tch0
     integer(INTL):: dil_mem,l0
     integer(INTD):: i0,i1,i2,i3,errc,tens_rank,tens_dims(MAX_TENSOR_RANK),tens_bases(MAX_TENSOR_RANK)
     integer(INTD):: ddims(MAX_TENSOR_RANK),ldims(MAX_TENSOR_RANK),rdims(MAX_TENSOR_RANK)
     integer(INTD):: dbase(MAX_TENSOR_RANK),lbase(MAX_TENSOR_RANK),rbase(MAX_TENSOR_RANK)
#endif
     real(realk):: r0
!}

      if(s == 1 .and. .not. present(o2tens))then
         call lsquit("ERROR(combine_and_transform_sigma): for scheme 1 we need the o2tens argument",-1)
      endif

      rest_o2_occ = .false.
      if(present(rest_occ_om2 ))rest_o2_occ = rest_occ_om2
      no2 = no
      if(present(act_no))no2 = act_no
      rest_sio4 = .true.
      if(present(sio4_ilej))rest_sio4 = sio4_ilej
      qu = .false.
      if(present(query))qu = query

      o2v2 = int((i8*no)*no*nv*nv,kind=long)
#ifdef VAR_MPI
      mode = MPI_MODE_NOCHECK
#endif

      scaleitby=1.0E0_realk
      if(present(scal)) scaleitby = scal
      second_trafo_step=.false.
      lsa=fa+la-1
      lsg=fg+lg-1

      !building the sigma matrices is a complicated matter, due to non-uniform
      !batchsizes. A distiction between the different possible cases has to be
      !made and the elements copied respectively
      case_sel=0
      if(fa>=fg)then
         if(lsa>=lsg)then
            !case 0 is a triagular submatrix:
            !                                full2=tlen
            ! ---|                          |----|
            ! \  |                          | \  |  full1=tlen
            !  \ |   build full from that   |  \ |
            !   \|                          |---\|
            case_sel=1
            full1=tlen
            full2=tlen
            full1T=tlen
            full2T=tlen
            l1=0
            l2=0
         else if(lsa<lsg)then
            !case 1 is a  submatrix of the type:
            !                                   full2=lg-goffs
            ! ------|                          |------|
            ! \     |                          | \    |
            !  \    |   build full from that   |  \   |full1=la
            !   \---|                          |___\__|
            !                                  |   |0 |
            !            full1T=lg-goffs-tlen  |___|__|
            case_sel=2                         !full2T=la 
            full1=la
            full2=lg-goffs
            full1T=lg-goffs-tlen
            full2T=la
            l1=la
            l2=0
            second_trafo_step=.true.
         endif
      else
         if(lsa>=lsg)then
            !case 2 is a  submatrix of the type:
            !                                 full2=lg
            ! |---|                         |------|
            ! |   |                         |0 |   |
            ! |   |  build full from that   |__|   | full1=aoffs+tlen
            !  \  |                         |  |\  |
            !   \ |                         |  | \ |
            !    \|     full1T=lg           |__|__\|
            !                               full2T=aoffs
            case_sel=3
            full1=aoffs+tlen
            full2=tlen
            full1T=tlen
            full2T=aoffs
            l1=aoffs
            l2=aoffs
            second_trafo_step=.true.
         else if(lsa<lsg)then
            if(lsa>=fg)then
               !case 3 is a  submatrix of the type: 
               !                                       full2=lg
               ! |-----|                         |--------|
               ! |     |                         |0 |     |
               ! |     |  build full from that   |__|     |
               !  \    |                         |  |\    |  full1 = la
               !   \   |                         |  | \   |
               !    \--|                         |  |__\__|
               !                     full1T=lg   |     | 0|
               !                                 |_____|__|
               !                                    full2T = la
               case_sel=4
               full1=la
               full2=lg
               full1T=lg
               full2T=la
               l1=aoffs
               l2=aoffs
               second_trafo_step=.true.
            else
               !case 4 is a  submatrix of dimensions la lg, and can thus be transformed in
               !full anyways
               case_sel=5
               full1=la
               full2=lg
               full1T=lg
               full2T=la
               l1=la
               l2=0
               second_trafo_step=.true.
            endif
         endif
      endif
      if(case_sel==0)call lsquit("ERROR(combine_and_transform_sigma):case not known",DECinfo%output)

      !print *,"-------------------------------------------------------------------------------------"
      !print *,"       case        dim1        dim2            la          lg        tlen        tred"
      !print *,case_sel,full1,full2,la,lg,tlen,tred
      !print *,"-------------------------------------------------------------------------------------"

      if( qu )then

         !ATTENTION: KEEP UP TO DATE
         !should always be zeroed outside, else there are stack allocation
         !problems
         !s0 = 0
         !s2 = 0
         !s3 = 0

         !w0:
         if(second_trafo_step)then
            s0 = max(s0,(i8*nor)*full1*full2+(i8*nor)*full1T*full2T)
         else
            s0 = max(s0,(i8*nor)*full1*full2)
         endif

         !w2:
         s2 = max(s2,(i8*nor)*full1*full2)
         s2 = max(s2,(i8*nv*nv)*nor)
         if(.not.rest_o2_occ)then
            s2 = max(s2,(i8*nv*nv)*no*no)
         endif
         if(second_trafo_step)then
            s2 = max(s2,(i8*full1T)*full2T*nor)
            s2 = max(s2,(i8*nv)*nv*nor)
         endif
         if( .not. rest_sio4 )then
            s2 = max(s2,(i8*no2*no2)*nor)
         endif

         !w3:
         s3 = max(s3,(i8*nor)*full1*full2)
         s3 = max(s3,(i8*nv*nor*full1)*full2)
         s3 = max(s3,(i8*no2)*nor*full1)
         if(second_trafo_step)then
            s3 = max(s3,(i8*nv*nor)*full1T)
        endif

      else

         !Zero the elements to update for testing, not needed in a performance
         !implementation
         pos=(i8*nor)*full1*full2
         if(second_trafo_step)pos=pos+(i8*nor)*full1T*full2T
         w0(1_long:pos)=0.0E0_realk

         !set required variables
         ttri      = tlen * (tlen+1)/2
         tsq       = tlen *  tlen
         pos       = 1
         occ       = 1
         dim_big   = (i8 * full1 ) * full2
         dim_small = (i8 * full1T) * full2T

         !$OMP PARALLEL DEFAULT(NONE)&
         !$OMP SHARED(w0,w3,case_sel,nor,goffs,lg,la,full1,full1T,ttri,tred,&
         !$OMP full2,full2T,tlen,l1,second_trafo_step,aoffs,dim_big,dim_small,l2)&
         !$OMP PRIVATE(occ,gamm,gamm_i_b,pos,nel2cp,pos2,jump,ft1,ft2,ncph,pos21,&
         !$OMP dims,i)
         !$OMP DO
         do occ=1,nor
            do gamm=1,lg-goffs
               gamm_i_b=gamm+goffs
               !SYMMETRIC COMBINATION OF THE SIGMAS

               !calculate the old position
               !**************************
               if(case_sel==3.or.case_sel==4)then
                  pos=1+(gamm      -1)*full1+(occ-1)*dim_big
               else
                  pos=1+(gamm+aoffs-1)*full1+(occ-1)*dim_big
               endif

               !calculate the new position
               !**************************
               if(gamm>tlen)then

                  !get the elements from the rectangular part of the batch

                  nel2cp = l1
                  pos2   = 1_long + ttri+(gamm-tlen-1)*(la-aoffs)+(occ-1)*((la-aoffs)*(lg-goffs-tlen)+ttri)

                  if(case_sel==4)then
                     nel2cp = nel2cp + tlen
                     pos2   = pos2   + tlen * aoffs + (gamm-tlen-1) * (la-tlen)+(occ-1)*(aoffs*lg)
                  endif

                  if(second_trafo_step)then
                     jump = full1T
                     ft1  = full1T
                     ft2  = full2T
                  else
                     jump = full1
                     ft1  = full1
                     ft2  = full2
                  endif

               else
                  !get the elements from the triangular part of the batch
                  pos2   = 1  + (gamm*(gamm-1)/2)+(gamm-1)*aoffs+(occ-1)*tred
                  nel2cp = l2 + gamm
                  jump   = full1
                  ft1    = full1
                  ft2    = full2
               endif

               w0(pos:pos+nel2cp-1) = w3(pos2:pos2+nel2cp-1)

               !get corresponding position in sigma- and add to output
               pos21=pos2+tred*nor
               w0(pos:pos+nel2cp-1) =w0(pos:pos+nel2cp-1) + w3(pos21:pos21+nel2cp-1)    

               !ANTI-SYMMETRIC COMBINATION OF THE SIGMAS
               pos = gamm+aoffs+(occ-1)*ft1*ft2

               if(second_trafo_step.and.gamm>tlen) pos = pos+full1*full2*nor-tlen

               if(case_sel==3.or.case_sel==4)then
                  !fill diagonal part
                  if(gamm>tlen)then
                     ncph=0
                  else
                     ncph=gamm
                  endif

                  !because of the intrinsic omp-parallelizaton of daxpy the following
                  !lines replace the daxpy calls
                  pos = aoffs+gamm+(occ-1)*dim_big
                  do i=0,ncph-1
                     w0(pos+i*full1)=w0(pos+i*full1) + w3(pos2 +aoffs+i)
                     w0(pos+i*full1)=w0(pos+i*full1) - w3(pos21+aoffs+i)
                  enddo

                  !fill small matrix
                  if(gamm>tlen)then
                     ncph=nel2cp
                  else
                     ncph=nel2cp-gamm
                  endif
                  pos = gamm+(occ-1)*full1T*full2T+dim_big*nor
                  do i=0,ncph-1
                     w0(pos+i*full1T) = w0(pos+i*full1T) + w3(pos2 +i)
                     w0(pos+i*full1T) = w0(pos+i*full1T) - w3(pos21+i)
                  enddo
               else
                  do i=0,nel2cp-1
                     w0(pos+i*jump) = w0(pos+i*jump) + w3(pos2 +i)
                     w0(pos+i*jump) = w0(pos+i*jump) - w3(pos21+i)
                  enddo
               endif
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL


         !Print the individual contributions
         !if(.false.)then
         !  do occ=1,nor
         !    print *,"symmetric"
         !    print *,w3(1+(occ-1)*tred:tred+(occ-1)*tred)
         !    !print *,"antisymmetric"
         !    !print *,w3(nor*tred+1+(occ-1)*tred:nor*tred+tred+(occ-1)*tred)
         !    print *,"a g -> 1. and l.",1+(occ-1)*full1*full2,full1*full2+(occ-1)*full1*full2
         !    do alpha=1,full1
         !      print *,(w0(alpha+(gamm-1)*full1+(occ-1)*full1*full2),gamm=1,full2)
         !    enddo
         !    print *,""
         !    if(second_trafo_step)then
         !      print *,"g a -> 1. and l.",dim_big*nor+1+(occ-1)*full1T*full2T,dim_big*nor+full1T*full2T+(occ-1)*full1T*full2T
         !      do gamm=1,full1T
         !        print *,(w0(dim_big*nor+gamm+(alpha-1)*full1T+(occ-1)*full1T*full2T),alpha=1,full2T)
         !      enddo
         !      print *,""
         !    endif
         !    print *,""
         !  enddo
         !endif

         !add up the contributions for the sigma [ alpha gamma ] contributions
         ! get the order sigma[ gamma i j alpha ]
         call mat_transpose(full1,full2*nor,1.0E0_realk,w0,0.0E0_realk,w2)
         !transform gamma -> b
         call dgemm('t','n',nv,nor*full1,full2,1.0E0_realk,xvirt(fg+goffs),nb,w2,full2,0.0E0_realk,w3,nv)
         !transform alpha -> a , order is now sigma [ a b i j]
         if(s == 1) then !`DIL: Tensor contraction 8
#ifdef DIL_ACTIVE
          if(DIL_DEBUG) then
           write(DIL_CONS_OUT,'("#DEBUG(DIL): Process ",i6,"[",i6,"] starting tensor contraction 8:")')&
            &infpar%lg_mynum,infpar%mynum
          endif
          tcs='D(a,b,i)+=L(u,a)*R(b,i,u)'
          call dil_clean_tens_contr(tch0)
          call dil_set_tens_contr_args(tch0,'d',errc,tens_distr=o2tens)
          if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC8: DA: ',infpar%lg_mynum,errc
          if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC8: Destination arg set failed!',-1)
          tens_rank=2; tens_dims(1:tens_rank)=(/nb,nv/)
          call dil_set_tens_contr_args(tch0,'l',errc,tens_rank,tens_dims,xvirt)
          if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC8: LA: ',infpar%lg_mynum,errc
          if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC8: Left arg set failed!',-1)
          tens_rank=3; tens_dims(1:tens_rank)=(/nv,nor,full1/); tens_bases(1:tens_rank)=(/0,0,fa-1/)
          call dil_set_tens_contr_args(tch0,'r',errc,tens_rank,tens_dims,w3,tens_bases)
          if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC8: RA: ',infpar%lg_mynum,errc
          if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC8: Right arg set failed!',-1)
          call dil_set_tens_contr_spec(tch0,tcs,errc,&
               &ldims=(/int(full1,INTD),int(nv,INTD)/),lbase=(/int(fa-1,INTD),0_INTD/),&
               &rdims=(/int(nv,INTD),int(nor,INTD),int(full1,INTD)/),rbase=(/0_INTD,0_INTD,int(fa-1,INTD)/),&
               &alpha=scaleitby)
          if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC8: CC: ',infpar%lg_mynum,errc
          if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC8: Contr spec set failed!',-1)
          dil_mem=dil_get_min_buf_size(tch0,errc)
          if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC8: BS: ',infpar%lg_mynum,errc,dil_mem
          if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC8: Buf size set failed!',-1)
          dil_mem=dil_buf_size*8_INTL; call dil_prepare_buffer(tch0,dil_mem,errc,dil_buffer)
          if(errc.ne.0) call lsquit('ERROR(ccsd_residual_integral_driven): TC8: Buf alloc failed!',-1)
          call dil_tensor_contract(tch0,DIL_TC_EACH,dil_mem,errc,locked=lock_outside)
          if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC8: TC: ',infpar%lg_mynum,errc
          if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC8: Tens contr failed!',-1)
#else
          call lsquit('ERROR(combine_and_transform_sigma): This part (7) of Scheme 1 requires DIL backend!',-1)
#endif
         else
            call dgemm('t','t',nv,nv*nor,full1,1.0E0_realk,xvirt(fa),nb,w3,nor*nv,0.0E0_realk,w2,nv)

            if(.not.rest_o2_occ)then
               !square up the contributions if the residual itself has no restricions
               !in the indices i and j
               call squareup_block_triangular_squarematrix(w2,nv,no,do_block_transpose = .true.)

               if(s==4.or.s==3)then
                  if( present(order) )then
                     call array_reorder_4d(scaleitby,w2,nv,nv,no,no,order,1.0E0_realk,omega%elm1)
                  else
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
                     call assign_in_subblocks(omega%elm1,'+',w2,o2v2,scal2=scaleitby)
#else
                     !$OMP WORKSHARE
                     omega%elm1(1_long:o2v2) = omega%elm1(1_long:o2v2) + scaleitby * w2(1_long:o2v2)
                     !$OMP END WORKSHARE
#endif
                  endif
               else if(s==2)then
#ifdef VAR_MPI
                  if( .not.alloc_in_dummy.and.lock_outside )call tensor_lock_wins(omega,'s',mode)
                  !$OMP WORKSHARE
                  w2(1_long:o2v2) = scaleitby*w2(1_long:o2v2)
                  !$OMP END WORKSHARE
                  call time_start_phase(PHASE_COMM, at=twork)
                  call tensor_add(omega,1.0E0_realk,w2,wrk=w3,iwrk=s3)
                  call time_start_phase(PHASE_WORK, at=tcomm)
#endif
               endif
            else
#ifdef VAR_LSDEBUG
               if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)not implemented for other schemes than 4",-1)
               if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)only implemented for pairs for now",-1)
#endif
               pos1 = 1
               call array_reorder_2d(scaleitby,w2(pos1:nv*nv+pos1-1),nv,nv,[2,1],1.0E0_realk,omega%elm1)

            endif
         endif

         !If the contributions are split in terms of the sigma matrix this adds the
         !sigma [gamma alpha] contributions
         if(second_trafo_step)then
            !get the order sigma [gamma i j alpha]
            pos2 = 1+full1*full2*nor
            if(case_sel==3.or.case_sel==4)then
               l1=fa
               l2=fg
            else
               l1=fa
               l2=fg+goffs+tlen
            endif
            call mat_transpose(full1T,full2T*nor,1.0E0_realk,w0(pos2:full1T*full2T*nor+pos2-1),0.0E0_realk,w2)

#ifdef  VAR_MPI
            if( lock_outside .and. (s==2.or.s==1) )then
               call time_start_phase(PHASE_COMM, at=twork)
               if( alloc_in_dummy )then
                  call lsmpi_win_flush(omega%wi(1),local=.true.)
               else
                  call tensor_unlock_wins(omega,.true.)
               endif
               call time_start_phase(PHASE_WORK, at=tcomm)
            endif
#endif 

            !transform gamma -> a
            call dgemm('t','n',nv,nor*full1T,full2T,1.0E0_realk,xvirt(l1),nb,w2,full2T,0.0E0_realk,w3,nv)
            !transform alpha -> , order is now sigma[b a i j]
            if(s == 1)then !`DIL: Tensor contraction 9
#ifdef DIL_ACTIVE
             if(DIL_DEBUG) then
              write(DIL_CONS_OUT,'("#DEBUG(DIL): Process ",i6,"[",i6,"] starting tensor contraction 9:")')&
               &infpar%lg_mynum,infpar%mynum
             endif
             tcs='D(a,b,i)+=L(u,a)*R(b,i,u)'
             call dil_clean_tens_contr(tch0)
             call dil_set_tens_contr_args(tch0,'d',errc,tens_distr=o2tens)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC9: DA: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC9: Destination arg set failed!',-1)
             tens_rank=2; tens_dims(1:tens_rank)=(/nb,nv/)
             call dil_set_tens_contr_args(tch0,'l',errc,tens_rank,tens_dims,xvirt)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC9: LA: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC9: Left arg set failed!',-1)
             tens_rank=3; tens_dims(1:tens_rank)=(/nv,nor,full1T/); tens_bases(1:tens_rank)=(/0,0,l2-1/)
             call dil_set_tens_contr_args(tch0,'r',errc,tens_rank,tens_dims,w3,tens_bases)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC9: RA: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC9: Right arg set failed!',-1)
             call dil_set_tens_contr_spec(tch0,tcs,errc,&
                  &ldims=(/int(full1T,INTD),int(nv,INTD)/),lbase=(/int(l2-1,INTD),0_INTD/),&
                  &rdims=(/int(nv,INTD),int(nor,INTD),int(full1T,INTD)/),rbase=(/0_INTD,0_INTD,int(l2-1,INTD)/),&
                  &alpha=scaleitby)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC9: CC: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC9: Contr spec set failed!',-1)
             dil_mem=dil_get_min_buf_size(tch0,errc)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC9: BS: ',infpar%lg_mynum,errc,dil_mem
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC9: Buf size set failed!',-1)
             dil_mem=dil_buf_size*8_INTL; call dil_prepare_buffer(tch0,dil_mem,errc,dil_buffer)
             if(errc.ne.0) call lsquit('ERROR(ccsd_residual_integral_driven): TC9: Buf alloc failed!',-1)
             call dil_tensor_contract(tch0,DIL_TC_EACH,dil_mem,errc,locked=lock_outside)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC9: TC: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC9: Tens contr failed!',-1)
#else
             call lsquit('ERROR(combine_and_transform_sigma): This part (8) of Scheme 1 requires DIL backend!',-1)
#endif
            else
               call dgemm('t','t',nv,nv*nor,full1T,1.0E0_realk,xvirt(l2),nb,w3,nor*nv,0.0E0_realk,w2,nv)

               if(.not.rest_o2_occ)then
                  call squareup_block_triangular_squarematrix(w2,nv,no,do_block_transpose = .true.)
                  if(s==4.or.s==3)then
                     if( present(order) )then
                        call array_reorder_4d(scaleitby,w2,nv,nv,no,no,order,1.0E0_realk,omega%elm1)
                     else
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
                        call assign_in_subblocks(omega%elm1,'+',w2,o2v2,scal2=scaleitby)
#else
                        !$OMP WORKSHARE
                        omega%elm1(1_long:o2v2) = omega%elm1(1_long:o2v2) + scaleitby * w2(1_long:o2v2)
                        !$OMP END WORKSHARE
#endif
                     endif
                  else if(s==2)then
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
                     call assign_in_subblocks(w2,'=',w2,o2v2,scal2=scaleitby)
#else
                     !$OMP WORKSHARE
                     w2(1_long:o2v2) = scaleitby*w2(1_long:o2v2)
                     !$OMP END WORKSHARE
#endif
                     call time_start_phase(PHASE_COMM, at=twork)
                     call tensor_add(omega,1.0E0_realk,w2,wrk=w3,iwrk=s3)
                     call time_start_phase(PHASE_WORK, at=tcomm)
                  endif
               else
#ifdef VAR_LSDEBUG
                  if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)not implemented&
                     & for other schemes than 4",-1)
                  if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)only implemented&
                     & for pairs for now",-1)
#endif
                  pos1 = 1
                  call array_reorder_2d(scaleitby,w2(pos1:nv*nv+pos1-1),nv,nv,[2,1],1.0E0_realk,omega%elm1)
               endif
            endif
         endif


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Construct the B2 term from the intermediates in w0
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !add up the contributions for the sigma [ alpha gamma ] contributions
         ! get the order sigma[ gamma i j alpha ]
         call mat_transpose(full1,full2*nor,1.0E0_realk,w0,0.0E0_realk,w2)
#ifdef VAR_MPI
         if(lock_outside.and.(s==2.or.s==1))then
            call time_start_phase(PHASE_COMM, at=twork)
            if( alloc_in_dummy )then
               call lsmpi_win_flush(omega%wi(1),local=.true.)
            else
               call tensor_unlock_wins(omega,.true.)
            endif
            call time_start_phase(PHASE_WORK, at=tcomm)
         endif
#endif
         !transform gamma -> l
         call dgemm('t','n',no2,nor*full1,full2,1.0E0_realk,xocc(fg+goffs),nb,w2,full2,0.0E0_realk,w3,no2)
         !transform alpha -> a , order is now sigma [ k l i j]
         if( rest_sio4 )then
            call dgemm('t','t',no2,no2*nor,full1,1.0E0_realk,xocc(fa),nb,w3,nor*no2,1.0E0_realk,sio4%elm1,no2)
         else
            select case(s)
            case(1)
#ifdef DIL_ACTIVE
             if(DIL_DEBUG) then !`DIL: Tensor contraction 10
              write(DIL_CONS_OUT,'("#DEBUG(DIL): Process ",i6,"[",i6,"] starting tensor contraction 10:")')&
               &infpar%lg_mynum,infpar%mynum
             endif
             tcs='D(k,l,i)+=L(u,k)*R(l,i,u)'
             call dil_clean_tens_contr(tch0)
             call dil_set_tens_contr_args(tch0,'d',errc,tens_distr=sio4)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC10: DA: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC10: Destination arg set failed!',-1)
             tens_rank=2; tens_dims(1:tens_rank)=(/nb,no2/)
             call dil_set_tens_contr_args(tch0,'l',errc,tens_rank,tens_dims,xocc)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC10: LA: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC10: Left arg set failed!',-1)
             tens_rank=3; tens_dims(1:tens_rank)=(/no2,nor,full1/); tens_bases(1:tens_rank)=(/0,0,fa-1/)
             call dil_set_tens_contr_args(tch0,'r',errc,tens_rank,tens_dims,w3,tens_bases)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC10: RA: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC10: Right arg set failed!',-1)
             call dil_set_tens_contr_spec(tch0,tcs,errc,&
                  &ldims=(/int(full1,INTD),int(no2,INTD)/),lbase=(/int(fa-1,INTD),0_INTD/),&
                  &rdims=(/int(no2,INTD),int(nor,INTD),int(full1,INTD)/),rbase=(/0_INTD,0_INTD,int(fa-1,INTD)/),&
                  &alpha=1E0_realk)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC10: CC: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC10: Contr spec set failed!',-1)
             dil_mem=dil_get_min_buf_size(tch0,errc)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC10: BS: ',infpar%lg_mynum,errc,dil_mem
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC10: Buf size set failed!',-1)
             dil_mem=dil_buf_size*8_INTL; call dil_prepare_buffer(tch0,dil_mem,errc,dil_buffer)
             if(errc.ne.0) call lsquit('ERROR(ccsd_residual_integral_driven): TC10: Buf alloc failed!',-1)
             call dil_tensor_contract(tch0,DIL_TC_EACH,dil_mem,errc,locked=lock_outside)
             if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC10: TC: ',infpar%lg_mynum,errc
             if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC10: Tens contr failed!',-1)
#else
             call lsquit('ERROR(combine_and_transform_sigma): This part (9) of Scheme 1 requires DIL backend!',-1)
#endif
            case(2)
               call dgemm('t','t',no2,no2*nor,full1,1.0E0_realk,xocc(fa),nb,w3,nor*no2,0.0E0_realk,w2,no2)

               call squareup_block_triangular_squarematrix(w2,no,no,do_block_transpose = .true.)
#ifdef VAR_MPI
               if( lock_outside .and..not. alloc_in_dummy )call tensor_lock_wins(sio4,'s',mode)
               call time_start_phase(PHASE_COMM, at=twork)
               call tensor_add(sio4,1.0E0_realk,w2,wrk=w3,iwrk=s3)
               if( alloc_in_dummy )then
                  call lsmpi_win_flush(sio4%wi(1),local=.true.)
               else
                  if( lock_outside )call tensor_unlock_wins(sio4,.true.)
               endif
               call time_start_phase(PHASE_WORK, at=tcomm)
#endif
            case(3,4)
               call dgemm('t','t',no2,no2*nor,full1,1.0E0_realk,xocc(fa),nb,w3,nor*no2,0.0E0_realk,w2,no2)
#ifdef VAR_PTR_RESHAPE
               t1(1:no2,1:no2,1:nor)     => w2
               !this is a workaround for PGI 15.5 and 15.7
               h => sio4%elm1
               h1(1:no,1:no,1:no2,1:no2) => h
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
               call c_f_pointer(c_loc(w2(1)),t1,[no2,no2,nor])
               call c_f_pointer(c_loc(sio4%elm1(1)),h1,[no,no,no2,no2])
#else
               call lsquit('ERROR(combine_and_transform_sigma): unable to reshape pointers!',-1)
#endif
               do j=no,1,-1
                  do i=j,1,-1
                     call array_reorder_2d(1.0E0_realk,t1(:,:,i+j*(j-1)/2),no2,no2,[2,1],1.0E0_realk,h1(i,j,:,:))
                     if(i /= j)then
                        h1(j,i,:,:) = h1(j,i,:,:) +  t1(:,:,i+j*(j-1)/2)
                     endif
                  enddo
               enddo
            end select
         endif


         !If the contributions are split in terms of the sigma matrix this adds the
         !sigma [gamma alpha i j] contributions
         if(second_trafo_step)then
            !get the order sigma [gamma i j alpha]

            pos2 = 1+full1*full2*nor
            if(case_sel==3.or.case_sel==4)then
               l1=fa
               l2=fg
            else
               l1=fa
               l2=fg+goffs+tlen
            endif
            call mat_transpose(full1T,full2T*nor,1.0E0_realk,w0(pos2:full1T*full2T*nor+pos2-1),0.0E0_realk,w2)
            !transform gamma -> l
            call dgemm('t','n',no2,nor*full1T,full2T,1.0E0_realk,xocc(l1),nb,w2,full2T,0.0E0_realk,w3,no2)

            !transform alpha -> k, order is now sigma[k l i j]
            if( rest_sio4 )then
               call dgemm('t','t',no2,no2*nor,full1T,1.0E0_realk,xocc(l2),nb,w3,nor*no2,1.0E0_realk,sio4%elm1,no2)
            else
               select case(s)
               case(1)
#ifdef DIL_ACTIVE
                if(DIL_DEBUG) then !`DIL: Tensor contraction 11
                 write(DIL_CONS_OUT,'("#DEBUG(DIL): Process ",i6,"[",i6,"] starting tensor contraction 11:")')&
                  &infpar%lg_mynum,infpar%mynum
                endif
                tcs='D(k,l,i)+=L(u,k)*R(l,i,u)'
                call dil_clean_tens_contr(tch0)
                call dil_set_tens_contr_args(tch0,'d',errc,tens_distr=sio4)
                if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC11: DA: ',infpar%lg_mynum,errc
                if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC11: Destination arg set failed!',-1)
                tens_rank=2; tens_dims(1:tens_rank)=(/nb,no2/)
                call dil_set_tens_contr_args(tch0,'l',errc,tens_rank,tens_dims,xocc)
                if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC11: LA: ',infpar%lg_mynum,errc
                if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC11: Left arg set failed!',-1)
                tens_rank=3; tens_dims(1:tens_rank)=(/no2,nor,full1T/); tens_bases(1:tens_rank)=(/0,0,l2-1/)
                call dil_set_tens_contr_args(tch0,'r',errc,tens_rank,tens_dims,w3,tens_bases)
                if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC11: RA: ',infpar%lg_mynum,errc
                if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC11: Right arg set failed!',-1)
                call dil_set_tens_contr_spec(tch0,tcs,errc,&
                     &ldims=(/int(full1T,INTD),int(no2,INTD)/),lbase=(/int(l2-1,INTD),0_INTD/),&
                     &rdims=(/int(no2,INTD),int(nor,INTD),int(full1T,INTD)/),rbase=(/0_INTD,0_INTD,int(l2-1,INTD)/),&
                     &alpha=1E0_realk)
                if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC11: CC: ',infpar%lg_mynum,errc
                if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC11: Contr spec set failed!',-1)
                dil_mem=dil_get_min_buf_size(tch0,errc)
                if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC11: BS: ',infpar%lg_mynum,errc,dil_mem
                if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC11: Buf size set failed!',-1)
                dil_mem=dil_buf_size*8_INTL; call dil_prepare_buffer(tch0,dil_mem,errc,dil_buffer)
                if(errc.ne.0) call lsquit('ERROR(ccsd_residual_integral_driven): TC11: Buf alloc failed!',-1)
                call dil_tensor_contract(tch0,DIL_TC_EACH,dil_mem,errc,locked=lock_outside)
                if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC11: TC: ',infpar%lg_mynum,errc
                if(errc.ne.0) call lsquit('ERROR(combine_and_transform_sigma): TC11: Tens contr failed!',-1)
#else
                call lsquit('ERROR(combine_and_transform_sigma): This part (10) of Scheme 1 requires DIL backend!',-1)
#endif
               case(2)
                  call dgemm('t','t',no2,no2*nor,full1T,1.0E0_realk,xocc(l2),nb,w3,nor*no2,0.0E0_realk,w2,no2)
                  call squareup_block_triangular_squarematrix(w2,no,no,do_block_transpose = .true.)
#ifdef VAR_MPI
                  call time_start_phase(PHASE_COMM, at=twork)
                  if( lock_outside .and..not. alloc_in_dummy )call tensor_lock_wins(sio4,'s',mode)
                  call tensor_add(sio4,1.0E0_realk,w2,wrk=w3,iwrk=s3)
                  if(alloc_in_dummy)then
                     call lsmpi_win_flush(sio4%wi(1),local=.true.)
                  else
                     if( lock_outside )call tensor_unlock_wins(sio4,.true.)
                  endif
                  call time_start_phase(PHASE_WORK, at=tcomm)
#endif
               case(3,4)
                  call dgemm('t','t',no2,no2*nor,full1T,1.0E0_realk,xocc(l2),nb,w3,nor*no2,0.0E0_realk,w2,no2)
#ifdef VAR_PTR_RESHAPE
                  t1(1:no2,1:no2,1:nor)     => w2
                  !this is a workaround for PGI 15.5 and 15.7
                  h => sio4%elm1
                  h1(1:no,1:no,1:no2,1:no2) => h
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                  call c_f_pointer(c_loc(w2(1)),t1,[no2,no2,nor])
                  call c_f_pointer(c_loc(sio4%elm1(1)),h1,[no,no,no2,no2])
#else
                  call lsquit('ERROR(combine_and_transform_sigma): unable to reshape pointers!',-1)
#endif
                  do j=no,1,-1
                     do i=j,1,-1
                        call array_reorder_2d(1.0E0_realk,t1(:,:,i+j*(j-1)/2),no2,no2,[2,1],1.0E0_realk,h1(i,j,:,:))
                        if(i /= j)then
                           h1(j,i,:,:) = h1(j,i,:,:) + t1(:,:,i+j*(j-1)/2)
                        endif
                     enddo
                  enddo
               end select
            endif
         endif
      endif


   end subroutine combine_and_transform_sigma


   !> \brief Construct integral matrix which uses symmetries in the virtual
   !> indices
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine get_I_cged(Int_out,Int_in,m,nv,acc_h)
      implicit none
      !> leading dimension m and virtual dimension
      integer,intent(in)::m,nv
      !> integral output with indices reduced to c>=d
      real(realk),intent(inout)::Int_out(m*(nv*(nv+1))/2)
      !>full integral input m*c*d
      real(realk),intent(in) :: Int_in(m*nv*nv)
      integer ::d,pos,pos2,a,b,c,cged,dc
      logical :: doit
#ifdef VAR_OPENACC
      integer(kind=acc_handle_kind),intent(in),optional :: acc_h
      integer(kind=acc_handle_kind) :: ac_
#ifdef VAR_PGF90
      logical :: stuff_here
      stuff_here = acc_is_present(Int_out,m*(nv*(nv+1))/2*8) .and. acc_is_present(Int_in,m*nv*nv*8)
#else
      logical, parameter :: stuff_here = .true.
#endif
#else
      integer,intent(in),optional :: acc_h
      integer :: ac_
      logical, parameter :: stuff_here = .false.
#endif
      ac_ = 0
#ifdef VAR_OPENACC
      ac_ = acc_async_sync
#endif
      if(present(acc_h)) ac_ = acc_h

      if(stuff_here)then
         !$acc parallel loop present(Int_out,Int_in) default(none) copyin(nv,m) private(pos,pos2,d,dc) async(ac_)
         do d=1,nv

            !calculate target position in array Int_out
            pos =1
            do dc=1,d-1
               pos=pos+m*(nv-dc+1)
            enddo

            !calculate origin position in array Int_in
            pos2=1+(d-1)*m+(d-1)*nv*m

            !do the copy
            Int_out(pos:pos+m*(nv-d+1)-1) = Int_in(pos2:pos2+m*(nv-d+1)-1)

         enddo
         !$acc end parallel loop
      else
         do d=1,nv

            !calculate target position in array Int_out
            pos =1
            do dc=1,d-1
               pos=pos+m*(nv-dc+1)
            enddo

            !calculate origin position in array Int_in
            pos2=1+(d-1)*m+(d-1)*nv*m

            !do the copy
            Int_out(pos:pos+m*(nv-d+1)-1) = Int_in(pos2:pos2+m*(nv-d+1)-1)

         enddo
      endif

   end subroutine get_I_cged

   !Even though the following two functions are only needed inside
   !get_I_plusminus_le, they were moved here since PGI decided that internal
   !procedures may not be targets of function pointers
   function a_plus_b(a,b) result(c)
      implicit none
      real(realk), intent(in)  :: a, b
      real(realk) :: c
      c = a + b
   end function a_plus_b
   function a_minus_b(a,b) result(c)
      implicit none
      real(realk), intent(in)  :: a, b
      real(realk) :: c
      c = a - b
   end function a_minus_b

   !> \brief Construct symmetric and antisymmentric combinations of an itegral matrix 
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine get_I_plusminus_le(w2,op,fa,fg,la,lg,nb,tlen,tred,goffs,s2,faleg,laleg)
      implicit none
      integer(kind=8), intent(in) :: s2
      !> blank workspace
      real(realk),intent(inout),target :: w2(s2)
      !> integer specifying the first element in alpha and gamma batch
      integer,intent(in) :: fa,fg
      !> integer specifying the length of alpha and gamma batches
      integer,intent(in) :: la,lg
      !> number of upper triangular elements in the current batch
      integer,intent(in) :: tred
      !> number of elements in triangular part
      integer,intent(in) :: tlen
      !> number of ao indices
      integer,intent(in) :: nb
      !> offset in the gamma batch
      integer,intent(in) :: goffs
      !> character specifying which operation to use
      character, intent(in) :: op
      !> integers to specify the batch of a<=g that has to be calculated, faleg
      !is the first element of the batch, laleg is the length of the bathc
      integer,intent(in) :: faleg,laleg
      
      integer :: i,alpha,beta,gamm,delta,cagc,actr,bs,bctr
      integer :: alpha_b,beta_b,gamm_b,delta_b,elements,aleg
      integer :: globg, globa, loca,eldiag,elsqre,nbnb,nrnb
      real(realk) ::chk,chk2,el
      real(realk),pointer :: trick(:,:,:)
      logical :: modb
      procedure(ab_eq_c), pointer :: a_op_b => null()

      select case(op)
      case ('+')
         a_op_b => a_plus_b
      case ('-')
         a_op_b => a_minus_b
      case default
         call lsquit("ERROR(get_I_plusminus_le): wrong op on input",-1)
      end select

      bs=int(sqrt(((8.0E6_realk)/1.6E1_realk)))
      !bs=5
      !print *,"block size",bs,(bs*bs*8)/1024.0E0_realk
      nbnb=(nb/bs)*bs
      modb=(mod(nb,bs)>0)
      bctr = bs-1

#ifdef VAR_PTR_RESHAPE
      trick(1:nb,1:nb,1:laleg) => w2
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
      call c_f_pointer(c_loc(w2(1)),trick,[nb,nb,laleg])
#else
      call lsquit('ERROR(get_I_plusminus_le): unable to reshape pointers!',-1)
#endif
      aleg=0
      actr=0

      do gamm=0,lg-1
         do alpha=0,la-1
            if(fa+alpha<=fg+gamm.and.aleg+1>=faleg.and.aleg+1<faleg+laleg)then
               !aleg = (alpha+(gamm*(gamm+1))/2) 
               eldiag = actr*nb*nb
               elsqre = alpha*nb*nb+gamm*nb*nb*la
               !print *,alpha,gamm,1+eldiag,nb*nb+eldiag,1+elsqre,nb*nb+elsqre,aleg,laleg,nb*nb
               if(fa+alpha==fg+gamm) call dscal(nb*nb,0.5E0_realk,w2(1+elsqre),1)
               call dcopy(nb*nb,w2(1+elsqre),1,w2(1+eldiag),1)
               !$OMP PARALLEL PRIVATE(el,delta_b,beta_b,beta,delta)&
               !$OMP SHARED(bs,bctr,trick,nb,actr,nbnb,modb,a_op_b)&
               !$OMP DEFAULT(NONE)
               if(nbnb>0)then
                  !$OMP DO
                  do delta_b=1,nbnb,bs
                     do beta_b=delta_b+bs,nbnb,bs
                        do delta=0,bctr
                           do beta=0,bctr
                              el=trick(beta+beta_b,delta_b+delta,actr+1)
                              trick(beta+beta_b,delta_b+delta,actr+1)=&
                                 &a_op_b(trick(delta_b+delta,beta+beta_b,actr+1),trick(beta+beta_b,delta_b+delta,actr+1))
                              trick(delta_b+delta,beta+beta_b,actr+1)=a_op_b(el,trick(delta_b+delta,beta+beta_b,actr+1))
                           enddo
                        enddo
                     enddo
                  enddo
                  !$OMP END DO NOWAIT
               endif
               if(nbnb>0.and.modb)then
                  !$OMP DO
                  do delta_b=1,nbnb,bs
                     do delta=0,bctr
                        do beta=nbnb+1,nb
                           el=trick(beta,delta_b+delta,actr+1)
                           trick(beta,delta_b+delta,actr+1)=&
                              &a_op_b(trick(delta_b+delta,beta,actr+1),trick(beta,delta_b+delta,actr+1))
                           trick(delta_b+delta,beta,actr+1)=a_op_b(el,trick(delta_b+delta,beta,actr+1))
                        enddo
                     enddo
                  enddo
                  !$OMP END DO NOWAIT
               endif
               if(nbnb>0)then
                  !$OMP DO
                  do delta_b=1,nbnb,bs
                     do delta=0,bctr
                        do beta=delta+1,bctr
                           el=trick(beta+delta_b,delta_b+delta,actr+1)
                           trick(beta+delta_b,delta_b+delta,actr+1)=&
                              &a_op_b(trick(delta_b+delta,beta+delta_b,actr+1),trick(beta+delta_b,delta_b+delta,actr+1))
                           trick(delta_b+delta,beta+delta_b,actr+1)=a_op_b(el,trick(delta_b+delta,beta+delta_b,actr+1))
                        enddo
                        trick(delta+delta_b,delta_b+delta,actr+1)=&
                           &a_op_b(trick(delta_b+delta,delta+delta_b,actr+1),trick(delta+delta_b,delta_b+delta,actr+1))
                     enddo
                  enddo
                  !$OMP END DO NOWAIT
               endif
               !$OMP END PARALLEL 
               if(modb)then
                  do delta=nbnb+1,nb
                     do beta=delta+1,nb
                        el=trick(beta,delta,actr+1)
                        trick(beta,delta,actr+1)=a_op_b(trick(delta,beta,actr+1),trick(beta,delta,actr+1))
                        trick(delta,beta,actr+1)=a_op_b(el,trick(delta,beta,actr+1))
                     enddo
                     trick(delta,delta,actr+1)=a_op_b(trick(delta,delta,actr+1),trick(delta,delta,actr+1))
                  enddo
               endif
               actr=actr+1
            endif

            !count aleg
            if(fa+alpha<=fg+gamm)then
               aleg=aleg+1
            endif
         enddo
      enddo

      trick => null()


   end subroutine get_I_plusminus_le


   !> \brief subroutine to add contributions to the sio4 matrix which enters the
   !B2.2 term in the "non"-parallel region
   !> \author Patrick Ettenhuber
   !> \Date September 2012
   subroutine add_int_to_sio4(w0,w2,w3,nor,no,nv,nb,fa,fg,la,lg,xo,sio4,act_no)
      implicit none
      !> workspace containing the paritially transformed integrals ordered as I(i
      !gamma alpha j)
      real(realk),pointer :: w0(:)
      !> arbitrary workspace of correct size
      real(realk),pointer :: w2(:),w3(:)
      !> number of occupied, virutal and ao indices
      integer, intent(in) :: nor,no,nv,nb
      !> first alpha and first gamma indices of the current loop
      integer, intent(in) :: fa,fg
      !> lengths of the alpha ang gamma batches in the currnet loop
      integer, intent(in) :: la,lg
      !> transformation matrix for t1 transformed integrals "Lambda p"
      real(realk),intent(in) :: xo(nb*no)
      integer, intent(in), optional :: act_no
      !> sio4 storage space to update during the batched loops
      real(realk),pointer :: sio4(:)
      integer :: pos,i,j,no2
      integer(kind=8) :: pos1, pos2

      no2 = no
      if(present(act_no))no2 = act_no

      ! (w3):I[alpha gamma i j] <- (w0):I[i gamma alpha j]
      call array_reorder_4d(1.0E0_realk,w0,no,lg,la,no,[2,3,1,4],0.0E0_realk,w2)
      ! (w2):I[alpha gamma i <= j] <- (w3):I[alpha gamma i j]
      do j=1,no
         do i=1,j
            pos1=1_long+((i+j*(j-1)/2)-1)*la*(lg*i8)
            pos2=1_long+(i-1)*la*lg+(j-1)*la*lg*(no*i8)
            !call dcopy(la*lg,w2(1+(i-1)*la*lg+(j-1)*la*lg*no),1,w3(pos),1)
            w3(pos1:pos1+la*lg-1) = w2(pos2:pos2+la*lg-1)
         enddo
      enddo
      ! (w3):I[ alpha i <= j gamma] <- (w2):I[alpha gamma i <= j]
      call array_reorder_3d(1.0E0_realk,w3,lg,la,nor,[2,3,1],0.0E0_realk,w2)
      ! (w2):I[ l i <= j gamma] <- (w3):Lambda^p [alpha l ]^T I[alpha i <= j gamma]
      call dgemm('t','n',no2,nor*lg,la,1.0E0_realk,xo(fa),nb,w2,la,0.0E0_realk,w3,no2)
      ! (sio4):I[ k l i <= j] <-+ (w2):Lambda^p [gamma k ]^T I[l i <= j , gamma]^T
      call dgemm('t','t',no2,nor*no2,lg,1.0E0_realk,xo(fg),nb,w3,nor*no2,1.0E0_realk,sio4,no2)

   end subroutine add_int_to_sio4

   !> \brief Get the b2.2 contribution constructed in the kobayashi scheme after
   !the loop to avoid steep scaling ste  !> \author Patrick Ettenhuber
   !> \date December 2012
   subroutine get_B22_contrib_mo(sio4,t2,w1,w2,no,nv,om2,s,lock_outside,tw,tc,no_par,order,tmp_tens)
      implicit none
      !> the sio4 matrix from the kobayashi terms on input
      type(tensor), intent(in) :: sio4
      !> amplitudes
      !real(realk), intent(in) :: t2(*)
      type(tensor), intent(inout) :: t2
      !> some workspave
      real(realk), intent(inout) :: w1(:)
      real(realk), pointer :: w2(:)
      !> number of occupied, virutal and ao indices
      integer, intent(in) :: no,nv
      !> residual to be updated
      !real(realk), intent(inout) :: om2(*)
      type(tensor), intent(inout) :: om2
      !> integer specifying the calc-scheme
      integer, intent(in) :: s
      logical, intent(in) :: lock_outside
      !> work and communication time in B2 term
      real(realk), intent(inout) :: tw,tc
      logical, intent(in),optional :: no_par
      integer, intent(in),optional :: order(4)
      type(tensor), intent(inout), optional:: tmp_tens !already contains some previous contributions
      integer :: nor
      integer :: ml,l,tl,fai,lai
      integer :: tri,fri
      integer(kind=ls_mpik) :: nod,me,nnod,massa,mode
      real(realk) :: nrm1,nrm2,nrm3,nrm4
      integer ::  mv((nv*nv)/2),st
      integer(kind=8) :: o2v2,pos1,pos2,i,j,pos
      logical :: traf,np
      integer :: o(4)
!{`DIL:
#ifdef DIL_ACTIVE
     character(256):: tcs
     type(dil_tens_contr_t):: tch
     integer(INTL):: dil_mem,l0
     integer(INTD):: i0,i1,i2,i3,errc,tens_rank,tens_dims(MAX_TENSOR_RANK),tens_bases(MAX_TENSOR_RANK)
     integer(INTD):: ddims(MAX_TENSOR_RANK),ldims(MAX_TENSOR_RANK),rdims(MAX_TENSOR_RANK)
     integer(INTD):: dbase(MAX_TENSOR_RANK),lbase(MAX_TENSOR_RANK),rbase(MAX_TENSOR_RANK)
#endif
     real(realk):: r0
!}

      call time_start_phase(PHASE_WORK)

      np = .false.
      if(present(no_par)) np = no_par
      o = [1,2,3,4]
      if(present(order)) o  = order

      if(s==1.and.(.not.present(tmp_tens))) call lsquit('ERROR(get_B22_contrib_mo): Scheme 1 requires <tmp_tens> argument!',-1)

      me    = 0
      massa = 0
      nnod  = 1
#ifdef VAR_MPI
      massa = infpar%master
      nnod  = infpar%lg_nodtot
      me    = infpar%lg_mynum
      mode  = int(MPI_MODE_NOCHECK,kind=ls_mpik)
#endif
      o2v2=(i8*no)*no*nv*nv

      !Setting transformation variables for each rank
      !**********************************************
      if(np)then
         fai = 1
         tl  = nv*nv
      else
         call mo_work_dist(nv*nv,fai,tl,traf)
      endif

      if(DECinfo%PL>3.and.me==0)then
         write(DECinfo%output,'("Trafolength in striped B2:",I5)')tl
      endif

      nor=no*(no+1)/2

      ! do contraction
      if((s==4.or.s==3).and.traf)then


         w1=0.0E0_realk
         call dgemm('n','n',tl,nor,no*no,0.5E0_realk,t2%elm1(fai),nv*nv,sio4%elm1,no*no,0.0E0_realk,w1(fai),nv*nv)
         
         call squareup_block_triangular_squarematrix(w1,nv,no,do_block_transpose = .true.)

         if(present(order))then
            call array_reorder_4d(1.0E0_realk,w1,nv,nv,no,no,order,1.0E0_realk,om2%elm1)
         else
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
            call assign_in_subblocks(om2%elm1,'+',w1,o2v2)
#else
            !$OMP WORKSHARE
            om2%elm1(1:o2v2) = om2%elm1(1:o2v2) + w1(1:o2v2)
            !$OMP END WORKSHARE
#endif
         endif


      else if(s==2.and.traf)then

#ifdef VAR_MPI
         o = [1,2,3,4]
         call tensor_contract(0.5E0_realk,t2,sio4,[3,4],[1,2],2,1.0E0_realk,om2,o,force_sync=.true.)
#endif

      else if(s==1)then

#ifdef DIL_ACTIVE
         if(DIL_DEBUG) then !`DIL: Tensor contraction 12
          write(DIL_CONS_OUT,'("#DEBUG(DIL): Process ",i6,"[",i6,"] starting tensor contraction 12:",3(1x,i7))')&
          &infpar%lg_mynum,infpar%mynum
         endif
         tcs='D(a,b,s)+=L(a,b,i,j)*R(i,j,s)'
         call dil_clean_tens_contr(tch)
         call dil_set_tens_contr_args(tch,'d',errc,tens_distr=tmp_tens)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC12: DA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_B22_contrib_mo): TC12: Destination arg set failed!',-1)
         call dil_set_tens_contr_args(tch,'l',errc,tens_distr=t2)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC12: LA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_B22_contrib_mo): TC12: Left arg set failed!',-1)
         call dil_set_tens_contr_args(tch,'r',errc,tens_distr=sio4)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC12: RA: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_B22_contrib_mo): TC12: Right arg set failed!',-1)
         call dil_set_tens_contr_spec(tch,tcs,errc,alpha=0.5E0_realk)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC12: CC: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_B22_contrib_mo): TC12: Contr spec set failed!',-1)
         dil_mem=dil_get_min_buf_size(tch,errc)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC12: BS: ',infpar%lg_mynum,errc,dil_mem
         if(errc.ne.0) call lsquit('ERROR(get_B22_contrib_mo): TC12: Buf size set failed!',-1)
         dil_mem=dil_buf_size*8_INTL; call dil_prepare_buffer(tch,dil_mem,errc,dil_buffer)
         if(errc.ne.0) call lsquit('ERROR(ccsd_residual_integral_driven): TC12: Buf alloc failed!',-1)
         call dil_tensor_contract(tch,DIL_TC_ALL,dil_mem,errc,locked=.false.)
         if(DIL_DEBUG) write(DIL_CONS_OUT,*)'#DIL: TC12: TC: ',infpar%lg_mynum,errc
         if(errc.ne.0) call lsquit('ERROR(get_B22_contrib_mo): TC12: Tens contr failed!',-1)
         call lsmpi_barrier(infpar%lg_comm)
         call dil_update_abij_with_abc(om2,tmp_tens,errc,locked=.false.) !omega2 is locked, tmp_tens not
         if(errc.ne.0) call lsquit('ERROR(get_B22_contrib_mo): TC12: Sym update failed!',-1)
         call lsmpi_barrier(infpar%lg_comm)
#endif

      else
         call lsquit('ERROR(get_B22_contrib_mo): Unknown IF-THEN branch!',-1)
      endif

      call time_start_phase(PHASE_COMM, at = tw )
   end subroutine get_B22_contrib_mo



   subroutine print_tensor_unfolding_with_labels(mat,d1,l1,m1,d2,l2,m2,label)
      implicit none
      integer, intent(in) :: m1,m2
      integer, intent(in) :: d1(m1),d2(m2)
      character,intent(in) :: l1(m1),l2(m2)
      real(realk), intent(in) :: mat(*)
      character*(*), intent(in) :: label
      integer :: i,a,b,id1(m1),id2(m2),cd1,cd2

      print *,label(1:len(label))
      write (*,'("(")',advance='no') 
      cd1 = 1
      do i=1,m1
         if(i>1)write (*,'(X)',advance='no')
         write (*,'(A)',advance='no') l1(i)
         cd1 = cd1 * d1(i)
      enddo
      write (*,'(")  /  (")',advance='no') 
      cd2 = 1
      do i=1,m2
         if(i>1)write (*,'(X)',advance='no')
         write (*,'(A)',advance='no') l2(i)
         cd2 = cd2 * d2(i)
      enddo
      write (*,'(")  ")',advance='no') 

      do b = 1, cd2
         call get_midx(b,id2,d2,m2)
         write (*,'(12X,"(")',advance='no')
         do i = 1, m2
            if(i>1)write (*,'("/")',advance='no')
            write (*,'(I4)',advance='no') id2(i)
         enddo
         write (*,'(")")',advance='no')
      enddo
      write(*,*)

      do a = 1, cd1
         call get_midx(a,id1,d1,m1)
         write (*,'(" (")',advance='no') 
         do i = 1, m1
            if(i>1)write (*,'("/")',advance='no')
            write (*,'(I4)',advance='no') id1(i)
         enddo
         write (*,'(") ")',advance='no') 
         do b = 1, cd2
            write (*,'(g18.5)',advance='no') mat(a+(b-1)*cd1)
            if(b<cd2)write (*,'(",",X)',advance='no')
         enddo
         if(a<cd1)then
            write (*,'(";")')
         else
            write (*,*)
         endif
      enddo



   end subroutine print_tensor_unfolding_with_labels

   subroutine squareup_block_triangular_squarematrix(BTM,blockidx,sidx,do_block_transpose)
      implicit none
      integer ,intent(in) :: blockidx,sidx
      real(realk), intent(inout) :: BTM(blockidx**2*sidx**2)
      logical, intent(in), optional :: do_block_transpose
      logical :: dbt
      integer :: blocksize
      integer :: mv((blockidx**2)/2),st,i,j,pos1,pos2

      dbt = .false.
      if(present(do_block_transpose))dbt=do_block_transpose

      blocksize = blockidx**2

      do j=sidx,1,-1
         do i=j,1,-1
            pos1=1+((i+j*(j-1)/2)-1)*blocksize
            pos2=1+(i-1)*blocksize+(j-1)*sidx*blocksize
            if(j/=1) BTM(pos2:pos2+blocksize-1) = BTM(pos1:pos1+blocksize-1)
         enddo
      enddo
      do j=sidx,1,-1
         do i=j,1,-1
            pos1=1+(i-1)*blocksize+(j-1)*sidx*blocksize
            pos2=1+(j-1)*blocksize+(i-1)*sidx*blocksize
            if(i/=j) BTM(pos2:pos2+blocksize-1) = BTM(pos1:pos1+blocksize-1)
         enddo
      enddo

      if(dbt)then
         do j=sidx,1,-1
            do i=j,1,-1
               pos1=1+(i-1)*blocksize+(j-1)*sidx*blocksize
               call alg513(BTM(pos1:blocksize+pos1-1),blockidx,blockidx,blocksize,mv,blocksize/2,st)
            enddo
         enddo
      endif


   end subroutine squareup_block_triangular_squarematrix


   !> \brief: determine i and j from ij
   !> \author: Janus Juul Eriksen
   !> \date: july 2013
   subroutine calc_i_leq_j(ij,full,i,j)

      implicit none

      !> composite ij index
      integer, intent(in) :: ij,full
      !> i and j
      integer, intent(inout) :: i,j
      !> integers
      integer :: gauss_sum,gauss_sum_old,series

      ! in a N x N lower triangular matrix, there is a total of (N**2 + N)/2 elements.
      ! in column 1, there are N rows, in column 2, there are (N-1) rows, ...,  in
      ! column (N-1), there are 2 rows, and in column N, there are 1 row.
      ! this is a gauss sum of 1 + 2 + 3 + ... + N-2 + N-1 + N
      ! for a given value of i, the value of ij can thus at max be (i**2 + i)/2 (gauss_sum).
      ! if gauss_sum for a given i (series) is smaller than ij, we loop.
      ! when gauss_sum is greater than ij, we use the value of i for the present loop cycle
      ! and calculate the value of j from the present ij and previous gauss_sum values.
      ! when gauss_sum is equal to ij, we are on the diagonal and i == j (== series).

      do series = 1,full

         gauss_sum = int((series**2 + series)/2)

         if (gauss_sum .lt. ij) then

            gauss_sum_old = gauss_sum

            cycle

         else if (gauss_sum .eq. ij) then

            j = series
            i = series

            exit

         else

            j = ij - gauss_sum_old
            i = series

            exit

         end if

      end do

   end subroutine calc_i_leq_j

   subroutine calc_i_geq_j(ij,full,i,j)

      implicit none

      !> composite ij index
      integer, intent(in) :: ij,full
      !> i and j
      integer, intent(inout) :: i,j
      !> integers
      integer :: igeqj,series,gauss_sum,gauss_sum_old

      do series = 1,full

         gauss_sum = series + (series-1) * full - ((series-1)*series)/2

         if (gauss_sum .lt. ij) then

            gauss_sum_old = gauss_sum

            cycle

         else if (gauss_sum .eq. ij) then

            j = series
            i = series

            exit

         else

            j = series - 1
            i = j + ij - gauss_sum_old

            exit

         end if

      end do

   end subroutine calc_i_geq_j

   function get_nbuffs_scheme_0() result (nbuffs)
      implicit none
      integer :: nbuffs
      nbuffs = 5
   end function get_nbuffs_scheme_0

   function get_split_scheme_0(full) result (split)
      implicit none
      integer,intent(in) :: full
      integer :: split
      integer :: nnod

      nnod = 0
#ifdef VAR_MPI
      nnod = infpar%lg_nodtot
#endif

      split = full/int(sqrt(float(nnod))+2)

   end function get_split_scheme_0

   subroutine simulate_intloop_and_get_worksize(maxint,nb,nbg,nba,bs,intspec,setting)
      implicit none
      integer(kind=long), intent(out) :: maxint
      integer,intent(in) :: nb,nbg,nba,bs
      type(lssetting),intent(inout) :: setting
      Character,intent(in) :: intspec(5)
      integer :: gb,ab,fg,lg,fa,la
      integer :: ntiles
      integer :: starts(4), tdim(4), dims(4), ntpm(4)
      integer :: ndimA,  ndimB,  ndimC,  ndimD
      integer :: ndimAs, ndimBs, ndimCs, ndimDs
      integer :: startA, startB, startC, startD
      integer :: nbtchsg, nbtchsa
      integer :: i,nba1,nbg1,as,gs

      nbtchsa = nb / nba
      if( mod( nb, nba ) > 0 ) nbtchsa = nbtchsa + 1
      nbtchsg = nb / nbg
      if( mod( nb, nbg ) > 0 ) nbtchsg = nbtchsg + 1

      nba1 = nba
      nbg1 = nbg

      maxint = 0
      !simulate loop
      do gb=1,nbtchsg

         !short hand notation
         fg = 1 + (gb-1)*nbg
         lg = nb - fg + 1
         if( lg >= nbg )then
            lg = nbg
         endif

         if( lg > bs )then
            gs = bs
         else
            gs = lg
         endif


         do ab=1,nbtchsa
            !short hand notation
            fa = 1 + (ab-1)*nba
            la = nb - fa + 1
            if( la >= nba )then
               la = nba
            endif

            if( la > bs )then
               as = bs
            else
               as = la
            endif

            dims = [nb,nb,la,lg]
            tdim = [bs,bs,as,gs]

            call tensor_get_ntpm(dims,tdim,4,ntpm,ntiles = ntiles)

            do i = 1, ntiles

               call get_midx(i,starts,ntpm,4)

               ndimA  = dims(1)
               ndimB  = dims(2)
               ndimC  = dims(3)
               ndimD  = dims(4)

               startA = 1  + (starts(1)-1)*tdim(1)
               startB = 1  + (starts(2)-1)*tdim(2)
               startC = fa + (starts(3)-1)*tdim(3)
               startD = fg + (starts(4)-1)*tdim(4)

               call II_GET_ERI_INTEGRALBLOCK_INQUIRE(DECinfo%output,DECinfo%output,setting,&
                  & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
                  & ndimAs,ndimBs,ndimCs,ndimDs,INTSPEC)

               maxint = max(maxint,(i8*ndimAs*ndimBs)*ndimCs*ndimDs)

            enddo
         enddo
      enddo
   end subroutine simulate_intloop_and_get_worksize


   subroutine successive_4ao_mo_trafo(ao,WXYZ,WW,w,XX,x,YY,y,ZZ,z,WRKWXYZ)
      implicit none
      integer, intent(in) :: ao,w,x,y,z
      real(realk), intent(inout) :: WXYZ(ao*ao*ao*ao),WRKWXYZ(w*ao*ao*ao)
      real(realk), intent(in) :: WW(ao,w),XX(ao,x),YY(ao,y),ZZ(ao,z)
      !WXYZ(ao,ao ao ao)^T WW(ao,w)   -> WRKWXYZ (ao ao ao,w)
      call dgemm('t','n',ao*ao*ao,w,ao,1.0E0_realk,WXYZ,ao,WW,ao,0.0E0_realk,WRKWXYZ,ao*ao*ao)
      ! WRKWXYZ(ao,ao ao w)^T XX(ao,x)   -> WXYZ (ao ao w, x)
      call dgemm('t','n',ao*ao*w,x,ao,1.0E0_realk,WRKWXYZ,ao,XX,ao,0.0E0_realk,WXYZ,ao*ao*w)
      ! WXYZ(ao, ao w x)^T YY(ao,y)   -> WRKYXYX (ao w x,y)
      call dgemm('t','n',ao*w*x,y,ao,1.0E0_realk,WXYZ,ao,YY,ao,0.0E0_realk,WRKWXYZ,ao*w*x)
      ! WRKWXYZ(ao, w x y)^T ZZ(ao,z)^T   -> WXYZ (wxyz)
      call dgemm('t','n',w*x*y,z,ao,1.0E0_realk,WRKWXYZ,ao,ZZ,ao,0.0E0_realk,WXYZ,w*x*y)
   end subroutine successive_4ao_mo_trafo

   subroutine solver_energy_full_mod(no,nv,nfrags,offset,t2,t1,integral,occ_orbitals, &
         & virt_orbitals,FragOccEner,FragVirtEner,tmp_fragener,rmax,DTOA,DTVA,atom)

      implicit none

      !> solver doubles amplitudes and VOVO integrals (ordered as (a,b,i,j))
      type(tensor), intent(inout) :: t2, integral
      !> solver singles amplitudes
      type(tensor), intent(inout) :: t1
      !> dimensions
      integer, intent(in) :: no, nv, nfrags, offset
      !> occupied orbital information
      type(decorbital), dimension(no+offset), intent(inout) :: occ_orbitals
      !> virtual orbital information
      type(decorbital), dimension(nv), intent(inout) :: virt_orbitals
      !> Fragment energies array:
      real(realk), dimension(nfrags,nfrags), intent(inout) :: FragOccEner
      real(realk), dimension(nfrags,nfrags), intent(inout) :: FragVirtEner
      real(realk), dimension(nfrags,nfrags), intent(inout) :: tmp_fragener
      integer, intent(in) :: atom
      real(realk), intent(in) :: rmax, DTOA(:,:),DTVA(:,:)
      !> integers
      integer :: i,j,a,b
      integer :: atomI,atomJ,atomA,atomB
      !> energy reals
      real(realk) :: energy_tmp_1, energy_tmp_2
      real(realk), pointer :: t1p(:,:), t2p(:,:,:,:), inp(:,:,:,:)

      ! Pointer to avoid OMP problems:
      inp => integral%elm4(:,:,:,:)
      t2p => t2%elm4(:,:,:,:)
      t1p => t1%elm2(:,:)

      ! Get occupied partitioning energy:
      if (.not.DECinfo%OnlyVirtPart) then
         tmp_fragener=0.0e0_realk
         energy_tmp_1=0.0e0_realk
         energy_tmp_2=0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:FragOccEner),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,occ_orbitals,offset,DECinfo,DTVA,atom,rmax)
         do j=1,no
            atomJ = occ_orbitals(j+offset)%CentralAtom
            do i=1,no
               atomI = occ_orbitals(i+offset)%CentralAtom
               if(atomI == atomJ .and. atomI == atom)then

                  do b=1,nv
                     if(DTVA(b,atom)<=rmax)then
                        do a=1,nv
                           if(DTVA(a,atom)<=rmax)then

                              energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                              if(DECinfo%use_singles)then
                                 energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                              else
                                 energy_tmp_2 = 0.0E0_realk
                              endif
                              FragOccEner(AtomI,AtomJ) = FragOccEner(AtomI,AtomJ) &
                                 & + energy_tmp_1 + energy_tmp_2

                           endif
                        end do
                     endif
                  end do
               endif

            end do
         end do
         !$END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:tmp_fragener),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,occ_orbitals,offset,DECinfo,DTVA,atom,rmax)
         do j=1,no
            atomJ = occ_orbitals(j+offset)%CentralAtom
            do i=1,no
               atomI = occ_orbitals(i+offset)%CentralAtom
               if(atomI == atomJ .and. atomI == atom)then

                  do b=1,nv
                     if(DTVA(b,atom)<=rmax)then
                        do a=1,nv
                           if(DTVA(a,atom)<=rmax)then

                              energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                              if(DECinfo%use_singles)then
                                 energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                              else
                                 energy_tmp_2 = 0.0E0_realk
                              endif
                              tmp_fragener(AtomI,AtomJ) = tmp_fragener(AtomI,AtomJ) &
                                 & + energy_tmp_1 + energy_tmp_2

                           endif
                        end do
                     endif
                  end do
               endif

            end do
         end do
         !$END PARALLEL DO

         FragOccEner(:,:) = 2.0E0_realk * FragOccEner(:,:) - tmp_fragener

         do AtomI=1,nfrags
            do AtomJ=AtomI+1,nfrags
               FragOccEner(AtomI,AtomJ) = FragOccEner(AtomI,AtomJ) + FragOccEner(AtomJ,AtomI)
               FragOccEner(AtomJ,AtomI) = FragOccEner(AtomI,AtomJ)
            end do
         end do

         ! reorder from (a,b,j,i) to (a,b,i,j)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

      end if
      if (.not.DECinfo%OnlyOccPart) then
         tmp_fragener=0.0e0_realk
         energy_tmp_1=0.0e0_realk
         energy_tmp_2=0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomA,j,atomB,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:FragVirtEner),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,virt_orbitals,offset,DECinfo,DTOA,atom,rmax)
         do b=1,nv
            atomB = virt_orbitals(b)%CentralAtom
            do a=1,nv
               atomA = virt_orbitals(a)%CentralAtom
               if(atomA == atomB .and. atomA == atom)then

                  do j=1,no
                     if(DTOA(j+offset,atom)<=rmax)then
                        do i=1,no
                           if(DTOA(i+offset,atom)<=rmax)then

                              energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                              if(DECinfo%use_singles)then
                                 energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                              else
                                 energy_tmp_2 = 0.0E0_realk
                              endif
                              FragVirtEner(atomA,atomB) = FragVirtEner(atomA,atomB) &
                                 & + energy_tmp_1 + energy_tmp_2
                           endif

                        end do
                     endif
                  end do

               endif
            end do
         end do
         !$END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomA,j,atomB,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:tmp_fragener),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,virt_orbitals,offset,DECinfo,DTOA,atom,rmax)
         do b=1,nv
            atomB = virt_orbitals(b)%CentralAtom
            do a=1,nv
               atomA = virt_orbitals(a)%CentralAtom

               if(atomA == atomB .and. atomA == atom)then
                  do j=1,no
                     if(DTOA(j+offset,atom)<=rmax)then
                        do i=1,no
                           if(DTOA(i+offset,atom)<=rmax)then

                              energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                              if(DECinfo%use_singles)then
                                 energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                              else
                                 energy_tmp_2 = 0.0E0_realk
                              endif
                              tmp_fragener(atomA,atomB) = tmp_fragener(atomA,AtomB) &
                                 & + energy_tmp_1 + energy_tmp_2

                           endif
                        end do
                     endif
                  end do
               endif

            end do
         end do
         !$END PARALLEL DO

         FragVirtEner(:,:) = 2.0E0_realk * FragVirtEner(:,:) - tmp_fragener

         do AtomI=1,nfrags
            do AtomJ=AtomI+1,nfrags
               FragVirtEner(AtomI,AtomJ) = FragVirtEner(AtomI,AtomJ) + FragVirtEner(AtomJ,AtomI)
               FragVirtEner(AtomJ,AtomI) = FragVirtEner(AtomI,AtomJ)
            end do
         end do

         ! reorder from (a,b,j,i) to (a,b,i,j) in case of later use
         call tensor_reorder(integral,[1,2,4,3])

      end if

      inp => null()
      t2p => null()
      t1p => null()

   end subroutine solver_energy_full_mod

   !> \brief: calculate atomic and pair fragment contributions to solver (CCSD, MP2 ...)
   !> correlation energy for full molecule calculation.
   !> \author: Janus Juul Eriksen
   !> \date: February 2013
   subroutine solver_energy_full(no,nv,nfrags,offset,t2,t1,integral,occ_orbitals, &
         & virt_orbitals,FragOccEner,FragVirtEner,tmp_fragener)

      implicit none

      !> solver doubles amplitudes and VOVO integrals (ordered as (a,b,i,j))
      type(tensor), intent(inout) :: t2, integral
      !> solver singles amplitudes
      type(tensor), intent(inout) :: t1
      !> dimensions
      integer, intent(in) :: no, nv, nfrags, offset
      !> occupied orbital information
      type(decorbital), dimension(no+offset), intent(inout) :: occ_orbitals
      !> virtual orbital information
      type(decorbital), dimension(nv), intent(inout) :: virt_orbitals
      !> Fragment energies array:
      real(realk), dimension(nfrags,nfrags), intent(inout) :: FragOccEner
      real(realk), dimension(nfrags,nfrags), intent(inout) :: FragVirtEner
      real(realk), dimension(nfrags,nfrags), intent(inout) :: tmp_fragener
      !> integers
      integer :: i,j,a,b,atomI,atomJ
      !> energy reals
      real(realk) :: energy_tmp_1, energy_tmp_2
      real(realk), pointer :: t1p(:,:), t2p(:,:,:,:), inp(:,:,:,:)

      ! Pointer to avoid OMP problems:
      inp => integral%elm4(:,:,:,:)
      t2p => t2%elm4(:,:,:,:)
      t1p => t1%elm2(:,:)

      ! Get occupied partitioning energy:
      if (.not.DECinfo%OnlyVirtPart) then
         tmp_fragener=0.0e0_realk
         energy_tmp_1=0.0e0_realk
         energy_tmp_2=0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:FragOccEner),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,occ_orbitals,offset,DECinfo)
         do j=1,no
            atomJ = occ_orbitals(j+offset)%CentralAtom
            do i=1,no
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nv
                  do a=1,nv

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     FragOccEner(AtomI,AtomJ) = FragOccEner(AtomI,AtomJ) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:tmp_fragener),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,occ_orbitals,offset,DECinfo)
         do j=1,no
            atomJ = occ_orbitals(j+offset)%CentralAtom
            do i=1,no
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nv
                  do a=1,nv

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     tmp_fragener(AtomI,AtomJ) = tmp_fragener(AtomI,AtomJ) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         FragOccEner(:,:) = 2.0E0_realk * FragOccEner(:,:) - tmp_fragener

         do AtomI=1,nfrags
            do AtomJ=AtomI+1,nfrags
               FragOccEner(AtomI,AtomJ) = FragOccEner(AtomI,AtomJ) + FragOccEner(AtomJ,AtomI)
               FragOccEner(AtomJ,AtomI) = FragOccEner(AtomI,AtomJ)
            end do
         end do

         ! reorder from (a,b,j,i) to (a,b,i,j)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

      end if
      if (.not.DECinfo%OnlyOccPart) then
         tmp_fragener=0.0e0_realk
         energy_tmp_1=0.0e0_realk
         energy_tmp_2=0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:FragVirtEner),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,virt_orbitals,DECinfo)
         do b=1,nv
            atomJ = virt_orbitals(b)%CentralAtom
            do a=1,nv
               atomI = virt_orbitals(a)%CentralAtom

               do j=1,no
                  do i=1,no

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     FragVirtEner(AtomI,AtomJ) = FragVirtEner(AtomI,AtomJ) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:tmp_fragener),&
         !$OMP SHARED(t2p,t1p,inp,no,nv,virt_orbitals,DECinfo)
         do b=1,nv
            atomJ = virt_orbitals(b)%CentralAtom
            do a=1,nv
               atomI = virt_orbitals(a)%CentralAtom

               do j=1,no
                  do i=1,no

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     tmp_fragener(AtomI,AtomJ) = tmp_fragener(AtomI,AtomJ) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         FragVirtEner(:,:) = 2.0E0_realk * FragVirtEner(:,:) - tmp_fragener

         do AtomI=1,nfrags
            do AtomJ=AtomI+1,nfrags
               FragVirtEner(AtomI,AtomJ) = FragVirtEner(AtomI,AtomJ) + FragVirtEner(AtomJ,AtomI)
               FragVirtEner(AtomJ,AtomI) = FragVirtEner(AtomI,AtomJ)
            end do
         end do

         ! reorder from (a,b,j,i) to (a,b,i,j) in case of later use
         call tensor_reorder(integral,[1,2,4,3])

      end if

      inp => null()
      t2p => null()
      t1p => null()

   end subroutine solver_energy_full


   !> Purpose: calculate atomic fragment contributions to CCSD-like
   !           correlation energy for full molecule calculation.
   !           works for MP2, RPA, CC2, CCD, CCSD...
   !
   !> Author:  Pablo Baudin (Based on Janus's routine)
   !> Date:    Feb. 2015
   subroutine solver_decnp_full(nocc,nvirt,nfrags,offset,t2,t1,integral,occ_orbitals,&
            & virt_orbitals,FragOccEner,FragVirtEner,tmp_fragener)

      implicit none

      !> ccsd doubles amplitudes and VOVO integrals (ordered as (a,b,i,j))
      type(tensor), intent(inout) :: t2, integral
      !> ccsd singles amplitudes
      type(tensor), intent(inout) :: t1
      !> dimensions
      integer, intent(in) :: nocc, nvirt, nfrags, offset
      !> occupied orbital information
      type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
      !> virtual orbital information
      type(decorbital), dimension(nvirt), intent(inout) :: virt_orbitals
      !> Fragment energies array:
      real(realk), dimension(nfrags,nfrags), intent(inout) :: FragOccEner
      real(realk), dimension(nfrags,nfrags), intent(inout) :: FragVirtEner
      real(realk), dimension(nfrags,nfrags), intent(inout) :: tmp_fragener
      !> integers
      integer :: i,j,a,b,atomI
      !> energy reals
      real(realk) :: energy_tmp_1, energy_tmp_2
      real(realk), pointer :: t1p(:,:), t2p(:,:,:,:), inp(:,:,:,:)

      ! Pointer to avoid OMP problems:
      inp => integral%elm4(:,:,:,:)
      t2p => t2%elm4(:,:,:,:)
      t1p => t1%elm2(:,:)

      ! Get occupied partitioning energy:
      if (.not.DECinfo%OnlyVirtPart) then
         tmp_fragener=0.0e0_realk
         energy_tmp_1=0.0e0_realk
         energy_tmp_2=0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:FragOccEner),&
         !$OMP SHARED(t2p,t1p,inp,nocc,nvirt,occ_orbitals,offset,DECinfo)
         do j=1,nocc
            do i=1,nocc
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nvirt
                  do a=1,nvirt

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     FragOccEner(AtomI,AtomI) = FragOccEner(AtomI,AtomI) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:tmp_fragener),&
         !$OMP SHARED(t2p,t1p,inp,nocc,nvirt,occ_orbitals,offset,DECinfo)
         do j=1,nocc
            do i=1,nocc
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nvirt
                  do a=1,nvirt

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     tmp_fragener(AtomI,AtomI) = tmp_fragener(AtomI,AtomI) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! get total fourth--order energy contribution
         FragOccEner = 2.0E0_realk * FragOccEner - tmp_fragener

         ! reorder from (a,b,j,i) to (a,b,i,j)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)
      end if
      if (.not.DECinfo%OnlyOccPart) then
         tmp_fragener=0.0e0_realk
         energy_tmp_1=0.0e0_realk
         energy_tmp_2=0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:FragVirtEner),&
         !$OMP SHARED(t2p,t1p,inp,nocc,nvirt,virt_orbitals,DECinfo)
         do b=1,nvirt
            do a=1,nvirt
               atomI = virt_orbitals(a)%CentralAtom

               do j=1,nocc
                  do i=1,nocc

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     FragVirtEner(AtomI,AtomI) = FragVirtEner(AtomI,AtomI) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(integral,[1,2,4,3])
         inp => integral%elm4(:,:,:,:)

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,a,b,energy_tmp_1,energy_tmp_2),&
         !$OMP REDUCTION(+:tmp_fragener),&
         !$OMP SHARED(t2p,t1p,inp,nocc,nvirt,virt_orbitals,DECinfo)
         do b=1,nvirt
            do a=1,nvirt
               atomI = virt_orbitals(a)%CentralAtom

               do j=1,nocc
                  do i=1,nocc

                     energy_tmp_1 = t2p(a,b,i,j) * inp(a,b,i,j)
                     if(DECinfo%use_singles)then
                        energy_tmp_2 = t1p(a,i) * t1p(b,j) * inp(a,b,i,j)
                     else
                        energy_tmp_2 = 0.0E0_realk
                     endif
                     tmp_fragener(AtomI,AtomI) = tmp_fragener(AtomI,AtomI) &
                        & + energy_tmp_1 + energy_tmp_2

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         FragVirtEner(:,:) = 2.0E0_realk * FragVirtEner(:,:) - tmp_fragener

         ! reorder from (a,b,j,i) to (a,b,i,j) in case of later use
         call tensor_reorder(integral,[1,2,4,3])

      end if

      inp => null()
      t2p => null()
      t1p => null()

   end subroutine solver_decnp_full


   !> \brief: calculate E[4] contribution to ccsd(t) energy correction for full molecule.
   !> \author: Janus Juul Eriksen
   !> \date: February 2013
   subroutine ccsdpt_energy_e4_full(nocc,nvirt,nfrags,offset,ccsd_doubles,ccsdpt_doubles, &
         & occ_orbitals,virt_orbitals,eccsdpt_cou_occ,eccsdpt_cou_virt,eccsdpt_exc,ccsdpt_e4)

      implicit none

      !> ccsd and ccsd(t) doubles amplitudes
      type(tensor), intent(inout) :: ccsd_doubles, ccsdpt_doubles
      !> dimensions
      integer, intent(in) :: nocc, nvirt, nfrags, offset
      !> occupied orbital information
      type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
      type(decorbital), dimension(nvirt), intent(inout) :: virt_orbitals
      !> etot
      real(realk), intent(inout) :: ccsdpt_e4
      real(realk), dimension(nfrags,nfrags), intent(inout) :: eccsdpt_cou_occ
      real(realk), dimension(nfrags,nfrags), intent(inout) :: eccsdpt_cou_virt
      real(realk), dimension(nfrags,nfrags), intent(inout) :: eccsdpt_exc
      !> integers
      integer :: i,j,a,b,atomI,atomJ
      !> energy reals
      real(realk) :: energy_tmp, energy_res_cou, energy_res_exc

      ! *************************************************************
      ! ************** do energy for full molecule ******************
      ! *************************************************************

      ! ***********************
      !   do E[4] energy part
      ! ***********************

      ! ***note: we only run over nval (which might be equal to nocc_tot if frozencore = .false.)
      ! so we only assign orbitals for the space in which the core orbitals (the offset) are omited

      if (.not.DECinfo%OnlyVirtPart) then
         eccsdpt_exc = 0.0_realk
         energy_res_cou = 0.0E0_realk
         energy_res_exc = 0.0E0_realk
         ccsdpt_e4 = 0.0E0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_cou,eccsdpt_cou_occ), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,occ_orbitals,offset)
         do j=1,nocc
            atomJ = occ_orbitals(j+offset)%CentralAtom
            do i=1,nocc
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nvirt
                  do a=1,nvirt

                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_cou_occ(AtomI,AtomJ) = eccsdpt_cou_occ(AtomI,AtomJ) + energy_tmp
                     energy_res_cou = energy_res_cou + energy_tmp

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(ccsd_doubles,[1,2,4,3])

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_exc,eccsdpt_exc), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,occ_orbitals,offset)
         do j=1,nocc
            atomJ = occ_orbitals(j+offset)%CentralAtom
            do i=1,nocc
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nvirt
                  do a=1,nvirt

                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_exc(AtomI,AtomJ) = eccsdpt_exc(AtomI,AtomJ) + energy_tmp
                     energy_res_exc = energy_res_exc + energy_tmp

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! get total fourth--order energy contribution
         eccsdpt_cou_occ(:,:) = 4.0E0_realk * eccsdpt_cou_occ(:,:) &
            & - 2.0E0_realk * eccsdpt_exc
         ccsdpt_e4 = 4.0E0_realk * energy_res_cou - 2.0E0_realk * energy_res_exc

         ! for the e4 pair fragment energy matrix,
         ! we put the pair energy Delta E_IJ into both entry (I,J) and (J,I)

         do AtomJ=1,nfrags
            do AtomI=AtomJ+1,nfrags
               eccsdpt_cou_occ(AtomI,AtomJ) = eccsdpt_cou_occ(AtomI,AtomJ) &
                  & + eccsdpt_cou_occ(AtomJ,AtomI)
               eccsdpt_cou_occ(AtomJ,AtomI) =  eccsdpt_cou_occ(AtomI,AtomJ)
            end do
         end do

         ! reorder from (a,b,j,i) to (a,b,i,j) in case of later use
         call tensor_reorder(ccsd_doubles,[1,2,4,3])

      end if
      if (.not.DECinfo%OnlyOccPart) then
         eccsdpt_exc = 0.0_realk
         energy_res_cou = 0.0E0_realk
         energy_res_exc = 0.0E0_realk
         ccsdpt_e4 = 0.0E0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_cou,eccsdpt_cou_virt), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,virt_orbitals)
         do b=1,nvirt
            atomJ = virt_orbitals(b)%CentralAtom
            do a=1,nvirt
               atomI = virt_orbitals(a)%CentralAtom

               do j=1,nocc
                  do i=1,nocc

                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_cou_virt(AtomI,AtomJ) = eccsdpt_cou_virt(AtomI,AtomJ) + energy_tmp
                     energy_res_cou = energy_res_cou + energy_tmp

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(ccsd_doubles,[1,2,4,3])

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_exc,eccsdpt_exc), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,virt_orbitals)
         do b=1,nvirt
            atomJ = virt_orbitals(b)%CentralAtom
            do a=1,nvirt
               atomI = virt_orbitals(a)%CentralAtom

               do j=1,nocc
                  do i=1,nocc

                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_exc(AtomI,AtomJ) = eccsdpt_exc(AtomI,AtomJ) + energy_tmp
                     energy_res_exc = energy_res_exc + energy_tmp

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! get total fourth--order energy contribution
         eccsdpt_cou_virt(:,:) = 4.0E0_realk * eccsdpt_cou_virt(:,:) &
            & - 2.0E0_realk * eccsdpt_exc
         ccsdpt_e4 = 4.0E0_realk * energy_res_cou - 2.0E0_realk * energy_res_exc

         ! for the e4 pair fragment energy matrix,
         ! we put the pair energy Delta E_IJ into both entry (I,J) and (J,I)

         do AtomJ=1,nfrags
            do AtomI=AtomJ+1,nfrags
               eccsdpt_cou_virt(AtomI,AtomJ) = eccsdpt_cou_virt(AtomI,AtomJ) &
                  & + eccsdpt_cou_virt(AtomJ,AtomI)
               eccsdpt_cou_virt(AtomJ,AtomI) =  eccsdpt_cou_virt(AtomI,AtomJ)
            end do
         end do

         ! reorder from (a,b,j,i) to (a,b,i,j) in case of later use
         call tensor_reorder(ccsd_doubles,[1,2,4,3])
      end if



      ! ******************************************************************
      ! ************** done w/ energy for full molecule ******************
      ! ******************************************************************

   end subroutine ccsdpt_energy_e4_full


   !> \brief: calculate E[5] contribution to ccsd(t) energy correction for full molecule.
   !> \author: Janus Juul Eriksen
   !> \date: February 2013
   subroutine ccsdpt_energy_e5_full(nocc,nvirt,nfrags,offset,ccsd_singles,ccsdpt_singles,&
         & occ_orbitals,unocc_orbitals,e5_occ,e5_virt,ccsdpt_e5)

      implicit none

      !> ccsd and ccsd(t) singles amplitudes
      type(tensor), intent(inout) :: ccsd_singles, ccsdpt_singles
      !> dimensions
      integer, intent(in) :: nocc, nvirt, nfrags, offset
      !> occupied orbital information
      type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
      !> virtual orbital information
      type(decorbital), dimension(nvirt), intent(inout) :: unocc_orbitals
      !> etot
      real(realk), intent(inout) :: ccsdpt_e5
      real(realk), dimension(nfrags,nfrags), intent(inout) :: e5_occ
      real(realk), dimension(nfrags,nfrags), intent(inout) :: e5_virt
      !> integers
      integer :: i,a,AtomI,AtomA
      !> tmp energy real
      real(realk) :: energy_tmp

      ! ***********************
      !   do E[5] energy part
      ! ***********************

      ccsdpt_e5 = 0.0_realk
      energy_tmp = 0.0e0_realk
      !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,a,energy_tmp,AtomI,AtomA),&
      !$OMP SHARED(ccsd_singles,ccsdpt_singles,nocc,nvirt,offset,occ_orbitals,unocc_orbitals),&
      !$OMP REDUCTION(+:ccsdpt_e5),REDUCTION(+:e5_occ)
      do i=1,nocc
         AtomI = occ_orbitals(i+offset)%secondaryatom
         do a=1,nvirt
            AtomA = unocc_orbitals(a)%secondaryatom

            energy_tmp = ccsd_singles%elm2(a,i) * ccsdpt_singles%elm2(a,i)
            e5_occ(AtomA,AtomI) = e5_occ(AtomA,AtomI) + energy_tmp
            ! Important to update both (AtomI,AtomA) and (AtomA,AtomI) 
            e5_occ(AtomI,AtomA) = e5_occ(AtomA,AtomI)
            ccsdpt_e5 = ccsdpt_e5 + energy_tmp

         end do
      end do
      !$OMP END PARALLEL DO

      ! get total fifth-order energy correction
      e5_occ(:,:) = 2.0E0_realk * e5_occ(:,:)
      ccsdpt_e5 = 2.0E0_realk * ccsdpt_e5

      ! virtual partioning is the same as occupied for [5]:
      e5_virt(:,:) = e5_occ(:,:)


      ! ******************************
      !   done with E[5] energy part
      ! ******************************

   end subroutine ccsdpt_energy_e5_full


   !> \brief: calculate E[4] contribution to ccsd(t) energy correction for full molecule.
   !> \author: Pablo Baudin (from Janus Juul Eriksen)
   !> \date: Mar 2015
   subroutine ccsdpt_decnp_e4_full(nocc,nvirt,nfrags,offset,ccsd_doubles,ccsdpt_doubles, &
         & occ_orbitals,virt_orbitals,eccsdpt_cou_occ,eccsdpt_cou_virt,eccsdpt_exc,ccsdpt_e4)

      implicit none

      !> ccsd and ccsd(t) doubles amplitudes
      type(tensor), intent(inout) :: ccsd_doubles, ccsdpt_doubles
      !> dimensions
      integer, intent(in) :: nocc, nvirt, nfrags, offset
      !> occupied orbital information
      type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
      type(decorbital), dimension(nvirt), intent(inout) :: virt_orbitals
      !> etot
      real(realk), intent(inout) :: ccsdpt_e4
      real(realk), dimension(nfrags,nfrags), intent(inout) :: eccsdpt_cou_occ
      real(realk), dimension(nfrags,nfrags), intent(inout) :: eccsdpt_cou_virt
      real(realk), dimension(nfrags,nfrags), intent(inout) :: eccsdpt_exc
      !> integers
      integer :: i,j,a,b,atomI,atomJ
      !> energy reals
      real(realk) :: energy_tmp, energy_res_cou, energy_res_exc

      ! *************************************************************
      ! ************** do energy for full molecule ******************
      ! *************************************************************

      ! ***********************
      !   do E[4] energy part
      ! ***********************

      ! ***note: we only run over nval (which might be equal to nocc_tot if frozencore = .false.)
      ! so we only assign orbitals for the space in which the core orbitals (the offset) are omited

      if (.not.DECinfo%OnlyVirtPart) then
         eccsdpt_exc = 0.0_realk
         energy_res_cou = 0.0E0_realk
         energy_res_exc = 0.0E0_realk
         ccsdpt_e4 = 0.0E0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_cou,eccsdpt_cou_occ), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,occ_orbitals,offset)
         do j=1,nocc
            do i=1,nocc
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nvirt
                  do a=1,nvirt

                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_cou_occ(AtomI,AtomI) = eccsdpt_cou_occ(AtomI,AtomI) + energy_tmp
                     energy_res_cou = energy_res_cou + energy_tmp

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(ccsd_doubles,[1,2,4,3])

         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_exc,eccsdpt_exc), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,occ_orbitals,offset)
         do j=1,nocc
            do i=1,nocc
               atomI = occ_orbitals(i+offset)%CentralAtom

               do b=1,nvirt
                  do a=1,nvirt

                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_exc(AtomI,AtomI) = eccsdpt_exc(AtomI,AtomI) + energy_tmp
                     energy_res_exc = energy_res_exc + energy_tmp

                  end do
               end do

            end do
         end do
         !$OMP END PARALLEL DO

         ! get total fourth--order energy contribution
         eccsdpt_cou_occ(:,:) = 4.0E0_realk * eccsdpt_cou_occ(:,:) &
            & - 2.0E0_realk * eccsdpt_exc
         ccsdpt_e4 = 4.0E0_realk * energy_res_cou - 2.0E0_realk * energy_res_exc

         ! reorder from (a,b,j,i) to (a,b,i,j) in case of later use
         call tensor_reorder(ccsd_doubles,[1,2,4,3])

      end if
      if (.not.DECinfo%OnlyOccPart) then
         eccsdpt_exc = 0.0_realk
         energy_res_cou = 0.0E0_realk
         energy_res_exc = 0.0E0_realk
         ccsdpt_e4 = 0.0E0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_cou,eccsdpt_cou_virt), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,virt_orbitals)
         do b=1,nvirt
            do a=1,nvirt
               atomI = virt_orbitals(a)%CentralAtom
       
               do j=1,nocc
                  do i=1,nocc
       
                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_cou_virt(AtomI,AtomI) = eccsdpt_cou_virt(AtomI,AtomI) + energy_tmp
                     energy_res_cou = energy_res_cou + energy_tmp
       
                  end do
               end do
       
            end do
         end do
         !$OMP END PARALLEL DO
       
         ! reorder from (a,b,i,j) to (a,b,j,i)
         call tensor_reorder(ccsd_doubles,[1,2,4,3])
       
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,a,b,energy_tmp), &
         !$OMP REDUCTION(+:energy_res_exc,eccsdpt_exc), &
         !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,virt_orbitals)
         do b=1,nvirt
            do a=1,nvirt
               atomI = virt_orbitals(a)%CentralAtom
       
               do j=1,nocc
                  do i=1,nocc
       
                     energy_tmp = ccsd_doubles%elm4(a,b,i,j) * ccsdpt_doubles%elm4(a,b,i,j)
                     eccsdpt_exc(AtomI,AtomI) = eccsdpt_exc(AtomI,AtomI) + energy_tmp
                     energy_res_exc = energy_res_exc + energy_tmp
       
                  end do
               end do
       
            end do
         end do
         !$OMP END PARALLEL DO
       
         ! get total fourth--order energy contribution
         eccsdpt_cou_virt(:,:) = 4.0E0_realk * eccsdpt_cou_virt(:,:) &
            & - 2.0E0_realk * eccsdpt_exc
         ccsdpt_e4 = 4.0E0_realk * energy_res_cou - 2.0E0_realk * energy_res_exc
       
         ! reorder from (a,b,j,i) to (a,b,i,j) in case of later use
         call tensor_reorder(ccsd_doubles,[1,2,4,3])
      end if

      ! ******************************************************************
      ! ************** done w/ energy for full molecule ******************
      ! ******************************************************************

   end subroutine ccsdpt_decnp_e4_full


   !> \brief: calculate E[5] contribution to ccsd(t) energy correction for full molecule.
   !> \author: Janus Juul Eriksen
   !> \date: February 2013
   subroutine ccsdpt_decnp_e5_full(nocc,nvirt,nfrags,offset,ccsd_singles,ccsdpt_singles,&
         & occ_orbitals,unocc_orbitals,e5_occ,e5_virt,ccsdpt_e5)

      implicit none

      !> ccsd and ccsd(t) singles amplitudes
      type(tensor), intent(inout) :: ccsd_singles, ccsdpt_singles
      !> dimensions
      integer, intent(in) :: nocc, nvirt, nfrags, offset
      !> occupied orbital information
      type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
      !> virtual orbital information
      type(decorbital), dimension(nvirt), intent(inout) :: unocc_orbitals
      !> etot
      real(realk), intent(inout) :: ccsdpt_e5
      real(realk), dimension(nfrags,nfrags), intent(inout) :: e5_occ
      real(realk), dimension(nfrags,nfrags), intent(inout) :: e5_virt
      !> integers
      integer :: i,a,AtomI,AtomA
      !> tmp energy real
      real(realk) :: energy_tmp

      ! ***********************
      !   do E[5] energy part
      ! ***********************

      if (.not.DECinfo%OnlyVirtPart) then
         ccsdpt_e5 = 0.0_realk
         energy_tmp = 0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,a,energy_tmp,AtomI),&
         !$OMP SHARED(ccsd_singles,ccsdpt_singles,nocc,nvirt,offset,occ_orbitals),&
         !$OMP REDUCTION(+:ccsdpt_e5),REDUCTION(+:e5_occ)
         do i=1,nocc
            AtomI = occ_orbitals(i+offset)%secondaryatom
            do a=1,nvirt

               energy_tmp = ccsd_singles%elm2(a,i) * ccsdpt_singles%elm2(a,i)
               e5_occ(AtomI,AtomI) = e5_occ(AtomI,AtomI) + energy_tmp
               ccsdpt_e5 = ccsdpt_e5 + energy_tmp

            end do
         end do
         !$OMP END PARALLEL DO

         ! get total fifth-order energy correction
         e5_occ(:,:) = 2.0E0_realk * e5_occ(:,:)
         ccsdpt_e5 = 2.0E0_realk * ccsdpt_e5
      end if
      if (.not.DECinfo%OnlyOccPart) then
         ccsdpt_e5 = 0.0_realk
         energy_tmp = 0.0e0_realk
         !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,a,energy_tmp,AtomA),&
         !$OMP SHARED(ccsd_singles,ccsdpt_singles,nocc,nvirt,unocc_orbitals),&
         !$OMP REDUCTION(+:ccsdpt_e5),REDUCTION(+:e5_virt)
         do i=1,nocc
            do a=1,nvirt
               AtomA = unocc_orbitals(a)%secondaryatom

               energy_tmp = ccsd_singles%elm2(a,i) * ccsdpt_singles%elm2(a,i)
               e5_virt(AtomA,AtomA) = e5_virt(AtomA,AtomA) + energy_tmp
               ccsdpt_e5 = ccsdpt_e5 + energy_tmp

            end do
         end do
         !$OMP END PARALLEL DO

         ! get total fifth-order energy correction
         e5_virt(:,:) = 2.0E0_realk * e5_virt(:,:)
         ccsdpt_e5 = 2.0E0_realk * ccsdpt_e5
      end if


      ! ******************************
      !   done with E[5] energy part
      ! ******************************

   end subroutine ccsdpt_decnp_e5_full


   subroutine get_t1_matrices(MyLsitem,t1,Co,Cv,xo,yo,xv,yv,fock,t1fock,sync)
     implicit none
     type(lsitem),intent(inout) :: MyLsItem
     type(tensor), intent(inout) :: t1, Co, Cv
     type(tensor), intent(inout) :: xo,xv,yo,yv
     type(tensor), intent(in)    :: fock
     type(tensor), intent(inout) :: t1fock
     logical, intent(in) :: sync
     integer :: ord2(2), nb,no,nv
     real(realk), pointer :: w1(:)

     nv=t1%dims(1)
     no=t1%dims(2)
     nb=Co%dims(1)

     ! synchronize singles data on slaves
     if(sync)call tensor_sync_replicated(t1)

     ! get the T1 transformation matrices
     call tensor_cp_data(Cv,yv)
     call tensor_cp_data(Cv,xv)
     ord2 = [1,2]
     call tensor_contract(-1.0E0_realk,Co,t1,[2],[2],1,1.0E0_realk,xv,ord2)

     call tensor_cp_data(Co,yo)
     call tensor_cp_data(Co,xo)
     call tensor_contract(1.0E0_realk,Cv,t1,[2],[1],1,1.0E0_realk,yo,ord2)

     !ONLY USE T1 PART OF THE DENSITY MATRIX AND THE FOCK 
     call mem_alloc(w1,nb**2)
     call dgemm('n','n',nb,no,nv,1.0E0_realk,yv%elm1,nb,t1%elm1,nv,0.0E0_realk,t1fock%elm1,nb)
     call dgemm('n','t',nb,nb,no,1.0E0_realk,t1fock%elm1,nb,xo%elm1,nb,0.0E0_realk,w1,nb)
     call II_get_fock_mat_full(DECinfo%output,DECinfo%output,MyLsItem%setting,nb,w1,.false.,t1fock%elm1)
     call daxpy(nb**2,1.0E0_realk,fock%elm1,1,t1fock%elm1,1)
     call mem_dealloc(w1)
   end subroutine get_t1_matrices


   end module cc_tools_module
