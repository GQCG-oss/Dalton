!> @file
!> Operations for general arrays
!> \author Patrick Ettenhuber

module tensor_tester_module

  use tensor_parameters_and_counters
  use lspdm_tensor_operations_module
  use reorder_tester_module
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use tensor_interface_module

  contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY TESTCASES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> \brief Test the array structure 
  !> \author Patrick Ettenhuber
  !> \date September 2012
  subroutine test_tensor_struct(fu_out)
    implicit none

    integer, intent(in) :: fu_out

    type(tensor) :: test1,test2,test3
    real(tensor_dp),pointer :: dummy1(:),tileget(:),dummy2(:)
    real(tensor_dp),pointer :: tileget2(:),buf1(:), buf2(:), buf3(:)
    real(tensor_dp) :: normher,ref,ref2,ref3
    integer(kind=long) :: testint
    logical :: master
    integer :: no,nv,nb,na,i,j,succ,to_get_from,ti,midx(4),output
    integer(kind=tensor_mpi_kind) :: sender, recver, nnod, rnk, me
    character(len=7) :: teststatus
    character(tensor_MSG_LEN) :: msg
    integer,pointer :: ord(:)
    master = .true.
    nnod   = 1_tensor_mpi_kind
    me     = 0
#ifdef VAR_MPI
    me = infpar%lg_mynum
    if(me /= 0) then
      master =.false.
    endif
    nnod = infpar%lg_nodtot
    if(nnod < 3) print*,"WARNING(test_tensor_struct): not enough MPI processes to test1 all features"
    output = fu_out
#endif
    nb =  21
    nv =  18
    no =  12
    na =  7
 

#ifdef VAR_MPI
    if(master)then
      write(output,*)"TESTING PDM TT_TILED ARRAY ALLOCATIONS"
      write(output,'(" Using",f8.3," GB of mem for the testarray")')&
      &(nv*no*(nv+nb)*8.0E0_tensor_dp)/(1024.E0_tensor_dp*1024.E0_tensor_dp*1024.E0_tensor_dp)
      testint=2

      call mem_alloc(dummy1,nb*na*nv*no)
      call mem_alloc(dummy2,nb*na*nv*no)
      call random_number(dummy1)
      call print_norm(dummy1,int(nb*na*nv*no,kind=8),ref)
      write(output,'("REFERENCE NORM:",f19.12)')ref
     

      !DIFFERENT ALLOCATION AND DEALLOCATION STEPS
      write (output,*) ""
      write (output,*) ""
      write (output,*) "TESTING SIMPLE ARRAY FUNCTIONS - MASTER DIRECTED"
      write (output,*) ""
      write (output,*)"ALLOC-DEALLOC TESTS"
      print *,"alloc dealloc tests"
      teststatus="SUCCESS"
      call tensor_init(test1,[nv,na,nv,nb],4,TT_TILED_DIST,AT_MASTER_ACCESS,[nv,no-1,1,2])
      call tensor_init(test2,[na,nb,nv,no],4,TT_TILED_DIST,AT_MASTER_ACCESS,[nv,no-1,1,2])
      call tensor_free(test2)
      call tensor_init(test2,[no,no+1,no-1,no+1],4,TT_TILED_DIST,AT_MASTER_ACCESS,[no,no-1,nv,nb])
      call tensor_free(test1)
      call tensor_free(test2)
      call tensor_print_mem_info(output,.true.,.false.,succ)
      if(succ/=0)teststatus=" FAILED"
      call tensor_init(test2,[nb,no,nv,no+1],4,TT_TILED_DIST,AT_MASTER_ACCESS,[nb,2,3,4])
      write (output,'(" ALLOC-DEALLOC TESTS: ",A7)')teststatus  
      print *,"DIFFERENT ALLOCATION AND DEALLOCATION STEPS: ",teststatus

      !ALLOCATING A FULL MATRIX AND PUT IT TO DISTRIBUTED MEMORY
      !check for errors via norm
      write(output,*)""
      write(output,*)""
      teststatus="SUCCESS"
      call tensor_init(test1,[nb,na,nv,no],4,TT_TILED_DIST,AT_MASTER_ACCESS,[nb,na-1,3,no/2])
      write (output,*) "CONVERT PREVIOUS ARRAY TO PDM TT_TILED" 
      call tensor_convert(dummy1,test1,[1,2,3,4])
      call print_norm(test1,normher)
      write(output,'("NORM OF PDM ARRAY  : ",f20.15)')normher
      if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
      write (output,'("CNVRT: NORM, TEST STATUS:",f19.10," : ",A7)')normher,teststatus
      print *,"ALLOCATING A FULL MATRIX AND PUT IT TO DISTRIBUTED MEMORY: ",normher,teststatus


      !GET A TILE OF A PDM ARRAY
      !calculate how many elements are in the desired tile, and allocate the
      !respective amount of memory in a fortran array
      write(output,*)""
      write(output,*)""
      write(output,*)"TESTING MPI_GET"
      do ti = 1, test1%ntiles
        call get_residence_of_tile(to_get_from,int(ti,kind=tensor_standard_int),test1)
        if(to_get_from /= me .or. nnod==1)then
         testint = ti
         exit
        endif 
      enddo
      call get_tile_dim(j,test1,testint)
      print *,"trying to get",testint," with size", j
      call mem_alloc(tileget,j)
      call random_number(tileget)
      !initiatilize with some weird number here 10 and after get compare the
      !norms of the tile and the local fortran array
      teststatus="SUCCESS"
      if(nnod>1) call tensor_print_tile_norm(test1,ti,ref)
      write(output,'("NORM OF TILE IN ARRAY   : ",f20.15)')ref
      if(nnod>1) call print_norm(tileget,int(j,kind=8),normher)
      write(output,'("NORM OF FORT BEFORE GET : ",f20.15)')normher
      if(nnod>1) call tensor_get_tile(test1,int(ti,kind=tensor_standard_int),tileget,j)
      if(nnod>1) call print_norm(tileget,int(j,kind=8),normher)
      write(output,'("NORM OF FORT AFTER GET  : ",f20.15)')normher
      if(nnod>1)then
        if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
      else
        print *,"GET A TILE OF A PDM ARRAY has been skipped"
      endif
      write (output,'("GET: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      call mem_dealloc(tileget)
      print *,"GET A TILE OF A PDM ARRAY: ",normher,teststatus

      write(output,*)""
      write(output,*)""
      write(output,*)"TESTING MPI_PUT"
      teststatus="SUCCESS"
      do ti = test1%ntiles,1, -1
        call get_residence_of_tile(to_get_from,int(ti,kind=tensor_standard_int),test1)
        if(to_get_from /= me .or. nnod==1)then
         testint = ti
         exit
        endif 
      enddo
      call get_tile_dim(j,test1,testint)
      call mem_alloc(tileget,j)
      call random_number(tileget)
      !initiatilize with some weird number here 10 and after get compare the
      !norms of the tile and the local fortran array
      if(nnod>1) call tensor_print_tile_norm(test1,ti,normher)
      write(output,'("NORM OF TILE BEFORE PUT : ",f20.15)')normher
      if(nnod>1) call print_norm(tileget,int(j,kind=8),ref)
      write(output,'("NORM OF FORT TO PUT     : ",f20.15)')ref
      if(nnod>1) call tensor_put_tile(test1,ti,tileget,j)
      if(nnod>1) call tensor_print_tile_norm(test1,ti,normher)
      write(output,'("NORM OF TILE AFTER PUT  : ",f20.15)')normher
      if(nnod>1)then 
        if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
      else
        print *,"PUT A TILE OF A PDM ARRAY has been skipped"
      endif
      write (output,'("PUT: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      print *,"TESTING MPI_PUT",normher,teststatus
     

      !GET A TILE FROM YOURSELF AS CHECK
      !testint=1
      !call mem_dealloc(tileget)
      !j=get_tileinfo_nels(test1,testint)
      !call mem_alloc(tileget,j)
      !call tensor_get_tile(test1,1,tileget,j)
      !call print_norm(tileget,j)
      !call tensor_print_tile_norm(test1,1)
      !call sleep(4)


      !CHECK MPI PUT AND ACCUMULATE IN THE SAME WAY
      !add three to the current tile and print the norm
      write(output,*)""
      write(output,*)""
      write(output,*)"TESTING MPI_ACCUMULATE"
      teststatus="SUCCESS"
      if(nnod>1)then
        do i=1,j
          tileget(i)=tileget(i)+3.0E0_tensor_dp
        enddo
        call print_norm(tileget,int(j,kind=8),ref)
      endif
      write(output,'("NORM LOCAL ACCUMULATION : ",f20.15)')ref
      !initialize the local tile with 3 and accumulate it --> compare norm
      tileget=3.0E0_tensor_dp
      if(nnod>1) call print_norm(tileget,int(j,kind=8),normher)
      write(output,'("NORM OF FORT TO ADD:      ",f20.15)')normher
      if(nnod>1) call tensor_accumulate_tile(test1,ti,tileget,j)
      if(nnod>1) call tensor_print_tile_norm(test1,ti,normher)
      write(output,'("NORM REMOTE ACCUMULATION: ",f20.15)')normher
      !use the tile with three in it, print its norm put and compare norms
      if(nnod>1)then 
        if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
      else
        print *,"ACCUMULATE A TILE OF A PDM ARRAY has been skipped"
      endif
      write (output,'("ACC: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      print *,"TESTING MPI_ACCUMULATE: ",normher,teststatus


      call tensor_free(test1)
      call tensor_init(test1,[nb,na,nv,no],4,TT_TILED_DIST,AT_MASTER_ACCESS,[0,0,0,0])
      write(output,*)""
      write(output,*)""
      write(output,*)"TESTING CONVERSION TO FORT"
      call tensor_convert(dummy1,test1)
      teststatus="SUCCESS"
      call print_norm(dummy1,int(nb*na*nv*no,kind=8),ref)
      write(output,'("NORM OF TT_DENSE ARRAY:      ",f20.15)')ref
      call print_norm(dummy1,int(no*nv*na*nb,kind=8),normher)
      write(output,'("NORM OF PDM ARRAY :       ",f20.15)')normher
      dummy2=1.0E13_tensor_dp
      call tensor_convert(test1,dummy2)
      call print_norm(dummy2,int(no*nv*na*nb,kind=8),normher)
      write(output,'("NORM OF CONTRACTED ARRAY: ",f20.15)')normher
      if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
      write (output,'("CTR: NORM, TEST STATUS:  ",f20.15," : ",A7)')normher,teststatus
      teststatus="SUCCESS"
      do i=1,no*nv*na*nb
        if(abs(dummy1(i)-dummy2(i))>1.0E-12)then
          print *,"element",i,dummy1(i),dummy2(i)
          teststatus=" FAILED"
        endif
      enddo
      write (output,'("ORDER: TEST STATUS:                              ",A7)')teststatus
      print *,"TESTING CONVERSION TO FORT: ",teststatus



      call mem_dealloc(dummy1)
      call mem_dealloc(dummy2)
      call mem_dealloc(tileget)
      !call tensor_print_mem_info(output,.true.,.false.)
      call tensor_free(test1)
      call tensor_free(test2)

      teststatus="SUCCESS"
      call tensor_print_mem_info(output,.true.,.false.,succ)
      !call tensor_print_mem_info(output,.true.,.false.)
      if(succ/=0)teststatus=" FAILED"
       write (output,'("FIRST HALF ALLOCATION: ",A7)')teststatus  
    endif
    !get the slaves into this routine
    if(master)then
      print *,"MASTER GETTING SLAVES"
      call ls_mpibcast(ARRAYTEST,infpar%master,infpar%lg_comm)
      write (output,*)""
      write (output,*)""
      write (output,*)"TESTING PARALLEL ACCESS TO THE SAME ROUTINES"
      write (output,*)""
    else
      print *,"SLAVE ARRIVED",infpar%lg_mynum
    endif
    call ls_mpibcast(output,infpar%master,MPI_COMM_LSDALTON)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!ALL OF THE SLAVES WILL BE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    testint=2

    !initialize a matrix
    teststatus="SUCCESS"
    if(master) write (output,*)"ALLOC-DEALLOC TESTS"
    call tensor_init(test1,[nb,nb+2,nb+3,nb+4],4,TT_TILED_DIST,AT_ALL_ACCESS,[nb,nb+2,40,2])
    call tensor_init(test2,[no+3,no+2,no+1,no],4,TT_TILED_DIST,AT_ALL_ACCESS,[no,40,40,10])
    call tensor_free(test1)
    call tensor_free(test2)
    call tensor_init(test2,[nb,na,nv,no],4,TT_TILED_DIST,AT_ALL_ACCESS,[nb,nv-1,1,2])
    call tensor_free(test2)
    call tensor_print_mem_info(output,.true.,.true.,succ)
    if(succ/=0)teststatus=" FAILED"
    call tensor_init(test2,[nb,no,nv,no+1],4,TT_TILED_DIST,AT_ALL_ACCESS,[nb,2,3,4])
    if(master) write (output,'(" ALLOC-DEALLOC TESTS: ",A7)')teststatus  
    if(master) write (output,*)"DONE -- NOW COMMUNICATION"
    if(master) write(output,*)""
    if(master) write(output,*)""
    print *,"ALL-INIT ALLOC-DEALLOC TESTS",teststatus
    !call lsmpi_barrier(infpar%lg_comm)

    !IF MY RANK IS NNOD-1, PUT A MATRIX CONTAINING 10 the first tile not on the
    !current rank
    teststatus="SUCCESS"
    rnk = nnod - 1
    do ti = test1%ntiles,1, -1
      call get_residence_of_tile(to_get_from,int(ti,kind=tensor_standard_int),test1)
      if(to_get_from /= nnod-1 .and. to_get_from/=nnod-2)then
       testint = ti
       exit
      endif 
    enddo


    if((infpar%lg_mynum==rnk.or.master).and. nnod > 2)then

      recver=rnk

      if(.not.master)then
        call get_tile_dim(j,test2,ti)
        call mem_alloc(tileget,j)
        tileget = 1.0E1_tensor_dp
        call tensor_put_tile(test2,ti,tileget,j)
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        call mem_dealloc(tileget)

      else
        call ls_mpisendrecv(ref,infpar%lg_comm,recver,infpar%master)
        write(output,'("NORM PARALLEL 3LPN: ",f20.15)')ref
      endif

    else
      if(master) print*,"WARNING: skipping test1 NORM PARALLEL 3LPN, not enough nodes"
    endif


    !BEFORE rank NNOD - 2  CAN GET THE TILE
    rnk = nnod - 2

    
    call lsmpi_barrier(infpar%lg_comm)
    if((infpar%lg_mynum==rnk.or.master).and.nnod>2)then
      recver=rnk
      if(.not.master)then
        call get_tile_dim(j,test2,ti)
        call mem_alloc(tileget,j)
        call tensor_get_tile(test2,int(ti,kind=tensor_standard_int),tileget,j)
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        do i=1,j
          tileget(i) = tileget(i) + 2.4E0_tensor_dp
        enddo
        call print_norm(tileget,int(j,kind=8),normher)
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        tileget = 2.4E0_tensor_dp
        call get_midx(ti,midx,test2%ntpm,test2%mode)
        call tensor_accumulate_tile(test2,midx,tileget,j)
        call mem_dealloc(tileget)
      else
        teststatus="SUCCESS"
        call ls_mpisendrecv(normher,infpar%lg_comm,recver,infpar%master)
        write(output,'("NORM PARALLEL 2LGN: ",f20.15)')normher
        if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
        write (output,'("PUT-GET: NORM, TEST STATUS: ",f19.10," : ",A7)')normher,teststatus

        call ls_mpisendrecv(ref,infpar%lg_comm,recver,infpar%master)
        write(output,'("NORM PARALLEL 2LAC: ",f20.15)')ref
      endif
    endif

    !BE CAREFUL ABOUT WHETER THE INFORMATION IS ALREADY TRANSMITTED --> AT
    !CRITICAL POINTS INSERT BARRIER STATEMENTS TO SYNCHONIZE THE NODES 
    call lsmpi_barrier(infpar%lg_comm)

    if(nnod>2)call tensor_print_tile_norm(test2,ti,normher)
    call tensor_free(test2)
    if(master)then
       teststatus="SUCCESS"
       write(output,'("NORM PARALLEL WORK: ",f20.15)')normher
       if(nnod>2)then
         if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
       else
         print *,"AS TEST WAS SKIPPED WE DO NOT CHECK FOR THE RESULT"
       endif
       write (output,'("ACC2    : NORM, TEST STATUS: ",f19.10," : ",A7)')ref,teststatus
       print *,"ACC2    : NORM, TEST STATUS:",teststatus
    endif
    if(master) write (output,*)""

  
    !test1 extracting a tile with a different ordering than the dense matrix, put
    !that into pdm, get these tiles on each node and put them in reversed
    !reordering back into the full array, check norms and order
    teststatus="SUCCESS"
    call lsmpi_barrier(infpar%lg_comm)
    call tensor_init(test2,[no-4,nv+3,nv/7,no],4,TT_TILED_DIST,AT_ALL_ACCESS,[no-4,nv+3,5,2])
    call tensor_init(test1,[nv/7,nv+3,no,no-4],4,TT_TILED_DIST,AT_ALL_ACCESS)
    call memory_allocate_tensor_dense(test1,.false.)
    call random_number(test1%elm1)
    call lsmpi_allreduce(test1%elm1,test1%nelms,infpar%lg_comm)
    if(infpar%lg_mynum==0)then
      write (msg,*)"local test1 norm master"
      call print_norm(test1%elm1,test1%nelms,msg)
    endif
    if(infpar%lg_nodtot>1)then
      rnk = 1
    else
      rnk = 0
    endif
    if(me==rnk)then
      write (msg,*)"local test1 norm slave"
      call print_norm(test1%elm1,test1%nelms,msg)
    endif
    call print_norm(test1%elm1,test1%nelms,ref)
    call tensor_convert(test1%elm1,test2,[4,2,1,3])
    call lsmpi_barrier(infpar%lg_comm)
    call print_norm(test2,normher)
    print *,"convert",ref,normher
    call tensor_mv_dense2tiled(test1,.false.)
    call memory_allocate_tensor_dense(test2,.false.)
    call lsmpi_barrier(infpar%lg_comm)
    test2%elm1=0.0E0_tensor_dp
    do i=1,test2%ntiles
      call get_tile_dim(j,test2,i)
      call mem_alloc(tileget,j)
      call tensor_get_tile(test2,int(i,kind=tensor_standard_int),tileget,j)
      call tile_in_fort(1.0E0_tensor_dp,tileget,i,int(test2%tdim),0.0E0_tensor_dp,&
                        &test2%elm1,test1%dims,4,[3,2,4,1])
      call mem_dealloc(tileget)
    enddo
    call lsmpi_barrier(infpar%lg_comm)
    call tensor_cp_tiled2dense(test1,.true.)
    call lsmpi_barrier(infpar%lg_comm)
    if(infpar%lg_mynum==0)then
      write (msg,*)"local test1 2 norm master"
      call print_norm(test2%elm1,test2%nelms,msg)
    endif
    if(infpar%lg_nodtot>1)then
      rnk = 1
    else
      rnk = 0
    endif
    if(infpar%lg_mynum==rnk)then
      write (msg,*)"local test1 2 norm slave"
      call print_norm(test2%elm1,test2%nelms,msg)
    endif
    call print_norm(test2%elm1,test2%nelms,normher)
    do i=1,test1%nelms
      if(abs(test1%elm1(i)-test2%elm1(i))>1.0E-12)then
        teststatus=" FAILED"
      endif
    enddo
    call tensor_deallocate_dense(test1)
    call tensor_deallocate_dense(test2)
    test1%itype=TT_TILED_DIST
    test2%itype=TT_TILED_DIST
    call tensor_free(test1)
    call tensor_free(test2)
    if(master)then
       write(output,'("PDM REORDERINGS: ",f20.15)')normher
       if(abs(normher-ref)>1.0E-12_tensor_dp)teststatus=" FAILED"
       write (output,'("PDMR    : NORM, TEST STATUS: ",f19.10," : ",A7)')ref,teststatus
    endif
    if(master) write (output,*)""

    !Test parallel tensor contractions
    teststatus="SUCCESS"
    call tensor_ainit(test1,[no,nv,no,nv],4,tdims=[no/2,nv/2,no/2,nv/2],atype="TDPD")
    call tensor_ainit(test2,[nv,no,nv],   3,tdims=[nv/2,no/2,nv/2],     atype="TDPD")
    call tensor_ainit(test3,[nv,no,nv],   3,tdims=[nv/2,no/2,nv/2],     atype="TDAR")
    call tensor_zero(test3)
    !do the local reference calculation on the master
    call mem_alloc(buf3,nv*no*nv)
    if(master)then
       call random_number(test1%elm1)
       call random_number(test2%elm1)
       call print_norm(test1%elm1,int(no*nv*no*nv,kind=8))
       call print_norm(test2%elm1,int(nv*no*nv,kind=8))
       call mem_alloc(buf1,no*nv*no*nv)
       call mem_alloc(buf2,nv*no*nv)
       call array_reorder_4d(1.0E0_tensor_dp,test1%elm1,no,nv,no,nv,[1,4,2,3],0.0E0_tensor_dp,buf1)
       call array_reorder_3d(1.0E0_tensor_dp,test2%elm1,nv,no,nv,     [3,2,1],0.0E0_tensor_dp,buf2)
       call dgemm('n','n',no*nv,nv,no*nv,1.0E0_tensor_dp,buf1,no*nv,buf2,no*nv,0.0E0_tensor_dp,buf3,no*nv)
       call print_norm(buf3,int(no*nv*nv,kind=8),ref)
       call mem_dealloc(buf1)
       call mem_dealloc(buf2)
    endif
    call ls_mpibcast(test1%elm1,test1%nelms,infpar%master,infpar%lg_comm)
    call ls_mpibcast(test2%elm1,test2%nelms,infpar%master,infpar%lg_comm)
    call lsmpi_barrier(infpar%lg_comm)
    call tensor_mv_dense2tiled(test1,.true.)
    call tensor_mv_dense2tiled(test2,.true.)
    call lsmpi_barrier(infpar%lg_comm)
    call print_norm(test1)
    call print_norm(test2)

    call mem_alloc(ord,3)
    call mem_alloc(buf2,nv*no*nv)

    ord = [2,1,3]
    buf2 = 0.0E0_tensor_dp

    call tensor_lock_local_wins(test3,'e')
    call tensor_contract(1.0E0_tensor_dp,test1,test2,[2,3],[3,2],2,0.0E0_tensor_dp,test3,ord)
    call tensor_unlock_local_wins(test3)

    call mem_dealloc(ord)

    call tensor_gather(1.0E0_tensor_dp,test3,0.0E0_tensor_dp,buf2,test3%nelms,oo=[2,1,3])
    call print_norm(buf2,test3%nelms,normher)
    if(master)then
       write(output,'("PDM CONTRACTION: ",f20.15)') normher
       if(abs(normher-ref)>1.0E-10_tensor_dp) teststatus=" FAILED"
       write (output,'("PDCONTR : NORM, TEST STATUS: ",f19.10," : ",A7)') ref,teststatus
    endif
    if(master)then
       do i=1,test3%nelms
          if(abs(buf2(i)-buf3(i))>1.0E-10_tensor_dp)then
             teststatus=" FAILED"
             if(infpar%lg_mynum==0)print *,i,buf2(i),buf3(i)
          endif
       end do
       write (output,'("PDCWORD :       TEST STATUS:                     : ",A7)') teststatus
    endif 

    call mem_dealloc(buf2)
    call mem_dealloc(buf3)

    call tensor_free(test1)
    call tensor_free(test2)
    call tensor_free(test3)

#endif

  end subroutine test_tensor_struct 
end module tensor_tester_module


#ifdef VAR_MPI
subroutine get_slaves_to_tensor_test()
  use tensor_tester_module,only:test_tensor_struct
  implicit none
  call test_tensor_struct(6)
end subroutine get_slaves_to_tensor_test
#endif
