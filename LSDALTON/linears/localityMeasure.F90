module localityMeasureMod
  use TYPEDEF,only: count_ncore
  use precision
  use orbspread_utilMod
  use kurtosis
  use matrix_util
  use loc_utils
  use typedeftype
  use matrix_module, only: matrix
  use matrix_operations 
  use matrix_util!, only: matrix_exponential
  use memory_handling
  use decompMod

  private 
  public :: LocalityMeasure

CONTAINS

  subroutine LocalityMeasure(CFG,ls,cmo,ncore,nval,nvirt)
    implicit none
    type(RedSpaceItem) :: CFG
    type(orbspread_data)  :: orbspread_input
    type(lsitem) :: ls
    type(matrix) :: cmo,CMOblock
    integer :: i,ncore,nval,nvirt,nbas,indx(1) 
    real(realk), pointer :: kurtvec(:),tmp(:)
    logical :: DoNotAllocateTmpM,Memreduced,MajorMemreduced
    DoNotAllocateTmpM = .TRUE.
    Memreduced = .TRUE.
    MajorMemreduced = .TRUE.

    nbas = cmo%nrow

    ! *** COMPUTE FOURTH MOMENT FOR ORBITALS ***
    CFG%PFM_input%crossterms=.true.
    CFG%PFM_input%m=1
    IF(.NOT.MajorMemreduced)call kurt_initAO(CFG%PFM_input,ls,cmo%nrow)
    call orbspread_propint(orbspread_input,ls)
    if (ncore > 0) then
       !######### CORE #############
       CFG%PFM_input%norb=ncore
       !get core block
       call mat_init(CMOblock,nbas,ncore)
       call mem_alloc(tmp,nbas*ncore)
       call mat_retrieve_block(CMO,tmp,nbas,ncore,1,1)
       call mat_set_from_full(tmp,1d0,CMOblock)
       call mem_dealloc(tmp)
       IF(MajorMemreduced)THEN
          call compute_memreduced_noAOinit_omega(CFG%PFM_input,CMOblock,nbas,ls)
       ELSE
          call kurt_initMO(CFG%PFM_input,CMOblock,Memreduced)
       ENDIF
       call orbspread_init(orbspread_input,1,ncore,DoNotAllocateTmpM)
       call orbspread_update(orbspread_input,CMOblock)
       call orbspread_free(orbspread_input,DoNotAllocateTmpM)
       write(CFG%lupri,*)
       write(ls%lupri,'(a)') '  %LOC%  %%%%%%%%%%%%%%% CORE LOCALITY  %%%%%%%%%%%%%%%'
       call LocalityMeasure_print(CFG,orbspread_input,ncore,0)
       write(CFG%lupri,*)
       IF(DoNotAllocateTmpM) call mem_dealloc(orbspread_input%spread2) 
       IF(MajorMemreduced)THEN
          call mem_dealloc(CFG%PFM_input%omega)
       ELSE
          call kurt_freeMO(CFG%PFM_input,Memreduced)
       ENDIF
       call  mat_free(CMOblock)
    end if
    !######### VALENCE #############
    if (nval > 0) then
       CFG%PFM_input%norb = nval
       call mat_init(CMOblock,nbas,nval)
       call mem_alloc(tmp,nbas*nval)
       call mat_retrieve_block(CMO,tmp,nbas,nval,1,ncore+1)
       call mat_set_from_full(tmp,1d0,CMOblock)
       call mem_dealloc(tmp)
       IF(MajorMemreduced)THEN
          call compute_memreduced_noAOinit_omega(CFG%PFM_input,CMOblock,nbas,ls)
       ELSE
          call kurt_initMO(CFG%PFM_input,CMOblock,Memreduced)
       ENDIF
       call orbspread_init(orbspread_input,1,nval,DoNotAllocateTmpM)
       call orbspread_update(orbspread_input,CMOblock)
       call orbspread_free(orbspread_input,DoNotAllocateTmpM)
       indx=maxloc(orbspread_input%spread2)
       CFG%leastl_occ  = indx(1)+ncore
       indx=minloc(orbspread_input%spread2)
       CFG%mostl_occ = indx(1)+ncore
       write(CFG%lupri,*)
       write(CFG%lupri,'(a)') '  %LOC%  %%%%%%%%%%%%%% VALENCE LOCALITY  %%%%%%%%%%%%%%'
       call LocalityMeasure_print(CFG,orbspread_input,nval,ncore)
       write(CFG%lupri,*)
       IF(DoNotAllocateTmpM) call mem_dealloc(orbspread_input%spread2) 
       IF(MajorMemreduced)THEN
          call mem_dealloc(CFG%PFM_input%omega)
       ELSE
          call kurt_freeMO(CFG%PFM_input,Memreduced)
       ENDIF
       call  mat_free(CMOblock)
    end if

    !######### VIRTUAL #############
    if (nvirt > 0) then
       CFG%PFM_input%norb = nvirt
       call mat_init(CMOblock,nbas,nvirt)
       call mem_alloc(tmp,nbas*nvirt)
       call mat_retrieve_block(CMO,tmp,nbas,nvirt,1,ncore+nval+1)
       call mat_set_from_full(tmp,1d0,CMOblock)
       call mem_dealloc(tmp)
       IF(MajorMemreduced)THEN
          call compute_memreduced_noAOinit_omega(CFG%PFM_input,CMOblock,nbas,ls)
       ELSE
          call kurt_initMO(CFG%PFM_input,CMOblock,Memreduced)
       ENDIF
       call orbspread_init(orbspread_input,1,nvirt,DoNotAllocateTmpM)
       call orbspread_update(orbspread_input,CMOblock)
       call orbspread_free(orbspread_input,DoNotAllocateTmpM)
       indx=maxloc(orbspread_input%spread2)
       CFG%leastl_virt = indx(1)+ncore+nval
       indx=minloc(orbspread_input%spread2)
       CFG%mostl_virt  = indx(1)+ncore+nval
       write(CFG%lupri,*)
       write(CFG%lupri,'(a)') '  %LOC%  %%%%%%%%%%%%%% VIRTUAL LOCALITY  %%%%%%%%%%%%%%'
       call LocalityMeasure_print(CFG,orbspread_input,nvirt,ncore+nval)
       IF(DoNotAllocateTmpM) call mem_dealloc(orbspread_input%spread2) 
       write(CFG%lupri,*)
       IF(MajorMemreduced)THEN
          call mem_dealloc(CFG%PFM_input%omega)
       ELSE
          call kurt_freeMO(CFG%PFM_input,Memreduced)
       ENDIF
       call  mat_free(CMOblock)
    end if


    call orbspread_propint_free(orbspread_input)
    IF(.NOT.MajorMemreduced)THEN
       call kurt_freeAO(CFG%PFM_input)
    ENDIF
  end subroutine LocalityMeasure

  subroutine LocalityMeasure_print(CFG,inp,ndim,offset)
    implicit none
    type(RedSpaceItem) :: CFG
    type(orbspread_data) :: inp
    integer, intent(in) :: ndim,offset
    real(realk),pointer :: kurtvec(:)  
    integer :: i

    call mem_alloc(kurtvec,ndim)
    if (CFG%orb_debug .or. CFG%all_orb_locality) then
       do i=1,ndim
          write(CFG%lupri, '(a,i5,f15.3,f15.3)') 'Orbital number: sigma_2,   sigma_4 :'&
               &, i+offset,dsqrt(inp%spread2(i)),&
               &dsqrt(dsqrt(CFG%PFM_input%omega(i)))
       end do
       write(CFG%lupri,*)
       write(CFG%lupri,*) '--------------------------------------------------'
       write(CFG%lupri,*)
    end if
    write(CFG%lupri,'(a,f7.2,i5)') '  %LOC%   Max. sigma_2 and orb.number: ',&
         &dsqrt(maxval(inp%spread2)),maxloc(inp%spread2)+offset
    kurtvec= dsqrt(dsqrt(CFG%PFM_input%omega))
    write(CFG%lupri,'(a,f7.2,i5)') '  %LOC%   Max. sigma_4 and orb.number: ',&
         &maxval(kurtvec),maxloc(kurtvec)+offset
    write(CFG%lupri,'(a,f7.2,i5)') '  %LOC%   Min. sigma_2 and orb.number: ',&
         &dsqrt(minval(inp%spread2)),minloc(inp%spread2)+offset
    kurtvec= dsqrt(dsqrt(CFG%PFM_input%omega))
    write(CFG%lupri,'(a,f7.2,i5)') '  %LOC%   Min. sigma_4 and orb.number: ',&
         &minval(kurtvec),minloc(kurtvec)+offset
    call mem_dealloc(kurtvec)

  end subroutine LocalityMeasure_print

end module localityMeasureMod
