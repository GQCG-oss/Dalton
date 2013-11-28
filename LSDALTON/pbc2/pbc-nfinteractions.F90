#ifdef MOD_UNRELEASED
MODULE pbc_interactions
  use precision
  use fundamental
  USE TYPEDEF
  use memory_handling
  use integralinterfaceMOD
  USE matrix_module
  use matrix_operations
  !USE multipole_pbc
  use pbc_MSC
  USE pbc_matrix_operations
!  USE harmonics_pbc
  USE lattice_type
  USE lattice_vectors
  use lstiming
  CONTAINS


!Finds the cutoff distance to the reference cell for one-op
SUBROUTINE find_cutoff_onep(lupri,luerr,setting,nbast,lattice, &
                       latt_cell,refcell,numvecs)
  implicit none
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri, luerr,nbast,numvecs ! nlayer 
  TYPE(moleculeinfo),INTENT(INout),target :: latt_cell(numvecs)
  TYPE(moleculeinfo),intent(inout),target :: refcell
  TYPE(lvec_list_t),intent(INOUT) :: lattice
  TYPE(LSSETTING),intent(inout):: SETTING 

  ! local variables
  INTEGER(short),allocatable        :: GAB1(:,:)!(naobatch,naobatch)
  TYPE(MATRIX)  :: S!,Gab1
  REAL(realk) :: overlap,origin(3)!,distance
  INTEGER :: j,k,nbatch1,nbatch2
  INTEGER ::  idx
  INTEGER :: l1,l2,l3
  INTEGER, DIMENSION(3) :: fdim
  LOGICAL :: ntimes
!  DOUBLE PRECISION latcoord(3), origin(3),latt_vec_std(3)
  REAL(realk) :: latt_vec_std(3)!,distance

 ! write(lupri,*) 'nbast', nbast
  call mat_init(S,nbast,nbast)
  call mat_zero(S)
  !call mat_init(Gab1,nbast,nbast)
  !call mat_zero(Gab1)
  origin(1:3)=1.0
  latt_vec_std=0.0
  ntimes = .true.
!  write(*,*) 'debug 1'
!  write(*,*) setting%Output%resultTensor%nbast
!  write(*,*) 'debug 2'
  call set_lstime_print(.false.)
  fdim(:)= 0
  if (lattice%ldef%is_active(1)) fdim(1) = 1
  if (lattice%ldef%is_active(2)) fdim(2) = 1
  if (lattice%ldef%is_active(3)) fdim(3) = 1

  DO idx=1,size(lattice%lvec)
     
     !call find_latt_index(idx,l1,l2,l3,fdim,lattice,lattice%max_layer)
     call find_latt_vectors(idx,l1,l2,l3,fdim,lattice)

     !call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2)
     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
     !setting%samemol(1,2)=.false.
     !setting%samemol(2,1)=.false.
     call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast, lattice%lvec(idx)%maxGab)

     !call mat_print(Gab1,1,nbast,1,nbast,lupri)
     !setting%samemol(1,2)=.true.
     !setting%samemol(2,1)=.true.
     !write(lupri,*) 'screening matrix'
     !call mat_print(Gab1,1,nbast,1,nbast,lupri)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3)
     !setting%samemol(1,3)=.false.
     !setting%samemol(3,1)=.false.

    ! call II_get_overlap(lupri,luerr,setting,S)
!     CALL mat_print(S,1,S%nrow,1,S%ncol,6) 
     !setting%samemol(1,3)=.true.
     !setting%samemol(3,1)=.true.
    ! call find_latt_index(idx,-l1,l2,l3,fdim,lattice,lattice%max_layer)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast,lattice%lvec(idx)%maxGab)
    ! call find_latt_index(idx,l1,-l2,l3,fdim,lattice,lattice%max_layer)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast,lattice%lvec(idx)%maxGab)
    ! call find_latt_index(idx,l1,l2,-l3,fdim,lattice,lattice%max_layer)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast,lattice%lvec(idx)%maxGab)
    ! call find_latt_index(idx,l1,-l2,-l3,fdim,lattice,lattice%max_layer)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast,lattice%lvec(idx)%maxGab)
    ! call find_latt_index(idx,-l1,l2,-l3,fdim,lattice,lattice%max_layer)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast,lattice%lvec(idx)%maxGab)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call find_latt_index(idx,-l1,-l2,l3,fdim,lattice,lattice%max_layer)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast,lattice%lvec(idx)%maxGab)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call find_latt_index(idx,-l1,-l2,-l3,fdim,lattice,lattice%max_layer)
    ! call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),2,refcell,3,latt_cell(idx),4)
    ! call II_get_maxGabelm_ScreenMat(lupri,luerr,setting,nbast,lattice%lvec(idx)%maxGab)

    !  overlap=0d0
    !  DO j=1,nbast*nbast
    !    overlap=overlap+abs(S%elms(j))
    !  ENDDO
 
      !if(overlap .lt. 1E-10 .and. ntimes) then
      !  k=max(l3,l2)
      !  lattice%nneighbour = max(k,l1)
      !  ntimes=.false.
      !  write(lupri,*) 'lattice%nneighbour =', lattice%nneighbour
      !  write(*,*) 'lattice%nneighbour =', lattice%nneighbour
      !  !call mat_print(s,1,s%nrow,1,s%ncol,6)
      !endif 
      !if(.not. ntimes) exit

  !  END DO
  ! END DO
  ENDDO


!  deallocate(latt_cell)
  call set_lstime_print(.true.)
  call mat_free(S)
  !call lsquit('Quiting cutoffonepart,reason: debug',lupri)
  !call mat_free(Gab1)
!  CALL mat_print(S,1,S%nrow,1,S%ncol,6) 
!  write(lupri,*) 'tester S: ', S%elms
!  STOP
   

END SUBROUTINE find_cutoff_onep


!Finds the cutoff distance to the reference cell for two-op
SUBROUTINE find_cutoff_twop(lupri,luerr,setting,nbast,lattice, &
                       latt_cell,refcell,numvecs,nfdensity)
  implicit none
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri, luerr,nbast,numvecs ! nlayer 
  TYPE(moleculeinfo),INTENT(IN) :: latt_cell(numvecs)
  TYPE(moleculeinfo),intent(inOUT) :: refcell
  TYPE(matrix),intent(in) :: nfdensity(numvecs)
  TYPE(lvec_list_t),intent(INOUT) ::lattice

  ! local variables
  TYPE(LSSETTING),intent(inout):: SETTING 
  TYPE(MATRIX),pointer  ::  Kx(:), K_tmp(:)
  REAL(realk) :: exact_xch,origin(3),distance
  INTEGER :: j,k,checknf1,checknf2,checknf3,checknf
  INTEGER ::  idx,index2,index3,newcell
  INTEGER :: l1,l2,l3,il1,il2,il3
  INTEGER :: maxl1,maxl2,maxl3
  INTEGER :: l11,l12,l13,il21,il22,il23,il31,il32,il33
  INTEGER, DIMENSION(3) :: fdim
  LOGICAL :: ntimes
!  DOUBLE PRECISION latcoord(3), origin(3),latt_vec_std(3)
  REAL(realk) :: latt_vec_std(3)!,distance

  call set_lstime_print(.false.)
  write(lupri,*) 'Finding nlayer for Kop',numvecs, lattice%ldef%is_active(1),&
  lattice%max_layer

  call mem_alloc(K_tmp,1)
  call mem_alloc(Kx,1)

  call mat_init(Kx(1),nbast,nbast)
  call mat_zero(Kx(1))
  call mat_init(K_tmp(1),nbast,nbast)
  call mat_zero(K_tmp(1))


  origin(1:3)=1.0
  latt_vec_std=0.0
  ntimes = .false.

  fdim(:)= 0
  if (lattice%ldef%is_active(1)) fdim(1) = 1
  if (lattice%ldef%is_active(2)) fdim(2) = 1
  if (lattice%ldef%is_active(3)) fdim(3) = 1

  DO l3=0,lattice%max_layer*fdim(3)
   if(ntimes) exit
   DO l2=0,lattice%max_layer*fdim(2)
    if(ntimes) exit
    DO l1=0,lattice%max_layer*fdim(1)
    if(ntimes) exit
     
    call find_latt_index(idx,l1,l2,l3,fdim,lattice,lattice%max_layer)

    call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3)

     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.




  DO index2=1,numvecs

     call find_latt_vectors(index2,il21,il22,il23,fdim,lattice)
     !if(abs(il21) .gt.15) CYCLE!lattice%nneighbour) 
     !if(abs(il22) .gt.15) CYCLE!lattice%nneighbour) 
     !if(abs(il23) .gt.15) CYCLE!lattice%nneighbour) 

     !Skal kun ha avstand til origo ikke latt_std_vec her
     !call calc_distance(distance,lattice%lvec(index2)%std_coord,latt_vec_std)


   DO index3=1,numvecs

    !write(lupri,*) 'electron rep get'
    !CALL mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,6)
    !CALL mat_print(D1(1),1,F(1)%nrow,1,F(1)%ncol,6)

    call find_latt_vectors(index3,il31,il32,il33,fdim,lattice)


!   ToDo Remove: base on CS screening at some point
    l11=il21+il31
    if( abs(l11) .gt.lattice%max_layer) CYCLE
    checknf1=l11-il21
    !l1 blir større enn max slik at newcell blir negativ
    !sjekk dette, er ikke sikker på om l1 skal være slik.
    !Likedan man l2 og l3
    if(abs(il31) .gt. lattice%ndmat) CYCLE! then
    
    l12=il22+il32
    if( abs(l12) .gt.lattice%max_layer) CYCLE
    checknf2=l12-il22
    if(abs(il32) .gt. lattice%ndmat) CYCLE !then 
    
    l13=il23+il33
    if( abs(l13) .gt.lattice%max_layer) CYCLE
    checknf3=l13-il23
    if(abs(il33) .gt. lattice%ndmat) CYCLE !then 

    

    call find_latt_index(newcell,l11,l12,l13,fdim,lattice,lattice%max_layer)
    call find_latt_index(checknf,il31,il32,il33,&
                          fdim,lattice,lattice%max_layer)

   ! write(lupri,*) index1,index2,newcell,l1,l2,l3

    setting%samefrag=.false.
    setting%samemol=.false.
    ! ToDo Should be turned on at some point

    call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3,latt_cell(index2),2,latt_cell(newcell),4)

    call mat_zero(K_tmp(1))
    call II_get_exchange_mat(lupri,luerr,setting,nfdensity(checknf:checknf),&
         &                   1,.false.,K_tmp)

    call mat_daxpy(0.5D0,K_tmp(1),Kx(1))

   ENDDO
  ENDDO

   exact_xch=0.0_realk
   DO j=1,nbast*nbast
     exact_xch=exact_xch+abs(Kx(1)%elms(j))
   ENDDO
 
   if(exact_xch .lt. 1E-10 .and. .not. ntimes) then
     k = max(l3,l2)
     k = max(k,l1)
     lattice%kx1=k
     lattice%kx2=k
     lattice%kx3=k
     ntimes=.true.
     write(lupri,*) 'nlayer Kx =', k
     write(lupri,*) 'lattice%kx1 =', lattice%kx1
     write(lupri,*) 'lattice%kx2 =', lattice%kx3
     write(lupri,*) 'lattice%kx2 =', lattice%kx3
   endif 

  CALL mat_zero(Kx(1))

  if(ntimes) exit

    END DO
   END DO
  ENDDO

!  deallocate(latt_cell)
   call mat_free(K_tmp(1))
   call mat_free(Kx(1))
   call mem_dealloc(K_tmp)
   call mem_dealloc(Kx)
  call set_lstime_print(.true.)
   write(lupri,*) 'Finished with find_cutoff_two_p'
!  CALL mat_print(S,1,S%nrow,1,S%ncol,6) 
!  write(lupri,*) 'tester S: ', S%elms
!  STOP

END SUBROUTINE find_cutoff_twop


SUBROUTINE screening_ovl(basinfo)
IMPLICIT NONE
TYPE(basissetinfo) :: basinfo
INTEGER :: natomtypes
INTEGER :: i,j,k
!REAL(realk) :: PI=3.14159265
REAL(realk), pointer :: minexp(:)
REAL(realk) :: distance

  natomtypes=basinfo%natomtypes
  
  k=0
  DO i=1,natomtypes
  DO j=1,basinfo%atomtype(i)%nangmom
  k=k+1
  ENDDO
  ENDDO
  call mem_alloc(minexp,k)
  k=0
  DO i=1,natomtypes
  DO j=1,basinfo%atomtype(i)%nangmom
  k=k+1
  minexp(k)=minval(basinfo%atomtype(i)%shell(j)%segment(2)%exponents)
  ENDDO
  ENDDO
  
  !write(lupri,*) 'expon ', minval(minexp)
  distance=sqrt(1./minval(minexp)*log((pi/(2.*minval(minexp)))**3E0_realk*10E0_realk**20))
  !write(lupri,*) 'distance ',distance

!  distance=sqrt(1./0.2979640)*sqrt(log((pi/(2.*0.2979640))**(3)*10**20))
!  write(lupri,*) 'distance ',distance
   call mem_dealloc(minexp)


END SUBROUTINE screening_ovl

!Finds distance between two vectors
SUBROUTINE calc_distance(distance,tvec,latstdvec)
  IMPLICIT NONE
  REAL(realk), INTENT(IN) :: tvec(3),latstdvec(3)
  REAL(realk), INTENT(OUT) :: distance

  !a=tvec
  distance = (tvec(1)-latstdvec(1))**2 + (tvec(2)-latstdvec(2))**2 &
             +(tvec(3)-latstdvec(3)**2)

  distance=sqrt(distance)
  !write(lupri,*) 'debug: ', distance

END SUBROUTINE calc_distance



!SUBROUTINE for computing overlap integrals over different cells
!including for the reference cell.
SUBROUTINE pbc_overlap_k(lupri,luerr,setting,molecule,nbast,lattice,&
 latt_cell,refcell,numvecs,ovl)

  implicit none
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri, luerr, nbast,numvecs
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(moleculeinfo),intent(inout)          :: refcell
  !TYPE(BASISSETINFO) :: basinfo
  !COMPLEX(COMPLEXK), INTENT(INOUT) :: overlap(sizeov*sizeov)
  TYPE(lvec_list_t),intent(INOUT) ::lattice
  TYPE(matrix),target :: ovl(numvecs)

  ! local variables
  REAL(realk) :: kvec(3)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(moleculeinfo),intent(in) :: latt_cell(numvecs)
  TYPE(MATRIX)  :: S,S1
  INTEGER :: i,j,m,il1,il2,il3
  INTEGER :: idx,refindex
  INTEGER :: natoms,iunit,maxl1,maxl2,maxl3
  INTEGER, DIMENSION(3) :: fdim
  integer(short) :: gab1
  REAL(realk) :: latt_vec_std(3),phase1,phase2,phase3,origin(3)
  COMPLEX(COMPLEXK) :: phase
!  REAL(realk) :: PI=3.14159265358979323846D0
  INTEGER,save :: iter=0
  CHARACTER(len=10) :: numstr1,numstr2,numstr3
  CHARACTER(len=40) :: filename
  iter=iter+1



  natoms=molecule%natoms!number of atoms

  call set_lstime_print(.false.)

  call mat_init(S,nbast,nbast)
  call mat_zero(S)

  call mat_init(S1,nbast,nbast)
  call mat_zero(S1)

  origin(1:3)=1
  latt_vec_std=0.0

  write(lupri,*) 'Number of lattice vectors ', numvecs

  i=1
  !finds the lattice index for the reference cell.
  !i.e for l1=l2=l3=0
  call find_latt_index(refindex,0,0,0,fdim,lattice,lattice%max_layer)
  !loop over lattice cells
  maxl1=0
  maxl2=0
  maxl3=0
  DO idx=1,numvecs
     
     lattice%lvec(idx)%ovl_computed=.false.

  !Doing translations
     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3)
     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.

     call find_latt_vectors(idx,il1,il2,il3,fdim,lattice)
     !So that we do not consider negligible integrals
     !if(abs(il1) .gt. lattice%nneighbour) CYCLE
     !if(abs(il2) .gt. lattice%nneighbour) CYCLE
     !if(abs(il3) .gt. lattice%nneighbour) CYCLE
     gab1=lattice%lvec(idx)%maxgab
     if(gab1 .ge. -12) then
       maxl1=max(abs(il1),maxl1)
       maxl2=max(abs(il2),maxl2)
       maxl3=max(abs(il3),maxl3)
       lattice%lvec(idx)%ovl_computed=.true.
       !input%Basis%regular%atomtype(:)%shell(:)%segment(1)%exponents(:)

       if(.not. lattice%store_mats) then
         call mat_init(ovl(idx),nbast,nbast)
       endif
       call mat_init(lattice%lvec(idx)%oper(1),nbast,nbast)
       call mat_zero(lattice%lvec(idx)%oper(1))


   !g  et the overlap matrix for cell between reference 
   !a  nd cell l
   !    call II_get_overlap(lupri,luerr,setting,S1)
       call II_get_overlap(lupri,luerr,setting,lattice%lvec(idx)%oper(1))

       !S%elms(:)=S1%elms(:)*exp(phase)
       !call mat_add(1.D0,S,1.D0,S1,S)!phase kan gå inn her istedet
    !  call mat_daxpy(1.D0,lattice%lvec(idx)%oper(1),lattice%lvec(idx)%oper(3))!phase kan gå inn her istedet
  !     CALL mat_print(S1,1,S1%nrow,1,S1%ncol,lupri)
       !call mat_to_full(S1,1d0,lattice%lvec(idx)%Sl_mat)
       !DO k=1,nbast*nbast
       !   lattice%lvec(idx)%Sl_mat=S1%elms(k)
       !ENDDO

#ifdef DEBUGPBC
       write(lupri,*) 'OVERLAP MAT',il1
       call mat_print(lattice%lvec(idx)%oper(1),1,nbast,1,nbast,lupri)
#endif

       if(.not. lattice%store_mats) then
         call mat_copy(1D0,lattice%lvec(idx)%oper(1),ovl(idx))
       endif

       if(lattice%store_mats) then
         call pbc_get_file_and_write(lattice,nbast,nbast,idx,1,1,'            ')!1 refers to overlap
       endif

       
       !write(lupri,*) 'S%elms', S%elms!,nbast
       !For debugging only
       if(lattice%compare_elmnts .and. iter .eq. 1) then
         !compare integrals with the old pbc code

       !  write(lupri,*) 'comparing overlap elements with old pbc code'
       !  write(lupri,*) il1,il2,il3
         
         iunit=-1
         write(numstr1,'(I5)')  il1
         write(numstr2,'(I5)')  il2
         write(numstr3,'(I5)')  il3
         numstr1=adjustl(numstr1)
         numstr2=adjustl(numstr2)
         numstr3=adjustl(numstr3)
         filename='minovl'//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'

         CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
         DO j=1,nbast
            write(iunit,*) (lattice%lvec(idx)%oper(1)%elms(i+(j-1)*nbast),i=1,nbast)
         ENDDO
         call lsclose(iunit,'KEEP')

       endif

       call mat_free(lattice%lvec(idx)%oper(1))
     endif
       !if(phase .ne. cmplx(0.d0,0d0)) then
     !  write(lupri,*) 'phase= ',phase
     !  call write_matrix(kvec,3,3)
     !  stop
     !endif
       !CALL mat_print(S1,1,S1%nrow,1,S1%ncol,6)

   !should write to disk in stead.
  END DO
  lattice%oneop1=maxl1
  lattice%oneop2=maxl2
  lattice%oneop3=maxl3

  call mat_free(S1)
  call mat_free(S)
  call set_lstime_print(.true.)
  write(lupri,*) 'finished pbc_overlap_k'
!  STOP

END SUBROUTINE pbc_overlap_k

!
SUBROUTINE pbc_kinetic_k(lupri,luerr,setting,molecule,nbast,lattice,&
latt_cell,refcell,numvecs,nfdensity,f_1,E_kin)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr,nbast,numvecs
  !REAL(realk), INTENT(IN) :: kvec(3,3)
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(moleculeinfo),intent(inout)  :: refcell
  TYPE(moleculeinfo), DIMENSION(numvecs),INTENT(IN) :: latt_cell
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  REAL(realk) :: E_kin
  TYPE(lvec_list_t),intent(INOUT) ::lattice
  TYPE(matrix),INTENT(IN) :: nfdensity(numvecs)
  TYPE(matrix),INTENT(INOUT),target :: f_1(numvecs)

  ! local variables
  REAL(realk) :: kvec(3)
  REAL(realk), DIMENSION(12) :: Tvec
  TYPE(MATRIX)  :: kin
  REAL(realk) :: std_vec_length,E_h1,kinetic2,kineticst
  INTEGER :: i,j,k,m,il1,il2,il3,s,t,st
  INTEGER :: idx, refindex,indred
  INTEGER :: num_latvectors, natoms,iunit
  INTEGER(short) :: gab1
  INTEGER, DIMENSION(3) :: fdim
  REAL(realk) :: latt_vec_std(3), origin(3)
  REAL(realk) :: phase1,phase2,phase3
!  REAL(realk) :: PI=3.14159265358979323846D0
  COMPLEX(complexk) :: phase
  INTEGER,SAVE :: iter=0
  CHARACTER(len=10) :: numstr1,numstr2,numstr3
  CHARACTER(len=40) :: filename

  write(lupri,*) 'Starting pbc_kinetic_k'
  call set_lstime_print(.false.)
  natoms=molecule%natoms
  

  call mat_init(kin,nbast,nbast)
  call mat_zero(kin)

  iter=iter+1




  origin(1:3)=1D0
  latt_vec_std=0.0

  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)

  fdim(1:3) = 0

  call find_latt_index(refindex,0,0,0,fdim,lattice,lattice%max_layer)
  write(lupri,*) 'reference cell index=', refindex,lattice%max_layer
  tvec(1:12)=0.0
  i=1
  E_h1=0
  E_kin=0._realk

  DO idx=1,num_latvectors
     


    lattice%lvec(idx)%f1_computed=.false.
    if(.not. lattice%lvec(idx)%is_redundant) then
      call find_latt_vectors(idx,il1,il2,il3,fdim,lattice)
      call find_latt_index(indred,-il1,-il2,-il3,fdim,lattice,lattice%max_layer)
      if(il1**2+il2**2+il3**2 .gt. 0) lattice%lvec(indred)%is_redundant =.true.
      call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(idx),3)
      setting%samemol(1,3)=.false.
      setting%samemol(3,1)=.false.

      !call calc_distance(distance,lattice%lvec(idx)%std_coord,latt_vec_std)
      !if(abs(il1) .gt. lattice%nneighbour) CYCLE
      !if(abs(il2) .gt. lattice%nneighbour) CYCLE
      !if(abs(il3) .gt. lattice%nneighbour) CYCLE
    
      gab1=lattice%lvec(idx)%maxgab
      if(gab1 .ge. -12) then
        lattice%lvec(idx)%f1_computed=.true.
        lattice%lvec(indred)%f1_computed=.true.
        call mat_init(lattice%lvec(idx)%oper(1),nbast,nbast)
        call mat_zero(lattice%lvec(idx)%oper(1))
        call mat_init(f_1(idx),nbast,nbast)
        call mat_zero(f_1(idx))
        !call mat_init(lattice%lvec(idx)%oper(3),nbast,nbast)
        !call mat_zero(lattice%lvec(idx)%oper(3))


        !call II_get_kinetic(lupri,luerr,setting,kin)
        
        call II_get_kinetic(lupri,luerr,setting,lattice%lvec(idx)%oper(1))

        !call mat_daxpy(1.D0,lattice%lvec(idx)%oper(1),lattice%lvec(idx)%oper(3))
        call mat_copy(1._realk,lattice%lvec(idx)%oper(1),f_1(idx))

        call pbc_get_file_and_write(lattice,nbast,nbast,idx,2,1,'            ')!2 refers to kin

#ifdef DEBUGPBC
        write(lupri,*) 'KINETIC MAT',il1
        call mat_print(lattice%lvec(idx)%oper(1),1,nbast,1,nbast,lupri)
#endif

        if(nfdensity(idx)%init_magic_tag.EQ.mat_init_magic_value) THEN
        if(il1**2+il2**2+il3**2 .gt. 0) then
          E_kin=E_kin+2.*mat_dotproduct(lattice%lvec(idx)%oper(1),nfdensity(idx))
        else
          E_kin=E_kin+mat_dotproduct(lattice%lvec(idx)%oper(1),nfdensity(idx))
        endif
        endif

!         DO k=1,nbast*nbast
!            lattice%lvec(idx)%fck_vec(k)=lattice%lvec(idx)%oper(1)%elms(k)
!         ENDDO

#ifdef DEBUGPBC
!       write(*,*) 'kinetic DEBUGPBC'
        st=0
        kinetic2=0._realk
        kineticst=0._realk
        DO s=1,nbast
         DO t=1,nbast
          st=st+1
          kinetic2=kinetic2+lattice%lvec(idx)%oper(1)%elms(st)**2
          kineticst=kineticst+st*lattice%lvec(idx)%oper(1)%elms(st)**2
         ENDDO
        ENDDO
        write(lupri,'(A29)') 'DEBUGPBC kinetic contribution'
        write(lupri,'(A10,X,I3)') 'iteration:',iter
        write(lupri,'(A9,X,I3,x,I3,X,I3)') 't1 t2 t3:',il1,il2,il3
        write(lupri,'(A4,X,E16.8)') 'T^2:',kinetic2
        write(lupri,'(A10,X,E16.8)') 'sum ijT^2:',kineticst
#else

        if(lattice%compare_elmnts) then
          !compare integrals with the old pbc code

          !write(lupri,*) 'comparing overlap elements with old pbc code'
          !write(*,*) 'comparing overlap elements with old pbc code'
          
          iunit=-1
          write(numstr1,'(I5)')  il1
          write(numstr2,'(I5)')  il2
          write(numstr3,'(I5)')  il3
          numstr1=adjustl(numstr1)
          numstr2=adjustl(numstr2)
          numstr3=adjustl(numstr3)
          filename='minkin'//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'
!!          call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
          CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
!          write(iunit,*) nbast
          DO j=1,nbast
             write(iunit,*) (lattice%lvec(idx)%oper(1)%elms(i+(j-1)*nbast),i=1,nbast)
          ENDDO
          call lsclose(iunit,'KEEP')

        endif
#endif
        write(lupri,*) 'Kinetic mat finished for', il1,il2,il3

        call mat_free(lattice%lvec(idx)%oper(1))
       
      endif ! maxgab

      i=1
      j=0
    endif !is_redundant
      

       !write(lupri,*) 'kpnt, phase and lvec',bz%kpnt(m)%lambda(1),phase,lattice%lvec(idx)%lat_coord(1)
  END DO
  
    ! do m=1,bz%nk
    !  write(lupri,*) 'kpoint ',bz%kpnt(m)%lambda(1)
    !  !write(lupri,*) kdep(m)%kfockvec(k)
    !  call write_zmatrix(kdep(m)%kfockvec,nbast,nbast)
    ! enddo

#ifdef DEBUGPBC
 write(lupri,*) 'Kinetic energy', E_kin
 write(*,*) 'Kinetic energy', E_kin
#endif

  call mat_free(kin)
  call set_lstime_print(.true.)
  write(lupri,*) 'Finished pbc_kinetic_k'

!!
END SUBROUTINE pbc_kinetic_k
!
!
!SUBROUTINE pbc_nucattrc_k(lupri,luerr,setting,molecule,nbast,fock_mtx,lattice,latt_cell,refcell,numvecs,kvec)
SUBROUTINE pbc_nucattrc_k(lupri,luerr,setting,molecule,nbast,lattice,&
latt_cell,refcell,numvecs,nfdensity,f_1,E_en)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr,nbast,numvecs
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule,refcell
  !REAL(realk), INTENT(IN) :: kvec(3,3)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(moleculeinfo), DIMENSION(numvecs),INTENT(IN) :: latt_cell
  TYPE(matrix),intent(in) :: nfdensity(numvecs)
  TYPE(matrix),intent(inout),target :: f_1(numvecs)
!  TYPE(DALTONINPUT) :: INPUT
  TYPE(MATRIX)  :: H,H1
  REAL(realk),INTENT(INOUT) :: E_en

  ! local variables
  REAL(realk) :: kvec(3),nucatr2,nucatrst
  REAL(realk), DIMENSION(12) :: Tvec
  REAL(realk) :: std_vec_length,E_h1
  INTEGER :: i,j,k,m,il11,il12,il13,il21,il22,il23
  INTEGER :: index1,index2, refindex,s,t,st
  INTEGER :: num_latvectors, natoms,iunit,indred
  INTEGER(short) :: gab1
  INTEGER, DIMENSION(3) :: fdim
  REAL(realk):: latt_vec_std(3)
  REAL(realk) :: phase1,phase2,phase3,origin(3)
  COMPLEX(COMPLEXK) :: phase
  !COMPLEX(COMPLEXK) :: nucel_k(bz%nk,nbast*nbast)
  TYPE(lvec_list_t),intent(INOUT) ::lattice
  INTEGER,SAVE :: iter=0
  CHARACTER(len=10) :: numstr1,numstr2,numstr3
  CHARACTER(len=40) :: filename

  natoms=molecule%natoms
  call set_lstime_print(.false.)

  iter=iter+1

  !call mem_alloc
  call mat_init(H,nbast,nbast)
  call mat_zero(H)
  call mat_init(H1,nbast,nbast)
  call mat_zero(H1)

  origin(1:3)=1D0
  latt_vec_std=0.0


  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)
  
  fdim(1:3) = 0
  if (lattice%ldef%is_active(1)) fdim(1) = 1
  if (lattice%ldef%is_active(2)) fdim(2) = 1
  if (lattice%ldef%is_active(3)) fdim(3) = 1


  call find_latt_index(refindex,0,0,0,fdim,lattice,lattice%max_layer)
  write(lupri,*) 'reference cell index=', refindex

  tvec(1:12)=0.0
 
  E_en=0._realk

  DO index1=1,num_latvectors



!  call calc_distance(distance,lattice%lvec(index1)%std_coord,latt_vec_std)
  call find_latt_vectors(index1,il11,il12,il13,fdim,lattice)
  !if(abs(il11) .gt. lattice%nneighbour) CYCLE
  !if(abs(il12) .gt. lattice%nneighbour) CYCLE
  !if(abs(il13) .gt. lattice%nneighbour) CYCLE

  lattice%lvec(index1)%Vz_computed=.false.
  if(.not. lattice%lvec(index1)%is_redundant) then
    call find_latt_index(indred,-il11,-il12,-il13,fdim,lattice,lattice%max_layer)
    if(il11**2+il12**2+il13**2 .gt. 0) lattice%lvec(indred)%is_redundant =.true.
    gab1=lattice%lvec(index1)%maxgab
    if(gab1 .ge. -12) then
      lattice%lvec(index1)%Vz_computed=.true.
      lattice%lvec(indred)%Vz_computed=.true.
      call mat_init(lattice%lvec(index1)%oper(1),nbast,nbast)
      call mat_zero(lattice%lvec(index1)%oper(1))

      !phase=k1*il1+k2*il2+k3*il3
      !phase1 = dot_product(kvec(:,1),lattice%ldef%avec(:,1)*il11)
      !phase2 = dot_product(kvec(:,2),lattice%ldef%avec(:,2)*il12)
      !phase3 = dot_product(kvec(:,3),lattice%ldef%avec(:,3)*il13)
      !phase=CMPLX(0.,phase1+phase2+phase3)

      !nucel_k=0.0d0
      DO index2=1,num_latvectors


         call find_latt_vectors(index2,il21,il22,il23,fdim,lattice)
         if(abs(il21) .gt. lattice%nf) CYCLE 
         if(abs(il22) .gt. lattice%nf) CYCLE
         if(abs(il23) .gt. lattice%nf) CYCLE

        ! call calc_distance(distance,lattice%lvec(index2)%std_coord,latt_vec_std)

         call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1),2,latt_cell(index2),3)

         !setting%samemol(1,3)=.false.
         !setting%samemol(3,1)=.false.
         !setting%samemol(1,2)=.false.
         !setting%samemol(2,1)=.false.
         !setting%samemol(3,2)=.false.
         !setting%samemol(2,3)=.false.
         !These two screenings should be turned on at some point

         call II_get_nucel_mat(lupri,luerr,setting,H1)
         !ToDo: change mat_add to mat_daxpy
         !call mat_daxpy(1.D0,H1,H)
         call mat_daxpy(1.D0,H1,lattice%lvec(index1)%oper(1))
         !call mat_add(1.D0,H,1.D0,H1,H)
             


      ENDDO

!         i=1
!         j=0
!         DO k=1,nbast*nbast
!            j=j+1
!            if(j .gt. nbast) Then
!              j=1
!              i=i+1
!            ENDIF
!            fock_tmp(j,i)=H%elms(k)
!         ENDDO
         !fock_tmp=H%elms!*exp(phase)
#ifdef DEBUGPBC
      write(lupri,*) 'nucattrc mat',il11
      call mat_print(lattice%lvec(index1)%oper(1),1,H%nrow,1,H%ncol,lupri)
#endif
      !call write_matrix(fock_tmp,2,2)
      !STOP
      i=0

!     DO k=1,nbast*nbast
!        lattice%lvec(index1)%fck_vec(k)=lattice%lvec(index1)%fck_vec(k)&
!          +lattice%lvec(index1)%oper(1)%elms(k)
!     ENDDO

#ifdef DEBUGPBC
      !write(lupri,*) 'Hamilton',il11
      !call mat_print(f_1(index1),1,H%nrow,1,H%ncol,lupri)
      st=0
      nucatr2=0._realk
      nucatrst=0._realk
      DO s=1,nbast
       DO t=1,nbast
          st=st+1
          nucatr2=nucatr2+lattice%lvec(index1)%oper(1)%elms(st)**2
          nucatrst=nucatrst+st*lattice%lvec(index1)%oper(1)%elms(st)**2
       ENDDO
      ENDDO
      write(lupri,'(A27)') 'DEBUGPBC nuclear attraction'
      write(lupri,'(A10,X,I3)') 'iteration:',iter
      write(lupri,'(A9,X,I3,x,I3,X,I3)') 't1 t2 t3:',il11,il12,il13
      write(lupri,'(A4,X,E16.8)') 'Z^2:',nucatr2
      write(lupri,'(A10,X,E16.8)') 'sum ijZ^2:',nucatrst
#else

      if(lattice%compare_elmnts .and. iter .eq. 1) then
      !compare integrals with the old pbc code

      !  write(lupri,*) 'comparing overlap elements with old pbc code'
      !  write(lupri,*) il1,il2,il3
        
      iunit=-1
      write(numstr1,'(I5)')  il11
      write(numstr2,'(I5)')  il12
      write(numstr3,'(I5)')  il13
      numstr1=adjustl(numstr1)
      numstr2=adjustl(numstr2)
      numstr3=adjustl(numstr3)
      filename='minham'//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'
!!    call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
      CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
      !DO j=1,nbast
      !   write(iunit,*) (lattice%lvec(index1)%fck_vec(i+(j-1)*nbast),i=1,nbast)
      !ENDDO
      call mat_print(f_1(index1),1,H%nrow,1,H%ncol,iunit)
      call lsclose(iunit,'KEEP')

      filename='minnuc'//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'
!!    call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
      CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
      DO j=1,nbast
         write(iunit,*) (lattice%lvec(index1)%oper(1)%elms(i+(j-1)*nbast),i=1,nbast)
      ENDDO
      call lsclose(iunit,'KEEP')


      endif
#endif

      ! phase=phase*2.*pi
      !ToDo Change fock_mtx to matrix type
      ! call mem_dealloc(htemp)
      !call mat_print(H,1,H%nrow,1,H%ncol,6)
!     write(lupri,*) 'H%elms', H%elms
!     if(lattice%compare_elmnts) then
!     !comparing integrals with the old pbc code
!     write(lupri,*) 'comparing nucleus-electron elements with old pbc code'
!     matris(:,:) =0.0
!
!     write(lupri,*) il11,il12,il13

      call mat_daxpy(1.D0,lattice%lvec(index1)%oper(1),f_1(index1))
!
      if(nfdensity(index1)%init_magic_tag.EQ.mat_init_magic_value) THEN
        if(il11**2+il12**2+il13**2 .gt. 0) then
          E_en=E_en+2.*mat_dotproduct(lattice%lvec(index1)%oper(1),nfdensity(index1))
        else
          E_en=E_en+mat_dotproduct(lattice%lvec(index1)%oper(1),nfdensity(index1))
        endif
      endif
#ifdef DEBUGPBC
      write(lupri,*) 'Hamilton',il11
      call mat_print(f_1(index1),1,H%nrow,1,H%ncol,lupri)
#endif

!     call pbc_readopmat2(il11,il12,il13,matris,nbast,'NUCVNF2',.true.,.false.)
!     call write_matrix(matris,nbast,nbast)
!     write(lupri,*) ''
!     write(lupri,*) 'To compare with the matrix below'
!     CALL mat_print(H,1,H%nrow,1,H%ncol,6)
!     endif
!     write(*,*) 'fock(nlat)', index
!     write(*,*)  
!     call write_matrix(lattice%lvec(index1)%fck_vec,nbast,nbast)
!     call mat_daxpy(1.D0,lattice%lvec(index1)%oper(1),lattice%lvec(index1)%oper(3))
      call pbc_get_file_and_write(lattice,nbast,nbast,index1,3,1,'            ')!3 refers to nuc-el

      call mat_free(lattice%lvec(index1)%oper(1))
      call mat_zero(H)
  endif
endif ! is_redundant

  ENDDO

#ifdef DEBUGPBC
 write(lupri,*) 'Nuclear attraction energy', E_en
 write(*,*) 'Nuclear attraction energy', E_en
#endif

!     do m=1,bz%nk
!     write(*,*) 'fock(k)', m
!     write(*,*) 
!     DO j=1,nbast
!       write(*,*) (nucel_k(m,(i+(j-1)*nbast)),i=1,nbast)
!     ENDDO
!     enddo
!     stop
    call mat_free(H)
    call mat_free(H1)

  call set_lstime_print(.true.)
  write(lupri,*) 'finished with pbc_nucattrc_k'

END SUBROUTINE pbc_nucattrc_k
!
!
SUBROUTINE pbc_electron_rep_k(lupri,luerr,setting,molecule,nbast,&
 lattice,latt_cell,refcell,numvecs,nfdensity,g_2,E_J)
  IMPLICIT NONE

  TYPE(matrix),intent(in),DIMENSION(numvecs) :: nfdensity
  TYPE(matrix),target,intent(inout) :: g_2(numvecs)
  INTEGER, INTENT(IN) :: lupri, luerr ,nbast,numvecs  
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule,refcell
  TYPE(moleculeinfo), INTENT(IN), DIMENSION(numvecs) :: latt_cell
  !local
  REAL(realk) :: kvec(3)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(MATRIX),pointer  :: F(:),F_tmp(:)
  REAL(realk),INTENT(INOUT) :: E_J
  ! local variables
  REAL(realk) :: std_vec_length
  INTEGER :: i,j,k,il21,il22,il23,m,iunit,s,t,st
  INTEGER :: il31,il32,il33, newcell,il1,il2,il3,ir1,ir2,ir3
  INTEGER :: index1,index2,index3,checknf,maxl1,maxl2,maxl3
  INTEGER :: l1,l2,l3, num_latvectors, natoms,indred
  INTEGER, DIMENSION(3) :: fdim
  INTEGER(short) :: gab1,gab2,gabmaxsum
  REAL(realk) :: latt_vec_std(3),origin(3),coulombn2,coulombnst
  REAL(realk) :: phase1,phase2,phase3,valmax,E_nfC
  TYPE(lvec_list_t),intent(inout) ::lattice
  COMPLEX(complexk) :: phase
  INTEGER,SAVE :: iter=0
  INTEGER(short) :: valm1, nf1,nf2,nf3
  CHARACTER(len=10) :: numstr0,numstr1,numstr2,numstr3
  CHARACTER(len=40) :: filename
  CHARACTER(len=12) :: diis,stiter
  real(realk) :: TS,TE
  LOGICAL :: compute_int

  write(lupri,*) 'Starting pbc_coul_k'
  iter=iter+1
  call set_lstime_print(.false.)


  natoms=molecule%natoms

  if(iter .gt. 1) then
    write(stiter,'(I5)') iter-1
    stiter=adjustl(stiter)
    diis='diis_'//trim(stiter)//'_'
  endif
  
  call mem_alloc(F_tmp,1)
  call mem_alloc(F,1)

  !write(lupri,*) 'check, nr.'
  compute_int=.false.

  call mat_init(F(1),nbast,nbast)
  call mat_zero(F(1))
  call mat_init(F_tmp(1),nbast,nbast)
  call mat_zero(F_tmp(1))

  origin(1:3)=1D0
  latt_vec_std=0.0
 

  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)

!  write(*,*) 'COULOMB mat computing 0'

  maxl1=0
  maxl2=0
  maxl3=0
  E_nfC=0._realk
  DO index1=1,num_latvectors

!  write(*,*) 'COULOMB mat computing 1'
  !if(.not. lattice%lvec(index1)%is_redundant) then
  !call find_latt_vectors(index1,il1,il2,il3,fdim,lattice)
  il1=int(lattice%lvec(index1)%lat_coord(1))
  il2=int(lattice%lvec(index1)%lat_coord(2))
  il3=int(lattice%lvec(index1)%lat_coord(3))
  
!  CALL LSTIMER('START',TS,TE,LUPRI)
  lattice%lvec(index1)%g2_computed=.false.
  lattice%lvec(index1)%J_computed=.false.
  gab1=lattice%lvec(index1)%maxgab

  if(gab1  .lt. -18) CYCLE
  if(lattice%lvec(index1)%is_redundant) then
    call find_latt_index(indred,-il1,-il2,-il3,fdim,lattice,lattice%max_layer)
    if(lattice%lvec(indred)%g2_computed) then
      lattice%lvec(index1)%g2_computed=.true.
      lattice%lvec(index1)%J_computed=.true.
    endif
  endif

  if(.not. lattice%lvec(index1)%is_redundant) then
    call find_latt_index(indred,-il1,-il2,-il3,fdim,lattice,lattice%max_layer)
    if(il1**2+il2**2+il3**2 .gt. 0) lattice%lvec(indred)%is_redundant =.true.

    DO index2=1,num_latvectors

       !call find_latt_vectors(index2,il21,il22,il23,fdim,lattice)
       il21=int(lattice%lvec(index2)%lat_coord(1))
       il22=int(lattice%lvec(index2)%lat_coord(2))
       il23=int(lattice%lvec(index2)%lat_coord(3))

       gab2=lattice%lvec(index2)%maxgab

       if(gab2  .lt. -18) CYCLE

      if(abs(il21) .gt. lattice%ndmat) CYCLE !then 
      if( abs(il22) .gt.lattice%ndmat) CYCLE
      if(abs(il23) .gt. lattice%ndmat) CYCLE !then 
      if(nfdensity(index2)%init_magic_tag.NE.mat_init_magic_value) CYCLE

      !if(abs(il21) .le. lattice%nf .and. abs(il22) .le. lattice%nf) THEN
      !  if(abs(il23) .le. lattice%nf) THEN
!      write(*,*) 'COULOMB mat computing 2'


     DO index3=1,num_latvectors


      !call find_latt_vectors(index3,il31,il32,il33,fdim,lattice)
      il31=int(lattice%lvec(index3)%lat_coord(1))
      il32=int(lattice%lvec(index3)%lat_coord(2))
      il33=int(lattice%lvec(index3)%lat_coord(3))

      if(abs(il31) .le. lattice%nf .and. abs(il32) .le. lattice%nf) THEN
        if(abs(il33) .le. lattice%nf) THEN

      !if(abs(il31) .gt. lattice%ndmat) CYCLE !then 
      !if( abs(il32) .gt.lattice%ndmat) CYCLE
      !if(abs(il33) .gt. lattice%ndmat) CYCLE !then 
      !if(nfdensity(index3)%init_magic_tag.NE.mat_init_magic_value) CYCLE

!     ToDo Remove: base on CS screening at some point

      l1=il21+il31
      if( abs(l1) .gt.lattice%max_layer) CYCLE
      
      l2=il22+il32
      if( abs(l2) .gt.lattice%max_layer) CYCLE
      
      l3=il23+il33
      if( abs(l3) .gt.lattice%max_layer) CYCLE
      
      call find_latt_index(newcell,l1,l2,l3,fdim,lattice,lattice%max_layer)

!      write(*,*) 'COULOMB mat computing 3'
      !call find_latt_index(checknf,il31,il32,il33,&
                            !fdim,lattice,lattice%max_layer)

      !gab2=lattice%lvec(index3)%maxgab

      !if(gab2  .lt. -18) CYCLE
      
      ! ToDo Should be turned on at some point
      
      !nf1=lattice%nf
      !nf2=lattice%nf
      !nf3=lattice%nf
      !if(il1 .lt. 0)nf1=-lattice%nf
      !if(il2 .lt. 0)nf2=-lattice%nf
      !if(il3 .lt. 0)nf3=-lattice%nf

      !if(il21+il1 .gt. lattice%nf .and. il1 .gt. 0 ) cycle
      !if(il21 .lt. nf1 .and. il1 .lt. 0) CYCLE
      !if(il21 .lt. nf1

         gabmaxsum=gab1+gab2
         call mat_abs_max_elm(nfdensity(index2),valmax)
         if(valmax .gt. 0._realk) valm1=int(log10(valmax),kind=short)

!         write(*,*) 'COULOMB mat computing 4'
         !if(valm1  .lt. -18) CYCLE
         if(valm1+gabmaxsum  .ge. -12) Then
         !if(gabmaxsum +valm1 .ge. -10) Then
!          write(lupri,*) 'computing coul before gabmaxsum',il1,il2,il3,gabmaxsum
         !if(gabmaxsum  .ge. -12) Then
           compute_int=.true.
           lattice%lvec(index1)%g2_computed=.true.
           lattice%lvec(index1)%J_computed=.true.
           !lattice%lvec(indred)%is_redundant=.true.
           !lattice%lvec(indred)%g2_computed=.true.
          !call TYPEDEF_setmolecules(setting,latt_cell(index2),3,latt_cell(newcell),4)
          !call  II_get_maxGabelm_ScreenMat(lupri,luerr,setting,Gab2)
           if(lattice%lvec(index1)%oper(2)%init_magic_tag.NE.mat_init_magic_value) then
              call mat_init(lattice%lvec(index1)%oper(2),nbast,nbast)
              call mat_zero(lattice%lvec(index1)%oper(2))
           endif
 !          if(lattice%lvec(indred)%oper(2)%init_magic_tag.NE.mat_init_magic_value) then
 !             call mat_init(lattice%lvec(indred)%oper(2),nbast,nbast)
 !             call mat_zero(lattice%lvec(indred)%oper(2))
 !          endif
 !          lattice%lvec(indred)%g2_computed=.true.
 !          lattice%lvec(indred)%J_computed=.true.


          !if(lattice%lvec(indred)%oper(2)%init_magic_tag.NE.mat_init_magic_value) then
          !    call mat_init(lattice%lvec(indred)%oper(2),nbast,nbast)
          !    call mat_zero(lattice%lvec(indred)%oper(2))
          !endif

          !setting%samefrag=.true.
          maxl1=max(il1,maxl1)
          maxl2=max(il2,maxl2)
          maxl3=max(il3,maxl3)
#ifdef DEBUGPBC
          write(lupri,*) 'computing coul mat for',il1,il2,il3
#endif

          call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1),&
               2,latt_cell(index3),3,latt_cell(newcell),4)

#ifdef DEBUGPBC
write(lupri,*) 'TYPEDEF coul mat for',il1,il2,il3
#endif


          call mat_zero(F_tmp(1))
          call II_get_coulomb_mat(lupri,luerr,setting,nfdensity(index2:index2),&
               &F_tmp,1)
#ifdef DEBUGPBC
write(lupri,*) 'computed coul mat for',il1,il2,il3
#endif

          !call mat_daxpy(0.5_realk,F_tmp(1),F(1))
          call mat_daxpy(0.5_realk,F_tmp(1),lattice%lvec(index1)%oper(2))
         endif
       endif
     endif

#ifdef DEBUGPBC
      write(lupri,*) 'computing coul mat for l3',il31
#endif

     ENDDO
#ifdef DEBUGPBC
      write(lupri,*) 'computing coul mat for l2',il21
#endif
!      endif
!    endif
    ENDDO

 !   if(il1 .ne. 0 .or. il2 .ne. 0 .or. il3 .ne. 0) then
 !     call mat_trans(lattice%lvec(index1)%oper(2),lattice%lvec(indred)%oper(2))
 !     lattice%lvec(indred)%is_redundant=.true.
 !   endif

 ! endif !is_redundant


#ifdef DEBUGPBC
  if(lattice%lvec(index1)%g2_computed) write(lupri,*) 'COULOMB MAT',il1
  if(lattice%lvec(index1)%g2_computed) call mat_print(lattice%lvec(index1)%oper(2),1,nbast,1,nbast,lupri)
#endif
   

#ifdef DEBUGPBC
     if(lattice%lvec(index1)%g2_computed) then
       st=0
       Coulombn2=0._realk
       Coulombnst=0._realk
        DO s=1,nbast
         DO t=1,nbast
          st=st+1
          Coulombn2=Coulombn2+lattice%lvec(index1)%oper(2)%elms(st)**2
          Coulombnst=Coulombnst+st*lattice%lvec(index1)%oper(2)%elms(st)**2
         ENDDO
        ENDDO
        write(lupri,'(A26)') 'DEBUGPBC nearfield coulomb'
        write(lupri,'(A10,X,I3)') 'iteration:',iter
        write(lupri,'(A9,X,I3,x,I3,X,I3)') 't1 t2 t3:',il1,il2,il3
        write(lupri,'(A4,X,E16.8)') 'Z^2:',coulombn2
        write(lupri,'(A10,X,E16.8)') 'sum ijZ^2:',coulombnst
     endif
#else

     if(lattice%compare_elmnts ) then
       !compare integrals with the old pbc code

     !  write(lupri,*) 'comparing overlap elements with old pbc code'
     !  write(lupri,*) il1,il2,il3
     if(lattice%lvec(index1)%g2_computed) then
       
       iunit=-1
       write(numstr0,'(I5)')  iter
       write(numstr1,'(I5)')  il1
       write(numstr2,'(I5)')  il2
       write(numstr3,'(I5)')  il3
       numstr0=adjustl(numstr0)
       numstr1=adjustl(numstr1)
       numstr2=adjustl(numstr2)
       numstr3=adjustl(numstr3)
       filename='minjop'//trim(numstr0)//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'
!!       call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
       CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
       DO j=1,nbast
          write(iunit,*) (lattice%lvec(index1)%oper(2)%elms(i+(j-1)*nbast),i=1,nbast)
       ENDDO
       call lsclose(iunit,'KEEP')
      endif

     endif !comparing
#endif

     if(lattice%lvec(index1)%J_computed) Then
       if(.not. lattice%store_mats) then
         call mat_init(g_2(index1),nbast,nbast)
         call mat_copy(1._realk,lattice%lvec(index1)%oper(2),g_2(index1))
       endif
      ! DO k=1,nbast*nbast
      !    lattice%lvec(index1)%fck_vec(k)=lattice%lvec(index1)%fck_vec(k)&
      !    +lattice%lvec(index1)%oper(2)%elms(k)!+0.5d0*Kx(1)%elms(k)
      ! ENDDO
     ! ToDo mat_to_full

!    call mat_daxpy(1.0_realk,lattice%lvec(index1)%oper(2),lattice%lvec(index1)%oper(3))

      if((abs(il1) .le. lattice%ndmat .and. abs(il2) .le. lattice%ndmat)&
          .and.abs(il3) .le. lattice%ndmat)then
        !call mat_max_elm(lattice%lvec(index1)%oper(2),valmax)
        !write(*,*) 'max JOP element for lx',valmax,il1
       if(nfdensity(index1)%init_magic_tag.EQ.mat_init_magic_value) THEN
         if(il1**2+il2**2+il3**2 .gt. 0) then
           E_nfC=E_nfC+&
           1.0*mat_dotproduct(lattice%lvec(index1)%oper(2),nfdensity(index1))
         else
           E_nfC=E_nfC+&
           0.5*mat_dotproduct(lattice%lvec(index1)%oper(2),nfdensity(index1))
         endif
       endif
      endif
     endif


     !endif
     !if(lattice%lvec(index1)%g2_computed) write(lupri,*) 'Coulomb mat for layer',il1,il2,il3,index1
     compute_int=.false.
     endif !is_redundant
     !lattice%lvec(index1)%is_redundant=.false.
  !  if(iter .eq. 1) then
  !    call pbc_get_file_and_write(lattice,nbast,nbast,index1,4,2)
   ! else
   !   call pbc_get_file_and_write(lattice,nbast,nbast,index1,4,diis)
  !  endif
  
    !CALL mat_print(F(1),1,F_tmp(1)%nrow,1,F_tmp(1)%ncol,6)
    call mat_zero(F(1))

!    write(lupri,*) 'index1',index1, num_latvectors
  lattice%lvec(index1)%is_redundant=.false.
  ENDDO
  lattice%col1=maxl1
  lattice%col2=maxl2
  lattice%col3=maxl3

  call mat_free(F(1))
  call mat_free(F_tmp(1))
 ! do l1=1,nfsze
 ! call mat_free(D(l1))
 ! ENDDO
  call mem_dealloc(F)
  call mem_dealloc(F_tmp)
  E_J=E_nfC
#ifdef DEBUGPBC
 write(lupri,*) 'Electronic repulsion nearfield energy', E_J
 write(*,*) 'Electronic repulsion nearfield energy', E_J
#endif
  !deallocate(D)
  write(lupri,*) 'Check, nr. lattice vectors',num_latvectors
  call set_lstime_print(.true.)
  write(lupri,*) 'finished with pbc_elrep_k'
  !if(iter .eq. 3) STOP

END SUBROUTINE pbc_electron_rep_k
!
!
!
SUBROUTINE pbc_exact_xc_k(lupri,luerr,setting,molecule,nbast,&
lattice,latt_cell,refcell,numvecs,nfdensity,g_2,E_K)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr ,nbast,numvecs  
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule,refcell
  TYPE(matrix),intent(in) :: nfdensity(numvecs)
  TYPE(matrix),intent(inout) :: g_2(numvecs)
!  TYPE(lvec_data_t),intent(inout) :: nfdensity(nfsze)
  TYPE(moleculeinfo), INTENT(IN), DIMENSION(numvecs) :: latt_cell
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(MATRIX),pointer  ::  K_tmp(:)
  TYPE(MATRIX)  ::  Kx
  TYPE(lvec_list_t),intent(inout) ::lattice
  REAL(realk),INTENT(INOUT) :: E_K

  ! local variables
  !TYPE(moleculeinfo), pointer :: lattice_cell(:)
  !TYPE(lattice_cell_info_t),pointer :: lat_cells(:)
  REAL(realk) :: std_vec_length,exchang2,exchangst
  INTEGER :: i,j,k,il21,il22,il23,m
  INTEGER :: il31,il32,il33, newcell,il1,il2,il3,checknf1,checknf2,checknf3
  INTEGER :: index1,index2,index3,checknf,gabind
  INTEGER :: l1,l2,l3, num_latvectors, natoms,iunit,maxl1,maxl2,maxl3
  INTEGER :: gabl1,gabl2,gabl3,s,t,st,indred
  INTEGER, DIMENSION(3) :: fdim
  integer(short) :: gab1,gab2,maxgabsum
  REAL(realk) :: latt_vec_std(3),origin(3)!,debug_mat(nbast,nbast)
  REAL(realk) :: phase1,phase2,phase3,kvec(3),E_hfx
!  real(realk) :: pi=3.14159265
  COMPLEX(complexk) :: phase
  INTEGER, SAVE :: iter=0
  integer(short) :: valm1
  real(realk) :: TS,TE,valmax,val
  CHARACTER(len=10) :: numstr0,numstr1,numstr2,numstr3
  CHARACTER(len=40) :: filename
  CHARACTER(len=12) :: diis,stiter
  LOGICAL :: compute_int
  iter=iter+1

  natoms=molecule%natoms
  call set_lstime_print(.false.)
  write(lupri,*) 'check, nr.'

  if(iter .gt. 1) then
    write(stiter,'(I5)') iter-1
    stiter=adjustl(stiter)
    diis='diis_'//trim(stiter)//'_'
  endif

  call mem_alloc(K_tmp,1)

  call mat_init(Kx,nbast,nbast)
  call mat_zero(Kx)
  call mat_init(K_tmp(1),nbast,nbast)
  call mat_zero(K_tmp(1))

  origin(1:3)=1D0
  latt_vec_std=0.0
 

  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)

  maxl1=0
  maxl2=0
  maxl3=0

  compute_int=.false.
 
  E_hfx=0._realk
  DO index1=1,num_latvectors

  !call find_latt_vectors(index1,il1,il2,il3,fdim,lattice)
  il1=int(lattice%lvec(index1)%lat_coord(1))
  il2=int(lattice%lvec(index1)%lat_coord(2))
  il3=int(lattice%lvec(index1)%lat_coord(3))
  !if(abs(il1) .gt. lattice%kx1) CYCLE
  !if(abs(il2) .gt. lattice%kx2) CYCLE
  !if(abs(il3) .gt. lattice%kx3) CYCLE
!  CALL LSTIMER('START',TS,TE,LUPRI)
  !call mat_init(lattice%lvec(index1)%oper(1),nbast,nbast)
  !call mat_zero(lattice%lvec(index1)%oper(1))
  !
  !lattice%lvec(index1)%Kx_computed=.false.
  lattice%lvec(index1)%Kx_computed=.false.

  if(lattice%lvec(index1)%is_redundant) then
    call find_latt_index(indred,-il1,-il2,-il3,fdim,lattice,lattice%max_layer)
    if(lattice%lvec(indred)%Kx_computed) then
      lattice%lvec(index1)%g2_computed=.true.
      lattice%lvec(index1)%Kx_computed=.true.
    endif
  endif

  if(.not. lattice%lvec(index1)%is_redundant) then
    call find_latt_index(indred,-il1,-il2,-il3,fdim,lattice,lattice%max_layer)

    DO index2=1,num_latvectors

       !call find_latt_vectors(index2,il21,il22,il23,fdim,lattice)
       il21=int(lattice%lvec(index2)%lat_coord(1))
       il22=int(lattice%lvec(index2)%lat_coord(2))
       il23=int(lattice%lvec(index2)%lat_coord(3))
       gab1=lattice%lvec(index2)%maxgab
       if(gab1 .lt. -18) CYCLE

       !Skal kun ha avstand til origo ikke latt_std_vec her
       !call calc_distance(distance,lattice%lvec(index2)%std_coord,latt_vec_std)


     DO index3=1,num_latvectors


      !call find_latt_vectors(index3,il31,il32,il33,fdim,lattice)
      il31=int(lattice%lvec(index3)%lat_coord(1))
      il32=int(lattice%lvec(index3)%lat_coord(2))
      il33=int(lattice%lvec(index3)%lat_coord(3))

      if(abs(il31) .gt. lattice%ndmat) CYCLE! then
      if(abs(il32) .gt. lattice%ndmat) CYCLE !then 
      if(abs(il33) .gt. lattice%ndmat) CYCLE !then 
      if(nfdensity(index3)%init_magic_tag.NE.mat_init_magic_value) CYCLE

!     ToDo Remove: base on CS screening at some point
      l1=il21+il31
      if( abs(l1) .gt.lattice%max_layer) CYCLE
      !l1 blir større enn max slik at newcell blir negativ
      !sjekk dette, er ikke sikker på om l1 skal være slik.
      !Likedan man l2 og l3
      
      l2=il22+il32
      if( abs(l2) .gt.lattice%max_layer) CYCLE
      
      l3=il23+il33
      if( abs(l3) .gt.lattice%max_layer) CYCLE
      gabl1=l1-il1
      gabl2=l2-il2
      gabl3=l3-il3
      if( abs(gabl1) .gt.lattice%max_layer) CYCLE
      if( abs(gabl2) .gt.lattice%max_layer) CYCLE
      if( abs(gabl3) .gt.lattice%max_layer) CYCLE


      call find_latt_index(newcell,l1,l2,l3,fdim,lattice,lattice%max_layer)
      !call find_latt_index(checknf,il31,il32,il33,&
      !                      fdim,lattice,lattice%max_layer)

      call find_latt_index(gabind,gabl1,gabl2,gabl3,fdim,lattice,lattice%max_layer)

      gab2=lattice%lvec(gabind)%maxGab
      if(gab2 .lt. -18) CYCLE
      maxgabsum=gab1+gab2
      call mat_abs_max_elm(nfdensity(index3),valmax)
      if(valmax.gt.0._realk) valm1=int(log10(valmax),kind=short)


      !if(maxgabsum+valm1 .ge. -10) Then
      if(maxgabsum+valm1 .ge. -10) Then
        !write(lupri,*) 'maxgabsum',maxgabsum,valm1,valmax,il31
        if(lattice%lvec(index1)%oper(1)%init_magic_tag.NE.mat_init_magic_value) then
          call mat_init(lattice%lvec(index1)%oper(1),nbast,nbast)
          call mat_zero(lattice%lvec(index1)%oper(1))
        endif
        !write(lupri,*) 'gab1,gab2',gab1,gab2
        compute_int=.true.
        lattice%lvec(index1)%g2_computed=.true.
        lattice%lvec(index1)%Kx_computed=.true.
        maxl1=max(abs(il1),maxl1)
        maxl2=max(abs(il2),maxl2)
        maxl3=max(abs(il3),maxl3)

        call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1),3,latt_cell(index2),2,latt_cell(newcell),4)

        call mat_zero(K_tmp(1))
        call II_get_exchange_mat(lupri,luerr,setting,nfdensity(index3:index3),&
             &                   1,.false.,K_tmp)

        call mat_daxpy(1.D0,K_tmp(1),Kx)
        call mat_daxpy(0.5D0,K_tmp(1),lattice%lvec(index1)%oper(1))

        if(lattice%lvec(indred)%oper(1)%init_magic_tag.NE.mat_init_magic_value) then
          call mat_init(lattice%lvec(indred)%oper(1),nbast,nbast)
          call mat_zero(lattice%lvec(indred)%oper(1))
        endif
      endif

     ENDDO
    ENDDO

    if((il1**2 + il2**2 +il3**2 .gt. 0) .and. lattice%lvec(index1)%Kx_computed) then
      call mat_trans(lattice%lvec(index1)%oper(1),lattice%lvec(indred)%oper(1))
      lattice%lvec(indred)%is_redundant=.true.
    endif

    if(lattice%lvec(index1)%Kx_computed) write(lupri,*) 'Kx computed for', il1,il2,il3

  endif !is_redundant
 !    


!  CALL LSTIMER('xch one cell',TS,TE,LUPRI)

     if(lattice%compare_elmnts ) then
       !compare integrals with the old pbc code

       if(lattice%lvec(index1)%Kx_computed)then
       
       iunit=-1
       write(numstr0,'(I5)')  iter
       write(numstr1,'(I5)')  il1
       write(numstr2,'(I5)')  il2
       write(numstr3,'(I5)')  il3
       numstr0=adjustl(numstr0)
       numstr1=adjustl(numstr1)
       numstr2=adjustl(numstr2)
       numstr3=adjustl(numstr3)
       filename='minxch'//trim(numstr0)//trim(numstr1)//trim(numstr2)//trim(numstr3)//'.dat'
!!       call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
       CALL lsOPEN(IUNIT,filename,'unknown','FORMATTED')
       DO j=1,nbast
          write(iunit,*) (lattice%lvec(index1)%oper(1)%elms(i+(j-1)*nbast),i=1,nbast)
       ENDDO
       call lsclose(iunit,'KEEP')
     endif
   !if(checknf1 .eq. 0 .and. checknf2 .eq. 0 .and. checknf3 .eq. 0) then
     !call mat_print(nfdensity(checknf),1,nbast,1,nbast,6)
   !endif

     endif !comparing

  !   DO k=1,nbast*nbast
  !      lattice%lvec(index1)%fck_vec(k)=lattice%lvec(index1)%fck_vec(k)&
  !      +0.5*Kx%elms(k)
  !   ENDDO
     ! ToDo mat_to_full

     ! FOR debug only
#ifdef DEBUGPBC
     write(lupri,*) 'Exchange mat', il1
     if(lattice%lvec(index1)%oper(1)%init_magic_tag.EQ.mat_init_magic_value) THEN
       call mat_print(lattice%lvec(index1)%oper(1),1,nbast,1,nbast,lupri)
     endif
#endif


    !if(iter .eq. 1) then
    if(lattice%lvec(index1)%Kx_computed) then
#ifdef DEBUGPBC
       st=0
       exchang2=0._realk
       exchangst=0._realk
        DO s=1,nbast
         DO t=1,nbast
          st=st+1
          exchang2=exchang2+lattice%lvec(index1)%oper(1)%elms(st)**2
          exchangst=exchangst+st*lattice%lvec(index1)%oper(1)%elms(st)**2
         ENDDO
        ENDDO
        write(lupri,'(A23)') 'DEBUGPBC exact exchange'
        write(lupri,'(A10,X,I3)') 'iteration:',iter
        write(lupri,'(A9,X,I3,x,I3,X,I3)') 't1 t2 t3:',il1,il2,il3
        write(lupri,'(A4,X,E16.8)') 'Z^2:',exchang2
        write(lupri,'(A10,X,E16.8)') 'sum ijZ^2:',exchangst
#endif
      !call mat_abs_max_elm(lattice%lvec(index1)%oper(1), val) 
      !if(val .lt. 1E-10) write(*,*)'max K0l',val,il1
      if(lattice%store_mats) then
        call pbc_get_file_and_write(lattice,nbast,nbast,index1,5,1,'            ')!5 refers to
      else
        if(g_2(index1)%init_magic_tag.NE.mat_init_magic_value) then
          call mat_init(g_2(index1),nbast,nbast)
          call mat_zero(g_2(index1))
        endif
        call mat_daxpy(1.D0,lattice%lvec(index1)%oper(1),g_2(index1))
      endif
      if((abs(il1) .le. lattice%ndmat .and. abs(il2) .le. lattice%ndmat)&
        .and. abs(il3) .le. lattice%ndmat)then
        if(nfdensity(index1)%init_magic_tag.EQ.mat_init_magic_value) THEN
          E_hfx=E_hfx+0.5*mat_dotproduct(lattice%lvec(index1)%oper(1),nfdensity(index1))
        endif
      endif

    endif
    compute_int=.false.
    !lattice%lvec(index1)%is_redundant=.false.
    !else                                            !exact exhange 
    !  call pbc_get_file_and_write(lattice,nbast,nbast,index1,5,1,diis)
    !endif

    !CALL mat_print(lattice%lvec(index1)%oper(1),1,nbast,1,nbast,lupri)
    if(lattice%lvec(index1)%oper(1)%init_magic_tag.EQ.mat_init_magic_value) then
      call mat_free(lattice%lvec(index1)%oper(1))
    endif
    !CALL mat_print(F(1),1,F_tmp(1)%nrow,1,F_tmp(1)%ncol,6)
    call mat_zero(Kx)

!    write(lupri,*) 'index1',index1, num_latvectors
  ENDDO
  lattice%Kx1=maxl1
  lattice%Kx2=maxl2
  lattice%Kx3=maxl3

  call mat_free(Kx)
  call mat_free(K_tmp(1))
 ! do l1=1,nfsze
 ! call mat_free(D(l1))
 ! ENDDO
!  call mem_dealloc(Kx)
  call mem_dealloc(K_tmp)

  E_K=E_hfx

#ifdef DEBUGPBC
 write(lupri,*) 'Exact exchange energy', E_K
 write(*,*) 'Exact exchange energy', E_K
#endif

  write(lupri,*) 'Check, nr. lattice vectors',num_latvectors
  call set_lstime_print(.true.)
  write(lupri,*) 'finished with pbc_exact_xc_k'
!  STOP

END SUBROUTINE pbc_exact_xc_k


!Maybe not needed
SUBROUTINE pbc_overlap_int(lupri,luerr,setting,molecule,nbast,lattice,latt_cell,refcell,numvecs)
  implicit none
  INTEGER, INTENT(IN) :: lupri, luerr, nbast,numvecs
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(moleculeinfo),intent(inout)          :: refcell
  TYPE(moleculeinfo),intent(in) :: latt_cell(numvecs)

  ! local variables
  TYPE(MATRIX)  :: S
  INTEGER :: i,il1,il2,il3
  INTEGER ::  index,refindex
  INTEGER ::  natoms,origin(3)
  INTEGER, DIMENSION(3) :: fdim
  REAL(realk) :: latt_vec_std(3)
  TYPE(lvec_list_t),intent(INOUT) ::lattice


  natoms=molecule%natoms

  ! II_get_sphmom

  call mat_init(S,nbast,nbast)
  call mat_zero(S)
  origin(1:3)=1.0
  latt_vec_std=0.0

  write(lupri,*) 'Number of lattice vectors ', numvecs

  i=1
  call find_latt_index(refindex,0,0,0,fdim,lattice,lattice%max_layer)
  DO index=1,numvecs
     

     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index),3)
     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.

     call find_latt_vectors(index,il1,il2,il3,fdim,lattice)
     if(abs(il1) .gt. lattice%nneighbour) CYCLE
     if(abs(il2) .gt. lattice%nneighbour) CYCLE
     if(abs(il3) .gt. lattice%nneighbour) CYCLE
     !input%Basis%regular%atomtype(:)%shell(:)%segment(1)%exponents(:)


     call II_get_overlap(lupri,luerr,setting,S)
     !write(lupri,*) 'S%elms', S%elms!,nbast
     !if(lattice%compare_elmnts) then
     !  !compare integrals with the old pbc code

     !  write(lupri,*) 'comparing overlap elements with old pbc code'
     !  write(lupri,*) il1,il2,il3
!    ! write(lupri,*) 'molecule%p%atom%center', setting%molecule(1)%p%atom(1)%center(1)
!    ! write(lupri,*) 'setting%p%atom%center', setting%molecule(3)%p%atom(1)%center(1)
!
!    !   matris(:,:) =0.0
!    !   call pbc_readopmat2(il1,il2,il3,matris,nbast,'OVERLAP',.true.,.false.)
!    !   call write_matrix(matris,nbast,nbast)
!
!    !   write(lupri,*) ''
!    !   write(lupri,*) 'To compare with the matrix below'
     !  CALL mat_print(S,1,S%nrow,1,S%ncol,6)
     !endif

!     i=1
!     j=0
!     DO k=1,nbast*nbast
!        j=j+1
!        if(j .gt. nbast) Then
!          j=1
!          i=i+1
!        ENDIF
!        tmp_overl(j,i)=S%elms(k)
!     ENDDO
!!     !write(lupri,*) fock_tmp
!     i=0
!     DO k=refindex*nbast-nbast+1,refindex*nbast 
!       i=i+1
!       j=0
!     DO m=index*nbast-nbast+1,index*nbast 
!       j=j+1
!       write(lupri,*) 'k,m,i,j',k,m,i,j
!       overlap(k,m)=tmp_overl(i,j)!*coeff(k,m)
!       overlap(m,k)=overlap(k,m)!dagger if complex
!     ENDDO
!     ENDDO

!     CALL mat_print(S,1,S%nrow,1,S%ncol,6)
!     STOP
     !store the overlap somewhere
  END DO
  call mat_free(S)

  write(lupri,*) 'finished pbc_overlap_int'
!  STOP

END SUBROUTINE pbc_overlap_int

!
SUBROUTINE pbc_kinetic_int(lupri,luerr,setting,molecule,nbast,fock_mtx,sizef,lattice,latt_cell,refcell,numvecs)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr,nbast,sizef,numvecs
  COMPLEX(complexk), INTENT(INOUT) :: fock_mtx(sizef*sizef)
  TYPE(moleculeinfo),INTENT(IN) :: molecule
  TYPE(LSSETTING),intent(inout)   :: SETTING 

  ! local variables
  REAL(realk), DIMENSION(12) :: Tvec
  TYPE(moleculeinfo)  :: refcell
  TYPE(moleculeinfo), DIMENSION(numvecs),INTENT(IN) :: latt_cell
  TYPE(MATRIX)  :: kin
  REAL(realk) :: std_vec_length
  INTEGER :: i,j,il1,il2,il3
  INTEGER ::  index, refindex
  INTEGER :: num_latvectors, natoms
  INTEGER, DIMENSION(3) :: fdim
  REAL(realk) :: latt_vec_std(3),origin(3)
  REAL(realk) :: fock_tmp(nbast*nbast)
  TYPE(lvec_list_t),intent(INOUT) ::lattice

  write(lupri,*) 'Starting pbc_kinetic_int'
  natoms=molecule%natoms


  origin(1:3)=1.0
  latt_vec_std=0.0

  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)

  fdim(1:3) = 0

  call find_latt_index(refindex,0,0,0,fdim,lattice,lattice%max_layer)
  write(lupri,*) 'reference cell index=', refindex,lattice%max_layer
  fock_tmp= 0.0
  tvec(1:12)=0.0
  i=1

  call mat_init(kin,nbast,nbast)
  call mat_zero(kin)
  write(lupri,*) 'nbast',nbast,kin%ncol,kin%nrow

  DO index=1,num_latvectors
     


     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index),3)

     !call calc_distance(distance,lattice%lvec(index)%std_coord,latt_vec_std)
     call find_latt_vectors(index,il1,il2,il3,fdim,lattice)
     if(abs(il1) .gt. lattice%nneighbour) CYCLE
     if(abs(il2) .gt. lattice%nneighbour) CYCLE
     if(abs(il3) .gt. lattice%nneighbour) CYCLE
     !phase=cmplx(0,k1*il1+k2*il2+k3*il3)
    

     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.
     call II_get_kinetic(lupri,luerr,setting,kin)

!     if(lattice%compare_elmnts) then
!       !compare integrals with the old pbc code
!       write(lupri,*) 'comparing kinetic elements with old pbc code'
!       matris(:,:) =0.0
!       write(lupri,*) il1,il2,il3
!       call pbc_readopmat2(il1,il2,il3,matris,nbast,'KINETIC',.true.,.false.)
!       call write_matrix(matris,nbast,nbast)
!!       call pbc_readopmat2(0,0,0,matris,2,'OVERLAP',.true.,.false.)
!       write(lupri,*) ''
!       write(lupri,*) 'To compare with the matrix below'
!       CALL mat_print(kin1,1,kin1%nrow,1,kin1%ncol,6)
!     endif

     i=1
     j=0
!     call mat_print(kin1,1,kin1%nrow,1,kin1%ncol,6)
!     write(lupri,*) kin1%elms
!     STOP
  !   DO k=1,nbast*nbast
  !      j=j+1
  !      if(j .gt. nbast) Then
  !        j=1
  !        i=i+1
  !      ENDIF
  !      fock_tmp(j,i)=kin1%elms(k)
    ! ENDDO
    fock_tmp=kin%elms!*exp(phase)
     write(lupri,*) fock_tmp
     i=0
   !  DO k=refindex*nbast-nbast+1,refindex*nbast 
   !    i=i+1
   !    j=
   !  DO m=index*nbast-nbast+1,index*nbast 
   !    j=j+1
   !    fock_mtx(k,m)=fock_tmp(i,j)!*coeff(k,m)
   !    fock_mtx(m,k)=fock_mtx(k,m)!dagger if complex
   !  ENDDO
   !  ENDDO
        fock_mtx=fock_tmp

  END DO

  call mat_free(kin)
  write(lupri,*) 'Finished pbc_kinetic_int'

!!
END SUBROUTINE pbc_kinetic_int
!
!
SUBROUTINE pbc_nucattrc_int(lupri,luerr,setting,molecule,nbast,fock_mtx,sizef,lattice,latt_cell,refcell,numvecs)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr,nbast,sizef,numvecs
  TYPE(moleculeinfo),INTENT(IN) :: molecule,refcell
  COMPLEX(complexk), INTENT(INOUT) :: fock_mtx(sizef*sizef)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
!  TYPE(DALTONINPUT) :: INPUT
  TYPE(MATRIX)  :: H,H1

  ! local variables
  REAL(realk), DIMENSION(12) :: Tvec
!  TYPE(moleculeinfo), pointer :: lattice_cell(:)
  TYPE(moleculeinfo), DIMENSION(numvecs),INTENT(IN) :: latt_cell
  REAL(realk) :: std_vec_length
  INTEGER :: il11,il12,il13,il21,il22,il23
  INTEGER ::  index1,index2, refindex
  INTEGER ::  num_latvectors, natoms
  INTEGER, DIMENSION(3) :: fdim
  REAL(realk)::  latt_vec_std(3),origin(3)
  REAL(realk) :: fock_tmp(nbast*nbast)
  TYPE(lvec_list_t),intent(INOUT) ::lattice

  natoms=molecule%natoms


  call mat_init(H,nbast,nbast)
  call mat_zero(H)
  call mat_init(H1,nbast,nbast)
  call mat_zero(H1)

  origin(1:3)=1.0
  latt_vec_std=0.0


  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)
  
  fdim(1:3) = 0
  if (lattice%ldef%is_active(1)) fdim(1) = 1
  if (lattice%ldef%is_active(2)) fdim(2) = 1
  if (lattice%ldef%is_active(3)) fdim(3) = 1


  call find_latt_index(refindex,0,0,0,fdim,lattice,lattice%max_layer)
  write(lupri,*) 'reference cell index=', refindex

  tvec(1:12)=0.0
 

  DO index1=1,num_latvectors



!  call calc_distance(distance,lattice%lvec(index1)%std_coord,latt_vec_std)
  call find_latt_vectors(index1,il11,il12,il13,fdim,lattice)
  if(abs(il11) .gt. lattice%nneighbour) CYCLE
  if(abs(il12) .gt. lattice%nneighbour) CYCLE
  if(abs(il13) .gt. lattice%nneighbour) CYCLE
  !phase=(0,k1*il1+k2*il2+k3*il3)


  DO index2=1,num_latvectors


     call find_latt_vectors(index2,il21,il22,il23,fdim,lattice)
     if(abs(il21) .gt. lattice%nneighbour) CYCLE !uncomment this
     if(abs(il22) .gt. lattice%nneighbour) CYCLE
     if(abs(il23) .gt. lattice%nneighbour) CYCLE

    ! call calc_distance(distance,lattice%lvec(index2)%std_coord,latt_vec_std)

     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1),2,latt_cell(index2),3)

         setting%samemol(1,3)=.false.
         setting%samemol(3,1)=.false.
         setting%samemol(1,2)=.false.
         setting%samemol(2,1)=.false.
         setting%samemol(3,2)=.false.
         setting%samemol(2,3)=.false.

         call II_get_nucel_mat(lupri,luerr,setting,H1)
         call mat_add(1E0_realk,H,1E0_realk,H1,H)
         
!         write(lupri,*) 'H%elms', H%elms


  ENDDO

!     i=1
!     j=0
!     DO k=1,nbast*nbast
!        j=j+1
!        if(j .gt. nbast) Then
!          j=1
!          i=i+1
!        ENDIF
!        fock_tmp(j,i)=H%elms(k)
!     ENDDO
     !call mat_print(H,1,H%nrow,1,H%ncol,6)
     !call write_matrix(fock_tmp,2,2)
     !STOP
     fock_tmp=H1%elms
!     i=0
!     DO k=refindex*nbast-nbast+1,refindex*nbast 
!       i=i+1
!       j=0
!     DO m=index1*nbast-nbast+1,index1*nbast 
!       j=j+1
!       fock_mtx(k,m)=fock_mtx(k,m)+fock_tmp(i,j)!*coeff(k,m)
!       fock_mtx(m,k)=fock_mtx(k,m)!dagger if complex
!     ENDDO
!     ENDDO

     fock_mtx=fock_mtx+fock_tmp!*coeff(k,m)
     call mat_print(H1,1,H1%nrow,1,H1%ncol,6)
!    write(lupri,*) 'H%elms', H%elms
!  if(lattice%compare_elmnts) then
!    !comparing integrals with the old pbc code
!    write(lupri,*) 'comparing nucleus-electron elements with old pbc code'
!    matris(:,:) =0.0
!
!    write(lupri,*) il11,il12,il13
!
!    call pbc_readopmat2(il11,il12,il13,matris,nbast,'NUCVNF2',.true.,.false.)
!    call write_matrix(matris,nbast,nbast)
!    write(lupri,*) ''
!    write(lupri,*) 'To compare with the matrix below'
!    CALL mat_print(H,1,H%nrow,1,H%ncol,6)
!  endif

    call mat_zero(H1)

  ENDDO

  call mat_free(H)
  call mat_free(H1)
  write(lupri,*) 'finished with pbc_nucattrc_int'

END SUBROUTINE pbc_nucattrc_int
!
!
!SUBROUTINE for computing the nuclear repulsions
SUBROUTINE pbc_nucpot(lupri,luerr,setting,molecule,lattice,&
 latt_cell,refcell,numvecs,E_nn)

  implicit none
  ! input and output arguments
  INTEGER, INTENT(IN) :: lupri, luerr, numvecs
  TYPE(moleculeinfo),INTENT(INOUT) :: molecule
  TYPE(moleculeinfo),intent(inout) :: refcell
  TYPE(lvec_list_t),intent(INOUT) ::lattice
  Real(realk),INTENT(INOUT) :: E_nn

  ! local variables
  REAL(realk) :: kvec(3)
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(moleculeinfo),intent(in) :: latt_cell(numvecs)
  INTEGER :: i,j,m,il1,il2,il3
  INTEGER :: index,refindex
  INTEGER :: natoms,iunit,maxl1,maxl2,maxl3
  INTEGER, DIMENSION(3) :: fdim
  REAL(realk) :: latt_vec_std(3)
  REAL(realk) :: nucpot



  natoms=molecule%natoms!number of atoms

  call set_lstime_print(.false.)

  latt_vec_std=0.0

  write(lupri,*) 'Number of lattice vectors ', numvecs

  !finds the lattice index for the reference cell.
  !i.e for l1=l2=l3=0
  call find_latt_index(refindex,0,0,0,fdim,lattice,lattice%max_layer)

  E_nn=0._realk
  !loop over lattice cells

  DO index=1,numvecs
     

  !Doing translations
     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index),3)
     setting%samemol(1,3)=.false.
     setting%samemol(3,1)=.false.

     call find_latt_vectors(index,il1,il2,il3,fdim,lattice)
     !So that we do not consider negligible integrals
     if(abs(il1) .gt. lattice%nf) CYCLE
     if(abs(il2) .gt. lattice%nf) CYCLE
     if(abs(il3) .gt. lattice%nf) CYCLE

     if(index .eq. refindex) then

       call II_get_nucpot(lupri,luerr,setting,nucpot)
     else
       call pbc_get_nucpot(lupri,luerr,setting,nucpot)
     endif
     E_nn=E_nn+nucpot

  END DO
  !E_nn=E_nn/2.

 ! write(*,*) 'Debug 2', E_nn

END SUBROUTINE pbc_nucpot
!> \brief Calculates the nuclear repulsion energy contribution
!> \author S. Reine and T. Kjaergaard
!> \modified by J. Rekkedal
!> \date 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nucpot the nuclear repulsion energy contribution
SUBROUTINE pbc_get_nucpot(LUPRI,LUERR,SETTING,NUCPOT)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
integer             :: usemat
INTEGER               :: LUPRI,LUERR
REAL(realk)           :: nucpot
Integer               :: I,J
real(realk)           :: pq(3),distance
logical               :: NOBQBQ
call time_II_operations1()
NOBQBQ = SETTING%SCHEME%NOBQBQ
NUCPOT=0.0E0_realk
DO I=1,SETTING%MOLECULE(1)%p%Natoms
 IF(SETTING%MOLECULE(1)%p%ATOM(I)%phantom)CYCLE
 DO J=1,SETTING%MOLECULE(1)%p%Natoms
  IF(SETTING%MOLECULE(1)%p%ATOM(J)%phantom)CYCLE
  
  if(setting%molecule(1)%p%ATOM(I)%Pointcharge .and. &
    &setting%molecule(1)%p%ATOM(J)%Pointcharge .and. NOBQBQ) cycle

  pq(1) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)-SETTING%MOLECULE(3)%p%ATOM(J)%CENTER(1)
  pq(2) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)-SETTING%MOLECULE(3)%p%ATOM(J)%CENTER(2)
  pq(3) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)-SETTING%MOLECULE(3)%p%ATOM(J)%CENTER(3)
  Distance = sqrt(pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3))
  NUCPOT=NUCPOT+SETTING%MOLECULE(1)%p%ATOM(I)%Charge*SETTING%MOLECULE(3)%p%ATOM(J)%Charge/Distance
 ENDDO
ENDDO
NUCPOT=NUCPOT/2.
call time_II_operations2(JOB_II_get_nucpot)
END SUBROUTINE pbc_get_nucpot

SUBROUTINE pbc_electron_rep(lupri,luerr,setting,molecule,&
 nbast,lattice,latt_cell,refcell,numvecs,nfdensity,nfsze)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri, luerr ,nbast,numvecs,nfsze
  TYPE(moleculeinfo),INTENT(IN) :: molecule,refcell
  TYPE(LSSETTING),intent(inout)   :: SETTING 
  TYPE(matrix),intent(inout),DIMENSION(nfsze) :: nfdensity
  TYPE(MATRIX)  :: F(1),F_tmp(1)

  ! local variables
  REAL(realk), DIMENSION(12) :: Tvec
  !TYPE(moleculeinfo), pointer :: lattice_cell(:)
  TYPE(moleculeinfo), INTENT(IN), DIMENSION(numvecs) :: latt_cell
  !TYPE(lattice_cell_info_t),pointer :: lat_cells(:)
  REAL(realk) :: std_vec_length
  INTEGER :: il21,il22,il23
  INTEGER :: il31,il32,il33, newcell,il1,il2,il3
  INTEGER ::  index1,index2,index3
  INTEGER :: l1,l2,l3, num_latvectors, natoms
  INTEGER :: checknf,checknf1,checknf2,checknf3
  INTEGER, DIMENSION(3) :: fdim
  REAL(realk):: latt_vec_std(3),origin(3)
  TYPE(lvec_list_t),intent(inout) ::lattice

  natoms=molecule%natoms

!  call mat_init(D(1),nbast,nbast)
!  call mat_zero(D(1))
  call mat_init(F(1),nbast,nbast)
  call mat_zero(F(1))
  call mat_init(F_tmp(1),nbast,nbast)
  call mat_zero(F_tmp(1))

  origin(1:3)=1.0
  latt_vec_std=0.0
 

  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)

  std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  num_latvectors=size(lattice%lvec)


  tvec(1:12)=0.0
 

  DO index1=1,num_latvectors

  !call calc_distance(distance,lattice%lvec(index1)%std_coord,latt_vec_std)
  call find_latt_vectors(index1,il1,il2,il3,fdim,lattice)
!  if(abs(il1) .gt. lattice%nneighbour) CYCLE
!  if(abs(il2) .gt. lattice%nneighbour) CYCLE
!  if(abs(il3) .gt. lattice%nneighbour) CYCLE

  DO index2=1,num_latvectors

     call find_latt_vectors(index2,il21,il22,il23,fdim,lattice)
  if(abs(il21) .gt. lattice%nneighbour) CYCLE
  if(abs(il22) .gt. lattice%nneighbour) CYCLE
  if(abs(il23) .gt. lattice%nneighbour) CYCLE

     !Skal kun ha avstand til origo ikke latt_std_vec her
!     call calc_distance(distance,lattice%lvec(index2)%std_coord,latt_vec_std)


   DO index3=1,num_latvectors

    !write(lupri,*) 'electron rep get'
    !CALL mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,6)
    !CALL mat_print(D1(1),1,F(1)%nrow,1,F(1)%ncol,6)

    call find_latt_vectors(index3,il31,il32,il33,fdim,lattice)
!  if(abs(il31) .gt. lattice%nneighbour) CYCLE
!  if(abs(il32) .gt. lattice%nneighbour) CYCLE
!  if(abs(il33) .gt. lattice%nneighbour) CYCLE

    l1=il21+il31
    if(abs(l1) .gt. lattice%max_layer) CYCLE! then
    
    l2=il22+il32
    if(abs(l2) .gt. lattice%max_layer) CYCLE !then 
    
    l3=il23+il33
    if(abs(l3) .gt. lattice%max_layer) CYCLE !then 
    
    checknf1=l1-il21
    if(abs(checknf1) .gt. n_neighbour) CYCLE! then
    
    checknf2=l2-il22
    if(abs(checknf2) .gt. n_neighbour) CYCLE !then 
    
    checknf3=l3-il23
    if(abs(checknf3) .gt. n_neighbour) CYCLE !then 

    
    call find_latt_index(newcell,l1,l2,l3,fdim,lattice,lattice%max_layer)

    call find_latt_index(checknf,checknf1,checknf2,checknf3,fdim,lattice,&
     lattice%nneighbour)

    !Changes setting to point at different lattice cells
    call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index1),2,&
     latt_cell(index2),3,latt_cell(newcell),4)

    setting%samemol(1,2)=.false.
    setting%samemol(2,1)=.false.
    setting%samemol(1,3)=.false.
    setting%samemol(3,1)=.false.
    setting%samemol(1,4)=.false.
    setting%samemol(4,1)=.false.
    setting%samemol(2,3)=.false.
    setting%samemol(2,3)=.false.
    setting%samemol(3,4)=.false.
    setting%samemol(4,2)=.false.
    setting%samemol(3,4)=.false.
    setting%samemol(4,3)=.false.

 !   setting%samefrag=.false.
    
    call II_get_coulomb_mat(lupri,luerr,setting,nfdensity(checknf:checknf),F_tmp,1)
    !call II_get_exchange_mat(lupri,luerr,setting,D,F_tmp,1)
    
    call mat_add(1E0_realk,F(1),1E0_realk,F_tmp(1),F(1))
    
    !call exchange()
!    write(lupri,*) 'electron rep fin rep',nbast,F(1)%nrow,F(1)%ncol
!    CALL mat_print(F(1),1,F_tmp(1)%nrow,1,F_tmp(1)%ncol,6)
!     STOP
   ENDDO
  ENDDO
  !store the matrix here

!    if(lattice%compare_elmnts) then
!      
!      !compare integrals with the old pbc code
!      write(lupri,*) 'comparing electron-electron elements with old pbc code'
!      matris(:,:) =0.0
!
!      write(lupri,*) il1,il2,il3
!
!      call pbc_readopmat2(il1,il2,il3,matris,nbast,'COULOMB',.true.,.false.)
!      call write_matrix(matris,nbast,nbast)
!
!      write(lupri,*) ''
!      write(lupri,*) 'To compare with the matrix below'
!      CALL mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,6)
!
!    endif

    CALL mat_print(F(1),1,F_tmp(1)%nrow,1,F_tmp(1)%ncol,6)
    call mat_zero(F(1))

    write(lupri,*) 'index1',index1, num_latvectors
  ENDDO
  call mat_free(F(1))
  call mat_free(F_tmp(1))

  write(lupri,*) 'finished with pbc_elrep_int'
  call LSQuit('finished with pbc_elrep_int',lupri)

END SUBROUTINE pbc_electron_rep


SUBROUTINE pbc_complete_Fock_mtx(lupri,nbast,fock_mtx,sizef,cut,lattice)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lupri,nbast,sizef
  REAL(realk), INTENT(INOUT) :: fock_mtx(sizef,sizef)
  REAL(realk), INTENT(IN) :: cut
  ! Local variables
  REAL(realk) :: distance
  TYPE(lvec_list_t), intent(in) ::lattice
  INTEGER :: num_latvectors
  INTEGER :: cellij,i,j
  INTEGER :: diff1,diff2,diff3
  INTEGER :: fdim(3),il1,il2,il3,jl1,jl2,jl3
  REAL(realk) ::  latt_vec_std(3), origin(3)
  INTEGER :: refcell

!  call build_lvec_list(lattice)
  !call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)
  
!  write(lupri,*) 'debug: ',lattice%lvec(2)%std_coord!, latstdvec 

  num_latvectors=size(lattice%lvec)

  origin(1:3)=1.0
  latt_vec_std=0.0


  call latt_2_std_coord(origin,latt_vec_std,lattice%ldef%avec)
  !std_vec_length= sqrt(latt_vec_std(1)**2+latt_vec_std(2)**2+latt_vec_std(3)**2)

  call find_latt_index(refcell,0,0,0,fdim,lattice,lattice%max_layer)
  write(lupri,*) 'refcelle', refcell

  DO i=1,num_latvectors
   IF(i .eq. refcell) CYCLE
   fock_mtx(nbast*i+1-nbast:i*nbast,nbast*i+1-nbast:i*nbast)=&
   fock_mtx(refcell*nbast-nbast+1:refcell*nbast,refcell*nbast-nbast+1:refcell*nbast)
   call find_latt_vectors(i,il1,il2,il3,fdim,lattice)
   DO j=1,i
      IF(j .eq. refcell) CYCLE
      IF(j .eq. i) CYCLE
      call find_latt_vectors(j,jl1,jl2,jl3,fdim,lattice)
      !write(lupri,*) 'SEGMENTATION FAULT',j
!      write(lupri,*) 'refcelle', refcell,j
      diff1=il1-jl1
      if(abs(diff1) .gt. lattice%max_layer) CYCLE !should be max_layer
      if(abs(diff1) .gt. lattice%nneighbour) CYCLE !should be threshold
      diff2=il2-jl2
      if(abs(diff2) .gt.lattice%max_layer) CYCLE !should be max_layer
      if(abs(diff2) .gt.lattice%nneighbour) CYCLE !should be threshold
      diff3=il3-jl3
      if(abs(diff3) .gt. lattice%max_layer) CYCLE !should be max_layer
      if(abs(diff3) .gt. lattice%nneighbour) CYCLE !should be threshold

      call find_latt_index(cellij,diff1,diff2,diff3,fdim,lattice,lattice%max_layer)

      call calc_distance(distance,lattice%lvec(cellij)%std_coord,latt_vec_std)

      IF(distance .ge. cut) CYCLE
      fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)=&
      fock_mtx(refcell*nbast-nbast+1:refcell*nbast,cellij*nbast-nbast+1:cellij*nbast)
      fock_mtx(nbast*j+1-nbast:j*nbast,nbast*i+1-nbast:i*nbast)=&
      fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)
 !     IF(i /= j) THEN
 !       diff= i-j
 !       fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)=&
 !       fock_mtx(refcell*nbast-nbast+1:refcell*nbast,diff*nbast-nbast+1:diff*nbast)
 !       fock_mtx(nbast*j+1-nbast:j*nbast,nbast*i+1-nbast:i*nbast)=&
 !       fock_mtx(nbast*i+1-nbast:i*nbast,nbast*j+1-nbast:j*nbast)
 !     ENDIF
 !   ! DO k=refindex*nbast-nbast+1,refindex*nbast 
 !   !   i=i+1
 !   !   j=0
 !   !  DO m=index1*nbast-nbast+1,index1*nbast 
 !   !    fock_mtx(k,m)=fock_mtx(k,m)+fock_tmp(i,j)!*coeff(k,m)
 !   !    fock_mtx(m,k)=fock_mtx(k,m)
 !   !  ENDDO
 !   ! ENDDO
   ENDDO
  ENDDO
  write(lupri,*) 'FOCK MATRIX'
  call write_matrix(fock_mtx,sizef,sizef)

END SUBROUTINE pbc_complete_Fock_mtx

END MODULE pbc_interactions

#endif
