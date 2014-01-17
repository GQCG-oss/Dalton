#ifdef MOD_UNRELEASED
MODULE multipole_pbc
use precision
use harmonics_pbc
USE matrix_module
USE pbc_MSC
USE TYPEDEF
USE lattice_type
USE lattice_vectors
USE pbc_matrix_operations

IMPLICIT NONE
!integer, parameter :: COMPLEXK = KIND( (1D0,1D0) )
TYPE, PUBLIC :: orbital_order_t
  COMPLEX(complexk), pointer :: celms(:)
  REAL(realk),  pointer   :: elms(:)
END TYPE
CONTAINS


SUBROUTINE multipole_pbc_initial(orbital_ctrl,totdim,cor)
IMPLICIT NONE
TYPE(orbital_order_t),INTENT(INOUT),target ::orbital_ctrl
INTEGER,INTENT(IN) :: totdim
CHARACTER(len=1),INTENT(IN) :: cor

if(cor .eq. 'c' .or. cor .eq. 'C') THEN
  allocate(orbital_ctrl%celms(totdim))
elseif(cor .eq. 'r' .or. cor .eq. 'R') THEN
  allocate(orbital_ctrl%elms(totdim))
else
  write(*,*) 'Third input in routine multipole_pbc_initial must be,&
 & r, R, c or C'
  STOP
endif

END SUBROUTINE multipole_pbc_initial

!The below subroutine computes the local expansion at the the points relative
!to origin, It is not optimal yet, I have to fix it as well. 
SUBROUTINE pbc_local_expan_mult(setting,nbast,ll,latt_cell,&
                            &refcell,numvecs,maxmultmom)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: numvecs,maxmultmom,nbast
  TYPE(lvec_list_t),intent(in) :: ll
  TYPE(LSSETTING)   :: SETTING 
  TYPE(moleculeinfo),intent(inout) :: refcell
  TYPE(moleculeinfo),intent(in) :: latt_cell(numvecs)

  INTEGER, DIMENSION(3) :: fdim
  Integer :: refindex,index,il1,il2,il3,nSpherMom,l,m
  TYPE(matrix),pointer :: spherMom(:)
  COMPLEX(complexk) :: irreg_solid,irreg_solidt
  CHARACTER(len=20) ::filename
  CHARACTER(len=10) :: numtostring1,numtostring2,numtostring3
  INTEGER :: wunit

  wunit=9920
  nSpherMom=0

  DO l=0,maxmultmom
      nSpherMom=nSpherMom+2*l+1
  ENDDO
  
  write(*,*) 'PBC_LOCALMOMENT, maxMultmom,nspherMom', maxMultmom,nspherMom
  write(*,*) 'DEBUG INSIDE PBC_LOCALMOMENT'
  
  call mem_alloc(spherMom,nSpherMom)  
  DO l=1,nSpherMom
     call Mat_init(spherMom(l),nbast,nbast)
     call Mat_zero(spherMom(l))
  ENDDO
  call find_latt_index(refindex,0,0,0,fdim,ll,ll%max_layer)
     filename='Loclm'
  irreg_solidt=0.0
  irreg_solid =0.0
  !call LSOPEN(wunit,filename,'unknown','formatted')
  !DO l=0,2*maxmultmom
  !DO m=-l,l
  DO index=1,numvecs
     
     !phase1 = kvec(1)*ll%lvec(index)%std_coord(1)
     !phase2 = kvec(2)*ll%lvec(index)%std_coord(2)
     !phase3 = kvec(3)*ll%lvec(index)%std_coord(3)
     !phase=CMPLX(0.,phase1+phase2+phase3)
     !phase=exp(phase)

     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index),2)

     call find_latt_vectors(index,il1,il2,il3,fdim,ll)
     if(abs(il1) .gt. ll%nneighbour+1)CYCLE
     if(abs(il2) .gt. ll%nneighbour+1)CYCLE
     if(abs(il3) .gt. ll%nneighbour+1)CYCLE
     if((abs(il1) .ne. ll%nneighbour+1 .and. abs(il2) .ne. ll%nneighbour+1) &
     & .and. abs(il3) .ne. ll%nneighbour+1) CYCLE
    ! if(abs(il1) .le. abs(ll%nneighbour) .and. abs(il1) &
    !       &.gt. abs(ll%nneighbour+1)) CYCLE
    ! if(abs(il2) .le. abs(ll%nneighbour) .and. abs(il2) &
    !       &.gt. abs(ll%nneighbour+1)) CYCLE
    ! if(abs(il3) .le. abs(ll%nneighbour) .and. abs(il3) &
    !       &.gt. abs(ll%nneighbour+1)) CYCLE
     !input%Basis%regular%atomtype(:)%shell(:)%segment(1)%exponents(:)

     setting%scheme%PS_SCREEN = .FALSE.
     setting%scheme%CS_SCREEN = .FALSE.
     
     write(numtostring1,'(I5)')  il1
     write(numtostring2,'(I5)')  il2
     write(numtostring3,'(I5)')  il3

     numtostring1=adjustl(numtostring1)
     numtostring2=adjustl(numtostring2)
     numtostring3=adjustl(numtostring3)

     filename='Loclm'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'

     call LSOPEN(wunit,filename,'unknown','formatted')
   !  call pbc_T_tensor(ll%lvec(index)%std_coord,l,m,irreg_solid)
   !  irreg_solid=irreg_solid*factorial(l+abs(m))
   !  irreg_solidt=irreg_solidt+irreg_solid
      write(*,*) filename
      DO l=0,2*maxmultmom
        DO m=-l,l
           !write(*,*) ll%lvec(index)%std_coord, il1,il2,il3,index
           !STOP
           
           call pbc_T_tensor(ll%lvec(index)%std_coord,l,m,irreg_solid)
           irreg_solid=irreg_solid!factorial(l+abs(m))
           write(wunit,*) irreg_solid!,l,m,factorial(l+abs(m)),maxmultmom
        ENDDO
     ENDDO
     call LSCLOSE(wunit,'KEEP')

  

 ! ENDDO
 ! write(wunit,*) irreg_solidt ,l,m
 ! irreg_solidt=0.0
 ! ENDDO
  write(*,*) 'DEBUG LOCALMOMENT after first loop',l
  ENDDO


  DO l=1,nSpherMom
     call Mat_free(spherMom(l))
  ENDDO
  call mem_dealloc(spherMom)
     
  !call LSCLOSE(wunit,'KEEP')



END SUBROUTINE pbc_local_expan_mult

!This routine computes multipole moments for nearest neighbours, it will
!probably not be used in the program, but it is maybe nice to keep it for later
!convenience.
SUBROUTINE pbc_multipol_nn(setting,nbast,ll,latt_cell,&
                            &refcell,numvecs,maxmultmom)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: numvecs,maxmultmom,nbast
  TYPE(lvec_list_t),intent(in) :: ll
  TYPE(LSSETTING)   :: SETTING 
  TYPE(moleculeinfo),intent(inout) :: refcell
  TYPE(moleculeinfo),intent(in) :: latt_cell(numvecs)

  INTEGER, DIMENSION(3) :: fdim
  REAL(realk),DIMENSION(3) :: lvec
  Integer :: refindex,index,il1,il2,il3,nSpherMom,l,m
  TYPE(matrix),allocatable :: spherMom(:)
  COMPLEX(complexk) :: cart_solid

  nSpherMom=0
  DO l=0,maxmultmom
      nSpherMom=nSpherMom+2*l+1
  ENDDO
  
  write(*,*) 'PBC_MULTIPOL_NN, maxMultmom,nspherMom', maxMultmom,nspherMom
  write(*,*) 'DEBUG INSIDE PBC_MULTIPOL_NN'
  
  Allocate(spherMom(nSpherMom))  
  DO l=1,nSpherMom
     call Mat_init(spherMom(l),nbast,nbast)
     call Mat_zero(spherMom(l))
  ENDDO
  call find_latt_index(refindex,0,0,0,fdim,ll,ll%max_layer)
  DO index=1,numvecs
     
     !phase1 = kvec(1)*ll%lvec(index)%std_coord(1)
     !phase2 = kvec(2)*ll%lvec(index)%std_coord(2)
     !phase3 = kvec(3)*ll%lvec(index)%std_coord(3)
     !phase=CMPLX(0.,phase1+phase2+phase3)
     !phase=exp(phase)

     call TYPEDEF_setmolecules(setting,refcell,1,latt_cell(index),2)

     call find_latt_vectors(index,il1,il2,il3,fdim,ll)
     if(abs(il1) .gt. abs(ll%nneighbour) ) CYCLE
     if(abs(il2) .gt. abs(ll%nneighbour)) CYCLE
     if(abs(il3) .gt. abs(ll%nneighbour)) CYCLE
     !input%Basis%regular%atomtype(:)%shell(:)%segment(1)%exponents(:)

     setting%scheme%PS_SCREEN = .FALSE.
     setting%scheme%CS_SCREEN = .FALSE.
     lvec=ll%lvec(index)%std_coord
     DO l=0,maxmultmom
        DO m=-l,l
           call cartesian_solid_harmonic&
           &(lvec(1),lvec(2),lvec(3),l,m,cart_solid,.true.)
        ENDDO
     ENDDO


  ENDDO

  DO l=1,nSpherMom
     call Mat_free(spherMom(l))
  ENDDO

END SUBROUTINE pbc_multipol_nn


!!Computes the irregular solid, for the points in the vector
!lvec, I have to verify the coefficients.
SUBROUTINE pbc_T_tensor(lvec,l,m,irreg_solid)
  IMPLICIT NONE
  REAL(realk),intent(in),DIMENSION(3) :: lvec
  INTEGER,INTENT(IN) :: l,m
  COMPLEX(complexk),INTENT(OUT) :: irreg_solid
  !LOCAL VARIABLES
!  real(realk) :: pi=3.14159265
  COMPLEX(complexk) :: cart_solid
  INTEGER ::fildum
  !ll%lvec(j)%std_coord(123)
  !binfac=binomial(2*l+2*lambda+1,2*lambda)**0.5
  !call clebsch_g(lambda,mu,l,m,clebschfac)
  call cartesian_irreg_solid_harmonic(lvec(1),lvec(2),lvec(3),l,m,cart_solid,.false.)
  fildum=9930
  OPEN(fildum,FILE='solidcart')
  !write(fildum,*) cart_solid
  irreg_solid=cart_solid!/(lvec(1)**2+lvec(2)**2+lvec(3)**2)**(l+0.5)
  
  !irreg_solid=sqrt(factorial(l+m)/factorial(l-m))*sqrt(4.*pi/(2.*l+1.))*irreg_solid
  write(fildum,*) irreg_solid,(lvec(1)**2+lvec(2)**2+lvec(3)**2)**(l+0.5),l,m

END SUBROUTINE pbc_T_tensor


!This routine will do the iteration technique discussed in article
!J.Chem.Phys. Vol. 121, No.7 15 August 2004 by K.N. Kudin and G.E. Scuseria.
SUBROUTINE pbc_iterated_localmom(ll,niter,maxmultmom,filename1,filename2,&
                                  &numvecs,multpmom,nbast,siterlmt)
  IMPLICIT NONE
  TYPE(lvec_list_t),intent(in) :: ll
  INTEGER, INTENT(IN) :: niter,maxmultmom,numvecs
  INTEGER,INTENT(IN) :: nbast
  TYPE(lattice_cell_info_t),intent(INOUT)  :: multpmom!(:)
  TYPE(orbital_order_t), pointer ::  Siterlmt(:)
  !TYPE(lattice_cell_info_t),pointer,intent(INOUT)  :: multpmom!(:)
  !LOCAL VARIABLES
  !TYPE(orbital_order_t),allocatable :: taylornn(:)
  COMPLEX(complexk),allocatable :: taylornn(:)
  TYPE(orbital_order_t),allocatable ::  Siterlm(:)
  INTEGER :: iteration,l,m,j,k,l1,l2,orbit,i,l3,m3
  INTEGER :: wunit,wunit2,numnn,totalmom
  CHARACTER(len=13),INTENT(IN) :: filename1
  CHARACTER(len=12),INTENT(IN) :: filename2

  wunit = 9920
  wunit2=9930

  totalmom= 4*maxmultmom**2+4*maxmultmom+1
  Allocate(taylornn(4*maxmultmom**2+4*maxmultmom+1))
  Allocate(Siterlm(4*maxmultmom**2+4*maxmultmom+1))
  Allocate(Siterlmt(4*maxmultmom**2+4*maxmultmom+1))

  !multpmom should be the size of amount of 'valid vectors' 
  !number of naerest neighbour cells

  DO i=1,4*maxmultmom**2+4*maxmultmom+1
     call multipole_pbc_initial(siterlm(i),nbast*nbast,'c')
     call multipole_pbc_initial(siterlmt(i),nbast*nbast,'c')
     !call multipole_pbc_initial(taylornn(i),nbast*nbast,'c')
  ENDDO

!  DO i=1,2*maxmultmom**2+2*maxmultmom+1
!     call multipole_pbc_initial(siterlmt(i),nbast*nbast,'c')
!  ENDDO

  call LSOPEN(wunit,filename1,'old','formatted')
  !call LSOPEN(wunit2,'kopi-200.dat','unknown','formatted')
  orbit=0
  DO l=0,2*maxmultmom
   DO m=-l,l
   orbit=orbit+1
    READ(wunit,*) taylornn(orbit)
  !write(wunit2,*) taylornn(orbit)
   ENDDO
  ENDDO
  call LSCLOSE(wunit,'KEEP')
  write(*,*) 'Done reading first file: ', filename1,orbit,maxmultmom**2+2*maxmultmom+1
  !STOP
  !DO l=1,orbit
  !write(wunit2,*) taylornn(orbit)
  !ENDDO
  !call LSCLOSE(wunit2,'KEEP')
  !STOP

  l1=0
  DO l=0,maxmultmom
   DO m=-l,l
      l1=l1+1
      Siterlm(l1)%celms(:)=taylornn(l1)
   ENDDO
 ENDDO

 !call READ_multipole_files(ll,numvecs,maxmultmom,l1,multpmom,nbast,numnn)
 
 !S0=taylor0
 ! call pbc_local_expan_mult(lupri,luerr,setting,molecule,nbast,ll,latt_cell,&
 !                           &refcell,numvecs,basinfo,maxmultmom,kvec)
 ! call pbc_multipol_nn(lupri,luerr,setting,molecule,nbast,ll,latt_cell,&
 !                           &refcell,numvecs,basinfo,maxmultmom,kvec)
  DO iteration=0,niter-1
  ! DO i=1,numnn
  l1=0
      DO l=0,maxmultmom
       DO m=-l,l
       l1=l1+1
       l2=0
        DO j=0,maxmultmom
         DO k=-j,j
         l2=l2+1
         l3=l+j
         m3=m+k
         if(abs(m3) .gt. 2*maxmultmom) CYCLE
         if(l3 .gt. 2*maxmultmom) CYCLE

         
         if(angmom_mapping_to(l3,m3) .gt. 961) write(*,*) 'bug',l3,m3,angmom_mapping_to(l3,m3)

       !  Siterlmt(l1)%celms(:)=taylornn(l1)+&
       !   &Siterlm(angmom_mapping_to(l3,m3))%celms(:)/(3**(l+j+1))&
       !   &*numnn*multpmom%getmultipole(l2)%elms(:)+Siterlmt(l1)%celms(:)
          !&Siterlm(angmom_mapping_to(l3,m3))%celms(:)/(3**(l+j+1))*multpmom(i)%getmultipole(l2)%elms(:)+Siterlmt(l1)%celms(:)

         ENDDO
        ENDDO
       ENDDO
      ENDDO
      ENDDO
      DO l=1,totalmom
      siterlm(l)%celms=Siterlmt(l)%celms
      siterlmt(l)%celms= CMPLX(0.0D0,0.0D0)
  !   ENDDO
  ENDDO

  DO l=1,maxmultmom**2+2*maxmultmom+1
   DO i=1,nbast*nbast
      DO k=1,nbast*nbast
         Siterlmt(l)%celms(i)=Siterlm(l)%celms(i)*&
         &multpmom%getmultipole(l)%elms(k)+Siterlmt(l)%celms(i)          
      ENDDO                                 !This is not completely wrong       
   ENDDO                                    !I have to make sure that       
  ENDDO                                     !k =all inactive states.          

END SUBROUTINE pbc_iterated_localmom

SUBROUTINE READ_multipole_files(il1,il2,il3,maxmultmom,sphermom,nbast,lupri)
  IMPLICIT NONE
  !TYPE(lvec_list_t), INTENT(IN) :: ll
  INTEGER, INTENT(IN) :: il1,il2,il3, maxmultmom,nbast,lupri
  TYPE(lattice_cell_info_t), INTENT(INOUT) :: sphermom!(:)
  !LOCAL VARIABLES
  INTEGER :: fileunit,totalmom !For file detection
  INTEGER :: i,j,stat
  INTEGER :: sphermoms,ii, debdum
  CHARACTER(len=10) :: numtostring1,numtostring2,numtostring3
  CHARACTER(LEN=20) :: filename
  REAL(realk) :: testread
  !For debugging
  !REAL(realk) :: debugsumT

  fileunit=9950

  
  totalmom=(maxmultmom+1)**2
  !numnn=0
!  DO dummy=1,numvecs

!     call find_latt_vectors(dummy,il1,il2,il3,fdim,ll)

     if(abs(il1) .gt. n_neighbour) call lsquit('ERROR in READ_multipole_files&
     & vector not in NF',lupri)
     if(abs(il2) .gt. n_neighbour) call lsquit('ERROR in READ_multipole_files&
     & vector not in NF',lupri)
     if(abs(il3) .gt. n_neighbour) call lsquit('ERROR in READ_multipole_files&
     & vector not in NF',lupri)
!     numnn=numnn+1

!  ENDDO
  !allocate(sphermom(numnn))

  !DO dummy=1,numnn
  !   allocate(sphermom(dummy)%getmultipole(totalmom))
  !   DO ii=1,totalmom
  !   call Mat_init(sphermom(dummy)%getmultipole(ii),nbast,nbast)
  !   call Mat_zero(sphermom(dummy)%getmultipole(ii))
  !   ENDDO
  !ENDDO
  call mem_alloc(sphermom%getmultipole,totalmom)
  DO ii=1,totalmom
     call Mat_init(sphermom%getmultipole(ii),nbast,nbast)
     call Mat_zero(sphermom%getmultipole(ii))
  ENDDO

  !numnn=0
  !DO dummy=1,numvecs

     !call find_latt_vectors(dummy,il1,il2,il3,fdim,ll)

!     if(abs(il1) .gt. ll%nneighbour) CYCLE
!     if(abs(il2) .gt. ll%nneighbour) CYCLE
!     if(abs(il3) .gt. ll%nneighbour) CYCLE
     !numnn=numnn+1
     write(numtostring1,'(I5)')  il1
     write(numtostring2,'(I5)')  il2
     write(numtostring3,'(I5)')  il3
     numtostring1=adjustl(numtostring1)
     numtostring2=adjustl(numtostring2)
     numtostring3=adjustl(numtostring3)
  !ENDDO
     !filename='Qcdlm.dat'!&
     filename='Qcdlm'//trim(numtostring1)//trim(numtostring2)//trim(numtostring3)//'.dat'
     !filename2=&
     !&'kopi'//trim(numtostring1)//trim(numtostring2)//trim(numtostring2)//'.dat'

     call LSOPEN(fileunit,filename, 'old','unformatted')
!     read(fileunit,IOSTAT=stat) testread
     !write(*,*) testread,stat
!    BACKSPACE(fileunit)
!     if(stat .ne. 0) then
       rewind(fileunit) 
!     endif
     !OPEN(UNIT=9998,FILE=filename2)
     debdum=0
     !write(lupri,*) 'Debug 1: ', fileunit
     DO sphermoms=1,totalmom
        
        !DO ii=1,nbast*nbast
        debdum=debdum+1

        call mat_read_from_disk(fileunit,sphermom%getmultipole(sphermoms),.true.)
        !READ(fileunit,*) sphermom(numnn)%getmultipole(sphermoms)%elms!(ii)  
        !do j=1,nbast
        !READ(fileunit) (sphermom%getmultipole(sphermoms)%elms(i+(j-1)*nbast),&
        !i=1,nbast)
        !enddo

        !debugsumT=0d0
        !do j=1,nbast*nbast
        !  debugsumT=debugsumT+sphermom%getmultipole(sphermoms)%elms(j)
        !enddo
        !if(debugsumT .gt. 0d0) then
        !  write(*,*) 'debugsumT > 0d0',il1
        !  write(*,*) debugsumT , sphermoms
        !  write(*,*) sphermom%getmultipole(sphermoms)%elms
        !  stop
        !endif
        !write(9998,*) sphermom(dummy)%getmultipole(sphermoms)%elms!(ii)
        !write(*,*) 'Debug 2: ', dummy, ii,nbast*nbast
        !ENDDO
     ENDDO
        !write(*,*) 'Debug 3: ', filename
     !stop
     call LSCLOSE(fileunit,'KEEP') 

  !ENDDO

!  write(*,*) 'Debug 4: ', numnn
!  stop
  !DO ii=1,totalmom
  !   call Mat_free(sphermom%getmultipole(ii))
  !enddo




END SUBROUTINE READ_multipole_files


INTEGER FUNCTION ANGMOM_MAPPING_TO(l,m)
IMPLICIT NONE
INTEGER, INTENT(IN) :: l,m
!Local variables
INTEGER :: ll,ml,l1
!MAKE A BETTER ALGORITHM IT TAKES TO MUCH TIME
l1=0

DO ll=0,l
 DO ml=-ll,ll
 l1=l1+1
 IF(ml .ne. m) CYCLE
 IF(ll .ne. l ) CYCLE
 if(l1 .gt. 961) write(*,*) l,m,ml,ll,l1
 angmom_mapping_to=l1
 ENDDO
ENDDO

!if(l .eq. 0 .and. m .eq. 0) angmom_mapping_to=1

END FUNCTION ANGMOM_MAPPING_TO
END MODULE multipole_pbc

#endif
