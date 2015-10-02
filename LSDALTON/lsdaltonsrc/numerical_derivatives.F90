!> @file 
!> Contains the trilevel and atoms starting guess in addition to the construction of the GCbasis se PCCP 2009, 11, 5805-5813 
!> Trilevel and Atoms module
!> \author Katinka Dankel and Simen Reine
!> \date 2012-02-03
module Numerical_Hessian
use precision
use matrix_module, only: matrix
use configurationType, only: configitem
use typedeftype, only: lsitem
use molecule_typetype, only: moleculeinfo
use molecule_type, only: free_moleculeinfo
use memory_handling, only: mem_alloc, mem_dealloc
use typedef, only: copy_molecule
use matrix_operations, only :mat_init,mat_free
use energy_and_deriv, only: get_energy, get_gradient
use ls_util, only: ls_print_gradient
private
public :: get_numerical_hessian
contains

!> \brief 
!> \author Katinka Dankel and Simen Reine 
!> \date 2012-03-03
!> \param lupri logical unit number for the output file
subroutine get_numerical_hessian(lupri,luerr,ls,nbast,config,doNumHess,doNumGrad,doNumGradHess)
implicit none
integer, intent(in)     :: lupri,luerr, nbast
type(lsitem) :: ls
type(configItem)    :: config
!
integer :: i, j, k, m, DegFree,nAtoms
Type(Matrix) :: H1,S,D(1),F(1),C
real(realk) :: E(1), x, h,EigValues(3*ls%INPUT%MOLECULE%nAtoms)
real(realk),pointer :: gradient_min(:,:),gradient_plus(:,:),symmetric_hessian(:,:),Emin(:,:),Eplus(:,:),gradient(:,:), & 
     & analytical_gradient(:,:),analytical_gradient_min(:,:),analytical_gradient_plus(:,:), &
     & numerical_gradient(:,:),numerical_gradient_min(:,:),numerical_gradient_plus(:,:), &
     & numerical_hessian(:,:),semi_analytical_hessian(:,:)
Type(moleculeinfo),pointer :: unperturbed_molecule
character(len=1) :: V,L
real(realk),allocatable :: WORK(:)
integer :: LWORK,IERR
logical :: doNumHess,doNumGrad,doNumGradHess,debugAnaGrad
real(realk) :: Eerr

debugAnaGrad = .FALSE.

CALL mat_init(H1,nbast,nbast)
CALL mat_init(S,nbast,nbast)
CALL mat_init(D(1),nbast,nbast)
CALL mat_init(F(1),nbast,nbast)
CALL mat_init(C,nbast,nbast)

nullify(unperturbed_molecule)
allocate(unperturbed_molecule)
LWORK = 3*ls%INPUT%MOLECULE%nAtoms*3-1
allocate(WORK(LWORK))
nAtoms=ls%INPUT%MOLECULE%nAtoms
h = 1.0E-5_realk !1.0E-7_realk
call copy_molecule(ls%INPUT%MOLECULE,unperturbed_molecule,lupri)


call mem_alloc(numerical_hessian,3*ls%INPUT%MOLECULE%nAtoms,3*ls%INPUT%MOLECULE%nAtoms)
call mem_alloc(semi_analytical_hessian,3*ls%INPUT%MOLECULE%nAtoms,3*ls%INPUT%MOLECULE%nAtoms)
call mem_alloc(symmetric_hessian,3*ls%INPUT%MOLECULE%nAtoms,3*ls%INPUT%MOLECULE%nAtoms)
call mem_alloc(gradient,1,3*ls%INPUT%MOLECULE%nAtoms)
call mem_alloc(analytical_gradient,3,nAtoms)
call mem_alloc(analytical_gradient_plus,3,nAtoms)
call mem_alloc(analytical_gradient_min,3,nAtoms)
call mem_alloc(numerical_gradient_plus,ls%INPUT%MOLECULE%nAtoms,3)
call mem_alloc(numerical_gradient_min,ls%INPUT%MOLECULE%nAtoms,3)
call mem_alloc(numerical_gradient,ls%INPUT%MOLECULE%nAtoms,3)

if(doNumGrad) then
   IF(debugAnaGrad) THEN
       call get_gradient(E(1),Eerr,lupri,ls%INPUT%MOLECULE%nAtoms,S,F(1),D(1),ls,config,C,analytical_gradient)
   ENDIF
   call get_numerical_gradient(h,E(1),lupri,luerr,nbast,ls,H1,S,F(1),D(1),C,config,unperturbed_molecule,numerical_gradient)
   call LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,TRANSPOSE(numerical_gradient),ls%SETTING%MOLECULE(1)%p%nAtoms,'NUMGR')

   IF(debugAnaGrad) THEN
        call LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,analytical_gradient,ls%INPUT%MOLECULE%nAtoms,'Ana molgrad')
        call LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,TRANSPOSE(numerical_gradient),ls%INPUT%MOLECULE%nAtoms,'Num molgrad')
        DO i = 1,ls%INPUT%MOLECULE%nAtoms
             analytical_gradient(:,i) = numerical_gradient(i,:) - analytical_gradient(:,i)
        ENDDO
        write (*,*) "Difference: (Ana - Num) molgrad"
        write (lupri,*) "Difference: (Ana - Num) molgrad"
        call LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,analytical_gradient,ls%INPUT%MOLECULE%nAtoms,'Ana-Num molgrad')
   ENDIF
endif

!Checking whether the molecule is linear to determine vibrational degrees of freedom.
DegFree = 6
k=1
!do j = 1,3
!   do i=1,ls%INPUT%MOLECULE%nAtoms-1
!      if (ls%INPUT%MOLECULE%ATOM(i)%CENTER(j) == ls%INPUT%MOLECULE%ATOM(i+1)%CENTER(j)) then
!         k=k+1
!      endif
!   enddo
   
!   if (k == (ls%INPUT%MOLECULE%nAtoms)) then
!       DegFree = 5
!   endif
!enddo


!Analytical gradient used to calculate Hessian
if (doNumHess) then
   call get_hessian_from_analytical_gradient(ls,H1,S,D,F,C,E,ls%INPUT%MOLECULE%nAtoms,nbast,lupri,luerr,config, &   
        & unperturbed_molecule,analytical_gradient_min,analytical_gradient_plus,semi_analytical_hessian,h,symmetric_hessian)
   CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
   call print_hessian_and_vibrations(E(1),ls,symmetric_hessian,EigValues,DegFree,nAtoms,WORK,LWORK,IERR,lupri,luerr)
endif

!Numerical gradient used to calculated Hessian
if (doNumGradHess) then
   call get_hessian_from_numerical_gradient(ls,H1,S,D,F,C,E,nAtoms,nbast,lupri,luerr,&
     & config,unperturbed_molecule,numerical_gradient_min,numerical_gradient_plus,numerical_hessian,symmetric_hessian,h,&
     & analytical_gradient_min,analytical_gradient_plus)
   CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
   call print_hessian_and_vibrations(E(1),ls,symmetric_hessian,EigValues,DegFree,nAtoms,WORK,LWORK,IERR,lupri,luerr)
endif

call free_Moleculeinfo(unperturbed_molecule)
CALL mat_free(H1)
CALL mat_free(S)
CALL mat_free(D(1))
CALL mat_free(F(1))
CALL mat_free(C)
deallocate(unperturbed_molecule)
!nullify(unperturbed_molecule)
call mem_dealloc(numerical_hessian)
call mem_dealloc(semi_analytical_hessian)
call mem_dealloc(symmetric_hessian)
call mem_dealloc(gradient)
call mem_dealloc(analytical_gradient)
call mem_dealloc(analytical_gradient_min)
call mem_dealloc(analytical_gradient_plus)
call mem_dealloc(numerical_gradient_min)
call mem_dealloc(numerical_gradient_plus)
call mem_dealloc(numerical_gradient)

deallocate(WORK)

end subroutine get_numerical_hessian





!*******************************************************************************************************************!
!*******************************************************************************************************************!


!> \brief                                                                                                                                                                                                                                        
!> \author Katinka Dankel and Simen Reine
subroutine get_numerical_gradient(h,E,lupri,luerr,nbast,ls,H1,S,F,D,C,config,unperturbed_molecule,numerical_gradient)

implicit none
integer, intent(in)     :: lupri,nbast,luerr
type(lsitem)        :: ls
Type(Matrix)        :: H1,S,D(1),F(1),C
type(configItem)    :: config
real(realk), pointer :: numerical_gradient(:,:)
real(realk) :: E(1),Emin,Eplus,h
integer :: i, j
Type(moleculeinfo),pointer :: unperturbed_molecule
real(realk) :: Eerr

do i=1,ls%INPUT%MOLECULE%nAtoms
   do j=1, 3
   write(lupri,*) 'atom: ', i, 'out of ', ls%INPUT%MOLECULE%nAtoms
   write(lupri,*) 'coordinate: ', j, 'out of 3'
   write(*,*) 'atom: ', i, 'out of ', ls%INPUT%MOLECULE%nAtoms
   write(*,*) 'coordinate: ', j, 'out of 3'
      ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)-h 
      CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,ls%INPUT%MOLECULE%nAtoms,lupri,luerr)
      Emin=E(1)
      !call copy_molecule(unperturbed_molecule,ls%INPUT%MOLECULE,lupri)
      !call copy_molecule(unperturbed_molecule,Config%Molecule,lupri)
      
      ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)+(2*h)
      CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,ls%INPUT%MOLECULE%nAtoms,lupri,luerr)
      Eplus=E(1)
      ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)-h
      
      !call copy_molecule(unperturbed_molecule,ls%INPUT%MOLECULE,lupri)
      !call copy_molecule(unperturbed_molecule,Config%Molecule,lupri)
      numerical_gradient(i,j)=(Eplus-Emin)/(2*h)
        
    enddo
enddo



end subroutine get_numerical_gradient







!*******************************************************************************************************************!
!*******************************************************************************************************************!


!> \brief                                                                                                                                                                                                                                        
!> \author Katinka Dankel and Simen Reine
subroutine get_hessian_from_analytical_gradient(ls,H1,S,D,F,C,E,nAtoms,nbast,lupri,luerr,config,&
     & unperturbed_molecule,analytical_gradient_min,analytical_gradient_plus,semi_analytical_hessian,h,symmetric_hessian)
implicit none
integer, intent(in)     :: lupri,nbast,nAtoms,luerr
type(lsitem)        :: ls
Type(Matrix)        :: H1,S,D(1),F(1),C
type(configItem)    :: config
real(realk), pointer :: semi_analytical_hessian(:,:),analytical_gradient_min(:,:),analytical_gradient_plus(:,:),&
     &symmetric_hessian(:,:)
real(realk) :: E(1),Eplus,Emin,h
integer :: i,j,k,l,m
Type(moleculeinfo),pointer :: unperturbed_molecule
real(realk) :: Eerr

do i=1,nAtoms
   do j=1,3 

            !subtracting an increment h to coordinate j on atom i, calculating the gradient matix for this geometry 
            ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)-h 
            
            CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
            CALL get_gradient(E(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,analytical_gradient_min)
            call copy_molecule(ls%INPUT%MOLECULE,Config%Molecule,lupri)
            !Removing the perturbation
            call copy_molecule(unperturbed_molecule,ls%INPUT%MOLECULE,lupri)
            call copy_molecule(unperturbed_molecule,Config%Molecule,lupri)
           
            !adding an increment h to coordinate j on atom i, calculating the gradient for this geometry 
            ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)+h
            call copy_molecule(ls%INPUT%MOLECULE,Config%Molecule,lupri)
            CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
            CALL get_gradient(E(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,analytical_gradient_plus)
            !Removing the perturbation
            call copy_molecule(unperturbed_molecule,ls%INPUT%MOLECULE,lupri)
            call copy_molecule(unperturbed_molecule,Config%Molecule,lupri)
            
            
            !print*, 'ANALYTISKE GRADIENTEN'
            !do m=1, ls%INPUT%MOLECULE%nAtoms
            !   print*,     analytical_gradient_min(1,m),&
            !        &  analytical_gradient_min(2,m),&
            !        &  analytical_gradient_min(3,m)
            !enddo
            
            !do m=1, ls%INPUT%MOLECULE%nAtoms
            !   print*,     analytical_gradient_plus(1,m),&
             !       &  analytical_gradient_plus(2,m),&
              !      &  analytical_gradient_plus(3,m)
            !enddo
            
            !using analytical gradient to calculate Hessian (row by row)
            do k=1,nAtoms
               do  l=1,3
                  !hessian(alle radene,alle kolonner)
                  semi_analytical_hessian(3*(i-1)+j,3*(k-1)+l)=((analytical_gradient_plus(l,k)-analytical_gradient_min(l,k))/(2*h))
              enddo
            enddo
     enddo
enddo




!do i=1,3*ls%INPUT%MOLECULE%nAtoms
 !  do j=1,3*ls%INPUT%MOLECULE%nAtoms
  !    print*, semi_analytical_hessian(i,j)
 !  enddo
!enddo


do i=1,3*ls%INPUT%MOLECULE%nAtoms
   do j=1,i
      symmetric_hessian(i,j)= (semi_analytical_hessian(i,j)+semi_analytical_hessian(j,i))/2
   enddo
enddo


do i=1,3*ls%INPUT%MOLECULE%nAtoms
   do j=i,3*ls%INPUT%MOLECULE%nAtoms
      symmetric_hessian(i,j)= symmetric_hessian(j,i)
   enddo
enddo

end subroutine get_hessian_from_analytical_gradient

!*******************************************************************************************************************!
!*******************************************************************************************************************!



!> \brief                                                                                                                                                                                                                                        
!> \author Katinka Dankel and Simen Reine
subroutine get_hessian_from_numerical_gradient(ls,H1,S,D,F,C,E,nAtoms,nbast,lupri,luerr,config,&
     & unperturbed_molecule,numerical_gradient_min,numerical_gradient_plus,numerical_hessian,symmetric_hessian,h,&
     & analytical_gradient_min,analytical_gradient_plus)
implicit none
integer, intent(in)     :: lupri,nbast,nAtoms,luerr
type(lsitem)        :: ls
Type(Matrix)        :: H1,S,D(1),F(1),C
type(configItem)    :: config
real(realk), pointer :: numerical_hessian(:,:),numerical_gradient_min(:,:), numerical_gradient_plus(:,:),symmetric_hessian(:,:), &
     & analytical_gradient_min(:,:),analytical_gradient_plus(:,:)
real(realk) :: E(1),Eplus,Emin,h
integer :: i,j,k,l,m
Type(moleculeinfo),pointer :: unperturbed_molecule
real(realk) :: Eerr


do i=1,nAtoms
   do j=1,3 

            !subtracting an increment h to coordinate j on atom i, calculating the gradient matix for this geometry 
            ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)-h 
            call get_numerical_gradient(h,E,lupri,luerr,nbast,ls,H1,S,F,D,C,config,&
                       & unperturbed_molecule,numerical_gradient_min)
            CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
            CALL get_gradient(E(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,analytical_gradient_plus)
            
            !Removing the perturbation
            call copy_molecule(unperturbed_molecule,ls%INPUT%MOLECULE,lupri)
            call copy_molecule(unperturbed_molecule,Config%Molecule,lupri)
           
            !adding an increment h to coordinate j on atom i, calculating the gradient for this geometry 
            ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)=ls%INPUT%MOLECULE%ATOM(i)%CENTER(j)+h
            call get_numerical_gradient(h,E,lupri,luerr,nbast,ls,H1,S,F,D,C,config,&
                       &unperturbed_molecule,numerical_gradient_plus)
            CALL get_energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
            
            CALL get_gradient(E(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,analytical_gradient_plus)
            !Removing the perturbation
            call copy_molecule(unperturbed_molecule,ls%INPUT%MOLECULE,lupri)
            call copy_molecule(unperturbed_molecule,Config%Molecule,lupri)
            
            !write(lupri,*)
            !write(lupri,*) 'NUMERISKE GRADIENTEN', i, j
            !do m=1, ls%INPUT%MOLECULE%nAtoms
            !    write(lupri,*)     numerical_gradient_min(1,m),&
            !        &  numerical_gradient_min(2,m),&
            !        &  numerical_gradient_min(3,m)
            !enddo
           
            !write(lupri,*)
            !write(lupri,*) 'DEN ANALYTISKE GRADIENTEN', i,j
           ! do m=1, ls%INPUT%MOLECULE%nAtoms
           !     write(lupri,*)     analytical_gradient_min(1,m),&
           !         &  analytical_gradient_min(2,m),&
           !         &  analytical_gradient_min(3,m)
           ! enddo
           !  write(lupri,*)
            
           !  write(lupri,*) 'DEN NUMERISKE GRADIENTEN',i,j
           ! do m=1, ls%INPUT%MOLECULE%nAtoms
           !    write(lupri,*)     numerical_gradient_plus(1,m),&
           !         &  numerical_gradient_plus(2,m),&
           !         &  numerical_gradient_plus(3,m)
           ! enddo
           !  write(lupri,*)
            
           ! write(lupri,*) 'ANALYTISKE GRADIENTEN', i,j
           ! do m=1, ls%INPUT%MOLECULE%nAtoms
           !    write(lupri,*)     analytical_gradient_plus(1,m),&
           !         &  analytical_gradient_plus(2,m),&
           !         &  analytical_gradient_plus(3,m)
           ! enddo
           ! write(lupri,*)

            !using analytical gradient to calculate Hessian (row by row)
            do k=1,nAtoms
               do  l=1,3
                  !hessian(alle radene,alle kolonner)
                  numerical_hessian(3*(i-1)+j,3*(k-1)+l)=((numerical_gradient_plus(l,k)-numerical_gradient_min(l,k))/(2*h))
              enddo
            enddo
     enddo
enddo



do i=1,3*ls%INPUT%MOLECULE%nAtoms
   do j=1,i
      symmetric_hessian(i,j)= (numerical_hessian(i,j)+numerical_hessian(j,i))/2
   enddo
enddo

do i=1,3*ls%INPUT%MOLECULE%nAtoms
   do j=i,3*ls%INPUT%MOLECULE%nAtoms
      symmetric_hessian(i,j)= symmetric_hessian(j,i)
   enddo
enddo


!write(lupri,*)


!call ls_output(numerical_hessian,1,3*ls%INPUT%MOLECULE%nAtoms,1,3*ls%INPUT%MOLECULE%nAtoms,&
!     & 3*ls%INPUT%MOLECULE%nAtoms,3*ls%INPUT%MOLECULE%nAtoms,1,lupri)

end subroutine get_hessian_from_numerical_gradient


!*******************************************************************************************************************!
!*******************************************************************************************************************!



!> \brief                                                                                                                                                                                                                                        
!> \author Katinka Dankel and Simen Reine
subroutine print_hessian_and_vibrations(E,ls,symmetric_hessian,EigValues,DegFree,nAtoms,WORK,LWORK,IERR,lupri,luerr)
implicit none
integer, intent(in)     :: lupri,luerr,nAtoms
type(lsitem)        :: ls
type(configItem)    :: config
real(realk), pointer :: symmetric_hessian(:,:)
real(realk) :: E,Eplus,Emin,h,EigValues(:),Vibrational_Frequencies(3*nAtoms)
real(realk) :: ZPVE,Hvib,Hvib_Hartree,Svib,Svib_Hartree
real(realk) :: hc_kBT,hc,kB,T,Stot,Gibbs,CmEh,kBT_Hartree
integer :: i,j,k,LWORK,IERR,DegFree
character(len=1) :: V,L
real(realk),allocatable :: WORK(:)
IERR=0





!Writing the Hessian to file
write(lupri,*) 
write(lupri,*) 
write(lupri,*) 
write(lupri,*) 
write(lupri,'(8X,A)') '*************************************************************'
write(lupri,'(8X,A)') '*            MOLECULAR HESSIAN RESULTS                      *'
write(lupri,'(8X,A)') '*************************************************************'
write(lupri,*) 



write(lupri,*)
write(lupri,*)
write(lupri,*) '                            Molecular Hessian (au)'
write(lupri,*) '                            ----------------------'
write(lupri,*) 'The rows and colums represent the spacial coordiantes in the following order'
do i = 1,nAtoms
   write(lupri,*) ls%INPUT%MOLECULE%ATOM(i)%NAME, '(x,y,z) '
enddo
 
write(lupri,*)
call ls_output(symmetric_hessian,1,3*nAtoms,1,3*nAtoms,3*nAtoms,3*nAtoms,1,lupri)


!Mass Weighting the Hessian
 do i=1,3*nAtoms
    do  j=1,nAtoms
       do k=1,3
          symmetric_hessian(i,3*(j-1)+k)=symmetric_hessian(i,3*(j-1)+k)/(sqrt(ls%input%Molecule%Atom(j)%Mass))
          symmetric_hessian(3*(j-1)+k,i)=symmetric_hessian(3*(j-1)+k,i)/sqrt(ls%input%Molecule%Atom(j)%Mass)
       enddo
    enddo
 enddo

!Writing the Mass Weighted Hessian to file
write(lupri,*)
write(lupri,*)
write(lupri,*)
write(lupri,*) '                       Mass Weighted Molecular Hessian (au)'
write(lupri,*) '                       ------------------------------------'
write(lupri,*) 'The rows and colums represent the spacial coordiantes in the following order'
do i = 1,ls%INPUT%MOLECULE%nAtoms
   write(lupri,*) ls%INPUT%MOLECULE%ATOM(i)%NAME, '(x,y,z) '
enddo
write(lupri,*)
call ls_output(symmetric_hessian,1,3*nAtoms,1,3*nAtoms,3*nAtoms,3*nAtoms,1,lupri)


!Diagonalizing
call DSYEV('V','L',nAtoms*3,symmetric_hessian,nAtoms*3,EigValues,WORK,LWORK,IERR)
if (IERR .ne. 0) STOP 'Something went wrong in DSYEV in Solve_EigVal_RedSpace'
!deallocate(WORK)

       

!Printing the eigenvalues, the corresponding vibrational frequencies and the normal modes
write(lupri,*)
write(lupri,*) '                           Eigenvalues of Hessian (au)'
write(lupri,*) '                           --------------------------'
write(lupri,*) 
do i=1, 3*nAtoms
   write(lupri,*) EigValues(i)
enddo

write(lupri,*)
write(lupri,*) 
write(lupri,*) '                     Vibrational Frequencies (resiprocal cm)'
write(lupri,*) '                     ---------------------------------------'
write(lupri,*) 
do i=DegFree+1, 3*nAtoms
   Vibrational_frequencies(i) = EigValues(i)*26424608
   Vibrational_frequencies(i) = sqrt(Vibrational_frequencies(i))
   write(lupri,*) Vibrational_frequencies(i), 'cm^-1'
enddo


write(lupri,*)
write(lupri,*)
write(lupri,*) '                                Normal Modes'
write(lupri,*) '                               -------------'
write(lupri,*)

call ls_output(symmetric_hessian,1,3*nAtoms,DegFree+1,3*nAtoms,&
     & 3*nAtoms,3*nAtoms,1,lupri)
write(lupri,*)
write(lupri,*) 
write(lupri,*)
write(lupri,*) 



write(lupri,*)
write(lupri,*)
write(lupri,*) '                           ZPVE, Hvib, Svib'
write(lupri,*) '                          -----------------'
write(lupri,*)
write(lupri,*)

CmEh = 4.55633525E-6_realk !conversion factor from reciprocal cm to Hartree


ZPVE=0_realk
do i= DegFree+1, 3*nAtoms
   ZPVE = ZPVE + 0.5*Vibrational_frequencies(i)
enddo
write(lupri,*) ZPVE, ' cm-1'

ZPVE = ZPVE*CmEh !Converting ZPVE in reciprocal cm to Hartree
write(lupri,*) ZPVE, ' Hartree'
write(lupri,*)
write(lupri,*)
write(lupri,*)



Hvib=0_realk
hc_kBT = 4.825680E-3_realk !(h*c)/(kB*T) in reciprocal cm
hc=1.9864455684E-23_realk ! (h*c) in J*s*cm
do i= DegFree+1, 3*nAtoms
   Hvib = Hvib + Vibrational_frequencies(i)/(exp(hc_kBT* Vibrational_frequencies(i))-1)
   write(lupri,*) Hvib
enddo

write(lupri,*)
write(lupri,*) 'Hvib: ', Hvib, ' cm-1'!ok
Hvib_Hartree= Hvib*CmEh
write(lupri,*) 'Hvib: ', Hvib, ' Hartree'!ok
write(lupri,*)
write(lupri,*)

kBT_Hartree = 9.44185E-4_realk 
if (DegFree == 6) then
   write(lupri,*) 'hterm' , Hvib + 4*kBT_Hartree
   write(lupri,*) 'HTERM' , ZPVE + Hvib + 4*kBT_Hartree
   H = E + ZPVE + Hvib_Hartree + 4*kBT_Hartree
else 
   H = E + ZPVE + Hvib_Hartree + (7/2)*kBT_Hartree
endif

write(lupri,*) 'Htot: ', H, ' Hartree/particle'





!Svib = 0_realk
!kB = 1.38E-23_realk !Boltzmann
!write (lupri,*) kB
!do i= DegFree+1, 3*nAtoms
!   Svib = Svib + kB*((hc_kBT*Vibrational_frequencies(i)/(exp(hc_kBT*Vibrational_frequencies(i))-1)) - &
!        & log(1-exp(-hc_kBT*Vibrational_frequencies(i))))
!   write(lupri,*) Svib
!enddo


!write(lupri,*) 'Svib: ', Svib, ' J/K'
!Svib_Hartree= Svib/(4.35974394E-18)
!write(lupri,*) 'Svib: ', Svib_Hartree, ' Hartree'
!write(lupri,*) 'E: ', E, ' Hartree'



!Calculating the molecular entropy at room temperature
!if (DegFree == 6) then
!   Stot = Svib_Hartree + 3*kB/(4.35974394E-18)
!else 
!   Stot = Svib_Hartree + (5/2)*kB/(4.35974394E-18)
!endif

!write(lupri,*) 'S: ', Stot, ' Hartree'


!Calculating the molecular gibbs free energy at room temperature
!Gibbs = H + T*Stot
!write(lupri,*) 'Gtot: ', Gibbs
write(lupri,*)
!write(lupri,*) 'Gaussian: Energy -76.02187755 -76.0213409549'
write(lupri,*) 'Dalton: Energy  ', E
write(lupri,*)
!write(lupri,*) 'Gaussian: Zero-point correction 0.021538 0.021542'
write(lupri,*) 'Dalton: Zero-point correction  ', ZPVE
write(lupri,*)
!write(lupri,*) 'Gaussian: Thermal correction to enthalpy 0.024372 0.025320'
write(lupri,*) 'Dalton: Thermal correction to enthalpy ', ZPVE + Hvib_Hartree + 4*kBT_Hartree
write(lupri,*)
!write(lupri,*) 'Gaussian: Thermal correction to Gibbs free energy 0.003865 0.003868'
!write(lupri,*) 'Dalton: Thermal correction to Gibbs free energy  ', Stot*298.15+kBT_Hartree 
write(lupri,*)
write(lupri,*)
write(lupri,*)

end subroutine print_hessian_and_vibrations


end module Numerical_Hessian

