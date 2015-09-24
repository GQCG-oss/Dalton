module leastchange_module
use files
use precision
use decompMod
use TYPEDEF, only: getNbasis
use TYPEDEFTYPE, only: lssetting, lsitem
use IntegralInterfaceMOD
use matrix_util
use matrix_module
use matrix_operations
use LSparameters
use memory_handling
use TYPEDEF,only: count_ncore
contains

subroutine leastchange_to_oao_basis(decomp,A)
implicit none 
type(decompItem),intent(inout) :: decomp
type(Matrix) :: A, wrk

   call MAT_INIT(wrk,A%nrow,A%ncol)
   call mat_mul(decomp%U_inv,A,'n','n',1E0_realk,0E0_realk,wrk)
   call mat_mul(wrk,decomp%U_inv,'n','t',1E0_realk,0E0_realk,A)
   call mat_free(wrk)

end subroutine leastchange_to_oao_basis


function leastchange_propint(decomp,ls,lupri,luerr)
implicit none
type(decompItem),intent(inout) :: decomp
integer, parameter  :: nderiv=2, nMAT=10
TYPE(lsitem) :: ls
Real(realk),pointer :: leastchange_propint(:,:,:)
type(Matrix) :: propint(nMAT)
INTEGER               :: i, n,lupri,luerr
 
!     n   = ls%setting%BASIS(1)%p%REGULAR%nbast
     n = getNbasis(AORdefault,Contractedinttype,ls%SETTING%MOLECULE(1)%p,LUPRI)

     DO i=1,nMAT
        call mat_init(propint(I),n,n)
     ENDDO

     call II_get_carmom(lupri,luerr,ls%setting,propint,nMAT,nderiv,0E0_realk,0E0_realk,0E0_realk)

     call mat_free(propint(1))
     call mat_free(propint(6))
     call mat_free(propint(7))
     call mat_free(propint(9))

     call leastchange_to_oao_basis(decomp,propint(2))
     call leastchange_to_oao_basis(decomp,propint(3))
     call leastchange_to_oao_basis(decomp,propint(4))
     call leastchange_to_oao_basis(decomp,propint(5))
     call leastchange_to_oao_basis(decomp,propint(8))
     call leastchange_to_oao_basis(decomp,propint(10))

     allocate(leastchange_propint(n,n,6))

     do i=1,4
      call mat_to_full(propint(i+1),1E0_realk,leastchange_propint(:,:,i))
      call mat_free(propint(i+1))
     enddo
     
     call mat_to_full(propint(8),1E0_realk,leastchange_propint(:,:,5))
     call mat_free(propint(8))
     
     call mat_to_full(propint(10),1E0_realk,leastchange_propint(:,:,6))
     call mat_free(propint(10))

     
  return
end function leastchange_propint


subroutine leastchange_blocktransf(T,CMO,nbas,nocc,ncore)
implicit none
integer     :: nbas,ncore,nocc
real(realk) :: T(nbas,nbas), CMO(nbas,nbas)
integer     :: nval, nvirt,LWORK,INFO
real(realk),pointer :: U(:,:), S(:), VT(:,:),WORK(:), F(:,:),CB(:,:)
character*27 :: msg
INFO=0

 nvirt = nbas - nocc
 nval  = nocc - ncore

!core block
if (ncore.gt.0) then

 call mem_alloc(S,ncore)
 call mem_alloc(U,ncore,ncore)
 call mem_alloc(VT,ncore,ncore)
 call mem_alloc(CB,ncore,ncore)

 CB = CMO(1:ncore,1:ncore)

!optimal memory
 LWORK = -1
 call mem_alloc(WORK,5)
 call dgesvd('A','A',ncore,ncore,CB,ncore,S,U,ncore,VT,ncore,WORK,LWORK,INFO)
 LWORK=nint(WORK(1))
 call mem_dealloc(WORK)
 call mem_alloc(WORK,LWORK)
 call dgesvd('A','A',ncore,ncore,CB,ncore,S,U,ncore,VT,ncore,WORK,LWORK,INFO)

 call mem_dealloc(CB)
 call mem_dealloc(S)
 call mem_dealloc(WORK)

 if (INFO.ne. 0) then
   write(msg,'(A,I3)') 'dgesvd failed with INFO=',INFO
   call lsquit(msg,-1)
 endif

 
 call mem_alloc(F,ncore,ncore)

 call dgemm('T','T',ncore,ncore,ncore,1E0_realk,VT,ncore,U,ncore,0E0_realk,F,ncore)

 call mem_dealloc(U)
 call mem_dealloc(VT)

 call dgemm('N','N',nbas,ncore,ncore,1E0_realk,CMO(1,1),nbas,F,ncore,0E0_realk,T(1,1),nbas)

 call mem_dealloc(F)
endif
!valence block
!optimal memory

 call mem_alloc(S,nval)
 call mem_alloc(U,nval,nval)
 call mem_alloc(VT,nval,nval)
 call mem_alloc(CB,nval,nval)

 CB = CMO(ncore+1:ncore+nval,ncore+1:ncore+nval)

 call mem_alloc(WORK,5)
 LWORK = -1
 call dgesvd('A','A',nval,nval,CB,nval,S,U,nval,VT,nval,WORK,LWORK,INFO)
 LWORK=nint(WORK(1))
 call mem_dealloc(WORK)
 call mem_alloc(WORK,LWORK)
 call dgesvd('A','A',nval,nval,CB,nval,S,U,nval,VT,nval,WORK,LWORK,INFO)

 call mem_dealloc(CB)
 call mem_dealloc(S)
 call mem_dealloc(WORK)

 if (INFO.ne. 0) then
   write(msg,'(A,I3)') 'dgesvd failed with INFO=',INFO
   call lsquit(msg,-1)
 endif

 
 call mem_alloc(F,nval,nval)

 call dgemm('T','T',nval,nval,nval,1E0_realk,VT,nval,U,nval,0E0_realk,F,nval)

 call mem_dealloc(U)
 call mem_dealloc(VT)

 call dgemm('N','N',nbas,nval,nval,1E0_realk,CMO(1,ncore+1),nbas,F,nval,0E0_realk,T(1,ncore+1),nbas)

 call mem_dealloc(F)


!virtual block
!optimal memory
if (nvirt .gt. 0) then
 
 call mem_alloc(S,nvirt)
 call mem_alloc(U,nvirt,nvirt)
 call mem_alloc(VT,nvirt,nvirt)
 call mem_alloc(CB,nvirt,nvirt)

 CB = CMO(nocc+1:nbas,nocc+1:nbas)

 call mem_alloc(WORK,5)
 LWORK=-1
 call dgesvd('A','A',nvirt,nvirt,CB,nvirt,S,U,nvirt,VT,nvirt,WORK,LWORK,INFO)
 LWORK=nint(WORK(1))
 call mem_dealloc(WORK)
 call mem_alloc(WORK,LWORK)
 call dgesvd('A','A',nvirt,nvirt,CB,nvirt,S,U,nvirt,VT,nvirt,WORK,LWORK,INFO)

 call mem_dealloc(CB)
 call mem_dealloc(S)
 call mem_dealloc(WORK)

 if (INFO.ne. 0) then
   write(msg,'(A,I3)') 'dgesvd failed with INFO=',INFO
   call lsquit(msg,-1)
 endif

 
 call mem_alloc(F,nvirt,nvirt)

 call dgemm('T','T',nvirt,nvirt,nvirt,1E0_realk,VT,nvirt,U,nvirt,0E0_realk,F,nvirt)
 
 call mem_dealloc(U)
 call mem_dealloc(VT)

 call dgemm('N','N',nbas,nvirt,nvirt,1E0_realk,CMO(1,nocc+1),nbas,F,nvirt,0E0_realk,T(1,nocc+1),nbas)

 call mem_dealloc(F)
end if
end subroutine leastchange_blocktransf

subroutine leastchange_orbspread(spread,uindex,n,&
                     &T,DIPX,DIPY,DIPZ,SECX,SECY,SECZ)
implicit none
integer :: n, uindex(2*n)
real(realk) :: spread(n),T(n,n),DIPX(n,n),DIPY(n,n),DIPZ(n,n)
real(realk) :: SECX(n,n),SECY(n,n),SECZ(n,n)
real(realk),pointer :: tmp(:), tmpT(:,:)
real(realk) :: sec,dip
integer     :: i, itmp
real(realk), external :: ddot


 call leastchange_iuindex(uindex(n+1),uindex,n)

 call mem_alloc(tmp,n)
 call mem_alloc(tmpT,n,n)

 call leastchange_rowreorder(tmpT,T,uindex(n+1),n)

 T = tmpT

 call mem_dealloc(tmpT)


 

 do i=1,n
     call dsymv('u',n,1E0_realk,SECX,n,T(1,i),1,0E0_realk,tmp,1)
     sec = ddot(n,T(1,i),1,tmp,1)
     call dsymv('u',n,1E0_realk,DIPX,n,T(1,i),1,0E0_realk,tmp,1)
     dip = ddot(n,T(1,i),1,tmp,1)
     spread(i) = sec - dip*dip


     call dsymv('u',n,1E0_realk,SECY,n,T(1,i),1,0E0_realk,tmp,1)
     sec = ddot(n,T(1,i),1,tmp,1)
     call dsymv('u',n,1E0_realk,DIPY,n,T(1,i),1,0E0_realk,tmp,1)
     dip = ddot(n,T(1,i),1,tmp,1)
     spread(i) = spread(i) + sec - dip*dip
     

     call dsymv('u',n,1E0_realk,SECZ,n,T(1,i),1,0E0_realk,tmp,1)
     sec = ddot(n,T(1,i),1,tmp,1)
     call dsymv('u',n,1E0_realk,DIPZ,n,T(1,i),1,0E0_realk,tmp,1)
     dip = ddot(n,T(1,i),1,tmp,1)
     spread(i) = spread(i) + sec - dip*dip

     spread(i) = sqrt(spread(i))
 enddo

 call mem_dealloc(tmp)
end subroutine leastchange_orbspread

function leastchange_idmin(n,vec)
integer :: leastchange_idmin, n, i
real(realk) :: vec(n)
  
   leastchange_idmin=1
   do i=2,n
       if ((vec(leastchange_idmin)-vec(i))>1E-6_realk) leastchange_idmin=i
   enddo

return
end function leastchange_idmin


subroutine leastchange_rowreorder(B,A,uindex,n)
implicit none
integer :: n, i
integer :: uindex(n)
real(realk) :: A(n,n), B(n,n)

   do i=1,n
      call dcopy(n,A(uindex(i),1),n,B(i,1),n)
   enddo

end subroutine leastchange_rowreorder

subroutine leastchange_iuindex(iuindex,uindex,n)
integer :: n, i
integer :: iuindex(n),uindex(n)
  
   do i=1,n
      iuindex(uindex(i))=i
   enddo

end subroutine leastchange_iuindex   


subroutine leastchange_rowsort_diag(uindex,CMO,nocc,n)
implicit none
integer :: n,nocc, uindex(n)
real(realk) :: CMO(n,n)
integer, external :: idamax
integer :: i,j,tmp

 do i=1,nocc
       j=idamax(n-i+1,CMO(i,i),1)
       j=j+i-1
       if (j.gt.i) then
        call dswap(n, CMO(j,1),n,CMO(i,1),n)   
        tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
       endif
 enddo
 
end subroutine leastchange_rowsort_diag

subroutine leastchange_rowsort_occup(uindex,CMO,ncore,nocc,n)
implicit none
integer :: n,ncore,nocc, uindex(n)
real(realk) :: CMO(n,n)
real(realk),pointer :: occnum(:)
integer, external     :: idamax
real(realk), external :: ddot
integer :: i,j,tmp
real(realk) :: rtmp



 ! 1. Order orbitals according to occupied/virtual occupation
 ! --> "occupied" OAOs before "virtual" OAOs
 call mem_alloc(occnum,n)
 do i=1,n
   occnum(i) = ddot(nocc,CMO(i,1),n,CMO(i,1),n) 
 enddo

 do i=1,nocc
       j=idamax(n-i+1,occnum(i),1)
       j=j+i-1
       if (j.gt.i) then
        call dswap(n, CMO(j,1),n,CMO(i,1),n)   
        tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
        rtmp=occnum(i); occnum(i)=occnum(j); occnum(j)=rtmp
       endif
 enddo
 call mem_dealloc(occnum)


 ! 2. Order orbitals according to core/valence occupation
 ! --> "core" OAOs before "valence" OAOs
 call mem_alloc(occnum,nocc)
 do i=1,nocc
   occnum(i) = ddot(ncore,CMO(i,1),n,CMO(i,1),n) 
 enddo

 do i=1,ncore
       j=idamax(nocc-i+1,occnum(i),1)
       j=j+i-1
       if (j.gt.i) then
        call dswap(n, CMO(j,1),n,CMO(i,1),n)   
        tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
        rtmp=occnum(i); occnum(i)=occnum(j); occnum(j)=rtmp
       endif
 enddo
 call mem_dealloc(occnum)
 
end subroutine leastchange_rowsort_occup



subroutine leastchange_rowsort_bf(decomp,uindex,T,CMO,ncore,nocc,n,ls,&
                                 &lupri,luerr)
implicit none
TYPE(lsitem) :: ls
type(decompItem),intent(inout) :: decomp
real(realk) :: T(n,n), CMO(n,n), spread(n)
real(realk),pointer :: mx(:)
real(realk),pointer :: DIPX(:,:), DIPY(:,:), DIPZ(:,:)
real(realk),pointer :: SECX(:,:), SECY(:,:), SECZ(:,:)
Real(realk),pointer :: propint(:,:,:)
integer :: uindex(2*n), tmp, changes, ncore,nocc,nvirt,n,idx
integer, external :: idamax
real(realk),external :: ddot
integer :: i,j,lupri,luerr,iter
data iter /0/


  propint => leastchange_propint(decomp,ls,lupri,luerr)

  DIPX => propint(:,:,1)
  DIPY => propint(:,:,2)
  DIPZ => propint(:,:,3)
  SECX => propint(:,:,4)
  SECY => propint(:,:,5)
  SECZ => propint(:,:,6)

  !call dumpmat('DIPX.bin',DIPX,8*n*n)
  !call dumpmat('DIPY.bin',DIPY,8*n*n)
  !call dumpmat('DIPZ.bin',DIPZ,8*n*n)
  !call dumpmat('SECX.bin',SECX,8*n*n)
  !call dumpmat('SECY.bin',SECY,8*n*n)
  !call dumpmat('SECZ.bin',SECZ,8*n*n)

  nvirt = n - nocc

  call mem_alloc(mx,nvirt+1)

 do 

  changes=0

  do i=ncore+1,nocc

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE(j,tmp,spread,T) &
!$OMP FIRSTPRIVATE(CMO,uindex) SCHEDULE(DYNAMIC,1) 
     do j=nocc,n

       if (j.gt.nocc) then
        call dswap(n, CMO(j,1),n,CMO(i,1),n)   
        tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
       endif

       call leastchange_blocktransf(T,CMO,n,nocc,ncore)
       call leastchange_orbspread(spread,uindex,n,T,&
           &DIPX,DIPY,DIPZ,SECX,SECY,SECZ)

       if (j.gt.nocc) then
        call dswap(n, CMO(j,1),n,CMO(i,1),n)   
        tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
       endif

       mx(j-nocc+1) = spread(idamax(n,spread,1))

     enddo
!$OMP END PARALLEL DO

     idx=leastchange_idmin(nvirt+1,mx) 
     j= idx + nocc - 1

     if (j.gt.nocc) then
       call dswap(n, CMO(i,1),n,CMO(j,1),n) 
       tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
       changes = changes+1
     endif

  enddo

 iter = iter + 1
 write(*,'(2X,I4,1X,A,I4,1X,A,F8.5)') &
 &     iter, 'Changes= ', changes, 'Max spread= ', mx(idx)

 if (changes.eq. 0) exit

 enddo

 call mem_dealloc(propint)
 call mem_dealloc(mx)

end subroutine leastchange_rowsort_bf

end module leastchange_module

subroutine leastchangeOrbspreadStandalone(mx,ls,CMO,lupri,luerr)
use precision
use TYPEDEFTYPE, only: lsitem
use IntegralInterfaceMOD
use matrix_operations
use matrix_module
use memory_handling
implicit none
real(realk)              :: mx
TYPE(lsitem)             :: ls
type(Matrix)             :: CMO
integer,intent(in)       :: lupri, luerr

integer                  :: n, i, nel, nocc,nvirt,occidx,virtidx
integer, parameter       :: nderiv=2, nMAT=10
type(Matrix)             :: propint(nMAT), PROPT
real(realk)              :: sec,dip, maxocc, maxvirt
real(realk), pointer :: T(:,:), PROPTf(:,:)
real(realk), pointer :: spread(:)
real(realk), external    :: ddot
integer, external        :: idamax


  n=CMO%nrow

  DO i=1,nMAT
     call mat_init(propint(I),n,n)
  ENDDO
     call II_get_carmom(lupri,luerr,ls%setting,propint,nMAT,nderiv,0E0_realk,0E0_realk,0E0_realk)

! Free what we do not need
     call mat_free(propint(1))
     call mat_free(propint(6))
     call mat_free(propint(7))
     call mat_free(propint(9))

! SEC = SECX + SECY + SECZ
     call mat_daxpy(1E0_realk,propint(8),propint(5))
     call mat_free(propint(8))
     call mat_daxpy(1E0_realk,propint(10),propint(5))
     call mat_free(propint(10))

! SEC, prepare T and PROPT=SEC*T 
     call mat_init(PROPT,n,n)
     call mat_mul(propint(5),CMO,'n','n',1E0_realk,0E0_realk,PROPT)
     call mat_free(propint(5))

     call mem_alloc(PROPTf,n,n)
     call mat_to_full(PROPT,1E0_realk,PROPTf)
     call mat_free(PROPT)

     call mem_alloc(T,n,n)
     call mat_to_full(CMO,1E0_realk,T)

     call mem_alloc(spread,n)

 do i=1,n
     sec = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = sec
 enddo


!DIPX
     call mat_init(PROPT,n,n)
     call mat_mul(propint(2),CMO,'n','n',1E0_realk,0E0_realk,PROPT)
     call mat_free(propint(2))

     call mat_to_full(PROPT,1E0_realk,PROPTf)
     call mat_free(PROPT)

 do i=1,n
     dip = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = spread(i) - dip*dip
 enddo

!DIPY
     call mat_init(PROPT,n,n)
     call mat_mul(propint(3),CMO,'n','n',1E0_realk,0E0_realk,PROPT)
     call mat_free(propint(3))

     call mat_to_full(PROPT,1E0_realk,PROPTf)
     call mat_free(PROPT)

 do i=1,n
     dip = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = spread(i) - dip*dip
 enddo

!DIPZ
     call mat_init(PROPT,n,n)
     call mat_mul(propint(4),CMO,'n','n',1E0_realk,0E0_realk,PROPT)
     call mat_free(propint(4))

     call mat_to_full(PROPT,1E0_realk,PROPTf)
     call mat_free(PROPT)

 do i=1,n
     dip = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = spread(i) - dip*dip
     
     spread(i) = sqrt(spread(i))
 enddo

 call mem_dealloc(T)
 call mem_dealloc(PROPTf)

 mx = spread(idamax(n,spread,1))

 ! Number of electrons
 nel = ls%input%molecule%nelectrons
 ! number of occupied and virtual orbitals (assuming even number of electrons)
 nocc = nel/2
 nvirt = n-nocc
 ! Indices for orbitals with max orb. spread
 occidx = idamax(nocc,spread(1:nocc),1)
 virtidx = nocc + idamax(nvirt,spread(nocc+1:n),1)
 ! Max occ and virt orbital spreads
 maxocc = spread(occidx)
 maxvirt = spread(virtidx)

 write(lupri,*)
 write(lupri,*) 'Number of electrons  in molecule = ', nel
 write(lupri,*) 'Number of occ. orb.  in molecule = ', nocc
 write(lupri,*) 'Number of virt. orb. in molecule = ', nvirt
 write(lupri,*)
 write(lupri,*)

 !write(lupri,*) 'Orbspreads for occupied orbitals:'
 !write(lupri,*) '*********************************'
 !do i=1,nocc
 !   write(lupri,*) i, spread(i)
 !end do

 !write(lupri,*)
 !write(lupri,*)
 !write(lupri,*) 'Orbspreads for virtual orbitals:'
 !write(lupri,*) '********************************'
 !do i=nocc+1,n
 !   write(lupri,*) i, spread(i)
 !end do
 !write(lupri,*)
 !write(lupri,*)
 !write(lupri,*) 'Maximum orb. spread (MOS) summary'
 !write(lupri,*) '*********************************'

 write(lupri,*) 'OCCUPIED: (Orbital index, MOS) = ', occidx, maxocc
 write(lupri,*) 'VIRTUAL : (Orbital index, MOS) = ', virtidx, maxvirt
 write(lupri,*)
 write(lupri,*)

 call mem_dealloc(spread)

end subroutine leastchangeOrbspreadStandalone



subroutine leastchange_lcv(decomp,CMO,nocc,ls)
use precision
use leastchange_module
use decompMod
use loc_utils
use matrix_module
use matrix_operations
use memory_handling
implicit none
type(decompItem) :: decomp
type(Matrix) :: CMO
type(lsitem) :: ls
type(Matrix) :: tmp
integer :: n, nocc, ncore, i
integer, allocatable :: uindex(:)
real(realk),pointer :: CMOf(:,:), Tf(:,:)

  print*,'inside leastchange_lcv'
  !convert CMO to orthonormal basis
  n = CMO%nrow
  call mat_init(tmp,n,n)
  call mat_mul(decomp%U,CMO,'n','n',1E0_realk,0E0_realk,tmp)

  call mem_alloc(CMOf,n,n)
  call mat_to_full(tmp,1E0_realk,CMOf)
  call mat_free(tmp)
  !get ncore
  ncore = count_ncore(ls)

  !sorting index
  allocate(uindex(2*n))
  do i=1,n
    uindex(i)=i
  enddo

  !occupation sorting
  call leastchange_rowsort_occup(uindex,CMOf,ncore,nocc,n)

  !brute force sorting
  call mem_alloc(Tf,n,n)
  if (decomp%cfg_lcvbf) &
 &call leastchange_rowsort_bf(decomp,uindex,Tf,CMOf,ncore,nocc,n,ls,&
                            &decomp%lupri,decomp%luerr)

  !lcv basis
  call leastchange_blocktransf(Tf,CMOf,n,nocc,ncore)

  !inverse sorting index
  call leastchange_iuindex(uindex(n+1),uindex,n)

  !reorder to original
  call leastchange_rowreorder(CMOf,Tf,uindex(n+1),n)


  !back to gcscf basis
  call mem_dealloc(Tf)
  deallocate(uindex)


  call mat_init(tmp,n,n)
  call mat_set_from_full(CMOf,1E0_realk,tmp)

  call mem_dealloc(CMOf)
 
  call mat_mul(decomp%U_inv,tmp,'n','n',1E0_realk,0E0_realk,CMO)
  
  call mat_free(tmp)


end subroutine leastchange_lcv

subroutine leastchange_lcm(decomp,CMO,nocc,ls)
use leastchange_module
use decompMod
use loc_utils
use precision
use matrix_module
use matrix_operations
use memory_handling
implicit none
type(decompItem),intent(in) :: decomp
type(Matrix) :: CMO
type(lsitem) :: ls
type(Matrix) :: tmp
integer :: n, nocc, ncore
real(realk),pointer :: CMOf(:,:), Tf(:,:)

  write(decomp%lupri,*)'inside leastchange_lcm'
  print*,'inside leastchange_lcm'
  !convert CMO to orthonormal basis
  n = CMO%nrow
  call mat_init(tmp,n,n)
  call mat_mul(decomp%U,CMO,'n','n',1E0_realk,0E0_realk,tmp)
  call mem_alloc(CMOf,n,n)
  call mat_to_full(tmp,1E0_realk,CMOf)
  call mat_free(tmp)

  !get ncore
  ncore = count_ncore(ls)

  call mem_alloc(Tf,n,n)

  !lcm basis
  call leastchange_blocktransf(Tf,CMOf,n,nocc,ncore)

  call mem_dealloc(CMOf)

  !back to gcscf basis

  call mat_init(tmp,n,n)
  call mat_set_from_full(Tf,1E0_realk,tmp)

  call mem_dealloc(Tf)
 
  call mat_mul(decomp%U_inv,tmp,'n','n',1E0_realk,0E0_realk,CMO)
  
  call mat_free(tmp)

end subroutine leastchange_lcm



!> \brief Get normalized projected atomic orbitals (PAOs),
!> where the occupied space have been projected out.
!> I.e. the nbasis PAOs is a redundant set of
!> local orbitals which span the virtual space.
!> Only used for testing.
!> \author Kasper Kristensen
!> \date November 2011
subroutine get_PAOs(Cmo,S,MyLsitem,lupri)
  use matrix_util
  use files
  use TYPEDEFTYPE, only: lsitem
  use matrix_module
  use precision
  use matrix_operations
  use memory_handling

  implicit none
  !> Optimized MO coeffiecient (both occupied and virtual)
  type(matrix), intent(in) :: Cmo
  !> Overlap matrix in AO basis
  type(matrix), intent(in) :: S
  !> LSDALTON info
  type(lsitem), intent(inout) :: MyLsitem
  !> File unit number for DALTON.OUT
  integer, intent(in) :: lupri
  integer :: nelectrons
  type(matrix) :: DS, UnitMatrix, Q, QtSQ,D
  type(matrix) :: PAOs, PAO_overlap,Cocc,Cvirt
  type(matrix) :: SPCocc, SPCvirt
  integer :: nbasis, i,j, idx, nstart, nend, funit, nocc,nunocc, nunocc_end,lun
  real(realk), pointer :: Normalize_coeff(:), occ_vector(:), unocc_vector(:)
  real(realk), pointer :: eival(:)
  real(realk) :: coeff, tmp, occ_sum, unocc_sum,mx
  logical :: file_exist,OnMaster
  OnMaster=.TRUE.

  ! Initialize stuff
  ! ****************

  ! Number of electrons (assuming neutral molecule)
  nelectrons=0
  do i=1,MyLsitem%setting%MOLECULE(1)%p%nAtoms
     IF(MyLsitem%setting%MOLECULE(1)%p%Atom(i)%phantom)CYCLE
     nelectrons = nelectrons + MyLsitem%setting%MOLECULE(1)%p%Atom(i)%Charge
  enddo
  nbasis = cmo%nrow
  nocc = nelectrons/2
  nunocc = nbasis-nocc
  call mat_init(DS,nbasis,nbasis)
  call mat_init(UnitMatrix,nbasis,nbasis)
  call mat_init(Q,nbasis,nbasis)


  ! Occupied orbitals
  ! *****************
  call mat_init(Cocc,nbasis,nocc)
  nstart=1
  nend = nbasis*nocc
  Cocc%elms(nstart:nend) = cmo%elms(nstart:nend)


  ! Unoccupied orbitals
  ! *******************
  call mat_init(Cvirt,nbasis,nunocc)
  nstart = nbasis*nocc + 1
  nend = nbasis*nbasis
  nunocc_end = nbasis*nunocc
  Cvirt%elms(1:nunocc_end) = Cmo%elms(nstart:nend)


  ! Get density matrix from MOs
  ! ***************************
  ! Density matrix D = Cocc*Cocc^T 
  call mat_init(D,nbasis,nbasis)
  call mat_mul(Cocc,Cocc, 'n', 't', 1.d0, 0.d0,D)




  ! Project out occupied part of atomic orbitals
  ! ********************************************

  ! DS
  call mat_mul(D,S, 'n', 'n', 1.d0, 0.d0,DS)

  ! Unit matrix
  call mat_zero(UnitMatrix)
  do i=1,nbasis
     ! idx is the position in the matrix vector containing diagonal entry (i,i)
     idx = nbasis*(i-1) + i
     UnitMatrix%elms(idx) = 1.d0
  end do

  ! Q = 1-DS
  call mat_add(1.d0,UnitMatrix,-1.d0,DS,Q)
  call mat_free(UnitMatrix)
  call mat_free(DS)
  call mat_free(D)

  ! Q now contains the (non-normalized) expansion coeffisions for the PAOs



  ! Normalize PAOs
  ! **************
  call mat_init(QtSQ,nbasis,nbasis)

  ! Calculate Q^T S Q
  call util_AO_to_MO_2(S,Q,S,QtSQ,.true.)

  call mem_alloc(Normalize_coeff,nbasis)

  do i=1,nbasis
     ! idx is the position in the matrix vector containing diagonal entry (i,i)
     idx = nbasis*(i-1) + i
     ! Normalize coeffient for PAO_mu: 1/sqrt[ (S_SDS)_{mu,mu} ]
     if(QtSQ%elms(idx) < 0.0) then
        call lsquit('get_PAOs: Squared normalization factor smaller than zero,&
             & something is wrong!',-1)
     end if
     Normalize_coeff(i) = 1.d0/sqrt(QtSQ%elms(idx))
  end do
  call mat_free(QtSQ)


  ! For all orbitals, multiply orbital coefficients by normalization coefficients
  call mat_init(PAOs,nbasis,nbasis)
  call mat_zero(PAOs)
  do i=1,nbasis
     coeff = Normalize_coeff(i)
     do j=1,nbasis
        nstart  = (i-1)*nbasis + 1
        nend = i*nbasis
        PAOs%elms(nstart:nend) = coeff*Q%elms(nstart:nend) 
     end do
  end do
  call mat_free(Q)


  ! PAOs now contains the normalized projected atomic orbital coefficients.
  lun=-1
  CALL LSOPEN(lun,'pao','replace','UNFORMATTED')
  call mat_write_to_disk(lun,PAOs,OnMaster)
  call LSclose(lun,'KEEP')


  ! *******************************************************************************
  !                                Check PAOs                                     *
  ! *******************************************************************************

  ! Overlap matrix: PAO^T S PAO
  ! ***************************
  call mat_init(PAO_overlap,nbasis,nbasis)
  call util_AO_to_MO_2(S,PAOs,S,PAO_overlap,.true.)



  ! PAO/occupied overlap
  ! ********************

  ! Construct PAO/occupied overlap matrix: SPCocc = PAO^T S Cocc_can
  call mat_init(SPCocc,nbasis,nocc)
  call util_AO_to_MO_different_trans(PAOs, S, Cocc, SPCocc)
  call mat_free(Cocc)

  ! Construct PAO/virtual overlap matrix: SPCvirt = PAO^T S Cvirt_can
  call mat_init(SPCvirt,nbasis,nunocc)
  call util_AO_to_MO_different_trans(PAOs, S, Cvirt, SPCvirt)
  call mat_free(Cvirt)



  ! Vectors containing PAOs projected against occ and virt spaces
  ! *************************************************************

  call mem_alloc(occ_vector,nbasis)
  call mem_alloc(unocc_vector,nbasis)
  occ_vector=0.d0
  unocc_vector=0.d0


  ! Occupied vector
  ! '''''''''''''''
  do i=1,nbasis
     do j=1,nocc
        idx = (j-1)*nbasis + i
        tmp = SPCocc%elms(idx)
        ! occ_vector(i) = sum_j | <PAO(i) | Cocc_can(j)> |^2
        occ_vector(i) = occ_vector(i) + tmp**2
     end do
  end do
  call mat_free(SPCocc)


  ! Virtual vector
  ! ''''''''''''''
  do i=1,nbasis
     do j=1,nunocc
        idx = (j-1)*nbasis + i
        tmp = SPCvirt%elms(idx)
        ! unocc_vector(i) = sum_j | <PAO(i) | Cvirt_can(j)> |^2
        unocc_vector(i) = unocc_vector(i) + tmp**2
     end do
  end do
  call mat_free(SPCvirt)


  ! Maximum orbital spreads
  ! ***********************

  ! Get maximum orbital spreads for PAOs
  ! (only the maximum number in MaxOrbSpreads makes sense here 
  !  because we consider only the virtual space)
  write(lupri,*) 
  write(lupri,*) 'Calculating orbital spreads for projected atomic orbitals'
  write(lupri,*) 'NOTE: The orb spreads printed below ALL refer to PAOs!'
  write(lupri,*) '      - - The occupied/virtual partitioning is meaningless for the redundant PAOs'
  write(lupri,*) 
  call leastchangeOrbspreadStandalone(mx,MyLSitem,PAOs,lupri,lupri)
  call mat_free(PAOs)



  write(lupri,*) 
  write(lupri,*) 
  write(lupri,*) 
  write(lupri,'(5X,a)') '****************************************************************'
  write(lupri,'(5X,a)') '*             Projected atomic orbitals check                  * '
  write(lupri,'(5X,a)') '****************************************************************'
  write(lupri,*) 
  write(lupri,*) 


  write(lupri,'(5X,a)') 'Checking that PAOS are normalized'
  write(lupri,'(5X,a)') '---------------------------------'

  write(lupri,'(9X,a)') 'Orbital         Norm'
  do i=1,nbasis
     idx = nbasis*(i-1) + i
     write(lupri,'(5X,i8,4X,g18.8)') i, PAO_overlap%elms(idx)
  end do


  write(lupri,*) 
  write(lupri,*) 
  write(lupri,'(5X,a)') 'Checking that PAOS span virtual space'
  write(lupri,'(5X,a)') '-------------------------------------'

  write(lupri,'(9X,a,6X,a,10X,a)') 'Orbital', 'Occ. proj.', 'Virt. proj.'
  occ_sum=0.d0
  unocc_sum=0.d0
  do i=1,nbasis
     write(lupri,'(5X,i8,3X,g18.8,3X,g18.8)') i, occ_vector(i), unocc_vector(i)
     occ_sum = occ_sum + occ_vector(i)
     unocc_sum = unocc_sum + unocc_vector(i)
  end do
  write(lupri,'(20X,a)') "'''''''''''''''''''''''''''''''"
  write(lupri,'(10X,a,2X,g18.8,3X,g18.8)') 'SUM:',occ_sum, unocc_sum




  ! Free stuff
  call mem_dealloc(Normalize_coeff)
  call mem_dealloc(occ_vector)
  call mem_dealloc(unocc_vector)
  call mat_free(PAO_overlap)

end subroutine get_PAOs
