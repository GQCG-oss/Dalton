!----------------------------------------------------------------------
subroutine purify(lupri,F,S,nocc,method,Dnew)
!----------------------------------------------------------------------
!density matrix purification routine (Andreas Hesselmann, 2005) 
!circular dependencies and portability fixed by Pawel Salek.
!(based on dmm.py from PyQuante package)
use matrix_operations
  implicit none
  integer,intent(in)  :: lupri
  integer, intent(in) :: nocc, method
  type(Matrix), intent(in)  :: F,S
  type(Matrix), intent(inout) :: Dnew
  type(Matrix) :: L,Z,Ft,Id,D,D2,Df,Dp,Dp2,Dg,D3,T,Dbest
  integer, parameter :: maxiter=100
  real(realk), parameter :: thr=1E-12_realk,dev_factor=20E0_realk
  character(3), dimension(4), parameter :: methods=(/'TCP','TRS','MCW','PM '/)
  integer :: ndim,iter
  real(realk) :: emin,emax,Ne_curr,efermi,alpha,beta,gamma,Dsum,Dsumold, &
               & trf,trg,cn,trace_error,trace_error_best
  real(realk),external :: mat_sum

  write(lupri,'(1x,a)') "DM purification will be performed:"
  write(lupri,'(1x,a,a)') "Method: ",methods(method)
  ndim=S%nrow

  !transform F into orthonormal basis
  call mat_init(L,ndim,ndim)
  call mat_init(Z,ndim,ndim)
  call mat_init(Ft,ndim,ndim)
  call mat_cholesky(S,L)
  call mat_inverse_triang(L,Z)
  !call mat_to_minus_one_half(S,Z) 
  !call mat_simtran(F,Z,'N',Ft)
  call mat_mul(Z,   F,'N','N',1E0_realk,0E0_realk,Dnew)
  call mat_mul(Dnew,Z,'N','T',1E0_realk,0E0_realk,Ft)

  !initialise the density matrix
  call mat_init(Id,ndim,ndim)
  call mat_init(D,ndim,ndim)
  call mat_gershgorin_minmax(Ft,emin,emax)
  call mat_identity(Id)
  if(method.eq. 1.or.method.eq. 2) then
     call mat_add(emax/(emax-emin),Id,-1E0_realk/(emax-emin),Ft,D)
  elseif(method.eq. 3) then
     call Dinit_mcw(Ft,nocc,D)
  elseif(method.eq. 4) then
     efermi=mat_tr(Ft)/dble(ndim)
     beta=dble(nocc)/dble(ndim)
     alpha=min(dble(nocc)/(emax-efermi),dble(ndim-nocc)/(efermi-emin))
     alpha=alpha/dble(ndim)
     call mat_add(alpha*efermi,Id,-alpha,Ft,D)
     call mat_add(1E0_realk,D,beta,Id,D)
  else
     stop 'unknown method in purify!'
  endif

  !initialise matrices needed in the purification
  call mat_init(D2,ndim,ndim)
  call mat_init(Dbest,ndim,ndim)
  if(method.eq. 2) then
     call mat_init(Df,ndim,ndim)
     call mat_init(Dp,ndim,ndim)
     call mat_init(Dp2,ndim,ndim)
     call mat_init(Dg,ndim,ndim)
  elseif(method.eq. 3) then
     call mat_init(D3,ndim,ndim)
  elseif(method.eq. 4) then
     call mat_init(D3,ndim,ndim)
     call mat_init(T,ndim,ndim)
     Dsumold=mat_sum(D)
  endif

  !iterate on D updates (purify!)
  do iter=1,maxiter
     Ne_curr=mat_tr(D)
     trace_error=abs(Ne_curr-nocc)
     if(method.ne. 4) then
        write(lupri,'(1x,a,i4,3x,a,f16.10)') 'iter=',iter,'trace=',Ne_curr
        if(trace_error.lt.thr) exit
     endif
     call mat_mul(D,D,'N','N',1E0_realk,0E0_realk,D2)
     if(method.eq. 1) then
        if(Ne_curr.lt.nocc) then
           call mat_add(2E0_realk,D,-1E0_realk,D2,D)
        else
           call mat_copy(1E0_realk,D2,D)
        endif
     elseif(method.eq. 2) then
        call mat_add(4E0_realk,D,-3E0_realk,D2,Dp)
        call mat_mul(D2,Dp,'N','N',1E0_realk,0E0_realk,Df)
        trf=mat_tr(Df)
        call mat_add(1E0_realk,Id,-1E0_realk,D,Dp)
        call mat_mul(Dp,Dp,'N','N',1E0_realk,0E0_realk,Dp2)
        call mat_mul(D2,Dp2,'N','N',1E0_realk,0E0_realk,Dg)
        trg=mat_tr(Dg)
        gamma=dble(nocc-trf)/trg
        if(gamma.gt. 2E0_realk) then
           call mat_add(2E0_realk,D,-1E0_realk,D2,D)
        elseif(gamma.lt. 0E0_realk) then
           call mat_copy(1E0_realk,D2,D)
        else
           call mat_add(1E0_realk,Df,-gamma,Dg,D)
        endif
     elseif(method.eq. 3) then
        if(trace_error.gt. 0.5E0_realk) then
           !the McW purification can converge to wrong D if the occupation
           !numbers are out of bound (usually at beginning), so do TCP
           !iteration if far from converged D  
           if(Ne_curr.lt.nocc) then
              call mat_add(2E0_realk,D,-1E0_realk,D2,D)
           else
              call mat_copy(1E0_realk,D2,D)
           endif
        else
           call mat_mul(D,D2,'N','N',1E0_realk,0E0_realk,D3)
           call mat_add(3E0_realk,D2,-2E0_realk,D3,D)
        endif
     elseif(method.eq. 4) then
        call mat_mul(D,D2,'N','N',1E0_realk,0E0_realk,D3)
        call mat_add(1E0_realk,D2,-1E0_realk,D3,T)
        cn=mat_tr(T)
        call mat_add(1E0_realk,D,-1E0_realk,D2,T)
        cn=cn/mat_tr(T)
        if(cn.lt. 0.5E0_realk) then
           call mat_add(1E0_realk-2E0_realk*cn,D,1E0_realk+cn,D2,D)
           call mat_add(1E0_realk,D,-1E0_realk,D3,D)
           call mat_scal(1E0_realk/(1E0_realk-cn),D)
        else
           call mat_add((1E0_realk+cn)/cn,D2,-1E0_realk/cn,D3,D)
        endif
        Dsum=mat_sum(D)
        write(lupri,'(1x,a,i4,3x,a,f16.10)') &
             & 'iter=',iter,'sum[D(n+1)-D(n)]=',Dsum-Dsumold
        if((Dsum-Dsumold).lt.thr) exit
        Dsumold=Dsum
        cycle
     endif

     if(iter.eq. 1) then
        trace_error_best=trace_error
        call mat_copy(1E0_realk,D,Dbest)
     else
        if(trace_error.lt.trace_error_best) then
           !save this up to now best density matrix
           trace_error_best=trace_error
           call mat_copy(1E0_realk,D,Dbest)
        elseif(Ne_curr.lt. 0.or.Ne_curr.gt.dev_factor*nocc) then
           !go here when purification shows no convergency
           call mat_copy(1E0_realk,Dbest,D)
           exit
        endif
     endif
          
  enddo
  write(lupri,'(1x,a,i5,1x,a)') 'purification converged in',iter,'iterations'
  if(iter.eq.maxiter) then
     call mat_copy(1E0_realk,Dbest,D)
  endif
  call mat_free(Dbest)
  
  !transform density back to nonorthogonal AO basis
  !call mat_simtran(D,Z,'T',Dnew)
  call mat_mul(Z,D,'T','N',1E0_realk,0E0_realk,L)
  call mat_mul(L,Z,'N','N',1E0_realk,0E0_realk,Dnew)

  call mat_free(L)
  call mat_free(Z)
  call mat_free(Id)
  call mat_free(D)
  call mat_free(D2)
  call mat_free(Ft)
  if(method.eq. 2) then
     call mat_free(Df)
     call mat_free(Dp)
     call mat_free(Dp2)
     call mat_free(Dg)
  elseif(method.eq. 3) then
     call mat_free(D3)
  elseif(method.eq. 4) then
     call mat_free(D3)
     call mat_free(T)
  endif

  return
end subroutine purify


!----------------------------------------------------------------------
subroutine Dinit_mcw(F,nocc,D)
use matrix_operations
!----------------------------------------------------------------------
  implicit none
  type(Matrix), intent(in) :: F
  integer, intent(in)      :: nocc
  type(Matrix), intent(inout) :: D
  integer, parameter :: maxit=100
  real(realk), parameter :: tol=1E-7_realk,beta=0.5E0_realk
  integer :: ndim,i
  real(realk) :: emin,emax,elow,ehigh,de,alpha,nelow,nehigh,efermi, &
              &  nefermi
  type(Matrix) :: Id,C

  ndim=F%nrow
  call mat_init(Id,ndim,ndim)
  call mat_identity(Id)

  call mat_gershgorin_minmax(F,emin,emax)
  elow=emin
  ehigh=emax+20E0_realk
  de=emax-elow
  alpha=beta/de
  call mat_add(alpha*elow,Id,-alpha,F,D)
  call mat_add(1E0_realk,D,beta,Id,D)
  nelow=mat_tr(D)
  call mat_add(alpha*ehigh,Id,-alpha,F,D)
  call mat_add(1E0_realk,D,beta,Id,D)
  nehigh=mat_tr(D)
  
  do i=1,maxit
     efermi=0.5E0_realk*(elow+ehigh)
     call mat_add(alpha*efermi,Id,-alpha,F,D)
     call mat_add(1E0_realk,D,beta,Id,D)
     nefermi=mat_tr(D)
     if(abs(nocc-nefermi).lt.tol) exit
     if(nefermi.lt.nocc) then
        elow=efermi
        nelow=nefermi
     elseif(nefermi.gt.nocc) then
        ehigh=efermi
        nehigh=nefermi
     endif
  enddo
  alpha=min(beta/(emax-efermi),(1E0_realk-beta)/(efermi-emin))
  call mat_add(alpha*efermi,Id,-alpha,F,D)
  call mat_add(1E0_realk,D,beta,Id,D)
  
  call mat_free(Id)

  return
end subroutine Dinit_mcw

  
