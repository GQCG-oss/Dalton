!> @file 
!> Contains DSM extra term module.

!> Module for adding an extra, more expensive, term to the DSM function.
module dsm_x_term
   use av_utilities
   use matrix_util
   private
   public :: xterm_start,xterm_build_Gpara,&
        &    xterm_add_Eterm,xterm_end !,xterm_build_deriv, xterm_Eps
   !> G(Delta_par) = sum_i coef_i (F_i - h)
   type(matrix),save :: Gpara
   !> coefficients used in sum for Gpara
   real(realk),allocatable,save :: coef(:)
contains
!> \brief Initialize xterm module.
!> \author L. Thogersen
!> \date 2003
!> \param av Contains info about SCF averaging
!> \param ndim Number of basis functions
!> \param msize Number of previous F/D matrices
   subroutine xterm_start(av,ndim,msize)
   implicit none
     type(avItem),intent(in) :: av
     integer, intent(in) :: ndim,msize
     WRITE(av%lupri,*) 'xterm mode with ',msize,' vectors    grep'
     call mat_init(Gpara,ndim,ndim)
     ALLOCATE(coef(msize))
   end subroutine xterm_start

!> \brief G(Delta_par) = sum_i coef_i (F_i - h)
!> \param av Contains info about SCF averaging
   subroutine xterm_build_Gpara(av,msize,queue,Delta,S,H1)
      implicit none
      type(avItem),intent(inout) :: av
      integer, intent(in) :: msize
      TYPe(util_HistoryStore)  :: queue
      type(Matrix), intent(in) :: Delta,S,H1
      integer :: i
      call util_GET_PROJ_PART(av,queue,Delta,S,coef)
      call get_average_arr(av,'F',queue, av%dsm_history_size, coef, Gpara)
      !correction to F if sum_i c_i /= 1
      !and withdrawel of the one-electron part of F: G=F-h
      call mat_daxpy(-SUM(coef(1:msize)),H1,Gpara)
!
! Remember to free Gpara at the end!!
!
   end subroutine xterm_build_Gpara

!> \brief Add the term Tr (2Delta - Delta_parr) G(Delta_parr)  to DSMener
!> \param av Contains info about SCF averaging
   subroutine xterm_add_Eterm(av,msize,queue,Delta,Dav,S,DSMener)
      !use fock_evaluator
      !use diagonalization, only: get_didem 
      implicit none
      type(avItem),intent(inout) :: av
      integer, intent(in) :: msize
      TYPe(util_HistoryStore)  :: queue
      type(Matrix), intent(in) :: Delta,Dav,S
      real(realk), intent(inout) :: DSMener
      type(Matrix) :: Dpara,wrk
      real(realk) :: EHFcur,dGd
      integer :: ndim,i
      ndim = Delta%nrow

      call mat_init(Dpara,ndim,ndim)
      call mat_init(wrk,ndim,ndim)
      call get_average_arr(av,'D',queue, av%dsm_history_size, coef, Dpara)
      call mat_add(2E0_realk,Delta,-1E0_realk,Dpara,wrk)
      dGd = mat_dotproduct(wrk,Gpara)
      DSMener = DSMener + dGd
      !if(.false.) then
      !   call mat_add(1E0_realk,Delta,-1E0_realk,Dpara,wrk)
      !   WRITE(av%LUPRI,"('||d_orth||,||d_para||,||d||',3f10.4,'   grep')") &
      !        & util_Snorm(wrk,S),util_Snorm(Dpara,S), util_Snorm(Delta,S)
      !   if (ABS(util_Snorm(Delta,S)) > 1E-10_realk) then 
      !     WRITE(av%LUPRI,"('ratio: ||d_orth||/||d||',f10.4,'    grep')") util_Snorm(wrk,S)/util_Snorm(Delta,S)
      !   endif
      !   call GET_Didem(Dav,S,Dpara)  !Fcur used as scr
      !   !** Find new Fock matrix and SCF energy
      !   call fck_get_fock(Dpara,wrk,EHFcur)
      !   WRITE(av%LUPRI,"('EDSM,EDSM+TrdG(d),ESCF',3f16.10,' grep')") DSMener-dGd, DSMener,EHFcur
      !endif
      call mat_free(Dpara)
      call mat_free(wrk)

   end subroutine xterm_add_Eterm

!> \brief Hessian contribution for extra term.
!>  
!>  Some of the Hessian terms for the "extra" term Tr(2Delta - Delta_par)G(Delta_par) are evaluated here
!>  on the fly. \n
!> The Hes-terms are: 2Trd²Delta/dc_i dc_j G(Delta_para) - 2Trd²Delta_para/dc_i dc_j G(Delta_para)
!>                                 +2TrDelta d²G(Delta_para)/dc_i dc_j
!>
!      subroutine xterm_build_hes(av,msize,queue,ddDelta,Delta,S,H1,hes)
!       type(avItem),intent(inout) :: av
!       integer, intent(in) :: msize
!       type(util_HistoryStore)  :: queue
!       type(Matrix), intent(in) :: ddDelta, Delta, S, H1
!       real(realk), intent(out) :: hes
!       type(Matrix) :: wrk
!       real(Realk) :: do_dc(msize)
!       integer :: ndim
!
!       ndim = S%nrow
!       call mat_init(wrk,ndim,ndim)
!
!       hes = 2E0_realk*mat_dotproduct(ddDelta,Gpara) !2Trd²Delta/dc_i dc_j G(Delta_para)
!!
!! Evaluate d² o_i /d c_x d c_y for d² Delta_par / d c_x d c_y= sum_i [d² o_i /d c_x d c_y] D_i
!!
!       call util_GET_PROJ_PART(queue,ddDelta,S,do_dc)
!!
!! Evaluate   d² Delta_par / d c_x d c_y= sum_i [d² o_i /d c_x d c_y] D_i
!!
!       call get_average_arr(av,'D',queue, dsm_history_size, do_dc, wrk)
!       hes = hes - 2E0_realk*mat_dotproduct(wrk,Gpara) ! - 2Trd²Delta_para/dc_i dc_j G(Delta_para)
!!
!! Evaluate   d² G(Delta_par) / d c_x d c_y= sum_i [d² o_i /d c_x d c_y] (F_i - h)
!!
!       call get_average_arr('F',queue, dsm_history_size, do_dc, wrk)
!       call mat_daxpy(-SUM(do_dc(1:msize)),H1,wrk)
!       hes = hes + 2E0_realk*mat_dotproduct(Delta,wrk)  !  +2TrDelta d²G(Delta_para)/dc_i dc_j 
!       call mat_free(wrk)
!   
!   end subroutine xterm_build_hes

!!> \brief Construct derivatives of extra term.
!   subroutine xterm_build_deriv(msize,queue,Delta,dDelta,S,H1,grad,hes)
!      
!      implicit none
!      integer, intent(in) :: msize
!      type(util_HistoryStore)  :: queue
!      type(Matrix), intent(in) :: Delta, dDelta(msize), S, H1
!      real(realk), intent(inout) :: grad(msize),hes(msize,msize)
!      type(Matrix) :: wrk,dGpara(msize)
!      real(Realk) :: do_dc(msize,msize),grad_ps(msize)
!      integer :: ndim,i,j
!
!      ndim = S%nrow
!      call mat_init(wrk,ndim,ndim)
!!
!! Evaluate the gradient terms
!!
!      do i = 1,msize
!        call mat_init(dGpara(i),ndim,ndim)
!        grad(i) = grad(i) + 2E0_realk*mat_dotproduct(dDelta(i),Gpara)
!        call util_GET_PROJ_PART(queue,dDelta(i),S,do_dc(:,i))
!        call get_average_arr(av,'D',queue, dsm_history_size, do_dc(:,i), wrk)
!        grad(i) = grad(i) - 2E0_realk*mat_dotproduct(wrk,Gpara)
!        call get_average_arr(av,'F',queue, dsm_history_size, do_dc(:,i), dGpara(i))
!        call mat_daxpy(-SUM(do_dc(:,i)),H1,dGpara(i))
!        !grad_ps(i) = mat_Tr(dGpara(i))
!        grad(i) = grad(i) +2E0_realk*mat_dotproduct(Delta,dGpara(i))
!      enddo
!!
!! Evaluate the rest of the Hessian terms
!!
!      do i = 1,msize
!        call get_average_arr(av,'D',queue, dsm_history_size, do_dc(:,i), wrk)
!        do j = 1,msize
!          hes(i,j) = hes(i,j) + 2E0_realk*mat_dotproduct(dDelta(i),dGpara(j)) &
!                   &          + 2E0_realk*mat_dotproduct(dDelta(j),dGpara(i)) &
!                   &          - 2E0_realk*mat_dotproduct(wrk,dGpara(j))
!        enddo
!      enddo
!  
!      do i = 1,msize
!        call mat_free(dGpara(i))
!      enddo
!      call mat_free(wrk)
!   end subroutine xterm_build_deriv

!> \brief Shut down the extra term module.
   subroutine xterm_end()
   implicit none
       call mat_free(Gpara)
       DEALLOCATE(coef)
   end subroutine xterm_end

   !real(realk) function xterm_Eps(ndim,msize,queue,Delta,S,H1)

   !   integer, intent(in) :: ndim,msize
   !   type(Matrix), intent(in) :: Delta,S,H1
   !   type(util_HistoryStore), intent(in) :: queue
   !   type(Matrix) :: Dpara   
 
   !   call mat_init(Dpara,ndim,ndim)
   !   call xterm_build_Gpara(msize,queue,Delta,S,H1)
   !   call get_average_arr('D',queue, dsm_history_size, coef, Dpara)
   !   !xterm_Eps = mat_dotproduct(Dpara,Gpara)
   !   !xterm_Eps = mat_Tr(Dpara)
   !   xterm_Eps = mat_Tr(Gpara)
   !   call mat_free(Dpara)
   !end function xterm_Eps

end module dsm_x_term
