module orbspread_hess_prec_mod
!##########################################################
!#            GENERAL INTERFACE ROUTINES                  #
!# Below are routine that are kept outside modules.       #
!# They are interface routines to solvers (precond/       #
!# linear trans.) and lsdalton main program (optimloc).   #
!#                                                        #
!##########################################################
  use precision
  use matrix_module, only: matrix
  use matrix_operations 
  use matrix_operations_aux 
  use matrix_util!, only: matrix_exponential
  use matrix_util
  use loc_utils
  use typedeftype
  use LSTIMING
  use localitymeasureMod
  use decompMod
  private
  public :: orbspread_hesslin, orbspread_precond
  CONTAINS

!> \brief hessian linear transformation for power m orbital spread locality measure
!> \author B. Jansik
!> \date 2010
!> \param Hv linear transformation of H with trial vector V
!> \param V, trial vector
!> \param mu, level shift (H-muI)
!> \param norb, size of matrices (same as number of orbitals involved in localization)
!> \param orbspread_input structure holding orbspread related data
!!> \param m power
!!> \param spread2 vector of squared orbital spreads 
!!> \param R(3) matrix array holding XDIPLEN,YDIPLEN,ZDIPLEN components
!!> \param Q matrix holding sum XXSECMOM+YYSECMOM+ZZSECMOM
!!> \param tmpM(4) matrix array used as workspace
!> all parameters must be allocated prior to subroutine call
!> all type(matrix) arguments are of norb,norb size
    subroutine orbspread_hesslin(Hv,V,mu,norb,orbspread_input)!m,spread2,R,Q,tmpM)
      implicit none

      Type(Matrix), intent(inout) :: Hv
      Type(Matrix), intent(in)  :: V
      real(realk), intent(in)   :: mu
      integer, intent(in)       :: norb
      !  type(orbspread_data), intent(in), target :: orbspread_input
      type(orbspread_data), target :: orbspread_input

      integer       :: m
      !  Type(Matrix)  :: R(3), RV(3)
      Type(Matrix), pointer  ::  Q, G
      Type(Matrix) ::  QV ,tmpM
      real(realk), pointer   :: spread2(:)
      real(realk)  :: diagQV(norb), diagR(norb,3), diagRV(norb,3), tmp(norb)
      Type(Matrix) :: inptmpMx
      integer      :: x,y,i

      !pointer assignments
      !stupid f90, no support for pointer arrays!
      !  do i=1,3
      !     call mat_clone(R(i),orbspread_input%R(i))
      !     call mat_clone(RV(i),orbspread_input%tmpM(i))
      !  enddo
      Q       => orbspread_input%Q
!      QV      => orbspread_input%tmpM(4)
!      tmpM    => orbspread_input%tmpM(4) !it's ok, QV will not be needed at that point
      G       => orbspread_input%G

      spread2 => orbspread_input%spread2

      m       =  orbspread_input%m


      !job
      call mat_init(QV,Q%nrow,V%ncol)
      call mat_mul(Q,V,'n','n',1E0_realk,0E0_realk,QV)
      call mat_extract_diagonal(diagQV,QV)
      do i=1,norb
         tmp(i) = diagQV(i)*(spread2(i)**(m-2))
      enddo
      call mat_zero(Hv)
      call mat_dmul(tmp,Q,'n',-4E0_realk*m*(m-1),0E0_realk,Hv)
      do i=1,norb
         tmp(i) =  (spread2(i)**(m-1))
      enddo
      call mat_dmul(tmp,QV,'t',-2E0_realk*m,1E0_realk,Hv)
      call mat_dmul(tmp,QV,'n',-2E0_realk*m,1E0_realk,Hv)
      call mat_free(QV)


      call mat_init(inptmpMx,orbspread_input%R(1)%nrow,V%ncol)
      do x=1, 3
         call mat_mul(orbspread_input%R(x),V,'n','n',1E0_realk,0E0_realk,inptmpMx)
         call mat_extract_diagonal(diagRV(:,x),inptmpMx)
         call mat_extract_diagonal(diagR(:,x),orbspread_input%R(x))
      enddo
      do x=1, 3
         call mat_mul(orbspread_input%R(x),V,'n','n',1E0_realk,0E0_realk,inptmpMx)
         do i=1,norb
            tmp(i) = diagR(i,x)*diagQV(i)*(spread2(i)**(m-2))
         enddo
         call mat_dmul(tmp,orbspread_input%R(x),'n',8E0_realk*m*(m-1),1E0_realk,Hv)

         do i=1,norb
            tmp(i) = diagR(i,x)*diagRV(i,x)*(spread2(i)**(m-2))
         enddo
         call mat_dmul(tmp,Q,'n',8E0_realk*m*(m-1),1E0_realk,Hv)

         do y=1, 3
            do i=1,norb
               tmp(i) = diagR(i,x)*diagR(i,y)*diagRV(i,x)*(spread2(i)**(m-2))
            enddo
            call mat_dmul(tmp,orbspread_input%R(y),'n',-16E0_realk*m*(m-1),1E0_realk,Hv)
         enddo

         do i=1,norb
            tmp(i) = diagRV(i,x)*(spread2(i)**(m-1))
         enddo
         call mat_dmul(tmp,orbspread_input%R(x),'n',8E0_realk*m,1E0_realk,Hv)

         do i=1,norb
            tmp(i) = diagR(i,x)*(spread2(i)**(m-1))
         enddo
         call mat_dmul(tmp,inptmpMx,'t',4E0_realk*m,1E0_realk,Hv)
         
         call mat_dmul(tmp,inptmpMx,'n',4E0_realk*m,1E0_realk,Hv)
      enddo
      call mat_free(inptmpMx)

      call mat_mul(V,G,'n','n',0.5E0_realk,1E0_realk,Hv)
      !call mat_mul(V,G,'n','n',1E0_realk,1E0_realk,Hv)

      call mat_init(tmpM,Hv%ncol,Hv%nrow)
      call mat_trans(Hv,tmpM)
      call mat_daxpy(-1E0_realk,tmpM,Hv)
      call mat_free(tmpM)

      !call mat_scal(0.5E0_realk,Hv)
      if (dabs(mu) > 1.0E-8_realk) call mat_daxpy(-mu,V,Hv)

    end subroutine orbspread_hesslin

    subroutine orbspread_precond(Xout,X,mu,inp)
      implicit none
      type(Matrix), intent(inout) :: Xout
      type(Matrix) :: X
      real(realk),  intent(in) :: mu
      type(orbspread_data), intent(in) :: inp
      real(realk),pointer   :: tmp(:), tmpP(:)
      integer                   :: i, ne
      real(realk),pointer :: xvec(:),yvec(:)
      type(matrix)  :: xtemp

      ne = Xout%nrow*Xout%ncol

      call mat_zero(Xout)

      select case(matrix_type)
      case(mtype_dense)


         do i=1,ne
            if (dabs(inp%P%elms(i)- mu) > 1d-8) Xout%elms(i) = X%elms(i)/(inp%P%elms(i) - mu)
         enddo
      case(mtype_scalapack)
         call mat_copy(1d0,X,Xout)
         call mat_hdiv(Xout,inp%P,mu)
      case default

         call mem_alloc(tmp,X%nrow*X%ncol)
         call mat_to_full(X,1E0_realk,tmp)
         call mem_alloc(tmpP,inp%P%nrow*inp%P%ncol)
         call mat_to_full(inp%P,1E0_realk,tmpP)

         do i=1,ne
            if (dabs(tmpP(i) - mu)> 1d-8) tmp(i) = tmp(i)/(tmpP(i) - mu)
         enddo
 
         call mem_dealloc(tmpP)
         
         call mat_set_from_full(tmp,1E0_realk,Xout)

         call mem_dealloc(tmp)

      end select

    end subroutine orbspread_precond

end module orbspread_hess_prec_mod
