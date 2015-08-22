!> @file
!> Contains routines needed for CROP solver in CC driver
!> \author Marcin Ziolkowski

!> Set of subroutines for DIIS or CROP solver for linear or nonlinear set of
!> equations in CC methods.
module crop_tools_module

   use precision
   use tensor_interface_module
   
   
   ! DEC DEPENDENCIES (within deccc directory)   
   ! *****************************************
   use array2_simple_operations
   use array4_simple_operations
   
   ! Interface for CC2 or CCSD energies
   interface get_cc_energy
      module procedure get_cc_energy_arrold
      module procedure get_cc_energy_arrnew
   end interface
   ! Interface for MP2 and CCD energies
   interface get_mp2_energy
      module procedure get_mp2_energy_arrold
      module procedure get_mp2_energy_arrnew
   end interface

   integer,parameter :: SOLVE_AMPLITUDES      = 1
   integer,parameter :: SOLVE_AMPLITUDES_PNO  = 2
   integer,parameter :: SOLVE_MULTIPLIERS     = 3

   private
   
   public :: SolveLinearEquations, CalculateDIIScoefficients,&
        & CalculateDIIScoefficientsII, PrintMatrix,&
        & print_ccjob_header, print_ccjob_iterinfo,&
        & print_ccjob_summary, can_local_trans, local_can_trans,&
        & successive_wxyz_trafo, successive_xyxy_trafo,&
        & SOLVE_AMPLITUDES, SOLVE_AMPLITUDES_PNO,&
        & SOLVE_MULTIPLIERS, get_mp2_energy, get_cc_energy
   
   contains
   
   !> \brief Solve system of linear equations using DGESV lapack routine 
   !> \author Marcin Ziolkowski
   !> \param Amat Square matrix A in Ax=b
   !> \param Bvec RHS vector (b)
   !> \param n Size of the problem (A(n,n))
   subroutine SolveLinearEquations(Amat,Bvec,n)

      implicit none
      real(realk), dimension(n,n), intent(in) :: Amat
      real(realk), dimension(n), intent(inout) :: Bvec
      integer, dimension(n) :: ipiv
      integer :: n,infoLAPACK
      external dgesv

      infoLAPACK = 0

      call dgesv(n,1,Amat,n,ipiv,Bvec,n,infoLAPACK)

      if(infoLAPACK /= 0) then
         print *,"DGESV INFO:",infoLAPACK
         call lsquit('SolveLinearEquations: error in LAPACK DGESV routine',-1)
      end if

   end subroutine SolveLinearEquations

   !> \brief Solve CROP/DIIS system of equations to get mixing coefficients 
   !> \author Marcin Ziolkowski
   !> \param maxDIIS Maximum size of iterative subspace
   !> \param maxIter Maximum number of iterations
   !> \param iter Current iteration number
   !> \param Bmat B matrix in CROP/DIIS
   !> \param c Mixing coefficients
   !> \param verbose Verbose parameter
   subroutine CalculateDIIScoefficients(maxDIIS,maxIter,iter,Bmat,c,verbose)

      implicit none
      integer, intent(in) :: maxDIIS,maxIter,iter
      integer :: i,j,k,l
      real(realk), dimension(maxIter,maxIter), intent(in) :: Bmat
      real(realk), dimension(maxIter), intent(out) :: c
      logical, intent(in) :: verbose

      ! Maximal size of the subspace equations should be 4x4 (3 from vectors and 1
      ! from Lagrange coefficient). 
      real(realk), dimension(min(maxDIIS+1,iter+1),min(maxDIIS+1,iter+1)) :: bb
      real(realk), dimension(min(maxDIIS+1,iter+1)) :: cc

      c=0.0E0_realk

      bb=1.0E0_realk

      k=1
      do i=max(iter-maxDIIS+1,1),iter
         l=1
         do j=max(iter-maxDIIS+1,1),iter
            bb(k,l)=Bmat(i,j)
            l=l+1
         end do 
         k=k+1
      end do

      bb(min(maxDIIS+1,iter+1),min(maxDIIS+1,iter+1))=0E0_realk
      cc=0E0_realk; cc(min(maxDIIS+1,iter+1))=1E0_realk

      if(verbose) call PrintMatrix(bb,min(maxDIIS+1,iter+1),'B matrix  ')

      ! --- solve system of the linear equations ---   
      call SolveLinearEquations(bb,cc,min(maxDIIS+1,iter+1))
      ! ---

      k=1
      c=0E0_realk
      do i=max(iter-maxDIIS+1,1),iter
         c(i)=cc(k)
         k=k+1
      end do

      if(verbose) print *,"Coeffs: ",c(max(iter-maxDIIS+1,1):iter)
   end subroutine CalculateDIIScoefficients
   subroutine CalculateDIIScoefficientsII(nSS,Bmat,c,verbose)

      implicit none
      integer, intent(in) :: nSS
      integer :: i,j,k,l
      real(realk), dimension(nSS,nSS), intent(in) :: Bmat
      real(realk), dimension(nSS), intent(out) :: c
      logical, intent(in) :: verbose

      ! Maximal size of the subspace equations should be 4x4 (3 from vectors and 1
      ! from Lagrange coefficient). 
      real(realk), dimension(nSS+1,nSS+1) :: bb
      real(realk), dimension(nSS+1) :: cc

      c=0.0E0_realk

      cc=0E0_realk
      cc(nSS+1)=1E0_realk

      !build Bcrop in blocks
      bb(1:nSS,nSS+1) = 1.0E0_realk
      bb(nSS+1,1:nSS) = 1.0E0_realk
      bb(1:nSS,1:nSS) = Bmat
      bb(nSS+1,nSS+1) = 0E0_realk

      if(verbose) call PrintMatrix(bb,nSS+1,'B matrix  ')

      ! --- solve system of the linear equations ---   
      call SolveLinearEquations(bb,cc,nSS+1)
      ! ---

      k=1
      c = cc(1:nSS)

      if(verbose) print *,"Coeffs: ",c(1:nSS)

   end subroutine CalculateDIIScoefficientsII

   !> \brief Print simple fortran matrix in nice form, just for debuging purpose 
   !> \author Marcin Ziolkowski
   !> \param my_matrix Matrix to print
   !> \param m Size of the matrix
   !> \param title Title for the matrix
   subroutine PrintMatrix(my_matrix,m,title)

      implicit none
      real(realk), dimension(m,m) :: my_matrix
      integer :: m,i,j
      character(len=10), intent(in), optional :: title

      print *,''
      if(present(title)) print *,' *** ',title,' ***'
      do i=1,m
         print 100,(my_matrix(i,j),j=1,m)
      enddo
      print *,''

      100 format(15f10.5) 

      return
   end subroutine PrintMatrix

   !> \brief Standard mp2 correlation energy
   !> \author Kasper Kristensen
   !> \return Full molecular MP2 energy
   !> \param t2 Double amplitudes
   !> \param Lmo Two-electron integrals L_{bjai} = 2*g_{bjai} - g_{ajbi}
   function get_mp2_energy_arrold(t2,Lmo) result(ecorr)

      implicit none
      type(array4), intent(in) :: Lmo,t2
      real(realk) :: ecorr

      ! Ecorr = sum_{aibj} t2_{bjai}*Lmo_{bjai}
      Ecorr=t2*Lmo

   end function get_mp2_energy_arrold

   !> \brief Coupled-cluster correlation energy
   !> \author Marcin Ziolkowski
   !> \return Full molecular CC correlation energy
   !> \param t2 Single amplitudes
   !> \param t2 Double amplitudes
   !> \param gmo Two-electron integrals (ia|jb)
   !> \param nocc Number of occupied orbitals
   !> \param nvirt Number of unoccupied orbitals
   function get_cc_energy_arrold(t1,t2,gmo,nocc,nvirt) result(ecorr)

      implicit none
      type(array2), intent(in) :: t1
      type(array4), intent(in) :: gmo,t2
      integer, intent(in) :: nocc,nvirt
      real(realk) :: ecorr,ecorr_s,ecorr_d
      integer :: a,i,b,j

      ecorr = 0.0E0_realk
      ecorr_s = 0.0E0_realk
      ecorr_d = 0.0E0_realk

      do j=1,nocc
         do b=1,nvirt
            do i=1,nocc
               do a=1,nvirt
                  ecorr_d = ecorr_d + t2%val(a,i,b,j)* &
                     (2.0E0_realk*gmo%val(i,a,j,b)-gmo%val(i,b,j,a))
                  ecorr_s = ecorr_s + ( t1%val(a,i)*t1%val(b,j) ) * &
                     (2.0E0_realk*gmo%val(i,a,j,b)-gmo%val(i,b,j,a))
               end do
            end do
         end do
      end do

      if(DECinfo%cc_driver_debug) then
         print *,' Singles energy : ',ecorr_s
         print *,' Doubles energy : ',ecorr_d
      end if

      ecorr = ecorr_s + ecorr_d

      return
   end function get_cc_energy_arrold

   function get_cc_energy_arrnew(t1,t2,gmo,nocc,nvirt) result(ecorr)

      implicit none
      type(tensor), intent(inout) :: t1
      type(tensor), intent(in) :: t2
      type(tensor), intent(inout) :: gmo
      integer, intent(in) :: nocc,nvirt
      real(realk) :: ecorr,ecorr_s,ecorr_d
      integer :: a,i,b,j

      ecorr = 0.0E0_realk
      ecorr_s = 0.0E0_realk
      ecorr_d = 0.0E0_realk

      if(t2%itype==TT_DENSE.and.gmo%itype==TT_DENSE.and.(t1%itype==TT_DENSE.or.t1%itype==TT_REPLICATED))then

         do j=1,nocc
            do b=1,nvirt
               do i=1,nocc
                  do a=1,nvirt
                     ecorr_d = ecorr_d + t2%elm4(a,b,i,j)* &
                        (2.0E0_realk*gmo%elm4(i,a,j,b)-gmo%elm4(i,b,j,a))
                     ecorr_s = ecorr_s + ( t1%elm2(a,i)*t1%elm2(b,j) ) * &
                        (2.0E0_realk*gmo%elm4(i,a,j,b)-gmo%elm4(i,b,j,a))
                  end do
               end do
            end do
         end do

         if(DECinfo%cc_driver_debug) then
            print *,' Singles energy : ',ecorr_s
            print *,' Doubles energy : ',ecorr_d
         end if

         ecorr = ecorr_s + ecorr_d

      else if(t2%itype==TT_TILED_DIST.and.gmo%itype==TT_TILED_DIST)then

         t1%itype = TT_REPLICATED
         call tensor_sync_replicated(t1)
         ecorr    = get_cc_energy_parallel(t2,gmo,t1)
         t1%itype = TT_DENSE

      endif


   end function get_cc_energy_arrnew

   function get_mp2_energy_arrnew(t2,gmo,nocc,nvirt) result(ecorr)

      implicit none
      type(tensor), intent(in) :: t2
      type(tensor), intent(inout) :: gmo
      integer, intent(in) :: nocc,nvirt
      real(realk) :: ecorr,ecorr_d
      integer :: a,i,b,j

      ecorr = 0.0E0_realk
      ecorr_d = 0.0E0_realk

      if(t2%itype==TT_DENSE.and.gmo%itype==TT_DENSE)then
         !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) PRIVATE(i,a,j,b) SHARED(nocc,&
         !$OMP nvirt,t2,gmo) REDUCTION(+:ecorr_d)
         do j=1,nocc
            do b=1,nvirt
               do i=1,nocc
                  do a=1,nvirt
                     ecorr_d = ecorr_d + t2%elm4(a,b,i,j)* &
                        (2.0E0_realk*gmo%elm4(i,a,j,b)-gmo%elm4(i,b,j,a))
                  end do
               end do
            end do
         end do
         !$OMP END PARALLEL DO

         ecorr = ecorr_d

      else if(t2%itype==TT_TILED_DIST.and.gmo%itype==TT_TILED_DIST)then

         ecorr=get_cc_energy_parallel(t2,gmo)

      endif


   end function get_mp2_energy_arrnew



   
   !> \brief Get antisymmetrized double amplitudes
   !> \return Array4 structure with antisymmetrized double amplitudes
   function get_u(t2) result(u)

      implicit none
      type(array4), intent(inout) :: t2
      type(array4) :: u
      integer, dimension(4) :: dims
      integer :: a,i,b,j

      dims = t2%dims

#ifdef EXTRA_SIMPLE

      u = array4_init(dims)

      do j=1,dims(4)
         do b=1,dims(3)
            do i=1,dims(2)
               do a=1,dims(1)

                  u%val(a,i,b,j) = 2.0E0_realk*t2%val(a,i,b,j) - t2%val(a,j,b,i)

               end do
            end do
         end do
      end do

#else

      u = array4_duplicate(t2)
      call array4_scale(u,2.0E0_realk)
      call array4_reorder(t2,[1,4,3,2])
      call array4_add_to(u,-1.0E0_realk,t2)
      call array4_reorder(t2,[1,4,3,2])

#endif

   end function get_u
   
   !> \brief Print header and info about coupled-cluster job
   !> \author Marcin Ziolkowski
   !> \param ccPrintLevel Print level
   !> \param framgment_job Fragment job
   !> \param nbasis Number of basis functions
   !> \param nocc Number of occupied orbitals
   !> \param nvirt Number of unoccupied orbitals
   subroutine print_ccjob_header(ccmodel,ccPrintLevel,fj,multiplier_job,&
         &nbasis,nocc,nvirt,maxsub,restart,converged,old_iter)
      implicit none
      integer, intent(in) :: ccmodel,ccPrintLevel,nbasis,nocc,nvirt,maxsub,old_iter
      logical, intent(in) :: fj,multiplier_job,restart,converged

      if(ccPrintLevel > 0) then
         if(.not.fj) then
            if(multiplier_job)then
               write(DECinfo%output,'(/,a)') '-------------------------------'
               write(DECinfo%output,'(a)')   '  Coupled-cluster multipliers  '
               write(DECinfo%output,'(a,/)') '-------------------------------'
            else
               write(DECinfo%output,'(/,a)') '--------------------------'
               write(DECinfo%output,'(a)')   '  Coupled-cluster energy  '
               write(DECinfo%output,'(a,/)') '--------------------------'
            endif
            if(DECinfo%CCDhack)then
               write(DECinfo%output,'(a,a)')      'Wave function    = ','CCD'
            else
               write(DECinfo%output,'(a,a)')      'Wave function    = ',DECinfo%cc_models(ccModel)
            endif
            write(DECinfo%output,'(a,i4)')     'MaxIter          = ',DECinfo%ccMaxIter
            write(DECinfo%output,'(a,i4)')     'Num. b.f.        = ',nbasis
            write(DECinfo%output,'(a,i4)')     'Num. occ. orb.   = ',nocc
            write(DECinfo%output,'(a,i4)')     'Num. unocc. orb. = ',nvirt
            write(DECinfo%output,'(a,e8.1e2)') 'Convergence      = ',DECinfo%ccConvergenceThreshold
            write(DECinfo%output,'(a,l4)')     'Debug mode       = ',DECinfo%cc_driver_debug
            write(DECinfo%output,'(a,i4)')     'Print level      = ',ccPrintLevel
            write(DECinfo%output,'(a,l4)')     'Use CROP         = ',DECinfo%use_crop
            write(DECinfo%output,'(a,i4)')     'CROP subspace    = ',maxsub
            write(DECinfo%output,'(a,l4)')     'Preconditioner   = ',DECinfo%use_preconditioner
            write(DECinfo%output,'(a,l4)')     'Precond. B       = ',DECinfo%use_preconditioner_in_b
            write(DECinfo%output,'(a,l4)')     'Singles          = ',DECinfo%use_singles
         else
            write(DECinfo%output,'(/,a)') '  Coupled-cluster energy  -> Fragment job '
            write(DECinfo%output,'(a)')   '------------------------------------------'
            if(DECinfo%CCDhack)then
               write(DECinfo%output,'(a,a)')      'Wave function    = ','CCD'
            else
               write(DECinfo%output,'(a,a)')      'Wave function    = ',DECinfo%cc_models(ccModel)
            endif
            write(DECinfo%output,'(4x,a,l4)')     'Debug mode       = ',DECinfo%cc_driver_debug
            write(DECinfo%output,'(a,i4,$)')      'MaxIter          = ',DECinfo%ccMaxIter
            write(DECinfo%output,'(5x,a,e8.1e2)') 'Convergence      = ',DECinfo%ccConvergenceThreshold
            write(DECinfo%output,'(a,i4,$)')      'Num. b.f.        = ',nbasis
            write(DECinfo%output,'(5x,a,i4)')     'Print level      = ',ccPrintLevel
            write(DECinfo%output,'(a,i4,$)')      'Num. occ. orb.   = ',nocc
            write(DECinfo%output,'(5x,a,i4)')     'CROP subspace    = ',DECinfo%ccMaxDIIS
            write(DECinfo%output,'(a,i4,$)')      'Num. unocc. orb. = ',nvirt
            write(DECinfo%output,'(5x,a,l4)')     'Preconditioner   = ',DECinfo%use_preconditioner
         end if

         ! cc parameters
         !if(ccPrintLevel > 0) then
         !   if(fj) then
         !      write(DECinfo%output,'(/,a,a)') &
         !           '----  -------------   -------------   -------------   -------------   -------------     ------'
         !      write(DECinfo%output,'(a,a)') &
         !           'Iter   1norm(S)        1norm(D)        2norm(S+D)      Targ-N(S+D)     Targ-energy       time  '
         !      write(DECinfo%output,'(a,a)') &
         !           '----  -------------   -------------   -------------   -------------   -------------     ------'
         !   else
         !      write(DECinfo%output,'(/,a,a)') &
         !           '----  -------------   -------------   -------------   -------------   -------------     ------'
         !      write(DECinfo%output,'(a,a)') &
         !           'Iter   1norm(S)        1norm(D)        1norm(S+D)      2norm(S+D)      energy            time  '
         !      write(DECinfo%output,'(a,a)') &
         !           '----  -------------   -------------   -------------   -------------   -------------     ------'
         !   end if
         !end if
      end if

      if(.not.fj.or.DECinfo%PL>1)then
         if(.not.multiplier_job)then
            if(restart)then
               if( converged )then
                  print *
                  print '(X,a,I4)', '### Skipping CC calculation with &
                     &converged amplitudes of iteration',old_iter
                  print '(X,a)', '### --------------------------------&
                     &---------------------------------'

                  write(DECinfo%output,*)
                  write(DECinfo%output,'(X,a,I4)') '### Skipping CC calculation with &
                     &converged amplitudes of iteration',old_iter
                  write(DECinfo%output,*) '### --------------------------------&
                     &---------------------------------'
               else

                  print *
                  print '(X,a,I4)', '### Re-starting CC iterations &
                     &from amplitudes of iteration',old_iter
                  print '(X,a)', '### -----------------------------&
                     &-----------------------------'

                  write(DECinfo%output,*)
                  write(DECinfo%output,'(X,a,I4)') '### Re-starting CC iterations from amplitudes of iteration',old_iter
                  write(DECinfo%output,*) '### ----------------------------------------------------------'
               endif
            else
               print *
               print *, '### Starting CC iterations'
               print *, '### ----------------------'

               write(DECinfo%output,*)
               write(DECinfo%output,*) '### Starting CC iterations'
               write(DECinfo%output,*) '### ----------------------'
            endif
            print '(1X,a)',  '###  Iteration     Residual norm          CC energy'
            write(DECinfo%output,'(1X,a)')  '###  Iteration     Residual norm          CC energy'
         else
            print *
            print *, '### Starting  Lagrangian iterations'
            print *, '### -------------------------------'
            print '(1X,a)',  '###  Iteration     Residual norm'

            write(DECinfo%output,*)
            write(DECinfo%output,*) '### Starting Lagrangian iterations'
            write(DECinfo%output,*) '### ------------------------------'
            write(DECinfo%output,'(1X,a)')  '###  Iteration     Residual norm'
         endif
      endif

   end subroutine print_ccjob_header


   !it = iteration number
   !norm = norm of the total residual
   !ce = correlation energy
   !gm = multipliers yes or no?
   !fj = fragment job?
   subroutine print_ccjob_iterinfo(it,norm,ce,gm,fj)
      implicit none
      integer,intent(in)     :: it
      real(realk),intent(in) :: norm,ce
      logical, intent(in)    :: gm,fj
      if( .not. fj .or. DECinfo%PL>1 )then
         if( gm ) then
            print '(1X,a,2X,i4,5X,g19.9,4X)',  '### ',it, norm
            write(DECinfo%output,'(1X,a,2X,i4,5X,g19.9,4X)') &
               &   '### ',it, norm
         else
            print '(1X,a,2X,i4,5X,g19.9,4X,g19.9)',  '### ',it, norm,ce
            write(DECinfo%output,'(1X,a,2X,i4,5X,g19.9,4X,g19.9)') &
               &   '### ',it, norm,ce
         endif
      endif

   end subroutine print_ccjob_iterinfo



   ! bi is the logical for break_iterations => success of the crop procedure
   ! gm is the logical for get_mult(ipliers) => ccequations vs left-hand
   ! transformations
   ! fj is the logical for fragment_job => fj or not fj, that is the question
   ! li is the number of the last iteration
   ! ce is the correlation energy
   ! t* are timings
   ! nt*, nm* are the squared l_2 norms of the amplitudes and multipliers
   subroutine print_ccjob_summary(bi,gm,fj,li,us,ce,tew,tsw,tec,tsc,nt1,nt2,nm1,nm2)
      implicit none
      logical, intent(in)        :: bi,gm,fj,us
      integer, intent(in)        :: li
      real(realk), intent(in)    :: ce,tew,tsw,tec,tsc
      real(realk),intent(in) :: nt1
      real(realk),intent(in) :: nt2
      real(realk),intent(in),optional :: nm1
      real(realk),intent(in),optional :: nm2
      real(realk) :: snorm,dnorm,tnorm
      tnorm = 0.0E0_realk
      dnorm = 0.0E0_realk
      snorm = 0.0E0_realk

      if(gm)then
         if(us)then
            if(.not. present(nm1) ) call lsquit('ERROR(print_ccjob_summary) no singles multipliers norm present',-1)
            snorm = nm1
         endif
         if(.not. present(nm2) ) call lsquit('ERROR(print_ccjob_summary) no doubles multipliers norm present',-1)
         dnorm = nm2
      else
         if(us)snorm = nt1
         dnorm = nt2
      endif

      tnorm = sqrt(snorm+dnorm)
      if(us)snorm = sqrt(snorm)
      dnorm = sqrt(dnorm)


      if( .not. fj .or. DECinfo%PL>1)then
         write(DECinfo%output,*)
         write(DECinfo%output,'(/,a)') '-------------------------------'
         write(DECinfo%output,'(a)')   '  Coupled-cluster job summary  '
         write(DECinfo%output,'(a,/)') '-------------------------------'
         if(bi.and.li<DECinfo%ccMaxIter) then
            if(gm)then
               write(DECinfo%output,'(a)')     'Yeeehaw! left-transformations converged!'
            else
               write(DECinfo%output,'(a)')     'Hooray! CC equation is solved!'
            endif
         else
            write(DECinfo%output,'(a,i4,a)')  'Equation not solved in ', &
               & DECinfo%ccMaxIter, ' iterations!'
            call lsquit('CC equation not solved!',DECinfo%output)
         end if
         write(DECinfo%output,'(a,g10.3,a)') 'CCSOL: Total cpu time    = ',tec-tsc,' s'
         write(DECinfo%output,'(a,g10.3,a)') 'CCSOL: Total wall time   = ',tew-tsw,' s'

         !if(fj) then
         !   write(DECinfo%output,'(a,f16.10)')  'Frag. corr. energy = ',ce
         !else
            if(gm)then
               if(us)then
                  write(DECinfo%output,'(a,g14.7)')  'Singles multiplier norm  = ',snorm
               endif
               write(DECinfo%output,'(a,g14.7)')  'Doubles multiplier norm  = ',dnorm
               write(DECinfo%output,'(a,g14.7)')  'Total multiplier norm    = ',tnorm
            else
               if(us)then
                  write(DECinfo%output,'(a,g14.7)')  'Singles amplitudes norm  = ',snorm
               endif
               write(DECinfo%output,'(a,g14.7)')  'Doubles amplitudes norm  = ',dnorm
               write(DECinfo%output,'(a,g14.7)')  'Total amplitudes norm    = ',tnorm
               write(DECinfo%output,'(a,f17.10)')  'Corr. energy             = ',ce
            endif
         !end if
         write(DECinfo%output,'(a,i5)') 'Number of CC iterations  =', li
      else
         print '(1X,a,1X,i4,1X,a,1X,g19.9)', 'Fragment CC converged after',li,'iterations with norm',tnorm
         if(li==DECinfo%ccMaxIter)then
            call lsquit("ERROR(print_ccjob_summary): ccequations not converged in&
               & current fragment, the dec calculation becomes meaningless",-1)
         endif
      endif
   end subroutine print_ccjob_summary

   !> \brief: transform ccsd_doubles, ccsdpt_singles and ccsdpt_doubles from canonical to local basis
   !> \author: Patrick Ettenhuber adapted from Janus Juul Eriksen
   !> \date: April 2013
   !> \param: t2, gvovo, t1, no and nv are nocc and nvirt, respectively, 
   !<         and U_occ and U_virt are unitary matrices from canonical --> local basis
   subroutine can_local_trans(no,nv,nb,Uocc,Uvirt,vovo,vvoo,oovv,vo,ov,bo,bv,vv,oo)

      implicit none
      !> integers
      integer, intent(in) :: no, nv, nb
      !> general quantities with the respective sizes
      real(realk), intent(inout),optional :: vovo(nv*nv*no*no),vvoo(nv*nv*no*no),oovv(no*no*nv*nv)
      real(realk), intent(inout),optional :: bo(nb*no), bv(nb*nv)
      real(realk), intent(inout),optional :: vo(nv*no),ov(no*nv)
      real(realk), intent(inout),optional :: oo(no*no),vv(nv*nv)
      !> unitary transformation matrices - indices: (local,pseudo-canonical)
      real(realk), intent(in) :: Uocc(no*no), Uvirt(nv*nv)
      !> temp array2 and array4 structures
      real(realk),pointer :: tmp(:)
      integer(kind=8) :: wrksize
      logical :: bg

      wrksize = 0
      if(present(vovo)) wrksize = max(wrksize,(i8*nv**2)*no**2)
      if(present(vvoo)) wrksize = max(wrksize,(i8*nv**2)*no**2)
      if(present(oovv)) wrksize = max(wrksize,(i8*no**2)*nv**2)
      if(present(vo))   wrksize = max(wrksize,(i8*nv) * no)
      if(present(ov))   wrksize = max(wrksize,(i8*no) * nv)
      if(present(bo))   wrksize = max(wrksize,(i8*nb) * no)
      if(present(bv))   wrksize = max(wrksize,(i8*nb) * nv)
      if(present(oo))   wrksize = max(wrksize,(i8*no) * no)
      if(present(vv))   wrksize = max(wrksize,(i8*nv) * nv)

      bg = (mem_is_background_buf_init().and.mem_get_bg_buf_free()>=wrksize)

      if(bg)then
         call mem_pseudo_alloc(tmp,i8*wrksize)
      else
         call mem_alloc(tmp,wrksize)
      endif

      !successive transformation of vovo:
      if(present(vovo)) call successive_wxyz_trafo(nv,no,nv,no,vovo,Uvirt,Uocc,Uvirt,Uocc,tmp)
      !successive transformation of vvoo:
      if(present(vvoo)) call successive_wxyz_trafo(nv,nv,no,no,vvoo,Uvirt,Uvirt,Uocc,Uocc,tmp)
      !successive transformation of oovv:
      if(present(oovv)) call successive_wxyz_trafo(no,no,nv,nv,oovv,Uocc,Uocc,Uvirt,Uvirt,tmp)

      !if t1 trafo has to be done as well
      if(present(vo))then
         !U(a,A) t(AI)    -> t(aI)
         call dgemm('n','n',nv,no,nv,1.0E0_realk,Uvirt,nv,vo,nv,0.0E0_realk,tmp,nv)
         ! tmp(aI) U(i,I)^T   -> t(ai)
         call dgemm('n','t',nv,no,no,1.0E0_realk,tmp,nv,Uocc,no,0.0E0_realk,vo,nv)
      endif

      if(present(ov))then
         !U(i,I) t(IA)    -> t(iA)
         call dgemm('n','n',no,nv,no,1.0E0_realk,Uocc,no,ov,no,0.0E0_realk,tmp,no)
         ! tmp(iA) U(a,A)^T   -> t(ia)
         call dgemm('n','t',no,nv,nv,1.0E0_realk,tmp,no,Uvirt,nv,0.0E0_realk,ov,no)
      endif

      if(present(bo))then
         tmp(1:nb*no) = bo
         ! tmp(alpha,I) U(i,I)^T   -> Co(alpha,i)
         call dgemm('n','t',nb,no,no,1.0E0_realk,tmp,nb,Uocc,no,0.0E0_realk,bo,nb)
      endif

      if(present(bv))then
         tmp(1:nb*nv) = bv
         ! tmp(alpha,A) U(a,A)^T   -> Cv(alpha,a)
         call dgemm('n','t',nb,nv,nv,1.0E0_realk,tmp,nb,Uvirt,nv,0.0E0_realk,bv,nb)
      endif

      if(present(vv))then
         !U(a,A) t(AB)    -> t(aB)
         call dgemm('n','n',nv,nv,nv,1.0E0_realk,Uvirt,nv,vv,nv,0.0E0_realk,tmp,nv)
         ! tmp(aB) U(b,B)^T   -> t(ab)
         call dgemm('n','t',nv,nv,nv,1.0E0_realk,tmp,nv,Uvirt,nv,0.0E0_realk,vv,nv)
      endif

      if(present(oo))then
         !U(i,I) t(IJ)    -> t(iJ)
         call dgemm('n','n',no,no,no,1.0E0_realk,Uocc,no,oo,no,0.0E0_realk,tmp,no)
         ! tmp(iJ) U(j,J)^T   -> t(iJ)
         call dgemm('n','t',no,no,no,1.0E0_realk,tmp,no,Uocc,no,0.0E0_realk,oo,no)
      endif

      if(bg)then
         call mem_pseudo_dealloc(tmp)
      else
         call mem_dealloc(tmp)
      endif

   end subroutine can_local_trans

   subroutine local_can_trans(no,nv,nb,Uocc,Uvirt,vovo,vvoo,oovv,vo,ov,bo,bv,oo,vv)

      implicit none
      !> integers
      integer, intent(in) :: no, nv, nb
      !> general quantities with the respective sizes
      real(realk), intent(inout),optional :: vovo(nv*nv*no*no),vvoo(nv*nv*no*no),oovv(no*no*nv*nv)
      real(realk), intent(inout),optional :: bo(nb*no), bv(nb*nv)
      real(realk), intent(inout),optional :: vo(nv*no), ov(no*nv)
      real(realk), intent(inout),optional :: oo(no*no), vv(nv*nv)
      !> unitary transformation matrices
      !> unitary transformation matrices - indices: (local,pseudo-canonical)
      real(realk), intent(in) :: Uocc(no*no), Uvirt(nv*nv)
      real(realk) :: UoccT(no*no), UvirtT(nv*nv)
      !> temp array2 and array4 structures
      real(realk),pointer :: tmp(:)
      integer :: wrksize
      logical :: bg


      call mat_transpose(no,no,1.0E0_realk,Uocc,0.0E0_realk,UoccT)
      call mat_transpose(nv,nv,1.0E0_realk,Uvirt,0.0E0_realk,UvirtT)

      wrksize = 0
      if(present(vovo)) wrksize = max(wrksize,(i8*nv**2)*no**2)
      if(present(vvoo)) wrksize = max(wrksize,(i8*nv**2)*no**2)
      if(present(oovv)) wrksize = max(wrksize,(i8*no**2)*nv**2)
      if(present(vo))   wrksize = max(wrksize,(i8*nv) * no)
      if(present(ov))   wrksize = max(wrksize,(i8*no) * nv)
      if(present(bo))   wrksize = max(wrksize,(i8*nb) * no)
      if(present(bv))   wrksize = max(wrksize,(i8*nb) * nv)
      if(present(oo))   wrksize = max(wrksize,(i8*no) * no)
      if(present(vv))   wrksize = max(wrksize,(i8*nv) * nv)

      bg = (mem_is_background_buf_init().and.mem_get_bg_buf_free()>=wrksize)

      if(bg)then
         call mem_pseudo_alloc(tmp,i8*wrksize)
      else
         call mem_alloc(tmp,wrksize)
      endif

      !successive transformation of vovo:
      if(present(vovo))call successive_wxyz_trafo(nv,no,nv,no,vovo,UvirtT,UoccT,UvirtT,UoccT,tmp)
      !successive transformation of vvoo:
      if(present(vvoo))call successive_wxyz_trafo(nv,nv,no,no,vvoo,UvirtT,UvirtT,UoccT,UoccT,tmp)
      !successive transformation of oovv:
      if(present(oovv))call successive_wxyz_trafo(no,no,nv,nv,oovv,UoccT,UoccT,UvirtT,UvirtT,tmp)

      if(present(vo))then
         !U(a,A) t(AI)    -> t(aI)
         call dgemm('n','n',nv,no,nv,1.0E0_realk,UvirtT,nv,vo,nv,0.0E0_realk,tmp,nv)
         ! tmp(aI) U(i,I)^T   -> t(ai)
         call dgemm('n','n',nv,no,no,1.0E0_realk,tmp,nv,Uocc,no,0.0E0_realk,vo,nv)
      endif

      if(present(ov))then
         !U(i,I) t(IA)    -> t(iA)
         call dgemm('n','n',no,nv,no,1.0E0_realk,UoccT,no,ov,no,0.0E0_realk,tmp,no)
         ! tmp(iA) U(a,A)^T   -> t(ia)
         call dgemm('n','n',no,nv,nv,1.0E0_realk,tmp,no,Uvirt,nv,0.0E0_realk,ov,no)
      endif

      if(present(bo))then
         tmp(1:nb*no) = bo
         ! tmp(alpha,i) U(i,I)   -> Co(alpha,I)
         call dgemm('n','n',nb,no,no,1.0E0_realk,tmp,nb,Uocc,no,0.0E0_realk,bo,nb)
      endif

      if(present(bv))then
         tmp(1:nb*nv) = bv
         ! tmp(alpha,a) U(a,A)   -> Cv(alpha,A)
         call dgemm('n','n',nb,nv,nv,1.0E0_realk,tmp,nb,Uvirt,nv,0.0E0_realk,bv,nb)
      endif

      if(present(vv))then
         !U(a,A) t(AB)    -> t(aB)
         call dgemm('n','n',nv,nv,nv,1.0E0_realk,UvirtT,nv,vv,nv,0.0E0_realk,tmp,nv)
         ! tmp(aB) U(b,B)^T   -> t(ab)
         call dgemm('n','n',nv,nv,nv,1.0E0_realk,tmp,nv,Uvirt,nv,0.0E0_realk,vv,nv)
      endif

      if(present(oo))then
         !U(i,I) t(IJ)    -> t(iJ)
         call dgemm('n','n',no,no,no,1.0E0_realk,UoccT,no,oo,no,0.0E0_realk,tmp,no)
         ! tmp(iJ) U(j,J)^T   -> t(ij)
         call dgemm('n','n',no,no,no,1.0E0_realk,tmp,no,Uocc,no,0.0E0_realk,oo,no)
      endif

      if(bg)then
         call mem_pseudo_dealloc(tmp)
      else
         call mem_dealloc(tmp)
      endif
   end subroutine local_can_trans

   subroutine successive_wxyz_trafo(w,x,y,z,WXYZ,WW,XX,YY,ZZ,WRKWXYZ)
      implicit none
      integer, intent(in) :: w,x,y,z
      real(realk), intent(inout) :: WXYZ(w*x*y*z),WRKWXYZ(w*x*y*z)
      real(realk), intent(in) :: WW(w,W),XX(x,X),YY(y,Y),ZZ(z,Z)
      !WXYZ(W,XYZ)^T WW(w,W)^T   -> WRKWXYZ (XYZ,w)
      call dgemm('t','t',X*Y*Z,w,W,1.0E0_realk,WXYZ,W,WW,w,0.0E0_realk,WRKWXYZ,X*Y*Z)
      ! WRKWXYZ(X,YZw)^T XX(x,X)^T   -> WXYZ (YZw,x)
      call dgemm('t','t',Y*Z*w,x,X,1.0E0_realk,WRKWXYZ,X,XX,x,0.0E0_realk,WXYZ,Y*Z*w)
      ! WXYZ(Y,Zwx)^T YY(y,Y)^T   -> WRKWXYZ (Zwx,y)
      call dgemm('t','t',Z*w*x,y,Y,1.0E0_realk,WXYZ,Y,YY,y,0.0E0_realk,WRKWXYZ,Z*w*x)
      ! WRKWXYZ(Z,wxy)^T ZZ(z,Z)^T   -> WXYZ (wxyz)
      call dgemm('t','t',w*x*y,z,Z,1.0E0_realk,WRKWXYZ,z,ZZ,Z,0.0E0_realk,WXYZ,w*x*y)
   end subroutine successive_wxyz_trafo

   subroutine successive_xyxy_trafo(x,y,XYXY,XX,YY,WRKYXYX)
      implicit none
      integer, intent(in) :: x,y
      real(realk), intent(inout) :: XYXY(x*y*x*y),WRKYXYX(y*x*y*x)
      real(realk), intent(in) :: XX(x,x),YY(y,y)
      !XYXY(X,YXY)^T XX(x,X)^T   -> WRKYXYX (YXY,x)
      call dgemm('t','t',x*y*y,x,x,1.0E0_realk,XYXY,x,XX,x,0.0E0_realk,WRKYXYX,x*y*y)
      ! WRKYXYX(Y,XYx)^T YY(y,Y)^T   -> XYXY (XYx,y)
      call dgemm('t','t',x*x*y,y,y,1.0E0_realk,WRKYXYX,y,YY,y,0.0E0_realk,XYXY,x*x*y)
      ! XYXY(X,Yxy)^T XX(x,X)^T   -> WRKYXYX (Yxy,x)
      call dgemm('t','t',x*y*y,x,x,1.0E0_realk,XYXY,x,XX,x,0.0E0_realk,WRKYXYX,x*y*y)
      ! WRKYXYX(Y,xyx)^T YY(y,Y)^T   -> XYXY (xyxy)
      call dgemm('t','t',x*x*y,y,y,1.0E0_realk,WRKYXYX,y,YY,y,0.0E0_realk,XYXY,x*x*y)
   end subroutine successive_xyxy_trafo

end module crop_tools_module
