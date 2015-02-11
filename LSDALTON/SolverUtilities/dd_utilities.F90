!> @file
!> Contains Direct Density Optimization utilities module.

!> \brief Mainly routines for getting HOMO-LUMO gap and Hessian when using the parameterization for direct density optimization, D = D(X). 
!> \author S. Host
!> \date 2005
module direct_dens_util
use arhDensity
use decompMod
use files
use lstiming
use queue_ops
use queue_module
use typedeftype
use precision
use matrix_module
use matrix_operations
use matrix_operations_aux

!> \brief Contains stuff concerning HOMO-LUMO gap and Hessian eigenvalues
!> \author S. Host
!> \date 2005
type DDitem
    !> Eigenvectors for occupied orbitals
    real(realk),pointer   :: FUP_eivecs(:,:)
    !> Eigenvectors for virtual orbitals
    real(realk),pointer   :: FUQ_eivecs(:,:)
    !> Occupied orbital energies
    real(realk),pointer   :: FUP_eival(:)
    !> Virtual orbital energies
    real(realk),pointer   :: FUQ_eival(:)

    !> Eigenvectors for occupied alpha orbitals (if unrestricted) 
    real(realk),pointer   :: FUPalfa_eivecs(:,:)
    !> Eigenvectors for virtual alpha orbitals (if unrestricted) 
    real(realk),pointer   :: FUQalfa_eivecs(:,:)
    !> Occupied alpha orbital energies (if unrestricted)
    real(realk),pointer   :: FUPalfa_eival(:)
    !> Virtual alpha orbital energies (if unrestricted)
    real(realk),pointer   :: FUQalfa_eival(:)

    !> Eigenvectors for occupied beta orbitals (if unrestricted) 
    real(realk),pointer   :: FUPbeta_eivecs(:,:)
    !> Eigenvectors for virtual beta orbitals (if unrestricted) 
    real(realk),pointer   :: FUQbeta_eivecs(:,:)
    !> Occupied beta orbital energies (if unrestricted)
    real(realk),pointer   :: FUPbeta_eival(:)
    !> Virtual beta orbital energies (if unrestricted)
    real(realk),pointer   :: FUQbeta_eival(:)

    !> How many occupied orbital energies should be determined? 
    integer                  :: nocc
    !> How many virtual orbital energies should be determined?
    integer                  :: nvirt

    !> How many alpha occupied orbital energies should be determined? (if unrestricted) 
    integer                  :: nocca
    !> How many beta occupied orbital energies should be determined? (if unrestricted) 
    integer                  :: noccb
    !> How many alpha virtual orbital energies should be determined? (if unrestricted) 
    integer                  :: nvirta
    !> How many beta virtual orbital energies should be determined? (if unrestricted) 
    integer                  :: nvirtb

    !> Reduced Hessian for eigenvalue problem
    real(realk), pointer :: Ared(:,:)
    !> Overlap matrix for trial vectors for eigenvalue problem
    real(realk), pointer :: boverlaps(:,:)
    !> Logical unit number for file containing sigma vectors, eigenvalue problem
    integer                  :: lusigma
    !> Logical unit number for file containing trial vectors, eigenvalue problem
    integer                  :: lub
    !> Threshold for dd_Fock_solver
    real(realk)              :: fock_thresh 
    !> True if lowest ARH Hessian eival should be found instead of full Hessian eival
    logical                  :: debug_arh_hessian
    logical                  :: cfg_orbspread
    type(orbspread_data), pointer :: orbspread_input
end type DDitem

!> \brief Contains debug info (work in progress!)
!> \author S. Host
!> \date March 2010
type debugItem
      !> SCF iteration
      integer     :: scfit
      !> HOMO-LUMO gap obtained by diagonalization
      real(realk) :: diag_hlgap
      !> HOMO-LUMO gap obtained iteratively
      real(realk) :: iter_hlgap
      !> Lowest full 2nd order Hessian eigenvalue
      real(realk) :: heseival
      !> Lowest ARH Hessian eigenvalue
      real(realk) :: arheival
end type debugItem

contains

    !> \brief Initialize DDitem and determined number of needed orbital energies.
    !> \author S. Host
    !> \date 2005
    subroutine DD_INIT(decomp,DD)
    implicit none
       !> Contains settings for decomposition and OAO decomposed overlap
       type(decompItem), intent(inout) :: decomp
       !> Contains settings for HOMO-LUMO gap and Hessian eigenvalues
       type(DDitem),intent(inout) :: DD
       integer                    :: matdim, i
       integer                    :: nexci, nexci_a, nexci_b

       DD%cfg_orbspread = .false.
       !Only one vector is needed if
       !1) we only want the homo-lumo gap or
       !2) we make mo start guess for exc energies
       DD%nocc = 1 ; DD%nvirt = 1
       matdim = decomp%S%nrow
       DD%fock_thresh=1.0E-6_realk

       !Generate several initial guesses if several exc. energies are requested!

       !If declared in input, decide here how many start vectors should be used:
       if (decomp%cfg_startvectors) then
          DD%nocc = decomp%cfg_no_of_startvectors ; DD%nvirt = decomp%cfg_no_of_startvectors 
       else if (.not. decomp%cfg_rsp_mostart) then
          if (decomp%cfg_rsp_nexcit > 0) then
             !Test that there are actually as many exc energies as requested!
             if (decomp%cfg_unres) then
                nexci_a = decomp%nocca*(decomp%S%nrow - decomp%nocca)
                nexci_b = decomp%noccb*(decomp%S%nrow - decomp%noccb)
                if (nexci_a < decomp%cfg_rsp_nexcit) then
                   write(decomp%lupri,*)
                   write(decomp%lupri,'("You have requested too many excitation energies!")')
                   write(decomp%lupri,'("There are ", i3, " occupied alpha orbitals and ", i3, " virtual alpha orbitals")') &
     &               decomp%nocca, decomp%S%nrow - decomp%nocca
                   write(decomp%lupri,'("giving ", i3, " alpha excitation energies!")') nexci_a
                   write(decomp%lupri,'("Number of excitation energies is reduced to ", i4)') nexci_a
                   decomp%cfg_rsp_nexcit = nexci_a
                endif               
                if (nexci_b < decomp%cfg_rsp_nexcit) then
                   write(decomp%lupri,*)
                   write(decomp%lupri,'("You have requested too many excitation energies!")')
                   write(decomp%lupri,'("There are ", i3, " occupied beta orbitals and ", i3, " virtual beta orbitals")') &
     &               decomp%noccb, decomp%S%nrow - decomp%noccb
                   write(decomp%lupri,'("giving ", i3, " beta excitation energies!")') nexci_b
                   write(decomp%lupri,'("Number of excitation energies is reduced to ", i4)') nexci_b
                   decomp%cfg_rsp_nexcit = nexci_b                
                endif
             else
                nexci = decomp%nocc*(decomp%S%nrow - decomp%nocc)
                if (nexci < decomp%cfg_rsp_nexcit) then
                   write(decomp%lupri,*)
                   write(decomp%lupri,'("You have requested too many excitation energies!")')
                   write(decomp%lupri,'("There are ", i3, " occupied orbitals and ", i3, " virtual orbitals")') &
                   decomp%nocc, decomp%S%nrow - decomp%nocc
                   write(decomp%lupri,'("giving ", i3, " excitation energies!")') nexci
                   write(decomp%lupri,'("Number of excitation energies is reduced to ", i3)') nexci
                   decomp%cfg_rsp_nexcit = nexci
                endif
             endif  
             DD%nocc = decomp%cfg_rsp_nexcit ; DD%nvirt = decomp%cfg_rsp_nexcit
          else if (decomp%cfg_hessian_nvec > 1) then
             DD%nocc = decomp%cfg_hessian_nvec ; DD%nvirt = decomp%cfg_hessian_nvec
          endif
       endif

       if (.not. decomp%cfg_unres) then
          write(decomp%lupri,*)
          write(decomp%lupri,*) 'Number of occupied orbitals:', decomp%nocc
          write(decomp%lupri,*) 'Number of virtual orbitals:', decomp%S%nrow - decomp%nocc
          write(decomp%lupri,*)

          if (DD%nocc > decomp%nocc) DD%nocc = decomp%nocc
          if (DD%nvirt > (decomp%S%nrow - decomp%nocc)) DD%nvirt = decomp%S%nrow - decomp%nocc

          write(decomp%lupri,*) 'Number of occupied orbital energies to be found:', DD%nocc
          write(decomp%lupri,*) 'Number of virtual orbital energies to be found:', DD%nvirt
          write(decomp%lupri,*)
       else
          write(decomp%lupri,*)
          write(decomp%lupri,*) 'Number of occupied alpha orbitals:', decomp%nocca
          write(decomp%lupri,*) 'Number of occupied beta orbitals:', decomp%noccb
          write(decomp%lupri,*)
          write(decomp%lupri,*) 'Number of virtual alpha orbitals:', decomp%S%nrow - decomp%nocca
          write(decomp%lupri,*) 'Number of virtual beta orbitals:', decomp%S%nrow - decomp%noccb
          write(decomp%lupri,*)

          DD%nocca = DD%nocc ; DD%noccb = DD%nocc
          DD%nvirta = DD%nvirt ; DD%nvirtb = DD%nvirt

          if (DD%nocca > decomp%nocca) DD%nocca = decomp%nocca
          if (DD%noccb > decomp%noccb) DD%noccb = decomp%noccb

          if (DD%nvirta > (decomp%S%nrow - decomp%nocca)) DD%nvirta = decomp%S%nrow - decomp%nocca
          if (DD%nvirtb > (decomp%S%nrow - decomp%noccb)) DD%nvirtb = decomp%S%nrow - decomp%noccb

          write(decomp%lupri,*) 'Number of occupied alpha orbital energies to be found:', DD%nocca
          write(decomp%lupri,*) 'Number of occupied beta orbital energies to be found:', DD%noccb
          write(decomp%lupri,*)
          write(decomp%lupri,*) 'Number of virtual alpha orbital energies to be found:', DD%nvirta
          write(decomp%lupri,*) 'Number of virtual beta orbital energies to be found:', DD%nvirtb
          write(decomp%lupri,*)
       endif

       nullify(DD%FUP_eivecs)
       nullify(DD%FUQ_eivecs)
       nullify(DD%FUP_eival)
       nullify(DD%FUQ_eival)
       nullify(DD%FUPalfa_eivecs)
       nullify(DD%FUQalfa_eivecs)
       nullify(DD%FUPalfa_eival)
       nullify(DD%FUQalfa_eival)
       nullify(DD%FUPbeta_eivecs)
       nullify(DD%FUQbeta_eivecs)
       nullify(DD%FUPbeta_eival)
       nullify(DD%FUQbeta_eival)
       nullify(DD%Ared)
       nullify(DD%boverlaps)

       if (decomp%cfg_unres) then
          allocate(DD%FUPalfa_eivecs(matdim,DD%nocca))
          allocate(DD%FUQalfa_eivecs(matdim,DD%nvirta))
          allocate(DD%FUPalfa_eival(DD%nocca))
          allocate(DD%FUQalfa_eival(DD%nvirta))
          allocate(DD%FUPbeta_eivecs(matdim,DD%noccb))
          allocate(DD%FUQbeta_eivecs(matdim,DD%nvirtb))
          allocate(DD%FUPbeta_eival(DD%noccb))
          allocate(DD%FUQbeta_eival(DD%nvirtb))
       else
          allocate(DD%FUP_eivecs(matdim,DD%nocc))
          allocate(DD%FUQ_eivecs(matdim,DD%nvirt))
          allocate(DD%FUP_eival(DD%nocc))
          allocate(DD%FUQ_eival(DD%nvirt))
       endif
    end subroutine DD_INIT

    !> \brief Free matrices in DDitem.
    !> \author S. Host
    !> \date 2005
    subroutine DD_SHUTDOWN(decomp,DD)
    implicit none
       !> Contains settings for decomposition and OAO decomposed overlap
       type(decompItem), intent(in) :: decomp
       !> Contains settings for HOMO-LUMO gap and Hessian eigenvalues
       type(DDitem),intent(inout)   :: DD
       if (decomp%cfg_unres) then
          deallocate(DD%FUPalfa_eivecs)
          deallocate(DD%FUQalfa_eivecs)
          deallocate(DD%FUPalfa_eival)
          deallocate(DD%FUQalfa_eival)
          deallocate(DD%FUPbeta_eivecs)
          deallocate(DD%FUQbeta_eivecs)
          deallocate(DD%FUPbeta_eival)
          deallocate(DD%FUQbeta_eival)
       else
          deallocate(DD%FUP_eivecs)
          deallocate(DD%FUQ_eivecs)
          deallocate(DD%FUP_eival)
          deallocate(DD%FUQ_eival)
       endif
       nullify(DD%FUP_eivecs)
       nullify(DD%FUQ_eivecs)
       nullify(DD%FUP_eival)
       nullify(DD%FUQ_eival)
       nullify(DD%FUPalfa_eivecs)
       nullify(DD%FUQalfa_eivecs)
       nullify(DD%FUPalfa_eival)
       nullify(DD%FUQalfa_eival)
       nullify(DD%FUPbeta_eivecs)
       nullify(DD%FUQbeta_eivecs)
       nullify(DD%FUPbeta_eival)
       nullify(DD%FUQbeta_eival)
       nullify(DD%Ared)
       nullify(DD%boverlaps)
    end subroutine DD_SHUTDOWN

!> \brief Main driver for getting HOMO-LUMO gap and, if requested, lowest Hessian eigenvlaue(s). 
!> \author S. Host
!> \date 2005
!>
!> ==================================================================================== \n
!> Calculate final HOMO-LUMO gap, and, if requested, lowest Hessian eigenvalue \n
!> The homo-lumo gap is inexpensive and is always calculated. From the orbital
!> energies, we obtain starting guesses for both excitation energies and Hessian eigenvalues.
!> MAKE SURE THAT THE PROPER DU, FU, FUP, and FUQ ARE CREATED BEFORE CALLING THIS ROUTINE
!> (call to oao_get_transformed_matrices). \n 
!> Hessian eigenvalue and orbital energy solvers use different routines, since
!> one uses type(matrix) and the other uses real(realk) matrices.
!> ====================================================================================
!>
   subroutine DD_homolumo_and_heseigen(DD,decomp,debug,do_rsp_iniguess,howmany,rsp_iniguess,fifoqueue)
   use rspPrecond
   implicit none
          !> Contains settings for HOMO-LUMO gap and Hessian eigenvalues
          type(DDitem),intent(inout)  :: DD
          !> Contains settings for decomposition and OAO decomposed overlap
          type(decompItem),intent(inout) :: decomp
          !> Contains debug info
          type(debugItem),intent(inout) :: debug
          !> True if we are calculating response initial guess(es)
          logical, intent(in)      :: do_rsp_iniguess
          !> Number of rsp iniguesses to be calculated. Not referenced if do_rsp_iniguess=.false.
          integer,intent(in)       :: howmany
          !> Rsp initial guesses (output)
          type(matrix),intent(inout),optional :: rsp_iniguess(howmany)
          !> Contains subspace of previous density and Fock matrices
          type(modFIFO), intent(inout), optional :: fifoqueue
          type(matrix),allocatable :: iniguess(:)
          real(realk), allocatable :: iniguess_full(:,:), orb_energy_diff(:)
          integer, allocatable     :: FUQFUPindex(:,:)
          real(realk)              :: fac, hessian_eigenval(decomp%cfg_hessian_nvec), HOMOLUMOgap, T1, T2, norm
          !real(realk)              :: cpu1, cpu2, wall1, wall2
          integer                  :: mu, nu, matdim, i, k, j, posmin(1), posmax(1), fulldim, nguesses
          type(matrix)             :: scr, scr2
          type(modFIFO)            :: dummyfifo
          logical                  :: printgap
          integer                  :: KK1,KK2
          real(realk)              :: TMP1,TMP2

   printgap = .true.
   DD%debug_arh_hessian = .false.

   if (decomp%cfg_check_converged_solution .and. do_rsp_iniguess) then
      call lsquit('Currently not possible to do both Hessian eigenvalues and excitation energies',decomp%lupri)
   endif

   if (do_rsp_iniguess) then
      if (.not. present(rsp_iniguess)) then
         call lsquit('rsp_iniguess must be present when do_rsp_iniguess=.true.',decomp%lupri)
      endif
   endif

   matdim = decomp%FUP%nrow
   fulldim = matdim
   !finding_arh_eigenvalue = .false.

   CALL LSTIMER('START ',T1,T2,decomp%lupri)
   call time_II_operations1
   !call gettim(cpu1,wall1)
   write(decomp%lupri,*)
   write(decomp%lupri,'("Calculation of HOMO-LUMO gap")')
   write(decomp%lupri,'("============================")')

   call DD_Fock_eigenvalue(DD,decomp,'h',0.0E0_realk)
   if (.not. decomp%cfg_hlgap_converged) then
      write (decomp%lupri,*)
      write (decomp%lupri,'("Calculation of occupied orbital energies did not converge")')
      printgap = .false.
   else
      write (decomp%lupri,*)
      write (decomp%lupri,'("Calculation of occupied orbital energies converged in", i6, " iterations!")') &
      & decomp%cfg_hlgap_nit_conv
   endif
   if (DD%nvirt > 0) then !Fix for atoms
      call DD_Fock_eigenvalue(DD,decomp,'l',0.0E0_realk)
      if (.not. decomp%cfg_hlgap_converged) then
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of virtual orbital energies did not converge")')
         printgap = .false.
      else
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of virtual orbital energies converged in", i6, " iterations!")') &
         & decomp%cfg_hlgap_nit_conv
      endif
      if (printgap) then
         HOMOLUMOgap = DD%FUQ_eival(1)-DD%FUP_eival(1)
         write (decomp%lupri,*)
         write (decomp%lupri,'("    E(LUMO):                     ", F12.6, " au")') DD%FUQ_eival(1)
         write (decomp%lupri,'("   -E(HOMO):                     ", F12.6, " au")') DD%FUP_eival(1)
         write (decomp%lupri,'("   ------------------------------")')
         write (decomp%lupri,'("    HOMO-LUMO Gap (iteratively): ", F12.6, " au")') HOMOLUMOgap
         write (decomp%lupri,*)
         debug%iter_hlgap = HOMOLUMOgap
      else
         write (decomp%lupri,*)
         write (decomp%lupri,*) "HOMO-LUMO gap cannot be printed"
         write (decomp%lupri,*)
      endif
   else
      write (decomp%lupri,*)
      write (decomp%lupri,*) 'No virtual orbitals - atomic calculation?'
      write (decomp%lupri,*)
   endif

   !call gettim(cpu2,wall2)
   !call util_print_time('HL GAP',cpu1,cpu2,wall1,wall2)
   CALL LSTIMER('HL GAP',T1,T2,decomp%lupri)
   call time_II_operations2(job_homolumo)
   if (decomp%cfg_check_converged_solution .or. do_rsp_iniguess) then

      if (decomp%cfg_check_converged_solution) then
         if (decomp%cfg_startvectors) then
            nguesses = decomp%cfg_no_of_startvectors
         else
            nguesses = decomp%cfg_hessian_nvec
         endif
      else
         nguesses = howmany
      endif
      allocate(iniguess(nguesses))

      do i = 1, nguesses
         call mat_init(iniguess(i), matdim, matdim)
      enddo

      allocate(iniguess_full(fulldim,fulldim))
      allocate(orb_energy_diff(DD%nocc*DD%nvirt))
      allocate(FUQFUPindex(DD%nocc*DD%nvirt,2))

      if (decomp%cfg_check_converged_solution) then
         write(decomp%lupri,*)
         write(decomp%lupri,'("Checking stability of solution: Calculation of lowest Hessian eigenvalue")')
         write(decomp%lupri,'("========================================================================")')
      else
         write(decomp%lupri,*)
         write(decomp%lupri,'("Calculating initial guess(es) for excitation energies from orb. energies")')
         write(decomp%lupri,'("========================================================================")')
      endif

      if (decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) '       Occ. orbitals:' !; call flshfo(lupri)
         write(decomp%lupri,*) '       =============='
         CALL LS_OUTPUT(DD%FUP_eival,1,DD%nocc,1,1,DD%nocc,1,1,decomp%LUPRI)
         !CALL LS_OUTPUT(DD%FUP_eivecs,1,matdim,1,DD_no_of_vecs,matdim,DD_no_of_vecs,1,LUPRI)

         write(decomp%lupri,*)
         write(decomp%lupri,*) '       Virt. orbitals:' !; call flshfo(lupri)
         write(decomp%lupri,*) '       ==============='
         CALL LS_OUTPUT(DD%FUQ_eival,1,DD%nvirt,1,1,DD%nvirt,1,1,decomp%LUPRI)
         !CALL LS_OUTPUT(DD%FUQ_eivecs,1,matdim,1,DD_no_of_vecs,matdim,DD_no_of_vecs,1,LUPRI)
      endif

      k = 0
      !Calculate orbital energy differences obtained from FUP/FUP eigenvalues
      do i = 1, DD%nvirt !FUQeival
         do j = 1, DD%nocc !FUPeival
            k = k + 1
            orb_energy_diff(k) = DD%FUQ_eival(i) - DD%FUP_eival(j)
            FUQFUPindex(k,1) = j
            FUQFUPindex(k,2) = i
         enddo
      enddo

      if (decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) '  Orbital energy differences  (occ index, virt index):' !; call flshfo(lupri)
         write(decomp%lupri,*) '======================================================'
         do i = 1, DD%nocc*DD%nvirt
            write(decomp%lupri, '(i5, F12.7, i19, i11)') i, orb_energy_diff(i), FUQFUPindex(i,1), FUQFUPindex(i,2)
         enddo

         write(decomp%lupri,*)
         write(decomp%lupri,*) '    Sorted orbital energy differences:'
         write(decomp%lupri,*) '=========================================='
      endif

      do i = 1, nguesses
         posmin = MINLOC(orb_energy_diff)
         posmax = MAXLOC(orb_energy_diff)
         if (decomp%info_stability) then
            write (decomp%lupri, "('Root no.', i4, F12.7, '   Position:', i4)") i, orb_energy_diff(posmin(1)),  posmin(1)
         endif
         if (do_rsp_iniguess) then
            KK1 = FUQFUPindex(posmin(1),1)
            KK2 = FUQFUPindex(posmin(1),2)
            do nu = 1, fulldim
               TMP1 = DD%FUP_eivecs(nu,KK2)
               do mu = 1, fulldim
                  iniguess_full(mu,nu) = DD%FUQ_eivecs(mu,KK1)*TMP1
               enddo
            enddo
         else
            KK1 = FUQFUPindex(posmin(1),1)
            KK2 = FUQFUPindex(posmin(1),2)
            do nu = 1, fulldim
               TMP1 = DD%FUP_eivecs(nu,KK2)
               TMP2 = DD%FUQ_eivecs(nu,KK2)
               do mu = 1, fulldim
                  iniguess_full(mu,nu) = DD%FUQ_eivecs(mu,KK1)*TMP1- &
                       &  DD%FUP_eivecs(mu,KK1)*TMP2
               enddo
            enddo
         endif
         !write(lupri,*) 'full Hes starting guess no', i
         !call LS_OUTPUT(iniguess_full, 1, fulldim, 1, fulldim, fulldim, fulldim, 1, lupri)

         orb_energy_diff(posmin(1)) = orb_energy_diff(posmax(1))

         call mat_set_from_full(iniguess_full,1.0E0_realk, iniguess(i), 'iniguess')
         if (do_rsp_iniguess) then
            !Normalize X*S*X=1
            call mat_init(scr,matdim,matdim)
            call mat_init(scr2,matdim,matdim)

            call oao_rsp_lintrans(decomp,iniguess(i),scr2,scr) !scr = rho
            fac = abs(mat_dotproduct(scr,iniguess(i)))
            call mat_scal(1.0E0_realk/sqrt(fac),iniguess(i))
            call x_from_oao_basis(decomp,iniguess(i),rsp_iniguess(i))

            call mat_free(scr)
            call mat_free(scr2)
         else
            !Normalize X*X=1 
            fac = abs(mat_dotproduct(iniguess(i),iniguess(i)))
            call mat_scal(1.0E0_realk/sqrt(fac),iniguess(i))
            !write(lupri,*) 'Hes starting guess no', i
            !call mat_print(iniguess(i),1,matdim,1,matdim,lupri)
         endif
      enddo
      !call flshfo(lupri)

      if (decomp%cfg_check_converged_solution) then
      !if (.not. debug_no_hessian) then !A bit ugly, I know... for debugging arh hessian which is
                                       !changed when we do a rejection, as opposed to the 'real' Hessian 
         !write(lupri,*) 'Calling hes eival' ; call flshfo(lupri)
         call DD_Hessian_eigenvalue(decomp, DD, iniguess, nguesses, hessian_eigenval, dummyfifo)

         write(decomp%lupri,*)
         write (decomp%lupri,'("*======================================================================*")')
         do i = 1, decomp%cfg_hessian_nvec
            write (decomp%lupri,'("                Hessian eigenvalue no.",i3,": ",F12.6, "      ")') i, hessian_eigenval(i) 
         enddo
         write (decomp%lupri,'("*======================================================================*")')
         debug%heseival = hessian_eigenval(1)
      endif

      if (present(fifoqueue)) then
         DD%debug_arh_hessian = .true.
      else
         DD%debug_arh_hessian = .false.
      endif

      if (DD%debug_arh_hessian) then
         call DD_Hessian_eigenvalue(decomp, DD, iniguess, nguesses, hessian_eigenval, fifoqueue)

         write(decomp%lupri,*)
         write(decomp%lupri,'("*======================================================================*")')
         write(decomp%lupri,'("                ARH lowest Hessian eigenvalue: ",F12.6, "      ")') hessian_eigenval(1)
         write(decomp%lupri,'("*======================================================================*")')
         debug%arheival = hessian_eigenval(1)
      endif

      do i = 1, nguesses
         call mat_free(iniguess(i))
      enddo
      deallocate(iniguess_full)
      deallocate(orb_energy_diff)
      deallocate(FUQFUPindex)

      !call gettim(cpu2,wall2)
      !write(lupri,*) 'cputime, walltime:', cpu2-cpu1, wall2-wall1
      !call util_print_time('HESEIG',cpu1,cpu2,wall1,wall2)
      CALL LSTIMER('HESEIG',T1,T2,decomp%lupri)
   endif
   end subroutine DD_homolumo_and_heseigen

   !> \brief Find lowest or highest eigenvalue of n x n matrix M by inverse iteration and reduced space
   !> \author S. Host
   !> \date 2005
   subroutine DD_Fock_eigenvalue(DD, decomp, desired_eigenval, guess)
   implicit none
         !> Contains settings for HOMO-LUMO gap and Hessian eigenvalues
         type(DDitem),intent(inout)  :: DD
         !> Contains settings for decomposition and OAO decomposed overlap
         type(decompItem),intent(inout) :: decomp
         !> h = find highest eigenvalue(s), l = find lowest eigenvalue(s)
         character, intent(in)    :: desired_eigenval
         !> Should be close to the desired eigenvalue for efficient solution.
         real(realk), intent(in)  :: guess 
         integer                  :: i, j, ndim, nroots
         real(realk), allocatable, dimension(:,:) :: iniguess_full, M_full, P_full, Mdamp, x
         real(realk), allocatable, dimension(:) :: eigenval

      if (desired_eigenval == 'l') then
         nroots = DD%nvirt
      else if (desired_eigenval == 'h') then
         nroots = DD%nocc
      endif

      ndim = decomp%U%nrow

      allocate(M_full(ndim,ndim), Mdamp(ndim,ndim))
      allocate(P_full(ndim,ndim))
      allocate(iniguess_full(ndim,nroots),x(ndim,nroots))
      allocate(eigenval(nroots))

      !write(lupri,*) 'M:'
      !call mat_print(M, 1, M%nrow, 1, M%ncol, lupri)
      !write(lupri,*) 'P:'
      !call mat_print(P, 1, M%nrow, 1, M%ncol, lupri)

      M_full = 0.0E0_realk; P_full = 0.0E0_realk
      if (desired_eigenval == 'l') then
         call mat_to_full(decomp%FUQ, 1.0E0_realk, M_full)
         call mat_to_full(decomp%QU, 1.0E0_realk, P_full)
      else if (desired_eigenval == 'h') then
         call mat_to_full(decomp%FUP, 1.0E0_realk, M_full)
         call mat_to_full(decomp%DU, 1.0E0_realk, P_full)
      endif

      !WRITE(LUPRI,*)'M_full'
      !CALL LS_OUTPUT(M_full,1,ndim,1,ndim,ndim,ndim,1,LUPRI)
      !WRITE(LUPRI,*)'P_full'
      !CALL LS_OUTPUT(P_full,1,ndim,1,ndim,ndim,ndim,1,LUPRI)

      !Damp M:
      Mdamp = M_full
      do i = 1, ndim
         Mdamp(i,i) = Mdamp(i,i) - guess
      enddo 

      !Find nroots starting guesses:
      call DD_initial_guess(M_full,ndim,nroots,iniguess_full)

      !Project starting guesses: This is done in orthonormalize inside solver
      x = iniguess_full

      if (desired_eigenval == 'h') then !Switch sign on matrix
         M_full = -M_full
         call DD_Fock_solver(decomp, DD, ndim, nroots, Mdamp, M_full, x, P_full, eigenval)
         if (decomp%cfg_hlgap_converged) then
            !eigenval = -eigenval
            DD%FUP_eival = -eigenval
            !Save eigenvectors so they can be used to generate starting guesses for Hessian eigenvalue
            !and/or excitation energies
            do i = 1, DD%nocc
               call DD_Fock_project(ndim,P_full,x(1:ndim,i),DD%FUP_eivecs(1:ndim,i))
            enddo
         endif         !eigenval = -eigenval
      else if (desired_eigenval == 'l') then
         call DD_Fock_solver(decomp,DD,ndim,nroots, Mdamp, M_full, x, P_full, eigenval)
         if (decomp%cfg_hlgap_converged) then
	    DD%FUQ_eival = eigenval
            !Save eigenvectors so they can be used to generate starting guesses for Hessian eigenvalue
            !and/or excitation energies
            do i = 1, DD%nvirt
               call DD_Fock_project(ndim,P_full,x(1:ndim,i),DD%FUQ_eivecs(1:ndim,i))
            enddo
         endif
      else
         STOP 'Unknown type of eigenvalue (DD_eigenvalue)'
      endif

      deallocate(x)
      deallocate(M_full,Mdamp)
      deallocate(P_full)
      deallocate(iniguess_full)
      deallocate(eigenval)
   end subroutine DD_Fock_eigenvalue 

   !> \brief Wrapper for obtaining lowest Hessian eigenvalue(s).
   !> \author S. Host
   !> \date 2005
   subroutine DD_Hessian_eigenvalue(decomp, DD, iniguess, nstart, eigenval, fifoqueue)
   implicit none
         !> Contains settings for decomposition and OAO decomposed overlap
         type(decompItem),intent(inout) :: decomp
         !> Contains settings concerning HOMO-LUMO gap and Hessian eigenvalues
         type(DDitem),intent(inout)  :: DD
         !> Number of startvectors in iniguess (default is decomp%cfg_hessian_nvec)
         integer, intent(in)      :: nstart
         !> Initial guess(es) for Hessian eigenvalue iterative solver
         type(matrix), intent(in) :: iniguess(nstart)
         !> Hessian eigenvalue(s)
         real(realk), intent(out) :: eigenval(decomp%cfg_hessian_nvec)
         !> Contains density and Fock matrices from previous iterations (for debugging only)
         type(modFIFO), intent(inout) :: fifoqueue
         type(matrix)             :: x(nstart)
         real(realk)              :: norm
         integer                  :: ndim, i

      ndim = iniguess(1)%nrow

      do i = 1, nstart
         call mat_init(x(i), ndim, ndim)
      enddo

      do i = 1, nstart
         !norm = sqrt(mat_sqnorm2(iniguess(i)))
         !write(lupri,*) 'DD_Hessian_eigenvalue: norm of Hes starting guess no', i, norm 

         call project_oao_basis(decomp, iniguess(i), 2, x(i))
         !call DD_project('h',iniguess(i),x(i)) !x = projected iniguess
         !norm = sqrt(mat_sqnorm2(x(i)))
         !write(lupri,*) 'DD_Hessian_eigenvalue: norm of proj starting guess no', i, norm 
      enddo
      call DD_solver(decomp, DD, x, nstart, eigenval, fifoqueue)

      do i = 1, nstart
         !write (lupri,*) 'in loop, deallocating vector no', i ; call flshfo(lupri)
         call mat_free(x(i))
         !write (lupri,*) 'in loop, after deallocating vector no', i ; call flshfo(lupri)
      enddo
      !write (lupri,*) 'exit from DD_Hessian_eigenvalue' ; call flshfo(lupri)

   end subroutine DD_Hessian_eigenvalue 

   !> \brief Obtain initial guess(es) for orbital eigenvalue iterative solver.
   !> \author S. Host
   !> \date 2005
   subroutine DD_initial_guess(M_full,ndim,nroots,iniguess_full)
   implicit none
         !> Number of basis functions
         integer, intent(in)        :: ndim
         !> FUP (if occupied eigenvalue are requested) or FUQ (if virtual eigenvalue are requested) 
         real(realk), intent(in)    :: M_full(ndim,ndim)
         !> Number of requested orbital energies (occupied or virtual)
         integer, intent(in)        :: nroots
         !> The initial guess(es) for the orbital eigenvalue iterative solver
         real(realk), intent(out)   :: iniguess_full(ndim,nroots)
         integer, allocatable       :: positions(:)
         integer                    :: i, pos(1)
         real(realk), allocatable   :: values(:), diag(:)

      allocate(positions(nroots))
      allocate(values(nroots))
      allocate(diag(ndim))

      !WRITE(LUPRI,*)'M_full'
      !CALL LS_OUTPUT(M_full,1,ndim,1,ndim,ndim,ndim,1,LUPRI)

      do i = 1, ndim
         diag(i) = M_full(i,i)
      enddo
      diag = abs(diag)

      !WRITE(LUPRI,*)'diag'
      !CALL LS_OUTPUT(diag,1,ndim,1,1,ndim,1,1,LUPRI)

      do i = 1, nroots
         pos = MAXLOC(diag)
         positions(i) = pos(1)
         diag(pos(1)) = 0.0E0_realk
      enddo

      !write (lupri,*) 'Positions of largest (absolute) diagonal elements of FUP/FUQ:', positions
      
      iniguess_full = 0.0E0_realk
      do i = 1, nroots
         iniguess_full(positions(i),i) = 1.0E0_realk
         !WRITE(LUPRI,*)'iniguess no.', i 
         !CALL LS_OUTPUT(iniguess_full,1,ndim,i,i,ndim,nroots,1,LUPRI)
      enddo

      deallocate(positions)
      deallocate(values)
      deallocate(diag)
   end subroutine DD_initial_guess

   !> \brief Preconditioned Conjugate Gradient method
   !> \author S. Host
   !> \date March 2005. Revised to fit the inverse iteration routine March 2006.
   !>  
   !> Used for preconditioning when solving orbital eigenvalue equations.
   !> Iterative solution of the linear equations A*x=b
   !>
   subroutine DD_Fock_PCG(decomp, vec_number, ndim, A, b, proj, x)
   implicit none
         !> Contains settings for decomposition and OAO decomposed overlap
         type(decompItem),intent(in) :: decomp
         !> Currently preconditioning residual with this number (for printout only)
         integer, intent(in)      :: vec_number
         !> Number of basis functions
         integer, intent(in)      :: ndim
         !> Left hand side of linear equations A*x=b
         real(realk), intent(in)  :: A(ndim,ndim)
         !> Right hand side of linear equations A*x=b
         real(realk), intent(in)  :: b(ndim)
         !> After diagonal preconditioning, project with this matrix
         real(realk), intent(in)  :: proj(ndim,ndim)
         !> Solution to system of linear equations A*x=b
         real(realk), intent(out) :: x(ndim) 
         integer                  :: N, j, i, k !, vecdim
         real(realk)              :: alfa, beta, err, alfa_num, alfa_den, beta_num, beta_den, &
                                   & thresh, norm
         real(realk),allocatable, dimension(:) :: p, Hr, Hr_new, res, scr

   allocate(p(ndim))
   allocate(Hr(ndim))
   allocate(Hr_new(ndim))
   allocate(scr(ndim))
   allocate(res(ndim))

   !Threshold of convergence:
   !norm = sqrt(mat_sqnorm2(b))
   call VECNORM(b, ndim, norm)
   thresh = norm*1.0E-2_realk
   !Initialize iteration number and x:
   N=0
   !call MAT_ZERO(x)
   x = 0.0E0_realk

   !Solving:
   do
      if (N == 15 + 1) then
         if (decomp%info_stability) write(decomp%lupri,'("          Forced exit after ", i4, " PCG iterations")') N
         exit
      endif
      !If first iteration, initialize:
      if (N == 0) then
         call DD_Fock_lintrans(ndim,A,x,scr)
         res = b - scr
         call VECNORM(res, ndim, err)
         if (decomp%info_stability) then
            write (decomp%LUPRI, '("          VECTOR ", i3, ": Error of PCG iteration", i3, ":  ", E18.10)') vec_number, N, err
         endif
         call DD_Fock_precond(ndim,A,res,scr)
         call DD_Fock_project(ndim,Proj,scr,p)
         Hr = p
      else
         beta_den = dot_product(res, Hr)
         call DD_Fock_lintrans(ndim,A,p,scr)
         alfa_num = beta_den
         alfa_den = dot_product(p, scr)
         if (ABS(alfa_den) < 1.0E-12_realk) then
            if (decomp%info_stability) then
               write (decomp%lupri,*) "Warning: Exit from PCG loops due to denominator = 0!"
            endif
            if (N == 1) then
               x = b
            endif
            EXIT
         endif
         alfa = alfa_num/alfa_den
         x = x + alfa*p
         res = res - alfa*scr     !r(n+1) = r(n) - alfa*A*p(n)
         !Convergence of PCG?
         call VECNORM(res, ndim, err)
         if (decomp%info_stability) then
            write (decomp%LUPRI, '("          VECTOR ", i3, ": Error of PCG iteration", i3, ":  ", E18.10)') vec_number, N, err
         endif
         if (err < thresh) then
            if (decomp%info_stability) then
               write (decomp%LUPRI, '("          VECTOR ", i3, ": PCG converged in ", i3, " iterations! ")') vec_number, N
            endif
            EXIT
         endif
         call DD_Fock_precond(ndim, A,res,scr)
         call DD_Fock_project(ndim,Proj,scr,Hr_new)
         beta_num = dot_product(res, Hr_new)
         beta = beta_num/beta_den
         p = Hr_new + beta*p
         Hr = Hr_new
      endif
      N = N + 1 
   enddo

   deallocate(p)
   deallocate(Hr)
   deallocate(Hr_new)
   deallocate(scr)
   deallocate(res)
   end subroutine DD_Fock_PCG

   !> \brief Eigenvalue solver for obtaining orbital energies.
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_Fock_solver(decomp, DD, ndim, nroots, Mdamp, M, iniguess, Proj, mu)
   implicit none
          !> Contains settings for decomposition and OAO decomposed overlap
          type(decompItem),intent(inout) :: decomp
          !> Contains settings concerning HOMO-LUMO gap and Hessian eigenvalues
          type(DDItem),intent(inout)  :: DD
          !> Number of basis functions
          integer, intent(in)      :: ndim
          !> Number of orbital energies to be determined
          integer, intent(in)      :: nroots
          !> Obtain nroots lowest eigenvalues of this matrix
          real(realk), intent(in)  :: M(ndim,ndim)
          !> Matrix M damped by initial guess (inverse iteration algorithm)
          real(realk), intent(in)  :: Mdamp(ndim,ndim)
          !> Initial guess(es) for eigenvalue solver
          real(realk), intent(inout) :: iniguess(ndim,nroots)
          !> After preconditioning, project with this matrix
          real(realk), intent(in)  :: Proj(ndim,ndim)
          !> The final eigenvalue(s)
          real(realk), intent(out) :: mu(nroots)
          real(realk)              :: val, err(nroots), norm, thresh
          real(realk),allocatable,dimension(:,:) :: b_current, scrvec, x, res, sigma_cur
          integer                  :: max_it, pos, iter, i, j 
          integer                  :: n_on_disk, ndim_red, nb_curr, nnew, nconv
          logical                  :: err_went_up, conv(nroots)

    max_it = decomp%cfg_homolumo_maxit
    DD%lusigma = -1 ; DD%lub = -1
    CALL lsOPEN(DD%lusigma,'DD_sigmavecs','unknown','UNFORMATTED')
    CALL lsOPEN(DD%lub,'DD_bvecs','unknown','UNFORMATTED')
    allocate(DD%Ared(max_it,max_it), DD%boverlaps(max_it,max_it))
    DD%Ared = 0.0E0_realk
    conv = .false.
    nconv = 0

    call VECNORM(iniguess(1:ndim,1),ndim,norm)
    allocate(b_current(ndim,nroots))
    allocate(sigma_cur(ndim,nroots))
    allocate(scrvec(ndim,nroots))
    allocate(res(ndim,nroots))
    allocate(x(ndim,nroots))
    thresh = norm*DD%fock_thresh
    if (decomp%info_stability) then
       write (decomp%lupri,*) 'Threshold for iterative determination of orbital energies:', thresh
    endif

    !linear_dep = .false.

    !Orthonormalize nroots starting guesses against each other
    nnew = nroots
    n_on_disk = 0
    ndim_red  = 0
    ! IN : nnew is number of new non-orth vectors
    ! OUT : nnew is number of accepted orth vectors, written on disk
    ! IN (not changed) : n_on_disk = number of vecs on disk before call to DD_orthonormalize
    call DD_orthonormalize(decomp,nroots,iniguess,ndim,nnew,n_on_disk,conv,Proj,DD%lub,b_current)
    nb_curr = nnew  !no of vecs in b_current
    n_on_disk = n_on_disk + nnew

    do i = 1, nb_curr
       call DD_Fock_lintrans(ndim,M,b_current(1:ndim,i),sigma_cur(1:ndim,i))
       write(DD%lusigma) sigma_cur(1:ndim,i)
    enddo

    !call flshfo(decomp%lupri)
    iter = 0
    do  
       if (iter + nroots - 1 >= decomp%cfg_homolumo_maxit) then
          if (.not. decomp%cfg_hlgap_needed) then !just print a warning, don't quit
             WRITE(decomp%LUPRI,'(/A)') &
             &     'WARNING: Max. number of HOMO-LUMO iterations reached '
             write(decomp%lupri,*) 'You may try increasing .HLMAXIT in input. Currently set to', decomp%cfg_homolumo_maxit
             decomp%cfg_hlgap_converged = .false.
             exit
          else
             WRITE(decomp%LUPRI,'(/A)') &
             &     'Max. number of HOMO-LUMO iterations reached '
             write(decomp%lupri,*) 'Increase .HLMAXIT in input. Currently set to', decomp%cfg_homolumo_maxit
             CALL lsQUIT(' Max. number of HOMO-LUMO iterations reached ',decomp%lupri)
          endif
       endif
       iter = iter + 1
       call DD_Fock_red_space(DD,decomp,nroots,ndim,ndim_red,n_on_disk,nb_curr,b_current,sigma_cur,mu,x,res)
       !Convergence?
       conv = .false.
       nconv = 0
       do i = 1, nroots 
          call vecnorm(res(1:ndim,i),ndim,err(i))
          if (err(i) < thresh) then
             conv(i) = .true.
             nconv = nconv + 1
          endif
          if (decomp%info_stability) then
                write (decomp%LUPRI, '("    FUP/FUQ eival: Residual for vector ", &
                & i3, ":  ", E18.10, " and mu: ",E18.10, "  It = ", i3, "  Conv = ", L2)') &
                & i, err(i), mu(i), iter, conv(i)
          endif
       enddo
       if (ALL(conv)) then
          if (decomp%info_stability) then
             write (decomp%lupri,*) '   Convergence of all vectors!'
             do i = 1, nroots 
                write (decomp%lupri, '("   Eigenvalue for vector", i3, " is:", E14.6)') i, mu(i)
             enddo
          endif
          decomp%cfg_hlgap_converged = .true.
          decomp%cfg_hlgap_nit_conv = iter
          exit
       endif
       !Preconditioning of residual = resP and b(iter+1) = resP:
       do i = 1, nroots !Precond each residual individually
          if (.not. conv(i)) then !Only do precond if residual not converged
             call DD_Fock_PCG(decomp, i, ndim, Mdamp, res(1:ndim,i), proj, scrvec(1:ndim,i))
          endif
       enddo
       !Orthonormalize b_current against set of previous b vectors and save:
       nnew = nroots
       ! IN : nnew is number of new non-orth vectors
       ! OUT : nnew is total number of orth vector written on disk
       ! IN (not changed) : n_on_disk = number of vecs on disk before call to DD_orthonormalize
       call DD_orthonormalize(decomp,nroots,scrvec,ndim,nnew,n_on_disk,conv,Proj,DD%lub,b_current)
       nb_curr = nnew !no of vecs in b_current
       n_on_disk = n_on_disk + nnew
       do i = 1, nb_curr
          call DD_Fock_lintrans(ndim,M,b_current(1:ndim,i),sigma_cur(1:ndim,i))
          write(DD%lusigma) sigma_cur(1:ndim,i)
       enddo
       !if (linear_dep) then
       !   write (lupri,*) '   Linear dependencies, exiting'
       !   write (lupri,*) '   Current approximation to lowest eigenvalue, mu =', mu
       !   exit
       !endif
    enddo

    iniguess = x

    deallocate(b_current)
    deallocate(scrvec)
    deallocate(res)
    deallocate(x)
    CALL lsCLOSE(DD%lusigma,'DELETE')
    CALL lsCLOSE(DD%lub,'DELETE')
    deallocate(DD%Ared, DD%boverlaps)
    end subroutine DD_Fock_solver

   !> \brief Eigenvalue solver for obtaining Hessian eigenvalues.
   !> \author S. Host
   !> \date March 2005.
    subroutine DD_solver(decomp, DD, iniguess, nstart, mu, fifoqueue)
    implicit none
          !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
          type(decompItem),intent(inout) :: decomp
          !> Contains settings concerning HOMO-LUMO gap and Hessian eigenvalues
          type(DDitem),intent(inout)  :: DD
          !> Number of start vectors in iniguess (default is decomp%cfg_hessian_nvec)
          integer, intent(in)         :: nstart
          !> Initial guess(es) for Hessian eigenvalue(s)
          type(matrix), intent(inout) :: iniguess(nstart)
          !> The final Hessian eigenvalue(s)
          real(realk), intent(out) :: mu(decomp%cfg_hessian_nvec)
          !> Contains density and Fock matrices from previous iterations (for debugging only)
          type(modFIFO), intent(inout) :: fifoqueue
          real(realk)              :: val, err(decomp%cfg_hessian_nvec), norm, thresh, T1, T2
          type(matrix), allocatable :: b_current(:), scrvec(:), sigma_cur(:), x(:), res(:) 
          integer                  :: max_it, pos, iter, i, j, rowdim, coldim, &
                                    & idum, ldum, heseig_lun, restart_lun, iteration
          logical                  :: err_went_up, restart_from_disk, convergence
          logical, allocatable     :: conv(:)
          integer                  :: nnew, n_on_disk, ndim_red, nb_curr, maxsize
          logical :: OnMaster

    OnMaster=.TRUE.
    max_it = decomp%cfg_check_maxit
    iter = 0      !Size of reduced space. Differs from 'iteration' if no of start vectors /= 1
    iteration = 1
    rowdim = iniguess(1)%nrow
    coldim = iniguess(1)%ncol
    norm = sqrt(mat_sqnorm2(iniguess(1)))
    !restart_from_disk = .false.
    !INQUIRE(file='heseig.restart',EXIST=restart_from_disk)
    !nnew = decomp%cfg_hessian_nvec 
    nnew = nstart
    n_on_disk = 0
    ndim_red  = 0    

    allocate(b_current(nstart))
    allocate(sigma_cur(nstart))
    allocate(res(nstart))
    allocate(x(nstart))
    allocate(scrvec(nstart))
    allocate(conv(nstart))

    conv = .false.

    do i = 1, nstart
       call mat_init(b_current(i),rowdim,coldim)
       call mat_init(sigma_cur(i),rowdim,coldim)
    enddo
    do i = 1, decomp%cfg_hessian_nvec
       call mat_init(res(i),rowdim,coldim)
       call mat_init(x(i),rowdim,coldim)
       call mat_init(scrvec(i),rowdim,coldim)
    enddo

    ! The absolute threshold cannot be larger than this, otherwise it is very
    ! likely that we converge to a solution which is not the lowest eigenvalue
    thresh = 1.0E-6_realk

    DD%lusigma = -1 ; DD%lub = -1
    CALL lsOPEN(DD%lusigma,'hessigma','unknown','UNFORMATTED')
    CALL lsOPEN(DD%lub,'hesbvec','unknown','UNFORMATTED')

    maxsize = max_it + nstart - 1
    allocate(DD%Ared(maxsize,maxsize), DD%boverlaps(maxsize,maxsize))
    DD%Ared = 0.0E0_realk
      !Lowest eigenvalue(s)
    call DD_Hes_orthonormalize(decomp,iniguess,nnew,n_on_disk,conv,DD%lub,b_current)

    nb_curr = nnew  !no of vecs in b_current
    n_on_disk = n_on_disk + nnew

    !write(lupri,*) 'WARNING: hessian restart option removed!' ; call flshfo(lupri)
    !if (restart_from_disk .and. .false.) then
    !   !if (DD_no_of_vecs /= 1) then
    !   !  WRITE(LUPRI,*) 'DD_no_of_vecs not compatible with restart_from_disk '
    !   !  CALL lsQUIT('DD_no_of_vecs not compatible with restart_from_disk',decomp%lupri)
    !   !endif
    !   restart_lun = -1  !initialization
    !   call lsopen(restart_lun,'heseig.restart','OLD', &
    !               & 'UNFORMATTED')
    !   rewind restart_lun
    !   call mat_read_from_disk(restart_lun,scrvec)
    !   call lsclose(restart_lun,'KEEP')
    !   WRITE(LUPRI,*)
    !   WRITE(LUPRI,*) '*** RESTART HES EIGENVAL FROM VECTOR ON DISK - READ FROM heseig.restart  ***'
    !   WRITE(LUPRI,*)
    !   call normalize(scrvec)
    !   call DD_project(scrvec,b_current)
    !   call normalize(b_current)

    !   call mat_write_to_disk(DD%lub,b_current)
    !   call DD_lintrans(0.0E0_realk,scrvec,b_current,scrvec)
    !   call mat_write_to_disk(DD%lusigma,scrvec)
    !   call DD_red_space(decomp,DD,iter,b_current,mu,x,res)
    !else
    !   do i = DD_no_of_vecs, 1, -1 !Do backwords such that b_current corresponds to HOMO-LUMO guess
    !      iter = iter + 1
    !      scrvec = iniguess(i)
    !      call normalize(scrvec)
    !      call DD_project(scrvec,b_current)
    !!write(lupri,*) 'iniguess projected'
    !!call mat_print(b_current,1,ndim,1,ndim,lupri)
    !      call normalize(b_current)
    !!write(lupri,*) 'iniguess projected, normalized'
    !!call mat_print(b_current,1,ndim,1,ndim,lupri)
    !      call mat_write_to_disk(DD%lub,b_current)
    !      call DD_lintrans(0.0E0_realk,scrvec,b_current,scrvec)
    !!write(lupri,*) 'iniguess projected, linearly transformed'
    !!call mat_print(scrvec,1,ndim,1,ndim,lupri)

    !      call mat_write_to_disk(DD%lusigma,scrvec)
    !      call DD_red_space(decomp,DD,iter,b_current,mu,x,res)
    !   enddo
    !!endif

    do i = 1, nb_curr !linear transformation of start vectors
       call DD_lintrans(DD,decomp,0.0E0_realk,b_current(i),sigma_cur(i),fifoqueue)
       call mat_write_to_disk(DD%lusigma,sigma_cur(i),OnMaster)
    enddo

    iter = 0
    do 
       iter = iter + 1
       if (decomp%info_stability) then
          call LSTIMER('START ', T1, T2, decomp%lupri) 
       endif
       if (iter == max_it) then
             WRITE(decomp%LUPRI,'(/A)') &
             &     'Max. number of Hessian iterations reached '
             write(decomp%lupri,*) 'Increase .STAB MAXIT in input. Currently set to', decomp%cfg_check_maxit
             CALL lsQUIT(' Max. number of Hessian iterations reached ',decomp%lupri)
       endif
       call DD_red_space(decomp,DD,ndim_red,n_on_disk,nb_curr,b_current,sigma_cur,mu,x,res)

       !We are now done with the initial guess. If number of start vectors is
       !larger than number of required solutions, we can free the extra allocated
       !matrices:
       if (iter == 1 .and. nstart > decomp%cfg_hessian_nvec) then
          do i = nstart, decomp%cfg_hessian_nvec+1, -1
             call mat_free(b_current(i))
             call mat_free(sigma_cur(i))
             deallocate(conv) ; allocate(conv(decomp%cfg_hessian_nvec))
          enddo
       endif

       !Convergence?
       conv = .false.
       do i = 1, decomp%cfg_hessian_nvec
          err(i) = SQRT(mat_sqnorm2(res(i)))
          if (err(i) < thresh) then
             conv(i) = .true.
          endif
          if (decomp%info_stability) then
                write (decomp%LUPRI, '("    Hessian eival: Residual for vector ", &
                & i3, ":  ", E18.10, " and eival: ",E18.10, "  It = ", i5, "  Conv = ", L2)') &
                & i, err(i), mu(i), iter, conv(i)
          endif
       enddo
       convergence = .true.
       do i = 1, decomp%cfg_hessian_nvec
          if (.not. conv(i)) then
             convergence = .false.
          endif
       enddo
       if (convergence) then
          if (decomp%info_stability) then
             write (decomp%lupri,*) '   Convergence of all vectors!'
             do j = 1, decomp%cfg_hessian_nvec
                write (decomp%lupri, *) '   Eigenvalue for vector ', j, ' is: ', mu(j)
             enddo
          endif
          exit
       endif
       !Preconditioning of residual = resP and b(iter+1) = resP:
       do i = 1, decomp%cfg_hessian_nvec !Precond each residual individually
          if (.not. conv(i)) then        !Only do precond if residual not converged
             call DD_precond(decomp, 0.0E0_realk,res(i),scrvec(i))
          endif
       enddo
       !Restart possibility for Hessian eigenvalue
       !   heseig_lun = -1
       !   call lsopen(heseig_lun,'heseig.restart','UNKNOWN','UNFORMATTED')
       !   rewind heseig_lun
       !   call mat_write_to_disk(heseig_lun,x)
       !   call lsclose(heseig_lun,'KEEP')
       !End restart 
       !Orthonormalize b_current against set of previous b vectors and save:
       nnew = decomp%cfg_hessian_nvec
       call DD_Hes_orthonormalize(decomp,scrvec,nnew,n_on_disk,conv,DD%lub,b_current)
       nb_curr = nnew !no of vecs in b_current
       n_on_disk = n_on_disk + nnew
       !if (linear_dep) then
       !   write (lupri,*) '   Linear dependencies, exiting'
       !   write (lupri,*) '   Current approximation to lowest eigenvalue, mu =', mu
       !   exit
       !endif
       !Calculate and save new sigma vector:
       do i = 1, nb_curr
          call DD_lintrans(DD,decomp,0.0E0_realk,b_current(i),sigma_cur(i),fifoqueue)
          call mat_write_to_disk(DD%lusigma,sigma_cur(i),OnMaster)
       enddo
       if (decomp%info_stability) then
          call LSTIMER('HES IT ', T1, T2, decomp%lupri)
       endif
    enddo

    CALL lsCLOSE(DD%lusigma,'DELETE')
    CALL lsCLOSE(DD%lub,'DELETE')
    deallocate(DD%Ared, DD%boverlaps)
    do i = 1, decomp%cfg_hessian_nvec
       call mat_free(b_current(i))
       call mat_free(sigma_cur(i))
       call mat_free(scrvec(i))
       call mat_free(res(i))
       call mat_free(x(i))
    enddo
    deallocate(b_current)
    deallocate(sigma_cur)
    deallocate(res)
    deallocate(x)
    deallocate(scrvec)
    deallocate(conv)
   end subroutine DD_solver

   !> \brief Solve in reduced space for Hessian eigenvalues.
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_red_space(decomp,DD,ndim_red,n_on_disk,nb_curr,b_current,sigma_cur,mu,x,res)
      implicit none
      !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
      type(decompItem),intent(in) :: decomp
      !> Contains settings concerning HOMO-LUMO gap and Hessian eigenvalues
      type(DDitem),intent(inout)  :: DD
      !> Dimension of red space - is updated
      integer, intent(inout)   :: ndim_red  
      !> Number of vectors on disk
      integer, intent(in)      :: n_on_disk
      !> Number of new b-vectors in b_current
      integer, intent(in)      :: nb_curr
      !> Contains new trial vector(s)
      type(matrix), intent(in) :: b_current(nb_curr)
      !> Contains new sigma vector(s)
      type(matrix), intent(in) :: sigma_cur(nb_curr)  
      !> Lowest eigenvalue(s) obtained in reduced space
      real(realk), intent(out) :: mu(:) 
      !> Full space solution vector x obtained from solution in reduced space (output)
      type(matrix),intent(inout) :: x(decomp%cfg_hessian_nvec)
      !> Full space residual (output)
      type(matrix),intent(inout) :: res(decomp%cfg_hessian_nvec)
      integer                  :: rowdim, coldim, i, j, IERR, k
      integer, allocatable, dimension(:)       :: IPIV
      type(matrix)                             :: sigma_n, bn, scrmat
      real(realk), allocatable, dimension(:)   :: eigenval, FV1, FV2, RHS, scrvec
      real(realk),allocatable, dimension(:,:)  :: A, eigenvec, xred
      real(realk)                              :: value
      logical :: OnMaster
      OnMaster = .TRUE.
   rowdim = b_current(1)%nrow
   coldim = b_current(1)%ncol
   allocate(A(n_on_disk,n_on_disk), eigenvec(n_on_disk,n_on_disk))
   allocate(eigenval(n_on_disk), FV1(n_on_disk), FV2(n_on_disk))
   allocate(RHS(n_on_disk), IPIV(n_on_disk), scrvec(n_on_disk))
   allocate(xred(n_on_disk,decomp%cfg_hessian_nvec))

   call mat_init(sigma_n,rowdim,coldim)
   call mat_init(bn,rowdim,coldim)
   call mat_init(scrmat,rowdim,coldim)

   !setup 'diagonal' elements
   do i = 1, nb_curr !loop over b
      do j = 1, nb_curr !loop over sigma
         value = mat_dotproduct(b_current(i),sigma_cur(j))
         DD%Ared(ndim_red+i,ndim_red+j) = value
         DD%boverlaps(ndim_red+i,ndim_red+j) = mat_dotproduct(b_current(i),b_current(j))
         DD%boverlaps(ndim_red+j,ndim_red+i) = DD%boverlaps(ndim_red+i,ndim_red+j)
      enddo
   enddo
   rewind(DD%lusigma)
   do j = 1, nb_curr
      do i = 1, n_on_disk  !Setup lower half of reduced matrix
         call mat_read_from_disk(DD%lusigma,sigma_n,OnMaster)
         value = mat_dotproduct(b_current(j),sigma_n)
         DD%Ared(ndim_red+j,i) = value
      enddo 
      rewind(DD%lusigma)
   enddo
   rewind(DD%lub)
   do j = 1, nb_curr 
      do i = 1, n_on_disk   !Explicitly calculate upper half:
         call mat_read_from_disk(DD%lub,bn,OnMaster)
         value = mat_dotproduct(bn,sigma_cur(j))
         DD%Ared(i,ndim_red+j) = value
         !Set up matrix of b vector overlaps (just for testing):
         DD%boverlaps(ndim_red+j,i) = mat_dotproduct(b_current(j),bn)
         DD%boverlaps(i,ndim_red+j) = DD%boverlaps(ndim_red+j,i)
      enddo
      rewind(DD%lub)
   enddo

   ndim_red = n_on_disk

   if (decomp%INFO_stability_REDSPACE) then
      write (decomp%lupri,*) 'b vector overlaps, DD_fock solver:'
      call LS_OUTPUT(DD%boverlaps, 1, ndim_red, 1, ndim_red, decomp%cfg_check_maxit, decomp%cfg_check_maxit, 1, decomp%lupri)

      write (decomp%lupri,*) 'Reduced E2, DD_fock solver:'
      call LS_OUTPUT(DD%Ared, 1, ndim_red, 1, ndim_red, decomp%cfg_check_maxit, decomp%cfg_check_maxit, 1, decomp%lupri)
   endif

   !Setup reduced matrix with proper dimensionguess
   A(1:ndim_red,1:ndim_red) = DD%Ared(1:ndim_red,1:ndim_red)

   !Solve eigenvalue problem A*X = mu*X in reduced space:
   call RS(ndim_red,ndim_red,A,eigenval,1,eigenvec,FV1,FV2,IERR)
   if (IERR /= 0) then
      write (decomp%LUPRI,*) 'Problem in RS (called from DD_Fock_red_space), IERR =', IERR
      STOP 'Problem in RS (called from DD_Fock_red_space)'
   endif
   do i = 1, decomp%cfg_hessian_nvec
      mu(i) = eigenval(i)
      !...and corresponding eigenvector:
      xred(1:ndim_red, i) = eigenvec(1:ndim_red,i)
   enddo
   !write (lupri,*) 'Solution vector, DD solver:'
   !call LS_OUTPUT(xred, 1, iter, 1, 1, iter, 1, 1, lupri)

   !Now obtain residual and x in real space by
   !res = Ax - mu*x = xred(1)*sigma1 + xred(2)*sigma2 + ...
   !                 - mu*xred(1)*b1  - mu*xred(2)*b2  - ...
   !x = xred(1)*b1 + xred(2)*b2  + ...
   do i = 1, decomp%cfg_hessian_nvec
      call mat_zero(res(i)) ; call mat_zero(x(i))
   enddo
   rewind(DD%lusigma) ; rewind(DD%lub)
   do j = 1, decomp%cfg_hessian_nvec
      do i = 1, ndim_red
         call mat_read_from_disk(DD%lusigma,sigma_n,OnMaster)
         call mat_read_from_disk(DD%lub,bn,OnMaster)

         call mat_add(xred(i,j), sigma_n, -xred(i,j)*mu(j), bn, scrmat)

         call mat_daxpy(1.0E0_realk, scrmat, res(j))
         call mat_daxpy(xred(i,j), bn, x(j))
      enddo
     !Don't rewind last time, file pointer should be at end of file or new sigmavec overwrites old ones!
      if (j /= decomp%cfg_hessian_nvec) then  
         rewind(DD%lusigma) ; rewind(DD%lub)
      endif
   enddo

   deallocate(A,eigenval,eigenvec,FV1,FV2,RHS,IPIV,scrvec,xred)
   call mat_free(sigma_n)
   call mat_free(bn)
   call mat_free(scrmat)
   end subroutine DD_red_space

   !> \brief Solve in reduced space for orbital eigenvalues.
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_Fock_red_space(DD,decomp,nroots,ndim,ndim_red,n_on_disk,nb_curr,b_current,sigma_cur,mu,x,res)
      implicit none
      !> Contains settings concerning HOMO-LUMO gap and Hessian eigenvalues
      type(DDitem),intent(inout) :: DD
      !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
      type(decompItem),intent(in) :: decomp
      !> Total number of solutions to be found
      integer, intent(in)      :: nroots 
      !> Number of basis functions
      integer, intent(in)      :: ndim
      !> Dimension of red space - is updated
      integer, intent(inout)   :: ndim_red  
      !> Number of vectors on disk
      integer, intent(in)      :: n_on_disk
      !> Number of new b-vectors in b_current
      integer, intent(in)      :: nb_curr
      !> Current trial vector(s)
      real(realk), intent(in)  :: b_current(ndim, nroots)
      !> Current linearly transformed trial vector(s)
      real(realk), intent(in)  :: sigma_cur(ndim, nroots)  
      !> Lowest eigenvalue(s) in reduced space
      real(realk), intent(out) :: mu(:) 
      !> Full space trial vector obtained from solution in reduced space
      real(realk), intent(out) :: x(ndim, nroots)
      !> Full space residual
      real(realk), intent(out) :: res(ndim, nroots)
      integer                  :: i, j, IERR, k
      integer, allocatable, dimension(:)      :: IPIV
      real(realk), allocatable, dimension(:)  :: sigma_n, bn, scrmat
      real(realk), allocatable, dimension(:)  :: eigenval, FV1, FV2, RHS, scrvec
      real(realk),allocatable, dimension(:,:) :: A, eigenvec, xred
      real(realk)                             :: value

   allocate(A(n_on_disk,n_on_disk), eigenvec(n_on_disk,n_on_disk))
   allocate(eigenval(n_on_disk), FV1(n_on_disk), FV2(n_on_disk))
   allocate(RHS(n_on_disk), IPIV(n_on_disk), scrvec(n_on_disk))
   allocate(xred(n_on_disk,nroots))

   allocate(sigma_n(ndim))
   allocate(bn(ndim))
   allocate(scrmat(ndim))

   !setup 'diagonal' elements
   do i = 1, nb_curr !loop over b
      do j = 1, nb_curr !loop over sigma
         value = dot_product(b_current(1:ndim,i),sigma_cur(1:ndim,j))
         DD%Ared(ndim_red+i,ndim_red+j) = value
         DD%boverlaps(ndim_red+i,ndim_red+j) = dot_product(b_current(1:ndim,i),b_current(1:ndim,j))
         DD%boverlaps(ndim_red+j,ndim_red+i) = DD%boverlaps(ndim_red+i,ndim_red+j)
      enddo
   enddo
   rewind(DD%lusigma)
   do j = 1, nb_curr
      do i = 1, n_on_disk  !Setup lower half of reduced matrix
         read(DD%lusigma) sigma_n
         value = dot_product(b_current(1:ndim,j),sigma_n)
         DD%Ared(ndim_red+j,i) = value
      enddo 
      rewind(DD%lusigma)
   enddo
   rewind(DD%lub)
   do j = 1, nb_curr 
      do i = 1, n_on_disk   !Explicitly calculate upper half:
         read(DD%lub) bn
         value = dot_product(bn,sigma_cur(1:ndim,j))
         DD%Ared(i,ndim_red+j) = value
         !Set up matrix of b vector overlaps (just for testing):
         DD%boverlaps(ndim_red+j,i) = dot_product(b_current(1:ndim,j),bn)
         DD%boverlaps(i,ndim_red+j) = DD%boverlaps(ndim_red+j,i)
      enddo
      rewind(DD%lub)
   enddo

   ndim_red = n_on_disk

   if (decomp%INFO_stability_redspace) then
      write (decomp%lupri,*) 'b vector overlaps, DD_fock solver:'
      call LS_OUTPUT(DD%boverlaps, 1, ndim_red, 1, ndim_red, decomp%cfg_homolumo_maxit, decomp%cfg_homolumo_maxit, 1, decomp%lupri)

      write (decomp%lupri,*) 'Reduced E2, DD_fock solver:'
      call LS_OUTPUT(DD%Ared, 1, ndim_red, 1, ndim_red, decomp%cfg_homolumo_maxit, decomp%cfg_homolumo_maxit, 1, decomp%lupri)
   endif

   !Setup reduced matrix with proper dimensionguess
   A(1:ndim_red,1:ndim_red) = DD%Ared(1:ndim_red,1:ndim_red)

   !Solve eigenvalue problem A*X = mu*X in reduced space:
   call RS(ndim_red,ndim_red,A,eigenval,1,eigenvec,FV1,FV2,IERR)
   if (IERR /= 0) then
      write (decomp%LUPRI,*) 'Problem in RS (called from DD_Fock_red_space), IERR =', IERR
      STOP 'Problem in RS (called from DD_Fock_red_space)'
   endif
   do i = 1, nroots
      mu(i) = eigenval(i)
      !...and corresponding eigenvector:
      xred(1:ndim_red, i) = eigenvec(1:ndim_red,i)
   enddo
   !Now obtain residual and x in real space by
   !res = Ax - mu*x = xred(1)*sigma1 + xred(2)*sigma2 + ...
   !                 - mu*xred(1)*b1  - mu*xred(2)*b2  - ...
   !x = xred(1)*b1 + xred(2)*b2  + ...
   !call mat_zero(res) ; call mat_zero(x)
   res = 0.0E0_realk ; x = 0.0E0_realk
   rewind(DD%lusigma) ; rewind(DD%lub)
   do j = 1, nroots
      do i = 1, ndim_red
         read(DD%lub) bn
         read(DD%lusigma) sigma_n
         scrmat = xred(i,j)*sigma_n - xred(i,j)*mu(j)*bn
         res(1:ndim,j) = res(1:ndim,j) + scrmat
         x(1:ndim,j)   = x(1:ndim,j) + xred(i,j)*bn
      enddo
     !Don't rewind last time, file pointer should be at end of file or new sigmavec overwrites old ones!
      if (j /= nroots) then  
         rewind(DD%lusigma) ; rewind(DD%lub)
      endif
   enddo

   deallocate(A,eigenval,eigenvec,FV1,FV2,RHS,IPIV,scrvec,xred)
   deallocate(sigma_n)
   deallocate(bn)
   deallocate(scrmat)
   end subroutine DD_Fock_red_space

   !> \brief Linear transformation for Hessian eigenvalue solver.
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_lintrans(DD,decomp, shift, x, MX, fifoqueue)
   implicit none
          !> Contains settings concerning HOMO-LUMO gap and Hessian eigenvalues
          type(DDitem),intent(inout) :: DD
          !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
          type(decompItem),intent(in) :: decomp
          !> Level shift/damping on linear transformation
          real(realk)              :: shift
          !> Trial vector to be linearly transformed
          type(matrix), intent(in)    :: x
          !> Linearly transformed trial vector (output)
          type(matrix), intent(inout) :: Mx
          !> Contains density and Fock matrices from previous iterations (for debugging only)
          type(modFIFO), intent(inout) :: fifoqueue
          type(solverItem)            :: arh
          integer                     :: ndens

          call arh_set_default_config(arh)
          if (DD%cfg_orbspread) then
             arh%orbspread_input => DD%orbspread_input
             arh%cfg_orbspread = .true.
          else if (DD%debug_arh_hessian) then
             ndens = fifoqueue%offset
             if (ndens > 0) then 
                allocate(arh%fifometric(ndens,ndens))
                allocate(arh%inv_fifometric(ndens,ndens))
                allocate(arh%fifoM(ndens,ndens))
                call fifo_inverse_metric(arh,fifoqueue)
                call arh_get_M(arh,fifoqueue)
             endif
          else
             arh%set_do_2nd_order = .true.
             arh%set_arhterms     = .false.
             arh%lupri            = decomp%lupri
          endif
          call arh_lintrans(arh,decomp,x,2,shift,MX,fifoqueue)

   end subroutine DD_lintrans

   !> \brief Linear transformation for orbital eigenvalue solver.
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_Fock_lintrans(ndim, M, x, MX)
   implicit none
          !> Number of basis functions
          integer, intent(in)      :: ndim
          !> Linearly transform by multiplying with this vector
          real(realk), intent(in)  :: M(ndim,ndim)
          !> Trial vector to be linearly transformed
          real(realk), intent(in)  :: x(ndim)
          !> Linearly transformed trial vector
          real(realk), intent(out) :: Mx(ndim)

         Mx = matmul(M,x)
   end subroutine DD_Fock_lintrans

   !> \brief Branches out to chosen preconditioning for Hessian eigenvalue(s).
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_precond(decomp, omega, x, xprec)
   implicit none
          !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
          type(decompItem),intent(inout) :: decomp
          !> Level shift
          real(realk), intent(in)    :: omega
          !> Vector to be preconditioned
          type(matrix), intent(in)   :: x
          !> Preconditioned input vector (output)
          type(matrix),intent(inout) :: xprec
          type(matrix)               :: scr
          type(solverItem)           :: arh
          type(modFIFO)              :: dummyfifo

   if (decomp%cfg_NOPREC) then
      xprec = x
   else if (decomp%cfg_orbspread) then
      arh%orbspread_input => decomp%orbspread_input
      arh%cfg_orbspread = .true.
      call arh_precond(arh,decomp,x,2,omega,xprec)
   else !if (cfg_do_2nd_order) then
      call mat_init(scr,x%nrow,x%ncol)

      call arh_set_default_config(arh)
      arh%set_arhterms     = .false.
      arh%lupri            = decomp%lupri

      decomp%pcg_preconditioning = .true.
      
      call project_oao_basis(decomp, x, 2, scr)
      call arh_PCG(arh, decomp, scr, xprec, omega, 2, dummyfifo)

      decomp%pcg_preconditioning = .false.
      call mat_free(scr)
   !else   
   !   call mat_init(scr,x%nrow,x%nrow)

   !   scr = x
   !   call mat_ao_precond(2,omega,decomp%FUP,decomp%FUQ,decomp%DU,scr)
   !   call project_oao_basis(decomp, scr, 2, xprec)

   !   call mat_free(scr)      
   endif

   end subroutine DD_precond

   !> \brief Diagonal preconditioning of trial vector/residual for Hessian eigenvalue(s).
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_diag_precond(decomp, omega, x, xprec)
   implicit none
          !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
          type(decompItem),intent(in) :: decomp
          !> Level shift
          real(realk), intent(in)    :: omega
          !> Vector to be preconditioned
          type(matrix), intent(in)   :: x
          !> Preconditioned input vector (output)
          type(matrix),intent(inout) :: xprec
          type(matrix)               :: scr

   if (decomp%cfg_NOPREC) then
      xprec = x
   else   
      call mat_init(scr,x%nrow,x%nrow)

      scr = x
      call mat_ao_precond(2,omega,decomp%FUP,decomp%FUQ,decomp%DU,scr)
      call project_oao_basis(decomp, scr, 2, xprec)

      call mat_free(scr)      
   endif

   end subroutine DD_diag_precond

   !> \brief Preconditioning for orbital eigenvalue solver.
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_Fock_precond(ndim, M, x, xprec)
   implicit none
          !> Number of basis functions
          integer, intent(in)      :: ndim
          !> Precondition trial vector by diagonal elements of this matrix
          real(realk), intent(in)  :: M(ndim,ndim)
          !> Trial vector/residual to be preconditioned
          real(realk), intent(in)  :: x(ndim)   
          !> Preconditioned trial vector
          real(realk), intent(out) :: xprec(ndim)
          integer                  :: i

      xprec = x
      do i = 1, ndim
         if (abs(M(i,i)) > 1.0E-9_realk) then 
            xprec(i) = xprec(i)/M(i,i)
         endif
      enddo

   end subroutine DD_Fock_precond

   !> \brief Project out redundant components for orbital eigenvalue solver.
   !> \author S. Host
   !> \date March 2005.
   subroutine DD_Fock_project(ndim,Proj,x,xproj)
   implicit none
          !> Number of basis functions
          integer, intent(in)      :: ndim
          !> Project trial vector with this matrix
          real(realk), intent(in)  :: Proj(ndim,ndim)
          !> Trial vector to be projected
          real(realk), intent(in)  :: x(ndim)
          !> Projected trial vector
          real(realk), intent(inout) :: xproj(ndim)

!          xproj = matmul(proj,x)
          call DGEMV('N',ndim,ndim,1.0E0_realk,Proj,ndim,X,1,0.0E0_realk,xproj,1)

   end subroutine DD_Fock_project

   !> \brief Orthonormalize new trial vector(s) for orbital eigenvalue solver.
   !> \author S. Host
   !> \date March 2005.
   !>
   !>   Orthogonalize new b-vectors against all previous b-vectors
   !>   and among themselves, and renormalize.
   !>   The b-vectors have the form \n
   !>         ( Y_dia    Z   )       ! Z_mu_nu  mu < nu ; Y_dia = Y_mu_mu \n
   !>         (  Y     Y_dia )       ! Y_mu_nu  mu > nu     \n
   !>   Each b-vector is in (Z, Y) form, the (Y, Z) vector is obtained
   !>   by transposing.  Each of the new b-vectors is
   !>   orthogonalized both against (Z, Y) and (Y, Z) for each previous
   !>   b-vectors. Orthogonalization is performed twice if round-off is large.
   !> 
  subroutine DD_orthonormalize(decomp,nroots,Bvec_tmp,ndim,Nb_new,Nb_prev,conv,P,lub,b_current)
    implicit none
    !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
    type(decompItem),intent(in) :: decomp
    !> Number of solutions to be found
    integer, intent(in)     :: nroots
    !> New non-orthogonal vectors
    real(realk),intent(inout) :: Bvec_tmp(:,:)
    !> Number of basis functions
    integer, intent(in)     :: ndim
    !> Input: number of new non-orthogonal b-vectors in Bvec_tmp. Output: Number of accepted new orthogonal b-vectors (linear deps. removed)
    integer, intent(inout)  :: Nb_new
    !> Number of previous (orthogonal) b-vectors located on file lub
    integer, intent(in)     :: Nb_prev
    !> If conv(i) = true, i'th root has converged and will not be modified
    logical, intent(in)     :: conv(nroots)
    !> Project out redundancies using this matrix
    real(realk), intent(in) :: P(ndim,ndim)
    !> Logical unit number for file containing previous (orthogonal) b-vectors
    integer, intent(in)     :: lub
    !> New accepted orthonormalized b-vectors
    real(realk), intent(out)  :: b_current(ndim,nb_new)
!local
    integer :: irx,irow,i,iround,iturn,k,ibvec,jbvec,jrx,iB,no_of_new
    integer,allocatable :: lin_depend(:) !linear dependency index array
    real(realk), allocatable :: B_scr(:)
    real(realk) :: TT,T1,T2, norm

    ALLOCATE(lin_depend(Nb_new))
    allocate(B_scr(ndim))  

! STEP 1: check initial linear dependencies among new vectors
!         if the vector is modified it is also projected!!!!!

    if (decomp%info_dd_normalize) then
       WRITE (decomp%LUPRI,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
       write(decomp%lupri,*) '=========== BEGIN DD_orthonormalize ==========='
       write(decomp%lupri,*) '****** STEP 1: Projection'
    endif
    do i = 1,Nb_new
      lin_depend(i) = 1  !initialize to non-dependent
      if (.not. conv(i)) then
         call DD_Fock_project(ndim,P,bvec_tmp(1:ndim,i),b_scr)
         bvec_tmp(1:ndim,i) = b_scr
         call vecnorm(bvec_tmp(1:ndim,i),ndim,norm)
         if (decomp%info_dd_normalize) write(decomp%lupri,*) 'Norm of vector', i, 'after projection:', norm
      endif
    end do
!---
    iround = 0
    iturn = 0
!loop start
 1500 iturn = iturn + 1
!
! Orthogonalize new b-vectors agains previous (old) b-vectors
! (both (Z, Y) and (Y, Z))
!
    if (decomp%info_dd_normalize) then
      write(decomp%lupri,*) '****** STEP 2: orthogonalize new b-vectors agains previous b'
      write(decomp%lupri,*) 'NB_PREV, NB_NEW ', NB_PREV, NB_NEW
      if (nb_prev == 0) write (decomp%lupri,*) 'No old vectors, skip STEP 2'
    endif

    rewind(lub)
    do k = 1,Nb_prev
      read(lub) b_scr  !b_scr is bvecs(k) 
      do irx = 1,Nb_new
        if (.not. conv(irx)) then
           !if lin_depend(irx) == 0, the vector is skipped
           !because of linear dependencies
           if (lin_depend(irx) /= 0) then
             !Orthogonalize to (Z Y)_old
             if (decomp%info_dd_normalize) write(decomp%lupri,*) 'prev index, new index', k, irx
             TT = dot_product(b_scr,Bvec_tmp(1:ndim,irx))
             Bvec_tmp(1:ndim,irx) = Bvec_tmp(1:ndim,irx) - TT*b_scr
           endif
        endif
      enddo
    enddo

    if (decomp%info_dd_normalize) write(decomp%lupri,*) '****** STEP 3: &
       & Orthogonalize new b-vectors against each other '
!
! Orthogonalize new vectors against each other
!
    do ibvec = 1,Nb_new !index for current bvector
      if (.not. conv(ibvec)) then
         if (decomp%info_dd_normalize) write(decomp%lupri,*)'Current ibvec index is ', ibvec, lin_depend(ibvec)
   
         if (lin_depend(ibvec) /= 0) then
            jbvec = 1    !index for another bvector
            do jrx = 1,(ibvec-1)
               if (lin_depend(jrx) == 0) then
                  jbvec = jbvec + 1 !skip
               else
                  call VECNORM(Bvec_tmp(1:ndim,jbvec),ndim,T1)
                  if (decomp%info_dd_normalize) write(decomp%lupri,*) 'Norm of jbvec', jbvec, T1
                  T2 = dot_product(Bvec_tmp(1:ndim,jbvec),Bvec_tmp(1:ndim,ibvec))
                  if (decomp%info_dd_normalize) write(decomp%lupri,*) '<bjbvec|bibvec>', T2
                  TT = -T2/T1
                  if (decomp%info_dd_normalize) write(decomp%lupri,*) 'TT = -T2/T1', TT
                  Bvec_tmp(1:ndim,ibvec) = Bvec_tmp(1:ndim,ibvec) + TT*Bvec_tmp(1:ndim,jbvec)
                  jbvec = jbvec + 1
               endif
            enddo
         endif
!
! NORMALIZE VECTOR NUMBER Ibvec
!
         if (lin_depend(ibvec) /= 0) then !not linear dependent
           if (decomp%info_dd_normalize) write(decomp%lupri,*) 'DD_FOCK_SOLVER: Orthonormalize -> normalize'
           call DD_normalize(decomp,ndim,iturn,ibvec,iround,lin_depend(ibvec),Bvec_tmp(1:ndim,ibvec))
         endif
       endif    
    enddo !ibvec

!WHAT IS THIS ALL ABOUT?? could it be more clear what actually happens?
! Get rid of the GO TO stuff!!
    if (iround > 0) then
      iround = 0
      if (iturn == 1) GO TO 1500
      STOP 'orthonormalize sw-error; iround > 0 and iturn > 1'
    endif
!
! Add new vectors to file !!!!!!!!!!!!!!!!!!!!!!!!! (in memory together with old ones)
! ib counts how many "acceptable" vectors there are in total
!    
    no_of_new = 0
    iB = Nb_prev
    do irx = 1,Nb_new
      if (.not. conv(irx)) then
         if (lin_depend(irx) /= 0) then
           iB = iB + 1
           no_of_new = no_of_new + 1
           b_current(1:ndim,no_of_new) = Bvec_tmp(1:ndim,irx)
           write(lub) Bvec_tmp(1:ndim,irx)
           if (decomp%info_dd_normalize) write(decomp%lupri,*) 'DD_orthonormalize: new vector is added'
         endif
      endif
    enddo
!
!     Set NB_NEW to actual number of acceptable new trial vectors
!
    Nb_new = iB - Nb_prev

    if (decomp%info_dd_normalize) then
       write(decomp%lupri,*) '============ END DD_orthonormalize ============'
       write (decomp%lupri,*) 'Number of new vectors accepted and written to disk: ', nb_new
    endif

    DEALLOCATE(lin_depend)
    deallocate(B_scr)

  end subroutine DD_orthonormalize

   !> \brief Orthonormalize new trial vector(s) for Hessian eigenvalue solver.
   !> \author S. Host
   !> \date March 2005.
   !>
   !>   Orthogonalize new b-vectors against all previous b-vectors
   !>   and among themselves, and renormalize.
   !>   The b-vectors have the form \n
   !>         ( Y_dia    Z   )       ! Z_mu_nu  mu < nu ; Y_dia = Y_mu_mu \n
   !>         (  Y     Y_dia )       ! Y_mu_nu  mu > nu     \n
   !>   Each b-vector is in (Z, Y) form, the (Y, Z) vector is obtained
   !>   by transposing.  Each of the new b-vectors is
   !>   orthogonalized both against (Z, Y) and (Y, Z) for each previous
   !>   b-vectors. Orthogonalization is performed twice if round-off is large.
   !> 
  subroutine DD_Hes_orthonormalize(decomp,Bvec_tmp,Nb_new,Nb_prev,conv,lub,b_current)
    implicit none
    !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
    type(decompItem),intent(in) :: decomp
    !> New non-orthogonal vectors
    type(matrix),intent(inout) :: Bvec_tmp(Nb_new)
    !> Input: number of new non-orthogonal b-vectors in Bvec_tmp. Output: Number of accepted new orthogonal b-vectors (linear deps. removed)
    integer, intent(inout)  :: Nb_new
    !> Number of previous (orthogonal) b-vectors located on file lub
    integer, intent(in)     :: Nb_prev
    !> If conv(i) = true, i'th root has converged and will not be modified
    logical, intent(in)     :: conv(Nb_new)
    !> Logical unit number for file containing previous (orthogonal) b-vectors
    integer, intent(in)     :: lub
    !> New accepted orthonormalized b-vectors
    type(matrix), intent(inout) :: b_current(Nb_new) !output
!local
    integer :: rowdim,coldim,irx,irow,i,iround,iturn,k,ibvec,jbvec,jrx,iB,no_of_new
    integer,allocatable :: lin_depend(:) !linear dependency index array
    type(matrix) :: B_scr
    real(realk) :: TT,T1,T2, norm
    logical :: OnMaster
    OnMaster = .true.
    rowdim = Bvec_tmp(1)%nrow
    coldim = Bvec_tmp(1)%ncol
    ALLOCATE(lin_depend(Nb_new))
    call mat_init(B_scr,rowdim,coldim)  

! STEP 1: check initial linear dependencies among new vectors
!         if the vector is modified it is also projected!!!!!

    if (decomp%info_dd_normalize) then
       WRITE (decomp%LUPRI,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
       write(decomp%lupri,*) '=========== BEGIN DD_orthonormalize ==========='
       write(decomp%lupri,*) '****** STEP 1: Projection'
    endif
    do i = 1,Nb_new
      if (decomp%info_dd_normalize) write(decomp%lupri,*) 'Vector ', i, ': conv =', conv(i)
      if (.not. conv(i)) then
         lin_depend(i) = 1  !initialize to non-dependent
         call project_oao_basis(decomp, bvec_tmp(i), 2, b_scr)
         call mat_assign(bvec_tmp(i),b_scr)
         norm = sqrt(mat_sqnorm2(bvec_tmp(i)))
         if (decomp%info_dd_normalize) write(decomp%lupri,*) 'Norm of vector', i, 'after projection:', norm
      endif
    end do
!---
    iround = 0
    iturn = 0
!loop start
 1500 iturn = iturn + 1
!
! Orthogonalize new b-vectors agains previous (old) b-vectors
! (both (Z, Y) and (Y, Z))
!
    if (decomp%info_dd_normalize) then
      write(decomp%lupri,*) '****** STEP 2: orthogonalize new b-vectors agains previous b'
      write(decomp%lupri,*) 'NB_PREV, NB_NEW ', NB_PREV, NB_NEW
      if (nb_prev == 0) write (decomp%lupri,*) 'No old vectors, skip STEP 2'
    endif

    rewind(lub)
    do k = 1,Nb_prev
      call mat_read_from_disk(lub,b_scr,OnMaster)
      do irx = 1,Nb_new
        if (.not. conv(irx)) then
           !if lin_depend(irx) == 0, the vector is skipped
           !because of linear dependencies
           if (lin_depend(irx) /= 0) then
             !Orthogonalize to (Z Y)_old
             if (decomp%info_dd_normalize) write(decomp%lupri,*) 'prev index, new index', k, irx
             TT = mat_dotproduct(b_scr,Bvec_tmp(irx))
             call mat_daxpy(-TT,b_scr,Bvec_tmp(irx))
           endif
        endif
      enddo
    enddo

    if (decomp%info_dd_normalize) write(decomp%lupri,*) '****** STEP 3: &
       & Orthogonalize new b-vectors against each other '
!
! Orthogonalize new vectors against each other
!
    do ibvec = 1,Nb_new !index for current bvector
      if (.not. conv(ibvec)) then
         if (decomp%info_dd_normalize) write(decomp%lupri,*)'Current ibvec index is ', ibvec, lin_depend(ibvec)
   
         if (lin_depend(ibvec) /= 0) then
            jbvec = 1    !index for another bvector
            do jrx = 1,(ibvec-1)
               if (lin_depend(jrx) == 0) then
                  jbvec = jbvec + 1 !skip
               else
                  T1 = sqrt(mat_sqnorm2(Bvec_tmp(jbvec)))
                  if (decomp%info_dd_normalize) write(decomp%lupri,*) 'Norm of jbvec', jbvec, T1
                  T2 = mat_dotproduct(Bvec_tmp(jbvec),Bvec_tmp(ibvec))
                  if (decomp%info_dd_normalize) write(decomp%lupri,*) '<bjbvec|bibvec>', T2
                  TT = -T2/T1
                  if (decomp%info_dd_normalize) write(decomp%lupri,*) 'TT = -T2/T1', TT
                  call mat_daxpy(TT,Bvec_tmp(jbvec),Bvec_tmp(ibvec))
                  jbvec = jbvec + 1
               endif
            enddo
         endif
!
! NORMALIZE VECTOR NUMBER Ibvec
!
         if (lin_depend(ibvec) /= 0) then !not linear dependent
           if (decomp%info_dd_normalize) write(decomp%lupri,*) 'DD_SOLVER: Orthonormalize -> normalize'
           call DD_Hes_normalize(decomp,iturn,ibvec,iround,lin_depend(ibvec),Bvec_tmp(ibvec))
         endif
       endif    
    enddo !ibvec

!WHAT IS THIS ALL ABOUT?? could it be more clear what actually happens?
! Get rid of the GO TO stuff!!
    if (iround > 0) then
      iround = 0
      if (iturn == 1) GO TO 1500
      STOP 'orthonormalize sw-error; iround > 0 and iturn > 1'
    endif
!
! Add new vectors to file !!!!!!!!!!!!!!!!!!!!!!!!! (in memory together with old ones)
! ib counts how many "acceptable" vectors there are in total
!    
    no_of_new = 0
    iB = Nb_prev
    do irx = 1,Nb_new
      if (.not. conv(irx)) then
         if (lin_depend(irx) /= 0) then
           iB = iB + 1
           no_of_new = no_of_new + 1
           call mat_assign(b_current(no_of_new),Bvec_tmp(irx))
           call mat_write_to_disk(lub,Bvec_tmp(irx),OnMaster)
           if (decomp%info_dd_normalize) write(decomp%lupri,*) 'DD_orthonormalize: new vector is added'
         endif
      endif
    enddo
!
!     Set NB_NEW to actual number of acceptable new trial vectors
!
    Nb_new = iB - Nb_prev

    if (decomp%info_dd_normalize) then
       write(decomp%lupri,*) '============ END DD_orthonormalize ============'
       write (decomp%lupri,*) 'Number of new vectors accepted and written to disk: ', nb_new
    endif

    DEALLOCATE(lin_depend)
    call mat_free(B_scr)

  end subroutine DD_Hes_orthonormalize

  !> \brief Normalize new trial vector(s) for orbital eigenvalue solver.
  !> \author S. Host
  !> \date March 2005.
  subroutine DD_normalize(decomp,ndim,iturn,ibvec,iround,lin_depend_i,Bvec_tmp_i)
    implicit none
    !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
    type(decompItem),intent(in) :: decomp
    !> Number of basis functions
    integer, intent(in) :: ndim
    !> Keep track of how many times normalize has been called
    integer, intent(in) :: iturn
    !> Number of vector to be normalized (only for printout)
    integer, intent(in) :: ibvec
    !> Keep track of how many times normalize has been called
    integer, intent(inout) :: iround
    !> Set to zero if vector is linearly dependent
    integer, intent(inout) :: lin_depend_i
    !> Trial vector to be normalized
    real(realk), intent(inout) :: Bvec_tmp_i(ndim)
    real(realk) :: TT, norm

    !norm of current vector
    call vecnorm(Bvec_tmp_i,ndim,TT)
    if (decomp%info_dd_normalize) write(decomp%lupri,*) 'norm of vector BEFORE', TT
    if (TT < 1.0E-12_realk) then
      !first type of linear dependence
      if (decomp%info_dd_normalize) then
        write(decomp%lupri,*) 'Vector ', ibvec, 'corresponding to root nr.', lin_depend_i, &
        ' is linearly dependent, TT=', TT, 'and is removed from the chosen set'
      endif
      lin_depend_i = 0  !linear dependent
    elseif(TT < 1E-5_realk) then
      if (iturn == 1) then
        iround = iround + 1
      else
        !Second type of linear dependence
        if (decomp%info_dd_normalize) then
          write(decomp%lupri,*)'Vector ',ibvec,' is linearly dependent, TT=', TT
          write(decomp%lupri,*)'too small norm^2 in second Gram-Schmidt orthonorm.'
        endif
        lin_depend_i = 0
      endif
    END IF
    IF (lin_depend_i /= 0) THEN
      IF (TT < 1E-8_realk) THEN
         TT = 1E0_realk / TT     
         BVEC_tmp_i = TT*BVEC_tmp_i
         call vecnorm(Bvec_tmp_i,ndim,TT)
      END IF
      TT = 1E0_realk / TT
      BVEC_tmp_i = TT*BVEC_tmp_i
      if (decomp%info_dd_normalize) then
        call vecnorm(Bvec_tmp_i,ndim,norm)
        if (decomp%info_dd_normalize) write(decomp%lupri,*) 'norm of vector AFTER', norm  
      endif
    END IF
  end subroutine DD_normalize

  !> \brief Normalize new trial vector(s) for Hessian eigenvalue solver.
  !> \author S. Host
  !> \date March 2005.
  subroutine DD_Hes_normalize(decomp,iturn,ibvec,iround,lin_depend_i,Bvec_tmp_i)
    implicit none
    !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
    type(decompItem),intent(in) :: decomp
    !> Keep track of how many times normalize has been called
    integer, intent(in) :: iturn
    !> Number of vector to be normalized (only for printout)
    integer, intent(in) :: ibvec
    !> Keep track of how many times normalize has been called
    integer, intent(inout) :: iround
    !> Set to zero if vector is linearly dependent
    integer, intent(inout) :: lin_depend_i
    !> Trial vector to be normalized
    type(matrix), intent(inout) :: Bvec_tmp_i
    real(realk) :: TT, norm

    !norm of current vector
    TT = sqrt(mat_sqnorm2(Bvec_tmp_i))
    if (decomp%info_dd_normalize) write(decomp%lupri,*) 'norm of vector BEFORE', TT
    if (TT < 1.0E-12_realk) then
      !first type of linear dependence
      if (decomp%info_dd_normalize) then
        write(decomp%lupri,*) 'Vector ', ibvec, 'corresponding to root nr.', lin_depend_i, &
        ' is linearly dependent, TT=', TT, 'and is removed from the chosen set'
      endif
      lin_depend_i = 0  !linear dependent
    elseif(TT < 1E-5_realk) then
      if (iturn == 1) then
        iround = iround + 1
      else
        !Second type of linear dependence
        if (decomp%info_dd_normalize) then
          write(decomp%lupri,*)'Vector ',ibvec,' is linearly dependent, TT=', TT
          write(decomp%lupri,*)'too small norm^2 in second Gram-Schmidt orthonorm.'
        endif
        lin_depend_i = 0
      endif
    END IF
    IF (lin_depend_i /= 0) THEN
      IF (TT < 1E-8_realk) THEN
         TT = 1E0_realk / TT     
         call mat_scal(TT,BVEC_tmp_i)
         TT = sqrt(mat_sqnorm2(Bvec_tmp_i))
      END IF
      TT = 1E0_realk / TT
      call mat_scal(TT,BVEC_tmp_i)
      if (decomp%info_dd_normalize) then
        norm = sqrt(mat_sqnorm2(Bvec_tmp_i))
        if (decomp%info_dd_normalize) write(decomp%lupri,*) 'norm of vector AFTER', norm  
      endif
    END IF

   end subroutine DD_Hes_normalize

   !> \brief Diagonalize FUQ and FUP to get all virtual and occupied orbital eigenvalues.  
   !> \author S. Host
   !> \date March 2005 (Rewritten to calculate one eigenvalue at a time (and not eigenvectors) May 2012 SR)
   subroutine dd_debug_homolumo(decomp,diag_hlgap)
   implicit none
        !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
        type(decompItem),intent(in) :: decomp
        !> HOMO-LUMO gap obtained by diagonalization
        real(realk),intent(out)  :: diag_hlgap
        integer                  :: ndim
        type(matrix)             :: Mat
        real(realk)              :: HOMO, LUMO
        integer                  :: nocc

   if (decomp%cfg_unres) then
     call dd_debug_homolumo_unres(decomp)
     diag_hlgap = 0.0E0_realk
     write(decomp%lupri,*) 'WARNING: HOMO-LUMO gap not calculated for unrestricted!'
   else
     ndim = decomp%FUP%nrow
     nocc = decomp%nocc
     call mat_init(Mat,ndim,ndim)

!    If the homo-orbital is negative it must be the eigenvalue with index nocc
     call mat_assign(Mat,decomp%FUP)
     call mat_dsyevx(Mat,HOMO,nocc)

!    If the homo orbital is positive the eigenvalue with index nocc is zero
     IF (abs(HOMO) .LE. 1.0E-7_realk) then
!      If the homo orbital is positive it must be the last eigenvalue
       call mat_assign(Mat,decomp%FUP)
       call mat_dsyevx(Mat,HOMO,ndim)
     ENDIF

!    If all virtual orbitals are positive the lumo is the eigenvalue with index nocc+1
     call mat_assign(Mat,decomp%FUQ)
     call mat_dsyevx(Mat,LUMO,nocc+1)

!    If the lumo orbital is negative the eigenvalue with index nocc is zero
     IF (abs(LUMO) .LE. 1.0E-7_realk) then
!      If the lumo orbital is negative it is the eigenvalue with index 1
       call mat_assign(Mat,decomp%FUQ)
       call mat_dsyevx(Mat,LUMO,1)
     ENDIF

     call mat_free(Mat)

     write (decomp%lupri,*)
     write (decomp%lupri,'("    E(LUMO):                  ", F12.6, " au")') LUMO
     write (decomp%lupri,'("   -E(HOMO):                  ", F12.6, " au")') HOMO
     write (decomp%lupri,'("   -------------------------------------")')
     write (decomp%lupri,'("    HOMO-LUMO Gap (by diag.): ", F12.6, " au")') LUMO-HOMO
     write (decomp%lupri,*)
     
     diag_hlgap = LUMO-HOMO

   endif

   end subroutine dd_debug_homolumo

   !> \brief Diagonalize FUQ and FUP to get all virtual and occupied orbital eigenvalues (for unrestricted).  
   !> \author S. Host
   !> \date March 2005
   subroutine dd_debug_homolumo_unres(decomp)
   implicit none
        !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
        type(decompItem),intent(in) :: decomp
        !real(realk),intent(out)  :: diag_hlgap
        integer                  :: ndim, fulldim
        real(realk), allocatable :: fullmat(:,:), FUQfull(:,:), FUPfull(:,:), eigenval(:), eigenvec(:,:)
        real(realk), allocatable :: temp(:)
        integer, allocatable     :: itemp(:), IFAIL(:)
        real(realk)              :: VL, VU !, HOMO, LUMO
        integer                  :: IL, IU, neig, ltemp, INFO, m 
   
   INFO=0
   ndim = decomp%FUP%nrow
   fulldim = ndim*2
   VL = 0.0E0_realk
   VU = 0.0E0_realk
   Ltemp = 8*ndim 
   allocate(fullmat(fulldim,fulldim))
   allocate(FUQfull(ndim,ndim),FUPfull(ndim,ndim),eigenval(ndim))
   allocate(eigenvec(ndim,ndim),temp(Ltemp),Itemp(5*ndim),IFAIL(ndim))

   !1. Do occupied orbitals:
   call mat_to_full(decomp%FUP,1.0E0_realk,fullmat)
      !1a. Do alpha occupied eigenvalues:
      FUPfull(1:ndim,1:ndim) = fullmat(1:ndim,1:ndim)

      call DSYEVX('V', 'A', 'U', ndim, FUPfull, ndim, VL, VU, IL, IU, &
        &  0.0E0_realk, neig, eigenval, eigenvec, ndim, temp, Ltemp, Itemp, &
        &  IFAIL, INFO )

      if (info /= 0) STOP 'Problem in DSYEVX (dd_debug_homolumo_unres)'

      write(decomp%lupri,*) "Occupied alpha orbital energies:"
      do m = 1, neig
         write(decomp%lupri,'(i6,F15.7)') m, eigenval(m)
      enddo
      !do m = neig, 1, -1
      !   if (abs(eigenval(m)) > 1.0E-7_realk) then
      !      HOMO = eigenval(m)
      !      exit
      !   endif
      !enddo
      !write(lupri,*) 'Highest alpha occupied eigenvalue:', HOMO

      !1b. Do beta occupied eigenvalues:
      FUPfull(1:ndim,1:ndim) = fullmat(ndim+1:fulldim,ndim+1:fulldim)

      call DSYEVX('V', 'A', 'U', ndim, FUPfull, ndim, VL, VU, IL, IU, &
        &  0.0E0_realk, neig, eigenval, eigenvec, ndim, temp, Ltemp, Itemp, &
        &  IFAIL, INFO )

      if (info /= 0) STOP 'Problem in DSYEVX (dd_debug_homolumo_unres)'

      write(decomp%lupri,*) "Occupied beta orbital energies:"
      do m = 1, neig
         write(decomp%lupri,'(i6,F15.7)') m, eigenval(m)
      enddo
      !do m = neig, 1, -1
      !   if (abs(eigenval(m)) > 1.0E-7_realk) then
      !      HOMO = eigenval(m)
      !      exit
      !   endif
      !enddo
      !write(lupri,*) 'Highest beta occupied eigenvalue:', HOMO

   !2. Do virtual orbitals:
   call mat_to_full(decomp%FUQ,1.0E0_realk,fullmat)
      !2a. Do alpha virtual eigenvalues:
      FUQfull(1:ndim,1:ndim) = fullmat(1:ndim,1:ndim)

      call DSYEVX('V', 'A', 'U', ndim, FUQfull, ndim, VL, VU, IL, IU, &
        &                   0.0E0_realk, neig, eigenval, eigenvec, ndim, temp, Ltemp, Itemp, &
        &                   IFAIL, INFO )

      if (info /= 0) STOP 'Problem in DSYEVX (dd_debug_homolumo_unres)'

      !write(decomp%lupri,*) "Dim of FUQ and eigenvalues of FUQ, # found:", ndim, neig
      write(decomp%lupri,*) "Virtual alpha orbital energies:"
      do m = 1, neig
         write(decomp%lupri,'(i6,F15.7)') m, eigenval(m)
      enddo
      !do m = 1, neig
      !   if (abs(eigenval(m)) > 1.0E-7_realk) then
      !      LUMO = eigenval(m)
      !      exit
      !   endif
      !enddo
      !write(decomp%lupri,*) 'Lowest virtual eigenvalue:', LUMO

      !2b. Do beta virtual eigenvalues:
      FUQfull(1:ndim,1:ndim) = fullmat(ndim+1:fulldim,ndim+1:fulldim)

      call DSYEVX('V', 'A', 'U', ndim, FUQfull, ndim, VL, VU, IL, IU, &
        &                   0.0E0_realk, neig, eigenval, eigenvec, ndim, temp, Ltemp, Itemp, &
        &                   IFAIL, INFO )

      if (info /= 0) STOP 'Problem in DSYEVX (dd_debug_homolumo_unres)'

      !write(decomp%lupri,*) "Dim of FUQ and eigenvalues of FUQ, # found:", ndim, neig
      write(decomp%lupri,*) "Virtual beta orbital energies:"
      do m = 1, neig
         write(decomp%lupri,'(i6,F15.7)') m, eigenval(m)
      enddo
      !do m = 1, neig
      !   if (abs(eigenval(m)) > 1.0E-7_realk) then
      !      LUMO = eigenval(m)
      !      exit
      !   endif
      !enddo
      !write(decomp%lupri,*) 'Lowest virtual eigenvalue:', LUMO

      !write(decomp%lupri,*) 'HOMO-LUMO gap:', LUMO-HOMO

      !write (decomp%lupri,*)
      !write (decomp%lupri,'("    E(LUMO):                  ", F12.6, " au")') LUMO
      !write (decomp%lupri,'("   -E(HOMO):                  ", F12.6, " au")') HOMO
      !write (decomp%lupri,'("   -------------------------------------")')
      !write (decomp%lupri,'("    HOMO-LUMO Gap (by diag.): ", F12.6, " au")') LUMO-HOMO
      !write (decomp%lupri,*)
   
   !diag_hlgap = LUMO-HOMO
   deallocate(FUPfull)
   deallocate(FUQfull)
   deallocate(eigenval)
   deallocate(eigenvec)
   deallocate(temp)
   deallocate(Itemp)
   deallocate(IFAIL)

   end subroutine dd_debug_homolumo_unres

end module direct_dens_util
