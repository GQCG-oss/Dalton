!> @file
!> Contains Direct Density Optimization utilities module for unrestricted calculations.
module direct_dens_util_unres
use precision
use matrix_module
use matrix_operations
use direct_dens_util
use matrix_operations_unres_dense
use queue_module
use rspPrecond
!use queue_ops
contains


!> \brief Main driver for getting HOMO-LUMO gap and, if requested, lowest Hessian eigenvlaue(s). 
!> \author S. Host
!> \date 2005
!>
!> ==================================================================================== \n
!>  For unrestricted SCF!!! 
!>  It is really not pretty with these 'not behind the curtain routines' but 
!>  necessary since mat operations like mat_to_full cannot be behind the curtain
!>  because of the difference in dimensions. \n
!>  \n
!>  Calculate final HOMO-LUMO gap, and, if requested, lowest Hessian eigenvalue \n
!>  \n
!>  The homo-lumo gap is inexpensive and is always calculated. From the orbital
!   energies, we obtain starting guesses for both excitation energies and Hessian eigenvalues.
!>  MAKE SURE THAT THE PROPER DU, FU, FUP, and FUQ ARE CREATED BEFORE CALLING THIS ROUTINE
!>  (call to oao_get_transformed_matrices). \n 
!> ====================================================================================
!>
   subroutine DD_homolumo_and_heseigen_unres(DD,decomp,debug,do_rsp_iniguess,howmany,rsp_iniguess)
   use rspPrecond
   implicit none
          !> Contains settings for HOMO-LUMO gap and Hessian eigenvalues
          type(DDitem),intent(inout)  :: DD
          !> Contains settings for decomposition and OAO decomposed overlap
          type(decompItem),intent(inout) :: decomp
          !> Contains debug info
          type(debugItem),intent(inout)  :: debug
          !> True if we are calculating response initial guess(es)
          logical, intent(in)      :: do_rsp_iniguess
          !> Number of rsp iniguesses to be calculated. Not referenced if do_rsp_iniguess=.false.
          integer,intent(in)       :: howmany
          !> Rsp initial guesses (output)
          type(matrix),intent(inout),optional :: rsp_iniguess(howmany)
          type(matrix),allocatable :: iniguess(:)
          real(realk), allocatable :: iniguess_full(:,:), alfa_orb_energy_diff(:), beta_orb_energy_diff(:) 
          integer, allocatable     :: alfaFUQFUPindex(:,:), betaFUQFUPindex(:,:)
          real(realk)              :: fac, hessian_eigenval(decomp%cfg_hessian_nvec), HOMOLUMOgap_a, HOMOLUMOgap_b, TOTALgap, &
                                    & T1, T2, norm, alfamin, betamin !, cpu1, cpu2, wall1, wall2
          integer                  :: mu, nu, matdim, i, k, j, posmin(1), posmax(1), fulldim, nguesses
          real(realk), allocatable :: fullmat1(:,:), fullmat2(:,:), halfmat(:,:), proj(:,:)
          logical                  :: alfa, NoAlfa, NoBeta, alfaocc_conv, alfavirt_conv, betaocc_conv, betavirt_conv
          type(matrix)             :: scr, scr2
          real(realk),external     :: frob_norm
          type(modFIFO)            :: dummyfifo
          integer                  :: KK1,KK2
          real(realk)              :: TMP1,TMP2

   alfaocc_conv  = .false.
   alfavirt_conv = .false.
   betaocc_conv  = .false.
   betavirt_conv = .false.

   DD%debug_arh_hessian = .false.

   if (decomp%cfg_check_converged_solution .and. do_rsp_iniguess) then
      call lsquit('Currently not possible to do both Hessian eigenvalues and excitation energies',decomp%lupri)
   endif

   if (do_rsp_iniguess) then
      if (.not. present(rsp_iniguess)) then
         call lsquit('rsp_iniguess must be present when do_rsp_iniguess=.true.',decomp%lupri)
      endif
   endif

   NoAlfa = .false.
   NoBeta = .false.

   matdim = decomp%FUP%nrow
   fulldim = 2*matdim
   !finding_arh_eigenvalue = .false.
   allocate(fullmat1(fulldim,fulldim),fullmat2(fulldim,fulldim))
   allocate(halfmat(matdim,matdim),proj(matdim,matdim))

   call time_II_operations1
   CALL LSTIMER('START ',T1,T2, decomp%lupri)
   !call gettim(cpu1,wall1)
   if (decomp%info_stability) then
      write(decomp%lupri,*)
      write (decomp%lupri,'("Calculation of orbital energies")')
      write (decomp%lupri,'("===============================")')
   endif

!1. Find eigenvalues of alfa part of FUP
   if (decomp%info_stability) then 
      write (decomp%lupri,*)
      write (decomp%lupri,'("Alfa occupied orbitals:")')
      write (decomp%lupri,'("=======================")')
   endif
   fullmat1 = 0.0E0_realk ; fullmat2 = 0.0E0_realk
   call mat_to_full(decomp%FUP,1.0E0_realk,fullmat1)
   call mat_to_full(decomp%DU,1.0E0_realk,fullmat2)
   !Set halfmat to alfa part of FUP:
   halfmat(1:matdim,1:matdim) = fullmat1(1:matdim,1:matdim)
   !Set proj to alfa part of DUone:
   proj(1:matdim,1:matdim) = fullmat2(1:matdim,1:matdim)
       !write(lupri,*) 'FUP halfmat, alfa part **'
       !call LS_OUTPUT(halfmat, 1, matdim, 1, matdim, matdim, matdim, 1, lupri)
   !norm = frob_norm(halfmat,matdim,matdim)
   !if (norm < 1.0E-12_realk) then
   if (DD%nocca == 0) then
       NoAlfa = .true.
       write (decomp%lupri,'("There appears to be no alfa occupied orbitals, continuing...")')
   else
      call DD_Fock_eigenvalue_unres('a',DD,decomp,halfmat,proj,'h',0.0E0_realk)
      if (decomp%cfg_hlgap_converged) then
         alfaocc_conv = .true.
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of alpha occupied orbital energies converged in", &
           & i6, " iterations!")') decomp%cfg_hlgap_nit_conv
      else
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of alpha occupied orbital energies did not converge")')
      endif
   endif

!2. Find eigenvalues of beta part of FUP
   if (decomp%info_stability) then
      write (decomp%lupri,*)
      write (decomp%lupri,'("Beta occupied orbitals:")')
      write (decomp%lupri,'("=======================")')
   endif
   !Set halfmat to beta part of FUP:
   halfmat(1:matdim,1:matdim) = fullmat1(matdim+1:fulldim,matdim+1:fulldim)
   !Set proj to beta part of DUone:
   proj(1:matdim,1:matdim) = fullmat2(matdim+1:fulldim,matdim+1:fulldim)
       !write(lupri,*) 'FUP halfmat, beta part **'
       !call LS_OUTPUT(halfmat, 1, matdim, 1, matdim, matdim, matdim, 1, lupri)
   !norm = frob_norm(halfmat,matdim,matdim)
   !if (norm < 1.0E-12_realk) then
   if (DD%noccb == 0) then
       NoBeta = .true.
       write (decomp%lupri,'("There appears to be no beta occupied orbitals, continuing...")')
   else
      call DD_Fock_eigenvalue_unres('b',DD,decomp,halfmat,proj,'h',0.0E0_realk)
      if (decomp%cfg_hlgap_converged) then
         betaocc_conv = .true.
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of beta occupied orbital energies converged in", &
           & i6, " iterations!")') decomp%cfg_hlgap_nit_conv
      else
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of beta occupied orbital energies did not converge")')
      endif
   endif

!3. Find eigenvalues of alfa part of FUQ
   if (decomp%info_stability) then
      write (decomp%lupri,*)
      write (decomp%lupri,'("Alfa virtual orbitals:")')
      write (decomp%lupri,'("=====================")')
   endif
   fullmat1 = 0.0E0_realk ; fullmat2 = 0.0E0_realk
   call mat_to_full(decomp%FUQ,1.0E0_realk,fullmat1)
   call mat_to_full(decomp%QU,1.0E0_realk,fullmat2)
   !Set halfmat to alfa part of FUQ:
   halfmat(1:matdim,1:matdim) = fullmat1(1:matdim,1:matdim)
   !Set proj to alfa part of QU:
   proj(1:matdim,1:matdim) = fullmat2(1:matdim,1:matdim)
   if (DD%nvirta == 0) then
       NoAlfa = .true.
       write (decomp%lupri,'("There appears to be no alfa virtual orbitals, continuing...")')
   else
      call DD_Fock_eigenvalue_unres('a',DD,decomp,halfmat,proj,'l',0.0E0_realk)
      if (decomp%cfg_hlgap_converged) then
         alfavirt_conv = .true.
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of alpha virtual orbital energies converged in", &
           & i6, " iterations!")') decomp%cfg_hlgap_nit_conv
      else
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of alpha virtual orbital energies did not converge")')
      endif
   endif

!4. Find eigenvalues of beta part of FUQ
   if (decomp%info_stability) then
      write (decomp%lupri,*)
      write (decomp%lupri,'("Beta virtual orbitals:")')
      write (decomp%lupri,'("=====================")')
   endif
   !Set halfmat to alfa part of FUQ:
   halfmat(1:matdim,1:matdim) = fullmat1(matdim+1:fulldim,matdim+1:fulldim)
       !write(lupri,*) 'FUQ halfmat, beta part **'
       !call LS_OUTPUT(halfmat, 1, matdim, 1, matdim, matdim, matdim, 1, lupri)
   !Set proj to alfa part of QU:
   proj(1:matdim,1:matdim) = fullmat2(matdim+1:fulldim,matdim+1:fulldim)
   !norm = frob_norm(halfmat,matdim,matdim)
   !if (norm < 1.0E-12_realk) then
   if (DD%nvirtb == 0) then
       NoBeta = .true.
       write (decomp%lupri,'("There appears to be no beta virtual orbitals, continuing...")')
   else
      call DD_Fock_eigenvalue_unres('b',DD,decomp,halfmat,proj,'l',0.0E0_realk)
      if (decomp%cfg_hlgap_converged) then
         betavirt_conv = .true.
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of beta virtual orbital energies converged in", &
           & i6, " iterations!")') decomp%cfg_hlgap_nit_conv
      else
         write (decomp%lupri,*)
         write (decomp%lupri,'("Calculation of beta virtual orbital energies did not converge")')
      endif
   endif

   !do i = 1, DD_nocc !Very quick fix for debugging small systems, should be better
   !   if (abs(DD_FUPalfa_eival(i)) < 1.0E-12_realk) then
   !      DD_FUPalfa_eival(i) = -1.0E3_realk
   !   else if (abs(DD_FUPbeta_eival(i)) < 1.0E-12_realk) then
   !      DD_FUPbeta_eival(i) = -1.0E3_realk
   !   endif
   !enddo

   if (NoAlfa) then
      write (decomp%lupri,*)
      write (decomp%lupri,*) "There is no alpha gap!"
      write (decomp%lupri,*)
   else if (DD%nocca > 0) then
      if (alfaocc_conv .and. alfavirt_conv) then
         HOMOLUMOgap_a = DD%FUQalfa_eival(1)-DD%FUPalfa_eival(1)
         write (decomp%lupri,*)
         write (decomp%lupri,'("    E(LUMO), alpha:                      ", F12.6, " au")') DD%FUQalfa_eival(1)
         write (decomp%lupri,'("   -E(HOMO), alpha:                      ", F12.6, " au")') DD%FUPalfa_eival(1)
         write (decomp%lupri,'("   -------------------------------------")')
         write (decomp%lupri,'("    ""Alpha"" HOMO-LUMO Gap (iteratively): ", F12.6, " au")') HOMOLUMOgap_a
         write (decomp%lupri,*)
      else
         write (decomp%lupri,*)
         write (decomp%lupri,*) "Alpha HOMO-LUMO gap cannot be printed"
      endif
   else
      write (decomp%lupri,*)
      write (decomp%lupri,*) "No alpha electrons!"
      write (decomp%lupri,*)
   endif

   if (NoBeta) then
      write (decomp%lupri,*)
      write (decomp%lupri,*) "There is no beta gap!"
      write (decomp%lupri,*)
   else if (DD%noccb > 0) then
      if (betaocc_conv .and. betavirt_conv) then
         HOMOLUMOgap_b = DD%FUQbeta_eival(1)-DD%FUPbeta_eival(1)
         write (decomp%lupri,*)
         write (decomp%lupri,'("    E(LUMO), beta:                       ", F12.6, " au")') DD%FUQbeta_eival(1)
         write (decomp%lupri,'("   -E(HOMO), beta:                       ", F12.6, " au")') DD%FUPbeta_eival(1)
         write (decomp%lupri,'("   -------------------------------------")')
         write (decomp%lupri,'("    ""Beta"" HOMO-LUMO Gap (iteratively):  ", F12.6, " au")') HOMOLUMOgap_b
         write (decomp%lupri,*)
      else
         write (decomp%lupri,*)
         write (decomp%lupri,*) "Beta HOMO-LUMO gap cannot be printed"
      endif
   else
      write (decomp%lupri,*)
      write (decomp%lupri,*) "No beta electrons!"
      write (decomp%lupri,*)
   endif

   if (DD%noccb > 0 .and. DD%nocca > 0) then
      if (alfaocc_conv .and. alfavirt_conv .and. betaocc_conv .and. betavirt_conv) then
         if (.not. NoAlfa .and. .not. NoBeta) then
            TOTALgap = min(DD%FUQalfa_eival(1),DD%FUQbeta_eival(1)) - max(DD%FUPalfa_eival(1),DD%FUPbeta_eival(1))
            write (decomp%lupri,*)
            write (decomp%lupri,'("    E(LUMO):                             ", F12.6, " au")') &
                     & min(DD%FUQalfa_eival(1),DD%FUQbeta_eival(1))
            write (decomp%lupri,'("   -E(HOMO):                             ", F12.6, " au")') &
                     & max(DD%FUPalfa_eival(1),DD%FUPbeta_eival(1))
            write (decomp%lupri,'("   -------------------------------------")')                                       
            write (decomp%lupri,'("    ""Overall"" HOMO-LUMO Gap (iteratively): ", F10.6, " au")') TOTALgap           
            write (decomp%lupri,*)
         endif
      else
         write (decomp%lupri,*)
         write (decomp%lupri,*) "Overall HOMO-LUMO gap cannot be printed"
         write (decomp%lupri,*)
      endif
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

      allocate(iniguess_full(matdim,matdim))
      allocate(alfa_orb_energy_diff(DD%nocca*DD%nvirta))
      allocate(beta_orb_energy_diff(DD%noccb*DD%nvirtb))
      allocate(alfaFUQFUPindex(DD%nocca*DD%nvirta,2))
      allocate(betaFUQFUPindex(DD%noccb*DD%nvirtb,2))

      if (decomp%info_stability) then
         write(decomp%lupri,*)
         write (decomp%lupri,'("Checking stability of solution: Calculation of lowest Hessian eigenvalue")')
         write (decomp%lupri,'("========================================================================")')
      endif

      if (.not. NoAlfa .and. decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) '       Alpha occ. eigenvalues' !; call flshfo(decomp%lupri)
         write(decomp%lupri,*) '       ======================'
         CALL LS_OUTPUT(DD%FUPalfa_eival,1,DD%nocca,1,1,DD%nocca,1,1,decomp%LUPRI)
         !CALL LS_OUTPUT(DD_FUPalfa_eivecs,1,matdim,1,DD_nocc,matdim,DD_nocc,1,LUPRI)
      endif

      if (.not. NoBeta .and. decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) '       Beta occ. eigenvalues' !; call flshfo(decomp%lupri)
         write(decomp%lupri,*) '       ====================='
         CALL LS_OUTPUT(DD%FUPbeta_eival,1,DD%noccb,1,1,DD%noccb,1,1,decomp%LUPRI)
         !CALL LS_OUTPUT(DD_FUPbeta_eivecs,1,matdim,1,DD_nocc,matdim,DD_nocc,1,LUPRI)
      endif

      if (.not. NoAlfa .and. decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) '       Alpha virt. eigenvalues' !; call flshfo(lupri)
         write(decomp%lupri,*) '       ======================='
         CALL LS_OUTPUT(DD%FUQalfa_eival,1,DD%nvirta,1,1,DD%nvirta,1,1,decomp%LUPRI)
         !CALL LS_OUTPUT(DD_FUQalfa_eivecs,1,matdim,1,DD_nvirt,matdim,DD_nvirt,1,LUPRI)
      endif

      if (.not. NoBeta .and. decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) '       Beta virt. eigenvalues' !; call flshfo(lupri)
         write(decomp%lupri,*) '       ======================'
         CALL LS_OUTPUT(DD%FUQbeta_eival,1,DD%nvirtb,1,1,DD%nvirtb,1,1,decomp%LUPRI)
         !CALL LS_OUTPUT(DD_FUQbeta_eivecs,1,matdim,1,DD_nvirt,matdim,DD_nvirt,1,LUPRI)
      endif

      k = 0
      !Calculate alpha orbital energy differences obtained from FUP/FUP eigenvalues
      do i = 1, DD%nvirta !FUQeival
         do j = 1, DD%nocca !FUPeival
            k = k + 1
            alfa_orb_energy_diff(k) = DD%FUQalfa_eival(i) - DD%FUPalfa_eival(j)
            alfaFUQFUPindex(k,1) = i
            alfaFUQFUPindex(k,2) = j
         enddo
      enddo

      k = 0
      !Calculate beta orbital energy differences obtained from FUP/FUP eigenvalues
      do i = 1, DD%nvirtb !FUQeival
         do j = 1, DD%noccb !FUPeival
            k = k + 1
            beta_orb_energy_diff(k) = DD%FUQbeta_eival(i) - DD%FUPbeta_eival(j)
            betaFUQFUPindex(k,1) = i
            betaFUQFUPindex(k,2) = j
         enddo
      enddo

      if (.not. NoAlfa .and. decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) 'Alpha orbital energy differences  (occ index, virt index):' !; call flshfo(lupri)
         write(decomp%lupri,*) '=========================================================='
         do i = 1, DD%nvirta*DD%nocca
            write(decomp%lupri, '(i5, F18.7, i19, i11)') i, alfa_orb_energy_diff(i), alfaFUQFUPindex(i,1), alfaFUQFUPindex(i,2)
         enddo
      endif

      if (.not. NoBeta .and. decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) 'Beta orbital energy differences  (occ index, virt index):' !; call flshfo(lupri)
         write(decomp%lupri,*) '========================================================='
         do i = 1, DD%nvirtb*DD%noccb
            write(decomp%lupri, '(i5, F18.7, i19, i11)') i, beta_orb_energy_diff(i), betaFUQFUPindex(i,1), betaFUQFUPindex(i,2)
         enddo
      endif

      if (decomp%info_stability) then
         write(decomp%lupri,*)
         write(decomp%lupri,*) '    Sorted orbital energy differences:'
         write(decomp%lupri,*) '=========================================='
      endif

      do i = 1, nguesses
         alfa = .false.
         !First determine if smallest orb energy diff is in alfa or beta part:
         if (NoAlfa) then
            alfa = .false.
         else if (NoBeta) then
            alfa = .true.
         else
            alfamin = MINVAL(alfa_orb_energy_diff) ; betamin = MINVAL(beta_orb_energy_diff)
            if (alfamin < betamin .and. DD%nocca > 0) then
               alfa = .true.
            endif
         endif
         if (alfa) then
            posmin = MINLOC(alfa_orb_energy_diff)
            posmax = MAXLOC(alfa_orb_energy_diff)
            if (decomp%info_stability) write (decomp%lupri, "('Root no.', i4, F12.7, '   Position:', i4, ' (ALPHA)')") &
                                  & i, alfa_orb_energy_diff(posmin(1)),  posmin(1)
            if (do_rsp_iniguess) then
               KK1 = alfaFUQFUPindex(posmin(1),1)
               KK2 = alfaFUQFUPindex(posmin(1),2)
               do nu = 1, matdim
                  TMP1 = DD%FUPalfa_eivecs(nu,KK2)
                  do mu = 1, matdim
                     iniguess_full(mu,nu) = DD%FUQalfa_eivecs(mu,KK1)*TMP1
                  enddo
               enddo
            else
               KK1 = alfaFUQFUPindex(posmin(1),1)
               KK2 = alfaFUQFUPindex(posmin(1),2)
               do nu = 1, matdim
                  TMP1 = DD%FUPalfa_eivecs(nu,KK2)
                  TMP2 = DD%FUQalfa_eivecs(nu,KK2)
                  do mu = 1, matdim
                     iniguess_full(mu,nu) = DD%FUQalfa_eivecs(mu,KK1)*TMP1- &
                          &  DD%FUPalfa_eivecs(mu,KK1)*TMP2
                  enddo
               enddo
            endif
            !write(lupri,*) 'full Hes starting guess no', i
            !call LS_OUTPUT(iniguess_full, 1, matdim, 1, matdim, matdim, matdim, 1, lupri)

            alfa_orb_energy_diff(posmin(1)) = alfa_orb_energy_diff(posmax(1))
            call mat_unres_dense_part_from_full(iniguess_full,'a',iniguess(i))
            call mat_unres_dense_part_from_full(iniguess_full,'b',iniguess(i)) !Thomas, Stinne 19/2-10
            !iniguess%elms(1:matdim*matdim) = iniguess_full
            !call mat_set_from_full(iniguess_full,1.0E0_realk, iniguess(i))
         else
            posmin = MINLOC(beta_orb_energy_diff)
            posmax = MAXLOC(beta_orb_energy_diff)
            if (decomp%info_stability) write (decomp%lupri, "('Root no.', i4, F12.7, '   Position:', i4, ' (BETA)')") &
                                  & i, beta_orb_energy_diff(posmin(1)),  posmin(1)
            if (do_rsp_iniguess) then
               KK1 = betaFUQFUPindex(posmin(1),1)
               KK2 = betaFUQFUPindex(posmin(1),2)
               do nu = 1, matdim
                  TMP1 = DD%FUPbeta_eivecs(nu,KK2)
                  do mu = 1, matdim
                     iniguess_full(mu,nu) = DD%FUQbeta_eivecs(mu,KK1)*TMP1
                  enddo
               enddo
            else
               KK1 = betaFUQFUPindex(posmin(1),1)
               KK2 = betaFUQFUPindex(posmin(1),2)
               do nu = 1, matdim
                  TMP1 = DD%FUPbeta_eivecs(nu,KK2)
                  TMP2 = DD%FUQbeta_eivecs(nu,KK2)
                  do mu = 1, matdim
                     iniguess_full(mu,nu) = DD%FUQbeta_eivecs(mu,KK1)*TMP1 &
                          &  - DD%FUPbeta_eivecs(mu,KK1)*TMP2
                  enddo
               enddo
            endif
            !write(lupri,*) 'full Hes starting guess no', i
            !call LS_OUTPUT(iniguess_full, 1, matdim, 1, matdim, matdim, matdim, 1, lupri)
            beta_orb_energy_diff(posmin(1)) = beta_orb_energy_diff(posmax(1))
            call mat_unres_dense_part_from_full(iniguess_full,'a',iniguess(i)) !Thomas, Stinne 17/2-10
            call mat_unres_dense_part_from_full(iniguess_full,'b',iniguess(i))
            !iniguess%elmsb(1:matdim*matdim) = iniguess_full
            !call mat_set_from_full(iniguess_full,1.0E0_realk, iniguess(i))
         endif

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
      !Construct starting guess from eigenvalues from FUP/FUQ eigenvectors

      !do mu = 1, matdim
      !   do nu = 1, matdim
      !      iniguess_full(mu,nu) = DD_FUQ_eivecs(mu,1)*DD_FUP_eivecs(nu,1) - &  !We use lowest solution
      !                          &  DD_FUP_eivecs(mu,1)*DD_FUQ_eivecs(nu,1)
      !   enddo
      !enddo

      !call mat_set_from_full(iniguess_full,1.0E0_realk, iniguess, 'iniguess')
      !!Normalize X*X=1 
      !fac = abs(mat_dotproduct(iniguess,iniguess))
      !call mat_scal(1.0E0_realk/sqrt(fac),iniguess)

      if (decomp%cfg_check_converged_solution) then
      !if (.not. debug_no_hessian) then !A bit ugly, I know... for debugging arh hessian which is
                                       !changed when we do a rejection, as opposed to the 'real' Hessian 
         !write(lupri,*) 'Calling hes eival' ; call flshfo(lupri)
         !cfg_do_2nd_order = .true. 
!HACK
!debug_dd_lintra = .true.
!END HACK
         call DD_Hessian_eigenvalue(decomp, DD, iniguess, nguesses, hessian_eigenval,dummyfifo)
         !cfg_do_2nd_order = .false.

         write(decomp%lupri,*)
         write (decomp%lupri,'("*======================================================================*")')
         do i = 1, decomp%cfg_hessian_nvec
            write (decomp%lupri,'("                Hessian eigenvalue no.",i3,": ",F12.6, "      ")') i, hessian_eigenval(i) 
         enddo
         write (decomp%lupri,'("*======================================================================*")')
         !call flshfo(decomp%lupri)
         debug%heseival = hessian_eigenval(1)
      !endif
      endif

      !if (cfg_density_method == cfg_f2d_arh .and. debug_arh_hessian) then
      !   finding_arh_eigenvalue = .true.
      !   if (.true.) then !First find HL gap w arh linear transformation
      !      cfg_arhterms = .false.
      !      call DD_Hessian_eigenvalue(iniguess, hessian_eigenval,fifoqueue)
      !      cfg_arhterms = .true.
      !      write (lupri,'("*======================================================================*")')
      !      write (lupri,'("             HOMO-LUMO gap using ARH linear transf: ",F12.6, "      ")') hessian_eigenval
      !      write (lupri,'("*======================================================================*")')

      !   endif
      !   call DD_Hessian_eigenvalue(iniguess, hessian_eigenval,fifoqueue)
      !   finding_arh_eigenvalue = .false.
      !   write(lupri,*)
      !   write (lupri,'("*======================================================================*")')
      !   write (lupri,'("                ARH lowest Hessian eigenvalue: ",F12.6, "      ")') hessian_eigenval 
      !   write (lupri,'("*======================================================================*")')
      !   debug_arh_arheival = hessian_eigenval(1)
      !endif
      !call flshfo(lupri)

      do i = 1, nguesses
         call mat_free(iniguess(i))
      enddo
      deallocate(iniguess_full)
      deallocate(alfa_orb_energy_diff, beta_orb_energy_diff)
      deallocate(alfaFUQFUPindex, betaFUQFUPindex)
      !call gettim(cpu2, wall2)
      !call util_print_time('HESEIG',cpu1,cpu2,wall1,wall2)
      CALL LSTIMER('HESEIG',T1,T2,decomp%lupri)
   endif
   end subroutine DD_homolumo_and_heseigen_unres

   !> \brief Find lowest or highest eigenvalue of n x n matrix M by inverse iteration and reduced space
   !> \author S. Host
   !> \date 2005
   subroutine DD_Fock_eigenvalue_unres(spintype, DD, decomp, M, P, desired_eigenval, guess)
   implicit none
         !> Find orbital energies of alpha part (spintype='a') or beta part (spintype='b')
         character, intent(in)       :: spintype
         !> Contains settings for HOMO-LUMO gap and Hessian eigenvalues
         type(DDitem),intent(inout)  :: DD
         !> Contains settings for decomposition and OAO decomposed overlap
         type(decompItem),intent(inout) :: decomp
         !> Find highest or lowest eigenvalue(s) of this matrix
         real(realk), intent(inout)  :: M(:,:)
         !> After preconditioning, project with this matrix
         real(realk), intent(inout)  :: P(:,:)                
         !> h = find highest eigenvalue(s), l = find lowest eigenvalue(s)
         character, intent(in)       :: desired_eigenval
         !> Should be close to the desired eigenvalue for efficient solution.
         real(realk), intent(in)     :: guess 
         integer                  :: i, j, nroots, nocc, nvirt, ndim
         real(realk), allocatable, dimension(:,:) :: iniguess_full, M_full, P_full, Mdamp, x
         real(realk), allocatable, dimension(:) :: eigenval

      ndim = decomp%U%nrow

      if (spintype == 'a') then
         nocc = DD%nocca ; nvirt = DD%nvirta
      else if (spintype == 'b') then
         nocc = DD%noccb ; nvirt = DD%nvirtb
      else
         stop 'Moron!' 
      endif

      if (desired_eigenval == 'h') then
         nroots = nocc
      else if (desired_eigenval == 'l') then
         nroots = nvirt
      else
         stop 'Moron!' 
      endif

      allocate(Mdamp(ndim,ndim))
      allocate(iniguess_full(ndim,nroots),x(ndim,nroots))
      allocate(eigenval(nroots))

      !write(lupri,*) 'M:'
      !call mat_print(M, 1, M%nrow, 1, M%ncol, lupri)
      !write(lupri,*) 'P:'
      !call mat_print(P, 1, M%nrow, 1, M%ncol, lupri)

      !Damp M:
      Mdamp = M
      do i = 1, ndim
         Mdamp(i,i) = Mdamp(i,i) - guess
      enddo 

      !Find nroots starting guesses:
      call DD_initial_guess(M,ndim,nroots,iniguess_full)

      !Project starting guesses: This is done in orthonormalize inside solver
      x = iniguess_full

      if (desired_eigenval == 'h') then !Switch sign on matrix
         M = -M
         call DD_Fock_solver(decomp, DD, ndim, nroots, Mdamp, M, x, P, eigenval)
         !Save eigenvectors so they can be used to generate starting guesses for Hessian eigenvalue
         !and/or excitation energies
         !eigenval = -eigenval
         if (spintype == 'a') then
            DD%FUPalfa_eival = -eigenval
            do i = 1, nocc
               call DD_Fock_project(ndim,P,x(1:ndim,i),DD%FUPalfa_eivecs(1:ndim,i))
            enddo
         else if (spintype == 'b') then
            DD%FUPbeta_eival = -eigenval
            do i = 1, nocc
               call DD_Fock_project(ndim,P,x(1:ndim,i),DD%FUPbeta_eivecs(1:ndim,i))
            enddo
         endif
      else if (desired_eigenval == 'l') then
         call DD_Fock_solver(decomp, DD, ndim, nroots, Mdamp, M, x, P, eigenval)
         if (spintype == 'a') then
            DD%FUQalfa_eival = eigenval
            do i = 1, nvirt
               call DD_Fock_project(ndim,P,x(1:ndim,i),DD%FUQalfa_eivecs(1:ndim,i))
            enddo
         else if (spintype == 'b') then
            DD%FUQbeta_eival = eigenval
            do i = 1, nvirt
               call DD_Fock_project(ndim,P,x(1:ndim,i),DD%FUQbeta_eivecs(1:ndim,i))
            enddo
         endif
      else
         STOP 'Unknown type of eigenvalue (DD_eigenvalue)'
      endif

      deallocate(x)
      deallocate(Mdamp)
      deallocate(iniguess_full)
      deallocate(eigenval)
   end subroutine DD_Fock_eigenvalue_unres 

end module direct_dens_util_unres
