! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module rsp_equations

! ajt FIXME Add detailed description of contents here
!
! ajt jan10 Added place to store (cache) solutions of response
!           equations, solved_eqs(100), and avoid re-solving the same
!           equations later. Solved_eqs can also hold excitation densities,
!           used when solving for higher-order transition densities,
!           perturbed transition densities, as well as projected perturbed
!           densities.
!
! ajt feb10 Added solving for the non-resonant component of resonant
!           equations. Separated some utility routines from rsp_solve_eq.
!
!> Delivers perturbed density matrices (solutions of response equations).
!> Property integrals are obtained through module rsp_contribs.
!> Right-hand-sides are contracted here, and solutions are obtained
!> by invoking one of several response solvers, configured at compile
!> time by #define statements.

!  version history (response code, not just this file):
!  2010-04-21 identical (shared) files: dalton/branches/linsca-ng revision  7878
!                                       dirac/trunk               revision 10659


module lsdalton_rsp_equations

   use precision
   use lsdalton_matrix_defop  !matrix type and operators
   use lsdalton_rsp_contribs  !integrals and integral contractions
   use RSPsolver, only: rsp_molcfg
   implicit none

   public pert_fock    !contract perturbed Fock-matrices
   public pert_scf_eq  !contract perturbed scf equation residuals
   public pert_dens    !perturbed (or response-) densities and Fock-matrices
   public rsp_equations_debug  !turn debug printing on or off
   public rsp_eq_sol  !saved solutions
   public rsp_eq_sol_grow  !reallocate rsp_eq_sol larger
   public rsp_eq_sol_empty !empty and deallocate rsp_eq_sol
   public rsp_eq_truncate_order !above this order return zero instead of solving equations

   ! ajt Freqency derivative hack, not intended to be public
   public rsp_solve_eq
   public rsp_saved_sol

   !> turn on or off debugging in this file
   logical :: rsp_equations_debug = .false.

   integer :: rsp_eq_truncate_order = huge(1)

   !> Type for saving/caching the solution of a response equation,
   !> so to avoid re-solving the same equation later in the program.
   !> Note that D is not a vector, thus all fld%ncomp == 1
   type rsp_saved_sol
      !> molecule and settings for this equation
      type(rsp_molcfg), pointer :: mol
      !> Field lables, component indices and frequencies.
      !> Should be ordered according to:
      !> 1) decreasing label (index within field_list in rsp_contribs),
      !> 2) increasing component indices %comp,
      !> 2) decreasing abs(dreal(frequency)) (- before +),
      !> 4) decreasing abs(imag(frequency)) (- before +),
      !> function pert_order(fld(:)) in rsp_contribs gives proper ordering.
      type(pert_field), pointer :: fld(:)
      !> the density matrix solving the response equation
      !> ajt FIXME keep the 'actual' matrix here, but return aliases
      !>           of it from pert_dens, until the last use,
      !>           when the actual matrix is returned, and D is cleared.
      type(matrix) :: D
      !> reference count. Increment while rsp_eq_counting, compute on first
      !> call after counting, free when it hits zero.
      !> ajt FIXME counting is not implemented yet.
      integer :: refc
      !> ajt FIXME other matrices too: F, SD, FD, DFDp, FpDS, DpSD.
      !>           Also, to save memory, files can be used,
      !>           with filenames stored here.
   end type

   !> To keep collection of saved response equation solutions.
   !> Grows in units of 8 when needed
   type(rsp_saved_sol), allocatable :: rsp_eq_sol(:)

   ! ajt For some reason, if I put this 'private' above type cached_sol,
   !     Doxygen does not document the fields inside the type, dispite
   private !it being declared public at the top. Might be a bug in Doxygen.

contains


   !> Contract perturbed integrals with (an expansion of) perturbed density
   !> matrices to get the corresponding perturbed Fock matrix.
   !> Optionally, also return the corresponding perturbed overlap integrals
   subroutine pert_fock(mol, p, dimp, D, F, S, comp, freq)
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(rsp_molcfg), intent(inout) :: mol
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> dimension of each p, and thus also of F
      integer,           intent(in) :: dimp(:)
      !> density matrix expansion for p(:) unperturbed and perturbed
      !> density matrices. size(D) = product(1+dimp)
      type(matrix),      intent(in) :: D(:)
      !> output, perturbed Fock matrices. Deferred shape (*) to permit
      !> re-ranking. size(F)=product(dimp)
      type(matrix),   intent(inout) :: F(*)
      !--------------------------------------------------------------------
      !> optionally also output perturbed overlap matrices.
      !> Deferred shape (*) to permit re-ranking
      type(matrix),   optional, intent(inout) :: S(*)
      !> starting component within each p, default 1 1 .. 1
      integer,        optional, intent(in)    :: comp(:)
      !> frequency of each p default all zero
      complex(realk), optional, intent(in)    :: freq(:)
      !-----------------------------------------------
      integer        :: ccomp(size(p)), nf, nd
      complex(realk) :: ffreq(size(p))
      character(4)   :: UNP(0)
      logical        :: bas(size(p))
      if (size(dimp) /= size(p)) call lsquit('pert_fock: Differing numbers ' &
               // 'of perturbations and dimensions, size(dimp) /= size(p)',mol%lupri)
      nf = product(dimp)
      nd = product(1+dimp)
      if (size(D) /= nd) call lsquit('pert_fock: Wrong number of perturbed ' &
               // 'densities, size(D) /= product(dimp+1)',mol%lupri)
      ccomp = 1
      if (present(comp)) then
         ccomp = comp !ajt fixme Check number and bounds
      end if
      ffreq = 0
      if (present(freq)) then
         ffreq = freq !ajt fixme Check number
      end if
      ! fill F with perturbed one-electron integrals, and
      ! optionally S with perturbed overlap integrals
      call rsp_oneint(mol, D(1), p, dimp, F, S, comp=ccomp, freq=ffreq)
      ! add p-perturbed Coulomb-exchange+Kohn-Sham contracted against
      call rsp_twoint(mol, p, D(1:1), dimp, F, comp=ccomp) !unperturbed D
      ! add unperturbed Coulomb-exchange+Kohn-Sham contribution
      call rsp_twoint(mol, UNP, D, dimp, F)
      ! additional contributions for each basis-dependent subset of p
      bas = pert_basdep(p)
      if (size(p)==2 .and. bas(size(p))) then
         call rsp_twoint(mol, p(2:2), D(1:1+dimp(1)), &
              dimp, F, perm=(/2,1/), comp=ccomp(2:2))
      endif
      if (size(p)==2 .and. bas(1))then
         call rsp_twoint(mol, p(1:1), (/D(1),D(1+dimp(1)+1:1+dimp(1)+dimp(2))/), &
              dimp, F, comp=ccomp(1:1))
      endif
      if (size(p) > 2 .and. any(bas)) &
           call lsquit('pert_fock: general case not implemented',mol%lupri)
   end subroutine


   !> Computes the p'th derivative of idempotency relation DSD-D
   !> and SCF equation FDS-SDF, ie. the p-perturbed SCF equations.
   subroutine pert_scf_eq(mol, S0, p, dimp, Dp, Fp, FDSp, DSDp, comp, freq)
     implicit none
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(rsp_molcfg), intent(inout) :: mol
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> number of components within each p (dimension)
      integer,           intent(in) :: dimp(:)
      !> perturbation expantion of density matrix
      type(matrix),      intent(in) :: Dp(:)
      !> perturbation expantion of Fock matrix
      type(matrix),      intent(in) :: Fp(:)
      !-----------------------------------
      !> resulting perturbed SCF equation (or residual)
      !> Deferred shape (*) to permit any rank. size(FDSp) = product(dimp)
      type(matrix),   intent(inout) :: FDSp(*)
      !> resulting perturbed idempotency condition
      !> Deferred shape (*) to permit any rank, size(DSDp) = product(dimp)
      type(matrix),   intent(inout) :: DSDp(*)
      !-------------------------------------
      !> starting indices for each p, default 1
      integer,        optional, intent(in) :: comp(:)
      !> frequency of each p, default 0
      complex(realk), optional, intent(in) :: freq(:)
      !------------------------------------------
      type(matrix), dimension(size(Dp)) :: Sp, DSp, SDp, FDp, DFp
      integer        :: np, nd, i, j, i0, j0, ij, c(size(p)) !copy of comp or 1
      complex(realk) :: w(size(p)) !copy of freq or zero
      character(4)   :: UNP(0)
      logical        :: bas(size(p)), presnt
      i=0
      j=0
      ij=0
      np = product(dimp)
      nd = product(dimp+1)
      bas = pert_basdep(p)
      if (size(Dp) /= nd) call lsquit('pert_scf_eq: Wrong number of perturbed ' &
                // 'densities, size(Dp) /= product(dimp+1)',mol%lupri)
      if (size(Fp) /= nd) call lsquit('pert_scf_eq: Wrong number of perturbed ' &
                // 'densities, size(Dp) /= product(dimp+1)',mol%lupri)
      c = 1
      if (present(comp)) then
         c = comp !ajt fixme Check number and bounds
      end if
      w = 0
      if (present(freq)) then
         w = freq !ajt fixme Check number
      end if
      ! start by zeroing FDSp and DSDp
      do i = 1, np
         FDSp(i) = 0E0_realk * S0
         DSDp(i) = 0E0_realk * S0
      end do
      ! zero perturbed overlap Sp, and FDp/DFp (which multiply perturbed overlap)
      Sp(1)  = 0E0_realk * S0
      Sp(2:) = (/(Sp(1), i=2,size(Sp)) /)
      FDp(:) = (/(Sp(1), i=1,size(FDp))/)
      DFp(:) = (/(Sp(1), i=1,size(DFp))/)
      ! Contract -RHS's. The contraction formulas are derivatives of the
      ! rearranged first-order formula:
      ! d/da (F - i/2 d/dt S) D S - S D (F + i/2 tb\b S) =
      !     = (Fa - i/2 Sat) D S - S D (Fa + i/2 Sat)
      !     + (F D - i/2 St D - i S Dt) Sa - Sa (D F + i/2 D St + i Dt S)
      !     + (F Da - i/2 St Da - i/2 S Dat) S - S (Da F + i/2 Da St + i/2 Dat S)
      if (size(p)==0) then
         ! trivial case: Unperturbed equations
         FDSp(1) = Fp(1)*Dp(1)*S0 - S0*Dp(1)*Fp(1)
         DSDp(1) = Dp(1)*S0*Dp(1) - Dp(1)
      else if (size(p)==1) then
         presnt = .not.all((/(iszero(Dp(i)), i=2,1+dimp(1))/))
         if (presnt) then
            DSp(1) = Dp(1)*S0
            if (bas(1)) FDp(1) = Fp(1)*Dp(1)
            if (bas(1)) call rsp_oneint(mol, S0, p, dimp, S=Sp(2:1+dimp(1)), comp=c)
            do i0 = 1, dimp(1)
               i = 1+i0
               FDSp(i0) =            (Fp(i) - w(1)/2 * Sp(i)) * DSp(1)
               FDSp(i0) = FDSp(i0) - trans(DSp(1)) * (Fp(i) + w(1)/2 * Sp(i))
               FDSp(i0) = FDSp(i0) + FDp(1) * Sp(i) - Sp(i) * trans(FDp(1))
               FDSp(i0) = FDSp(i0) + (Fp(1) - w(1)/2 * S0) * Dp(i) * S0
               FDSp(i0) = FDSp(i0) - S0 * Dp(i) * (Fp(1) + w(1)/2 * S0)
               DSDp(i0) =            Dp(1) * Sp(i) * Dp(1) - Dp(i)
               DSDp(i0) = DSDp(i0) + Dp(i)*trans(DSp(1))
               DSDp(i0) = DSDp(i0) + DSp(1)*Dp(i)
               Sp(i)=0
            end do
            DSp(1)=0; FDp(1)=0
         end if
      else if (size(p)==2) then
         presnt = .not.all((/(iszero(Dp(ij)), ij = 2+dimp(1)+dimp(2), &
                                              (dimp(1)+1)*(dimp(2)+1))/))
         if (presnt) then
            DSp(1) = Dp(1) * S0
            if (all(bas)) FDp(1) = Fp(1) * Dp(1)
            if (all(bas))call rsp_oneint(mol, S0, p, dimp, S=Sp(2+dimp(1)+dimp(2) : &
                                         (dimp(1)+1)*(dimp(2)+1)), comp=c)
            do i0 = 1, dimp(1)*dimp(2)
               i = 1+dimp(1)+dimp(2)+i0
               FDSp(i0) =            (Fp(i) - (w(1)+w(2))/2 * Sp(i)) * DSp(1)
               FDSp(i0) = FDSp(i0) - trans(DSp(1)) * (Fp(i) + (w(1)+w(2))/2 * Sp(i))
               FDSp(i0) = FDSp(i0) + FDp(1) * Sp(i) - Sp(i) * trans(FDp(1))
               FDSp(i0) = FDSp(i0) + (Fp(1) - (w(1)+w(2))/2 * S0) * Dp(i) * S0
               FDSp(i0) = FDSp(i0) - S0 * Dp(i) * (Fp(1) + (w(1)+w(2))/2 * S0)
               DSDp(i0) =            Dp(1) * Sp(i) * Dp(1) - Dp(i)
               DSDp(i0) = DSDp(i0) + Dp(i)*trans(DSp(1))
               DSDp(i0) = DSDp(i0) + DSp(1)*Dp(i)
               Sp(i)=0
            end do
            DSp(1)=0; FDp(1)=0
         end if
         presnt = .not.all((/(iszero(Dp(i)), i=2,1+dimp(1))/)) .and. &
                  .not.all((/(iszero(Dp(j)), j=2+dimp(1),1+dimp(1)+dimp(2))/))
         if (presnt) then
            if (bas(1)) call rsp_oneint(mol, S0, p(1:1), dimp(1:1), &
                                        S=Sp(2:1+dimp(1)), comp=c(1:1))
            if (bas(2)) call rsp_oneint(mol, S0, p(2:2), dimp(2:2), S=Sp(2+dimp(1): &
                                        1+dimp(1)+dimp(2)), comp=c(2:2))
            do j0 = 0, dimp(2)-1
               j = 2+dimp(1)+j0
               DSp(j) = Dp(j)*S0 + Dp(1)*Sp(j)
               SDp(j) = S0*Dp(j) + Sp(j)*Dp(1)
               if (bas(1)) FDp(j) = (Fp(j) - w(2)/2 * Sp(j)) * Dp(1) &
                                  + (Fp(1) - w(2) * S0) * Dp(j)
               if (bas(1)) DFp(j) = Dp(1) * (Fp(j) + w(2)/2 * Sp(j)) &
                                  + Dp(j) * (Fp(1) + w(2) * S0)
               do i0 = 0, dimp(1)-1
                  i = 2+i0
                  ij = 1+i0+dimp(1)*j0
                  !               (Fa - i/2 Sat) D S - S D (Fa + i/2 Sat)
                  ! +     (F D - i/2 St D - i S Dt) Sa - Sa (D F + i/2 D St + i Dt S)
                  ! + (F Da - i/2 St Da - i/2 S Dat) S - S (Da F + i/2 Da St + i/2 Dat S)
                  FDSp(ij) = FDSp(ij) + (Fp(i) - w(1)/2 * Sp(i)) * DSp(j)
                  FDSp(ij) = FDSp(ij) - SDp(j) * (Fp(i) + w(1)/2 * Sp(i))
                  FDSp(ij) = FDSp(ij) + FDp(j)*Sp(i)
                  FDSp(ij) = FDSp(ij) - Sp(i)*DFp(j)
                  FDSp(ij) = FDSp(ij) + (Fp(j) - (w(1)+w(2))/2 * Sp(j)) * Dp(i) * S0
                  FDSp(ij) = FDSp(ij) - S0 * Dp(i) * (Fp(j) + (w(1)+w(2))/2 * Sp(j))
                  FDSp(ij) = FDSp(ij) + (Fp(1) - w(1)/2 * S0) * Dp(i) * Sp(j)
                  FDSp(ij) = FDSp(ij) - Sp(j) * Dp(i) * (Fp(1) + w(1)/2 * S0)
                  DSDp(ij) = DSDp(ij) + Dp(j) * Sp(i) * Dp(1)
                  DSDp(ij) = DSDp(ij) + Dp(1) * Sp(i) * Dp(j)
                  DSDp(ij) = DSDp(ij) + Dp(i)*SDp(j)
                  DSDp(ij) = DSDp(ij) + DSp(j)*Dp(i)
               end do
               DSp(j)=0; SDp(j)=0; FDp(j)=0; DFp(j)=0; Sp(j)=0
            end do
            Sp(2:1+dimp(1)) = 0
         end if
      else
         call lsquit('pert_scf_eq error: Too high order, size(p) > 2 not implemented',mol%lupri)
      end if
   end subroutine


   !> Solve the p-perturbed response equations for Dp and Fp.
   !> If the solution is already stored in rsp_eq_sol, the solution is fetched
   !> from there. Otherwise, the RHSs are computed before calling the solver
   !> ajt FIXME use type pert_field instead of p,dimp,comp,freq
   !> ajt FIXME remove vectorization to simplify higher-order generalization
   !> ajt FIXME make D,F,Dp,Fp optional, with default in rsp_eq_sol
   subroutine pert_dens(mol, S0, p, dimp, D, F, Dp, Fp, comp, freq)
     implicit none
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(rsp_molcfg), target, intent(inout) :: mol
      !> unperturbed overlap matrix
      type(matrix), intent(in) :: S0
      !> perturbation lables
      character(*), intent(in) :: p(:)
      !> dimension of property (each perturbation)
      integer,      intent(in) :: dimp(:)
      !> lower-than-p'th order density matrices
      type(matrix), intent(in) :: D(:)
      !> lower-than-p'th order Fock matrices
      type(matrix), intent(in) :: F(:)
      !-------------------------------
      !> resulting p-perturbed density. Deferred shape (*) to permit any rank
      type(matrix), intent(inout) :: Dp(*)
      !> resulting p-perturbed Fock. Deferred shape (*) to permit any rank
      type(matrix), intent(inout) :: Fp(*)
      !-----------------------------------
      !> starting component for each p. Default 1 1 ... 1
      integer,        optional, intent(in) :: comp(:)
      !> frequency of each p. Default all zero (static)
      complex(realk), optional, intent(in) :: freq(:)
      !------------------------------------------
      type(pert_field) :: fld(size(p))
      integer        :: ord(size(p)), ccomp(size(p)), pd, pd1, sym, nanti, i, j
      complex(realk) :: ffreq(size(p))
      type(matrix), pointer :: DSDp(:), FDSp(:), Sp(:)
      i=0
      allocate(DSDp(product(dimp)))
      allocate(FDSp(product(dimp)))
      allocate(Sp(product(dimp)))
      if (size(dimp) /= size(p)) call lsquit('pert_dens: Differing numbers ' &
               // 'of perturbations and dimensions, size(dimp) /= size(p)',mol%lupri)
      pd  = product(dimp)
      pd1 = product(dimp+1)
      if (size(D) /= pd1-pd) call lsquit('pert_dens: Wrong number of lower-order ' &
               // 'perturbed densities, size(D) /= product(dimp+1) - product(dimp)',mol%lupri)
      if (size(F) /= pd1-pd) call lsquit('pert_dens: Wrong number of lower-order ' &
               // 'perturbed densities, size(D) /= product(dimp+1) - product(dimp)',mol%lupri)
      ccomp = 1
      if (present(comp)) then
         ccomp = comp !ajt FIXME check number and bounds
      end if
      ffreq = 0
      if (present(freq)) then
         ffreq = freq !ajt FIXME check number
      end if
      ! determine the symmetry of the right-hand-sides
      ! ajt FIXME using sym, pert_scf_eq can be optimized
      nanti = count(pert_antisym(p))
      sym = merge(merge(1, -1, mod(nanti,2)==1), 0, all(ffreq==0) &
                  .or. (size(p)==1 .and. .not.any(pert_basdep(p))))
      ! search through rsp_eq_sol for an already calculated solution
      if (pd==1 .and. allocated(rsp_eq_sol)) then
         ! put fields in standard order, in order to compare with rsp_eq_sol
         fld = (/(pert_field(p(i),ffreq(i),ccomp(i),dimp(i)), &
                  i=1, size(p))/)
         ord = pert_order(fld)
         fld = (/(fld(ord(i)), i=1, size(p))/)
         ! search
         do i=1, size(rsp_eq_sol)
            if (.not.associated(rsp_eq_sol(i)%mol, mol)) cycle
            if (size(rsp_eq_sol(i)%fld) /= size(p)) cycle
            if (any(rsp_eq_sol(i)%fld%label /= fld%label)) cycle
            if (any(rsp_eq_sol(i)%fld%comp /= fld%comp)) cycle
            if (any(rsp_eq_sol(i)%fld%freq /= fld%freq) .and. &
                any(rsp_eq_sol(i)%fld%freq /=-fld%freq)) cycle
            if (all(rsp_eq_sol(i)%fld%freq == fld%freq)) then
               Dp(1) = rsp_eq_sol(i)%D
            else
               Dp(1) = (-1E0_realk)**nanti * trans(rsp_eq_sol(i)%D)
            end if
            ! calculate the corresponding Fock matrix, return
            call pert_fock(mol, p, dimp, (/D,Dp(1)/), Fp(:1), &
                           comp=ccomp, freq=ffreq)
            deallocate(DSDp)
            deallocate(FDSp)
            deallocate(Sp)
            return
         end do
      end if
      ! if an 'EXCI' equation was not caught in the previous
      ! loop, the corresponding excitation density has not
      ! been saved in solved_eqs(:), so quit
      if (size(p)==1 .and. p(1)=='EXCI') &
         call lsquit('rsp_equations/pert_dens error: EXCI specified but no ' &
                // 'excitation densities saved in rsp_equations',mol%lupri)
      ! zero highest-order D and F for RHS calculation
      do i = 1, pd
         Dp(i) = 0E0_realk * D(1)
         Fp(i) = 0E0_realk * F(1)
      end do
      ! if first order FREQ equation, result is zero (so just return)
      if (size(p)==1 .and. p(1)=='FREQ') return
      ! calculate all lower-order contributions to the minus-right-hand-sides
      ! DSDp and FDSp. Do this only when the equation won't be truncated one line below
      if (.not.(size(p) > rsp_eq_truncate_order)) &
         call pert_scf_eq(mol, S0, p, dimp, (/D,Dp(:pd)/), (/F,Fp(:pd)/), FDSp, DSDp, &
                          comp=ccomp, freq=ffreq)

      ! calculate all lower-order contributions to the p-perturbed Fock
      ! matrices, as well as the p-perturbed overlap matrices
      call pert_fock(mol, p, dimp, (/D,Dp(:pd)/), Fp(:pd), Sp, &
                     comp=ccomp, freq=ffreq)

      ! AJT/TK if truncate order is has been set and this equation is above that order,
      ! cleanup DSDp and FDSp and return the zeros already in Dp and the Fp already
      ! containing lower-order contribs
      if (size(p) > rsp_eq_truncate_order) then
         Sp(:) = 0
         deallocate(DSDp)
         deallocate(FDSp)
         deallocate(Sp)
         return
      end if

      ! top DSDp off with DSpD and FDSp with (F-w/2S)DSp-SpD(F+w/2S)
      do i = 1, pd
         !ajt FIXME This will not result in completely (anti-)symmetric
         !          DSDp and FDSp. Use +=sym*trans(..) to achieve this
         DSDp(i) = DSDp(i) + D(1)*Sp(i)*D(1)
         FDSp(i) = FDSp(i)  + (F(1) - sum(ffreq)/2 * S0)*D(1)*Sp(i) &
                 - Sp(i)*D(1)*(F(1) + sum(ffreq)/2 * S0)
      end do
      Sp(:) = 0 !free
      do i=1,pd
         ! KK/AJT quick fix. Avoid solving the same response equation twice for
         ! second-order response equations (i.e. when the two perturbations are identical).
         if (size(p) == 2) then
            ! if the two perturbations are identical
            if (dimp(1) == 3 .and. dimp(2) == 3 &
                  .and. p(1)==p(2) .and. ffreq(1) == ffreq(2)) then
               ! component 2 equals component 4 (similarly for 3,7 and 6,8)
               j = merge(2, merge(3, merge(6, -1, i==8), i==7), i==4)
               if (j/=-1) then
                  DSDp(i) = 0
                  FDSp(i) = 0
                  Dp(i) = Dp(j)
                  Fp(i) = Fp(j)
                  cycle
               end if
            end if
         end if
         ! if no symmetry can be exploited, calculate response parameter
         call rsp_solve_eq(mol, S0, D(1), F(1), sym, sum(ffreq), 1, &
                           DSDp(i:i), FDSp(i:i), Dp(i:i), Fp(i:i))
#ifdef SYS_AIX
         !ajt For some reason Xlf copy-in-copy-out-s Dp(i) and Fp(i), which
         !    will break Dp(i)%init_self_ptr => Dp(i), resulting in memory leak.
         !    Until one figures out how to avoid this, re-associate here
         call re_associate_self_ptrace(Dp(i))
         call re_associate_self_ptrace(Fp(i))
#endif
      end do
      ! if this is a 1-component equation, save the result
      ! ajt FIXME add proper ref counting
      ! ajt FIXME use subroutine for inserting
      if (pd==1) then
         ! make sure there is a vacant slot in rsp_eq_sol
         call rsp_eq_sol_grow(1)
         do i=1, size(rsp_eq_sol)
            if (.not.associated(rsp_eq_sol(i)%mol)) exit
         end do
         ! set reference to molecule + rsp config
         rsp_eq_sol(i)%mol => mol
         ! find standard order of fields, as required in rsp_eq_sol
         fld = (/(pert_field(p(i),ffreq(i),ccomp(i),dimp(i)), &
                  i=1, size(p))/)
         ord = pert_order(fld)
         ! set up field
         allocate(rsp_eq_sol(i)%fld(size(p)))
         rsp_eq_sol(i)%fld = (/(fld(ord(i)), i=1, size(p))/)
         ! copy density. ajt FIXME don't copy, swap and return alias
         rsp_eq_sol(i)%D = Dp(1)
      end if
      print* !ajt Blank line after prints from rsp_solve_eq
      deallocate(DSDp)
      deallocate(FDSp)
      deallocate(Sp)
   end subroutine

   subroutine re_associate_self_ptrace(A)
     type(matrix), intent(inout), target :: A
     A%init_self_ptr => A
   end subroutine

   ! if just moving with B=A, 'ownership' of data in A%D will
   ! not be moved to B%D, resulting in memory leak.
   ! This subroutine transfers ownership too
   subroutine move_eq_sol(A, B)
     type(rsp_saved_sol), target, intent(inout) :: A, B
     B = A
     if (associated(A%D%init_self_ptr, A%D)) then
        B%D%init_self_ptr => B%D
        A%D%init_self_ptr => B%D
     end if
     A%D = 0 !will mat_nullify it
   end subroutine


   !> Grow array storing saved solutions, rsp_eq_sol, so there are n vacancies.
   !> If rsp_eq_sol already has n vacancies, just return
   subroutine rsp_eq_sol_grow(n)
     implicit none
     integer, intent(in) :: n
     type(rsp_saved_sol), allocatable :: tmp(:)
     integer :: i, nvac, prevsiz, newsiz
     ! current size, number of vacancies, return if sufficient
     i=0
     nvac = 0
     prevsiz = 0
     if (allocated(rsp_eq_sol)) then
        nvac = count((/(.not.associated(rsp_eq_sol(i)%mol), &
                        i=1, size(rsp_eq_sol))/))
        prevsiz = size(rsp_eq_sol)
     end if
     if (nvac >= n) return
     ! if allocated, move to tmp
     if (prevsiz /= 0) then
        allocate(tmp(size(rsp_eq_sol)))
        do i = 1, prevsiz
           call move_eq_sol(rsp_eq_sol(i), tmp(i))
        end do
        deallocate(rsp_eq_sol)
     end if
     ! calculate new size, a multiple of 8, allocate
     newsiz = ((prevsiz - nvac + n + 7) / 8) * 8
     allocate(rsp_eq_sol(newsiz))
     ! transfer from tmp, nullify rest
     if (prevsiz /= 0) then
        do i = 1, prevsiz
           call move_eq_sol(tmp(i), rsp_eq_sol(i))
        end do
        deallocate(tmp)
     end if
     do i = prevsiz+1, newsiz
        nullify(rsp_eq_sol(i)%mol)
        nullify(rsp_eq_sol(i)%fld)
        rsp_eq_sol(i)%refc = 0
     end do
   end subroutine


   !> Free stored response equations
   !> ajt FIXME add optional rsp_molcfg argument for masking
   !> ajt FIXME add optional pert_field argument for masking
   subroutine rsp_eq_sol_empty()
     integer :: i
     if (.not.allocated(rsp_eq_sol)) return
     do i=1, size(rsp_eq_sol)
        if (.not.associated(rsp_eq_sol(i)%mol)) cycle
        nullify(rsp_eq_sol(i)%mol)
        deallocate(rsp_eq_sol(i)%fld)
        rsp_eq_sol(i)%D = 0  !free
     end do
     deallocate(rsp_eq_sol)
   end subroutine


   !> Solve an scf response equation by invoking the configured response
   !> solver on the residuals (-RHS'es) in DSDp and FDSp, obtaining
   !> solutions Dp (perturbed density) and Fp (perturbed Fock)
   !> ajt Should at some point also return the final residuals, at least
   !>     FDSp, for use in Sellers-type higher-precision formulas
   !> ajt aug09 Added argument sym, which specifies the symmetry of the minus-
   !>           right-hand-side FDSp: anti-symmetric=-1, symmetric=+1, general=0
   !> ajt sep09 Added argument neq (number of equations), so several can be
   !>           solved in one go
   !> ajt feb10 For linsca, added comparison with excitation energies,
   !>           and projection
   subroutine rsp_solve_eq(mol, S0, D0, F0, sym, freq, neq, DSDp, FDSp, Dp, Fp)
      !       in:          --------------------------------------------
      !       out:                                                      ------
      use RSPsolver, only: rsp_init, rsp_solver
      use RSPsymsolver, only: rsp_sym_init, rsp_sym_solver
      use COMPLEXSOLVER, only: rsp_complex_init, rsp_complex_solver
      use COMPLEXSYMSOLVER, only: rsp_sym_complex_init, rsp_sym_complex_solver
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(rsp_molcfg), intent(inout) :: mol
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> unperturbed density matrix
      type(matrix),      intent(in) :: D0
      !> unperturbed Fock matrix
      type(matrix),      intent(in) :: F0
      !> FDSp anti:-1, symm:+1, general:0
      integer,           intent(in) :: sym
      !> frequency of equation
      complex(realk),    intent(in) :: freq
      !> number of equations
      integer,           intent(in) :: neq
      !> DSD-D residual (-RHS)
      type(matrix),   intent(inout) :: DSDp(neq)
      !> FDS-SDF residual (-RHS)
      type(matrix),   intent(inout) :: FDSp(neq)
      !-----------------------------------------
      !> resulting perturbed density
      type(matrix),   intent(inout) :: Dp(neq)
      !> resulting perturbed Fock
      type(matrix),   intent(inout) :: Fp(neq)
      !-----------------------------------------
      type(matrix)   :: DS
      character*4    :: UNP(0)
      integer        :: i, j
      real(realk)    :: rhs_norm
      type(matrix)   :: Xph(1)
      logical, save  :: first = .true.
      real(realk)    :: freq1(1) !duplicate but array(1)
      ! for complex equations
      logical, save  :: first_complex = .true.
      real(realk)    :: gamma_saved
      type(matrix)   :: reXph(1), imXph(1)
      ! for resonant equations (real freq coincides with +-an excitation energy).
      ! Excitation energies and densities are stored in solved_eqs.
      ! Excitation and de-excitation should be projected out of the
      ! right-hand-sides, and de-/excitation vectors should be passed to the solver
      integer,pointer       :: iexci(:)     !indices in solved_eqs
      type(matrix),   allocatable :: Vexci(:)     !resonant excitation vector(s)
      complex(realk), allocatable :: tsmom(:,:,:) !transition moments

      ! KK fix
      nullify(iexci)
      ! ajt Since we are called from pert_dens only, this check is not really needed...
      if (sym /= -1 .and. sym /= 0 .and. sym /= +1) &
         call lsquit("rsp_solve_eq error: Argument 'sym' is not -1, 0 or 1",mol%lupri)
      ! prepare the right-hand-sides
      DS = D0*S0
      call scf_eq_prep_rhs(mol, D0, DS, sym, neq, FDSp, DSDp, Fp, Dp)
      ! check whether the frequency is resonant (+-excitation energy)
      call scf_eq_find_resonances(mol, freq, iexci)
      ! if resonant, project away the resonant component form FDSp,
      ! while storing the +freq and -freq transition moments in tsmom
      if (associated(iexci)) then
         allocate(Vexci(size(iexci)))
         allocate(tsmom(2, neq, size(iexci)))
         do i=1, size(iexci)
            call proj_resonant_rhs(S0, DS, sym, freq, neq,           &
                                   rsp_eq_sol(iexci(i))%fld(1)%freq, &
                                   rsp_eq_sol(iexci(i))%D,           &
                                   FDSp, Fp, Vexci(i), tsmom(:,:,i))
         end do
      end if
      ! call solver
      do i=1, neq
         rhs_norm = norm(FDSp(i))
         ! in case norm defined differently, be slightly more conservative than solver (/10)
         print *, 'before response solver: norm(RHS) = ', rhs_norm
         if (rhs_norm < mol%solver%rsp_thresh / 10) then
            print *, '=> skipping this equation'
            Xph(1) = 0E0_realk * S0
         else if (FDSp(i)%complex .or. abs(mol%solver%rsp_gamma) > 1e-20 &
                 & .or. abs(aimag(freq)) > 1e-20 ) then
            !on the first run, initialize the solvers
            if (first_complex) then
              if (mol%solver%rsp_stdnew) then
               call rsp_sym_init(2, 2, 2, 1, 2)
              endif
              if (mol%solver%rsp_cmplxnew) then
               call rsp_sym_complex_init(1, 1, 1, 1, 1)
              endif
               IF(mol%solver%UseExcitationVecs)then
                call rsp_init(2, 2, 2, 1, 2+mol%solver%rsp_eigenvecs)
               ELSE
                call rsp_init(2, 2, 2, 1, 2)
               ENDIF
                call rsp_complex_init(1, 1, 1, 1, 1)
               first = .false.
               first_complex = .false.
            end if
            ! if -RHS FDSp isn't complex, add a zero complex matrix to it
            if (.not.FDSp(i)%complex) &
               FDSp(i) = FDSp(i) + (0E0_realk,1E0_realk)*(0*FDSp(i))
            !configure the frequency and damping factor
            freq1(1) = real(freq)
            gamma_saved = mol%solver%rsp_gamma !save, for restore after solver
            mol%solver%rsp_gamma = mol%solver%rsp_gamma + aimag(freq)
            ! the solver currently doesn't seem to work with zero gamma
            if (abs(mol%solver%rsp_gamma) < 1d-6) &
               mol%solver%rsp_gamma = 2d-6
            ! allocate the solution matrix Xph (Xp/2)
            Xph(1) = 0*FDSp(i)
            call mat_ensure_alloc(Xph(1), only_alloc=.true.)
            !create aliases for the real and imaginary parts of FDSp and Xph
            reXph(:) = (/mat_get_part(Xph(1), imag=.false.)/)
            imXph(:) = (/mat_get_part(Xph(1), imag=.true. )/)
            ! solve for the real and imaginary part of the FDSp simlultaneously
            if (mol%solver%rsp_cmplxnew) then
               call rsp_sym_complex_solver(mol, F0, D0, S0, 1,          &
                               (/mat_get_part(FDSp(i), imag=.false.)/), &
                               freq1, 1, reXph(:), imXph(:), .true.,    &
                               gdi=(/mat_get_part(FDSp(i), imag=.true.)/))
            else
               call rsp_complex_solver(mol, F0, D0, S0, 1,               &
                               (/mat_get_part(FDSp(i), imag=.false.)/),  &
                               freq1, 1, reXph(1:1), imXph(1:1), .true., &
                               gdi=(/mat_get_part(FDSp(i), imag=.true.)/))
            endif 
            !clear re and im aliases of FDSp and Xph
            reXph = 0
            imXph = 0
            !restore gamma to avoid breaking other code
            mol%solver%rsp_gamma = gamma_saved
         else
            if (first) then
              if ((mol%solver%rsp_stdnew) .and. (abs(real(freq)) .gt.  1E-6_realk)) then
               call rsp_sym_init(1, 1, 1, 1, 1)
              endif
            !to be fixed when new solver contains projections
               call rsp_init(1, 1, 1, 1, 1)
               first = .false.
            end if
            freq1(1) = real(freq)
            Xph(1) = 0*Dp(i)
            call mat_ensure_alloc(Xph(1))
            ! if non-resonant, no excitation vectors to pass
            if (.not.associated(iexci)) then
               if (mol%solver%rsp_stdnew) then
                  if (abs(real(freq)) .gt.  1E-6_realk) then
                     call rsp_sym_solver(mol, D0, S0, F0, .true., 1, FDSp(i:i), &
                                               freq1(1:1), Xph(1))
                  else
                     call rsp_solver(mol, D0, S0, F0, .true., 1, FDSp(i:i), &
                                     freq1(1:1), Xph(1))
                  endif
               else
               call rsp_solver(mol, D0, S0, F0, .true., 1, FDSp(i:i), &
                               freq1(1:1), Xph(1))
               endif
            else !if resonant, pass excitation vectors Vexci
               call rsp_solver(mol, D0, S0, F0, .true., 1, FDSp(i:i), &
                               freq1(1:1), Xph(1), Xproject=Vexci)
            end if
         end if
         ! if (anti-)symmetric Dp (static with anti-/symmetric FDSp,DSDp),
         ! make sure Dp it comes out completely symmetric
         if (sym /= 0 .and. freq==0 .and. .not.mol%solver%rsp_complex) then
            Dp(i) = 1/2E0_realk * Dp(i) + 2E0_realk * DS*Xph(1)
            Dp(i) = Dp(i) - 1E0_realk*sym * trans(Dp(i))
         else
            Dp(i) = Dp(i) + 2E0_realk * (DS*Xph(1) - Xph(1)*trans(DS))
         end if
         Xph(1) = 0
         FDSp(i) = 0
      end do
      ! if this is a projected resonant equation, add contribution from
      ! the opposite (de-)excitation (then deallocate)
      if (associated(iexci)) then
         Vexci(:) = 0 !free excitation vectors
         deallocate(Vexci)
         deallocate(iexci)
         deallocate(tsmom)
      end if
      DS = 0
      ! add last contribution G(Dp) to Fp
      call rsp_twoint(mol, UNP, (/D0, Dp/), (/neq/), Fp)
#ifdef SYS_AIX
      !ajt Xlf for some reason does copy-in-copy-out of Dp and Fp, which will
      !    offset type(matrix)%init_self_ptr, resulting in memory-leak.
      !    So nullify here, and re-associate after return to pert_dens
      do i = 1,neq
         nullify(Dp(i)%init_self_ptr)
         nullify(Fp(i)%init_self_ptr)
      end do
#endif
   end subroutine


   !> Private subroutine for preparing right-hand-sides before calling the solver
   !> Calculates particular component of Dp from DSDp, adds contribution
   !> from that to FDSp, and projects out the particular component from FDSp
   !> (this is for numerical stability). Precalculated DS is taken as input.
   subroutine scf_eq_prep_rhs(mol, D0, DS, sym, neq, FDSp, DSDp, Fp, Dp)
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(rsp_molcfg),  intent(inout) :: mol
      type(matrix),      intent(in)    :: D0, DS, Fp(neq)
      integer,           intent(in)    :: sym, neq
      type(matrix),      intent(inout) :: FDSp(neq), DSDp(neq), Dp(neq)
      character*4 :: UNP(0)
      integer     :: i
      do i=1, neq
         ! if the minus-right-hand-sides are (anti-)symmetric, one gets away
         ! with half the number of multiplications (though twice as many lines)
         if (sym /= 0) then
            Dp(i) = 1/2E0_realk * DSDp(i) - DS*DSDp(i)
            Dp(i) = Dp(i) - sym * 1E0_realk * trans(Dp(i)) !symmetrize or anti-symmetrize
            FDSp(i) = FDSp(i)*DS
            FDSp(i) = FDSp(i) - sym * 1E0_realk * trans(FDSp(i)) + Fp(i)
            ! ajt Ensure that the perturbed Fock is precisely (anti-)symmetric
            !ajt wrong: Fp(i) = 1/2d0 * Fp(i) - 1/2d0*sym * trans(Fp(i)) !(anti-)symmetrize
         else
            Dp(i) = DSDp(i) - DSDp(i)*trans(DS) - DS*DSDp(i)
            FDSp(i) = FDSp(i)*DS - trans(DS)*FDSp(i) + Fp(i)
         end if
         DSDp(i) = 0
      end do
      ! contract Coulomb-exchange+Kohn-Sham with particular components Dp(;)
      call rsp_twoint(mol, UNP, (/D0, Dp/), (/neq/), FDSp)
      ! complete right-hand-sides
      do i=1, neq
         if (sym /= 0) then
            FDSp(i) = FDSp(i)*DS
            FDSp(i) = FDSp(i) + 1E0_realk * sym * trans(FDSp(i))
         else
            FDSp(i) = FDSp(i)*DS - trans(DS)*FDSp(i)
         end if
      end do
   end subroutine


   !> Look through rsp_eq_sol for whether freq is an excitation energy.
   !> If yes, allocate iexci(:), and fill it with indices in rsp_eq_sol
   subroutine scf_eq_find_resonances(mol, freq, iexci)
      !ajt FIXME use response_wrapper_type_module, only: RSPSOLVERinputitem
      type(rsp_molcfg),  target, intent(in)  :: mol
      complex(realk),            intent(in)  :: freq
      integer, pointer :: iexci(:)
      real(realk) :: dif, mindif, thr
      integer     :: deg, i, iter, imin
      if (.not.allocated(rsp_eq_sol)) return
      !DTHR is a threshold for when states are considered degenerate.
      !Here used slightly different, to determine when the rsp equation
      !is singular, but it makes sense to use same threshold, right. TK
      thr = mol%solver%degenerateTHR
      do iter = 1,2    !iter=1 counts degeneracy, then allocates iexci
         deg = 0       !iter=2 records indices in iexci
         mindif = huge(1E0_realk)
         do i=1, size(rsp_eq_sol)
            ! wrong molecule, skip
            if (.not.associated(rsp_eq_sol(i)%mol, mol)) cycle
            ! not first order, skip
            if (size(rsp_eq_sol(i)%fld) /= 1) cycle
            ! wrong equation, skip
            if (rsp_eq_sol(i)%fld(1)%label /= 'EXCI') cycle
            ! compute difference, keep minimum
            dif = min(abs(rsp_eq_sol(i)%fld(1)%freq - freq), &
                      abs(rsp_eq_sol(i)%fld(1)%freq + freq))
            if (dif < mindif) imin = i
            if (dif < mindif) mindif = dif
            ! not resonant, skip
            if (dif > thr) cycle
            deg = deg+1
            if (iter==2) iexci(deg) = i
         end do
         if (iter==1 .and. deg==0) exit !no resonant excitations
         ! several excitations are close to +-freq, but the solver/input
         ! does not indicate that any are degenerate.
         ! Then only pick the one closest to +-freq
         if (iter==1 .and. deg/=0 .and. &
              .not.mol%solver%degeneratestates) then
            allocate(iexci(1))
            iexci(1) = imin
            exit
         end if
         if (iter==1 .and. deg/=0) allocate(iexci(deg))
      end do
   end subroutine


   ! Project the excitation densities out of the -right-hand-sides FDSp,
   ! and return the corresponding excitation vector Vexci and
   ! transition moments tsmom. tsmom(1,:) are de-excitation moments (-exci),
   ! while tsmom(2,:) are excitation moments (+exci).
   ! ajt may10 Added forgotten projection of Fp
   subroutine proj_resonant_rhs(S0, DS, sym, freq, neq, &
                                exci, Dx, FDSp, Fp, Vx, tsmom)
      type(matrix),   intent(in)    :: S0, DS
      integer,        intent(in)    :: sym, neq
      complex(realk), intent(in)    :: freq, exci
      type(matrix),   intent(in)    :: Dx
      type(matrix),   intent(inout) :: FDSp(neq), Fp(neq), Vx
      complex(realk), intent(out)   :: tsmom(2,neq)
      type(matrix) :: SDxS, SVxS
      integer      :: i, j, k
      Vx = Dx*trans(DS) - DS*Dx !for dot/trace with FDSp
      SDxS = S0*Dx*S0           !for removal from FDSp
      SVxS = S0*Vx*S0           !for removal from Fp
      ! loop over -right-hand-sides (equations)
      do i=1, neq
         ! excitation transition moment (frequency +exci)
         tsmom(2,i) = dot(Vx, FDSp(i))
         ! de-excitation transition moment (frequency -exci)
         if (sym/=0) tsmom(1,i) = sym * tsmom(2,i)
         if (sym==0) tsmom(1,i) = trace(Vx, FDSp(i))
         ! remove resonant component from Fp, so Fp/Dp will
         ! solve the response equation FpDS+FDpS...-wSDpS=0
         if (abs(exci+freq) < abs(exci-freq)) then
            Fp(i) = Fp(i) + tsmom(1,i) * trans(SVxS)
         else
            Fp(i) = Fp(i) - tsmom(2,i) * SVxS
         end if
      end do
      SDxS = 0
      SVxS = 0
   end subroutine


end module
