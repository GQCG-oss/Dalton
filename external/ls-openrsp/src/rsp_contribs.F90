! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module rsp_contribs

!> This module contains routines for calculating contributions
!> to molecular properties (1st order, linear response, etc.),
!> and perturbed Fock matrices. 
module lsdalton_rsp_contribs
  use precision
  use integraloutput_type,   only: initIntegralOutputDims
  use ls_Integral_Interface, only: ls_getIntegrals, &
       ls_attachDmatToSetting, &
       ls_freeDmatFromSetting
  use TYPEDEFTYPE, only: LSSETTING
  use TYPEDEF, only: retrieve_output
  use GCtransMod, only: GCAO2AO_transform_matrixD2
  use Integralparameters
  use lsdalton_matrix_defop
  use RSPsolver, only: rsp_molcfg
  use dal_interface , only: di_GET_GbDs, di_GET_GbDs_and_XC_linrsp
  use matrix_module, only: Matrixp, matrix
  use matrix_operations, only: mat_init, mat_assign, mat_free
  use integralinterfaceMod
  use II_XC_interfaceModule
   implicit none
   public rsp_oneave
   public rsp_twoave
   public rsp_oneint
   public rsp_twoint
!   public rsp_molcfg
   public pert_basdep
   public pert_shape
   public pert_antisym
   public pert_field
   public pert_order
   private


   !> Struct denoting a perturbing field in a response function
   !> or response equation. A response equation (or density)
   !> corresponds to an array of pert_field.
   !> Analogously, a response function corresponds to an array
   !> of pert_field whose frequencies (freq) sum to zero.
   type pert_field
      sequence
      !> 4-char pert label
      character(4)   :: label
      !> frequency
      complex(realk) :: freq
      !> first component
      integer        :: comp
      !> number of components
      integer        :: ncomp
   end type


   !> private struct to collect properies of pertrubing fields
   type pert_field_info
      !> four-letter abbreviation
      character(4)  :: label
      !> long name
      character(64) :: name
      !> number of components (when known, 0 otherwise)
      integer       :: ncomp
      !> anti-symmetric (1,3,5th ord.) perturbed integrals
      logical       :: anti
      !> basis dependent (sa. GEO and MAG)
      logical       :: bas
      !> one-electron operator linear in field strength (EL)
      logical       :: lin
      !> one-electron operator quadratic in field strength (MAGO)
      logical       :: quad
   end type


   ! to compactify the table
   logical, parameter :: T = .true.
   logical, parameter :: F = .false.


   !> ajt jan10: EXCI is a ZERO (no) perturbation, and is introduced to
   !>            allow the same code to contract response functions and
   !>            "generalized transition moments".
   !> ajt may10: FREQ is also a ZERO (no) perturbation, and is introduced to
   !>            allow the same code to contract response functions and
   !>            frequency-differentiated response functions.
   type(pert_field_info) :: field_list(12) = &                         !nc an ba ln qu
      (/pert_field_info('EXCI', 'Generalized "excitation" field'      , 1, F, F, T, T), &
        pert_field_info('FREQ', 'Generalized "freqency" field'        , 1, F, F, T, T), &
        pert_field_info('EL  ', 'Electric field'                      , 3, F, F, T, F), &
        pert_field_info('VEL ', 'Velocity'                            , 3, T, F, T, F), &
        pert_field_info('MAGO', 'Magnetic field w/o. London orbitals' , 3, T, F, F, T), &
        pert_field_info('MAG ', 'Magnetic field with London orbitals' , 3, T, T, F, F), &
        pert_field_info('ELGR', 'Electric field gradient'             , 6, F, F, T, F), &
        pert_field_info('VIBM', 'Displacement along vibrational modes',-1, F, T, F, F), & !-1=mol-dep
        pert_field_info('GEO ', 'Nuclear coordinates'                 ,-1, F, T, F, F), & !-1=mol-dep
        pert_field_info('NUCM', 'Nuclear magnetic moment'             ,-1, T, F, F, T), & !-1=mol-dep
        pert_field_info('AOCC', 'AO contraction coefficients'         ,-1, F, T, F, F), & !-1=mol-dep
        pert_field_info('AOEX', 'AO exponents'                        ,-1, F, T, F, F)/)  !-1=mol-dep

   character(1), parameter :: xyz(3) = (/'X','Y','Z'/)

#ifndef PRG_DIRAC
   logical, external :: do_dft
#endif /* PRG_DIRAC */

contains


   function prefix_zeros(n, l)
      integer, intent(in) :: n, l !number, length
      character(l)        :: prefix_zeros !resulting n in ascii
      character(1), parameter :: char0to9(0:9) &
            = (/'0','1','2','3','4','5','6','7','8','9'/)
      integer :: i, k
      k = n
      do i = l, 1, -1
         prefix_zeros(i:i) = char0to9(mod(k,10))
         k = k / 10
      end do
      if (k /= 0) call lsquit('prefix_zeros error: Argument integer does not fit ' &
                           // 'in the specified number of ASCII caracters',-1)
   end function


   !> Contracts the 1-electron integrals perturbed by the perturbations p(:)
   !> with the perturbed density matrices in D(:) (e.g. D=(/Dxy/) for a 2nd
   !> order density), and ADDS the result to the property (rsp func) array
   !> E(:). Front for the private subroutine 'oneave' below, checking the
   !> arguments' dimensions, and doing permutations
   !> S0 is passed as argument only as reference to nuclei and basis set
   subroutine rsp_oneave(mol, S0, p, D, dime, E, perm, comp, freq, DFD)
     implicit none
      !> structure containing integral program settings
      type(rsp_molcfg),  intent(inout) :: mol
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> p(np) perturbation lables
      character(*),      intent(in) :: p(:)
      !> (un)perturbed density matrices to
      !> contract perturbed one-electron integrals with.
      !> If perm present, size(D) = product(dime(perm(np+1:np+nd))),
      !> if perm not present, size(D) = product(dime(np+1:np+nd))
      type(matrix),      intent(in) :: D(:)
      !> dime(np+nd) = shape(E), dimensions of property
      integer,           intent(in) :: dime(:)
      !-------------------------------------------------------------
      !> property contributions, works incrementally, thus
      !> contributions are ADDED to E(*). size(E) = product(dime)
      complex(realk),    intent(inout) :: E(*) 
      !-------------------------------------------------------------
      !> perm(np+nd), permutation of indices.
      !> For each dimension of p and D, the corresponding dimension in E.
      !> Default 1 2 ... np+nd (no permutation)
      integer,        optional, intent(in) :: perm(:)
      !> comp(np), starting component index for each p.
      !> Default 1 1 ... 1
      integer,         optional, intent(in) :: comp(:)
      !> freq(np), complex frequencies
      !> for each p, default all zero. Multiply the half-derivative
      !> overlap integrals, thus no contribution if basis independent of p
      complex(realk), optional, intent(in) :: freq(:)
      !> optional perturbed energy-weighted density matrices.
      !> Contracted against perturbed overlap integrals
      type(matrix),   optional, intent(in) :: DFD(:) 
      !-------------------------------------------------------------
      integer        :: idxp(size(p)), pperm(size(dime)), ccomp(size(p)), &
                        stepe(size(dime)), ddime(size(dime)), idxe(size(dime)), &
                        i, j, k, tmpi, nd
      character(4)   :: pp(size(p)), tmpp
      complex(realk) :: ffreq(size(p)), tmpf
      logical        :: zero, bas
      complex(realk), pointer :: Etmp(:)
      i=0
      allocate(Etmp(product(dime)))

      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero. If there is more than 1
      ! linear perturbation (like 'EL') or more than two quadratic (like 'MAGO')
      zero = .false.
      j = 0
      k = 0
      do i = 1, size(p)
         if (.not.field_list(idxp(i))%lin .and. &
             .not.field_list(idxp(i))%quad) cycle
         ! if both linear and quadratic, interpret as ZERO (like EXCI)
         if (field_list(idxp(i))%lin .and. &
             field_list(idxp(i))%quad) zero = .true.
         if (j /= 0) then !if second lin or quad perturbation
            ! if j was EL, then i is second EL (or first MAGO)
            if (field_list(idxp(j))%lin) zero = .true.
            ! if j was MAGO, and different from i, also zero
            if (field_list(idxp(j))%quad .and. &
               idxp(i) /= idxp(j)) zero = .true. !two different quadratic
            !third linear or quadratic perturbation
            if (k /= 0) zero = .true.
            k = j
         end if
         j = i
      end do
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dime) < size(p)) call lsquit('rsp_oneave argument error: ' &
               // 'More perturbations than dimensions of property, size(dime) < size(p)', mol%lupri)
      if (any(dime <= 0)) call lsquit('rsp_oneave argument error: ' &
               // 'Property has a zero or negative dimension, dime <= 0', mol%lupri)
      ! compute step lengths in E (cumulative products of dimensions)
      stepe(1) = 1
      do i = 2, size(dime)
         stepe(i) = stepe(i-1)*dime(i-1)
      end do
      ! reorder dimensions and step lengths in E according to permutation argument perm
      ddime = dime
      if (present(perm)) then
         if (size(perm) /= size(dime)) call lsquit('rsp_oneave argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dime)', mol%lupri)
         ! verify that perm is indeed a permutation
         do i = 1, size(dime)-1
            if (perm(i) <= 0 .or. perm(i) > size(dime) .or. &
                any(perm(i) == perm(i+1:size(dime)))) call lsquit('rsp_oneave ' &
                      // 'argument error: Permutation must contain each number exactly once', mol%lupri)
         end do
         ddime = (/( dime(perm(i)), i=1,size(dime))/)
         stepe = (/(stepe(perm(i)), i=1,size(dime))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call lsquit('rsp_oneave argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)', mol%lupri)
         if (any(comp <= 0)) call lsquit('rsp_oneave argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0', mol%lupri)
         ccomp = comp
      end if
      if (any(ccomp + ddime(:size(p)) - 1 > pert_shape(mol,p))) &
         call lsquit('rsp_oneave argument error: Lowest component index plus ' &
                // 'dimension exceeds dimension of perturbation, comp + dime > pert_shape(mol,p)', mol%lupri)
      ! check optional argument freq, default to zero
      ffreq = 0
      if (present(freq)) then
         if (size(freq) /= size(p)) call lsquit('rsp_oneave ' &
               // 'argument error: Wrong number of frequencies, size(freq) /= size(p)', mol%lupri)
         ffreq = freq
      end if
      ! if unperturbed density and anti-symmetric integral, also zero
      j = count((/(field_list(idxp(i))%anti, i=1,size(p))/))
      if (size(p)==size(dime) .and. mod(j,2)==1) then
         ! if not all bas-dep, or all bas-dep and zero frequencies
         if (.not.all((/(field_list(idxp(i))%bas,i=1,size(p))/)) &
             .or. all(ffreq==0)) zero = .true.
      end if
      ! sort perturbations p so that idxp is descending
      pp = p
      do i = 1, size(p)
         !find perturbation after i, with highest idxp (secondly highest ddime)
         j = i
         do k = j+1, size(p)
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddime(k) > ddime(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddime(i);  ddime(i) = ddime(j);  ddime(j) = tmpi
         tmpi = stepe(i);  stepe(i) = stepe(j);  stepe(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
         tmpf = ffreq(i);  ffreq(i) = ffreq(j);  ffreq(j) = tmpf
      end do
      ! verify that we have the correct number of perturbed densities
      nd = product(ddime(size(p)+1:size(dime)))
      if (size(D) /= nd) call lsquit('rsp_oneave error: Number of' &
               // 'perturbed densities D does not correspond to dime (and perm)', mol%lupri)
      ! verify number of DFD, and that all are defined
      bas = all((/(field_list(idxp(i))%bas, i=1,size(p))/))
      if (present(DFD)) then
         if (size(DFD) /= nd) call lsquit('rsp_oneave error: Number of' &
               // 'perturbed DFD differs from number of perturbed densities D', mol%lupri)
         ! if no basis perturbation (or zero perturbed integrals,
         ! DFD will not be used, so verify they are defined
         if (.not.bas .or. zero) then
            if (.not.all((/(isdef(DFD(i)), i=1,nd)/))) &
               call lsquit('rsp_oneave error: Undefined matrix in argument DFD(:)', mol%lupri)
         end if
      end if
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call lsquit('rsp_oneave error: Undefined matrix in argument D(:)', mol%lupri)
         return
      end if
      ! everything set up, so call core procedure oneave.
      ! Argument nd=0 is used when averaging over unperturbed density,
      ! in which case also perturbed nuclear attraction should be included
      nd = merge(0, nd, size(p)==size(dime))
      call oneave(mol, S0, size(p), pp, ccomp, ddime(:size(p)), ffreq, nd, D, Etmp, DFD)
      ! add oneavg property contribution in temporary array Etmp(:) to
      ! resulting property array E(:), while permuting indices according to pperm
      idxe = 0
      do j = 1, product(dime)
         i = 1 + sum(idxe * stepe)
         E(i) = E(i) + Etmp(j)
         do k = 1, size(dime)
            idxe(k) = idxe(k) + 1
            if (idxe(k) /= ddime(k)) exit
            idxe(k) = 0
         end do
      end do
      deallocate(Etmp)
   end subroutine


   subroutine oneave(mol, S0, np, p, c, dp, w, nd, D, E, DFD)
     implicit none
      !> structure containing integral program settings
      type(rsp_molcfg),  intent(inout):: mol 
      !> unperturbed overlap, to know its dimension
      type(matrix),      intent(in)  :: S0
      !> number of perturbations and order of density
      integer,           intent(in)  :: np
      !> perturbation lables
      character(*),      intent(in)  :: p(np)
      !> lowest component of each perturbation
      integer,           intent(in)  :: c(np)
      !> dimensions of property integrals
      integer,           intent(in)  :: dp(np)
      !> frequency of each p
      complex(realk),    intent(in)  :: w(np)
      !> dimensions of property integrals
      integer,           intent(in)  :: nd
      !> un-/perturbed density matrices,
      !> size(D) = product(1+de(np+1:np+nd))
      type(matrix),      intent(in)  :: D(max(1,nd))
      !-----------------------------------------------------
      !> resulting one-electron property contributions,
      !> size(E) = product(dp) * nd
      complex(realk),    intent(inout) :: E(*)
      !> un-/perturbed energy-weighted density matrices,
      !> Should have size(DFD) = size(D)
      type(matrix), optional, intent(in) :: DFD(:)
      !-----------------------------------------------------
      type(matrix), allocatable :: AA(:) !scratch matrices
      real(realk),  allocatable :: RR(:) !scratch
      type(matrix) :: A(9) !scratch matrices
      real(realk)  :: R(6) !scratch
      real(realk),allocatable :: expval(:)
      integer      :: i, j, k, l, ii, jj, kk, ll, na

      i=0
      na = mol%natoms
      A(1) = 0*S0
      A(2:) = (/(A(1), i=2,size(A))/) !scratch matrices
      if (np==0) then
         call lsquit('rsp_oneave error: unperturbed one-electron contribution requested', mol%lupri)
      else if (np==1 .and. p(1)=='EL  ') then
         ! contract -dipole integrals 'DIPLEN ' with densities
         ! direct calculation of oneel integrals in LSDALTON
         do i=1,3
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         call II_get_prop(mol%lupri, mol%luerr, mol%setting, A(:3), 3, 'DIPLEN ')
         do i=0, dp(1)-1
            do j=0, max(0,nd-1)
               E(1+i+dp(1)*j) = trace(A(c(1)+i),D(1+j))
            end do
         end do
         ! 'no' densities means D(1) contains unperturbed density matrix.
         ! Then also add nuclear attraction contribution to -dipole moment
         if (nd==0) then
            call prop_intifc_nuc(mol, p, R(:3))
            E(:dp(1)) = E(:dp(1)) + R(c(1):c(1)+dp(1)-1)
         end if
      else if (np==1 .and. p(1)=='MAGO' .and. nd/=0) then
         do i=1,3
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         call II_get_prop(mol%lupri, mol%luerr, mol%setting, A(:3), 3, 'ANGMOM ')
         do i=0, dp(1)-1
            do j=0, nd-1
               E(1+i+dp(1)*j) = trace(A(c(1)+i),D(1+j)) / 2
            end do
         end do
         A(:3) = 0
      else if (np==1 .and. p(1)=='MAG ' .and. nd/=0) then
         if (w(1)/=0) &
            call lsquit('rsp_oneave error: MAG with freq/=0 not implemented', mol%lupri)
         if (.not.present(DFD)) &
            call lsquit('rsp_oneave error: MAG must have matrix DFD present', mol%lupri)
         do i=1,6
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         call II_get_prop(mol%lupri, mol%luerr, mol%setting, A(:3), 3, 'MAGMOM ')
         ! ajt to Thomas: what we really need here is half-differentiated overlap
         ! Thomas to ajt: the call II_get_magderivOverlapR or II_get_magderivOverlapL 
         ! depending on ket or bra differentiation
         call II_get_magderivOverlap(A(4:6), mol%setting, mol%lupri, mol%luerr)
         do i=0, dp(1)-1
            do j=0, max(0,nd-1)
               E(1+i+dp(1)*j) = trace(A(c(1)+i),D(1+j)) - trace(A(c(1)+i+3),DFD(1+j))
            end do
         end do
         A(:6) = 0
      else if (np==1 .and. p(1)=='ELGR') then
         ! contract -quadrupole integrals 'THETA' (preferred) or 'SECMOM' with densities.
         do i=1,6
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         ! ajt-to-Thomas: Replace with direct call to II_get_prop, placing THETA/SECMOM
         !                in A(:6). (note that args A(:1) and R(:1) are not used in this
         !                call to prop_intifc_1el)
         ! Thomas-to-ajt: The correct call is
         ! call II_get_prop(mol%lupri, mol%luerr, mol%setting, A(1:6), 6, 'THETA  ')           
         call lsquit('Andreas need to put in some code.',mol%lupri)
!         call prop_intifc_1el(mol, p, A(:1), R(:1), A(:6))
!         do i=0, dp(1)-1 !ELGR
!            do j=0, max(0,nd-1) !D
!               E(1+i+dp(1)*j) = trace(A(c(1)+i),D(1+j))
!            end do
!         end do
         ! 'no' densities means D(1) contains unperturbed density matrix.
         ! Then also add nuclear attraction contribution to -quadrupole moment
!         if (nd==0) then
!            ! ajt-to-Thomas: Replace with direct call to II_get_* nuclear contribution to -quadrupole moment
!            call prop_intifc_nuc(mol, p, R(:3))
!            E(:dp(1)) = E(:dp(1)) + R(c(1):c(1)+dp(1)-1)
!         end if
!         Thomas-to-ajt: the call is 
!         call II_get_nuc_Quad(Q,mol%SETTING,mol%lupri,mol%luerr)
!         where Q is a (3,3) matrix of real(realk)
!         Note that there is a factor 2 difference between this 
!         Routine and the one in NUCQDR
!         Since I do not know the theory and only implemented according 
!         to you specifications you must figure out if I should include a factor
!         2 in my routine or what?  
         A(:6) = 0
      else if (np==1 .and. p(1)=='GEO ') then
         ! one-electron integral contribution
         allocate(RR(3*na))
         do j=0, max(0,nd-1)
            ! ajt-to-Thomas: Replce with direct call to get_oneElectron_gradient, placing
            !                gradient contribution in RR(:3*na), including reorthonormalization
            !                contribution: tr Hg D(1+j) - tr Sg DFD(1+j)
            ! Thomas-to-ajt: If you know what to do then we the fuck not do it
            call prop_intifc_1el(mol, p, (/D(1+j),DFD(1+j)/), RR, A(:1))
            E(1+dp(1)*j:dp(1)*(j+1)) = RR(c(1):c(1)+dp(1)-1)
         end do
         ! nuclear repulsion contribution to unperturbed:nd=0 gradient
         if (nd==0) then
            call prop_intifc_nuc(mol, p, RR)
            E(:dp(1)) = E(:dp(1)) + RR(c(1):c(1)+dp(1)-1)
         end if
      else if (np==1 .and. p(1)=='NUCM' .and. nd/=0) then
         allocate(expval(nd*3*mol%natoms))
         call II_get_prop_expval(mol%lupri,mol%luerr,mol%setting,expval,D(1:nd),nd,3*mol%natoms ,'PSO    ')
         do j=0, nd-1 !D
            do i=0, dp(1)-1 !NUCM
               E(1+i+dp(1)*j) = -2*expval(1+i+dp(1)*j)
            end do
         end do
         deallocate(expval)
      else if (np==2 .and. all(p==(/'MAG','EL '/)) .and. nd /= 0) then
         do i=1,9
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         call II_get_prop(mol%lupri,mol%luerr,mol%setting,A,9,'D-CM1   ') 
         do k = 0, nd-1
            do j = 0, dp(2)-1 !EL indices
               do i = 0, dp(1)-1 !MAG indices
                  E(1+i+dp(1)*(j+dp(2)*k)) = trace(A(c(1)+i+(c(2)+j-1)*3), D(k+1))
               enddo
            enddo
         enddo
         A(:9) = 0
      else if (np==2 .and. p(1)=='NUCM' .and. p(2)=='MAG ') then
         !Thomas-to-ajt : you properly want to generalize this to number of density matrices greater than 1. TK
         allocate(expval(9*mol%natoms))
         ! order (3,3,NATOM) (magnetic X,Y,Z, atomic X,Y,Z, for each Atom
         call II_get_prop_expval(mol%lupri,mol%luerr,mol%setting,expval,D(1),1,9*mol%natoms,'NST    ')
         do l=1,mol%natoms !natoms
            do j=1,3       !nuc
               do i=1,3    !mag
                  E(j+(l-1)*3+(i-1)*3*mol%nAtoms) = 2*expval(i + (j-1)*3 + (l-1)*3*3)
               enddo
            enddo
         enddo
         deallocate(expval)
      else
         print *,'rsp_oneave: no integrals for these perturbations: ', &
                 (p(i)//' ', i=1,np)
         call lsquit('rsp_oneave: no such integrals', mol%lupri)
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine



   !> Contract the 2-electron integrals perturbed by the perturbations p(:)
   !> with the perturbed density matrix expansion in D(:) (e.g. D=(/D,Dx,Dy,Dxy/)
   !> for a 2nd order expansion), and ADDS the result to the property (rsp func)
   !> array E(:). Front for the private subroutine 'twoave' below, checking the
   !> arguments' dimensions, and doing permutations
   !> D(1) serves as reference to nuclei, basis and model/functional
   subroutine rsp_twoave(mol, p, D, dime, E, perm, comp)
     implicit none
      !> structure containing integral program settings
      type(rsp_molcfg),  intent(in) :: mol
      !> p(np) perturbation lables
      character(*),      intent(in) :: p(:)
      !> (un)perturbed density matrices to contract perturbed one-electron
      !> integrals with. If perm present, size(D) = product(dime(perm(np+1:np+nd))),
      !> if perm not present, size(D) = product(dime(np+1:np+nd))
      type(matrix),      intent(in) :: D(:)
      !> dime(np+nd) = shape(E), dimensions of property
      integer,           intent(in) :: dime(:)
      !------------------------------------------------------------------------
      !> property contributions, works incrementally, thus contributions are
      !> ADDED to E(*). size(E) = product(dime)
      complex(realk),    intent(inout) :: E(*) !
      !------------------------------------------------------------------------
      !> perm(np+nd), permutation of indices. For each dimension of p and D,
      !> the corresponding dimension in E. Default 1 2 ... np+nd (no permutation)
      integer, optional, intent(in) :: perm(:)
      !> comp(np), starting component index for each p. Default 1 1 ... 1
      integer, optional, intent(in) :: comp(:)
      !------------------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dime)), ccomp(size(p)), &
                      stepe(size(dime)), ddime(size(dime)), &
                      idxe(size(dime)), i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      logical      :: zero
      complex(realk), pointer :: Etmp(:)
      integer      :: sizep
      i=0
      sizep = size(p)
      IF(size(p).GE.1)THEN
         IF(p(1).EQ.'NOOP') sizep = 0
      ENDIF
      allocate(Etmp(product(dime)))

      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,sizep)/)
      ! determine whether these integrals are zero.
      zero = any((/(field_list(idxp(i))%lin  .or. &
                    field_list(idxp(i))%quad .or. &
               .not.field_list(idxp(i))%bas,  i=1,sizep)/))
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dime) < sizep) call lsquit('rsp_twoave argument error: ' &
               // 'More perturbations than dimensions of property, size(dime) < sizep',-1)
      if (any(dime <= 0)) call lsquit('rsp_twoave argument error: ' &
               // 'Property has a zero or negative dimension, dime <= 0',-1)
      ! compute step lengths in E (cumulative products of dimensions)
      stepe(1) = 1
      do i = 2, size(dime)
         stepe(i) = stepe(i-1)*dime(i-1)
      end do
      ! reorder dimensions and step lengths in E according to permutation argument perm
      ddime = dime
      if (present(perm)) then
         if (size(perm) /= size(dime)) call lsquit('rsp_twoave argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dime)',-1)
         ! verify that perm is indeed a permutation
         do i = 1, size(dime)-1
            if (perm(i) <= 0 .or. perm(i) > size(dime) .or. &
                any(perm(i) == perm(i+1:size(dime)))) call lsquit('rsp_twoave ' &
                      // 'argument error: Permutation must contain each number exactly once',-1)
         end do
         ddime = (/( dime(perm(i)), i=1,size(dime))/)
         stepe = (/(stepe(perm(i)), i=1,size(dime))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= sizep) call lsquit('rsp_twoave argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)',-1)
         if (any(comp <= 0)) call lsquit('rsp_twoave argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0',-1)
         ccomp = comp
      end if
      if(sizeP.GT.0)then
         if (any(ccomp + ddime(:sizep) - 1 > pert_shape(mol,p))) &
              call lsquit('rsp_twoave argument error: Lowest component index plus ' &
              // 'dimension exceeds dimension of perturbation, comp + dime > pert_shape(mol,p)',-1)
      endif
      ! sort perturbations p so that idxp is descending
      pp = p
      do i = 1, sizep
         !find perturbation after i, with highest idxp (secondly highest ddime)
         j = i
         do k = j+1, sizep
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddime(k) > ddime(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddime(i);  ddime(i) = ddime(j);  ddime(j) = tmpi
         tmpi = stepe(i);  stepe(i) = stepe(j);  stepe(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
      end do
      ! verify that we have the correct number of perturbed densities
      nd = product(1+ddime(sizep+1:size(dime)))
      if (size(D) /= nd) call lsquit('rsp_twoave error: Number of' &
               // 'perturbed densities D does not correspond to dime (and perm)',-1)
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call lsquit('rsp_twoave error: Undefined matrix in argument D(:)',-1)
         deallocate(Etmp)
         return
      end if
      ! everything set up, so call core procedure oneave.
      ! Argument nd=0 is used when averaging over unperturbed density,
      ! in which case also perturbed nuclear attraction should be included
      call twoave(mol, sizep, size(dime)-sizep, pp, ccomp, ddime, D, Etmp)
      ! add oneavg property contribution in temporary array Etmp(:) to
      ! resulting property array E(:), while permuting indices according to pperm
      idxe = 0
      do j = 1, product(dime)
         i = 1 + sum(idxe * stepe)
         E(i) = E(i) + Etmp(j)
         do k = 1, size(dime)
            idxe(k) = idxe(k) + 1
            if (idxe(k) /= ddime(k)) exit
            idxe(k) = 0
         end do
      end do
      deallocate(Etmp)
   end subroutine



   subroutine twoave(mol, np, nd, p, c, de, D, E)
     implicit none
      !> structure containing integral program settings
      type(rsp_molcfg),  intent(in)  :: mol
      !> number of perturbations and order of density
      integer,           intent(in)  :: np, nd
      !> perturbation lables
      character(*),      intent(in)  :: p(np)
      !> lowest component of each perturbation
      integer,           intent(in)  :: c(np)
      !> dimensions of property (E)
      integer,           intent(in)  :: de(np+nd)
      !> un-/perturbed density matrices (expension),
      !> size(D) = product(1+de(np+1:np+nd))
      type(matrix),      intent(in)  :: D(*)
      !> where to ADD property contributions
      !> (works incrementally), size(E) = product(de)
      complex(realk),    intent(inout) :: E(*)
      !-------------------------------------
      real(realk), allocatable :: RR(:) !scratch
      real(realk)  :: R(6) !scratch
      integer      :: i, j, k, l, m, n, ii, jj, kk, ll, pd, pd1, na
      integer      :: lupri,luerr,nbast
      type(matrix) :: A(9) !scratch matrices
      logical      :: ok
      i=0
!!$      Interface 
!!$         SUBROUTINE II_get_xc_geoderiv_molgrad(LUPRI,LUERR,SETTING,nbast,D,grad,natoms)
!!$           use precision
!!$           use TYPEDEFTYPE, only: LSSETTING
!!$           use Matrix_module, only: matrix
!!$           IMPLICIT NONE
!!$           INTEGER               :: LUPRI,LUERR,nbast,natoms
!!$           TYPE(LSSETTING)       :: SETTING
!!$           TYPE(MATRIX)          :: D
!!$           REAL(REALK)           :: GRAD(3,natoms)
!!$         END SUBROUTINE II_GET_XC_GEODERIV_MOLGRAD
!!$      end Interface
!!$      interface 
!!$         subroutine II_get_xc_geoderiv_FxDgrad(LUPRI,LUERR,SETTING,nbast,D,b,grad,natoms)
!!$           use precision
!!$           use TYPEDEFTYPE, only: LSSETTING
!!$           use Matrix_module, only: matrix
!!$           IMPLICIT NONE
!!$           INTEGER               :: LUPRI,LUERR,nbast,natoms
!!$           TYPE(LSSETTING)       :: SETTING
!!$           TYPE(MATRIX)          :: D,b
!!$           REAL(REALK)           :: GRAD(3,natoms)
!!$         end subroutine II_get_xc_geoderiv_FxDgrad
!!$      end interface
!!$      interface
!!$         SUBROUTINE II_get_xc_geoderiv_GxDgrad(LUPRI,LUERR,SETTING,nbast,D,a,b,grad,natoms)
!!$           use precision
!!$           use TYPEDEFTYPE, only: LSSETTING
!!$           use Matrix_module, only: matrix
!!$           IMPLICIT NONE
!!$           INTEGER               :: LUPRI,LUERR,nbast,natoms
!!$           TYPE(LSSETTING)       :: SETTING
!!$           TYPE(MATRIX)          :: D,b,a
!!$           REAL(REALK)           :: GRAD(3,natoms)
!!$         END SUBROUTINE II_GET_XC_GEODERIV_GXDGRAD
!!$      end interface

      lupri = mol%lupri
      luerr = mol%luerr
      nbast = D(1)%nrow
      na = mol%natoms
      A(1) = 0 * D(1)
      A(2:) = (/(A(1), i=2,size(A))/) !scratch matrices
      pd  = product(de(np+1:np+nd))   !product of density dimensions
      pd1 = product(1+de(np+1:np+nd)) !size of D(*)
      if (np==0) then
         if (nd==0) call lsquit('rsp_twoave: unperturbed energy/integrals requested', &
                                mol%lupri)
         !The highest order Ds multiply the unperturbed F, which is unknown here,
         !so insist that no highest-order Ds are present
         if (.not.all((/(iszero(D(i)), i=pd1-pd+1,pd1)/))) &
            call lsquit('rsp_twoave: unperturbed Fock matrix requested', mol%lupri)
         !ajt fixme Stopping at nd=2, for now
         if (nd==2) then
            ! contract first density to Fock matrix, then trace with last
            do i = 0, de(1)-1
               ! Coulomb and exchange
               call twofck('  ', D(1), D(2+i), A(1:1))
               ! Kohn-Sham exchange-correlation
               call twofck_ks(mol,1, (/D(1),D(2+i)/), A(1))
               ! trace with first density matrix
               do j = 0, de(2)-1
                  E(1+i+de(1)*j) = trace(A(1),D(2+de(1)+j))
               end do
            end do
         else if (nd==3) then
            ! contract two of first-order densities to Fock, then trace with third
            do j = 0, de(2)-1
               do i = 0, de(1)-1
                  call twofck_ks(mol,2, (/D(1),D(2+i),D(2+de(1)+j)/), A(1))
                  do k = 0, de(3)-1
                     E(1+i+de(1)*(j+de(2)*k)) = trace(A(1),D(2+de(1)+de(2)+k))
                  end do
                  A(1) = 0 * A(1)
               end do
            end do
            ! contract first-order densities to Fock, then trace with second-order
            !ajt Fixme
            if (.not.all((/(iszero(D(i)), i=2+de(1)+de(2)+de(3),pd1)/))) &
               call lsquit('rsp_twoave: nd = 3 not fully implemented', mol%lupri)
         else if (nd==4) then
            !ajt FIXME dft-cubic contributions should go here
            ! verify that the contribution is zero (HF), otherwise quit
            j = 0
            n = 1
            m = de(1)*de(2)*de(3)*de(4)
            k = (de(1)+1)*(de(2)+1)*(de(3)+1)*(de(4)+1) - m
            ok = (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                  all((/(iszero(D(i)), i=k+1,k+m)/)))
            j=j+n; n=de(1); m=de(2)*de(3)*de(4); k=k-m
            ok = (ok .and. (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                            all((/(iszero(D(i)), i=k+1,k+m)/))))
            j=j+n; n=de(2); m=de(1)*de(3)*de(4); k=k-m
            ok = (ok .and. (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                            all((/(iszero(D(i)), i=k+1,k+m)/))))
            j=j+n; n=de(3); m=de(1)*de(2)*de(4); k=k-m
            ok = (ok .and. (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                            all((/(iszero(D(i)), i=k+1,k+m)/))))
            j=j+n; n=de(4); m=de(1)*de(2)*de(3); k=k-m
            ok = (ok .and. (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                            all((/(iszero(D(i)), i=k+1,k+m)/))))
            j=j+n; n=de(1)*de(2); m=de(3)*de(4); k=k-m
            ok = (ok .and. (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                            all((/(iszero(D(i)), i=k+1,k+m)/))))
            j=j+n; n=de(1)*de(3); m=de(2)*de(4); k=k-m
            ok = (ok .and. (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                            all((/(iszero(D(i)), i=k+1,k+m)/))))
            j=j+n; n=de(2)*de(3); m=de(1)*de(4); k=k-m
            ok = (ok .and. (all((/(iszero(D(i)), i=j+1,j+n)/)) .or. &
                            all((/(iszero(D(i)), i=k+1,k+m)/))))
            if (ok .and. .not.mol%setting%do_dft) then
               E(:de(1)*de(2)*de(3)*de(4)) = 0
            else if (mol%setting%do_dft) then
               call lsquit('rsp_twoave: XC-contribution to cubic (nd=4) not implemented', mol%lupri)
            else
               call lsquit('rsp_twoave: nd = 4 not fully implemented', mol%lupri)
            end if
         else if (nd > 4) then
            call lsquit('rsp_twoave: nd > 3 not implemented', mol%lupri)
         end if
      !London magnetic, no nd==0 because because MAG anti
      else if (np==1 .and. p(1)=='MAG' .and. nd /= 0) then
         ! II_get_MagDeriv_(F,G)xD_DFT takes initialized
         ! II also takes initiated matrices
         do i = 4, 6
            call mat_ensure_alloc(A(i))
         end do
         if (mol%setting%do_dft) then !scratch in A(4:6)
            do i = 7, 9
               call mat_ensure_alloc(A(i))
            end do
         end if
         if (all((/(iszero(D(i)), i=pd1-pd+1,pd1)/))) then
            E(:product(de)) = 0
         else
            ! contract unperturbed density to Fock, then trace with highest-order
            ! Coulomb-exchange - sym Dmat
            call II_get_magderivF(mol%lupri,mol%luerr,mol%SETTING,nbast,D(1:1),A(4:6))
            ! Kohn-Sham contribution
            if (mol%setting%do_dft) then
               call II_get_xc_magderiv_kohnsham_mat(mol%lupri,mol%luerr,mol%SETTING,nbast,D(1),A(7:9))
               do i = 1, 3
                  A(i) = A(3+i) + A(6+i)
               end do
            end if
            do j = 0, pd-1
               do i = 0, de(1)-1
                  E(1+i+de(1)*j) = trace(A(c(1)+i), D(pd1-pd+1+j))
               end do
            end do
         end if
         if (nd==1) then
            ! nothing more
         else if (nd==2) then
            ! contract first density to Fock matrix, then trace with second
            do j = 0, de(2)-1
               if (iszero(D(2+j))) cycle
               ! Coulomb-exchange with a possible non sym Dmat
               call II_get_magderivF(mol%lupri,mol%luerr,mol%SETTING,nbast,D(2+j:2+j),A(4:6))
               ! Kohn-Sham
               if (mol%setting%do_dft) then
                  if (D(2+j)%complex) then
                     call II_get_xc_magderiv_linrsp(mol%lupri, mol%luerr, mol%SETTING, nbast, &
                                       (/ mat_get_part(D(2+j), imag=.false.) /), D(1), A(7:9), 1)
                  else
                     call II_get_xc_magderiv_linrsp(mol%lupri, mol%luerr, mol%SETTING, nbast, &
                                                    D(2+j:2+j), D(1), A(7:9), 1)
                  end if
                  do i = 1, 3
                     A(i) = A(i+3) + (-1.d0)*A(6+i) !negative sign on KS contrib
                  end do
                  if (D(2+j)%complex) then
                     call II_get_xc_magderiv_linrsp(mol%lupri,mol%luerr,mol%SETTING,nbast, &
                                         (/ mat_get_part(D(2+j), imag=.true.) /), D(1), A(7:9),1)
                     do i = 1, 3
                        A(i) = A(3+i) - (0E0_realk,1E0_realk) * A(6+i) !negative sign on KS contrib
                     end do
                  end if
               else
                  do i = 1, 3
                     A(i) = A(i+3)
                  end do
               end if
               ! trace with second density matrix
               do k = 0, de(3)-1
                  do i = 0, de(1)-1
                     E(1+i+de(1)*(j+de(2)*k)) = E(1+i+de(1)*(j+de(2)*k)) &
                                              + trace(A(c(1)+i), D(2+de(2)+k))
                  end do
               end do
            end do
         else
            call lsquit('rsp_twoave: MAG, nd > 2 not implemented',-1)
         end if
         do i = 4, 6
            A(i)=0  
         end do
         if (mol%setting%do_dft) then 
            do i = 7, 9
               A(i)=0
            end do
         end if
      else if (np==1 .and. p(1)=='GEO ') then
         allocate(RR(3*na*2))
         ! highest-order contribution
         do j = 0, pd-1
            if (iszero(D(pd1-pd+1+j))) then
               E( 1+de(1)*j : de(1)*(j+1) ) = 0
               cycle
            end if
            ! Coulomb-exchange
            if (nd==0) call prop_intifc_2el(mol, p, D(1:1), RR(:3*na), A(1:1))
            if (nd/=0) call prop_intifc_2el(mol, p, (/D(1), D(pd1-pd+1+j)/), &
                                            RR(:3*na), A(1:1))
            ! for molgra (nd==0), factor 1/2 on these integrals
            if (nd==0) RR(:3*na) = RR(:3*na)/2
            if (mol%setting%do_dft) then
               if (nd==0) then
                  call II_get_xc_geoderiv_molgrad(LUPRI, LUERR, &
                                      mol%SETTING, nbast, D(1), &
                                      RR(3*na+1:3*na*2), na)
               else 
                  call II_get_xc_geoderiv_FxDgrad(LUPRI, LUERR, &
                                  mol%SETTING, nbast, D(1), D(pd1-pd+1+j), &
                                  RR(3*na+1:3*na*2), na)
               end if
               RR(:3*na) = RR(:3*na) + RR(3*na+1:3*na*2)
            end if
            E( 1+de(1)*j : de(1)*(j+1) ) = RR(c(1):c(1)+de(1)-1)
         end do
         if (nd==0 .or. nd==1) then
            ! nothing more
         else if (nd==2) then
            ! integrals over products of first order densities
            do k = 0, de(3)-1
               do j = 0, de(2)-1
                  ! Coulomb-exchange
                  call prop_intifc_2el(mol, p, (/D(2+j), D(2+de(2)+k)/), &
                                       RR(:3*na), A(1:1))
                  ! Kohn-Sham exchange-correlation
                  if (mol%setting%do_dft) then
                     call II_get_xc_geoderiv_GxDgrad(LUPRI, LUERR, &
                                         mol%SETTING, nbast, D(1), &
                                         D(2+j), D(2+de(2)+k), &
                                         RR(3*na+1:3*na*2), na)
                     RR(:3*na) = RR(:3*na) + RR(3*na+1:3*na*2)
                  end if
                  E( 1+de(1)*(j+de(2)*k) : de(1)*(1+j+de(2)*k) ) &
                            = E( 1+de(1)*(j+de(2)*k) : de(1)*(1+j+de(2)*k) ) &
                            + RR(c(1):c(1)+de(1)-1)
               end do
            end do
         else
            call lsquit('rsp_twoave: GEO, nd > 2 not implemented', mol%lupri)
         end if
         deallocate(RR)
         if (mol%setting%do_dft) print* !after all the "...integrated to nn electrons..." prints
      else
         print *,'rsp_twoave: no integrals for these perturbations: ', &
                 (p(i)//' ', i=1,np)
         call lsquit('rsp_twoave: no such integrals', mol%lupri)
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine



   !> Calculates the 1-electron integrals perturbed by the perturbations p(:)
   !> and ADDS the integrals to the perturbed Fock matrices F(:).
   !> Front for the private subroutine 'oneint' below, checking the
   !> arguments' dimensions, and doing permutations
   !> S0 is passed as argument only as reference to nuclei and basis set
   subroutine rsp_oneint(mol, S0, p, dimp, F, S, comp, freq)
     implicit none
      !> mol/basis data needed by integral program
      type(rsp_molcfg),  intent(inout) :: mol
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> shape(F), size(dimp) = size(p), dimensions of perturbed
      !> Fock matrices F(:)
      integer,           intent(in) :: dimp(:)
      !-----------------------------------------------------------
      !> perturbed Fock matrices to fill with perturbed integrals,
      !> size(F) = product(dimp)
      type(matrix),  optional, intent(inout) :: F(*)
      !> optionally return the corresponding perturbed overlap matrices
      type(matrix),  optional, intent(inout) :: S(*)
      !-----------------------------------------------------------
      !> comp(np), starting component index for each p. Default 1 1 ... 1
      integer,       optional, intent(in)    :: comp(:)
      !> freq(np), complex frequencies
      !> for each p, default all zero. These multiply the
      !> half-derivative overlap integrals, thus no contribution
      !> if basis independent of p
      complex(realk), optional, intent(in)   :: freq(:)
      !-----------------------------------------------------------
      integer        :: idxp(size(p)), pperm(size(dimp)), ccomp(size(p)), &
                        stepf(size(dimp)), ddimp(size(dimp)), idxf(size(dimp)), &
                        i, j, k, tmpi
      character(4)   :: pp(size(p)), tmpp
      complex(realk) :: ffreq(size(p)), tmpf
      logical        :: zero, bas
      type(matrix),pointer :: Ftmp(:), Stmp(:)
      i=0
      allocate(Ftmp(product(dimp)))
      allocate(Stmp(product(dimp)))

      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero. If there is more than 1
      ! linear perturbation (like 'EL') or more than two quadratic (like 'MAGO')
      zero = .false.
      j = 0
      k = 0
      do i = 1, size(p)
         if (.not.field_list(idxp(i))%lin .and. &
             .not.field_list(idxp(i))%quad) cycle
         ! if both linear and quadratic, interpret as ZERO (like EXCI)
         if (field_list(idxp(i))%lin .and. &
             field_list(idxp(i))%quad) zero = .true.
         if (j /= 0) then !if second lin or quad perturbation
            ! if j was EL, then i is second EL (or first MAGO)
            if (field_list(idxp(j))%lin) zero = .true.
            ! if j was MAGO, and different from i, also zero
            if (field_list(idxp(j))%quad .and. &
               idxp(i) /= idxp(j)) zero = .true. !two different quadratic
            !third linear or quadratic perturbation
            if (k /= 0) zero = .true.
            k = j
         end if
         j = i
      end do
      ! check dimensions argument dimp, verify that dimensions are positive
      if (size(dimp) /= size(p)) call lsquit('rsp_oneint argument error: ' &
               // 'Different number of perturbations and dimensions, size(dimp) /= size(p)', mol%lupri)
      if (any(dimp <= 0)) call lsquit('rsp_oneint argument error: ' &
               // 'Perturbations have a zero or negative dimension, dimp <= 0', mol%lupri)
      ! compute step lengths in F (cumulative product of dimp)
      stepf(1) = 1
      do i = 2, size(dimp)
         stepf(i) = stepf(i-1)*dimp(i-1)
      end do
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call lsquit('rsp_oneint argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)', mol%lupri)
         if (any(comp <= 0)) call lsquit('rsp_oneint argument error: ' &
                    // 'Lowest component indices must be positive, but comp <= 0', mol%lupri)
         ccomp = comp
      end if
      if (any(ccomp + dimp - 1 > pert_shape(mol,p))) &
         call lsquit('rsp_oneint argument error: Lowest component index plus ' &
                // 'dimension exceeds dimension of perturbation, comp + dimp - 1 > pert_shape(mol,p)', mol%lupri)
      ! check optional argument freq, default to zero
      ffreq = 0
      if (present(freq)) then
         if (size(freq) /= size(p)) call lsquit('rsp_oneave ' &
               // 'argument error: Wrong number of frequencies, size(freq) /= size(p)', mol%lupri)
         ffreq = freq
      end if
      ! sort perturbations p so that idxp is descending
      pp = p
      ddimp = dimp
      do i = 1, size(p)
         !find perturbation after i, with highest idxp (secondly highest ddimp)
         j = i
         do k = j+1, size(p)
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddimp(k) > ddimp(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddimp(i);  ddimp(i) = ddimp(j);  ddimp(j) = tmpi
         tmpi = stepf(i);  stepf(i) = stepf(j);  stepf(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
         tmpf = ffreq(i);  ffreq(i) = ffreq(j);  ffreq(j) = tmpf
      end do
      ! zero (deallocate) F and S before calling oneint
      do i = 1, product(dimp)
         if (present(F)) F(i) = 0*S0
         if (present(S)) S(i) = 0*S0
      end do
      ! early return if the result is zero
      if (zero) return
      ! everything set up, so call core procedure oneint
      call oneint(mol, S0, size(p), pp, ccomp, ddimp, ffreq, Ftmp, Stmp)
      ! add oneavg property contribution in temporary array Etmp(:) to
      ! resulting property array E(:), while permuting indices according to pperm
      bas = all(pert_basdep(p))
      idxf = 0
      do j = 1, product(dimp)
         i = 1 + sum(idxf * stepf)
         if (present(F)) then
            call mat_move(Ftmp(j), F(i))
         else
            Ftmp(j) = 0
         end if
         if (present(S) .and. bas) then
            call mat_move(Stmp(j), S(i))
         else
            Stmp(j) = 0
         end if
         do k = 1, size(dimp)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimp(k)) exit
            idxf(k) = 0
         end do
      end do
      deallocate(Ftmp)
      deallocate(Stmp)
   end subroutine

   subroutine oneint(mol, S0, np, p, c, dp, w, F, S)
      !> mol/basis data needed by integral program
      type(rsp_molcfg),  intent(inout) :: mol
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> number of perturbations
      integer,           intent(in) :: np
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> lowest component of each perturbation
      integer,           intent(in) :: c(np)
      !> dimensions of property (F and S)
      integer,           intent(in) :: dp(np)
      !> frequency of each p
      complex(realk),    intent(in) :: w(np)
      !--------------------------------------------
      !> perturbed Fock matrices
      type(matrix),   intent(inout) :: F(*)
      !> perturbed overlap matrices
      type(matrix),   intent(inout) :: S(*)
      !--------------------------------------------
      integer      :: i, j, k, ii, jj, kk
      type(matrix) :: A(6) !scratch matrices
      real(realk)  :: R(1) !dummy
      if (np==0) then
         call lsquit('rsp_oneint error:' &
                // ' Unperturbed one-electron integrals requested', mol%lupri)
      else if (np==1 .and. p(1)=='EL  ') then
         do i=1,3
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         call II_get_prop(mol%lupri, mol%luerr, mol%setting, A(:3), 3, 'DIPLEN ')
         do i=0, dp(1)-1
            F(1+i) = A(c(1)+i)
         end do
         A(:3) = 0
      else if (np==1 .and. p(1)=='MAGO') then
         do i=1,3
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         call II_get_prop(mol%lupri, mol%luerr, mol%setting, A(:3), 3, 'ANGMOM ')
         do i=0, dp(1)-1
            F(1+i) = (1/2d0) * A(c(1)+i) 
         end do
         A(:3) = 0
      else if (np==1 .and. p(1)=='MAG ') then
         do i=1,3
            A(i) = 0*S0
            call mat_ensure_alloc(A(i))
         end do
         call II_get_prop(mol%lupri, mol%luerr, mol%setting, A(:3), 3, 'MAGMOM ')
         do i=0, dp(1)-1
            F(1+i) = A(c(1)+i) !ajt-to-Thomas: Verify sign/scale by comparison with other results
         end do
         if (w(1)==0) then !no -i/2 Tb contribution
            call II_get_magderivOverlap(A(1:3), mol%setting, mol%lupri, mol%luerr)
            do i=0, dp(1)-1
               S(1+i) = A(c(1)+i) !ajt-to-Thomas: Verify sign/scale by comparison with other results
            end do
         else 
            call II_get_magderivOverlapR(A(1:3), mol%setting, mol%lupri, mol%luerr)
            do i=0, dp(1)-1
               S(1+i) = A(c(1)+i) 
            end do
            do i = 0, dp(1)-1
               F(i+1) = F(i+1) - w(1)/2 * S(i+1)
               F(i+1) = F(i+1) - w(1)/2 * trans(S(i+1))
               S(i+1) = S(i+1) - trans(S(i+1))
            enddo
         end if
         A(:3) = 0
      else
         print *,'rsp_oneint: No integrals for these perturbations:', &
                    (' ' // p(i),i=1,np)
         call lsquit('rsp_oneint: No such integrals', mol%lupri)
      end if
    end subroutine


   !> Contracts the 2-electron and Kohn-Sham integrals perturbed by the
   !> perturbations p(:) with the perturbed density matrix expansion in D(:)
   !> (e.g. D=(/D,Dx,Dy,Dxy/) for a 2nd order expansion), and ADD the resulting
   !> Fock matrix contibution to the array F(:).
   !> Front for the private subroutine 'twoave' below, checking the arguments'
   !> dimensions, and doing permutations
   !> D(1) serves as reference to nuclei, basis and model/functional
   subroutine rsp_twoint(mol, p, D, dimf, F, perm, comp)
     implicit none
      !> mol/basis data needed by integral program
      type(rsp_molcfg),  intent(in) :: mol
      !> p(np) perturbation lables
      character(*),      intent(in) :: p(:)
      !> (un)perturbed density matrices to contract perturbed
      !> one-electron integrals with.
      !> If perm present, size(D) = product(1+dime(perm(np+1:np+nd))),
      !> if not present, size(D) = product(1+dime(np+1:np+nd))
      type(matrix),      intent(in) :: D(:)
      !> dime(np+nd) = shape(F), dimensions of perturbed Fock matrices F(:)
      integer,           intent(in) :: dimf(:)
      !---------------------------------------------------------------
      !> Perturbed Fock matrices, works incrementally,
      !> thus contributions are ADDED to F(*). size(F) = product(dimf)
      type(matrix), intent(inout) :: F(*)
      !---------------------------------------------------------------
      !> perm(np+nd), permutation of indices.
      !> For each dimension of p and D, the corresponding dimension of F.
      !> Default 1 2 ... np+nd (no permutation)
      integer,      optional, intent(in) :: perm(:)
      !> comp(np), starting component index for each p. Default 1 1 ... 1
      integer,      optional, intent(in) :: comp(:)
      !---------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dimf)), ccomp(size(p)), &
                      & stepf(size(dimf)), ddimf(size(dimf)), idxf(size(dimf)), &
                      & i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      logical      :: zero
      type(matrix),pointer :: Ftmp(:)
      integer      :: sizep
      i=0
      allocate(Ftmp(product(dimf)))
      ! verify that all perturbation labels exist, find index of each
      sizep = size(p)
      IF(size(p).GE.1)THEN
         IF(p(1).EQ.'NOOP') sizep = 0
      ENDIF
      idxp = (/(idx(p(i)), i=1,sizep)/)
      ! determine whether these integrals are zero.
      zero = any((/(field_list(idxp(i))%lin  .or. &
                    field_list(idxp(i))%quad .or. &
               .not.field_list(idxp(i))%bas, i=1,sizep)/))
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dimf) < sizep) call lsquit('rsp_twoint argument error: ' &
               // 'More perturbations than dimensions of F(:), size(dimf) < size(p)',-1)
      if (any(dimf <= 0)) call lsquit('rsp_twoint argument error: ' &
               // 'Perturbed Fock F(:) has a zero or negative dimension, dimf <= 0',-1)
      ! compute step lengths in E (cumulative products of dimensions)
      stepf(1) = 1
      do i = 2, size(dimf)
         stepf(i) = stepf(i-1)*dimf(i-1)
      end do
      ! reorder dimensions and step lengths in F according to permutation argument perm
      ddimf = dimf
      if (present(perm)) then
         if (size(perm) /= size(dimf)) call lsquit('rsp_twoint argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dimf)',-1)
         ! verify that perm is indeed a permutation
         do i = 1, size(dimf)-1
            if (perm(i) <= 0 .or. perm(i) > size(dimf) .or. &
                any(perm(i) == perm(i+1:size(dimf)))) call lsquit('rsp_twoint ' &
                      // 'argument error: Permutation must contain each number exactly once',-1)
         end do
         ddimf = (/( dimf(perm(i)), i=1,size(dimf))/)
         stepf = (/(stepf(perm(i)), i=1,size(dimf))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= sizep) call lsquit('rsp_twoint argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)',-1)
         if (any(comp <= 0)) call lsquit('rsp_twoint argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0',-1)
         ccomp = comp
      end if
      if(sizeP.GT.0)then
         if (any(ccomp + ddimf(:sizep) - 1 > pert_shape(mol,p))) &
              call lsquit('rsp_twoint argument error: Lowest component index plus ' &
              // 'dimension exceeds dimension of perturbation, comp + dimf > pert_shape(mol,p)',-1)
      endif
      ! sort perturbations p so that idxp is descending
      pp = p
      do i = 1, sizep
         !find perturbation after i, with highest idxp (secondly highest ddime)
         j = i
         do k = j+1, sizep
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddimf(k) > ddimf(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddimf(i);  ddimf(i) = ddimf(j);  ddimf(j) = tmpi
         tmpi = stepf(i);  stepf(i) = stepf(j);  stepf(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
      end do
      ! verify that we have the correct number of perturbed densities
      nd = product(1+ddimf(sizep+1:size(dimf)))
      if (size(D) /= nd) call lsquit('rsp_twoint error: Number of' &
               // 'perturbed densities D does not correspond to dimf (and perm)',-1)
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call lsquit('rsp_twoint error: Undefined matrix in argument D(:)',-1)
         return
      end if
      ! permute perturbed Fock matrices in F(:) over to Ftmp(:)
      idxf = 0
      do j = 1, product(dimf)
         i = 1 + sum(idxf * stepf)
         call mat_move(F(i), Ftmp(j))
         do k = 1, size(dimf)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimf(k)) exit
            idxf(k) = 0
         end do
      end do
      ! everything set up, so call core procedure twoint
      call twoint(mol, sizep, size(dimf)-sizep, pp, ccomp, ddimf, D, Ftmp)
      ! 'un-permute' perturbed Fock matrices from Ftmp(:) back into F(:)
      idxf = 0
      do j = 1, product(dimf)
         i = 1 + sum(idxf * stepf)
         call mat_move(Ftmp(j), F(i))
         do k = 1, size(dimf)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimf(k)) exit
            idxf(k) = 0
         end do
      end do
      Ftmp = 0
      deallocate(Ftmp)
   end subroutine

   subroutine twoint(mol, np, nd, p, c, df, D, F)
      !> mol/basis data needed by integral program
      type(rsp_molcfg),  intent(in) :: mol
      !> number of perturbations and order of density
      integer,           intent(in) :: np, nd
      !> perturbation lables
      character(*),      intent(in) :: p(np)
      !> lowest component of each perturbation
      integer,           intent(in) :: c(np)
      !> dimensions of perturbed Fock matrix F
      integer,           intent(in) :: df(np+nd)
      !> un-/perturbed density matrices (expension),
      !> size(D) = product(1+df(np+1:np+nd))
      type(matrix),      intent(in) :: D(*)
      !--------------------------------------------------
      !> where to ADD property contributions
      !> (works incrementally), size(F) = product(df)
      type(matrix),   intent(inout) :: F(*)
      !--------------------------------------------------
      integer      :: i, j, k, l, ii, jj, kk, ll, pd, pd1,lupri,luerr,nbmat
      type(matrix) :: A(6) !scratch matrices
!!$      interface
!!$         SUBROUTINE II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,SETTING,nbast,D,F)
!!$           use precision
!!$           use TYPEDEFTYPE, only: LSSETTING
!!$           use Matrix_module, only: matrix
!!$           IMPLICIT NONE
!!$           INTEGER               :: LUPRI,LUERR,nbast
!!$           TYPE(LSSETTING)       :: SETTING
!!$           TYPE(MATRIX)          :: D,F(3)
!!$         end SUBROUTINE II_get_xc_magderiv_kohnsham_mat
!!$      end interface
      lupri = mol%lupri
      luerr = mol%luerr
      A(1) = 0*D(1)
      ! ajt: ifort 11.1 (at least) seems to miscompile this line, somehow, resulting in a leaving a 
      !      seemingly valid and initialized type(matrix) on the stack, where it will coincide with 
      !      local variable A(2) in a subsequent call to subroutine oneint, causing an error when
      !      attempting to change that matrix (deallocate previous content)
      ! A(2:) = (/(A(1), i=2,size(A))/) !scratch matrices
      do i = 2, size(A)
         A(i) = A(1)
      end do
      pd  = product(df(np+1:np+nd))   !product of density dimensions
      pd1 = product(1+df(np+1:np+nd)) !size of D(*)
      if (np==0) then
         if (nd==0) call lsquit('rsp_twoint error: Unperturbed ' &
                           // 'two-electron Fock matrix requested',-1)
         ! di_GET_GbDs and di_get_sigma expects A initialized (allocated)
         call mat_ensure_alloc(A(1))
         do i = 0, pd-1
            if (iszero(D(pd1-pd+1+i))) cycle
            ! Coulomb-exchange
            call twofck('  ', D(1), D(pd1-pd+1+i), A(1:1))
            F(i+1) = F(i+1) + A(1)
            ! Kohn-Sham exchange-correlation
            call twofck_ks(mol,1, (/D(1),D(pd1-pd+1+i)/), F(1+i))
         end do
         if (nd==0 .or. nd==1 .or. .not.mol%setting%do_dft) then
            ! nothing more
         else if (nd==2) then
            do j = 0, df(2)-1
               do i = 0, df(1)-1
                  call twofck_ks(mol,2, (/D(1),D(2+i),D(2+df(1)+j)/), F(1+i+df(1)*j))
               end do
            end do
         else
            call lsquit('rsp_twoint: nd > 2 not implemented with DFT', mol%lupri)
         end if
         A(1)=0
      else if (np==1 .and. p(1)=='MAG ') then
         do i = 1, 3
            call mat_ensure_alloc(A(i))
         end do
         IF(nd==0)THEN
            ! Construction of the magnetic derivative Fock matrix
            ! Fb = Gb(D(1)) + Fxcb  (without differentiation on the D mat and no 1electron contribution)
            !
            ! Coulomb-exchange matrix constructed from first order magnetic
            ! derivative integrals, meaning the Coulomb and Exchange matrix
            ! contribution to the magnetic derivative Fock matrix
            ! (without differentiation on the D mat)
            call II_get_magderivF(LUPRI,LUERR,MOL%SETTING,D(1)%nrow,D(1:1),A(1:3))
            !Fb = Fb + Gb(D(1))
            do i=0, df(1)-1
               F(i+1) = F(i+1) + A(c(1)+i)
            enddo
            if (mol%setting%do_dft) then !scratch in A(4:6)
               ! Calculates the Kohn-Sham exchange-correlation contribution to the 
               ! magnetic derivative kohn-sham matrix (without differentiation on the D mat)
               call II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,MOL%SETTING,D(1)%nrow,D(1),A(1:3))
               !Fb = Fb + Fxcb
               do i=0, df(1)-1
                  F(i+1) = F(i+1) + A(c(1)+i)
               enddo
            endif
         ELSEIF(nd==1)THEN
            ! Construction of the magnetic derivative linear response matrix 
            ! Fb = Gb(D(2+j)) + Gxcb(D(2+j))  (without differentiation on the D mat and no 1electron contribution)
            if (mol%setting%do_dft) then !scratch in A(4:6)
               do i = 4, 6
                  call mat_ensure_alloc(A(i))
               end do
            end if
            do j=0, df(2)-1
               ! Coulomb-exchange matrix constructed from first order magnetic
               ! derivative integrals, meaning the Coulomb and Exchange matrix
               ! contribution to the magnetic derivative Fock matrix
               ! (without differentiation on the D mat)
               call II_get_magderivF(LUPRI,LUERR,MOL%SETTING,D(1)%nrow,D(2+j:2+j),A(1:3))
               if (mol%setting%do_dft) then 
                  ! Calculates the Kohn-Sham exchange-correlation contribution to the 
                  ! magnetic derivative kohn-sham matrix (without differentiation on the D mat)
                  nbmat = 1
                  call II_get_xc_magderiv_linrsp(LUPRI,LUERR,MOL%SETTING,D(1)%nrow,D(3+j:3+j),D(1),A(4:6),nbmat)
                  A(1) = A(1) + A(4)
                  A(2) = A(2) + A(5)
                  A(3) = A(3) + A(6)
               endif
               do i=0, df(1)-1
                  F(1+i+df(1)*j) = F(1+i+df(1)*j) + A(c(1)+i)
               enddo
             enddo
             if (mol%setting%do_dft) then !scratch in A(4:6)
                do i = 4, 6
                   A(i)=0
                end do
             endif
         ELSEIF(nd.GE.2)THEN
            call lsquit('Case not implemented in twoint',-1)
         ENDIF
         ! if nd==0: A(:3) = Gb(D(1)) + Fxcb
         !           F(:df(1)) += A(c(1):c(1)+df(1)-1)
         
         ! if nd==1: do j=0, df(2)-1
         !               A(:3) = Gb(D(2+j)) + Gxcb(D(2+j))
         !               do i=0, df(1)-1
         !                  F(1+i+df(1)*j) += A(c(1)+i)
         ! if nd>=2: lsquit - not implemented
         A(:3) = 0
      else
         print *,'rsp_twoint: No integrals for these perturbations: ', &
                 (p(i)//' ', i=1,np)
         call lsquit('rsp_twoint: No such integrals', mol%lupri)
      end if
      A(:) = 0 !delete scratch matrices
    end subroutine twoint

   subroutine twofck(what, Dmat, Bmat, F, a, b)
   ! Add (un)perturbed 2e-contributions to Fock matrices
   ! what='  ' -> unperturbed
   ! what='G ' -> 3 geometry-perturbed wrt. nucleus a
   ! what='M ' -> 3 magnetic-field-perturbed (London)
   ! what='GG' -> 9 2nd-order geometry-perturbed wrt. nuclei a, b
   ! what='MM' -> 6 2nd-order magnetic-field-perturbed (London)
   ! what='GM' -> 9 2nd-order mixed geometry and magnetic
      character*2,       intent(in)    :: what
      type(matrix),      intent(in)    :: Dmat
      type(matrix),      intent(in)    :: Bmat(1)
      type(matrix),      intent(inout) :: F(:)
      integer, optional, intent(in)    :: a, b
      type(matrix) reimF(1)
      integer      nf, i
      type(matrix) :: Gxc(1)
      integer      :: nbast, ndmat, luoutout
      type(matrix) :: D_complex(1)
      nf = size(F)
      luoutout = 6
      !make usre F is properly allocated
      do i = 1, nf
         if (iszero(F(i)) .or. (Bmat(1)%complex .and. .not.F(i)%complex)) then
            F(i) = 0*Bmat(1) !will deallocate
            call mat_ensure_alloc(F(i))
         end if
      end do
      !process complex Bmat with two calls
	  nbast = Bmat(1)%nrow
	  call mat_init(D_complex(1),nbast,nbast)
      if (Bmat(1)%complex) then
         ! real part
         reimF(1:1) = (/mat_get_part(F(1), imag=.false.)/)
         !call di_GET_GbDs(6, 6, mat_get_part(Bmat(1), imag=.false.), reimF(1))
         call mat_assign(D_complex(1), mat_get_part(Bmat(1), imag=.false.))
         ndmat = 1
         call di_GET_GbDs_and_XC_linrsp(reimF, Gxc, luoutout, luoutout, D_complex, ndmat, nbast, Dmat, .false.)

         ! imaginary part
         reimF(1:1) = (/mat_get_part(F(1), imag=.true.)/)
         !call di_GET_GbDs(6, 6, mat_get_part(Bmat(1), imag=.true.), reimF(1))
         call mat_assign(D_complex(1), mat_get_part(Bmat(1), imag=.true.))
         ndmat = 1
         call di_GET_GbDs_and_XC_linrsp(reimF, Gxc, luoutout, luoutout, D_complex, ndmat, nbast, Dmat, .false.)
       
         
         reimF(1:1) = 0
      else
         !call di_GET_GbDs(6, 6, D(1), F(1))
         ndmat = 1
         call di_GET_GbDs_and_XC_linrsp(F, Gxc, luoutout, luoutout, Bmat, ndmat, nbast, Dmat, .false.)
      end if
      call mat_free(D_complex(1))
   end subroutine



   subroutine twofck_ks(mol, n, D, F)
     implicit none
      type(rsp_molcfg),  intent(in)    :: mol
      integer,           intent(in)    :: n
      type(matrix),      intent(in)    :: D(n+1)
      type(matrix),      intent(inout) :: F
      type(matrix) :: A(1)
      integer      :: i, lupri, luerr, nbast
      i=0
!!$      interface
!!$         SUBROUTINE II_get_xc_linrsp(LUPRI,LUERR,SETTING,nbast,b,D,G,nbmat)
!!$           use precision
!!$           use TYPEDEFTYPE, only: LSSETTING
!!$           use Matrix_module, only: matrix
!!$           IMPLICIT NONE
!!$           INTEGER               :: LUPRI,LUERR,nbast,nbmat
!!$           TYPE(LSSETTING)       :: SETTING
!!$           TYPE(MATRIX)          :: b(nbmat),D,G(nbmat)
!!$         end SUBROUTINE II_get_xc_linrsp
!!$      end interface
!!$      interface
!!$         SUBROUTINE II_get_xc_quadrsp(LUPRI,LUERR,SETTING,nbast,b,c,D,T)
!!$           use precision
!!$           use TYPEDEFTYPE, only: LSSETTING
!!$           use Matrix_module, only: matrix
!!$           IMPLICIT NONE
!!$           INTEGER               :: LUPRI,LUERR,nbast
!!$           TYPE(LSSETTING)       :: SETTING
!!$           TYPE(MATRIX)          :: b,c,D,T
!!$         end SUBROUTINE II_get_xc_quadrsp
!!$      end interface
      ! no density functional, then return
      if (.not.mol%setting%do_dft) return
      ! if any of the perturbed density matrices are zero, no contribution
      if (any((/(iszero(D(i)), i=2,n+1)/))) return
      lupri = mol%lupri
      luerr = mol%luerr
      nbast = mol%zeromat%nrow
      A(1) = 0*D(1)
      call mat_ensure_alloc(A(1))
      if (n==0) then
         call lsquit('rsp_contribs/twofck_ks error: Kohn-Sham contribution ' &
                // 'to unperturbed Fock matrix requested, but not implemented', mol%lupri)
      else if (n==1) then
         if (D(2)%complex) then
            call II_get_xc_linrsp(lupri, luerr, mol%setting, nbast, &
                                  (/mat_get_part(D(2), imag=.false.)/), &
                                  D(1), A, 1)
         else
            call II_get_xc_linrsp(lupri, luerr, mol%setting, nbast, &
                                  D(2:2), D(1), A, 1)
         end if
         F = F + A(1)
         if (D(2)%complex) then
            call II_get_xc_linrsp(lupri, luerr, mol%setting, nbast, &
                                  (/mat_get_part(D(2), imag=.true.)/), &
                                  D(1), A, 1)
            F = F + (0E0_realk,1E0_realk)*A(1)
         end if
      else if (n==2) then
         if (D(2)%complex .and. D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, &
                                   mat_get_part(D(2), imag=.false.), &
                                   mat_get_part(D(3), imag=.false.), D(1),A(1))
         if (D(2)%complex .and. .not.D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, &
                                   mat_get_part(D(2), imag=.false.), D(3), D(1), A(1))
         if (.not.D(2)%complex .and. D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, D(2), &
                                   mat_get_part(D(3), imag=.false.), D(1), A(1))
         if (.not.D(2)%complex .and. .not.D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, &
                                   D(2), D(3), D(1), A(1))
         F = F + A(1)
         if (D(2)%complex .and. D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, &
                                   mat_get_part(D(2), imag=.false.), &
                                   mat_get_part(D(3), imag=.true.), D(1), A(1))
         if (.not.D(2)%complex .and. D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, D(2), &
                                   mat_get_part(D(3), imag=.true.), D(1), A(1))
         if (D(3)%complex) F = F + (0E0_realk,1E0_realk)*A(1)
         if (D(2)%complex .and. D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, &
                                   mat_get_part(D(2), imag=.true.), &
                                   mat_get_part(D(3), imag=.false.), D(1), A(1))
         if (D(2)%complex .and. .not.D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, &
                                   mat_get_part(D(2), imag=.true.), D(3), D(1),A(1))
         if (D(2)%complex) F = F + (0E0_realk,1E0_realk)*A(1)
         if (D(2)%complex .and. D(3)%complex) &
            call II_get_xc_quadrsp(lupri, luerr, mol%setting, nbast, &
                                   mat_get_part(D(2), imag=.true.), &
                                   mat_get_part(D(3), imag=.true.), D(1), A(1))
         if (D(2)%complex .and. D(3)%complex) F = F - A(1)
      else
         call lsquit('rsp_contribs/twofck_ks error: n > 2 not implemented', mol%lupri)
      end if
      A(1) = 0
   end subroutine



   function idx(p)
      character(*) :: p
      integer      :: idx, i, j
      do idx = 1, size(field_list)
          if (field_list(idx)%label == p) return
      end do
      if (p(1:min(3,len(p))) == 'AUX') then
         read (p(min(4,len(p)):),'(i1)', iostat=j) i
         if (j /= 0 .or. i < 0 .or. i > 9) &
            call lsquit('rsp_contribs error: ' &
                     // 'Auxiliary label index should be within 0..9, ' // p,-1)
         idx = 1
         return
      end if
      call lsquit('Perturbation not found: ' // p,-1)
   end function


   function pert_antisym(p)
     implicit none
      character(*), intent(in) :: p(:)
      logical :: pert_antisym(size(p))
      integer :: i
      i=0
      pert_antisym = (/(field_list(idx(p(i)))%anti, i=1,size(p))/)
   end function


   !> shape (dimensions) of property p(:)
   function pert_shape(mol, p)
     implicit none
      !> structure containing the integral program settings
      type(rsp_molcfg),  intent(in) :: mol
      !> field lables
      character(*),      intent(in) :: p(:)
      integer :: pert_shape(size(p)), i,sizep
      i=0
      sizeP = size(p)
      IF(size(p).GE.1)THEN
         IF(p(1).EQ.'NOOP')THEN
            sizeP = 0
            pert_shape(1) = 1
         ENDIF
      ENDIF
      pert_shape = (/(field_list(idx(p(i)))%ncomp, i=1,sizeP)/)
      ! loop through mol-dependent
      do i=1, sizeP
         if (pert_shape(i) /= -1) then
            ! cycle
         else if (p(i) == 'GEO') then
            pert_shape(i) = 3 * mol%natoms
         else if (p(i) == 'GEO' .or. p(i) == 'NUCM') then
            pert_shape(i) = 3 * mol%natoms
         else
            call lsquit('pert_shape error: Number of comp. unknown for ' // p(i), mol%lupri)
         end if
      end do
   end function


   function pert_basdep(p)
     implicit none
      character(*), intent(in) :: p(:)
      logical :: pert_basdep(size(p))
      integer :: i
      i=0
      pert_basdep = (/(field_list(idx(p(i)))%bas, i=1,size(p))/)
   end function


   !> Find the reordering of fi(:) that puts it in "standard" order:
   !> sorted by 1) decreasing %label's index in field_list (GEO before EL),
   !>           2) increasing component index %comp,
   !>           3) decreasing absolute value of real part
   !>              of %freq (- before +)
   !>           4) decreasing absolute value of imaginary
   !>              part of %freq (- before +)
   function pert_order(fld) result(ord)
     implicit none
      type(pert_field), intent(in) :: fld(:)
      integer :: ord(size(fld)), i, j, k
      i=0
      ord = (/(i, i=1, size(fld))/)
      do i=1, size(fld)
         !find perturbation after i, with highest idxp (secondly highest ddime)
         j = i
         do k = j+1, size(fld)
            if (idx(fld(ord(k))%label) <  idx(fld(ord(j))%label)) cycle
            if (idx(fld(ord(k))%label) == idx(fld(ord(j))%label)) then
               if (fld(ord(k))%comp >  fld(ord(j))%comp) cycle
               if (fld(ord(k))%comp == fld(ord(j))%comp) then
                  if (   abs(real(fld(ord(k))%freq)) &
                      <  abs(real(fld(ord(j))%freq))) cycle
                  if (   abs(real(fld(ord(k))%freq)) &
                      == abs(real(fld(ord(j))%freq))) then
                     if (   abs(aimag(fld(ord(k))%freq)) &
                         <= abs(aimag(fld(ord(j))%freq))) cycle
                  end if
               end if
            end if
            j = k  !new minimum
         end do
         !swap entries i and j
         k = ord(i)
         ord(i) = ord(j)
         ord(j) = k
      end do
   end function



   !> Nuclei-nuclei or nuclei-field contributions to properties
   !> ajt This should be moved to separate module
   !>     prop_intifc, together with prop_intifc_2el,
   !>     prop_intifc_ksm and prop_intifc_1el
   subroutine prop_intifc_nuc(mol, flds, pnuc)
      !> structure containing the integral program settings
      type(rsp_molcfg),  intent(in)  :: mol
      !> field lables, must currently be either EL or ELGR
      character(4),      intent(in)  :: flds(:)
      !> resulting integral averages
      real(realk),       intent(inout) :: pnuc(:)
      !---------------------------------------
      if (size(flds) /= 1) &
         call lsquit('prop_intifc_nuc: size(flds) /= 1, but only ' &
                // 'EL and GEO are currently implemented', mol%lupri)
      ! verify label and size
      select case(flds(1))
      case ('EL')
         ! ajt-to-Thomas: insert direct call in rsp_oneave, and remove this
         if (size(pnuc) /= 3) &
            call lsquit('prop_intifc_nuc: For field EL, size(pnuc) must be 3', mol%lupri)
         call II_get_nucdip(mol%setting, pnuc)
         pnuc = -pnuc
      case ('GEO ')
         ! ajt-to-Thomas: insert direct call in rsp_oneave, and remove this
         if (size(pnuc) /= 3 * mol%natoms) &
            call lsquit('prop_intifc_nuc: For field ELGR, size(pnuc) must be 6', mol%lupri)
         call II_get_nn_gradient(pnuc, mol%setting, mol%lupri, mol%luerr)
      case default
         call lsquit('prop_intifc_nuc: Unimplemented or unknown field "' // flds(1) &
                // '". Only EL (-nucdip), ELGR (-qdrnuc) and GEO (gradnn) ' &
                // 'are implemented so far', mol%lupri)
      end select
   end subroutine



   !> Call one-electron integral program to calculate some
   !> integrals ints(:) and some averages avgs(:).
   !> ajt This should eventually be moved to separate module
   !>     prop_intifc, together with prop_intifc_2el,
   !>     prop_intifc_ksm and prop_intifc_nuc
   !> ajt 1011 This should be removed, as its functionality is supplied by
   !>     II_get_prop in module Integral_interface
   subroutine prop_intifc_1el(mol, flds, dens, avgs, ints)
      type lsint_arg
         integer       :: idens
         integer       :: dim1
         integer       :: dim2
         integer       :: dim5
         integer       :: ao2
         integer       :: ao3
         integer       :: nderiv
         integer       :: oper
         integer       :: spec
         integer       :: fac
      end type
      !> structure containing the integral program settings
      type(rsp_molcfg),  intent(in)    :: mol
      !> field lables, must currently be either (/'EL'/) or (/'ELGR'/)
      character(4),      intent(in)    :: flds(:)
      !> density matrices to contract integrals with
      type(matrix),      intent(in)    :: dens(:)
      !> resulting integral averages
      real(realk),       intent(inout)   :: avgs(:)
      !> resulting integral matrices
      type(matrix),      intent(inout) :: ints(:)
      !----------------------------------------------
      type(LSSETTING), pointer :: st
      type(lsint_arg) :: run(3)
      integer         :: nf, nb, na, nrun, nmat, i, j
      real(realk)     :: tmp(3, (size(avgs)+2)/3)
      type(matrix)    :: mat(10) !scratch including lower-order integral matrices
      type(matrix)    :: D_AO
      nf = size(flds)
      nb = mol%zeromat%nrow
      na = mol%natoms
      st => mol%setting
      ! verify field labels flds(:), while configuring run / num_runs
      nrun = 1  !default
      if (nf==0) then
         call lsquit('prop_intifc_1el: size(flds) == 0 (unperturbed) not implemented', mol%lupri)
      else if (nf==1 .and. flds(1)=='ELGR') then
         ! ajt-to-Thomas: insert direct call to II_get_prop from oneint and oneave, and remove this
         if (size(ints) /= 6) &
            call lsquit('prop_intifc_1el: For field ELGR, size(ints) must be 6', mol%lupri)
         run(1) = lsint_arg(-1,nb,nb,10,AORdefault,AOEmpty,2,CarmomOperator,RegularSpec,1)
      else if (nf==1 .and. flds(1)=='GEO ') then
         ! ajt-to-Thomas: insert direct call to II_get_oneElectron_gradient from oneave, and remove this
         if (size(dens) /= 2) &
            call lsquit('prop_intifc_1el: For field GEO, only size(dens)=2 implemented', mol%lupri)
         if (size(avgs) /= 3*na) &
            call lsquit('prop_intifc_1el: For field GEO, size(avgs) must be 3*natoms', mol%lupri)
         run(1) = lsint_arg(1,3,na,1,AORdefault,AONuclear,1,NucpotOperator, GradientSpec, 2)
         run(2) = lsint_arg(1,3,na,1,AOEmpty,  AORdefault,2,KineticOperator,GradientSpec, 2)
         run(3) = lsint_arg(2,3,na,1,AORdefault,AOEmpty,  1,OverlapOperator,GradientSpec,-2)
         nrun = 3
      else
         call lsquit("prop_intifc_1el: unknown field argument 'flds':" // flds(1), mol%lupri)
      end if
      ! if there are any integral runs, prepare matrices
      nmat = 0
      if (any(run(:nrun)%idens == -1)) then
         ! zero output integral matrices (before accumulating)
         do i=1, size(ints)
            ints(i) = 0 * mol%zeromat
         end do
         ! allocate temporary integral matrices
         nmat = maxval(run(:nrun)%dim5)
         do i=1, nmat
            mat(i) = 0*mol%zeromat
            call mat_ensure_alloc(mat(i))
         end do
      end if
      ! loop over calls to integral program LSint
      ! KK & AJT: Initialize avgs
      avgs=0
      do i=1, nrun
         ! if this is an average run (not an integral run), register density
         ! matrix to average over, in setting
         if (run(i)%idens /= -1)then
            IF(st%IntegralTransformGC)THEN
               D_AO = 0*mol%zeromat
               call mat_ensure_alloc(D_AO)
               call GCAO2AO_transform_matrixD2(dens(run(i)%idens),D_AO,st,mol%lupri)
               call ls_attachDmatToSetting(D_AO, 1, st,'LHS', 1, &
                    merge(2, 3, run(i)%ao2 == AORdefault), .TRUE.,mol%lupri)
            ELSE
               call ls_attachDmatToSetting(dens(run(i)%idens), 1, st,'LHS', 1, &
                    merge(2, 3, run(i)%ao2 == AORdefault), .TRUE.,mol%lupri)
            ENDIF
         endif
         ! set up dimensions of output
         call initIntegralOutputDims(st%output, run(i)%dim1, &
                                     run(i)%dim2, 1, 1, run(i)%dim5)
         ! run integral program
         st%scheme%cmorder = run(i)%nderiv
         st%scheme%intTHRESHOLD = st%SCHEME%THRESHOLD*st%SCHEME%ONEEL_THR
         call ls_getIntegrals(AORdefault, run(i)%ao2, run(i)%ao3, AOEmpty, &
                              run(i)%oper,run(i)%spec, Contractedinttype, st, &
                              mol%lupri, mol%luerr)
         st%scheme%cmorder = 0
         if (run(i)%idens /= -1)then
            IF(st%IntegralTransformGC)THEN
               D_AO=0
            ENDIF
         endif
         ! retrieve averages or integrals
         if (run(i)%idens /= -1) then
            ! retrieve and accumulate averages
            call retrieve_output(mol%lupri, st, tmp, st%IntegralTransformGC)
            avgs = avgs + run(i)%fac * (/tmp/)
            ! unregister averaged density matrix
            call ls_freeDmatFromSetting(st)
         else
            ! retrieve and accumulate integrals
            call retrieve_output(mol%lupri, st, mat(:run(i)%dim5), st%IntegralTransformGC)
            do j=1, size(ints)
               ints(j) = ints(j) + mat(j + run(i)%dim5 - size(ints))
            end do
         end if
      end do
      ! free temporary matrices
      mat(:nmat) = 0
   end subroutine



   !> Call 2-electron integral program.
   !> FIXME Should support several avg / fock contractions at once
   subroutine prop_intifc_2el(mol, flds, dens, avgs, fock)
      !> structure containing the integral program settings
      type(rsp_molcfg),     intent(in) :: mol
      !> field lables, must currently be (/'GEO'/)
      character(4),         intent(in) :: flds(:)
      !> density matrices to contract integrals with
      type(matrix), target, intent(in) :: dens(:)
      !> resulting integral averages
      real(realk),          intent(inout) :: avgs(:)
      !> resulting integral matrices
      type(matrix),      intent(inout) :: fock(:)
      !--------------------------------------------
      type(matrix), target :: densT
      type(LSSETTING), pointer :: setting
      integer :: nb, na
      nb = mol%zeromat%nrow
      na = mol%natoms
      setting => mol%setting
      if (size(flds) /= 1) &
         call lsquit('prop_intifc_2el: size(flds) /= 1, but only GEO implemented', mol%lupri)
      ! verify field labels flds(:), while configuring
      select case(flds(1))
      case ('GEO')
         if (size(dens) /= 1 .and. size(dens) /= 2) &
            call lsquit('prop_intifc_1el: For field GEO, only size(dens)=1,2 implemented', mol%lupri)
      case default
         call lsquit('prop_intifc_1el: Unimplemented or unknown field "' // flds(1) &
                // '". Only EL (DIPLEN), ELGR (SECMOM) and GEO are implemented', mol%lupri)
      end select
      !ajt FIXME Set up and call Coulomb and Exchange directly

      ! KK and AJT: Send in transposed matrix due to conventions twoelectron_gradient routine.
      densT = trans(dens(size(dens)))
      call II_get_twoElectron_gradient(avgs, na, (/Matrixp(dens(1))/), &
                                       (/Matrixp(densT)/), 1, 1, &
                                       setting, mol%lupri, mol%luerr)
      ! factor 8, it seems
      avgs = avgs * 8

      ! Free temporary transposed matrix
      densT=0

   end subroutine


end module
