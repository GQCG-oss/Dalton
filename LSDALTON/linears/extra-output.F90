!> @file
!> Contains the module extra_output
!
!> brief Routines that make extra output if requested
!> author C. Nygaard
!> date June 2012
!
!> Basically this module was made to contain the subroutine print_orbital_info2
!>  because I didn't know where else to put it
!> If anyone else knows where to put this routine, feel free to move it
!> If anyone else have routines that fit into this module, feel free to add them,
!>  just remember to change this documentation then

module extra_output

use files
use matrix_operations

contains

!> \brief Find and print MO energies and occupations, print MO's to file
!> \author C. Nygaard
!> \date June 12. 2012
!> \param Dao The density matrix in the AO basis
!> \param Fao The Fock matrix in the AO basis
!> \param S The AO overlap matrix
!> \param filename Name of the file the orbitals are printed to
!> \param unres True if calculation is unrestricted, false otherwise
!> \param ps True if the solutions is expected to be pure-state, false otherwise
!> \param lupri LUN for standard output
!
!> An aufbau solution is not assumed when the orbitals are created.
!==============================================================================
subroutine print_orbital_info2(Dao, Fao, S, filename, unres, ps, lupri)

implicit none

type(matrix), intent(in)     :: Dao, Fao, S
character(len=*), intent(in) :: filename
logical, intent(in)          :: unres, ps
integer, intent(in)          :: lupri
integer                      :: lud

logical :: OnMaster
OnMaster=.TRUE.

if (ps) then
  call print_orbital_info_ps (Dao, Fao, S, filename, unres, lupri)
else
  call print_orbital_info_es (Dao, Fao, S, filename, unres, lupri)
endif

!output DAO to file
lud = -1
call lsopen (lud, "dao.out", "UNKNOWN", "UNFORMATTED")
call mat_write_to_disk (lud, Dao, OnMaster)
call lsclose (lud, "KEEP")

end subroutine print_orbital_info2
!==============================================================================

!> \brief Find and print MO energies and occupations, print MO's to file
!> \author C. Nygaard
!> \date June 10. 2012
!> \param Dao The density matrix in the AO basis
!> \param Fao The Fock matrix in the AO basis
!> \param S The AO overlap matrix
!> \param filename Name of the file the orbitals are printed to
!> \param unres True if calculation is unrestricted, false otherwise
!> \param lupri LUN for standard output
!
!> Finds molecular orbitals by first diagonalising the AO density matrix,
!> then diagonalising the Fock matrix projected in the occupied space and 
!> in the virtual space each.
!==============================================================================
subroutine print_orbital_info_ps (Dao, Fao, S, filename, unres, lupri)

implicit none

type(matrix), intent(in)     :: Dao, Fao, S
character(len=*), intent(in) :: filename
logical, intent(in)          :: unres
integer, intent(in)          :: lupri

integer                      :: i, luc, luoccs
integer                      :: Nbast, Nocc, Nvirt, Ndim
type(matrix)                 :: SDS, Cmo, Cocc, Cvirt, Focc, Fvirt
real(realk), pointer         :: occ(:), orbE(:)
type(matrix)                 :: tmp, eivec
real(realk), pointer         :: eival(:)

logical :: OnMaster
OnMaster=.TRUE.

!initialisations
Nocc = nint(mat_TrAB (Dao, S))
Nbast = S%nrow
if (unres) then
  Nocc = Nocc/2
  Ndim = 2*Nbast
else
  Ndim = Nbast
endif
Nvirt = Nbast - Nocc
call mem_alloc (occ,Ndim)
call mem_alloc (orbE,Ndim)
call mat_init (Cmo, Nbast, Nbast)
call mat_init (SDS, Nbast, Nbast)
call mat_init (Cocc, Nbast, Nocc)
call mat_init (Focc, Nocc, Nocc)
call mat_init (Cvirt, Nbast, Nvirt)
call mat_init (Fvirt, Nvirt, Nvirt)

call mat_init (tmp, Nbast, Nbast)
call mat_mul (S, Dao, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
call mat_mul (tmp, S, 'n', 'n', -1.0E0_realk, 0.0E0_realk, SDS)
           !minus so the occupied orbitals to comes first
call mat_free (tmp)

call mat_diag_f (SDS, S, occ, Cmo)
occ = -occ !minus to get positive occupation numbers

!Check if the occupation numbers are correct (either zero or one)
do i=1,Nocc
  if (occ(i)-1.0E0_realk < -1.0E-8_realk .or. occ(i)-1.0E0_realk > 1.0E-8_realk) then
    write (lupri, *) 'Orbital should be occupied, but has occupation number', occ(i)
    call lsquit ('Something is wrong with the occupations', lupri)
  endif
  if (unres) then
    if (occ(i+Nbast)-1.0E0_realk < -1.0E-8_realk .or. occ(i+Nbast)-1.0E0_realk > 1.0E-8_realk) then
      write (lupri, *) 'Orbital should be occupied, but has occupation number', occ(i+Nbast)
      call lsquit ('Something is wrong with the occupations', lupri)
    endif
  endif
enddo
do i=Nocc+1,Nbast
  if (occ(i) < -1.0E-8_realk .or. occ(i) > 1.0E-8_realk) then
    write (lupri, *) 'Orbital should be empty, but has occupation number', occ(i)
    call lsquit ('Something is wrong with the occupations', lupri)
  endif
  if (unres) then
    if (occ(i+Nbast) < -1.0E-8_realk .or.  occ(i+Nbast) > 1.0E-8_realk) then
      write (lupri, *) 'Orbital should be empty, but has occupation number', occ(i+Nbast)
      call lsquit ('Something is wrong with the occupations', lupri)
    endif
  endif
enddo

!Make the occupied part of Fock matrix: Focc = Cocc^T Fao Cocc
call mat_init (tmp, Nocc, Nbast)
call mat_section (Cmo, 1, Nbast, 1, Nocc, Cocc)
call mat_mul (Cocc, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
call mat_mul (tmp, Cocc, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Focc)
call mat_free (tmp)
!Diagonalise Focc
call mat_init (tmp, Nbast, Nocc)
call mat_init (eivec, Nocc, Nocc)
if (unres) then
  call mem_alloc (eival,2*Nocc)
else
  call mem_alloc (eival,Nocc)
endif
call mat_assign(eivec,Focc)
call mat_dsyev (eivec, eival, Nocc)
call mat_mul (Cocc, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
call mat_assign(Cocc,tmp)
call mat_insert_section (Cocc, 1, Nbast, 1, Nocc, Cmo)
orbE(1:Nocc) = eival(1:Nocc)
if (unres) orbE(Nbast+1:Nbast+Nocc) = eival(Nocc+1:2*Nocc)
call mem_dealloc (eival)
call mat_free (eivec)
call mat_free (tmp)

!Make the unoccupied part of Fock matrix: Fvirt = Cvirt^T Fao Cvirt
call mat_init (tmp, Nvirt, Nbast)
call mat_section (Cmo, 1, Nbast, Nocc+1, Nbast, Cvirt)
call mat_mul (Cvirt, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
call mat_mul (tmp, Cvirt, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Fvirt)
call mat_free (tmp)
!Diagonalise Fvirt
call mat_init (tmp, Nbast, Nvirt)
call mat_init (eivec, Nvirt, Nvirt)
if (unres) then
  call mem_alloc (eival,2*Nvirt)
else
  call mem_alloc (eival,Nvirt)
endif
call mat_assign(eivec,Fvirt)
call mat_dsyev (eivec, eival, Nvirt)
call mat_mul (Cvirt, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
call mat_assign(Cvirt,tmp)
call mat_insert_section (Cvirt, 1, Nbast, Nocc+1, Nbast, Cmo)
orbE(Nocc+1:Nbast) = eival(1:Nvirt)
if (unres) orbE(Nbast+Nocc+1:2*Nbast) = eival(Nvirt+1:2*Nvirt)
call mem_dealloc (eival)
call mat_free (eivec)
call mat_free (tmp)

!output Cmo to file
luc = -1
call lsopen (luc, filename, "UNKNOWN", "UNFORMATTED")
call mat_write_to_disk (luc, Cmo, OnMaster)
call lsclose (luc, "KEEP")

!output occs and orbe to file
luoccs = -1
call lsopen (luoccs, "occs-orbe.out", "UNKNOWN", "FORMATTED")
write (luoccs, *) "#", "orbital", "occupation", "orbital energy"
do i=1,Ndim
  write (luoccs, *) i, occ(i), orbE(i)
enddo
call lsclose (luoccs, "KEEP")

!output orbital energies and occupations to standard output
if (unres) then
  write (lupri, *) 'Orbital energies of occupied alpha-orbitals:'
  write (lupri, *) orbE(1:Nocc)
  write (lupri, *) 'Orbital energies of occupied beta-orbitals:'
  write (lupri, *) orbE(Nbast+1:Nbast+Nocc)
  write (lupri, *) 'Orbital energies of virtual alpha-orbitals:'
  write (lupri, *) orbE(Nocc+1:Nbast)
  write (lupri, *) 'Orbital energies of virtual beta-orbitals:'
  write (lupri, *) orbE(Nbast+Nocc+1:2*Nbast)
else
  write (lupri, *) 'Orbital energies of occupied orbitals:'
  write (lupri, *) orbE(1:Nocc)
  write (lupri, *) 'Orbital energies of virtual orbitals:'
  write (lupri, *) orbE(Nocc+1:Nbast)
endif

!write (lupri, *)
!write (lupri, *) 'Cmo:'
!call mat_print (Cmo, 1, Nbast, 1, Nbast, lupri)
!write (lupri, *)

!closing down
call mat_free (Cmo)
call mat_free (SDS)
call mat_free (Cocc)
call mat_free (Focc)
call mat_free (Cvirt)
call mat_free (Fvirt)
call mem_dealloc (occ)
call mem_dealloc (orbE)

end subroutine print_orbital_info_ps
!==============================================================================

!> \brief Find and print MO energies and occupations, print MO's to file
!> \author C. Nygaard
!> \date June 10. 2012
!> \param Dao The density matrix in the AO basis
!> \param Fao The Fock matrix in the AO basis
!> \param S The AO overlap matrix
!> \param filename Name of the file the orbitals are printed to
!> \param unres True if calculation is unrestricted, false otherwise
!> \param lupri LUN for standard output
!
!> Finds molecular orbitals by first diagonalising the AO density matrix,
!> then diagonalising the Fock matrix projected in the occupied space and 
!> in the virtual space each.
!==============================================================================
subroutine print_orbital_info_es (Dao, Fao, S, filename, unres, lupri)

use matrix_operations_aux, only: mat_get_elm

implicit none

type(matrix), intent(in)     :: Dao, Fao, S
character(len=*), intent(in) :: filename
logical, intent(in)          :: unres
integer, intent(in)          :: lupri

integer                      :: i, j, luc, luoccs
integer                      :: Nbast, Nocc, Nact, Nvirt, Ndim
integer                      :: Noccb, Nactb, Nvirtb
type(matrix)                 :: SDS, Cmo, Cocc, Cvirt, Cact, Focc, Fact, Fvirt, SDSact
real(realk), pointer         :: occ(:), orbE(:)
type(matrix)                 :: tmp, eivec
real(realk), pointer         :: eival(:)
real(realk)                  :: efermi
real(realk)                  :: check

logical :: OnMaster
OnMaster=.TRUE.

!initialisations
Nocc = nint(mat_TrAB (Dao, S))
Nbast = S%nrow
if (unres) then
  Nocc = Nocc/2
  Ndim = 2*Nbast
else
  Ndim = Nbast
endif
Nvirt = Nbast - Nocc
call mem_alloc (occ,Ndim)
call mem_alloc (orbE,Ndim)
call mat_init (Cmo, Nbast, Nbast)
call mat_init (SDS, Nbast, Nbast)

call mat_init (tmp, Nbast, Nbast)
call mat_mul (S, Dao, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
call mat_mul (tmp, S, 'n', 'n', -1.0E0_realk, 0.0E0_realk, SDS)
           !minus so the occupied orbitals to comes first
call mat_free (tmp)

call mat_diag_f (SDS, S, occ, Cmo)
occ = -occ !minus to get positive occupation numbers

!Create spaces
Nocc = 0 ; Nact = 0 ; Nvirt = 0
Noccb = 0 ; Nactb = 0 ; Nvirtb = 0
do i=1,Nbast
  if (occ(i)-1.0E0_realk < 1.0E-8_realk .and. occ(i)-1.0E0_realk > -1.0E-8_realk) then
    !orbital is occupied
    Nocc = Nocc + 1
  elseif (occ(i) < 1.0E-8_realk .and. occ(i) > -1.0E-8_realk) then
    !orbital is active
    Nvirt = Nvirt + 1
  else
    !orbital is active
    Nact = Nact + 1
  endif
  if (unres) then
    if (occ(i+Nbast) < 1.0E-8_realk .and.  occ(i+Nbast) > -1.0E-8_realk) then
      !orbital is occupied
      Noccb = Nocc + 1
    elseif (occ(i+Nbast)-1.0E0_realk < 1.0E-8_realk .and. occ(i+Nbast)-1.0E0_realk > -1.0E-8_realk) then
      !orbital is virtual
      Nvirt = Nvirt + 1
    else
      !orbital is active
      Nact = Nact + 1
    endif
  endif
enddo
if (Nocc+Nact+Nvirt /= Nbast) then
  call lsquit ("Something is wrong in print_orbital_info_es", lupri)
endif
if (unres) then
  Nocc = min(Nocc, Noccb)
  Nvirt = min(Nvirt, Nvirtb)
  Nact = Nbast - Nocc - Nvirt
endif

if (Nact == 0) then
  call print_orbital_info_ps (Dao, Fao, S, filename, unres, lupri)
else

  call mat_init (Cocc, Nbast, Nocc)
  call mat_init (Focc, Nocc, Nocc)
  call mat_init (Cact, Nbast, Nact)
  call mat_init (Fact, Nact, Nact)
  call mat_init (Cvirt, Nbast, Nvirt)
  call mat_init (Fvirt, Nvirt, Nvirt)
  call mat_init (SDSact, Nact, Nact)
  
  !Make the occupied part of Fock matrix: Focc = Cocc^T Fao Cocc
  if (Nocc > 0) then
    call mat_init (tmp, Nocc, Nbast)
    call mat_section (Cmo, 1, Nbast, 1, Nocc, Cocc)
    call mat_mul (Cocc, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_mul (tmp, Cocc, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Focc)
    call mat_free (tmp)
    !Diagonalise Focc
    call mat_init (tmp, Nbast, Nocc)
    call mat_init (eivec, Nocc, Nocc)
    if (unres) then
      call mem_alloc (eival,2*Nocc)
    else
      call mem_alloc (eival,Nocc)
    endif
    call mat_assign(eivec,Focc)
    call mat_dsyev (eivec, eival, Nocc)
    call mat_mul (Cocc, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_assign(Cocc,tmp)
    call mat_insert_section (Cocc, 1, Nbast, 1, Nocc, Cmo)
    orbE(1:Nocc) = eival(1:Nocc)
    if (unres) orbE(Nbast+1:Nbast+Nocc) = eival(Nocc+1:2*Nocc)
    call mem_dealloc (eival)
    call mat_free (eivec)
    call mat_free (tmp)
  endif
  
  if (Nvirt > 0) then
    !Make the unoccupied part of Fock matrix: Fvirt = Cvirt^T Fao Cvirt
    call mat_init (tmp, Nvirt, Nbast)
    call mat_section (Cmo, 1, Nbast, Nocc+Nact+1, Nbast, Cvirt)
    call mat_mul (Cvirt, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_mul (tmp, Cvirt, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Fvirt)
    call mat_free (tmp)
    !Diagonalise Fvirt
    call mat_init (tmp, Nbast, Nvirt)
    call mat_init (eivec, Nvirt, Nvirt)
    if (unres) then
      call mem_alloc (eival,2*Nvirt)
    else
      call mem_alloc (eival,Nvirt)
    endif
    call mat_assign(eivec,Fvirt)
    call mat_dsyev (eivec, eival, Nvirt)
    call mat_mul (Cvirt, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_assign(Cvirt,tmp)
    call mat_insert_section (Cvirt, 1, Nbast, Nocc+Nact+1, Nbast, Cmo)
    orbE(Nocc+Nact+1:Nbast) = eival(1:Nvirt)
    if (unres) orbE(Nbast+Nocc+Nact+1:2*Nbast) = eival(Nvirt+1:2*Nvirt)
    call mem_dealloc (eival)
    call mat_free (eivec)
    call mat_free (tmp)
  endif
  
  if (Nact > 0) then
    !Make the active part of Fock matrix: Fact = Cact^T Fao Cact
    call mat_init (tmp, Nact, Nbast)
    call mat_section (Cmo, 1, Nbast, Nocc+1, Nocc+Nact, Cact)
    call mat_mul (Cact, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_mul (tmp, Cact, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Fact)
    call mat_free (tmp)
    !Diagonalise Fact
    call mat_init (tmp, Nbast, Nact)
    call mat_init (eivec, Nact, Nact)
    if (unres) then
      call mem_alloc (eival,2*Nact)
    else
      call mem_alloc (eival,Nact)
    endif
    call mat_assign(eivec,Fact)
    call mat_dsyev (eivec, eival, Nact)
    call mat_mul (Cact, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_assign(Cact,tmp)
    call mat_free (eivec)
    call mat_free (tmp)
    !Check if all active orbitals have same energy
    efermi = eival(1)
    do i=2,Nact
      if (efermi-eival(i) < -1.0E-5_realk .or. efermi-eival(i) > 1.0E-5_realk) then
        write (lupri, *) 'WARNING: Active orbitals do not have the same energy!'
      endif
    enddo
    orbE(Nocc+1:Nocc+Nact) = eival(1:Nact)
    if (unres) orbE(Nbast+Nocc+1:Nbast+Nocc+Nact) = eival(Nact+1:2*Nact)
    call mem_dealloc (eival)
    !Make the active part of SDS
    call mat_init (tmp, Nact, Nbast)
    call mat_mul (Cact, SDS, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_mul (tmp, Cact, 'n', 'n', 1.0E0_realk, 0.0E0_realk, SDSact)
    call mat_free (tmp)
    !Diagonalise SDSact
    call mat_init (tmp, Nbast, Nact)
    call mat_init (eivec, Nact, Nact)
    if (unres) then
      call mem_alloc (eival,2*Nact)
    else
      call mem_alloc (eival,Nact)
    endif
    call mat_assign(eivec,SDSact)
    call mat_dsyev (eivec, eival, Nact)
    call mat_mul (Cact, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
    call mat_assign(Cact,tmp)
    call mat_free (eivec)
    call mat_free (tmp)
    !Create active part of Cmo and occupations
    call mat_insert_section (Cact, 1, Nbast, Nocc+1, Nocc+Nact, Cmo)
    occ(Nocc+1:Nocc+Nact) = -eival(1:Nact)
    !debug
    if (.false.) then
      call mat_init (tmp, Nact, Nbast)
      call mat_mul (Cact, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
      call mat_mul (tmp, Cact, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Fact)
      call mat_free (tmp)
      do i=Nocc+1,Nocc+Nact
        call mat_get_elm (Fact, i, i, orbE(i))
        do j=Nocc+1,Nocc+Nact
          if (i /= j) then
            call mat_get_elm (Fact, i, j, check)
            write (lupri, *) "off-diagonal elements of Fact =", check
          endif
        enddo
      enddo
    endif
    !end debug
    if (unres) occ(Nbast+Nocc+1:Nbast+Nocc+Nact) = eival(Nact+1:2*Nact)
    call mem_dealloc (eival)
  endif
  
  !output Cmo to file
  luc = -1
  call lsopen (luc, filename, "UNKNOWN", "UNFORMATTED")
  call mat_write_to_disk (luc, Cmo, OnMaster)
  call lsclose (luc, "KEEP")
 
  !output occs and orbe to file
  luoccs = -1
  call lsopen (luoccs, "occs-orbe.out", "UNKNOWN", "FORMATTED")
  write (luoccs, *) "#", "orbital", "occupation", "orbital energy"
  do i=1,Ndim
    write (luoccs, *) i, occ(i), orbE(i)
  enddo
  call lsclose (luoccs, "KEEP")
 
  !output orbital energies and occupations to standard output
  if (unres) then
    write (lupri, *) 'Orbital energies of occupied alpha-orbitals:'
    write (lupri, *) orbE(1:Nocc)
    write (lupri, *) 'Orbital energies of occupied beta-orbitals:'
    write (lupri, *) orbE(Nbast+1:Nbast+Nocc)
    write (lupri, *) 'Orbital occupations and energies of active alpha-orbitals:'
    do i=Nocc+1, Nocc+Nact
      write (lupri, *) i, occ(i), orbE(i)
    enddo
    write (lupri, *) 'Orbital occupations and energies of active beta-orbitals:'
    do i=Nbast+Nocc+1, Nbast+Nocc+Nact
      write (lupri, *) i, occ(i), orbE(i)
    enddo
    write (lupri, *) 'Orbital energies of virtual alpha-orbitals:'
    write (lupri, *) orbE(Nocc+Nact+1:Nbast)
    write (lupri, *) 'Orbital energies of virtual beta-orbitals:'
    write (lupri, *) orbE(Nbast+Nocc+Nact+1:2*Nbast)
  else
    write (lupri, *) 'Orbital energies of occupied orbitals:'
    write (lupri, *) orbE(1:Nocc)
    write (lupri, *) 'Orbital occupations and energies of active orbitals:'
    do i=Nocc+1, Nocc+Nact
      write (lupri, *) i, occ(i), orbE(i)
    enddo
    write (lupri, *) 'Orbital energies of virtual orbitals:'
    write (lupri, *) orbE(Nocc+Nact+1:Nbast)
  endif
 
!  write (lupri, *)
!  write (lupri, *) 'Cmo:'
!  call mat_print (Cmo, 1, Nbast, 1, Nbast, lupri)
!  write (lupri, *)
 
  !closing down
  call mat_free (Cocc)
  call mat_free (Focc)
  call mat_free (Cact)
  call mat_free (Fact)
  call mat_free (Cvirt)
  call mat_free (Fvirt)
  call mat_free (SDSact)

endif

call mat_free (Cmo)
call mat_free (SDS)
call mem_dealloc (occ)
call mem_dealloc (orbE)

end subroutine print_orbital_info_es
!==============================================================================

end module extra_output
