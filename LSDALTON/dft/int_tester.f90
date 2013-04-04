! Program for testing XC integration performance and correctness.
!
! The aim is to load a benchmark molecule, set up the AO information
! and then directly proceed to XC integration.  By omitting the SCF
! setup, coulom integral evaluation, this will speed up considerably
! benchmarking the XC code for a range of molecules and settings,
! allowing for a systematic studies and maybe even development of
! automatically tuned code.
!
! To compile it, type "make int_tester" in the main dalton directory.
!
program int_tester
  implicit none
  double precision, allocatable :: work(:)
  double precision, allocatable :: dmat(:,:), fksm(:,:)
  integer :: lwork, iprfck, nbast
  double precision :: edfty

  lwork = 100000000
  nbast = 100 !TK: I set the variable. TK
  allocate(work(lwork))
  call herini
  call readin(WORK,LWORK,.false.)
  deallocate(work)
  call dftsetfunc("BP86")
  lwork = NBAST*NBAST*4
  iprfck=1

  allocate(work(NBAST*NBAST*4))
  allocate(dmat(nbast,nbast), fksm(nbast,nbast))
  call DFTKSMb(DMAT,FKSM,EDFTY,WORK,LWORK,IPRFCK)
  deallocate(work)
  print *, "Computed XC energy ", EDFTY
end program int_tester
