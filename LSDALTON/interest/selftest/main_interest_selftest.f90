! -------------------------------------------------------------------------
!
! Program:      Interest
!
! File:         main_interest_selftest.f90 
!
! Description:  Performing selftest calculations of all integrals implemented in src/
!
! Licensing:    This code is distributed under the GNU LGPL license
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Revisions:    
!
! -------------------------------------------------------------------------
!program interest_selftest
program main

  use module_interest_selftest 

  implicit none

  write(6,'(a)') 
  write(6,'(a)') 'Interest is performing a selftest ...'
  write(6,'(a)') 

  call interest_initialize()
  call interest_selftest_overlap()
  call interest_selftest_dipole()
  call interest_selftest_nuclear_attraction_point_nucleus()

end program
