!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
module mlcc3_init
!
!
!  mlcc3 initialiasion routines
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
   use mlcc3_active_spaces
   use mlcc_typedef
   use mlcc3_data
   use mlcc_block_import
   use mlcc_work
!
!   
contains
   subroutine sirius_reader
!
!  Sirius reader routine
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: read in data from SIRIUS etc. and store in arrays
!
      implicit none
!
      integer  :: lusifc = -1 
      integer  :: idummy = 0
      integer  :: n_symmetries, n_basis_sym, n_orbitals_sym
      integer  :: i
!      
!
!     
!     Open Sirius Fock file
!     ---------------------
!
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',idummy,'.FALSE.')
      rewind(lusifc)
!
!
!     Read in various stuff from Sirius Fock file. Things depending on symmetry is mostly
!     discarded at the end of the subroutine as we do not use symmetry. Information in 
!     file should be in Cholesky orbital format. If Cholesky orbitals has been generated,
!     SIRIFC will contain the data in the Cholesky basis
!
      call mollab('TRCCINT ',lusifc,lupri)
!      
      read(lusifc) n_symmetries, n_orbitals, n_basis, n_lambda, n_occ, &
      &            n_orbitals_sym, n_basis_sym, nuclear_potential, scf_energy
!
!      
      if (n_symmetries /= 1) call quit('error in mlcc3_init: MLCC not implemented with symmetry')
!      
!      
!     Calculate number of virtuals and amplitudes
!
      n_vir          = n_orbitals - n_occ
      n_t1am         = n_vir*n_occ
      n_t2am         = n_t1am*n_t1am
      n_t2am_pack    = n_t1am*(n_t1am+1)/2
      n_v_2          = n_vir**2
      n_v_3          = n_vir**3
      n_basis_2      = n_basis**2
      n_basis_2_pack = n_basis*(n_basis+1)/2
      n_ao_ints      = n_basis*n_basis_2_pack
      n_bas_orb      = n_basis*n_orbitals
!      
!     Allocate space for Fock diagonal and coefficients.
!
      call work_allocator(Fock_diagonal,n_orbitals)
!      
      call work_allocator(orb_coefficients,n_lambda)
!
!     Read in Fock diagonal and coefficients
!
      read(lusifc) (Fock_diagonal(i),i=1,n_orbitals)
      read(lusifc) (orb_coefficients(i),i=1,n_lambda)
!
!     Done with file
!
      call gpclose(lusifc,'KEEP')
!
   end subroutine sirius_reader
!
!
   subroutine mlcc3_lambda_ao
!
!  lambda and AO routine
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: read in one electron integrals and calculate lambda matrices
!
      implicit none
!
      integer  :: work_remains, i
!      
      real(dp), dimension(:), pointer :: work_point_end

!
!     Allocate space for Fock matrices, standard and t1-transformed
!     and read in one electron integrals
!
      call work_allocator(mo_fock_mat,n_basis_2)
      call work_allocator(mo_fock_mat_t1,n_basis_2)
      if(resp_option) then !C1 density
         call work_allocator(mo_fock_mat_c1,n_basis_2)
         call dzero(mo_fock_mat_c1,n_basis_2)
      end if
!      
      work_point_end => work_end()
      work_remains = work_free()
!      
      call ccrhs_oneao(mo_fock_mat,work_point_end,work_remains)
!
      call dcopy(n_basis_2,mo_fock_mat,1,mo_fock_mat_t1,1)
!
      if(print_mlcc3 .ge. 10) then
         call around('AO one electron integrals')
         call output(mo_fock_mat,1,n_basis,1,n_basis,n_basis,n_basis,1,lupri)
      end if
      if(print_mlcc3 .ge. 10) then
         call around('AO t1 one electron integrals')
         call output(mo_fock_mat_t1,1,n_basis,1,n_basis,n_basis,n_basis,1,lupri)
      end if
!      
!     Calculate Lambda matrices
!
      call work_allocator(lambda_hole,n_bas_orb)
      call work_allocator(lambda_part,n_bas_orb)
      call mlcc3_lambda
!     
      if(resp_option) then !response lambda matrices
!     
         call work_allocator(lambda_hole_resp,n_bas_orb)
         call work_allocator(lambda_part_resp,n_bas_orb)
         call mlcc3_lambda_resp
!     
      end if
!     
!     Calculate the ao densities
!
      call work_allocator(ao_density,n_basis_2)
      call work_allocator(ao_density_t1,n_basis_2)
      if(resp_option) then !C1 density
         call work_allocator(ao_density_c1,n_basis_2)
      end if
!
      call dgemm('N','T',n_basis,n_basis,n_occ,one,orb_coefficients,n_basis,orb_coefficients,&
     &            n_basis,zero,ao_density,n_basis)
!  
      call dgemm('N','T',n_basis,n_basis,n_occ,one,lambda_part,n_basis,lambda_hole,&
     &            n_basis,zero,ao_density_t1,n_basis)
!  
      if(resp_option) then !C1 density
         call dgemm('N','T',n_basis,n_basis,n_occ,one,lambda_part,n_basis,lambda_hole_resp,&
        &            n_basis,zero,ao_density_c1,n_basis)
!  
      end if
!
!
      if(print_mlcc3 .ge. 10) then
         call around('orbital coefficients')
         call output(orb_coefficients,1,n_basis,1,n_basis,n_basis,n_basis,1,lupri)
      end if
      if(print_mlcc3 .ge. 10) then
         call around('AO density')
         call output(ao_density,1,n_basis,1,n_basis,n_basis,n_basis,1,lupri)
      end if
      if(print_mlcc3 .ge. 10) then
         call around('AO density transformed')
         call output(ao_density_t1,1,n_basis,1,n_basis,n_basis,n_basis,1,lupri)
      end if
!
!
   end subroutine  mlcc3_lambda_ao
!   
!
   subroutine sirius_deallocator
!
!     Sirius deallocator
!     Authors Henrik Koch and Rolf H. Myhre
!     December 2014
!
!     Purpose: deallocate arrays allocated in sirius_reader. Easier to doit in correct order
!
      implicit none
!
      if(resp_option) then
         call work_deallocator(ao_density_c1)
      end if
      call work_deallocator(ao_density_t1)
      call work_deallocator(ao_density)
      if(resp_option) then
         call work_deallocator(lambda_part_resp)
         call work_deallocator(lambda_hole_resp)
      end if
      call work_deallocator(lambda_part)
      call work_deallocator(lambda_hole)
      if(resp_option) then
         call work_deallocator(mo_fock_mat_c1)
      end if
      call work_deallocator(mo_fock_mat_t1)
      call work_deallocator(mo_fock_mat)
      call work_deallocator(t2am)
      call work_deallocator(t1am)
      call work_deallocator(orb_coefficients)
      call work_deallocator(Fock_diagonal)
!      
   end subroutine sirius_deallocator      
!
!
   subroutine sirius_pt_deallocator
!
!     Sirius deallocator
!     Authors Henrik Koch and Rolf H. Myhre
!     June 2015
!
!     Purpose: deallocate arrays allocated in sirius_reader. MLCCD(T)
!
      implicit none
!
      call work_deallocator(t2am)
      call work_deallocator(t1am)
      call work_deallocator(orb_coefficients)
      call work_deallocator(Fock_diagonal)
!      
   end subroutine sirius_pt_deallocator      
!
!
   subroutine omega_setup(omega1_in,omega2_in,c1_in,c2_in)
!
!     Omega pointer setup
!     Authors Henrik Koch and Rolf H. Myhre
!     December 2014
!
!     Purpose: set up pointers to omega vectors
!
      implicit none
!      
      integer  :: iopt
!
      real(dp), dimension(n_t1am), intent(in), target       :: omega1_in
      real(dp), dimension(n_t2am_pack), intent(in), target  :: omega2_in
!      
      real(dp), dimension(n_t1am), intent(in), target       :: c1_in
      real(dp), dimension(n_t2am_pack), intent(in), target  :: c2_in
!      
!
!     Allocate amplitude arrays and read in
      call work_allocator(t1am,n_t1am)
      call work_allocator(t2am,n_t2am_pack)
!   
      iopt = 3
      call cc_rdrsp('R0',0,1,iopt,model,t1am,t2am)
!   
!     Initialise omega pointers
      omega1 => omega1_in(1:n_t1am)
      omega2 => omega2_in(1:n_t2am_pack)
!      
      if(resp_option) then
!
!        Initialise trial vector pointers
         c1am => c1_in(1:n_t1am)
         c2am => c2_in(1:n_t2am_pack)
!
      end if
!      
!      
   end subroutine omega_setup
!
!
   subroutine amplitude_setup
!
!     Amplitude set up
!     Authors Henrik Koch and Rolf H. Myhre
!     December 2014
!
!     Purpose: set up pointers to and read in X-amplitudes
!
      implicit none
!      
      integer  :: iopt
!
!
!     Allocate amplitude arrays and read in
      call work_allocator(t1am,n_t1am)
      call work_allocator(t2am,n_t2am_pack)
!   
      iopt = 3
      call cc_rdrsp('R0',0,1,iopt,model,t1am,t2am)
!   
   end subroutine amplitude_setup
!
!
   subroutine mlcc3_lambda
!
!     Lambda matrix calculator
!     Authors Henrik Koch and Rolf H. Myhre
!     December 2014
!
!     Purpose: calculate the lambda matrices
!
      implicit none
!      
      integer  :: vir_off
!
      vir_off = n_basis*n_occ + 1
!
      call dcopy(n_lambda,orb_coefficients,1,lambda_hole,1)
      call dcopy(n_lambda,orb_coefficients,1,lambda_part,1)
!      
      call dgemm('N','N',n_basis,n_occ,n_vir,one, &
     &           orb_coefficients(vir_off),n_basis,t1am,n_vir,one, &
     &           lambda_hole,n_basis)
!
      call dgemm('N','T',n_basis,n_vir,n_occ,-one, &
     &           orb_coefficients,n_basis,t1am,n_vir,one, &
     &           lambda_part(vir_off),n_basis)
!
!
!
   end subroutine mlcc3_lambda
!
   subroutine mlcc3_lambda_resp
!
!     Response lambda matrix calculator
!     Authors Henrik Koch and Rolf H. Myhre
!     March 2015
!
!     Purpose: calculate the lambda matrices used for response integrals
!
      implicit none
!      
      integer  :: vir_off
!
      vir_off = n_basis*n_occ + 1
!
      call dzero(lambda_hole_resp,n_lambda)
      call dzero(lambda_part_resp,n_lambda)
!
      call dgemm('N','N',n_basis,n_occ,n_vir,one, &
     &           lambda_hole(vir_off),n_basis,c1am,n_vir,zero, &
     &           lambda_hole_resp,n_basis)
!
      call dgemm('N','T',n_basis,n_vir,n_occ,-one, &
     &           lambda_part,n_basis,c1am,n_vir,zero, &
     &           lambda_part_resp(vir_off),n_basis)
!
   end subroutine mlcc3_lambda_resp
!
end module mlcc3_init
