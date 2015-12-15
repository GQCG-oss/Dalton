module mlcc3_data
!
!
!  mlcc3 types
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: Store variables required througout MLCC3
!
   use mlcc_typedef
!   
   implicit none
!
!
!  Model variable from old Dalton
   character(10)                    :: model ='CCSD      '
!   
!  Energy or response vector
   logical                          :: resp_option = .false.
   real(dp)                         :: freq
!   
!  Integer variables
   integer                          :: n_orbitals, n_basis, n_lambda, n_occ, n_vir
   integer                          :: n_t1am, n_t2am, n_t2am_pack
   integer                          :: n_v_2, n_v_3, n_basis_2, n_basis_2_pack, n_bas_orb
   integer                          :: n_ao_ints
!
!  Info read from file  
   real(dp)                         :: nuclear_potential,scf_energy
!   
   integer                          :: print_mlcc3 = 1
!   
!  Basic variable
   logical                          :: mlcc3_active
   logical                          :: mlcc3_nrg_spa
   logical                          :: mlcc3_nrg_gen
   integer                          :: n_active, n_general
   integer                          :: n_occ_inp,n_vir_inp
   integer                          :: n_gen_o_inp,n_gen_v_inp
   logical                          :: set = .false.
!
!  Various pointers
!
!  Packed Omega vectors from CCSD
   real(dp), dimension(:), pointer  :: omega1            => null()
   real(dp), dimension(:), pointer  :: omega2            => null()
!   
!  T1 and T2 amplitudes, T2 packed
   real(dp), dimension(:), pointer  :: t1am              => null()
   real(dp), dimension(:), pointer  :: t2am              => null()
!
!  T1 and T2 amplitudes, T2 packed
   real(dp), dimension(:), pointer  :: c1am              => null()
   real(dp), dimension(:), pointer  :: c2am              => null()
!
!  HF orbital coefficients and fock diagonal elements
   real(dp), dimension(:), pointer  :: orb_coefficients  => null()
   real(dp), dimension(:), pointer  :: Fock_diagonal     => null()
!   
!  The MO Fock matrix, standard and T1- and C1-transformed
   real(dp), dimension(:), pointer  :: mo_fock_mat       => null()
   real(dp), dimension(:), pointer  :: mo_fock_mat_t1    => null()
   real(dp), dimension(:), pointer  :: mo_fock_mat_c1    => null()
!   
!  Lambda matrices
   real(dp), dimension(:), pointer  :: lambda_hole       => null()
   real(dp), dimension(:), pointer  :: lambda_part       => null()
!   
   real(dp), dimension(:), pointer  :: lambda_hole_resp  => null()
   real(dp), dimension(:), pointer  :: lambda_part_resp  => null()
!   
!  AO density matrices, standard and T1- and C1-transformed
   real(dp), dimension(:), pointer  :: ao_density        => null()
   real(dp), dimension(:), pointer  :: ao_density_t1     => null()
   real(dp), dimension(:), pointer  :: ao_density_c1     => null()
!   
!  File names
!
!  MO integrals
   character(len=12)                :: bDck_file_name = "bDck_mo_ints"
   character(len=12)                :: Dbkc_file_name = 'Dbkc_mo_ints'
   character(len=12)                :: Ljck_file_name = 'Ljck_mo_ints'
   character(len=12)                :: jLkc_file_name = 'jLkc_mo_ints'
   character(len=12)                :: jbkc_file_name = 'jbkc_mo_ints'
!
!  response MO integrals
   character(len=12)                :: bDck_resp_name = "bDck_mo_resp"
   character(len=12)                :: Dbkc_resp_name = 'Dbkc_mo_resp'
   character(len=12)                :: Ljck_resp_name = 'Ljck_mo_resp'
   character(len=12)                :: jLkc_resp_name = 'jLkc_mo_resp'
!
!
end module mlcc3_data
