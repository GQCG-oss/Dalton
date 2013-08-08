module lsdalton_fock_module
use typedeftype, only: lsitem
use matrix_module
private
public :: lsint_fock_data_type, lsint_fock_data
Type lsint_fock_data_type
 type(lsitem),pointer  :: ls
 type(Matrix),pointer  :: H1
 integer               :: lupri, luerr
end Type lsint_fock_data_type

type(lsint_fock_data_type),save  ::  lsint_fock_data

end module lsdalton_fock_module

