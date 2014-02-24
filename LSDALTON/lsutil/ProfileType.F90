!> @file 
!> Contains 
module profile_type

  type profileinput
     logical :: doProf
     logical :: Coulomb
     logical :: Coulombecont
     logical :: Exchange
     logical :: ExchangemanyD
     logical :: Exchangegrad
     logical :: Exchangeecont
     logical :: XC
     logical :: XCFGRAD
     logical :: XCLINRSP
     logical :: XCENERGY
     logical :: Fock
     logical :: overlap
     logical :: NEGRAD
     logical :: IchorDEC
     logical :: Ichor
     logical :: IchorProfInputBasis
     CHARACTER(len=20) :: IchorProfInputBasisString(4)
     logical :: IchorProfDoThermite
     logical :: IchorProfDoIchor
  end type profileinput
  
CONTAINS
  subroutine prof_set_default_config(profinput)
    implicit none
    type(profileinput) :: profinput
    profinput%doProf = .FALSE.
    profinput%Coulomb = .FALSE.
    profinput%CoulombEcont = .FALSE.
    profinput%Exchange = .FALSE.
    profinput%ExchangeManyD = .FALSE.
    profinput%Exchangegrad = .FALSE.
    profinput%Exchangeecont = .FALSE.
    profinput%XC = .FALSE.
    profinput%XCenergy = .FALSE.
    profinput%XCfgrad = .FALSE.
    profinput%XClinrsp = .FALSE.
    profinput%Fock = .FALSE.
    profinput%Overlap = .FALSE.
    profinput%NEGRAD = .FALSE.

    profinput%Ichor = .FALSE.
    profinput%IchorDEC = .FALSE.
    profinput%IchorProfInputBasis = .FALSE.
    profinput%IchorProfInputBasisString(1) = '                    '
    profinput%IchorProfInputBasisString(2) = '                    '
    profinput%IchorProfInputBasisString(3) = '                    '
    profinput%IchorProfInputBasisString(4) = '                    '
    profinput%IchorProfDoThermite = .FALSE.
    profinput%IchorProfDoIchor = .FALSE.

  end subroutine prof_set_default_config
  
end module profile_type
