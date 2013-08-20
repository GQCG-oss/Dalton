!
! DALTON-side interface routines for the Polarizable Continuum Model
! Luca Frediani, Roberto Di Remigio 2011-2013
!
! We divide the interface routines into:
!   1. cavity-related routines;
!   2. solver-related routines.
!

!
!                   CAVITY-RELATED ROUTINES
!
! We shall provide the following data to the module:
!   1. the number of nuclei in the molecule. NUCLEI variable;
!   2. the atomic number of each nucleus. CHARGES vector;
!   3. the coordinates of each nucleus. CENTERS matrix.
! PCMSolver will provide to cavity formation depending
! on the selected mode:
!   1. Mode = Implicit means that the spheres are
!      centered according to the CENTERS matrix with
!      radii taken from the internal library;
!   2. Mode = Atoms means that we specify sphere centers
!      and radii for certain atoms in the molecule. All the
!      other spheres are created as in the Implicit mode.
!


module pcm_interface
  
   use iso_c_binding
   use pcm_integrals
   use pcm_write
   use pcmmod_cfg
   use pcm_utils

   implicit none 
 
   public collect_nctot
   public collect_atoms
   public energy_pcm_drv
   public oper_ao_pcm_drv
!   public pcm_energy_driver   
!   public pcm_fock_driver
!   public pcm_lintra_driver
!   public get_pcm_energy

   private

   real(8), allocatable :: charges(:)
   real(8), allocatable :: potentials(:)
   real(8), allocatable :: tess_cent(:, :)
   real(8)              :: pcm_energy
! nr_points is 32 bit to avoid mess in passing it from C++ to FORTRAN
   integer(4)           :: nr_points = -1

   contains 
      
      subroutine collect_nctot(nr_nuclei) bind(c, name='collect_nctot_')

#include "mxcent.h"
#include "nuclei.h"

      integer(c_int), intent(out) :: nr_nuclei
      
      nr_nuclei = nucdep

      end subroutine collect_nctot
      
      subroutine collect_atoms(charges, centers) bind(c, name='collect_atoms_')

#include "mxcent.h"
#include "nuclei.h"
      
      real(c_double), intent(out) :: charges(*)
      real(c_double), intent(out) :: centers(3,*)
      
      integer :: i, j, k 

! Get coordinates
      call getacord(centers)
! Get charges      
      i = 0
      do j = 1, nucind
         do k = 1, nucdeg(j)
            i = i + 1
            charges(i) = charge(j)
         enddo
      enddo
      
      end subroutine

      subroutine energy_pcm_drv(dcao, dvao, pol_ene, work, lfree)

#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "inftap.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "dftcom.h"

        real(8)      :: dcao(*)
        real(8)      :: dvao(*)                                                       
        real(8)      :: pol_ene
        real(8)      :: work(*)
        integer      :: lfree
! Local variables
        character(7) :: potName, chgName, chgName1, chgName2
        integer(4)   :: nts 
        integer      :: kda, kdb, kpot, kcent, kfree, lwork
        real(8), allocatable :: nuc_pot(:), nuc_pol_chg(:), ele_pot(:), ele_pol_chg(:)
        real(8)      :: factor
                                                                                      
        call get_cavity_size(nts)
        allocate(nuc_pot(nts))
        allocate(nuc_pol_chg(nts))
        allocate(ele_pot(nts))
        allocate(ele_pol_chg(nts))
        kda   = 1
        kdb   = kda   + nnbasx
        kpot  = kdb   + nnbasx
        kcent = kpot  + nts
        kfree = kcent + 3 * nts
        lwork = lfree - kfree + 1
        if (lwork .lt. 0) then 
                call errwrk('energy_pcm_drv', kfree, lfree) 
        end if
                                                                                      
! 1) Get tessera data
        call get_tess_centers(work(kcent))
        
        if ((nasht.gt.0) .and. .not. dftadd) then
           call dcopy(nnbasx, dcao, 1, work(kdb), 1)
           call daxpy(nnbasx, 1.0d0, dvao, 1, work(kdb), 1)
           call pksym1(work(kda), work(kdb), nbas, nsym, -1)
        else
           call pksym1(work(kda), dcao, nbas, nsym, -1)
        end if
                                                                                      
! 2) Compute potentials
! 3) Compute charges
        potName = 'NucPot'//CHAR(0)
        chgName = 'NucChg'//CHAR(0)
        !call nuc_pot_pcm(nts, work(kcent), work(kpot))
        call nuc_pot_pcm(nts, work(kcent), nuc_pot)
        call set_surface_function(nts, nuc_pot, potName)
        call comp_chg_pcm(potName, chgName)
        call get_surface_function(nts, nuc_pol_chg, chgName)
                                                                                      
        potName = 'ElePot'//CHAR(0)
        chgName = 'EleChg'//CHAR(0)
        call ele_pot_pcm(nts, work(kcent), ele_pot, work(kda), work(kfree), lfree)
        call set_surface_function(nts, ele_pot, potName)
        call comp_chg_pcm(potName, chgName)
        call get_surface_function(nts, ele_pol_chg, chgName)

        call comp_pol_ene_pcm(pol_ene, 0)
                                                                                      
        chgName  = 'TotChg'//CHAR(0)
        chgName1 = 'NucChg'//CHAR(0)
        chgName2 = 'EleChg'//CHAR(0)
                                                                                      
        call append_surf_func(chgName)
        call clear_surf_func(chgName)
        factor =  1.0
        call add_surface_function(chgName, factor, chgName1)
        call add_surface_function(chgName, factor, chgName2)

! Write to file MEP and ASC
        call pcm_write_file_separate(nts, nuc_pot, nuc_pol_chg, ele_pot, ele_pol_chg)

      end subroutine
      
!
! Calculate exp values of potentials on tesserae
! Input: symmetry packed Density matrix in AO basis
!        cavity points
! Output: expectation values of electrostatic potential on tesserae
!
      subroutine oper_ao_pcm_drv(oper, charge, work, lwork)

        real(8), intent(out) :: oper(*)
        real(8)              :: work(*)
        character            :: charge(*)
        integer              :: lwork

        integer(4)           :: nts
        integer              :: kcharge, kcenters, kfree, lfree

       
        call get_cavity_size(nts)                                                      
        kcharge = 1
        kcenters = kcharge + nts
        kfree = kcenters + 3 * nts
        lfree = lwork - kfree
                                                                                  
        if(lfree .le. 0) then 
                call quit('Not enough mem in oper_ao_pcm_drv')
        end if
        
        call get_surface_function(nts, work(kcharge), charge)
        call get_tess_centers(work(kcenters))
        call j1int_pcm(work(kcharge), nts, work(kcenters), .false.,            & 
                       oper, 1, .false., 'NPETES ', 1, work(kfree), lfree)

      end subroutine

end module
