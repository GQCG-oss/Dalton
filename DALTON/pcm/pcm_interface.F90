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

   public pcm_interface_initialize
   public pcm_interface_finalize

   public collect_nctot
   public collect_atoms
   public energy_pcm_drv
   public oper_ao_pcm_drv
!   public pcm_energy_driver   
!   public pcm_fock_driver
!   public pcm_lintra_driver
!   public get_pcm_energy

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   real(c_double), allocatable :: charges(:)
   real(c_double), allocatable :: potentials(:)
   real(c_double), allocatable :: tess_cent(:, :)
   real(c_double)              :: pcm_energy
   integer(c_int)              :: nr_points = -1

   contains 
      
      subroutine pcm_interface_initialize()                              
      
      call init_pcm
      call print_pcm

      call get_cavity_size(nr_points)

      allocate(tess_cent(nr_points, 3))
      call get_tess_centers(tess_cent)

      allocate(potentials(nr_points))
      allocate(charges(nr_points))

      pcm_energy = 0.0d0
              
      is_initialized = .true.
                                                                 
      end subroutine
                                                                    
      subroutine pcm_interface_finalize()
                                                                    
      is_initialized = .false.
                                                                    
      end subroutine
                                                                    
      subroutine check_if_interface_is_initialized()

      if (.not. is_initialized) then
         print *, 'error: you try to access pcm_interface'
         print *, '       but this interface is not initialized'
         stop 1
      end if

      end subroutine
                                                                    
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
      
      end subroutine collect_atoms

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
      real(8), allocatable :: nuc_pot(:), nuc_pol_chg(:)
      real(8), allocatable :: ele_pot(:), ele_pol_chg(:)
      integer              :: kda, kdb, kpot, kcent, kfree, lwork
      character(7)         :: potName, chgName, chgName1, chgName2
      real(8)              :: factor
                                                                                    
      allocate(nuc_pot(nr_points))
      allocate(nuc_pol_chg(nr_points))
      allocate(ele_pot(nr_points))
      allocate(ele_pol_chg(nr_points))

      kda   = 1
      kdb   = kda   + nnbasx
      kfree = kdb + nnbasx
      lwork = lfree - kfree + 1
      if (lwork .lt. 0) then 
              call errwrk('energy_pcm_drv', kfree, lfree) 
      end if
                                                                                    
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
      call nuc_pot_pcm(nr_points, tess_cent, nuc_pot)
      call set_surface_function(nr_points, nuc_pot, potName)
      call comp_chg_pcm(potName, chgName)
      call get_surface_function(nr_points, nuc_pol_chg, chgName)
                                                                                    
      potName = 'ElePot'//CHAR(0)
      chgName = 'EleChg'//CHAR(0)
      call ele_pot_pcm(nr_points, tess_cent, ele_pot, work(kda), work(kfree), lfree)
      call set_surface_function(nr_points, ele_pot, potName)
      call comp_chg_pcm(potName, chgName)
      call get_surface_function(nr_points, ele_pol_chg, chgName)

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
      call pcm_write_file_separate(nr_points, nuc_pot, nuc_pol_chg, ele_pot, ele_pol_chg)

      end subroutine energy_pcm_drv
      
      subroutine oper_ao_pcm_drv(oper, charge, work, lwork)
!
! Calculate exp values of potentials on tesserae
! Input: symmetry packed Density matrix in AO basis
!        cavity points
! Output: expectation values of electrostatic potential on tesserae
!

      real(8), intent(out) :: oper(*)
      real(8)              :: work(*)
      character            :: charge(*)
      integer              :: lwork

      integer              :: kcharge, kcenters, kfree, lfree

       
      kcharge = 1
      kfree = kcharge + nr_points
      lfree = lwork - kfree
                                                                                  
      if(lfree .le. 0) then 
              call quit('Not enough mem in oper_ao_pcm_drv')
      end if
        
      call get_surface_function(nr_points, work(kcharge), charge)
      call j1int_pcm(work(kcharge), nr_points, tess_cent, .false.,    & 
                     oper, 1, .false., 'NPETES ', 1, work(kfree), lfree)

      end subroutine oper_ao_pcm_drv

end module
