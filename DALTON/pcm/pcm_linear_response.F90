module pcm_linear_response
   
use, intrinsic :: iso_c_binding
use pcmmod_cfg
use pcm_write
use pcm_utils
use pcm_integrals

implicit none

public pcm_lr_initialize
public pcm_lr_finalize
public pcm_lin_res
public hello_world

private
! If false the interface will refuse to be accessed
logical                     :: is_initialized = .false.
real(c_double), allocatable :: tess_cent(:, :)
integer(c_int)              :: nr_points
integer(c_int)              :: nr_points_irr
! A (maybe clumsy) way of passing LUPRI without using common blocks   
integer                     :: global_print_unit
! A counter for the number of LR iterations
integer, save               :: LR_iteration_counter = 1

contains

   subroutine hello_world
   
   print *, "HELLO WORLD, iteration", LR_iteration_counter
           
   LR_iteration_counter = LR_iteration_counter + 1

   end subroutine hello_world
      
   subroutine pcm_lr_initialize(print_unit)
   
   use pcmmod_cfg, only: pcmmod_host_provides_input

   integer, intent(in) :: print_unit
   
   integer :: host_provides_input = 0

   global_print_unit = print_unit

   if (pcmmod_host_provides_input) host_provides_input = 1 

   call set_up_pcm(host_provides_input)
   call print_pcm

   call get_cavity_size(nr_points, nr_points_irr)

   allocate(tess_cent(3, nr_points))
   tess_cent = 0.0d0
   call get_tesserae(tess_cent)
           
   is_initialized = .true.

   end subroutine pcm_lr_initialize

   subroutine pcm_lr_finalize()

   if (.not. is_initialized) then
      print *, 'Error: pcm_linear_response was never initialized.'
      stop 1
   end if
   ! Free the memory taken from the free store both in Fortran and in C++
   deallocate(tess_cent)

   call tear_down_pcm
                                                                 
   is_initialized = .false.
                                                                 
   end subroutine pcm_lr_finalize

   subroutine check_if_interface_is_initialized()

   if (.not. is_initialized) then
      print *, 'Error: pcm_linear_response is not initialized.'
      stop 1
   end if

   end subroutine

   subroutine pcm_lin_res(trial_vec, density, cmo, work, lwork)                         
   
   real(8) :: trial_vec(*)
   real(8) :: density(*)
   real(8) :: cmo(*)
   real(8) :: work(*)
   integer :: lwork

   character(7) :: potName, chgName
   integer      :: kfree, lfree, irrep
                                                                                        
                                                                                        
   kfree = 1                                                                            
   lfree = lwork - kfree + 1                                                            
   if (lfree .lt. 0) then                                                               
           call errwrk('ieflno', kfree, lwork)                                          
   end if                                                                               
                                                                                        
   call pcm_v1(trial_vec, density, cmo, work(kfree), lfree)                             
   potName = 'PotOIT'//CHAR(0)                                                          
   chgName = 'ChgOIT'//CHAR(0)        
   ! This will have to be modified. When doing linear response we might be
   ! interested in some other irreducible representation
   irrep = 0
   call compute_asc(potName, chgName, irrep)
                                                                                        
   !call pcm_v1q0                                                                     
                                                                                        
   end subroutine pcm_lin_res                                                           

   subroutine pcm_v0q1(rsp_grad, rsp_dim, trial_vec, density, cmo, work, lwork)         

#include "priunit.h"                                                                 
#include "inforb.h"                                                                  
                                                                                           
 !                                                                                         
 ! 1) No symmetry yet (see rspief2 before impl symmetry)                                   
 !                                                                                         
                                                                                         
   integer :: rsp_dim, lwork                                                                 
   real(8) :: tessera(3), rsp_grad(rsp_dim, *)                                               
   real(8) :: work(*), density(*)                                                            
   real(8) :: trial_vec(*), cmo(*)                                                           
   real(8) :: expval                                                                         
                                                                                             
   character(7) ::  potname                                                                  
                                                                                             
   integer :: nosim = 1                                                                      
   integer :: isymtv = 1 ! NO SYMMETRY                                                       
   integer :: isympot = 1 ! NO SYMMETRY                                                      
   integer :: iprlocal = 0                                                                   
   logical :: trimat = .true.                                                                
   integer :: kvoit, kpotao, kpotmo, kpotmo2 
   integer :: ktemp, kucmo, kubo, kfree, lfree, nts, i
                                                                                            
   call get_cavity_size(nts)                                                                 
   kvoit = 1                                                                                 
   kpotao  = kvoit   + nts                                                                   
   kpotmo  = kpotao  + nnbasx * nsym                                                         
   kpotmo2 = kpotmo  + n2orbx                                                                
   ktemp   = kpotmo2 + n2orbx                                                                
   kucmo   = ktemp   + n2basx                                                                
   kubo    = kucmo   + norbt * nbast                                                         
   kfree   = kubo    + n2orbx * nosim                                                        
   lfree   = lwork   - kfree + 1                                                             
                                                                                             
   call upkcmo(cmo, work(kucmo))                                                             
   call rspzym(nosim, trial_vec, work(kubo))                                                 
   call dscal(nosim*n2orbx, -1.0d0, work(kubo), 1)                                           
                                                                                             
   if (kfree .gt. lwork) then 
           call errwrk('pcm_trans_pot', kfree, lwork)                          
   end if
   do i = 1, nts
! how many of these really have to be set to zero??                                  
      call dzero(work(kpotao), nnbasx)                                                  
      call dzero(work(kpotmo), n2orbx)                                                  
      call dzero(work(kpotmo2), n2orbx)                                                 
      call get_tesserae_centers(i, tessera)                                              
      call pot_int_tess(work(kpotao), tessera, trimat, work(kfree), lfree)              
      call uthu(work(kpotao), work(kpotmo), work(kucmo), work(kfree), nbast, norbt)     
      call dsptsi(norbt, work(kpotmo), work(kpotmo2))                                   
      call dzero(work(kpotmo), n2orbx)                                                  
      call oith1(isymtv, work(kubo), work(kpotmo2), work(kpotmo), isympot)              
      call melone(work(kpotmo), isympot, density, 1.0d0, expval, iprlocal, 'PotOIT')    
      work(kvoit + i - 1) = expval                                                      
      print *, 'Tr pot', i, expval                                                      
   end do                                                                               
   potname = "PotOIT"//CHAR(0)                                                          
   call set_surface_function(nts, work(kvoit), potname)                                 
   call print_surface_function(potname)                                                 
   
   end subroutine pcm_v0q1                                                              
                                                                                           
      !subroutine pcm_v1q0(rsp_grad, trial_vec, density, cmo, work, lwork)                 
      !                                                                                    
      !end subroutine pcm_v1q0                                                             
                                                                                           
                                                                                           
   subroutine pcm_v1(trial_vec, density, cmo, work, lwork)                             

#include "priunit.h"                                                                 
#include "inforb.h"                                                                  
                                                                                           
!                                                                                    
! 1) No symmetry yet (see rspief2 before impl symmetry)                              
!                                                                                    

   real(8) :: trial_vec(*), density(*), cmo(*), work(*)
   integer :: lwork

   real(8) :: tessera(3)
   real(8) :: expval
   logical :: trimat
   character(7) :: potname
   integer :: nosim, isymtv, isympot, iprlocal
   integer :: kvoit, kpotao, kpotmo, kpotmo2, ktemp, kucmo
   integer :: kubo, kfree, lfree, i, nts
                                                                                        
   nosim = 1 ! NO SYMMETRY                                                              
   isymtv = 1 ! NO SYMMETRY                                                             
   isympot = 1 ! NO SYMMETRY                                                            
   iprlocal = 0                                                                         
   trimat = .true.                                                                      
   call get_cavity_size(nts)                                                            
   kvoit   = 1                                                                          
   kpotao  = kvoit   + nts                                                              
   kpotmo  = kpotao  + nnbasx * nsym                                                    
   kpotmo2 = kpotmo  + n2orbx                                                           
   ktemp   = kpotmo2 + n2orbx                                                           
   kucmo   = ktemp   + n2basx                                                           
   kubo    = kucmo   + norbt * nbast                                                    
   kfree   = kubo    + n2orbx * nosim                                                   
   lfree   = lwork   - kfree + 1                                                        
                                                                                        
   call upkcmo(cmo, work(kucmo))                                                        
   call rspzym(nosim, trial_vec, work(kubo))                                            
   call dscal(nosim*n2orbx, -1.0d0, work(kubo), 1)                                      
                                                                                        
   if (kfree .gt. lwork) then 
           call errwrk('pcm_trans_pot', kfree, lwork)                     
   end if
   do i = 1, nts                                                                      
! how many of these really have to be set to zero??                                  
      call dzero(work(kpotao), nnbasx)                                                  
      call dzero(work(kpotmo), n2orbx)                                                  
      call dzero(work(kpotmo2), n2orbx)                                                 
      call get_tesserae_centers(i, tessera)                                              
      call pot_int_tess(work(kpotao), tessera, trimat, work(kfree), lfree)              
      call uthu(work(kpotao), work(kpotmo), work(kucmo), work(kfree), nbast, norbt)     
      call dsptsi(norbt, work(kpotmo), work(kpotmo2))                                   
      call dzero(work(kpotmo), n2orbx)                                                  
      call oith1(isymtv, work(kubo), work(kpotmo2), work(kpotmo), isympot)              
      call melone(work(kpotmo), isympot, density, 1.0d0, expval, iprlocal, 'PotOIT')    
      work(kvoit + i - 1) = expval                                                      
      print *, 'Tr pot', i, expval                                                      
   end do                                                                               
   
   potname = "PotOIT"//CHAR(0)                                                          
   call set_surface_function(nts, work(kvoit), potname)                                 
   call print_surface_function(potname)                                                 
   
   end subroutine pcm_v1                                                                
                                                                                           
!      subroutine pcm_pot_int(work, lwork)                                           
!#include <implicit.h>                                                               
!#include <priunit.h>                                                                
!#include <inforb.h>                                                                 
!                                                                                    
!      logical trimat                                                                
!      double precision tessera(3), work(*)                                          
!                                                                                    
!      trimat = .true.                                                               
!                                                                                    
!      call get_cavity_size(nts)                                                     
!      kcent = 1                                                                     
!      kpotao = kcent + 3 * nts                                                      
!      kpotmo = kpotao + nnbasx * nsym                                               
!      ktemp  = kpotmo + nnorbx                                                      
!      kucmo  = ktemp  + n2basx                                                      
!      kfree  = kucmo  + norbt * nbast                                               
!      lfree = lwork - kfree + 1                                                     
!      call dzero(work(kpotint), nnbasx * nsym)                                      
!      if (lfree .lt. 0) call errwrk('ieflno', kfree, lwork)                         
!      do i = 1, nts                                                                 
!         call get_tesserae_centers(i, tessera)                                       
!         call pot_int_tess(work(kpotao), tessera, trimat, work(kfree),              
!     &        lfree)                                                                
!      end do                                                                        
!      end subroutine pcm_pot_int

      subroutine pcm_oit_potential(vse1, nosim, udv, udvtr, jwopsy, cmo, bovecs, work, lwork)

      use pcm_integrals, only: j1x_pcm

#include "inforb.h"
#include "priunit.h"
#include "infrsp.h"
#include "maxorb.h"

      ! Parameters
      ! Passed variables
      real(8) :: vse1(*)
      integer :: nosim
      real(8) :: udv(*)
      real(8) :: udvtr(*) 
      integer :: jwopsy
      real(8) :: cmo(*)
      real(8) :: bovecs(*)
      real(8) :: work(*)
      integer :: lwork
      ! Local variables
      logical :: tofile
      integer :: iosim
      integer :: jubo, kfree, kubo, kucmo, lfree

      kubo = 1
      kucmo = kubo + nosim * n2orbx
      kfree = kucmo + norbt * nbast
      lfree = lwork - kfree + 1
      if (lfree .lt. 0) then 
         call quit('pcm_oit_potential: insufficient memory')
      end if
      call upkcmo(cmo,work(kucmo))      

      if (nosim.gt.0) then
         call rspzym(nosim,bovecs,work(kubo))
         call dscal(nosim*n2orbx,-1.0d0,work(kubo),1)
         if (iprrsp .ge. 55) then
            do iosim = 1, nosim
               jubo = kubo + (iosim - 1) * n2orbx
               write(lupri, '(a, i3, i3)') 'Orbital trial vector unpacked to matrix form (no.', iosim,' of', nosim,')'
               call output(work(jubo), 1, norbt, 1, norbt, norbt, norbt, 1, lupri)
            end do 
         end if
      end if

      call j1x_pcm(nr_points, nr_points_irr, tess_cent, nosim, vse1, work(kucmo), &
                   work(kubo), udv, udvtr, tofile, jwopsy, work(kfree), lfree)

      end subroutine pcm_oit_potential

end module pcm_linear_response
