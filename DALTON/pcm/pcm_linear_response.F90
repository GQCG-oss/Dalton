module pcm_linear_response
   
   use iso_c_binding
   use pcmmod_cfg
   use pcm_write
   use pcm_utils
   use pcm_integrals

   implicit none

   public pcm_lin_res
   public hello_world

   private

   integer, save :: lr_solver_iteration_counter = 1
   
   contains

      subroutine hello_world

      print *, "HELLO WORLD, iteration", lr_solver_iteration_counter

      lr_solver_iteration_counter = lr_solver_iteration_counter + 1

      end subroutine hello_world

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
                                                                                           
    !  call pcm_v1q0                                                                     
                                                                                           
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

end module pcm_linear_response
