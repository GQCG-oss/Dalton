module trustradius_mod
   use arhDensity, only: solveritem
   use TYPEDEFTYPE, only: lsitem
   use precision
   use davidson_settings
contains
!> \brief Update of trust-radius - see Molecular Electronic Structure Theory p. 615
!> \author S. Host
!> \date 2005
!>
!> To be used for arh and 2nd order optimization only.
!>
   subroutine update_trustradius(arh, ls, SCF_it, E, ndens)
   implicit none
      !> Contains setting for ARH solver
      type(SolverItem)               :: arh
      !> Contains settings for integrals
      type(lsitem),intent(inout)     :: ls
      !> Current SCF iteration
      integer, intent(in)            :: SCF_it 
      !> Current SCF energy
      real(realk)                    :: E
      !> Number of densities/Fock/KS matrices in subspace
      integer, intent(in)            :: ndens
      real(realk)                    :: r, tr_ratio
      character*8                    :: tr_criterion

   arh%debug_arh_scfit = SCF_it
   if (arh%set_optxelm) then
      tr_criterion = ' (xmax) '
   else
      tr_criterion = ' (xnorm)' 
   endif
   
   arh%step_accepted = .true.

   !The predicted SCF energy change is
   !Q = E(SCF)_n + (E1_red)T * c + 1/2*cT * (E2_red) * c
   !and the ratio between actual and predicted energy change is
   !    E(SCF)_n+1 - E(SCF)_n
   !r = _____________________
   !        Q - E(SCF)_n
 
   !The ratio denominator
   !Q - E(SCF)_n = (E1_red)T * c + 1/2*cT * (E2_red) * c
   !has been calculated in ARH module and is found in arh%denom  

   !write(arh%lupri,*) 'Denominator is:', arh%denom
   !write(arh%lupri,*) 'Numerator is:',  E - arh%old_energy
   !Now calculate the ratio:
   r = (E - arh%old_energy)/arh%denom
   if (r < 0) then
      arh%Nrejections = arh%Nrejections + 1
      if (arh%info_lineq) write(arh%lupri,*) 'Number of steps rejected:', arh%Nrejections
      arh%step_accepted = .false.
   endif
   !write(lupri,*) '% Old and new SCF energy:', ESCF_old, ESCF_new
   !write(lupri,*) '% SCF energy change:', ESCF_new - ESCF_old 
   !write(lupri,*) '% Predicted energy change:', denom 
   !write(lupri,*) '% Ratio between observed and predicted changes:', r 

  !if (.not. cfg_do_2nd_order) then
      if (arh%set_optxelm) then
      write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F14.4, F10.4, i10, i11, '    %%%')") &
           & arh%set_max_element, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, arh%D_para_D_tot, ndens, SCF_it-1
      else
      write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F14.4, F10.4, i10, i11, '    %%%')") &
           & arh%set_max_step, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, &
           & arh%D_para_D_tot, ndens, SCF_it-1
      endif
   !else
   !   if (arh%set_optxelm) then
   !   write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F10.4, F10.4, D14.4, D13.4, i5, '    %%%')") &
   !        & arh%set_max_element, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, &
   !        & arh%D_para_D_tot, arh%gradcontrib, arh%hescontrib, SCF_it-1
   !   else
   !   write (arh%lupri, "(F8.5, A8, F11.5, F11.5, F8.2, F10.4, F10.4, D14.4, D13.4, i5, '    %%%')") &
   !        & arh%set_max_step, tr_criterion, arh%maxelm, arh%xnorm, -arh%current_mu, r, &
   !        & arh%D_para_D_tot, arh%gradcontrib, arh%hescontrib, SCF_it-1
   !   endif
   !endif
 
!16/9 - 2010: Brano uses the solver for localization of orbital, and sometimes
!the ratio gets really large, which seems to cause problems. We therefore
!introduce:
! r > 0.95 double TR -> changed to 1.5E0_realk > r > 0.95
! 1.5 < r < 2.0 leave TR unchanged
! r > 2.0 contract TR
! /Stinne
!Removed, didn't work well! /Stinne 27/10-10

!if (.not. arh%lshift_by_hlgap) then
   if (arh%set_optxelm) then
      if (arh%info_lineq) write(arh%lupri,*) 'Trust radius based on max element!'
      tr_ratio = arh%maxelm/arh%set_max_element
      if (arh%maxelm > 5.0E-2_realk .and. arh%set_max_element > 1.0E-1_realk .and. &
          & tr_ratio < 0.9E0_realk .and. r > 0.0E0_realk) then !.and.r<1.0E0_realk.and. abs(arh%current_mu)>1.0E-2_realk)then
         if (arh%info_lineq) write(arh%lupri,*) '% Too large difference between trust-radius and actual step. &
                        & Reduce trust-radius'
         if (arh%set_max_element*arh%cfg_arh_contract > arh%maxelm) then
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', &
                & arh%set_max_element, arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_element = arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
         else
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', arh%set_max_element, arh%maxelm
            arh%set_max_element = arh%maxelm
            arh%set_max_step = arh%xnorm 
         endif
      else
      !Now take action - change the trust-radius in accordance with r:
         if (r > arh%cfg_arh_expand_crit .and. r < 1.5E0_realk) then  !we can take larger steps - expansion
            !arh%set_max_step = xnorm
            if (arh%trustradius_decreased) then
               if (arh%info_lineq) write(arh%lupri,*) "Trust-radius was decreased to obtain convergence, don't expand xmax!"
            !else if (abs(arh%current_mu) < 0.001 .and. arh%set_max_element >= cfg_trmaxelm) then
            !   write(lupri,*) "Mu ~ 0, trust region not expanded!"
            else if (r > 0.95E0_realk) then ! .and. r < 1.5E0_realk) then
               if (arh%info_lineq) write(arh%lupri,*) 'Large ratio, double trust-radius'
               arh%set_max_element = arh%set_max_element*2.0E0_realk
               arh%set_max_step = arh%set_max_step*2.0E0_realk
            !else if (tr_ratio < 0.9E0_realk) then
            !   write(lupri,*) '% Too large difference between trust-radius and actual step - &
            !            & do not expand trust-radius.'
            else
               if (arh%info_lineq) write(arh%lupri,*) '% Expand trust-radius by a factor',  arh%cfg_arh_expand, 'h(old), h(new):', &
                              & arh%set_max_element, arh%set_max_element*arh%cfg_arh_expand
               arh%set_max_element = arh%set_max_element*arh%cfg_arh_expand
               arh%set_max_step = arh%xnorm*arh%cfg_arh_expand
               if (arh%info_lineq) write(arh%lupri,*) '% Max ||X|| updated to', arh%set_max_step
               if (arh%info_lineq) write(arh%lupri,*) '% Max ||X|| updated, old, new:', arh%set_max_step, arh%xnorm !Keep it safe...
            endif
         else if (arh%cfg_arh_contract_crit < r .and. r < arh%cfg_arh_expand_crit) then ! .or. &
                !& (r > 1.5E0_realk .and. r < 2.0E0_realk)) then
            if (arh%info_lineq) write(arh%lupri,*) '% trust-radius ok, h = ', arh%set_max_element
            arh%set_max_step = arh%xnorm
            !if (r > 1.5E0_realk .and. r < 2.0E0_realk .and. arh%info_lineq) write(arh%lupri,*) '% 1.5 < r < 2.0, trust-radius unchanged'  
         else if (0.0E0_realk < r .and. r < arh%cfg_arh_contract_crit) then
         !else if ((0.0E0_realk < r .and. r < arh%cfg_arh_contract_crit) .or. &
         !        & r > 2.0E0_realk) then
            if (arh%info_lineq) write(arh%lupri,*) '% Reduce trust-radius, h(old), h(new):', arh%set_max_element, &
                 &arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_element = arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_step = arh%xnorm
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            !if (r > 2.0E0_realk .and. arh%info_lineq) write(arh%lupri,*) '% r > 2.0, trust-radius contracted'
         else if (r < 0.0E0_realk) then
            arh%set_max_element = arh%maxelm
            arh%set_max_step = arh%xnorm
            if (arh%info_lineq) WRITE(arh%lupri,*) &
                 &'Reject step and reduce trust-radius, h(old), h(new):',&
                 & arh%set_max_element, arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_element = arh%set_max_element*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            arh%step_accepted = .false.
         endif
      endif
   else
      if (arh%info_lineq) write(arh%lupri,*) 'Trust radius based on max norm!'
      tr_ratio = arh%xnorm/arh%set_max_step
      if (arh%xnorm > 5.0E-2_realk .and. arh%set_max_step > 1.0E-1_realk .and. &
       & tr_ratio < 0.9E0_realk .and. r > 0 .and. abs(arh%current_mu) > 1.0E-2_realk) then
         if (arh%info_lineq) write(arh%lupri,*) '% Too large difference between trust-radius and actual step. Reduce trust-radius'
         if (arh%set_max_step*arh%cfg_arh_contract > arh%xnorm) then
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', arh%set_max_step, arh%set_max_step*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
         else
            if (arh%info_lineq) write(arh%lupri,*) '% h(old), h(new):', arh%set_max_step, arh%xnorm
            arh%set_max_step = arh%xnorm
         endif
      else
         if (r > arh%cfg_arh_expand_crit) then ! .and. r < 1.5E0_realk) then  !we can take larger steps - expansion
            if (arh%trustradius_decreased) then
               if (arh%info_lineq) write(arh%lupri,*) "Trust-radius was decreased to obtain convergence, don't expand xnorm!"
            !else if (abs(arh%current_mu) < 0.001 .and. arh%set_max_step >= cfg_trmaxstep) then
            !   write(lupri,*) "Mu ~ 0, trust region not expanded!"
            else if (r > 0.95E0_realk) then ! .and. r < 1.5E0_realk) then
               if (arh%info_lineq) write(arh%lupri,*) 'Large ratio, double trust-radius'
               arh%set_max_step = arh%set_max_step*2.0E0_realk
            else
               if (arh%info_lineq) write(arh%lupri,*) '% Expand trust-radius by a factor',  arh%cfg_arh_expand, 'h(old), h(new):', &
                              & arh%set_max_step, arh%set_max_step*arh%cfg_arh_expand
               arh%set_max_step = arh%set_max_step*arh%cfg_arh_expand
            endif
         else if (arh%cfg_arh_contract_crit < r .and. r < arh%cfg_arh_expand_crit) then !  .or. &
               ! & (r > 1.5E0_realk .and. r < 2.0E0_realk)) then
            if (arh%info_lineq) then           
               !if (r > 1.5E0_realk .and. r < 2.0E0_realk) then
               !   write(arh%lupri,*) '% 1.5 < r < 2.0, trust-radius unchanged'
               !else
                  write(arh%lupri,*) '% trust-radius ok, h = ', arh%set_max_step
               !endif
            endif
         else if (0.0E0_realk < r .and. r < arh%cfg_arh_contract_crit) then
         !else if (0.0E0_realk < r .and. r < arh%cfg_arh_contract_crit .or. &
         !        & r > 2.0E0_realk) then
            if (arh%info_lineq) write(arh%lupri,*) '% Reduce trust-radius, h(old), h(new):', arh%set_max_step, &
                 &arh%set_max_step*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            !if (r > 2.0E0_realk .and. arh%info_lineq) write(arh%lupri,*) '% r > 2.0, trust-radius contracted'
         else if (r < 0.0E0_realk) then
            arh%set_max_step = arh%xnorm
            if (arh%info_lineq) write(arh%lupri,*) '% Reject step and reduce trust-radius, h(old), h(new):', &
                           & arh%set_max_step, arh%set_max_step*arh%cfg_arh_contract
            arh%set_max_step = arh%set_max_step*arh%cfg_arh_contract
            arh%step_accepted = .false.
         endif
      endif
   endif

   if (arh%set_max_element > arh%cfg_max_element) then
      if (arh%info_lineq) write (arh%lupri,*) 'Absolute max element allowed is', arh%cfg_max_element,', resetting trust radius!'
      arh%set_max_element = arh%cfg_max_element
   endif
   if (arh%set_max_step > arh%cfg_max_step) then
      if (arh%info_lineq) write (arh%lupri,*) 'Absolute max step allowed is', arh%cfg_max_step,' , resetting trust radius!'
      arh%set_max_step = arh%cfg_max_step
   endif

   !Removed 13/10/2010, it doesn't work properly. /Stinne
   !if (.not. arh%step_accepted .and. arh%Nrejections > 1 .and. arh%set_local) then
   !   !If there is more than one rejection in the semilocal/local area,
   !   !the problem is often related to integral accuracy. Set a keyword
   !   !to get higher integral accuracy: 
   !   write(arh%lupri,*) 'Warning: Too many rejections. Integral accuracy will be increased'
   !   ls%setting%scheme%threshold = ls%setting%scheme%threshold*1.0E-1_realk
   !   !This keyword is reset to false outside arh_solver
   !endif
!endif

   end subroutine update_trustradius


subroutine arhupdate_trustradius_david(CFG,r,ls)
implicit none
type(RedSpaceItem)      :: CFG
real(realk), intent(in) :: r
type(lsitem), intent(in)   :: ls

 CFG%use_max_element=.false.
 CFG%step_accepted = .true.

 if (r.lt.0d0 .and. CFG%r_denom < 0.0d0) then
    CFG%Stepsize = CFG%Stepsize/2.0d0
    CFG%step_accepted = .false.
    CFG%arh%step_accepted = CFG%step_accepted
    write (ls%lupri,*) 'Reject and contract *0.7'
    return
 endif

 if(r.gt.0.75d0) then
         CFG%Stepsize = min(CFG%Stepsize*1.2d0,CFG%max_stepsize)
         write(ls%lupri,*) 'Expand *1.5', r
 end if

 if(r.gt.0.25d0 .and. r.lt.0.75d0)  write(ls%lupri,*) 'Keep stepsize', r
 if (r.ge.0.0d0 .and. r.lt.0.25d0) then
    CFG%Stepsize = CFG%Stepsize*0.7d0
    write(ls%lupri,*) 'Contract *0.7', r
 endif

 CFG%arh%step_accepted = CFG%step_accepted

end subroutine arhupdate_trustradius_david

 end module trustradius_mod
