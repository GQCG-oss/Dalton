!> @file
!> Contains SCF statistics module.

!> \brief Statistics for SCF optimization.
!> \author L. Thogersen. Documented by S. Host.
!> \date 2003
MODULE scf_stats
   !FIXME:       Without the stat_tab information, the computation crashes
   use opttype
   public
   private :: scf_stats_print_table_header,stat_ao_grad,stat_oao_grad!,scf_stats_print_table !Stinne comment
   ! we do not try to collect data from more than 100 iterations.
   !INTEGER, PARAMETER :: max_iterations = 100
   INTEGER,parameter :: trilevel_max_linscf_iterations=100
   !> Current SCF iteration   
   INTEGER,save :: stat_current_iteration
   !> Table for info about SCF iterations
   real(realk),dimension(:,:),pointer,save :: stat_tab
   !> Table for AO and OAO gradient norms for each SCF iteration
   real(realk),dimension(:,:),pointer,save :: stat_gradnorm
   !> SCF energy for each SCF iteration
   real(realk),dimension(:),pointer,save   :: stat_energy
   integer,parameter :: stat_ao_grad=2,stat_oao_grad=1
   contains

   !> \brief Initialize statistics tables.
   !> \author L. Thogersen
   !> \date 2003
   subroutine scf_stats_init(opt)
      implicit none
      !> Contains general settings for SCF optimization
      type(OptItem), intent(in) :: opt
      integer :: maxit,i,j

      maxit = MAX(trilevel_max_linscf_iterations,opt%cfg_max_linscf_iterations)

      NULLIFY(stat_tab,stat_gradnorm,stat_energy)
      ALLOCATE(stat_tab(0:maxit+1,11),stat_gradnorm(0:maxit+1,2),stat_energy(0:maxit+1))

      do j = 1,11
        do i = 0,maxit+1
          stat_tab(i,j) = 0.0E0_realk
        enddo
      enddo

      do i = 0, maxit+1
         stat_gradnorm(i,1) = 0.0E0_realk
         stat_gradnorm(i,2) = 0.0E0_realk
      enddo
   end subroutine scf_stats_init

   !> \brief For debug. When doing ARH (in OAO basis), use this routine to print also AO gradient norm.
   !> \author S. Host
   !> \date February 2010
!!$   subroutine scf_stats_ao_gradnorm(iteration,gradnrm)
!!$      implicit none
!!$      !> Current SCF iteration
!!$      integer, intent(in) :: iteration
!!$      !> Gradient norm in AO basis
!!$      real(realk), intent(in) :: gradnrm
!!$
!!$      stat_current_iteration = iteration
!!$      stat_gradnorm(stat_current_iteration,stat_ao_grad) = gradnrm
!!$   end subroutine scf_stats_ao_gradnorm

   !> \brief Fills statistics tables for current SCF iteration. 
   !> \author L. Thogersen
   !> \date 2003
   subroutine scf_stats_update(iteration,gradnrm,E,opt)
      implicit none
      !> Current SCF iteration
      integer, intent(in) :: iteration
      !> Current SCF AO gradient norm
      real(realk), intent(in) :: gradnrm
      !> Current SCF energy
      real(realk), intent(in) :: E
      !> Contains general settings for SCF optimization
      type(optItem), intent(in) :: opt     
      integer :: j

      stat_current_iteration = iteration
      stat_energy(stat_current_iteration) = E
      IF(opt%cfg_oao_gradnrm)THEN
         stat_gradnorm(stat_current_iteration,stat_oao_grad) = gradnrm
      ELSE
         stat_gradnorm(stat_current_iteration,stat_ao_grad) = gradnrm
      ENDIF
      stat_tab(stat_current_iteration,1) = E
      if (iteration > 1) then
        stat_tab(iteration,2) = stat_tab(iteration,1)-stat_tab(iteration-1,1)
        call scf_stats_print_table(opt,iteration)
        !if (config%diag%DEBUG_RH_DSM_ECHANGE) then
        !  !** col9 for the SCF energy change in the RH-step
        !  stat_tab(iteration,9) = E - stat_tab(iteration,9)
        !endif
      else
        call scf_stats_print_table_header(opt)
        call scf_stats_print_table(opt,iteration)
        if (opt%do_trustregion) then
           call scf_stats_trustradius_header(opt%lupri)
        endif
      endif
   end subroutine scf_stats_update

   !> \brief Print the header for the statistics table.
   !> \author L. Thogersen
   !> \date 2003
   subroutine scf_stats_print_table_header(opt)
      implicit none
      !> Contains general settings for SCF optimization
      type(optItem), intent(in) :: opt

      !if (config%diag%DEBUG_RH_DSM_ECHANGE) then
      !  print"('********************************************************************************&
      !         &*****************************')"
      !  print"(' it        E(HF)           dEHF(RH)         dEHF(DSM)    DSMexit  DSM_alpha  RHshift  RHove&
      !         &rlap    gradient    ###')"
      !  print"('********************************************************************************&
      !         &*****************************')"
      !  WRITE(config%LUPRI,"('********************************************************************************&
      !       &*****************************###')")
      !  WRITE(config%LUPRI,"(' it        E(HF)           dEHF(RH)         dEHF(DSM)    DSMexit  DSM_alpha  RHshift  RHove&
      !         &rlap    gradient    ###')")
      !  WRITE(config%LUPRI,"('********************************************************************************&
      !         &*****************************###')")
      if (opt%cfg_density_method == opt%cfg_f2d_arh) then
        print*,'("************************************************************************************************")'
        IF(opt%cfg_oao_gradnrm)THEN
           print*,'it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift     OAO gradient'
        ELSE
           print*,' it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift     AO gradient'
        ENDIF
        print*,'("************************************************************************************************")'
        WRITE(opt%LUPRI,'("*************************************************************************************************###")')
        IF(opt%cfg_oao_gradnrm)THEN
           WRITE(opt%LUPRI,'(" it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift &
                &   OAO gradient     ###")')
        ELSE
           WRITE(opt%LUPRI,'(" it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift &
                &   AO gradient     ###")')
        ENDIF
        WRITE(opt%LUPRI,'("****************************************************************** &
                         & *******************************###")')
      else
        print*,'("*******************************************************************************************")'
        IF(opt%cfg_oao_gradnrm)THEN
           print*,' it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift  RHinfo  OAO gradient'
        ELSE
           print*,' it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift  RHinfo  AO gradient'
        ENDIF
        print*,'("*******************************************************************************************")'
        WRITE(opt%LUPRI,'("*********************************************************************************************###")')
        IF(opt%cfg_oao_gradnrm)THEN
           WRITE(opt%LUPRI,'(" it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift   RHinfo  OAO gradient    ###")')
        ELSE
           WRITE(opt%LUPRI,'(" it        E(HF)            dE(HF)      DSMexit  DSM_alpha RHshift   RHinfo  AO gradient     ###")')
        ENDIF
        WRITE(opt%LUPRI,'("*********************************************************************************************###")')
      endif
   end subroutine scf_stats_print_table_header

   !> \brief Print the header for trust-radius info (for ARH and 2nd order optimization only).
   !> \author S. Host
   !> \date 2007
   subroutine scf_stats_trustradius_header(lupri)
      implicit none
      !> Logical unit number for output file.
      integer,intent(in) :: lupri

        WRITE(LUPRI,'("*******************************************************************************************%%%")')
        WRITE(LUPRI,'(" Trust Radius   Max element     Norm     RHshift   Ratio  Dpar/Dtot  Ndens(FIFO)    SCF it %%%")')
        WRITE(LUPRI,'("*******************************************************************************************%%%")')
   end subroutine scf_stats_trustradius_header

   !> \brief Print the header for ARH debug info.
   !> \author S. Host
   !> \date 2007
   subroutine scf_stats_arh_header(lupri)
      implicit none
      !> Logical unit number for output file.
      integer,intent(in) :: lupri
                                                                                                                         
        WRITE(LUPRI,'("**********************************************************************///")')
        WRITE(LUPRI,'(" SCF it  GapDiag     GapIter   RedspLowEig  OnestEival  HessianEival  ///")')
        WRITE(LUPRI,'("**********************************************************************///")')
   end subroutine scf_stats_arh_header

   !> \brief Print the statistics info for given iteration.
   !> \author L. Thogersen
   !> \date 2003
   subroutine scf_stats_print_table(opt,iteration)
      implicit none
      !> Contains general settings for SCF optimization
      type(optItem), intent(in) :: opt
      !> Print statistics for this iteration
      integer, intent(in) :: iteration
      integer :: j,igrad

      !if (config%diag%DEBUG_RH_DSM_ECHANGE) then
      !  print"(i3,f18.10,2f17.11,f8.2,f13.5,f8.2,f13.7,e12.2)", iteration-1,stat_tab(iteration-1,1),&
      !         &stat_tab(iteration-1,9),stat_tab(iteration-1,10),stat_tab(iteration-1,6),&
      !         &stat_tab(iteration-1,11),(stat_tab(iteration-1,j),j=7,8),&
      !         &stat_gradnorm(iteration-1,1)
      !  WRITE(config%LUPRI,"(i3,f18.10,2f17.11,f8.2,f13.5,f8.2,f13.7,e12.2,'  ###')") iteration-1,stat_tab(iteration-1,1),&
      !         &stat_tab(iteration-1,9),stat_tab(iteration-1,10),stat_tab(iteration-1,6),&
      !         &stat_tab(iteration-1,11),(stat_tab(iteration-1,j),j=7,8),&
      !         &stat_gradnorm(iteration-1,1)
      IF(opt%cfg_oao_gradnrm)THEN
         igrad = stat_oao_grad
      ELSE
         igrad = stat_ao_grad
      ENDIF
      if (opt%cfg_density_method == opt%cfg_f2d_arh) then
        IF(iteration.EQ.0)THEN
           print '(A5,f18.10,f17.11,f8.2,f13.5,f8.2,1es12.2)', 'Atoms',stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),stat_tab(iteration,7),&
                &stat_gradnorm(iteration,igrad)
           WRITE(opt%LUPRI,'(A5,f18.10,f17.11,f8.2,f13.5,f8.2,1es12.2,"  ###")') 'Atoms',stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),stat_tab(iteration,7),&
                &stat_gradnorm(iteration,igrad)
        ELSE
           print '(i3,f18.10,f17.11,f8.2,f13.5,f8.2,1es12.2)', iteration,stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),stat_tab(iteration,7),&
                &stat_gradnorm(iteration,igrad)
           WRITE(opt%LUPRI,'(i3,f18.10,f17.11,f8.2,f13.5,f8.2,1es12.2,"  ###")') iteration,stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),stat_tab(iteration,7),&
                &stat_gradnorm(iteration,igrad)
        ENDIF
      else
        IF(iteration.EQ.0)THEN
           print '(A5,f18.10,f17.11,f8.2,f13.5,f8.2,f13.7,es12.2)', 'Atoms',stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),(stat_tab(iteration,j),j=7,8),&
                &stat_gradnorm(iteration,igrad)
           WRITE(opt%LUPRI,'(A5,f18.10,f17.11,f8.2,f13.5,f8.2,f13.7,es12.2,"  ###")') 'Atoms',stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),(stat_tab(iteration,j),j=7,8),&
                &stat_gradnorm(iteration,igrad)
        ELSE
           print '(i3,f18.10,f17.11,f8.2,f13.5,f8.2,f13.7,es12.2)', iteration,stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),(stat_tab(iteration,j),j=7,8),&
                &stat_gradnorm(iteration,igrad)
           WRITE(opt%LUPRI,'(i3,f18.10,f17.11,f8.2,f13.5,f8.2,f13.7,es12.2,"  ###")') iteration,stat_tab(iteration,1),&
                &stat_tab(iteration,2),stat_tab(iteration,6),stat_tab(iteration,11),(stat_tab(iteration,j),j=7,8),&
                &stat_gradnorm(iteration,igrad)
        ENDIF
      endif
      !  print '(i3,f18.10,f17.11,f8.2,f13.5,f8.2,f13.7,e12.2)', iteration-1,stat_tab(iteration-1,1),&
      !       &stat_tab(iteration-1,2),stat_tab(iteration-1,6),stat_tab(iteration-1,11),(stat_tab(iteration-1,j),j=7,8),&
      !       &stat_gradnorm(iteration-1)
      !  WRITE(LUPRI,'(i3,f18.10,f17.11,f8.2,f13.5,f8.2,f13.7,e12.2,"  ###")') iteration-1,stat_tab(iteration-1,1),&
      !       &stat_tab(iteration-1,2),stat_tab(iteration-1,6),stat_tab(iteration-1,11),(stat_tab(iteration-1,j),j=7,8),&
      !       &stat_gradnorm(iteration-1)
   end subroutine scf_stats_print_table

   !> \brief Print all statistics info after SCF optimization is finished.
   !> \author L. Thogersen
   !> \date 2003
   subroutine scf_stats_end_print(opt)
      implicit none
      !> Contains general settings for SCF optimization
      type(optItem), intent(in)  :: opt
      integer :: i,j,igrad

      if (stat_current_iteration>=opt%cfg_max_linscf_iterations) then
         stat_current_iteration = opt%cfg_max_linscf_iterations
      endif

      !if (config%diag%DEBUG_RH_DSM_ECHANGE) then
      !  WRITE(config%LUPRI,*)
      !  WRITE(config%LUPRI,"('*************************************************************************')")
      !  WRITE(config%LUPRI,"(' it       dE(RH)          dEHF(RH)         dE(DSM)         dEHF(DSM)')")
      !  WRITE(config%LUPRI,"('*************************************************************************')")
      !  do i=2,stat_current_iteration
      !     WRITE(config%LUPRI,"(i3,4f17.11)") i,stat_tab(i,4),stat_tab(i,9),stat_tab(i,5),stat_tab(i,10)
      !  enddo
      !else
      !  WRITE(opt%LUPRI,*)
      !  WRITE(opt%LUPRI,"('*************************************************************************')")
      !  WRITE(opt%LUPRI,"(' it       dE(HF)        dE(DSM)+dE(RH)      dE(RH)dE(DSM)  ')")
      !  WRITE(opt%LUPRI,"('*************************************************************************')")
      !  do i=1,stat_current_iteration
      !     WRITE(opt%LUPRI,"(i3,4f17.11)") i,(stat_tab(i,j),j=2,5)
      !  enddo
      !endif
      WRITE(opt%LUPRI,*)
      WRITE(opt%LUPRI,"('********************************************************')")
      WRITE(opt%LUPRI,"(' it       dE(HF)       DSMexit   RHshift    RHinfo ')")
      WRITE(opt%LUPRI,"('********************************************************')")
      do i=1,stat_current_iteration
         WRITE(opt%LUPRI,"(i3,f17.11,2f10.4,f13.7)") i,stat_tab(i,2),(stat_tab(i,j),j=6,8)
      enddo

      IF(opt%cfg_oao_gradnrm)THEN
         igrad = stat_oao_grad
      ELSE
         igrad = stat_ao_grad
      ENDIF
      WRITE(opt%LUPRI,*)
      WRITE(opt%LUPRI,'("======================================================================")')
      WRITE(opt%LUPRI,'("                   LINSCF ITERATIONS:")')
      IF(opt%cfg_oao_gradnrm)THEN
         WRITE(opt%LUPRI,'("  It.nr.               Energy                 OAO Gradient norm")')
      ELSE
         WRITE(opt%LUPRI,'("  It.nr.               Energy                 AO Gradient norm")')
      ENDIF
      WRITE(opt%LUPRI,'("======================================================================")')
      IF(stat_energy(0).LT.-1.0E-16_realk)THEN
         WRITE(opt%LUPRI,'(2x,A5,3x,f30.20,5x,d22.15)') 'Atoms',stat_energy(0),stat_gradnorm(0,igrad)
      ENDIF
      do i=1,stat_current_iteration
         WRITE(opt%LUPRI,'(2x,i3,5x,f30.20,5x,d22.15)') i,stat_energy(i),stat_gradnorm(i,igrad)
      enddo
      WRITE(opt%LUPRI,*)
      if (stat_current_iteration < opt%cfg_max_linscf_iterations) then
        WRITE(opt%LUPRI,'("      SCF converged !!!! ")')
        print*,"      SCF converged !!!! "
        !if (cfg_reduced_conv_scf) then
        !    write (lupri,*) 'NOTE: SCF converged to only', cfg_convergence_threshold
        !    write (lupri,*) 'because this is all that can be obtained within numerical accuracy'
        !    write (lupri,*) '(reduced Hessian is only symmetrical to this accuracy).'
        !endif
      else
         if(opt%cfg_max_linscf_iterations.EQ.1)THEN
            WRITE(opt%LUPRI,'("      SCF has failed to converge in ",i3," iterations")') stat_current_iteration
         elseif(.not. opt%opt_quit)THEN
            WRITE(opt%LUPRI,'("      SCF has failed to converge in ",i3," iterations")') stat_current_iteration
         ELSE
            WRITE(opt%LUPRI,'("      SCF has failed to converge in ",i3," iterations")') stat_current_iteration
            CALL lsQUIT('Computation terminated - convergence NOT obtained!!!!',opt%lupri)
         ENDIF
      endif
      WRITE(opt%LUPRI,'("         >>> Final results from LSDALTON <<<")')
      WRITE(opt%LUPRI,*)
      WRITE(opt%LUPRI,*) 
      if (opt%calctype == opt%dftcalc) then
         WRITE(opt%LUPRI,'("      Final DFT energy:              ",f20.12)') stat_energy(stat_current_iteration)
      else if (opt%calctype == opt%hfcalc) then
         WRITE(opt%LUPRI,'("      Final HF energy:               ",f20.12)') stat_energy(stat_current_iteration)
      else
         call lsquit('Calculation type has not been set',opt%lupri)
      endif
      WRITE(opt%LUPRI,'("      Nuclear repulsion:             ",f20.12)') opt%potnuc
      WRITE(opt%LUPRI,'("      Electronic energy:             ",f20.12)') stat_energy(stat_current_iteration)-opt%potnuc
      WRITE(opt%LUPRI,*)

   end subroutine scf_stats_end_print

   !> \brief Deallocate statistics tables.
   !> \author L. Thogersen
   !> \date 2003
   subroutine scf_stats_shutdown()

     DEALLOCATE(stat_tab,stat_gradnorm,stat_energy)
     NULLIFY(stat_tab,stat_gradnorm,stat_energy)
   end subroutine scf_stats_shutdown

END MODULE scf_stats

