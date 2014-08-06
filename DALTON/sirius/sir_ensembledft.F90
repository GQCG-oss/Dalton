!-----------------------------------------------------------------------------------------------------------------------
      subroutine get_sc_ensemble_energy(emydft_orig,EJCSR_orig,         &
     &                                  EJVSR_orig,EDSR_orig,           &
     &                                  EDFT_orig,EMYDFTAUX_orig,       &
     &                                  eci_orig,emy_ci_orig,           &
     &                                  energy_root_orig,eensemb_orig,  &
     &                                  cmo,srac_orig,nnashx,nnorbt,    &
     &                                  ncroot,istaci,weights,lupri,    &
     &                                  potnuc,                         &
     &                                  wrk,lwrk)

      !> brief: calculate the ensemble-DFT energy in a self-consistent manner (see notes by E. Fromager)
      !> author: Stefan Knecht and Manu Fromager
!-----------------------------------------------------------------------------------------------------------------------
      implicit none

      integer, intent(in)    :: ncroot,istaci,lupri,lwrk,nnashx,nnorbt
      real*8 , intent(in)    :: emydft_orig, EJCSR_orig, EJVSR_orig
      real*8 , intent(in)    :: EDFT_orig, EMYDFTAUX_orig, emy_ci_orig
      real*8 , intent(in)    :: eensemb_orig, EDSR_orig
      real*8 , intent(in)    :: eci_orig(ncroot)
      real*8 , intent(in)    :: energy_root_orig(ncroot)
      real*8 , intent(in)    :: potnuc
      real*8 , intent(in)    :: weights(ncroot), srac_orig(nnashx)
      real*8 , intent(in)    :: cmo(*)
      real*8 , intent(inout) :: wrk(lwrk)
!-----------------------------------------------------------------------------------------------------------------------
      integer                :: kw2, KW2A, KUDVREF
      integer                :: LW2A, ist
      real*8                 :: uESRDV, uEDSR, uEJVSR, uEJCSR
      real*8                 :: uEMYDFTAUX, uedft
      real*8                 :: UEJCVSR
      real*8                 :: VEENSEMB

      real*8 , allocatable   :: u_rho_ensemble(:)
      real*8 , allocatable   :: ufsr(:)
      real*8 , external      :: ddot

      kw2 = 1

      !> build new ensemble density matrix
      allocate(u_rho_ensemble(nnashx))
      u_rho_ensemble = 0

      do IST = 1,NCROOT
        KUDVREF = KW2
        KW2A    = KUDVREF + NNASHX
        LW2A    = LWRK   - KW2A
        CALL GETUDVREF('LUCITA   ',WRK(KUDVREF),IST,wrk,WRK(KW2A),LW2A)
        call daxpy(nnashx,weights(ist),WRK(KUDVREF),1,u_rho_ensemble,1)
      end do
      kw2 = 1

      allocate(ufsr(nnorbt))
      ufsr  = 0

      CALL SRFMAT(ufsr,cmo,u_rho_ensemble,                              &
     &            uEJCSR,uEJVSR,uEDSR,uEDFT,                            &
     &            uEMYDFTAUX,uEJCVSR,                                   &
     &            wrk(kw2),lwrk,0)

      uESRDV   = DDOT(NNASHX,u_rho_ensemble,1,srac_orig,1)

      VEENSEMB = eensemb_orig + uEDFT  - EDFT_orig - EDSR_orig + uEJCVSR&
     &         - uESRDV       - uEJVSR - EJVSR_orig

      write(lupri,*) 'contributions to the variational ensemble energy'
      write(lupri,*) '------------------------------------------------'
      write(lupri,*) 'eensemb_orig --> ',eensemb_orig
      write(lupri,*) 'uEDFT        --> ',uEDFT
      write(lupri,*) ' EDFT_orig   --> ',EDFT_orig
      write(lupri,*) ' EDSR_orig   --> ',EDSR_orig
      write(lupri,*) 'uEJCVSR      --> ',uEJCVSR
      write(lupri,*) 'uESRDV       --> ',uESRDV
      write(lupri,*) 'uEJVSR       --> ',uEJVSR
      write(lupri,*) ' EJVSR_orig  --> ',EJVSR_orig
      write(lupri,*) '------------------------------------------------'

      write(lupri,806) '*******************************'
      write(lupri,806) 'Final variational ensemble CI-srDFT energy:',   &
     &                  VEENSEMB
      write(lupri,806) '*******************************'
!         
  805   FORMAT( 1X,A,F25.12)
  806   FORMAT(/1X,A,F25.12)
  807   FORMAT(/1X,A,I2,A,F20.12)

      deallocate(u_rho_ensemble,ufsr)

      end subroutine get_sc_ensemble_energy
!-----------------------------------------------------------------------------------------------------------------------

