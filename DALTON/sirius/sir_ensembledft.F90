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
      integer                :: kw2, KDVREF, KW2A, KSRAC
      integer                :: LW2A, ist
      real*8                 :: EJTOT, ESRDV, EDSR, EJVSR, EJCSR, ETDFT
      real*8                 :: EMYDFT, ELRCI, ECIAUX, EMYDFTAUX
      real*8                 :: ECIAUX_gs, EENSEMB

      real*8 , allocatable   :: eci(:)
      real*8 , external      :: ddot

      allocate(eci(ncroot))
      eci = 0

      write(lupri,'(/A,I5)') 'SR 0-el. energy for input .STATE',        &
     &                          ISTACI
      write(lupri,'(A/)') '-------------------------------------'
      write(lupri,805) '  SR core    Hartree energy   :',EJCSR_orig
      write(lupri,805) '- SR valence Hartree energy   :',EJVSR_orig
      write(lupri,805) '+ SR Exchange-correlation     :',EDFT_orig
      write(lupri,805) '- SR Exchange-correlation pot.:',EDSR_orig
      write(lupri,806) '= Total eff. SR 0-el. energy  :',EMYDFT_orig

      kw2 = 1

      DO IST = 1,NCROOT
        KDVREF = KW2
        KW2A   = KDVREF + NNASHX
        LW2A   = LWRK   - KW2A
        CALL GETDVREF('LUCITA   ',WRK(KDVREF),IST,wrk,WRK(KW2A),LW2A)
        !> DVREF is equal to DV for state abs(ISTI), i.e. 1-el. energy can be calculated
        ESRDV = DDOT(NNASHX,WRK(KDVREF),1,srac_orig,1)
!       ESRDV = DDOT(NNASHX,WRK(KDVREF),1,WRK(KSRAC),1)

        WRITE(LUPRI,'(/A,I5)') 'CI-DFT energy for state no.',IST
        WRITE(LUPRI,'(A/)') '-------------------------------------'
        WRITE(LUPRI,805) 'SR eff. 1-el. energy          :',ESRDV
        EJTOT = ESRDV + EDSR + EJVSR + EJCSR
        WRITE(LUPRI,805) 'SR total Hartree energy       :',EJTOT
        ETDFT = EMYDFT + ESRDV
        WRITE(LUPRI,805) 'SR eff. total DFT energy      :',ETDFT
        ELRCI = ECI(IST) - ETDFT
        WRITE(LUPRI,805) 'LR total CI  energy           :',ELRCI
        WRITE(LUPRI,806) 'Total CI-DFT energy           :',ECI(IST)
        ECIAUX = ELRCI+EMYDFTAUX+ESRDV-POTNUC 
!
        if (IST .EQ. 1) then 
        ECIAUX_gs = ECIAUX ! saving the GS auxiliary energy 
        end if
!       write(lupri,*) 'eciaux_gs =', ECIAUX_gs
!
        write(lupri,'(/a)')                                             &
     &  ' decomposition of the auxiliary CI-srDFT energy'
        WRITE(LUPRI,805) 'ELRCI        :',ELRCI
        WRITE(LUPRI,805) 'EMYDFTAUX    :',EMYDFTAUX
        WRITE(LUPRI,805) 'ESRDV        :',ESRDV
        WRITE(LUPRI,805) 'POTNUC       :',POTNUC
        write(lupri,*) ''
        WRITE(LUPRI,807) 'Auxiliary CI-srDFT energy for root',          &
     &                    IST, ':   ', ECIAUX
        if (IST .GT. 1) then
!       WRITE(LUPRI,806) 'Auxiliary excitation energy   :',
!    &                    ECIAUX-ECIAUX_gs
        WRITE(LUPRI,807) 'Auxiliary excitation energy for root',        &
     &                    IST, ': ', ECIAUX-ECIAUX_gs
        endif
!       compute the ensemble CI-srDFT energy 
        EENSEMB = EENSEMB + weights(IST) * ECI(IST)
!       compute the ensemble ESRDV energy contribution 
!       ESRDV_ens = ESRDV_ens + weights(IST) * ESRDV
      END DO

      write(lupri,806) '*******************************'
      write(lupri,806) 'Final ensemble CI-srDFT energy:',EENSEMB
      write(lupri,806) '*******************************'
!         
  805   FORMAT( 1X,A,F25.12)
  806   FORMAT(/1X,A,F25.12)
  807   FORMAT(/1X,A,I2,A,F20.12)

      deallocate(eci)

      end subroutine get_sc_ensemble_energy
!-----------------------------------------------------------------------------------------------------------------------

