!     PARSOPPA = Parallel Second order Polarization Propagator Approximation  (this also includes RPA calculations)

      integer :: NUMDIS, KODCL1, KODCL2, KODBC1,                        & 
     &           KODBC2, KRDBC1, KRDBC2, KODPP1,                        &
     &           KODPP2, KRDPP1, KRDPP2, KCCFB1,                        &
     &           KINDXB, NTOT,   LCMO,                                  &
     &           KCMO  , KEND1 , LWORK1,                                &
     &           LTR1E , LTR1D , LRES1E, LRES1D,                        &
     &           LFOCK , LDENS , LBTR1E, LBTR1D,                        &
     &           LBTJ1E, LBTJ1D, KTR1E , KTR1D ,                        &
     &           KRES1E, KRES1D, KFOCK , KDENS ,                        &
     &           KBTR1E, KBTR1D, KBTJ1E, KBTJ1D,                        &
     &           KEND2,  LWORK2, KENDSV, LWORKSV,                       &
     &           LTR2E,  LTR2D,  LRES2E, LRES2D,                        &
     &           LSIGAI1, LSIGAI2, LSIGDA1, LSIGDA2,                    &
     &           LAIJ, LAAB, KTR2E, KTR2D, KRES2E,                      &
     &           KRES2D, KSIGAI1, KSIGAI2, KSIGDA1,                     &
     &           KSIGDA2, KAIJ, KAAB, isyres,                           &
     &           LT2AMP, parsoppalast,                                  &
     &           icdel1, copyldensai, copynit, copyisymtr,              &
     &           copyinewtr

      logical :: forceupdate
      logical :: getdensai
      character*5 :: copymodel
!  Pad using a dummy variable      
      character(len=3) :: padding_dummy

! Extra communicator needed in case some slaves are excluded from the calculation
      integer(mpi_integer_kind) :: soppa_comm_active

      common /parsoppa/
     & NTOT  , NUMDIS, KODCL1, KODCL2, KODBC1, LWORKSV,                 &
     & KODBC2, KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1,                  &
     & KRDPP2, KCCFB1, KINDXB, KCMO  , KEND1 , LWORK1,                  &
     & LTR1E , LTR1D , LRES1E, LRES1D, LFOCK , LDENS ,                  &
     & LBTR1E, LBTR1D, LBTJ1E, LBTJ1D, KTR1E , KTR1D ,                  &
     & KRES1E, KRES1D, KFOCK , KDENS , KBTR1E, KBTR1D,                  &
     & KBTJ1E, KBTJ1D, KEND2 , LWORK2, KENDSV,                          &
     & LTR2E,  LTR2D,  LRES2E, LRES2D, LSIGAI1,LSIGAI2,                 &
     & LSIGDA1,LSIGDA2,LAIJ,   LAAB,   KTR2E,  KTR2D,                   &
     & KRES2E, KRES2D, KSIGAI1,KSIGAI2,KSIGDA1,KSIGDA2,                 &
     & forceupdate,    LT2AMP, KAIJ,   KAAB,  isyres,                   &
     & LCMO,   icdel1, copyldensai ,                                    &
     & copynit,  copyisymtr, copyinewtr, getdensai,                     &
     & copymodel, padding_dummy


      common /parsoppa/ parsoppalast
      !  Very important !!!
      !  Always keep parsoppalast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

      save /parsoppa/

      common /soppa_communicators/ soppa_comm_active
