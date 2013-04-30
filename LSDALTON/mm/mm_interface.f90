! Notes: traditional dalton's WORK memory is used for some of the
! larger allocations before memory allocation in dalton is "fixed".
!
MODULE mm_interface_mod

   use files
   USE mm_global_paras_mod
   USE mm_stats_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_connect_interface,        &
             mm_disconnect_interface,     &
             mm_get_raw_data,             &
             mm_check_moms

   ! we use pointers here so we can allocate outside of the module
   TYPE(raw_mm_paras), POINTER, SAVE :: raw_paras(:)
   TYPE(J_index_type), POINTER, SAVE :: J_indices(:)
   REAL(REALK),        POINTER, SAVE :: dens(:)
   REAL(REALK),        POINTER, SAVE :: dens_in(:)
   REAL(REALK),        POINTER, SAVE :: densf_in(:)
   REAL(REALK),        POINTER, SAVE :: raw_qlm(:,:)
   REAL(REALK),        POINTER, SAVE :: raw_qlm_der(:,:,:)
   INTEGER,      POINTER, SAVE :: mom2atom(:,:)
   ! Number of multipole moments generated externally
   TYPE(mm_counters), SAVE :: n_mms 

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_connect_interface(scheme,WORK,KWRK,LWORK)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      INTEGER,      INTENT(IN)    :: LWORK
      INTEGER,      INTENT(INOUT) :: KWRK
      REAL(REALK),        INTENT(INOUT) :: WORK(LWORK)
      
      LOGICAL :: do_center
     
      NULLIFY (dens, raw_qlm, raw_paras, J_indices,dens_in,densf_in)
      IF ( stat_do_grad ) THEN   
         NULLIFY (raw_qlm_der)
         NULLIFY (mom2atom)
      ENDIF

      IF ( .NOT. MM_STANDALONE ) LUPRI = 6 ! Stinne 27/10-2010 CALL mm_GTUNIT(LUPRI)
      CALL mm_init_moments(scheme%raw_LMAX,WORK,KWRK,LWORK)
      CALL mm_read_in_raw_data(stat_do_grad)
!      CALL mm_check_moms(scheme%raw_LMAX,raw_qlm)
      CALL mm_translate_centres(scheme%system_size)

   END SUBROUTINE mm_connect_interface

!-------------------------------------------------------------------------------

   SUBROUTINE mm_disconnect_interface

      IMPLICIT NONE
      n_mms%elec = 0
      n_mms%nuc  = 0
      n_mms%tot  = 0
      n_mms%nbast = 0
      n_mms%nauxbas = 0
      CALL mm_free_moments

   END SUBROUTINE mm_disconnect_interface

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_moments(LMAX,WORK,KWRK,LWORK)
 
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: LMAX, LWORK
      INTEGER, INTENT(INOUT) :: KWRK
      REAL(REALK),   INTENT(INOUT) :: WORK(LWORK)

      CALL mm_get_n_mms_from_file(LMAX)
      IF(stat_do_grad) THEN 
        CALL mm_allocate_mms_arrays(LMAX,WORK,LWORK/(LMAX+1)**2/(6+1),6)
        KWRK = KWRK + ((LMAX+1)**2)*n_mms%tot*(6 + 1)
      ELSE
         CALL mm_allocate_mms_arrays(LMAX,WORK,LWORK/(LMAX+1)**2,0)
         KWRK = KWRK + ((LMAX+1)**2)*n_mms%tot
      ENDIF
   END SUBROUTINE mm_init_moments

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_n_mms_from_file(LMAX_in)
     USE density_fitting
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LMAX_in
      INTEGER :: LMAX, NBAST, ID, LUINTM,NAUXMOM
      INTEGER :: LMAX_DER,NBAST_DER,ID_DER, nnuc, nelec, nnuc_der, nelec_der, NAUXMOM_der
      REAL(REALK)   :: RD
      LOGICAL       :: ldumm1, ldumm2, LHS_AUX_der, RHS_AUX_der

      ldumm1=.false.
      ldumm2=.false.

      ! get fitting information
      LUINTM = 0
      CALL LSOPEN(LUINTM,INPUT_FILE0,'UNKNOWN','UNFORMATTED')
      REWIND (LUINTM)
      READ (LUINTM)  fit%LHS_AUX, fit%RHS_AUX
      CALL LSCLOSE(LUINTM,'KEEP')

      ! read undifferentiated moments
      LUINTM = 0
      CALL LSOPEN(LUINTM,INPUT_FILE2,'UNKNOWN','UNFORMATTED')
      REWIND (LUINTM)
      READ (LUINTM) LMAX, NBAST, nelec, nnuc, ID, &
     &           NAUXMOM, ldumm1, ldumm2
      CALL LSCLOSE(LUINTM,'KEEP')

      IF(stat_do_grad) THEN
         ! get differentiated info     
         LUINTM = 0
         CALL LSOPEN(LUINTM,INPUT_FILE9,'UNKNOWN','UNFORMATTED')
         REWIND (LUINTM)
         READ (LUINTM) LMAX_DER, NBAST_DER, nelec_der, nnuc_der, ID_der, &
     &              NAUXMOM_der, LHS_AUX_der, RHS_AUX_der
         CALL LSCLOSE(LUINTM,'KEEP')

         !consistency check
         IF (LMAX_DER    .NE. LMAX   ) CALL LSQUIT('ERROR 1 in get_n_mms_from_file',-1)
         IF (NBAST_DER   .NE. NBAST  ) CALL LSQUIT('ERROR 2 in get_n_mms_from_file',-1)
         IF (nelec_der   .NE. nelec  ) THEN
            print*,'nelec_der',nelec_der
            print*,'nelec',nelec
            CALL LSQUIT('ERROR 3 in get_n_mms_from_file',-1)
         ENDIF
         IF (nnuc_der    .NE. nnuc   ) CALL LSQUIT('ERROR 4 in get_n_mms_from_file',-1)
         IF (NAUXMOM_der .NE. NAUXMOM) CALL LSQUIT('ERROR 5 in get_n_mms_from_file',-1)
         fit%LHS_AUX   = LHS_AUX_der
         fit%RHS_AUX   = RHS_AUX_der
      END IF

      n_mms%elec  = nelec
      n_mms%nuc   = nnuc
      fit%NAUXMOM = NAUXMOM
      fit%AUX       = fit%LHS_AUX .OR. fit%RHS_AUX
      n_mms%tot     = n_mms%elec + n_mms%nuc
      n_mms%nbast   = NBAST
      n_mms%nauxbas = fit%NAUXMOM
      stat_raw_moms = n_mms%tot
      stat_nuc_moms = n_mms%nuc 

      IF (LMAX .NE. LMAX_in) CALL lsQUIT('LMAX inconsistency in MM interface!',-1)

   END SUBROUTINE mm_get_n_mms_from_file

!-------------------------------------------------------------------------------

   SUBROUTINE mm_allocate_mms_arrays(LMAX,WORK,maxmom,nucder)

      USE mm_memory_manager_mod, ONLY: mm_allocate, mm_deallocate

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LMAX, maxmom, nucder
      REAL(REALK),   TARGET     :: WORK((LMAX+1)**2,maxmom,nucder+1)
      INTEGER :: i

      CALL mm_allocate(MEM_RAW_QLM,raw_paras,n_mms%tot)

      ALLOCATE( J_indices(n_mms%tot) )
      ! initialise parameters
      DO i = 1, n_mms%tot 
         raw_paras(i)%cntr = zero
         raw_paras(i)%ext = zero
         raw_paras(i)%id = 0
         raw_paras(i)%batch = 0
         raw_paras(i)%Lmin = 0
         raw_paras(i)%map_up = 0
         raw_paras(i)%box = 0 
         raw_paras(i)%bra = 0
         raw_paras(i)%box_cntr = zero
         J_indices(i)%i_indx = 0
         J_indices(i)%j_indx = 0
      END DO

      CALL mm_allocate(MEM_RAW_QLM,dens,n_mms%tot)
      CALL mm_allocate(MEM_RAW_QLM,dens_in,n_mms%nbast*(n_mms%nbast-1)/2+n_mms%nbast)
      CALL mm_allocate(MEM_RAW_QLM,densf_in,n_mms%nauxbas)
      ! use work for the stuff below:
      CALL mm_allocate(MEM_RAW_QLM,raw_qlm,(LMAX+1)**2,n_mms%tot)
      IF ( stat_do_grad ) THEN
         CALL mm_allocate(MEM_RAW_QLM,raw_qlm_der,(LMAX+1)**2,n_mms%tot,nucder)
      ENDIF

      IF ( stat_do_grad ) THEN
         CALL mm_allocate(MEM_RAW_QLM,mom2atom,n_mms%tot,2)
      END IF

      raw_qlm(:,:) = zero  ! crucial since only non-zero written explicitly
      IF ( stat_do_grad ) THEN
         raw_qlm_der(:,:,:) = zero
      END IF

   END SUBROUTINE mm_allocate_mms_arrays

!------------------------------------------------------------------------------

   SUBROUTINE mm_free_moments
     USE mm_stats_mod
     USE mm_memory_manager_mod, ONLY: mm_deallocate
      IMPLICIT NONE

      IF (.NOT.ASSOCIATED(dens))      STOP "dens should be allocated!"
      IF (.NOT.ASSOCIATED(raw_qlm))   STOP "raw_qlm should be allocated!"
      IF (.NOT.ASSOCIATED(raw_paras)) STOP "raw_paras should be allocated!"
      IF (.NOT.ASSOCIATED(J_indices)) STOP "J_indices should be allocated!"
      IF (stat_do_grad) THEN
         IF (.NOT.ASSOCIATED(raw_qlm_der))   STOP "raw_qlm_der should be allocated!"
         IF (.NOT.ASSOCIATED(mom2atom))      STOP "mom2atom should be allocated!"
      END IF
      IF (ASSOCIATED(raw_qlm))THEN
         call mm_deallocate(raw_qlm)
      ENDIF
      IF (stat_do_grad) THEN
         IF (ASSOCIATED(raw_qlm_der))THEN
            call mm_deallocate(raw_qlm_der)
         ENDIF
      END IF
      IF (stat_do_grad) THEN
         IF (ASSOCIATED(mom2atom))THEN
            call mm_deallocate(mom2atom)
         ENDIF
      END IF
      IF (ASSOCIATED(dens))THEN
         call mm_deallocate(dens)
      ENDIF
      IF (ASSOCIATED(dens_in))THEN
         call mm_deallocate(dens_in)
      ENDIF
      IF (ASSOCIATED(densf_in))THEN
         call mm_deallocate(densf_in)
      ENDIF
      IF (ASSOCIATED(raw_paras))THEN
         call mm_deallocate(raw_paras)
      ENDIF
      IF (ASSOCIATED(J_indices))THEN
         DEALLOCATE(J_indices)
      ENDIF
      NULLIFY (raw_paras, dens, raw_qlm, J_indices)
      IF (stat_do_grad) THEN
         NULLIFY (raw_qlm_der,mom2atom)
      END IF

   END SUBROUTINE mm_free_moments

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_raw_data(scheme,LHS_mms,RHS_mms)
     USE mm_global_paras_mod
     USE mm_stats_mod
     USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(OUT) :: LHS_mms, RHS_mms
      INTEGER :: lo,hi,hil,lol, startreg, startaux, startnuc, endreg, endaux, endnuc

      if (stat_reorder_mom) then
         ! moments on file are normally ordered as regular, auxiliary, nuclear
         ! we changed however the order into regular, nuclear, auxiliary
         ! we have to account for this change here as well
         startreg = 1
         endreg   = n_mms%elec - fit%NAUXMOM
         startnuc = endreg + 1
         endnuc   = endreg + n_mms%nuc
         startaux = endnuc + 1
         endaux   = n_mms%tot
      else
         startreg = 1
         endreg   = n_mms%elec - fit%NAUXMOM
         startaux = endreg + 1
         endaux   = n_mms%elec
         startnuc = endaux + 1
         endnuc   = n_mms%tot
      end if

      ! get LHS range multipole data
      IF (scheme%LHS_mm_range == NUCLEAR_ONLY) THEN
         lo = startnuc
         hi = endnuc
         LHS_mms%startnuc  = 1
         LHS_mms%endnuc    = n_mms%nuc
      ELSEIF (scheme%LHS_mm_range == ELECTRONIC_ONLY) THEN
         lo = startreg
         hi = endreg
         LHS_mms%startnuc  = 1 + n_mms%elec - fit%NAUXMOM
         LHS_mms%endnuc    = LHS_mms%startnuc
         IF (fit%LHS_AUX) THEN
            lo = startaux
            hi = endaux
            LHS_mms%startnuc  = fit%NAUXMOM + 1
            LHS_mms%endnuc    = LHS_mms%startnuc
         END IF
      ELSEIF (scheme%LHS_mm_range == ALL_MOMENTS ) THEN
         lo = startreg
         hi = endnuc
         LHS_mms%startnuc  = endreg + 1
         LHS_mms%endnuc    = LHS_mms%startnuc + n_mms%nuc - 1
         IF (fit%LHS_AUX) THEN
            if (stat_reorder_mom) then
               lo = startnuc
               hi = endaux
               LHS_mms%startnuc  = 1
               LHS_mms%endnuc    = n_mms%nuc
            else
               lo = startaux
               hi = endnuc
               LHS_mms%startnuc  = 1 + fit%NAUXMOM
               LHS_mms%endnuc    = LHS_mms%startnuc + n_mms%nuc - 1
            end if
         END IF
      ELSE
         WRITE(*,*) 'ERROR in mm_get_raw_data: range:', scheme%LHS_mm_range,'LHS_AUX', fit%LHS_AUX
         CALL lsquit('unrecognized option in mm_get_raw_data',-1)
      END IF
      LHS_mms%qlm       => raw_qlm(:,lo:hi)
      LHS_mms%dens      => dens(lo:hi)
      LHS_mms%paras     => raw_paras(lo:hi)
      LHS_mms%J_indices => J_indices(lo:hi)
      IF (stat_do_grad) THEN
         lol = startnuc
         hil = endnuc
         raw_qlm_der(:,lol:hil,:) = raw_qlm_der(:,lol:hil,:)*0.5D0
         LHS_mms%qlm_der  => raw_qlm_der(:,lo:hi,:)
         LHS_mms%mom2atom => mom2atom(lo:hi,:)
      END IF

      ! get RHS range multipole data
      IF (scheme%RHS_mm_range == NUCLEAR_ONLY) THEN
         lo = startnuc
         hi = endnuc
      ELSEIF (scheme%RHS_mm_range == ELECTRONIC_ONLY) THEN
         lo = startreg
         hi = endreg
         IF (fit%RHS_AUX) THEN
            lo = startaux
            hi = endaux
         END IF
      ELSEIF (scheme%RHS_mm_range == ALL_MOMENTS) THEN
         lo = startreg
         hi = endnuc
         IF(stat_do_grad.and. .not. fit%RHS_AUX) THEN
            raw_qlm(:,startreg:endreg) = 0.5D0*raw_qlm(:,startreg:endreg)
         ENDIF
         IF (fit%RHS_AUX) THEN
            IF (stat_reorder_mom) THEN
               lo = startnuc
               hi = endaux
            ELSE
               lo = startaux
               hi = endnuc
            END IF
         END IF
      ELSEIF (scheme%RHS_mm_range == REG_MIN_AUX .AND. fit%RHS_AUX) THEN
          ! for the deriavtives we need the potential due to the difference in regular and auxiliary
          ! density |rho - rho_fit), therefore we here multiply the fitting moments by -1
          if(stat_reorder_mom) CALL lsQUIT('REG_MIN_AUX not with reordered moments',-1)
          if(.NOT. stat_do_grad) CALL lsQUIT('REG_MIN_AUX only for gradient calculations',-1)
          raw_qlm(:,startaux:endaux) = -1.0D0*raw_qlm(:,startaux:endaux)
          raw_qlm(:,startreg:endreg) = raw_qlm(:,startreg:endreg)*0.5D0
          lo = startreg
          hi = endaux
      ELSE
         WRITE(*,*) 'ERROR in mm_get_raw_data: range:', scheme%RHS_mm_range,'RHS_AUX', fit%RHS_AUX
         CALL lsquit('unrecognized option in mm_get_raw_data',-1)
      END IF
      RHS_mms%qlm       => raw_qlm(:,lo:hi)
      RHS_mms%dens      => dens(lo:hi)
      RHS_mms%paras     => raw_paras(lo:hi)
      RHS_mms%J_indices => J_indices(lo:hi)

!ANDREAS: do we every need the RHS_mms%qlm_der??
!     IF (stat_do_grad) RHS_mms%qlm_der => raw_qlm_der(:,lo:hi,:)

   END SUBROUTINE mm_get_raw_data

!-------------------------------------------------------------------------------
! Read in multipole moment data from interface file.
! In all this MM code we assume the order of moments is:
!   (0),(-1,0,1),(-2,-1,0,1,2)...

   SUBROUTINE mm_read_in_raw_data(do_grad)

     use mm_mem, only: mem_alloc_fmm,mem_dealloc_fmm
     USE density_fitting

      IMPLICIT NONE

      REAL(REALK)   :: EXTP, PX, PY, PZ, SPH
      REAL(REALK)   :: SPHDER, FAC, RAW_NUCDER, gradfactor
      INTEGER :: L,M, A,B, LMINA, LMINB, IBTCH, I, p
      INTEGER :: ID, LUINTM,MMBUFLEN, MAXBUFI,MAXBUFR,MAXBUFN, MXBUFI, MXBUFR
      INTEGER :: ISCOOR, J, IBUFFER, ISTART, IEND
      INTEGER :: ICENTA, ICENTB,II,IDONE
      INTEGER :: LUINTMI, LUINTMR, ALLOC_ERR, NREDI, NREDR, INI, INR
      INTEGER :: nreg, naux, nnuc, nele, n
      REAL(REALK),   ALLOCATABLE :: RBUF(:)
      INTEGER, ALLOCATABLE :: IBUF(:)
      LOGICAL :: USEMMBUF
      LOGICAL :: DO_GRAD

      gradfactor = 1.0E0_realk
      IF (do_grad) gradfactor = 0.5E0_realk ! CAREFUL: Factor introduced for closed shells

      nreg = n_mms%elec - fit%NAUXMOM  ! number of regular moments
      naux = fit%NAUXMOM               ! number of fitting moments
      nele = n_mms%elec                ! number of electronic (regular+fitting) moments
      nnuc = n_mms%nuc                 ! number of nuclear moments

      ! get info about use of buffer 
      CALL GETMMBUFINFO(USEMMBUF,MMBUFLEN,MAXBUFI,MAXBUFR,MAXBUFN)

      IF (USEMMBUF) THEN
         ! use buffer to read in from file 
         IF (MAXBUFI .NE.11) CALL lsQUIT('wrong MAXBUFI in mm_read_in_raw_data',-1)
         IF (MAXBUFR .NE. 5) CALL lsQUIT('wrong MAXBUFR in mm_read_in_raw_data',-1)
         IF (MAXBUFN .NE. 4) CALL lsQUIT('wrong MAXBUFN in mm_read_in_raw_data',-1)
         MXBUFI = MMBUFLEN*MAXBUFI
         MXBUFR = MMBUFLEN*MAXBUFR
         ! allocate memory for the buffer
         !call mem_alloc_fmm(IBUF,MXBUFI)
         ALLOCATE (IBUF(MXBUFI),STAT=ALLOC_ERR)
         IF (alloc_err .NE. 0) CALL lsQUIT('Memory allocation failed in mm_read_in_raw_data(I).',-1)
         !call mem_alloc_fmm(RBUF,MXBUFR)
         ALLOCATE (RBUF(MXBUFR),STAT=ALLOC_ERR)
         IF (alloc_err .NE. 0) CALL lsQUIT('Memory allocation failed in mm_read_in_raw_data(II).',-1)

         ! open MM data files
         LUINTMI = 0
         LUINTMR = 0
         CALL LSOPEN(LUINTMI,INPUT_FILE1,'UNKNOWN','UNFORMATTED')
         CALL LSOPEN(LUINTMR,INPUT_FILE8,'UNKNOWN','UNFORMATTED')
         REWIND (LUINTMI)
         REWIND (LUINTMR)
         ! read electronic multipole moments into buffer
         readloopA: DO
            READ (LUINTMI) NREDI
            READ (LUINTMR) NREDR
            IF (NREDI .EQ. -1) EXIT readloopA
            IF (NREDR/5 .NE. NREDI/11 ) CALL lsQUIT ("ERROR while reading the MM from file",-1)
            READ (LUINTMI) (IBUF(II),II=1,NREDI)
            READ (LUINTMR) (RBUF(II),II=1,NREDR)
            ! put buffer content to the corresponding places
            INR = 1
            DO INI = 1, NREDI-10,11
               ! integer buffer
               L      = IBUF(INI)
               M      = IBUF(INI+1)
               A      = IBUF(INI+2)
               B      = IBUF(INI+3)
               LMINA  = IBUF(INI+4)
               LMINB  = IBUF(INI+5)
               IBTCH  = IBUF(INI+6)
               I      = IBUF(INI+7)
               ICENTA = IBUF(INI+8)
               ICENTB = IBUF(INI+9)
               ISCOOR = IBUF(INI+10)
               ! real buffer
               EXTP  = RBUF(INR)
               PX    = RBUF(INR+1)
               PY    = RBUF(INR+2)
               PZ    = RBUF(INR+3)
               SPH   = RBUF(INR+4)
               INR = INR + 5
               IF (I > SIZE(raw_qlm,2))          CALL lsQUIT ('interface file error 1',-1)
               IF (((L+1)**2) > SIZE(raw_qlm,1)) CALL lsQUIT ('interface file error 2',-1)
                ! reorder data eventually:
                ! the moments on file are ordered as regular, auxiliary, nuclear
                ! we want the order regular, nuclear, auxiliary.
                ! therefore the auxiliary moments are shiftet \Andreas Krapp
               IF (stat_reorder_mom .AND. I .GT. nreg .AND. I .LE. nele) THEN
                  I = I + nnuc
                  IBTCH = IBTCH + nnuc
               END IF
                ! indices to map moments to orbitals in J-matrix
               J_indices(I)%i_indx = A
               J_indices(I)%j_indx = B
                ! index to batches of unique centre and extent
               raw_paras(I)%batch  = IBTCH
                ! minimum order for exact treatment (if uncontracted basis set)
               raw_paras(I)%Lmin = LMINA + LMINB
               raw_paras(I)%cntr = (/ PX, PY, PZ /)
               raw_paras(I)%ext  = EXTP
                ! see defn p.424 "Electronic Structure Theory", Helgaker et al
               p = L*(L+1) +M +1
                ! components (l,m) of MM without density factorised in 
               IF(do_grad) SPH = SPH*gradfactor
               raw_qlm(p,I) = SPH
            END DO
         END DO readloopA

         IF(NREDI .NE. -1) CALL lsQUIT("ERROR2 while reading MM from file",-1)

         ! next read nuclei data: charge, and location
         idone = 0
         n = n_mms%elec
         if (stat_reorder_mom) n = nreg
         readloopB: DO
            READ (LUINTMR) NREDR
            IF(NREDR.LE.0) EXIT readloopB
            READ (LUINTMR) (RBUF(II),II=1,NREDR)
            DO I = 1, NREDR/4
               p = n + I + idone
               II = (I-1)*4+1
               raw_qlm(1,p) = RBUF(II)
               raw_paras(p)%cntr(1)= RBUF(II+1)
               raw_paras(p)%cntr(2)= RBUF(II+2)
               raw_paras(p)%cntr(3)= RBUF(II+3)
               raw_paras(p)%batch = raw_paras(n)%batch + I + idone
               IF(do_grad) THEN
                  raw_qlm(1,p) = RBUF(II)*gradfactor 
                  raw_nucder=-raw_qlm(1,p)
                  raw_qlm_der(1,p,1) = raw_nucder
                  raw_qlm_der(1,p,2) = raw_nucder
                  raw_qlm_der(1,p,3) = raw_nucder
                  mom2atom(p,1)       = I
                  mom2atom(p,2)       = I
               END IF
               raw_paras(p)%Lmin   = 0     ! spherical charge
               ! the artificial extent of the nuclei is set to 1.0D0
               ! this has to match the settings in the routines 
               !    II_get_nn_gradient (LSint/IntegralInterface.f90) 
               !    getODcenter (LSint/ODbatch.f90) 
               raw_paras(p)%ext    = 1.0d0 ! point charges
               J_indices(p)%i_indx = 0     ! not relevant
               J_indices(p)%j_indx = 0     ! not relevant
               dens(p) = one
            END DO
            idone = idone + NREDR/4
         END DO readloopB

         IF (idone .NE. n_mms%nuc) CALL lsQUIT("ERROR3 while reading MM from file",-1)

         IF (NREDR.NE.-2) CALL lsQUIT("ERROR4 while reading MM from file",-1)

         ! close files
         CALL LSCLOSE(LUINTMI,'KEEP')
         CALL LSCLOSE(LUINTMR,'KEEP')

         ! open MM derivative data files
         IF(do_grad) THEN
            LUINTMI = 0
            LUINTMR = 0
            CALL LSOPEN(LUINTMI,INPUT_FILE10,'UNKNOWN','UNFORMATTED')
            CALL LSOPEN(LUINTMR,INPUT_FILE11,'UNKNOWN','UNFORMATTED')
            REWIND (LUINTMI)
            REWIND (LUINTMR)
            ! read electronic multipole moment derivatives into buffer
            readloopC: DO
               READ (LUINTMI) NREDI
               READ (LUINTMR) NREDR
               IF (NREDI .EQ. -1) EXIT readloopC
               IF (NREDR/5 .NE. NREDI/11 ) CALL lsQUIT ("ERROR while reading the MM from file",-1)
               READ (LUINTMI) (IBUF(II),II=1,NREDI)
               READ (LUINTMR) (RBUF(II),II=1,NREDR)
               ! put buffer content to the corresponding places
               INR = 1
               DO INI = 1, NREDI-10,11
                  ! integer buffer
                  L      = IBUF(INI)
                  M      = IBUF(INI+1)
                  A      = IBUF(INI+2)
                  B      = IBUF(INI+3)
                  LMINA  = IBUF(INI+4)
                  LMINB  = IBUF(INI+5)
                  IBTCH  = IBUF(INI+6)
                  I      = IBUF(INI+7)
                  ICENTA = IBUF(INI+8)
                  ICENTB = IBUF(INI+9)
                  ISCOOR = IBUF(INI+10)
                  ! real buffer
                  EXTP  = RBUF(INR)
                  PX    = RBUF(INR+1)
                  PY    = RBUF(INR+2)
                  PZ    = RBUF(INR+3)
                  SPHDER= RBUF(INR+4)
                  INR = INR + 5
                  IF (I > SIZE(raw_qlm,2))          CALL lsQUIT ('interface file error 1',-1)
                  IF (((L+1)**2) > SIZE(raw_qlm,1)) CALL lsQUIT ('interface file error 2',-1)
                   ! reorder data eventually:
                   ! the moments on file are ordered as regular, auxiliary, nuclear
                   ! we want the order regular, nuclear, auxiliary.
                   ! therefore the auxiliary moments are shiftet \Andreas Krapp
                  IF (stat_reorder_mom .AND. I .GT. nreg .AND. I .LE. nele) THEN
                     I = I + nnuc
                     IBTCH = IBTCH + nnuc
                  END IF
                   !
                  IF(A.NE.J_indices(I)%i_indx)   CALL LSQUIT('WRONG I_indx WHEN READING MM DER INFO FROM FILE',-1)
                  IF(B.NE.J_indices(I)%j_indx)   CALL LSQUIT('WRONG J_indx WHEN READING MM DER INFO FROM FILE',-1)
                  IF(IBTCH.NE.raw_paras(I)%batch)CALL LSQUIT('WRONG BATCH  WHEN READING MM DER INFO FROM FILE',-1)
                  IF(LMINA+LMINB.NE.raw_paras(I)%Lmin)&
                                              &  CALL LSQUIT('WRONG LMIN   WHEN READING MM DER INFO FROM FILE',-1)
                  IF(raw_paras(I)%cntr(1).NE.PX) CALL LSQUIT('WRONG CNTR-X WHEN READING MM DER INFO FROM FILE',-1)
                  IF(raw_paras(I)%cntr(2).NE.PY) CALL LSQUIT('WRONG CNTR-Y WHEN READING MM DER INFO FROM FILE',-1)
                  IF(raw_paras(I)%cntr(3).NE.PZ) CALL LSQUIT('WRONG CNTR-Z WHEN READING MM DER INFO FROM FILE',-1)
                  IF(EXTP.NE.raw_paras(I)%ext)   CALL LSQUIT('WRONG EXT    WHEN READING MM DER INFO FROM FILE',-1)
                  p = L*(L+1) +M +1
                   ! for the derivatives we need to know to which center and which direction the moments contribute
                   ! Ax = 1, Ay = 2, Az = 3, Bx = 4, By = 5, Bz = 6
                   ! we make the mapping of A and B onto the total list of atoms later 
                   ! using the mom2atom matrix
                  mom2atom(I,1) = ICENTA
                  mom2atom(I,2) = ICENTB
                  fac = 1.0D0 * gradfactor
                  if(icenta.eq.0.or.icentb.eq.0) fac = 2.0D0 * gradfactor
                  if(ISCOOR .ne. 0) raw_qlm_der(p,I,ISCOOR) = raw_qlm_der(p,I,ISCOOR) + SPHDER*fac
               END DO
            END DO readloopC

            ! close files
            CALL LSCLOSE(LUINTMI,'KEEP')
            CALL LSCLOSE(LUINTMR,'KEEP')
         END IF

         ! free memory
         !call mem_dealloc_fmm(RBUF)
         !call mem_dealloc_fmm(IBUF)
         DEALLOCATE (RBUF)
         DEALLOCATE (IBUF)

      ELSE
         ! this is the original routine, without using a buffer

         ! read electronic multipole moments into core

         LUINTM = 0
         CALL LSOPEN(LUINTM,INPUT_FILE1,'UNKNOWN','UNFORMATTED')
         REWIND (LUINTM)
         readloop3: DO
            READ (LUINTM) L,M, A,B, LMINA,LMINB,IBTCH,I,EXTP,PX,PY,PZ,SPH,ICENTA,ICENTB,ISCOOR
            IF (L .EQ. -1) EXIT readloop3
            iF (I > SIZE(raw_qlm,2))          CALL lsQUIT ('interface file error 1',-1)
            IF (((L+1)**2) > SIZE(raw_qlm,1)) CALL lsQUIT ('interface file error 2',-1)

             ! reorder data eventually:
             ! the moments on file are ordered as regular, auxiliary, nuclear
             ! we want the order regular, nuclear, auxiliary.
             ! therefore the auxiliary moments are shiftet \Andreas Krapp
            IF (stat_reorder_mom .AND. I .GT. nreg .AND. I .LE. nele) THEN
               I = I + nnuc
               IBTCH = IBTCH + nnuc
            END IF
             ! indices to map moments to orbitals in J-matrix
            J_indices(I)%i_indx = A
            J_indices(I)%j_indx = B
             ! index to batches of unique centre and extent
            raw_paras(I)%batch  = IBTCH
             ! minimum order for exact treatment (if uncontracted basis set)
            raw_paras(I)%Lmin = LMINA + LMINB
            raw_paras(I)%cntr = (/ PX, PY, PZ /)
            raw_paras(I)%ext  = EXTP
             ! see defn p.424 "Electronic Structure Theory", Helgaker et al
            p = L*(L+1) +M +1
             ! components (l,m) of MM and the MM of the derivatives wrt ISCOOR without density factorised in 
            IF(do_grad) SPH = SPH*gradfactor
            raw_qlm(p,I) = SPH
         END DO readloop3

         ! next read nuclei data: charge, and location
         n = n_mms%elec
         if (stat_reorder_mom) n = nreg
         DO I = 1, n_mms%nuc
            p = n + I
            READ (LUINTM) raw_qlm(1,p), raw_paras(p)%cntr(:)
            raw_qlm(1,p)=raw_qlm(1,p)
            raw_paras(p)%batch = raw_paras(n)%batch + I
            dens(p) = one
            IF(do_grad) THEN
               raw_qlm(1,p) = raw_qlm(1,p)*gradfactor
               raw_nucder=-raw_qlm(1,p)
               raw_qlm_der(1,p,1) = raw_nucder
               raw_qlm_der(1,p,2) = raw_nucder
               raw_qlm_der(1,p,3) = raw_nucder
               mom2atom(p,1)       = I
               mom2atom(p,2)       = I
               raw_paras(p)%Lmin   = 0     ! spherical charge
               ! the artificial extent of the nuclei is set to 1.0D0
               ! this has to match the settings in the routines 
               !    II_get_nn_gradient (LSint/IntegralInterface.f90) 
               !    getODcenter (LSint/ODbatch.f90) 
               raw_paras(p)%ext    = 1.0d0 ! point charges
               J_indices(p)%i_indx = 0     ! not relevant
               J_indices(p)%j_indx = 0     ! not relevant
            END IF
         END DO

         CALL LSCLOSE(LUINTM,'KEEP')

         ! read derivative electronic mutiplole moments into core
         IF(do_grad) THEN
            LUINTM = 0
            CALL LSOPEN(LUINTM,INPUT_FILE10,'UNKNOWN','UNFORMATTED')
            REWIND (LUINTM)
            readloop4: DO
               READ (LUINTM) L,M, A,B, LMINA,LMINB,IBTCH,I,EXTP,PX,PY,PZ,SPHDER,ICENTA,ICENTB,ISCOOR
               IF (L .EQ. -1) EXIT readloop4
               IF (I > SIZE(raw_qlm,2))          CALL lsQUIT ('interface file error 1',-1)
               IF (((L+1)**2) > SIZE(raw_qlm,1)) CALL lsQUIT ('interface file error 2',-1)

                ! reorder data eventually:
                ! the moments on file are ordered as regular, auxiliary, nuclear
                ! we want the order regular, nuclear, auxiliary.
                ! therefore the auxiliary moments are shiftet \Andreas Krapp
               IF (stat_reorder_mom .AND. I .GT. nreg .AND. I .LE. nele) THEN
                  I = I + nnuc
                  IBTCH = IBTCH + nnuc
               END IF
               IF(A .NE. J_indices(I)%i_indx)   CALL LSQUIT('WRONG I_indx WHEN READING MM DER INFO FROM FILE',-1)
               IF(B .NE. J_indices(I)%j_indx)   CALL LSQUIT('WRONG J_indx WHEN READING MM DER INFO FROM FILE',-1)
               IF(IBTCH .NE. raw_paras(I)%batch)CALL LSQUIT('WRONG BATCH  WHEN READING MM DER INFO FROM FILE',-1)
               IF(LMINA+LMINB.NE.raw_paras(I)%Lmin)&
                                            &   CALL LSQUIT('WRONG LMIN   WHEN READING MM DER INFO FROM FILE',-1)
               IF(raw_paras(I)%cntr(1) .NE. PX) CALL LSQUIT('WRONG CNTR-X WHEN READING MM DER INFO FROM FILE',-1)
               IF(raw_paras(I)%cntr(2) .NE. PY) CALL LSQUIT('WRONG CNTR-Y WHEN READING MM DER INFO FROM FILE',-1)
               IF(raw_paras(I)%cntr(3) .NE. PZ) CALL LSQUIT('WRONG CNTR-Z WHEN READING MM DER INFO FROM FILE',-1)
               IF(EXTP .NE. raw_paras(I)%ext)   CALL LSQUIT('WRONG EXT    WHEN READING MM DER INFO FROM FILE',-1)
               p = L*(L+1) +M +1
                ! for the derivatives we need to know to which center and which direction the moments contribute
                ! Ax = 1, Ay = 2, Az = 3, Bx = 4, By = 5, Bz = 6
                ! we make the mapping of A and B onto the total list of atoms later 
                ! using the mom2atom matrix
               mom2atom(I,1) = ICENTA
               mom2atom(I,2) = ICENTB
               fac = 1.0D0 * gradfactor
               IF( fit%AUX ) fac = 2.00D0 * gradfactor
               if (ISCOOR .ne.0) raw_qlm_der(p,I,ISCOOR) = raw_qlm_der(p,I,ISCOOR) + SPHDER*fac
            END DO readloop4
            CALL LSCLOSE(LUINTM,'KEEP')
         END IF
      END IF

      !
      ! read the density in
      !
      LUINTM = 0
      IF (fit%LHS_AUX .AND. .NOT. fit%RHS_AUX) THEN
         CALL LSOPEN(LUINTM,INPUT_FILE4,'UNKNOWN','UNFORMATTED')
      ELSE IF (fit%RHS_AUX) THEN
         CALL LSOPEN(LUINTM,INPUT_FILE5,'UNKNOWN','UNFORMATTED')
      ELSE IF (.NOT. fit%AUX) THEN
         CALL LSOPEN(LUINTM,INPUT_FILE3,'UNKNOWN','UNFORMATTED')
      ELSE
         CALL lsQUIT('what kind of density should be read in?',-1)
      END IF
      REWIND(LUINTM)
      READ (LUINTM) J
      IF (J .NE. n_mms%nbast) CALL lsQUIT('FMM-INTERFACE: Problem reading the density',-1)
      CALL ls_read(LUINTM,dens_in,J*(J+1)/2)
      CALL LSCLOSE(LUINTM,'KEEP')

      ! we read the "real" density, not the fitting density, so we can not take the n_mms%elec as such
      IEND = n_mms%elec - fit%NAUXMOM
      densloop: DO I = 1, IEND
         ! we saved only lower triangle, so we have to get the correct indices
          B=MIN(J_indices(I)%j_indx,J_indices(I)%i_indx)
          A=MAX(J_indices(I)%j_indx,J_indices(I)%i_indx)
         IBUFFER = (A-1)*(A-2)/2+(A-1) + B
         if(A.EQ.B) IBUFFER = A*(A-1)/2 + A
         dens(I) = dens_in(IBUFFER)
      END DO densloop

      !
      ! we read the fitting density if necessary
      !
      IF(fit%AUX) THEN
         LUINTM = 0
         IF (fit%LHS_AUX .AND. .NOT. fit%RHS_AUX) THEN
            CALL LSOPEN(LUINTM,INPUT_FILE6,'UNKNOWN','UNFORMATTED')
         ELSE IF (fit%RHS_AUX) THEN
            CALL LSOPEN(LUINTM,INPUT_FILE7,'UNKNOWN','UNFORMATTED')
         ELSE IF (fit%RHS_AUX .AND. fit%LHS_AUX) THEN
            CALL lsQUIT('RHSFIT .and. LHSFIT. mm_read_in_raw_data is confused!',-1)
         END IF
         REWIND(LUINTM)
         READ(LUINTM) J
         IF (J .NE. n_mms%nauxbas) CALL lsQUIT('FMM-INTERFACE: Problem reading the fitting-density',-1)
         CALL ls_read(LUINTM,densf_in,J)
         CALL LSCLOSE(LUINTM,'KEEP')
         ! we read the fitting density, so we can not take the n_mms%elec as such 
         ISTART = n_mms%elec - fit%NAUXMOM + 1
         IEND   = n_mms%elec
         IF(stat_reorder_mom) THEN
            ISTART = ISTART + n_mms%nuc
            IEND   = IEND + n_mms%nuc
         END IF
         densfloop: DO I = ISTART, IEND
            IF(fit%RHS_AUX) THEN
               IBUFFER = J_indices(I)%j_indx
            ELSE IF(fit%LHS_AUX .AND. .not. fit%RHS_AUX) THEN
               IBUFFER = J_indices(I)%j_indx
            END IF
           dens(I) = densf_in(IBUFFER)
         END DO densfloop
          dens((IEND+1):) = one
      END IF

  END SUBROUTINE mm_read_in_raw_data

!-------------------------------------------------------------------------------
! Routine to translate co-ordinates so all centres are at x,y,z >= 0
! Hence all box indices will be positive integers
!FIXME; is there an issue when we have ZERO SYSTEM SIZE???
   SUBROUTINE mm_translate_centres(system_size)

      IMPLICIT NONE
      REAL(REALK), INTENT(OUT) :: system_size 

      REAL(REALK)   :: sys_min(3), sys_max(3)
      INTEGER :: i

      ! determine system size
      sys_min = raw_paras(1)%cntr
      sys_max = raw_paras(1)%cntr
      DO i = 1, SIZE(raw_paras)
         sys_min(:) = MIN(sys_min(:),raw_paras(i)%cntr(:))
         sys_max(:) = MAX(sys_max(:),raw_paras(i)%cntr(:))
      END DO
      system_size = MAXVAL(sys_max - sys_min)

      DO i = 1, SIZE(raw_paras)
         raw_paras(i)%cntr = raw_paras(i)%cntr - sys_min
      END DO

   END SUBROUTINE mm_translate_centres
! this modified version centers the molecule in the boxes, so nothing is on the border, 
! this should reduce noise
!   SUBROUTINE mm_translate_centres(system_size,grain,do_center)
!
!      IMPLICIT NONE
!      REAL(REALK), INTENT(OUT) :: system_size
!      REAL(REALK), INTENT(IN)  :: grain
!
!      REAL(REALK)   :: sys_min(3), sys_max(3), rest(3)
!      INTEGER :: i, nbox(3)
!      LOGICAL       :: do_center

!      ! determine system size
!      sys_min = raw_paras(1)%cntr
!      sys_max = raw_paras(1)%cntr
!      DO i = 1, SIZE(raw_paras)
!         sys_min(:) = MIN(sys_min(:),raw_paras(i)%cntr(:))
!         sys_max(:) = MAX(sys_max(:),raw_paras(i)%cntr(:))
!      END DO
!      system_size = MAXVAL(sys_max - sys_min)

!      IF (do_center) THEN
!         ! center the system in the boxes, so that the border of the boxes is not in the atoms,
!         ! this should reduce noise in case of planar systems /Andreas Krapp 
!         DO i=1, 3
!            nbox(i) = ceiling((sys_max(i) - sys_min(i))/grain)
!            IF( nbox(i)*grain .eq. (sys_max(i) - sys_min(i))) THEN
!               ! if the systems length is exactly a multiple of the box sizes, we need one box more.
!               ! FIXME: some numerical issues here?
!               nbox(i) = nbox(i) + 1
!            END IF
!            rest(i) = nbox(i) * grain - (sys_max(i) - sys_min(i))
!         END DO
!         DO i = 1, SIZE(raw_paras)
!            raw_paras(i)%cntr = raw_paras(i)%cntr - sys_min + rest/2.0
!         END DO
!     ELSE
!         DO i = 1, SIZE(raw_paras)
!            raw_paras(i)%cntr = raw_paras(i)%cntr - sys_min 
!         END DO
!     END IF
!
!   END SUBROUTINE mm_translate_centres

!-------------------------------------------------------------------------------

   SUBROUTINE mm_check_moms(LMAX,qlm)

      IMPLICIT NONE

      REAL(REALK), INTENT(IN) :: qlm(:,:)
      INTEGER, INTENT(IN) :: LMAX
      INTEGER :: i,j, lm
      REAL(REALK) :: tmp

      write(lupri,*) "length of L,M array  =", SIZE(qlm,1)
      write(lupri,*) "length of mom array  =", SIZE(qlm,2)
      write(lupri,*) "length of dens array =", SIZE(dens)
      write(lupri,*) "number of electronic =", n_mms%elec

      IF ((1+LMAX)**2 /= SIZE(qlm,1)) CALL lsQUIT('error in mm_check_moms',-1)
      lm = 0
      DO i = 0, LMAX
         DO j = -i, i
            lm = lm+1
            tmp = (MAXVAL(ABS(qlm(lm,:))))
            IF ( tmp > zero ) write(*,'(2I4,2E25.15)')   &
                    i,j, MAXVAL(ABS(qlm(lm,:))), dens(lm)
         END DO
      END DO
      write(*,*)

!      DO i = 1, SIZE(qlm,2)
!         print *, i
!         print '(E25.15)', qlm(:,i)
!      END DO
!      write(*,*)

   END SUBROUTINE mm_check_moms

!-------------------------------------------------------------------------------

END MODULE mm_interface_mod

!===============================================================================

MODULE mm_qlm_processor

   USE mm_global_paras_mod
   USE mm_stats_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_renormalise_qlm,      &
             mm_factor_in_dens,       &
             mm_factor_in_dens_der,   &
             mm_get_T_sym_qlm,        &
             mm_pack_raw_moments,     &
             mm_pack_raw_parameters,  &
             mm_analyse_qlm_paras
            
CONTAINS

!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_renormalise_qlm(LMAX,qlm)

      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: LMAX
      REAL(REALK),   INTENT(INOUT) :: qlm(:,:)

      INTEGER :: i,L,M,p,pp
      REAL(REALK)   :: pref

      ! prefactor to symmetrize T-matrix
      DO i = 1, SIZE(qlm,2)
         DO L = 0, LMAX
            pp = L*(L+1) +1
            DO M = -L, -1
               pref = -one/(SQRT(two*FACTORIAL(l-m)*FACTORIAL(l+m)))
               p = pp+M
               qlm(p,i) = pref*qlm(p,i)
            END DO
            pref = one/FACTORIAL(L)
            p = pp  ! M=0
            qlm(p,i) = pref*qlm(p,i)
            DO M = 1, L
               pref = ((-1)**m)/SQRT(two*FACTORIAL(l-m)*FACTORIAL(l+m))
               p = pp+M
               qlm(p,i) = pref*qlm(p,i)
            END DO
         END DO
      END DO

   CONTAINS

      REAL(REALK) FUNCTION FACTORIAL(n)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: n 
         INTEGER :: i
         FACTORIAL = 1 
         DO i = n, 2, -1 
            FACTORIAL = FACTORIAL*i
         END DO
      END FUNCTION FACTORIAL

   END SUBROUTINE mm_renormalise_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE mm_factor_in_dens(dens,qlm)

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN)    :: dens(:)
      REAL(REALK),   INTENT(INOUT) :: qlm(:,:)

      INTEGER :: i

!FIXME: nicer way to do this with SPREAD??
      DO i = 1, SIZE(qlm,2)
         qlm(:,i) = qlm(:,i)*dens(i)
      END DO

   END SUBROUTINE mm_factor_in_dens

!-------------------------------------------------------------------------------

   SUBROUTINE mm_factor_in_dens_der(dens,qlm_der)

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN)    :: dens(:)
      REAL(REALK),   INTENT(INOUT) :: qlm_der(:,:,:)

      INTEGER :: i

!FIXME: nicer way to do this with SPREAD??
      DO i = 1, SIZE(qlm_der,2)
         qlm_der(:,i,:) = qlm_der(:,i,:)*dens(i)
      END DO

   END SUBROUTINE mm_factor_in_dens_der

!-------------------------------------------------------------------------------
! Prefactorising moments to symmetrize modified T-matrix

   SUBROUTINE mm_get_T_sym_qlm(LMAX,qlm_in,qlm_out)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: qlm_in(:,:)
      REAL(REALK),   INTENT(OUT) :: qlm_out(:,:)

      INTEGER :: i,L,u, hi,lo
      REAL(REALK)   :: pref

      DO i = 1, SIZE(qlm_in,2)
         DO L = 0, LMAX
            u = L*(L+1) +1     ! m=0
            hi = u+L
            lo = u-L
            pref = two*((-1)**L)
            qlm_out(lo:hi,i) = pref*qlm_in(lo:hi,i)
            qlm_out(u,i) = half*pref*qlm_in(u,i)
         END DO
      END DO

   END SUBROUTINE mm_get_T_sym_qlm

!-------------------------------------------------------------------------------
! get number of unique batches;
! also check raw_paras are sorted by batch (expected from input file);
! do this by checking the batch ID is always increasing

   SUBROUTINE get_nbatch(paras,nbatch)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)  :: paras(:)
      INTEGER,      INTENT(OUT) :: nbatch
      INTEGER :: i, ndim

      ndim = SIZE(paras)
      nbatch = 1
      DO i = 2, ndim
         IF ( paras(i)%batch < paras(i-1)%batch ) THEN
            CALL lsQUIT('batches input from source file not sorted!',-1)
         END IF
         ! some batch numbers are skipped so need to take care of this
         IF ( paras(i)%batch /= paras(i-1)%batch ) THEN
            nbatch = nbatch + 1
         END IF
      END DO

   END SUBROUTINE get_nbatch

!-------------------------------------------------------------------------------

   SUBROUTINE get_pkd_data(add_dens,raw_mms,pkd_paras,pkd_qlm)

      IMPLICIT NONE
      LOGICAL,            INTENT(IN)  :: add_dens
      TYPE(raw_mm_data),  INTENT(IN)  :: raw_mms
      TYPE(raw_mm_paras), INTENT(OUT) :: pkd_paras(:) 
      REAL(REALK),        INTENT(OUT) :: pkd_qlm(:,:)
      INTEGER :: i,j, last_batch

      j = 0
      last_batch = -1
      DO i = 1, SIZE(raw_mms%paras)
         IF ( raw_mms%paras(i)%batch == last_batch ) THEN
            ! element same batch as previous
            IF (add_dens) THEN
               pkd_qlm(:,j) = pkd_qlm(:,j) + raw_mms%qlm(:,i)*raw_mms%dens(i)
            ELSE
               pkd_qlm(:,j) = pkd_qlm(:,j) + raw_mms%qlm(:,i)
            END IF
            ! need to store most conservative Lmin (order of exact expansion)
            pkd_paras(j)%Lmin = MAX(pkd_paras(j)%Lmin, raw_mms%paras(i)%Lmin) 
         ELSE
            ! element in new batch
            j = j + 1
            pkd_paras(j) = raw_mms%paras(i)
            pkd_paras(j)%id = j
            IF (add_dens) THEN
               pkd_qlm(:,j) = raw_mms%qlm(:,i)*raw_mms%dens(i)
            ELSE
               pkd_qlm(:,j) = raw_mms%qlm(:,i)
            END IF
         END IF
         last_batch = raw_mms%paras(i)%batch
      END DO

   END SUBROUTINE get_pkd_data

!-------------------------------------------------------------------------------

   SUBROUTINE get_screened_nmom(screen_thr,qlm,skip,nskip)

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN)  :: screen_thr
      REAL(REALK),   INTENT(IN)  :: qlm(:,:)
      LOGICAL,       INTENT(OUT) :: skip(:)
      INTEGER, INTENT(OUT) :: nskip
      INTEGER :: i,j

      nskip = 0
      batches: DO i = 1, SIZE(qlm,2)
         skip(i) = .TRUE.
         nskip = nskip + 1
         lm_loop: DO j = 1, SIZE(qlm,1)
            IF ( ABS(qlm(j,i)) > screen_thr ) THEN
               skip(i) = .FALSE.
               nskip = nskip - 1
               EXIT lm_loop
            END IF
         END DO lm_loop
      END DO batches

   END SUBROUTINE get_screened_nmom

!-------------------------------------------------------------------------------
! here we squeeze all the significant batches of moments to be sequential
! at the top of the array, with all the insignificant moments overwritten;
! we use skip(:) as a logical mask to direct the skipping.

   SUBROUTINE squeeze_significant_batches(skip,paras,qlm)

      IMPLICIT NONE
      LOGICAL,            INTENT(IN)    :: skip(:)
      TYPE(raw_mm_paras), INTENT(INOUT) :: paras(:) 
      REAL(REALK),        INTENT(INOUT) :: qlm(:,:)
      INTEGER :: i, j

      IF (SIZE(paras) /= SIZE(qlm,2)) STOP 'paras and qlm should be same size!'
      IF (SIZE(paras) /= SIZE(skip)) STOP 'paras and skip should be same size!'

      j = 0
      DO i = 1, SIZE(paras)
         IF (skip(i)) CYCLE
         j = j + 1
         paras(j) = paras(i)
         qlm(:,j) = qlm(:,i)
      END DO

   END SUBROUTINE squeeze_significant_batches

!-------------------------------------------------------------------------------
! routine to drive the packing of a set of raw moments by batches,
! where members of a batch share a common centre and extent,
! including density factoring and screening if requested.
! note that no record is kept of the packing for later "unpacking".

   SUBROUTINE mm_pack_raw_moments(raw_mms,dens,dens_thr,pkd_paras,pkd_qlm,ndim)

      USE mm_memory_manager_mod, ONLY: mm_allocate, mm_deallocate

      IMPLICIT NONE
      TYPE(raw_mm_data),  INTENT(IN)  :: raw_mms
      LOGICAL,            INTENT(IN)  :: dens
      REAL(REALK),        INTENT(IN)  :: dens_thr
      TYPE(raw_mm_paras), POINTER     :: pkd_paras(:) 
      REAL(REALK),        POINTER     :: pkd_qlm(:,:)
      INTEGER,      INTENT(OUT) :: ndim
      INTEGER :: nbatch, nskip
      LOGICAL,POINTER :: skip(:)

      ! get number of unique batches;
      CALL get_nbatch(raw_mms%paras,nbatch)

      ! now build packed data by summing raw data in the same batch
      CALL mm_allocate(MEM_RAW_QLM,pkd_paras,nbatch)
      CALL mm_allocate(MEM_RAW_QLM,pkd_qlm,INT(SIZE(raw_mms%qlm,1)),nbatch)
      CALL get_pkd_data(dens,raw_mms,pkd_paras,pkd_qlm)
      ndim = nbatch

      IF ( dens .AND. (dens_thr > DENS_SCREEN_CUT) ) THEN
         !WRITE(LUPRI,'(A,E9.1)') ' MM RHS density screening on:', dens_thr
         ! perform density-based screening;
         ! note we do not reallocate whole array, but push all the
         ! significant terms to the top, and only this array section
         ! of significant moments is then pointed to.
!         ALLOCATE(skip(nbatch))
         call mem_alloc_fmm(skip,nbatch)
         CALL get_screened_nmom(dens_thr,pkd_qlm,skip,nskip)
         ndim = nbatch - nskip  ! number of significant moments
         CALL squeeze_significant_batches(skip,pkd_paras,pkd_qlm)
         call mem_dealloc_fmm(skip)
!         DEALLOCATE(skip)
      END IF

      ! this is just for statistics
      stat_pkd_moms_RHS = nbatch
      stat_screened_moms_RHS = ndim

   END SUBROUTINE mm_pack_raw_moments

!-------------------------------------------------------------------------------
! routine to drive the packing of a set of raw mm parameters by batches,
! where members of a batch share a common centre and extent,
! including the build of a mapping between the packed paras and the
! original raw paras.  This map can be used for later "unpacking".

   SUBROUTINE mm_pack_raw_parameters(raw_paras,pkd_paras,batch_map)

      USE mm_memory_manager_mod, ONLY: mm_allocate, mm_deallocate

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN) :: raw_paras(:) 
      TYPE(raw_mm_paras), POINTER    :: pkd_paras(:) 
      TYPE(id_list),      POINTER    :: batch_map(:) 
      INTEGER :: i,j, nbatch, last_batch

      ! get number of unique batches;
      CALL get_nbatch(raw_paras,nbatch)
      stat_pkd_moms_LHS = nbatch

      ! initialise packed paras and batch map
      CALL mm_allocate(MEM_RAW_QLM,pkd_paras,nbatch)
      ALLOCATE(batch_map(nbatch))
      DO i = 1, nbatch 
         batch_map(i)%occ = 0
         NULLIFY( batch_map(i)%head )
      END DO

      ! now build packed paras by compressing raw paras in same batch
      j = 0
      last_batch = -1
      DO i = 1, SIZE(raw_paras)
         IF ( raw_paras(i)%batch == last_batch ) THEN
            ! add raw parameter mapping to existing linked-list for this batch
            CALL add_batch_item(batch_map(j),i)
            ! need to store most conservative Lmin (order of exact expansion)
            pkd_paras(j)%Lmin = MAX(pkd_paras(j)%Lmin, raw_paras(i)%Lmin) 
         ELSE
            ! element in new batch
            j = j + 1
            pkd_paras(j) = raw_paras(i)
            pkd_paras(j)%id = j
            ! linked-list for this batch is empty, so start one
            batch_map(j)%occ = 1
            ALLOCATE( batch_map(j)%head )
            batch_map(j)%head%id = i           ! map back to raw paras here
            NULLIFY( batch_map(j)%head%next )  ! rest of list is empty
         END IF
         last_batch = raw_paras(i)%batch
      END DO

   CONTAINS

      SUBROUTINE add_batch_item(batch_list,raw_id)

         IMPLICIT NONE
         TYPE(id_list), INTENT(INOUT) :: batch_list
         INTEGER, INTENT(IN)    :: raw_id
         TYPE(id_node), POINTER :: new_node

         batch_list%occ = batch_list%occ + 1
         ALLOCATE( new_node )
         new_node%id = raw_id
         IF (ASSOCIATED(batch_list%head%next)) THEN
            ! more than one entry in list (including head)
            ! so point new_node to old second entry
            new_node%next => batch_list%head%next
            ! point head to new_node
            NULLIFY(batch_list%head%next)
            batch_list%head%next => new_node
         ELSE
            ! only head so far; make new_node our second entry
            batch_list%head%next => new_node
            NULLIFY(new_node%next)   ! end of list
         END IF

      END SUBROUTINE add_batch_item

   END SUBROUTINE mm_pack_raw_parameters

!-------------------------------------------------------------------------------

   SUBROUTINE mm_analyse_qlm_paras(paras)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN) :: paras(:) 
      REAL(REALK)   :: centre(3)
      INTEGER :: i,j, ncen

      ncen = 1
      DO i = 2, SIZE(paras) 
         ncen = ncen +1
         centre = paras(i)%cntr
         jloop: DO j = i-1, 1, -1
            IF ( centre(1) == paras(j)%cntr(1) ) THEN
               IF ( centre(2) == paras(j)%cntr(2) ) THEN
                  IF ( centre(3) == paras(j)%cntr(3) ) THEN
                     ncen = ncen -1
                     EXIT jloop
                  END IF
               END IF
            END IF
         END DO jloop
      END DO
      WRITE(LUPRI,'(A,I8)') ' distinct centres:', ncen

   END SUBROUTINE mm_analyse_qlm_paras

!-------------------------------------------------------------------------------

END MODULE mm_qlm_processor

!===============================================================================

MODULE mm_aux_qlm_builder

   USE mm_global_paras_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_qlm_builder,    &
             mm_free_qlm_builder,    &
             mm_get_aux_qlm
            
   ! we use pointers here so we can allocate outside of the module
   TYPE(raw_mm_paras), POINTER, SAVE :: LHS_paras(:) 
   TYPE(id_list),      POINTER, SAVE :: LHS_batch_map(:) 
   TYPE(raw_mm_paras), POINTER, SAVE :: RHS_paras(:) 
   REAL(REALK),        POINTER, SAVE :: RHS_qlm_T(:,:)
   REAL(REALK),        POINTER, SAVE :: RHS_qlm_W(:,:)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_qlm_builder

      IMPLICIT NONE

      NULLIFY(LHS_paras, LHS_batch_map)
      NULLIFY(RHS_paras, RHS_qlm_T, RHS_qlm_W)

   END SUBROUTINE mm_init_qlm_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_qlm_builder
     USE mm_memory_manager_mod, ONLY: mm_deallocate

      IMPLICIT NONE
      INTEGER :: i

      IF (ASSOCIATED(LHS_batch_map)) THEN
         DO i = 1, SIZE(LHS_batch_map)
            CALL free_batch_map(LHS_batch_map(i)%head)
         END DO
      END IF

      IF (ASSOCIATED(LHS_paras)) THEN
         call mm_deallocate(LHS_paras)
         DEALLOCATE(LHS_batch_map)
      ENDIF
      IF (ASSOCIATED(RHS_paras))THEN
         call mm_deallocate(RHS_paras)
      ENDIF
      IF (ASSOCIATED(RHS_qlm_T))THEN
         call mm_deallocate(RHS_qlm_T)
      ENDIF
      IF (ASSOCIATED(RHS_qlm_W))THEN
         call mm_deallocate(RHS_qlm_W)
      ENDIF
      NULLIFY(LHS_paras, LHS_batch_map)
      NULLIFY(RHS_paras, RHS_qlm_T, RHS_qlm_W)

   CONTAINS 

      RECURSIVE SUBROUTINE free_batch_map(node)

         IMPLICIT NONE
         TYPE(id_node), POINTER :: node

         IF (ASSOCIATED(node%next)) THEN
            CALL free_batch_map(node%next)
         END IF
         DEALLOCATE(node)
         NULLIFY(node)

      END SUBROUTINE free_batch_map

   END SUBROUTINE mm_free_qlm_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_aux_qlm(scheme,LHS_mms,RHS_mms)
      USE mm_stats_mod

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: LHS_mms, RHS_mms

       ! we use the SCALED solid harmonic formulation throughout
      IF (USE_UNSCALED_HARMONICS) CALL lsQUIT ('unscaled harmonics not coded!',-1)
      IF(stat_do_grad) THEN
         CALL get_normalised_qlm_der(scheme,LHS_mms%qlm,RHS_mms%qlm,LHS_mms%qlm_der,RHS_mms%qlm_der)
      ELSE
         CALL get_normalised_qlm(scheme,LHS_mms%qlm,RHS_mms%qlm)
      ENDIF

       ! get RHS parameter mappings and preconditioned W and T matrices first
      CALL get_RHS_data(scheme,RHS_mms)
       ! get LHS parameter mappings and preconditioned W and T matrices next;
       ! note we get LHS after RHS because LHS moments may
       ! alter the raw moments subsequently used by RHS
      CALL get_LHS_data(scheme,LHS_mms)

      ! all relevant data should now be held in %qlm_T and %qlm_W
      NULLIFY (LHS_mms%qlm, RHS_mms%qlm)
      IF(stat_do_grad) THEN
         NULLIFY (LHS_mms%qlm_der, RHS_mms%qlm_der)
      ENDIF

   END SUBROUTINE mm_get_aux_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE get_normalised_qlm(scheme,LHS_qlm,RHS_qlm)
      USE density_fitting
      USE mm_qlm_processor, ONLY: mm_renormalise_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(INOUT) :: LHS_qlm(:,:), RHS_qlm(:,:)
      INTEGER :: LMAX, i, qlm_dim

      LMAX = scheme%raw_LMAX

      IF (scheme%LHS_mm_range == ALL_MOMENTS) THEN
         ! all data renormalised via LHS
         CALL mm_renormalise_qlm(LMAX,LHS_qlm) 
      ELSE IF (scheme%RHS_mm_range == ALL_MOMENTS) THEN
         ! all data renormalised via RHS
         CALL mm_renormalise_qlm(LMAX,RHS_qlm) 
      ELSE
         CALL mm_renormalise_qlm(LMAX,LHS_qlm) 
!CAREFUL affects the results if renormalized twice!!!!!!!!!!!
!        IF ((scheme%RHS_mm_range /= scheme%LHS_mm_range) .OR. fit%AUX) THEN
         IF ((scheme%RHS_mm_range /= scheme%LHS_mm_range) .OR.  &
              & (fit%LHS_AUX .AND. .NOT. fit%RHS_AUX)     .OR.  & 
              & (fit%RHS_AUX .AND. .NOT. fit%LHS_AUX) ) THEN
!         IF (scheme%RHS_mm_range /= scheme%LHS_mm_range) THEN 
            ! LHS and RHS are disjoint (must avoid scaling some twice!)
            CALL mm_renormalise_qlm(LMAX,RHS_qlm)
         END IF
      END IF

   END SUBROUTINE get_normalised_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE get_normalised_qlm_der(scheme,LHS_qlm,RHS_qlm,LHS_qlm_der,RHS_qlm_der)
      USE density_fitting
      USE mm_stats_mod
      USE mm_qlm_processor, ONLY: mm_renormalise_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK), INTENT(INOUT) :: LHS_qlm_der(:,:,:), RHS_qlm_der(:,:,:)
      REAL(REALK), INTENT(INOUT) :: LHS_qlm(:,:), RHS_qlm(:,:)
      INTEGER :: LMAX, i, qlm_dim, iscoor

      LMAX = scheme%raw_LMAX

      IF (scheme%LHS_mm_range == ELECTRONIC_ONLY .AND. &
       &  scheme%RHS_mm_range == ELECTRONIC_ONLY ) THEN

         ! normalise moments via  RHS
         CALL mm_renormalise_qlm(LMAX,RHS_qlm)

         ! normalize derivative moments via LHS
         DO iscoor = 1, 6
            CALL mm_renormalise_qlm(LMAX,LHS_qlm_der(:,:,iscoor))
         END DO

      ELSE
         ! this is specially needed for the nuclear attraction 
         ! (nuc'|nuc+e)

         ! renormalise RHS
         CALL mm_renormalise_qlm(LMAX,RHS_qlm)

         ! renormalise LHS (derivatives)
         DO iscoor = 1, 6
            CALL mm_renormalise_qlm(LMAX,LHS_qlm_der(:,:,iscoor))
         END DO

      END IF

   END SUBROUTINE get_normalised_qlm_der


!-------------------------------------------------------------------------------

   SUBROUTINE get_LHS_data(scheme,LHS_mms)

      USE mm_memory_manager_mod, ONLY: mm_allocate, mm_deallocate
      USE mm_qlm_processor,  ONLY: mm_factor_in_dens,       &
                                   mm_pack_raw_parameters,  & 
                                   mm_factor_in_dens_der
      USE density_fitting
!ANDREAS
      USE mm_stats_mod
!

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: LHS_mms
      INTEGER :: i, qlm_dim


       ! get LHS parameters from raw data
      IF (scheme%pack_LHS) THEN
         CALL mm_pack_raw_parameters(LHS_mms%paras,LHS_paras,LHS_batch_map)
         qlm_dim = SIZE(LHS_paras)
         LHS_mms%batch_map => LHS_batch_map(:)
      ELSE
         qlm_dim = SIZE(LHS_mms%qlm,2)
         CALL mm_allocate(MEM_RAW_QLM,LHS_paras,qlm_dim)
         LHS_paras = LHS_mms%paras(:)
      END IF

       ! set up (packed) LHS parameters and paras:moment mappings
      ! the compact expression is:
      ! LHS_paras(:)%id = (/ (i,i=1,qlm_dim) /)
      ! but intel fortran compiler 7.1-generated code may deadlock
      ! for certain vector lengths (in particular, for insulin molecule).
      ! let's play safe here.
      DO i = 1, qlm_dim
         LHS_paras(i)%id = i
      END DO
      LHS_mms%paras => LHS_paras(:)

       ! LHS T-matrix
      SELECT CASE (scheme%T_con%LHS_mm_type)
      CASE (USE_RAW_QLM)
         LHS_mms%qlm_T => LHS_mms%qlm(:,:)
         IF(stat_do_grad) THEN 
            LHS_mms%qlm_der_T => LHS_mms%qlm_der(:,:,:)
         END IF
      CASE DEFAULT
         CALL lsQUIT('cannot reconcile LHS_mm_type',-1) 
      END SELECT

       ! factorise in density if required
      IF ( scheme%LHS_dens ) THEN
         CALL mm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         IF (stat_do_grad) THEN 
            CALL mm_factor_in_dens_der(LHS_mms%dens,LHS_mms%qlm_der_T)
         END IF
         NULLIFY(LHS_mms%dens)  ! should no-longer be used
      END IF

       ! LHS W-matrix
      NULLIFY(LHS_mms%qlm_W) ! not actually used

   END SUBROUTINE get_LHS_data


!-------------------------------------------------------------------------------

   SUBROUTINE get_RHS_data(scheme,RHS_mms)

      USE mm_memory_manager_mod, ONLY: mm_allocate, mm_deallocate
      USE mm_qlm_processor,  ONLY: mm_get_T_sym_qlm,     &
                                   mm_pack_raw_moments,  &
                                   mm_factor_in_dens,    &
                                   mm_analyse_qlm_paras
      USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: RHS_mms
      INTEGER :: LMAX, i, qlm_dim
      LOGICAL :: dens ! flags if density is factored in
      REAL(REALK) :: thr

      LMAX = scheme%raw_LMAX

       ! get RHS parameters and moments from raw data
      IF ( scheme%pack_RHS ) THEN
         ! we now pack all moments belonging to
         ! the same batch (same centre,extent);
         ! also factor in density and perform
         ! density-based screening if flagged
         dens = scheme%RHS_dens 
         thr = scheme%dens_screen_thr
         CALL mm_pack_raw_moments(RHS_mms,dens,thr,RHS_paras,RHS_qlm_W,qlm_dim)
      ELSE
         ! use unpacked moments and parameters (as above for LHS)
         qlm_dim = SIZE(RHS_mms%qlm,2)
         CALL mm_allocate(MEM_RAW_QLM,RHS_paras,qlm_dim)
         RHS_paras = RHS_mms%paras(:)
         ! RHS W-matrix
         CALL mm_allocate(MEM_RAW_QLM,RHS_qlm_W,(1+LMAX)**2,qlm_dim)
         RHS_qlm_W = RHS_mms%qlm(:,:)
         IF ( scheme%RHS_dens ) THEN
            CALL mm_factor_in_dens(RHS_mms%dens,RHS_qlm_W)
            NULLIFY(RHS_mms%dens)  ! should no-longer be used
         END IF
      END IF

      ! note we must take care to only point to the required
      ! array section since screening may mean that the array is
      ! larger than the number of significant moments we want 
      
       ! set up (packed) RHS parameters and paras:moment mappings
      DO i = 1, qlm_dim
         RHS_paras(i)%id = i
      END DO
      RHS_mms%paras => RHS_paras(:qlm_dim)
      RHS_mms%qlm_W => RHS_qlm_W(:,:qlm_dim)

       ! precondition RHS T-matrix
      SELECT CASE (scheme%T_con%RHS_mm_type)
      CASE (USE_RAW_QLM)
         RHS_mms%qlm_T => RHS_qlm_W(:,:qlm_dim)
      CASE (USE_T_SYM_QLM)
         CALL mm_allocate(MEM_RAW_QLM,RHS_qlm_T,(1+LMAX)**2,qlm_dim)
          ! build %qlm_T by rescaling significant %qlm_W
         CALL mm_get_T_sym_qlm(LMAX,RHS_mms%qlm_W,RHS_qlm_T)
         RHS_mms%qlm_T => RHS_qlm_T(:,:)
      CASE DEFAULT
         CALL lsQUIT('cannot reconcile RHS_mm_type',-1) 
      END SELECT

!      CALL mm_analyse_qlm_paras(RHS_mms%paras)

   END SUBROUTINE get_RHS_data

!-------------------------------------------------------------------------------

END MODULE mm_aux_qlm_builder

! Globally accessible subroutine that sets threshold for
! screening away contributions in T_contraction
SUBROUTINE mm_set_screen_threshold(threshold)
use mm_global_paras_mod
use mm_T_worker
use mm_T_contractors
implicit none
REAL(realk) :: threshold
  call mm_set_T_contractor_threshold(threshold)
  call mm_set_T_worker_threshold(threshold)
END SUBROUTINE mm_set_screen_threshold
   
