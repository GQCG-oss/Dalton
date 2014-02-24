MODULE mm_Vff_processor

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_get_E_from_Vff,        &
             mm_get_E_from_pkd_Vff,    &
             mm_get_grad_from_Vff,     &
             mm_get_grad_from_pkd_Vff, &
             mm_get_J_from_Vff,        &
             mm_get_J_from_pkd_Vff

CONTAINS

!-------------------------------------------------------------------------------
! Trivial check that the number of moments matches the number of potentials

   SUBROUTINE mm_verify_Vff_input(scheme,LHS_mms,Vff,J_or_E)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      TYPE(raw_mm_data),  INTENT(IN) :: LHS_mms
      REAL(REALK),        INTENT(IN) :: Vff(:,:) 
      CHARACTER(1),       INTENT(IN) :: J_or_E
      LOGICAL :: A,B,C

      A = (SIZE(LHS_mms%paras) /= SIZE(Vff,2))
      IF (A) CALL LSQUIT('incompatible SIZE of Vff and LHS moments!',-1)

      IF (J_or_E == 'J') THEN
         A = (scheme%LHS_mm_range == NUCLEAR_ONLY)
         B = (scheme%RHS_mm_range == NUCLEAR_ONLY)
         C = (scheme%LHS_mm_range == ALL_MOMENTS )
         IF ( (A.AND.B) .OR. C ) CALL LSQUIT('mm_ranges invalid for J_matrix!',-1)
      END IF

   END SUBROUTINE mm_verify_Vff_input

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_E_with_text(scheme,energy,text)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(INOUT) :: energy
      CHARACTER(*),       INTENT(OUT)   :: text
      LOGICAL :: A,B,C,D

      A = (scheme%LHS_mm_range == ELECTRONIC_ONLY)
      B = (scheme%RHS_mm_range == ELECTRONIC_ONLY)
      C = (scheme%LHS_mm_range == NUCLEAR_ONLY)
      D = (scheme%RHS_mm_range == NUCLEAR_ONLY)

      IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
         energy = half*energy
         text = "total classical Coulomb energy"
         IF (A) text = "classical Coulomb electronic energy"
         IF (C) text = "classical Coulomb nuclear repulsion"
      ELSE IF (A.OR.B) THEN
         IF (C.OR.D) THEN
            text = "classical Coulomb nuclear attraction"
         ELSE
            text = "e-n + 2*(e-e) energy"
         END IF
      ELSE
         ! range is different and neither is ELECTRONIC_ONLY
         text = "e-n + 2*(n-n) energy"
      END IF

   END SUBROUTINE mm_get_E_with_text

!-------------------------------------------------------------------------------
! Get energy from full contraction of LHS moments and the potentials
!-------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! (so that a simple contraction with the LHS moments will double count
!  interactions if the LHS and RHS moment ranges overlap).

   SUBROUTINE mm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,text)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(IN)  :: LHS_mms
      REAL(REALK),        INTENT(IN)  :: Vff(:,:) 
      REAL(REALK),        INTENT(OUT) :: energy
      CHARACTER(*),       INTENT(OUT) :: text

      REAL(REALK)   :: g 
      INTEGER :: u,v, lm_max
      energy = 0.0E0_realk 
      CALL mm_verify_Vff_input(scheme,LHS_mms,Vff,'E')

      IF (scheme%dynamic_LMAX_on) THEN
         DO u = 1, SIZE(LHS_mms%paras)
            ! lm_max should never exceed raw_LMAX, and the dims
            ! of Vff also never exceed (1+raw_LMAX)**2
            lm_max = (1 + LHS_mms%paras(u)%Lmin + scheme%LEXTRA)**2
            lm_max = MIN(lm_max, INT(SIZE(Vff,1)) )
            v = LHS_mms%paras(u)%id
            g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
            energy = energy + g
         END DO
      ELSE
         ! although Vff should be the same size, we test for generality
         IF (SIZE(LHS_mms%qlm_T,1) /= SIZE(Vff,1)) STOP 'mm_get_E_from_Vff:2'
         lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
         DO u = 1, SIZE(LHS_mms%paras)
            v = LHS_mms%paras(u)%id
            g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
            energy = energy + g
         END DO
      END IF

      CALL mm_get_E_with_text(scheme,energy,text)

   END SUBROUTINE mm_get_E_from_Vff

!-------------------------------------------------------------------------------
! Get energy from full contraction of LHS moments and the potentials
! on the basis that the LHS parameters and Vff may have been packed
! into batches (and must be "expanded").
!-------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! (so that a simple contraction with the LHS moments will double count
!  interactions if the LHS and RHS moment ranges overlap).

   SUBROUTINE mm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,text)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(IN)  :: LHS_mms
      REAL(REALK),        INTENT(IN)  :: Vff(:,:) 
      REAL(REALK),        INTENT(OUT) :: energy
      CHARACTER(*),       INTENT(OUT) :: text

      REAL(REALK)   :: g 
      INTEGER :: i,u,v,w, lm_max
      TYPE(id_node), POINTER :: batch_map

      CALL mm_verify_Vff_input(scheme,LHS_mms,Vff,'E')

      IF (scheme%dynamic_LMAX_on) THEN

         packed_loop: DO u = 1, SIZE(LHS_mms%paras)

            ! lm_max should never exceed raw_LMAX, and the dims
            ! of Vff also never exceed (1+raw_LMAX)**2
            lm_max = (1 + LHS_mms%paras(u)%Lmin + scheme%LEXTRA)**2
            lm_max = MIN(lm_max, INT(SIZE(Vff,1)) )
            v = LHS_mms%paras(u)%id  ! LHS packed (batch) moments
            batch_map => LHS_mms%batch_map(v)%head

            batch_members: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
               energy = energy + g
                ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members
               batch_map => batch_map%next
            END DO batch_members 

         END DO packed_loop

      ELSE

         lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
         packed_loop2: DO u = 1, SIZE(LHS_mms%paras)

            v = LHS_mms%paras(u)%id  ! LHS packed (batch) parameters 
            batch_map => LHS_mms%batch_map(v)%head

            batch_members2: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
               energy = energy + g
                ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members2
               batch_map => batch_map%next
            END DO batch_members2 

         END DO packed_loop2

      END IF

      CALL mm_get_E_with_text(scheme,energy,text)

   END SUBROUTINE mm_get_E_from_pkd_Vff

!--------------------------------------------------------------------------------
! Get gradient from full contraction of LHS moment derivatives and the potentials
!--------------------------------------------------------------------------------

   SUBROUTINE mm_get_grad_from_Vff(scheme,LHS_mms,Vff,grad,energy,numnuc)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(IN)    :: LHS_mms
      REAL(REALK),        INTENT(IN)    :: Vff(:,:)
      REAL(REALK),        INTENT(INOUT) :: grad(3*numnuc),energy
      INTEGER,      INTENT(IN)    :: numnuc

      REAL(REALK)   :: g, fac
      INTEGER :: u,v, lm_max, iscoor, icenta,icentb,iscool,v2
      INTEGER :: ix, ixb
      LOGICAL :: doa, dob

      CALL mm_verify_Vff_input(scheme,LHS_mms,Vff,'E')
   IF(scheme%LHS_mm_range ==  NUCLEAR_ONLY ) THEN
      IF (scheme%dynamic_LMAX_on) THEN
         CALL lsquit('dynamic_LMAX not with NUCLEAR_ONLY at LHS for mm_get_grad_from_Vff',-1)
      ELSE
         IF (SIZE(LHS_mms%qlm_der_T,1) /= SIZE(Vff,1)) STOP 'mm_get_grad_from_Vff:2'
         lm_max = MIN(SIZE(LHS_mms%qlm_der_T,1),SIZE(Vff,1))
         DO u = 1, SIZE(LHS_mms%paras)
            v = LHS_mms%paras(u)%id
            icenta=LHS_mms%mom2atom(v,1)
            ! nuclear contribution
            do ix = 1, 3
               if (icenta .gt. 0) then
                  iscoor = (icenta-1)*3 + ix
                  fac = 0.5D0
                  if (ix == 1) ixb = 4
                  if (ix == 2) ixb = 2
                  if (ix == 3)then
                     ixb = 3
                     fac = -1.0D0
                  end if
                  grad(iscoor) = grad(iscoor) + LHS_mms%qlm_der_T(1,v,ix) * Vff(ixb,v)*fac
               end if
             end do
         END DO
      END IF
   ELSEIF(scheme%LHS_mm_range ==  ALL_MOMENTS ) THEN
      IF (scheme%dynamic_LMAX_on) THEN
         CALL lsquit('dynamic_LMAX not with ALL_MOMENTS at LHS for mm_get_grad_from_Vff',-1)
      ELSE
         ! although Vff should be the same size, we test for generality
         IF (SIZE(LHS_mms%qlm_der_T,1) /= SIZE(Vff,1)) STOP 'mm_get_grad_from_Vff:2'
         lm_max = MIN(SIZE(LHS_mms%qlm_der_T,1),SIZE(Vff,1))
         DO u = 1, SIZE(LHS_mms%paras)
            v = LHS_mms%paras(u)%id
            icenta=LHS_mms%mom2atom(v,1)
            IF (v .GE. LHS_mms%startnuc .AND. v .LE. LHS_mms%endnuc) THEN
              ! nuclear contribution
               do ix = 1, 3
                  if (icenta .gt. 0) then
                     iscoor = (icenta-1)*3 + ix
                     fac = 0.5D0
                     if (ix == 1) ixb = 4
                     if (ix == 2) ixb = 2
                     if (ix == 3)then
                        ixb = 3
                        fac = -1.0D0
                     end if
                     grad(iscoor) = grad(iscoor) + LHS_mms%qlm_der_T(1,v,ix) * Vff(ixb,v)*fac
                  end if
               end do
            ELSE
               icentb=LHS_mms%mom2atom(v,2)
               doa = .false.
               dob = .false.
               if(icenta .gt. 0 ) doa = .true.
               if(icentb .gt. 0 ) dob = .true.
               do iscool = 1, 6
                  iscoor = 0
                  if (iscool .eq. 1 .and. doa) iscoor = (icenta-1)*3 + 1
                  if (iscool .eq. 2 .and. doa) iscoor = (icenta-1)*3 + 2
                  if (iscool .eq. 3 .and. doa) iscoor = (icenta-1)*3 + 3
                  if (iscool .eq. 4 .and. dob) iscoor = (icentb-1)*3 + 1
                  if (iscool .eq. 5 .and. dob) iscoor = (icentb-1)*3 + 2
                  if (iscool .eq. 6 .and. dob) iscoor = (icentb-1)*3 + 3
                  if (iscoor .gt. 0 ) &
                   & grad(iscoor) = grad(iscoor) + DOT_PRODUCT(LHS_mms%qlm_der_T(:lm_max,v,iscool),Vff(:lm_max,v))
               end do
               energy = energy + DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
            END IF
         END DO
      END IF
   ELSE
      IF (scheme%dynamic_LMAX_on) THEN
         DO u = 1, SIZE(LHS_mms%paras)
            ! lm_max should never exceed raw_LMAX, and the dims
            ! of Vff also never exceed (1+raw_LMAX)**2
            lm_max = (1 + LHS_mms%paras(u)%Lmin + scheme%LEXTRA)**2
            lm_max = MIN(lm_max, INT(SIZE(Vff,1)) )
            v = LHS_mms%paras(u)%id
            icenta=LHS_mms%mom2atom(v,1)
            icentb=LHS_mms%mom2atom(v,2)
            doa = .false.
            dob = .false.
            if(icenta .gt. 0 ) doa = .true.
            if(icentb .gt. 0 ) dob = .true.
            do iscool = 1, 6
               iscoor = 0
               if (iscool .eq. 1 .and. doa) iscoor = (icenta-1)*3 + 1
               if (iscool .eq. 2 .and. doa) iscoor = (icenta-1)*3 + 2
               if (iscool .eq. 3 .and. doa) iscoor = (icenta-1)*3 + 3
               if (iscool .eq. 4 .and. dob) iscoor = (icentb-1)*3 + 1
               if (iscool .eq. 5 .and. dob) iscoor = (icentb-1)*3 + 2
               if (iscool .eq. 6 .and. dob) iscoor = (icentb-1)*3 + 3
               if (iscoor .gt. 0) &
                   & grad(iscoor) = grad(iscoor) + DOT_PRODUCT(LHS_mms%qlm_der_T(:lm_max,v,iscool),Vff(:lm_max,v))
            end do
            energy = energy + DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
         END DO
      ELSE
         ! although Vff should be the same size, we test for generality
         IF (SIZE(LHS_mms%qlm_der_T,1) /= SIZE(Vff,1)) STOP 'mm_get_grad_from_Vff:2'
         lm_max = MIN(SIZE(LHS_mms%qlm_der_T,1),SIZE(Vff,1))
         DO u = 1, SIZE(LHS_mms%paras)
            v = LHS_mms%paras(u)%id
            icenta=LHS_mms%mom2atom(v,1)
            icentb=LHS_mms%mom2atom(v,2)
            doa = .false.
            dob = .false.
            if(icenta .gt. 0 ) doa = .true.
            if(icentb .gt. 0 ) dob = .true.
            do iscool = 1, 6
               iscoor = 0
               if (iscool .eq. 1 .and. doa) iscoor = (icenta-1)*3 + 1
               if (iscool .eq. 2 .and. doa) iscoor = (icenta-1)*3 + 2
               if (iscool .eq. 3 .and. doa) iscoor = (icenta-1)*3 + 3
               if (iscool .eq. 4 .and. dob) iscoor = (icentb-1)*3 + 1
               if (iscool .eq. 5 .and. dob) iscoor = (icentb-1)*3 + 2
               if (iscool .eq. 6 .and. dob) iscoor = (icentb-1)*3 + 3
               if (iscoor .gt. 0 ) &
                   & grad(iscoor) = grad(iscoor) + DOT_PRODUCT(LHS_mms%qlm_der_T(:lm_max,v,iscool),Vff(:lm_max,v))
            end do
            energy = energy + DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
         END DO
      END IF
   END IF

   END SUBROUTINE mm_get_grad_from_Vff

!----------------------------------------------------------------------------------
! Get gradient from full contraction of LHS moments derivatives and the potentials
!----------------------------------------------------------------------------------

   SUBROUTINE mm_get_grad_from_pkd_Vff(scheme,LHS_mms,Vff,grad,energy,numnuc)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(IN)  :: LHS_mms
      REAL(REALK),        INTENT(IN)  :: Vff(:,:)
      REAL(REALK),        INTENT(INOUT) :: grad(numnuc*3),energy
      INTEGER,      INTENT(IN)  :: numnuc 

      REAL(REALK)   :: g,fac
      INTEGER :: i,u,v,w, lm_max, iscoor, il, icenta, icentb, iscool,ix,ixb, w2
      TYPE(id_node), POINTER :: batch_map
      LOGICAL :: doa, dob

      CALL mm_verify_Vff_input(scheme,LHS_mms,Vff,'E')

   IF(scheme%LHS_mm_range == ALL_MOMENTS) THEN
      IF (scheme%dynamic_LMAX_on) THEN
         CALL lsquit('dynamic_LMAX not with ALL_MOMENTS at LHS for mm_get_grad_from_pkd_Vff',-1)
      ELSE
         lm_max = MIN(SIZE(LHS_mms%qlm_der_T,1),SIZE(Vff,1))
         packed_loop7: DO u = 1, SIZE(LHS_mms%paras)
            v = LHS_mms%paras(u)%id  ! LHS packed (batch) parameters 
            batch_map => LHS_mms%batch_map(v)%head
            batch_members7: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               icenta=LHS_mms%mom2atom(w,1)
               IF (w .GE. LHS_mms%startnuc .AND. w .LE. LHS_mms%endnuc ) THEN
               ! nuclear contribution
                  do ix = 1, 3
                     if (icenta .gt. 0) then
                        iscoor = (icenta-1)*3 + ix
                        fac = 0.5D0
                        if (ix == 1) ixb = 4
                        if (ix == 2) ixb = 2
                        if (ix == 3)then
                           ixb = 3
                           fac = -1.0D0
                        end if
                        grad(iscoor) = grad(iscoor) + LHS_mms%qlm_der_T(1,w,ix) * Vff(ixb,v)*fac
                     end if
                  end do
               ELSE
                  icentb=LHS_mms%mom2atom(w,2)
                  doa = .false.
                  dob = .false.
                  if(icenta .gt. 0 ) doa = .true.
                  if(icentb .gt. 0 ) dob = .true.
                  do iscool = 1, 6
                     iscoor = 0
                     if (iscool .eq. 1 .and. doa) iscoor = (icenta-1)*3 + 1
                     if (iscool .eq. 2 .and. doa) iscoor = (icenta-1)*3 + 2
                     if (iscool .eq. 3 .and. doa) iscoor = (icenta-1)*3 + 3
                     if (iscool .eq. 4 .and. dob) iscoor = (icentb-1)*3 + 1
                     if (iscool .eq. 5 .and. dob) iscoor = (icentb-1)*3 + 2
                     if (iscool .eq. 6 .and. dob) iscoor = (icentb-1)*3 + 3
                     if (iscoor .gt. 0) &
                  & grad(iscoor) = grad(iscoor) + DOT_PRODUCT(LHS_mms%qlm_der_T(:lm_max,w,iscool),Vff(:lm_max,v))
                  end do
                  energy = energy + DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
               END IF
               ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members7
               batch_map => batch_map%next
            END DO batch_members7
         END DO packed_loop7
     END IF
   ELSEIF(scheme%LHS_mm_range == NUCLEAR_ONLY) THEN
      IF (scheme%dynamic_LMAX_on) THEN
         CALL lsquit('dynamic_LMAX not with NUCLEAR_ONLY at LHS for mm_get_grad_from_pkd_Vff',-1)
      ELSE
         lm_max = MIN(SIZE(LHS_mms%qlm_der_T,1),SIZE(Vff,1))
         packed_loop6: DO u = 1, SIZE(LHS_mms%paras)
            v = LHS_mms%paras(u)%id  ! LHS packed (batch) parameters 
            batch_map => LHS_mms%batch_map(v)%head
            batch_members6: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               icenta=LHS_mms%mom2atom(w,1)
               do ix = 1, 3
                  if (icenta .gt. 0) then
                     iscoor = (icenta-1)*3 + ix
                     fac = 0.5D0
                     if (ix == 1) ixb = 4
                     if (ix == 2) ixb = 2
                     if (ix == 3)then
                        ixb = 3
                        fac = -1.0D0
                     end if
                     grad(iscoor) = grad(iscoor) + LHS_mms%qlm_der_T(1,w,ix) * Vff(ixb,v)*fac
                  end if
               end do
               ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members6
               batch_map => batch_map%next
            END DO batch_members6
         END DO packed_loop6
      END IF
   ELSE
      IF (scheme%dynamic_LMAX_on) THEN
         packed_loop: DO u = 1, SIZE(LHS_mms%paras)
            ! lm_max should never exceed raw_LMAX, and the dims
            ! of Vff also never exceed (1+raw_LMAX)**2
            lm_max = (1 + LHS_mms%paras(u)%Lmin + scheme%LEXTRA)**2
            lm_max = MIN(lm_max, INT(SIZE(Vff,1)) )
            v = LHS_mms%paras(u)%id  ! LHS packed (batch) moments
            batch_map => LHS_mms%batch_map(v)%head
            batch_members: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               icenta=LHS_mms%mom2atom(w,1)
               icentb=LHS_mms%mom2atom(w,2)
               doa = .false.
               dob = .false.
               if(icenta .gt. 0 ) doa = .true.
               if(icentb .gt. 0 ) dob = .true.
               do iscool = 1, 6
                  iscoor = 0
                  if (iscool .eq. 1 .and. doa) iscoor = (icenta-1)*3 + 1
                  if (iscool .eq. 2 .and. doa) iscoor = (icenta-1)*3 + 2
                  if (iscool .eq. 3 .and. doa) iscoor = (icenta-1)*3 + 3
                  if (iscool .eq. 4 .and. dob) iscoor = (icentb-1)*3 + 1
                  if (iscool .eq. 5 .and. dob) iscoor = (icentb-1)*3 + 2
                  if (iscool .eq. 6 .and. dob) iscoor = (icentb-1)*3 + 3
                  if (iscoor .gt. 0) &
                   & grad(iscoor) = grad(iscoor) + DOT_PRODUCT(LHS_mms%qlm_der_T(:lm_max,w,iscool),Vff(:lm_max,v))
               end do
               energy = energy + DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
               ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members
               batch_map => batch_map%next
            END DO batch_members

         END DO packed_loop

      ELSE
         lm_max = MIN(SIZE(LHS_mms%qlm_der_T,1),SIZE(Vff,1))

         packed_loop2: DO u = 1, SIZE(LHS_mms%paras)

            v = LHS_mms%paras(u)%id  ! LHS packed (batch) parameters 
            batch_map => LHS_mms%batch_map(v)%head

            batch_members2: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               icenta=LHS_mms%mom2atom(w,1)
               icentb=LHS_mms%mom2atom(w,2)
               doa = .false.
               dob = .false.
               if(icenta .gt. 0 ) doa = .true.
               if(icentb .gt. 0 ) dob = .true.
               do iscool = 1, 6
                  iscoor = 0
                  if (iscool .eq. 1 .and. doa) iscoor = (icenta-1)*3 + 1
                  if (iscool .eq. 2 .and. doa) iscoor = (icenta-1)*3 + 2
                  if (iscool .eq. 3 .and. doa) iscoor = (icenta-1)*3 + 3
                  if (iscool .eq. 4 .and. dob) iscoor = (icentb-1)*3 + 1
                  if (iscool .eq. 5 .and. dob) iscoor = (icentb-1)*3 + 2
                  if (iscool .eq. 6 .and. dob) iscoor = (icentb-1)*3 + 3
                  if (iscoor .gt. 0) &
                  & grad(iscoor) = grad(iscoor) + DOT_PRODUCT(LHS_mms%qlm_der_T(:lm_max,w,iscool),Vff(:lm_max,v))
               end do
               energy = energy + DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
               ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members2
               batch_map => batch_map%next
            END DO batch_members2

         END DO packed_loop2

      END IF
   END IF

   END SUBROUTINE mm_get_grad_from_pkd_Vff


!-------------------------------------------------------------------------------
! Build J_matrix components from contraction of the LHS moments and potentials 
!----------------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! This choice of Vff is fine here since the e-e interactions require a
! factor of 2 in the J-matrix, but NOT the e-n interactions, and if we
! let LHS = ELECTRONIC, and RHS = ELECTRONIC or NUCLEAR this is exactly
! what we get.

   SUBROUTINE mm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)

     USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(IN)    :: LHS_mms
      REAL(REALK),        INTENT(IN)    :: Vff(:,:) 
      REAL(REALK),        INTENT(INOUT) :: J_matrix(:,:) 

      REAL(REALK)   :: g
      INTEGER :: u,v, i,j, lm_max

      CALL mm_verify_Vff_input(scheme,LHS_mms,Vff,'J')

       ! note we put the IF test outside the loop for efficiency
      IF (scheme%dynamic_LMAX_on) THEN

         DO u = 1, SIZE(LHS_mms%paras)
            ! lm_max should never exceed raw_LMAX, and the dims
            ! of Vff also never exceed (1+raw_LMAX)**2
            lm_max = (1 + LHS_mms%paras(u)%Lmin + scheme%LEXTRA)**2
            lm_max = MIN(lm_max, INT(SIZE(Vff,1)) )
            v = LHS_mms%paras(u)%id
            g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
            i = LHS_mms%J_indices(v)%i_indx
            j = LHS_mms%J_indices(v)%j_indx
            IF (fit%LHS_AUX) i = 1
            J_matrix(i,j) = J_matrix(i,j) + g
            IF (i/=j .AND. (.NOT. fit%LHS_AUX)) J_matrix(j,i) = J_matrix(j,i) + g
         END DO

      ELSE

         ! although Vff should be the same size, we test for generality
         IF (SIZE(LHS_mms%qlm_T,1) /= SIZE(Vff,1)) STOP 'mm_get_J_from_Vff:2'
         lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
         DO u = 1, SIZE(LHS_mms%paras)
            v = LHS_mms%paras(u)%id
            g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
            i = LHS_mms%J_indices(v)%i_indx
            j = LHS_mms%J_indices(v)%j_indx
            IF (fit%LHS_AUX) i = 1
            J_matrix(i,j) = J_matrix(i,j) + g
            IF (i/=j .AND. (.NOT. fit%LHS_AUX)) J_matrix(j,i) = J_matrix(j,i) + g
         END DO

      END IF

   END SUBROUTINE mm_get_J_from_Vff

!-------------------------------------------------------------------------------
! Build J_matrix components from contracion of the LHS moments and potentials 
! This routine recognises that the LHS parameters may have been packed
! and thus expands the potential over LHS batches.
!----------------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments

   SUBROUTINE mm_get_J_from_pkd_Vff(scheme,LHS_mms,Vff,J_matrix)

     USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(IN)    :: LHS_mms
      REAL(REALK),        INTENT(IN)    :: Vff(:,:) 
      REAL(REALK),        INTENT(INOUT) :: J_matrix(:,:) 

      REAL(REALK)   :: g
      INTEGER :: u,v,w, i,j, lm_max
      TYPE(id_node), POINTER :: batch_map

      CALL mm_verify_Vff_input(scheme,LHS_mms,Vff,'J')

      ! Vff should now be the same size as the raw LHS parameters;
      ! However, the size of LHS %qlm_T will be larger;
      ! We need to "expand" Vff using the batch mapping.

       ! note we put the IF test outside the loops for efficiency
      IF (scheme%dynamic_LMAX_on) THEN

         packed_loop: DO u = 1, SIZE(LHS_mms%paras)

            ! lm_max should never exceed raw_LMAX, and the dims
            ! of Vff also never exceed (1+raw_LMAX)**2
            lm_max = (1 + LHS_mms%paras(u)%Lmin + scheme%LEXTRA)**2
            lm_max = MIN(lm_max, INT(SIZE(Vff,1)) )
            v = LHS_mms%paras(u)%id  ! LHS packed (batch) moments
            batch_map => LHS_mms%batch_map(v)%head

            batch_members: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
               i = LHS_mms%J_indices(w)%i_indx
               j = LHS_mms%J_indices(w)%j_indx
               IF (fit%LHS_AUX) i = 1
               J_matrix(i,j) = J_matrix(i,j) + g
               IF (i/=j .AND. (.NOT. fit%LHS_AUX)) J_matrix(j,i) = J_matrix(j,i) + g

                ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members
               batch_map => batch_map%next
            END DO batch_members 

         END DO packed_loop

      ELSE ! contract Vff with LHS basically using common raw_LMAX

         lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
         packed_loop2: DO u = 1, SIZE(LHS_mms%paras)

            v = LHS_mms%paras(u)%id  ! LHS packed (batch) moments
            batch_map => LHS_mms%batch_map(v)%head

            batch_members2: DO ! over batch list until pointer disassociated
               w = batch_map%id   ! raw LHS moment ID
               g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
               i = LHS_mms%J_indices(w)%i_indx
               j = LHS_mms%J_indices(w)%j_indx
               IF (fit%LHS_AUX) i = 1
               J_matrix(i,j) = J_matrix(i,j) + g
               IF (i/=j .AND. (.NOT. fit%LHS_AUX)) J_matrix(j,i) = J_matrix(j,i) + g

                ! only do next raw item in batch list if it exists
               IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members2
               batch_map => batch_map%next
            END DO batch_members2 

         END DO packed_loop2

      END IF

   END SUBROUTINE mm_get_J_from_pkd_Vff

!-------------------------------------------------------------------------------

END MODULE mm_Vff_processor

