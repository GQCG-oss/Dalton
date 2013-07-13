!===========================================================================
! dal_interface module contains all the routines that allow for
! exchanging information with the remaining part of dalton which
! may have different memory layout etc.

MODULE dal_interface
  use LStiming
  use precision
  use matrix_module, only: matrix, matrixp
  USE Matrix_operations, only: mat_trab, matrix_type, mtype_unres_dense,&
       & mat_sqnorm2, mat_dotproduct, mtype_dense, mat_to_full, mat_init, &
       & mat_zero, mat_free, mat_add, mat_print, mat_read_from_disk,mat_add,&
       & mat_abs_max_elm, mat_daxpy, mat_trans, mat_mul, mat_inv,&
       & mat_set_from_full, mat_assign, mat_scal
  use matrix_util,only: mat_get_isym, mcweeney_purify
   use TYPEDEFTYPE, only: lsitem, lssetting
   use TYPEDEF, only: typedef_setlist_valence2full,getNbasis, &
		& GCAO2AO_transform_matrixD2
   use dec_typedef_module, only: batchTOorb
   use Integralparameters
   use AO_TypeType, only: AOITEM
   use files, only: lsclose, lsopen
   use BUILDAOBATCH, only: build_batchesOfAOs, build_ao, &
        & determine_maxBatchOrbitalsize
   use daltoninfo, only: ls_free, build_ccfragmentlsitem
   use AO_Type, only: free_aoitem
   use lowdin_module, only: lowdin_diag
   use IntegralInterfaceDEC,only: II_precalc_decscreenmat, &
        & II_getBatchOrbitalScreenK, &
        & ii_get_decpacked4center_k_eri
   use memory_handling,only: mem_alloc, mem_dealloc
   use ks_settings, only: incremental_scheme, SaveF0andD0
   use lsdalton_fock_module, only: lsint_fock_data
   use screen_mod, only: DECscreenITEM, free_decscreen
   use IntegralInterfaceMOD
   use II_XC_interfaceModule
   use IIDFTINT, only: II_DFTsetFunc
#ifdef BUILD_GEN1INT
   ! debug GEN1INT
   use gen1int_host
#endif
#ifdef BUILD_CGTODIFF
   ! debug cgto_diff_eri_host
   use cgto_diff_eri_host_interface, only: cgto_diff_eri_DGD_Econt
#endif
INTERFACE di_GET_GbDs
	MODULE PROCEDURE di_GET_GbDsSingle, di_GET_GbDsArray
END INTERFACE
	
INTERFACE di_GET_GbDs_and_XC_linrsp
	MODULE PROCEDURE di_GET_GbDs_and_XC_linrsp_Single, di_GET_GbDs_and_XC_linrsp_Array
END INTERFACE

INTERFACE di_get_sigma_xc_cont_lsdalton
  MODULE PROCEDURE di_get_sigma_xc_contSingle_lsd,&
       & di_get_sigma_xc_contArray_lsd
END INTERFACE

   PUBLIC::&
        & di_debug_general, &
        & di_debug_general2, &
        & di_debug_4center_eri, &
        & di_debug_4center, &
        & di_decpacked, &
        & di_debug_ccfragment, &
        & di_get_overlap_and_H1, &
        & di_get_fock_LSDALTON, &
        & di_GET_GbDs, &
        & di_GET_GbDs_and_XC_linrsp, &
        & di_get_sigma_xc_cont_lsdalton,&
        & di_GET_GAOL_lsdalton, &
        & di_SCF_EnergyCont, &
        & di_get_dipole, &
        & fockenergy_F
   
   PRIVATE

   ! local pointer to H1 matrix allocated in linscf() subroutine. This allows
   ! access to H1 within this module, avoiding re-reading H1 from disk and
   ! avoiding extensive modifications to existing subroutine interfaces.
   type(MATRIX), pointer, save  :: lH1
CONTAINS
  subroutine di_debug_general(lupri,luerr,ls,nbast,S,D,debugProp)
      implicit none
      TYPE(lsitem),target :: ls
      integer,intent(in)    :: lupri,luerr,nbast
      type(matrix)          :: S,D
      real(realk) :: Etest,E
      logical :: debugProp
      IF(ls%input%dalton%DEBUGDECPACKED)then
         call di_decpacked(lupri,luerr,ls,nbast,D)
      ENDIF
      IF(ls%setting%scheme%interest)then
         print*,'di_debug_4center_eri_interest'
         call di_debug_4center_eri_interest(lupri,luerr,ls,nbast)
         print*,'di_debug_4center_eri_interest',ls%setting%scheme%interest
      endif
      if (ls%input%dalton%DEBUGMAGDERIV) then
         !Test geometric derivative overlap integrals
         call di_debug_magderiv_4center_eri(lupri, luerr,ls,nbast,D)
      endif
      IF (debugProp) THEN 
         call di_debug_PropertyIntegrals(lupri,luerr,ls%setting,nbast,S,D)
      ENDIF
      if (ls%input%dalton%DEBUGGEN1INT) then
#ifdef BUILD_GEN1INT
         !Test general 1 electron integrals by Bin Gao
         call gen1int_host_test(ls%setting,D,lupri)
         call di_gen1int_host_test(ls%setting,S,D,lupri)
#else
         call lsquit('.DEBUGGEN1INT requires OpenRSP -DBUILD_GEN1INT',lupri)
#endif
      endif
      if (ls%input%dalton%DEBUGCGTODIFF) then
#ifdef BUILD_CGTODIFF
         call cgto_diff_eri_DGD_Econt(D,Etest,ls%setting, lupri)
         WRITE(lupri,'(A,F14.9)')'Energi Contribution From cgto_diff_eri_DGD_Econt',Etest
         call II_get_Econt(LUPRI,LUERR,ls%SETTING,D,E,1,1.0E0_realk)
         WRITE(lupri,'(A,F14.9)')'Energi Contribution From Thermite Econt',E
         IF(ABS(E-Etest).GT.1.0E-14)THEN
            CALL LSQUIT('CGTO_DIFF_ERI E ERROR',lupri)
         ELSE
            WRITE(lupri,'(A)')'CGTO Energy Cont Succesfull'
         ENDIF
#else
         call lsquit('.DEBUGCGTODIFF requires OpenRSP -DBUILD_CGTODIFF',lupri)
#endif
      endif
    end subroutine di_debug_general

  subroutine di_debug_general2(lupri,luerr,ls,nbast,S,D)
      implicit none
      TYPE(lsitem),target :: ls
      integer,intent(in)    :: lupri,luerr,nbast
      type(matrix)          :: S,D
      
      if (ls%input%dalton%DEBUGEP) then
         !Test electrostatic potential
         call di_ep_test(ls%setting, D, lupri, luerr)
      endif
      if (ls%input%dalton%DEBUGGEODERIVOVERLAP) then
         !Test geometric derivative overlap integrals
         call di_geoderivOverlap_test(ls, D, lupri, luerr)
      endif
      if (ls%input%dalton%DEBUGGEODERIVKINETIC) then
         !Test geometric derivative kinetic integrals
         call di_geoderivKinetic_test(ls, D, lupri, luerr)
      endif
      if (ls%input%dalton%DEBUGGEODERIVKINETIC) then
         !Test geometric derivative nuclear-electron attraction integrals
         call di_geoderivNucel_test(ls, D, lupri, luerr)
      endif
      if (ls%input%dalton%DEBUGGEODERIVEXCHANGE) then
         !Test geometric derivative exchange integrals
         call di_geoderivExchange_test(ls, D, lupri, luerr)
      endif
      if (ls%input%dalton%DEBUGGEODERIVCOULOMB) then
         !Test geometric derivative Coulomb integrals
         call di_geoderivCoulomb_test(ls, D, lupri, luerr)
      endif
      if (ls%input%dalton%DEBUGMAGDERIVOVERLAP) then
         !Test magnetic derivative overlap integrals
         call di_magderivOverlap_test(ls, D, lupri, luerr)
      endif      
      if (ls%input%dalton%DEBUGscreen) then
         !Test screening integrals
         call di_screen_test(ls%setting, nbast, lupri, luerr)
      endif

    end subroutine di_debug_general2

#ifdef BUILD_GEN1INT
    subroutine di_gen1int_host_test(setting,S,D,lupri)
      implicit none
      type(lssetting)                :: setting
      integer                        :: lupri
      type(matrix),target,intent(in)        :: S,D
      !
      type(matrix) :: genS(1),h1,genh1(1)
      type(matrix),pointer :: Sx(:),genSx(:),h1x(:),genh1x(:),Tx(:)
      integer :: nbast,luerr,natoms,I,ndmat
      real(realk) :: THRESH
      type(matrixp) :: Dmat(1)
      real(realk),pointer :: grad(:,:),tmpgrad(:,:)
      type(matrixp)                  :: Dp(1)
      real(realk),pointer            :: reOrtho(:,:),gengradSx(:,:)
      Dmat(1)%p => D
      ndmat = 1
      THRESH=1.0E-15
      nbast = S%nrow
      luerr = lupri
      call mat_init(genS(1),nbast,nbast)
      call gen1int_host_get_overlap(setting, genS, lupri)
      call VerifyMatrices(MAT1=S,MAT2=genS(1),STRING='genS',THR=THRESH,lu=lupri)
      call mat_free(genS(1))
      WRITE(lupri,'(A)') ' genS Test Complet Successful' 
      nAtoms = setting%molecule(1)%p%nAtoms
      call mem_alloc(Sx,3*nAtoms)
      do I=1,3*nAtoms
         call mat_init(Sx(I),nbast,nbast)
      enddo
      call II_get_geoderivOverlap(Sx,natoms,setting,lupri,luerr)      

      call mem_alloc(genSx,3*nAtoms)
      do I=1,3*nAtoms
         call mat_init(genSx(I),nbast,nbast)
      enddo
      call gen1int_host_get_first_geoderiv_overlap(setting,genSx,nAtoms,lupri)
      do I=1,3*nAtoms
         call VerifyMatrices(Sx(I),genSx(I),'genSx',THRESH,lupri)
         call mat_free(Sx(I))
         call mat_free(genSx(I))
      enddo
      call mem_dealloc(Sx)
      call mem_dealloc(genSx)
      WRITE(lupri,'(A)') ' genSx Test Complet Successful' 

      call mem_alloc(gengradSx,3*nAtoms,1)
      gengradSx = 0.0E0_realk
      call gen1int_host_get_first_geoderiv_overlap_expval(setting, D, gengradSx, natoms,lupri)
      call mem_alloc(reOrtho,3,natoms)
      reOrtho = 0.0E0_realk
      Dp(1)%p => D
      call II_get_reorthoNormalization(reOrtho,Dp,1,setting,lupri,lupri)
      do I=1,Natoms
         IF(ABS(gengradSx(3*(I-1)+1,1)-reOrtho(1,I)).GT.1E-12_realk)THEN
            CALL LSQUIT('ERROR IN X di_gen1int_host_test Sx expval',lupri)
         ENDIF
         IF(ABS(gengradSx(3*(I-1)+2,1)-reOrtho(2,I)).GT.1E-12_realk)THEN
            CALL LSQUIT('ERROR IN Y di_gen1int_host_test Sx expval',lupri)
         ENDIF
         IF(ABS(gengradSx(3*(I-1)+3,1)-reOrtho(3,I)).GT.1E-12_realk)THEN
            CALL LSQUIT('ERROR IN Z di_gen1int_host_test Sx expval',lupri)
         ENDIF
      enddo
      call mem_dealloc(gengradSx)
      call mem_dealloc(reOrtho)
      WRITE(lupri,'(A)') ' expvalgenSx Test Complet Successful' 

      call mat_init(h1,nbast,nbast)
      call II_get_h1(lupri,luerr,setting,H1)
      call mat_init(genh1(1),nbast,nbast)
      call gen1int_host_get_h1(setting, genh1, lupri)
      call VerifyMatrices(h1,genh1(1),'genh1',THRESH,lupri)
      call mat_free(h1)
      call mat_free(genh1(1))
      WRITE(lupri,'(A)') ' genh1 Test Complet Successful' 
      

      call mem_alloc(genh1x,3*nAtoms)
      do I=1,3*nAtoms
         call mat_init(genh1x(I),nbast,nbast)
      enddo
      call gen1int_host_get_first_geoderiv_h1(setting,genh1x,nAtoms,lupri)
            
      call mem_alloc(h1x,3*nAtoms)
      do I=1,3*nAtoms
         call mat_init(h1x(I),nbast,nbast)
      enddo
      call mem_alloc(Tx,3*nAtoms)
      do I=1,3*nAtoms
         call mat_init(Tx(I),nbast,nbast)
      enddo
      call II_get_geoderivKinetic(Tx,Natoms,setting,lupri,luerr)
      call II_get_geoderivnucel(h1x,Natoms,setting,lupri,luerr)
      do I=1,3*nAtoms
         call mat_daxpy(1.0E0_realk, Tx(I), h1x(I))
         call VerifyMatrices(h1x(I),genh1x(I),'genh1xMatrix',THRESH,lupri)
         call mat_free(h1x(I))
         call mat_free(Tx(I))
      enddo
      call mem_dealloc(h1x)
      call mem_dealloc(Tx)
      WRITE(lupri,'(A)') ' genh1x Matrix Test Complete Successful' 

      call mem_alloc(Grad,3,nAtoms)
      GRAD=0.0E0_realk
      call mem_alloc(tmpGrad,3,nAtoms)
      CALL II_get_ne_gradient(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
      CALL DAXPY(3*natoms,1E0_realk,tmpGrad,1,Grad,1)
      CALL II_get_kinetic_gradient(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
      CALL DAXPY(3*natoms,1E0_realk,tmpGrad,1,Grad,1)

!      WRITE(lupri,*)'The h1grad'
!      call output(Grad,1,3,1,Natoms,3,natoms,1,lupri)

      do I=1,Natoms
         IF(ABS(mat_trAB(genh1x(3*(I-1)+1),D)-grad(1,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(genh1x(3*(I-1)+1),D)',mat_trAB(genh1x(3*(I-1)+1),D)
            print*,'grad(1,I)',grad(1,I)
            print*,'DIFF:',ABS(mat_trAB(genh1x(3*(I-1)+1),D)-grad(1,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN X di_gen1int_host_test_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(genh1x(3*(I-1)+2),D)-grad(2,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(genh1x(3*(I-1)+1),D)',mat_trAB(genh1x(3*(I-1)+2),D)
            print*,'grad(1,I)',grad(2,I)
            print*,'DIFF:',ABS(mat_trAB(genh1x(3*(I-1)+2),D)-grad(2,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN Y di_gen1int_host_test_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(genh1x(3*(I-1)+3),D)-grad(3,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(genh1x(3*(I-1)+1),D)',mat_trAB(genh1x(3*(I-1)+3),D)
            print*,'grad(1,I)',grad(3,I)
            print*,'DIFF:',ABS(mat_trAB(genh1x(3*(I-1)+3),D)-grad(3,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN Z di_gen1int_host_test_test',lupri)
         ENDIF
         call mat_free(genh1x(3*(I-1)+1))
         call mat_free(genh1x(3*(I-1)+2))
         call mat_free(genh1x(3*(I-1)+3))
      enddo

      call mem_dealloc(genh1x)
      call mem_dealloc(TMPGrad)
      WRITE(lupri,'(A)') ' genh1x Test Complet Successful' 
      !we use grad to test some more
      call mem_alloc(gengradSx,3*nAtoms,1)
      gengradSx = 0.0E0_realk
      call gen1int_host_get_first_geoderiv_h1_expval(setting, D, gengradSx, natoms,lupri)
      do I=1,Natoms
         IF(ABS(gengradSx(3*(I-1)+1,1)-grad(1,I)).GT.1E-8_realk)THEN
            print*,'gengradhx(3*(I-1)+1,1)',gengradSx(3*(I-1)+1,1)
            print*,'grad(1,I)          ',grad(1,I)
            CALL LSQUIT('ERROR IN X di_gen1int_host_test hx expval',lupri)
         ENDIF
         IF(ABS(gengradSx(3*(I-1)+2,1)-grad(2,I)).GT.1E-8_realk)THEN
            print*,'gengradhx(3*(I-1)+2,1)',gengradSx(3*(I-1)+2,1)
            print*,'grad(2,I)          ',grad(2,I)
            CALL LSQUIT('ERROR IN Y di_gen1int_host_test hx expval',lupri)
         ENDIF
         IF(ABS(gengradSx(3*(I-1)+3,1)-grad(3,I)).GT.1E-8_realk)THEN
            print*,'gengradhx(3*(I-1)+3,1)',gengradSx(3*(I-1)+3,1)
            print*,'grad(3,I)          ',grad(3,I)
            CALL LSQUIT('ERROR IN Z di_gen1int_host_test hx expval',lupri)
         ENDIF
      enddo
      call mem_dealloc(gengradSx)
      WRITE(lupri,'(A)') ' expvalgenhx Test Complet Successful' 
      call mem_dealloc(Grad)


!      NSIZE = (3*NATOMS)*(3*NATOMS+1)/2 !this only works for certain NATOMS values
      call mem_alloc(gengradSx,9*nAtoms*nAtoms,1)
      gengradSx = 0.0E0_realk
      call gen1int_host_get_second_geoderiv_overlap_expval(&
           & setting, D, gengradSx, natoms,lupri)
      call mem_alloc(genSx,9*nAtoms*nAtoms)
      do I=1,9*nAtoms*nAtoms
         call mat_init(genSx(I),nbast,nbast)
      enddo
      call gen1int_host_get_second_geoderiv_overlap(setting,genSx,nAtoms,lupri)
      do I=1,9*nAtoms*nAtoms
         IF(ABS(gengradSx(I,1)-mat_TrAB(genSx(I),D)).GT.1E-12_realk)THEN
            print*,'I=',I
            print*,'gengradSx(I,1)',gengradSx(I,1)
            print*,'mat_TrAB(genSx(I),D)',mat_TrAB(genSx(I),D)
            CALL LSQUIT('ERROR IN second_geoderiv_overlap',lupri)
         endif
      enddo
      WRITE(lupri,'(A)') ' second_geoderiv_overlap Completed Successful' 
      do I=1,9*nAtoms*nAtoms
         call mat_free(genSx(I))
      enddo
      call mem_dealloc(genSx)
      call mem_dealloc(gengradSx)

    end  subroutine di_gen1int_host_test
#endif

  subroutine di_debug_PropertyIntegrals(lupri,luerr,setting,nbast,S,Dmat)
      implicit none      
      type(lssetting)                :: setting
      integer                        :: lupri,luerr,nbast
      type(matrix),intent(in)        :: S,Dmat
      !
      integer                        :: propfile,X,NATOMS,B,ISYM
      type(matrix),pointer           :: matarray(:),OLDINT(:)
      type(matrix)                   :: tempm1,NONsym
      logical                        :: ReadOldProp,OnMaster
      real(realk)                    :: maxelm
      real(realk),pointer            :: expval(:)
      real(realk),parameter          :: thresh=1.0E-14
      character(len=8)   :: Label
      character(len=1)   :: CHRXYZ(-3:3)
      DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
      character(len=2)   :: ROTXYZ(6)
      DATA ROTXYZ /'XX','XY','XZ','YY','YZ','ZZ'/

      call mat_init(tempm1,S%nrow,S%ncol)
      call mat_init(NONsym,S%nrow,S%ncol)      

      call mat_trans(Dmat,NONsym)

      ReadOldProp = .FALSE.
      IF(ReadOldProp)THEN
         propfile=-1
         call lsopen(propfile,'propfile','OLD','UNFORMATTED')
         allocate(OLDINT(3))
         DO X=1,3
            call mat_init(OLDINT(X),nbast,nbast)
         ENDDO
      ENDIF

      allocate(matarray(3))
      DO X=1,3
         call mat_init(matarray(X),nbast,nbast)
      ENDDO
      NATOMS = setting%MOLECULE(1)%p%nAtoms

      !CALC ANGMOM
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,3,'ANGMOM ')

      !build NONsym mat
      call mat_daxpy(-0.349232E0_realk, matarray(1), NONsym)
      call mat_daxpy(+0.149232E0_realk, matarray(2), NONsym)
      call mat_daxpy(-0.249232E0_realk, matarray(3), NONsym)

      !Test
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP ANGMOM',CHRXYZ(X),' = ',ISYM
         if(ISYM.EQ.1)THEN
            WRITE(lupri,'(A,A,A)')'PROP ANGMOM',CHRXYZ(X),' IS SYMMETRIC'
            CALL LSQUIT('ANGMOM IS SYMMETRIC',lupri)
         elseif(isym.EQ.2.OR.isym.EQ.4)then
            WRITE(lupri,'(A,A,A)')'PROP ANGMOM',CHRXYZ(X),' IS ANTISYMMETRIC'
         else
            WRITE(lupri,'(A,A,A)')'PROP ANGMOM',CHRXYZ(X),' IS NON-SYMMETRIC'
            CALL LSQUIT('ANGMOM IS NON-SYMMETRIC',lupri)
         endif
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of ANGMOM',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.14)')'Maximum Element of ANGMOM',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read ANGMOM and verify
         OnMaster = .TRUE.
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 ANGMOM',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.8)')'PROPTEST4 ANGMOM',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF
      !CALC ANGLON
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,3,'ANGLON ')
      !Test
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh) !1 sym, 2 antisym, 3 nosym, 4 zero
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP ANGLON',CHRXYZ(X),' = ',ISYM
         if(ISYM.EQ.1)THEN
            WRITE(lupri,'(A,A,A)')'PROP ANGLON',CHRXYZ(X),' IS SYMMETRIC'
            CALL LSQUIT('ANGMOM IS SYMMETRIC',lupri)
         elseif(isym.EQ.2)then
            WRITE(lupri,'(A,A,A)')'PROP ANGLON',CHRXYZ(X),' IS ANTISYMMETRIC'
            CALL LSQUIT('ANGMOM IS ANTISYMMETRIC',lupri)
         else
            WRITE(lupri,'(A,A,A)')'PROP ANGLON',CHRXYZ(X),' IS NON-SYMMETRIC OR ZERO'
         endif
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of ANGLON',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of ANGLON',CHRXYZ(X),'=',maxelm
      ENDDO

      IF(ReadOldProp)THEN
         !read ANGLON and verify
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 ANGLON',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.14)')'PROPTEST4 ANGLON',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF

      !CALC 1ELPOT
      call II_get_prop(LUPRI,LUERR,SETTING,matarray(1:1),1,'1ELPOT ')
      ISYM = mat_get_isym(matarray(1),thresh) 
      WRITE(lupri,'(A,I2)')'SYMMETRY OF PROP 1ELPOT = ',ISYM
      IF(ReadOldProp)THEN
         !read and verify
         call mat_read_from_disk(propfile,OLDINT(1),OnMaster)
         call mat_add(1E0_realk,matarray(1),-1E0_realk,OLDINT(1),tempm1)
         WRITE(lupri,'(A,I2)')'PROPTEST3 1ELPOT=',mat_get_isym(tempm1,thresh)
         WRITE(lupri,'(A,F18.8)')'PROPTEST4 1ELPOT=',ABS(mat_trab(tempm1,tempm1))
         call II_get_nucel_mat(LUPRI,LUERR,SETTING,matarray(1))
         call mat_add(-1E0_realk,matarray(1),-1E0_realk,OLDINT(1),tempm1)
         WRITE(lupri,'(A,I2)')'PROPTEST1 1ELPOT=',mat_get_isym(tempm1,thresh)
         WRITE(lupri,'(A,F18.14)')'PROPTEST2 1ELPOT=',ABS(mat_trab(tempm1,tempm1))
      ELSE
         call II_get_nucel_mat(LUPRI,LUERR,SETTING,matarray(2))
         call mat_add(-1E0_realk,matarray(2),-1E0_realk,matarray(1),tempm1)
         WRITE(lupri,'(A,I2)')'PROPTEST5 1ELPOT=',mat_get_isym(tempm1,thresh)
         WRITE(lupri,'(A,F18.14)')'PROPTEST6 1ELPOT=',ABS(mat_trab(tempm1,tempm1))
      ENDIF

      !CALC LONMOM1
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,3,'LONMOM1')
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh) 
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP LONMOM1',CHRXYZ(X),' = ',ISYM
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of LONMOM1',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of LONMOM1',CHRXYZ(X),'=',maxelm
      ENDDO      
      
      !CALC LONMOM2
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,3,'LONMOM2')
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh) 
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP LONMOM2',CHRXYZ(X),' = ',ISYM
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of LONMOM2',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of LONMOM2',CHRXYZ(X),'=',maxelm
      ENDDO
      !CALC LONMOM
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,3,'LONMOM ')
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh) 
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP LONMOM',CHRXYZ(X),' = ',ISYM
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of LONMOM',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of LONMOM',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read and verify LONMOM
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 LONMOM',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.14)')'PROPTEST4 LONMOM',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF

      !CALC MAGMOM
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,3,'MAGMOM ')
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh) 
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP MAGMOM',CHRXYZ(X),' = ',ISYM
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of MAGMOM',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of MAGMOM',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read MAGMOM anf verify
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 MAGMOM',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.8)')'PROPTEST4 MAGMOM',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF

      !CALC S1MAG
      call II_get_magderivOverlap(matarray,setting,lupri,luerr)
      !Test
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP S1MAG',CHRXYZ(X),' = ',ISYM
         if(ISYM.EQ.1)THEN
            WRITE(lupri,'(A,A,A)')'PROP S1MAG',CHRXYZ(X),' IS SYMMETRIC'
            CALL LSQUIT('S1MAG IS SYMMETRIC',lupri)
         elseif(isym.EQ.2.OR.isym.EQ.4)then
            WRITE(lupri,'(A,A,A)')'PROP S1MAG',CHRXYZ(X),' IS ANTISYMMETRIC'
         else
            WRITE(lupri,'(A,A,A)')'PROP S1MAG',CHRXYZ(X),' IS NON-SYMMETRIC'
            CALL LSQUIT('S1MAG IS NON-SYMMETRIC',lupri)
         endif
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of S1MAG',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of S1MAG',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read S1MAG
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 S1MAG',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.14)')'PROPTEST4 S1MAG',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF
      !CALC S1MAGR
      call II_get_magderivOverlapR(matarray,setting,lupri,luerr)
      !Test S1MAGR
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP S1MAGR',CHRXYZ(X),' = ',ISYM
         if(ISYM.EQ.1)THEN
            WRITE(lupri,'(A,A,A)')'PROP S1MAGR',CHRXYZ(X),' IS SYMMETRIC'
!            CALL LSQUIT('S1MAG IS SYMMETRIC',lupri)
         elseif(isym.EQ.2.OR.isym.EQ.4)then
            WRITE(lupri,'(A,A,A)')'PROP S1MAGR',CHRXYZ(X),' IS ANTISYMMETRIC'
         else
            WRITE(lupri,'(A,A,A)')'PROP S1MAGR',CHRXYZ(X),' IS NON-SYMMETRIC'
!            CALL LSQUIT('S1MAG IS NON-SYMMETRIC',lupri)
         endif
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of S1MAGR',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of S1MAGR',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read S1MAGR
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 S1MAGR',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.14)')'PROPTEST4 S1MAGR',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
            IF(ABS(mat_trab(tempm1,tempm1)).GT.1.0E-10)THEN
               WRITE(lupri,*)'OLD S1MAGR',CHRXYZ(X)
               call mat_print(OLDINT(X),1,nbast,1,nbast,lupri)
               WRITE(lupri,*)'NEW S1MAGR',CHRXYZ(X)
               call mat_print(matarray(X),1,nbast,1,nbast,lupri)
            ENDIF
         ENDDO
      ENDIF
      !CALC S1MAGL
      call II_get_magderivOverlapL(matarray,setting,lupri,luerr)
      !Test S1MAGL
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP S1MAGL',CHRXYZ(X),' = ',ISYM
         if(ISYM.EQ.1)THEN
            WRITE(lupri,'(A,A,A)')'PROP S1MAGL',CHRXYZ(X),' IS SYMMETRIC'
!            CALL LSQUIT('S1MAG IS SYMMETRIC',lupri)
         elseif(isym.EQ.2.OR.isym.EQ.4)then
            WRITE(lupri,'(A,A,A)')'PROP S1MAGL',CHRXYZ(X),' IS ANTISYMMETRIC'
         else
            WRITE(lupri,'(A,A,A)')'PROP S1MAGL',CHRXYZ(X),' IS NON-SYMMETRIC'
!            CALL LSQUIT('S1MAG IS NON-SYMMETRIC',lupri)
         endif
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of S1MAGL',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of S1MAGL',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read S1MAGL
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 S1MAGL',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.14)')'PROPTEST4 S1MAGL',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
            IF(ABS(mat_trab(tempm1,tempm1)).GT.1.0E-10)THEN
               WRITE(lupri,*)'OLD S1MAGL',CHRXYZ(X)
               call mat_print(OLDINT(X),1,nbast,1,nbast,lupri)
               WRITE(lupri,*)'NEW S1MAGL',CHRXYZ(X)
               call mat_print(matarray(X),1,nbast,1,nbast,lupri)
            ENDIF
         ENDDO
      ENDIF

      IF(ReadOldProp)THEN
         DO X=1,3
            call mat_free(OLDINT(X))
         ENDDO
         deallocate(OLDINT)
         allocate(OLDINT(3*NATOMS))
         DO X=1,3*NATOMS
            call mat_init(OLDINT(X),nbast,nbast)
         ENDDO
      ENDIF
      DO X=1,3
         call mat_free(matarray(X))
      ENDDO
      deallocate(matarray)
      allocate(matarray(3*NATOMS))
      DO X=1,3*NATOMS
         call mat_init(matarray(X),nbast,nbast)         
      ENDDO
      !CALC PSO 
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,3*NATOMS,'PSO    ')
      DO X=1,3*NATOMS
         Label = 'PSO '//Char(X/100+48)//Char(mod(X,100)/10+48)&
              &//Char(mod(mod(X,100),10)+48)//' '         
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A8,A,I2)')'SYMMETRY OF PROP ',Label,' = ',ISYM
         WRITE(lupri,'(A,A8,A1,F18.9)')'Norm of ',Label,'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A8,A1,F18.8)')'Maximum Element of ',Label,'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read PSO and verify
         DO X=1,3*NATOMS
            Label = 'PSO '//Char(X/100+48)//Char(mod(X,100)/10+48)&
                 &//Char(mod(mod(X,100),10)+48)//' '         
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A8,A1,I2)')'PROPTEST3 ',Label,'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A8,A1,F18.14)')'PROPTEST4 ',Label,'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF
      ! TEST PSO EXPECTATION VALUE
      call mem_alloc(expval,3*NATOMS)
      call II_get_prop_expval(LUPRI,LUERR,SETTING,expval,(/ NONsym /),1,3*NATOMS,'PSO    ')
      DO X=1,3*NATOMS
         IF(ABS(expval(X)-mat_dotproduct(matarray(X),NONsym)).GT.1.0E-12_realk)THEN
            write(lupri,*)'X',X
            write(lupri,*)'expval(X)',expval(X)
            write(lupri,*)'mat_trab(matarray(X),Dmat)',mat_dotproduct(matarray(X),NONsym)
            call lsquit('error in expectation value PSO in di_debug_PropertyIntegrals',lupri)
         ENDIF
      ENDDO
      call mem_dealloc(expval)

      IF(ReadOldProp)THEN
         DO X=1,3*NATOMS
            call mat_free(OLDINT(X))
         ENDDO
         deallocate(OLDINT)
         allocate(OLDINT(9*NATOMS))
         DO X=1,9*NATOMS
            call mat_init(OLDINT(X),nbast,nbast)
         ENDDO
      ENDIF
      DO X=1,3*NATOMS
         call mat_free(matarray(X))
      ENDDO
      deallocate(matarray)
      allocate(matarray(9*NATOMS))
      DO X=1,9*NATOMS
         call mat_init(matarray(X),nbast,nbast)         
      ENDDO

      !CALC NSTLON
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,9*NATOMS,'NSTLON ')
      DO X=1,3*NATOMS
         DO B=1,3
            Label = Char(X/100+48)//Char(mod(X,100)/10+48)&
                 &//Char(mod(mod(X,100),10)+48)//'NSLO'//CHRXYZ(B)
            WRITE(lupri,'(A,A8,A,I2)')'SYMMETRY OF PROP ',Label,' = ',ISYM
            WRITE(lupri,'(A,A8,A1,F18.9)')'Norm of ',Label,'=',&
                 &sqrt(mat_sqnorm2(matarray(B+(X-1)*3))/nbast)
            call mat_abs_max_elm(matarray(B+(X-1)*3), maxelm) 
            WRITE(lupri,'(A,A8,A1,F18.8)')'Maximum Element of ',Label,'=',maxelm
         ENDDO
      ENDDO
      IF(ReadOldProp)THEN
         !read NSTLON
         DO X=1,3*NATOMS
            DO B=1,3
               Label = Char(X/100+48)//Char(mod(X,100)/10+48)&
                    &//Char(mod(mod(X,100),10)+48)//'NSLO'//CHRXYZ(B)
               call mat_read_from_disk(propfile,OLDINT(B+(X-1)*3),OnMaster)
               call mat_add(1E0_realk,matarray(B+(X-1)*3),-1E0_realk,OLDINT(B+(X-1)*3),tempm1)
               WRITE(lupri,'(A,A8,A1,I2)')'PROPTEST3 ',Label,'=',mat_get_isym(tempm1,thresh)
               WRITE(lupri,'(A,A8,A1,F18.14)')'PROPTEST4 ',Label,'=',ABS(mat_trab(tempm1,tempm1))
            ENDDO
         ENDDO
      ENDIF
      ! TEST NSTLON EXPECTATION VALUE
      call mem_alloc(expval,9*NATOMS)
      call II_get_prop_expval(LUPRI,LUERR,SETTING,expval,(/ NONsym /),1,9*NATOMS,'NSTLON ')
      DO X=1,9*NATOMS
         IF(ABS(expval(X)-mat_dotproduct(matarray(X),NONsym)).GT.1.0E-12_realk)THEN
            write(lupri,*)'X',X
            write(lupri,*)'expval(X)',expval(X)
            write(lupri,*)'mat_trab(matarray(X),Dmat)',mat_dotproduct(matarray(X),NONsym)
            call lsquit('error in expectation value NSTLON in di_debug_PropertyIntegrals',lupri)
         ENDIF
      ENDDO
      call mem_dealloc(expval)

      !CALC NSTNOL
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,9*NATOMS,'NSTNOL ')
      DO X=1,3*NATOMS
         DO B=1,3
            Label = Char(X/100+48)//Char(mod(X,100)/10+48)&
                 &//Char(mod(mod(X,100),10)+48)//'NSNL'//CHRXYZ(B)
            WRITE(lupri,'(A,A8,A,I2)')'SYMMETRY OF PROP ',Label,' = ',ISYM
            WRITE(lupri,'(A,A8,A1,F18.9)')'Norm of ',Label,'=',&
                 &sqrt(mat_sqnorm2(matarray(B+(X-1)*3))/nbast)
            call mat_abs_max_elm(matarray(B+(X-1)*3), maxelm) 
            WRITE(lupri,'(A,A8,A1,F18.8)')'Maximum Element of ',Label,'=',maxelm
         ENDDO
      ENDDO
      IF(ReadOldProp)THEN
         !read NSTNOL
         DO X=1,3*NATOMS
            DO B=1,3
               call mat_read_from_disk(propfile,OLDINT(B+(X-1)*3),OnMaster)
               Label = Char(X/100+48)//Char(mod(X,100)/10+48)&
                    &//Char(mod(mod(X,100),10)+48)//'NSNL'//CHRXYZ(B)
               call mat_add(1E0_realk,matarray(B+(X-1)*3),-1E0_realk,OLDINT(B+(X-1)*3),tempm1)
               WRITE(lupri,'(A,A8,A1,I2)')'PROPTEST3 ',Label,'=',mat_get_isym(tempm1,thresh)
               WRITE(lupri,'(A,A8,A1,F18.14)')'PROPTEST4 ',Label,'=',ABS(mat_trab(tempm1,tempm1))
            ENDDO
         ENDDO
      ENDIF
      !TEST NSTNOL expectation value
      call mem_alloc(expval,9*NATOMS)
      call II_get_prop_expval(LUPRI,LUERR,SETTING,expval,(/ NONsym /),1,9*NATOMS,'NSTNOL ')
      DO X=1,9*NATOMS
         IF(ABS(expval(X)-mat_dotproduct(matarray(X),NONsym)).GT.1.0E-12_realk)THEN
            write(lupri,*)'expval(X)',expval(X)
            write(lupri,*)'mat_trab(matarray(X),Dmat)',mat_dotproduct(matarray(X),NONsym)
            call lsquit('error in expectation value NSTNOL in di_debug_PropertyIntegrals',lupri)
         ENDIF
      ENDDO
      call mem_dealloc(expval)

      !CALC NST
      call II_get_prop(LUPRI,LUERR,SETTING,matarray,9*NATOMS,'NST    ')
      DO X=1,3*NATOMS
         DO B=1,3
            Label = Char(X/100+48)//Char(mod(X,100)/10+48)&
                 &//Char(mod(mod(X,100),10)+48)//' NST'//CHRXYZ(B)
            WRITE(lupri,'(A,A8,A,I2)')'SYMMETRY OF PROP ',Label,' = ',ISYM
            WRITE(lupri,'(A,A8,A1,F18.9)')'Norm of ',Label,'=',&
                 &sqrt(mat_sqnorm2(matarray(B+(X-1)*3))/nbast)
            call mat_abs_max_elm(matarray(B+(X-1)*3), maxelm) 
            WRITE(lupri,'(A,A8,A1,F18.8)')'Maximum Element of ',Label,'=',maxelm
         ENDDO
      ENDDO
      IF(ReadOldProp)THEN
         !read NST
         DO X=1,3*NATOMS
            DO B=1,3
               call mat_read_from_disk(propfile,OLDINT(B+(X-1)*3),OnMaster)
               Label = Char(X/100+48)//Char(mod(X,100)/10+48)&
                    &//Char(mod(mod(X,100),10)+48)//' NST'//CHRXYZ(B)
               call mat_add(1E0_realk,matarray(B+(X-1)*3),-1E0_realk,OLDINT(B+(X-1)*3),tempm1)
               WRITE(lupri,'(A,A8,A1,I2)')'PROPTEST3 ',Label,'=',mat_get_isym(tempm1,thresh)
               WRITE(lupri,'(A,A8,A1,F18.14)')'PROPTEST4 ',Label,'=',ABS(mat_trab(tempm1,tempm1))
            ENDDO
         ENDDO
      ENDIF
      !TEST NST expectation value
      call mem_alloc(expval,9*NATOMS)
      call II_get_prop_expval(LUPRI,LUERR,SETTING,expval,(/ NONsym /),1,9*NATOMS,'NST    ')
      DO X=1,9*NATOMS
         IF(ABS(expval(X)-mat_dotproduct(matarray(X),NONsym)).GT.1.0E-12_realk)THEN
            write(lupri,*)'expval(X)',expval(X)
            write(lupri,*)'mat_trab(matarray(X),Dmat)',mat_dotproduct(matarray(X),NONsym)
            call lsquit('error in expectation value NST in di_debug_PropertyIntegrals',lupri)
         ENDIF
      ENDDO
      call mem_dealloc(expval)

      !CALC DIPVEL
      call II_get_prop(LUPRI,LUERR,SETTING,matarray(1:3),3,'DIPVEL ')
      !Test
      DO X=1,3
         WRITE(lupri,'(A,A1)')'PROP DIPVEL',CHRXYZ(X)
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP DIPVEL',CHRXYZ(X),' = ',ISYM
         if(ISYM.EQ.1)THEN
            WRITE(lupri,'(A,A,A)')'PROP DIPVEL',CHRXYZ(X),' IS SYMMETRIC'
            CALL LSQUIT('DIPVEL IS SYMMETRIC',lupri)
         elseif(isym.EQ.2.OR.isym.EQ.4)then
            WRITE(lupri,'(A,A,A)')'PROP DIPVEL',CHRXYZ(X),' IS ANTISYMMETRIC'
         else
            WRITE(lupri,'(A,A,A)')'PROP DIPVEL',CHRXYZ(X),' IS NON-SYMMETRIC'
            CALL LSQUIT('DIPVEL IS NON-SYMMETRIC',lupri)
         endif
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of DIPVEL',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of DIPVEL',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read DIPVEL and verify
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 DIPVEL',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.14)')'PROPTEST4 DIPVEL',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF
      !CALC ROTSTR
      call II_get_prop(LUPRI,LUERR,SETTING,matarray(1:6),6,'ROTSTR ')
      DO X=1,6
         Label = ROTXYZ(X)//'ROTSTR'
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A8,A,I2)')'SYMMETRY OF PROP ',Label,' = ',ISYM
         WRITE(lupri,'(A,A8,A1,F18.9)')'Norm of ',Label,'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A8,A1,F18.8)')'Maximum Element of ',Label,'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read ROTSTR and verify
         DO X=1,6
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            Label = ROTXYZ(X)//'ROTSTR'
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A8,A1,I2)')'PROPTEST3 ',Label,'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A8,A1,F18.14)')'PROPTEST4 ',Label,'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF
      !CALC THETA
      call II_get_prop(LUPRI,LUERR,SETTING,matarray(1:6),6,'THETA  ')
      DO X=1,6
         Label = ROTXYZ(X)//'THETA'
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A8,A,I2)')'SYMMETRY OF PROP ',Label,' = ',ISYM
         WRITE(lupri,'(A,A8,A1,F18.9)')'Norm of ',Label,'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A8,A1,F18.8)')'Maximum Element of ',Label,'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read THETA and verify
         DO X=1,6
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            Label = ROTXYZ(X)//'THETA'
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A8,A1,I2)')'PROPTEST3 ',Label,'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A8,A1,F18.14)')'PROPTEST4 ',Label,'=',ABS(mat_trab(tempm1,tempm1))
            IF(ABS(mat_trab(tempm1,tempm1)).GT.1.0E-10)THEN
               WRITE(lupri,*)'OLD LABEL:',Label
               call mat_print(OLDINT(X),1,nbast,1,nbast,lupri)
               WRITE(lupri,*)'NEW LABEL:',Label
               call mat_print(matarray(X),1,nbast,1,nbast,lupri)
            ENDIF
         ENDDO
      ENDIF
      !CALC DIPLEN
      call II_get_prop(LUPRI,LUERR,SETTING,matarray(1:3),3,'DIPLEN ')
      !Test
      DO X=1,3
         ISYM = mat_get_isym(matarray(X),thresh)
         WRITE(lupri,'(A,A1,A,I2)')'SYMMETRY OF PROP DIPLEN',CHRXYZ(X),' = ',ISYM
         if(ISYM.EQ.1.OR.isym.EQ.4)THEN
            WRITE(lupri,'(A,A,A)')'PROP DIPLEN',CHRXYZ(X),' IS SYMMETRIC OR ZERO'
         else
            WRITE(lupri,'(A,A,A)')'PROP DIPLEN',CHRXYZ(X),' IS NOT SYMMETRIC'
            call mat_print(matarray(X),1,nbast,1,nbast,lupri)
            CALL LSQUIT('DIPLEN IS NOT SYMMETRIC',lupri)
         endif
         WRITE(lupri,'(A,A1,A1,F18.9)')'Norm of DIPLEN',CHRXYZ(X),'=',&
              &sqrt(mat_sqnorm2(matarray(X))/nbast)
         call mat_abs_max_elm(matarray(X), maxelm) 
         WRITE(lupri,'(A,A1,A1,F18.8)')'Maximum Element of DIPLEN',CHRXYZ(X),'=',maxelm
      ENDDO
      IF(ReadOldProp)THEN
         !read DIPLEN and verify
         DO X=1,3
            call mat_read_from_disk(propfile,OLDINT(X),OnMaster)
            call mat_add(1E0_realk,matarray(X),-1E0_realk,OLDINT(X),tempm1)
            WRITE(lupri,'(A,A1,A1,I2)')'PROPTEST3 DIPLEN',CHRXYZ(X),'=',mat_get_isym(tempm1,thresh)
            WRITE(lupri,'(A,A1,A1,F18.14)')'PROPTEST4 DIPLEN',CHRXYZ(X),'=',ABS(mat_trab(tempm1,tempm1))
         ENDDO
      ENDIF
      !CALC OVERLAP
      call II_get_prop(LUPRI,LUERR,SETTING,matarray(1:1),1,'Overlap')
      !verify OVERLAP
      call mat_add(1E0_realk,matarray(1),-1E0_realk,S,tempm1)
      WRITE(lupri,*)'PROPTEST1 OVERLAP ',mat_get_isym(tempm1,thresh)
      WRITE(lupri,*)'PROPTEST2 OVERLAP ',ABS(mat_trab(tempm1,tempm1))
      DO X=1,9*NATOMS
         call mat_free(matarray(X))
      ENDDO
      deallocate(matarray)
      IF(ReadOldProp)THEN
         DO X=1,9*NATOMS
            call mat_free(OLDINT(X))
         ENDDO
         deallocate(OLDINT)
         call lsclose(propfile,'KEEP')
      ENDIF
      
      call mat_free(tempm1)
      call mat_free(NONsym)
    end subroutine di_debug_PropertyIntegrals

   SUBROUTINE di_debug_4center_eri(lupri,lu_err,ls,nbast)
     IMPLICIT NONE
     integer,intent(in)      :: lupri,lu_err,nbast
     type(lsitem),intent(inout) :: ls
     !
     real(realk),pointer   :: integrals(:,:,:,:)
     integer :: iB,iA,iC,iD
     !THIS IS A SPECIAL HARDCODED DEBUG ROUTINE FOR WATER (TESTCASE LS_DEBUG4CENTER)
     call mem_alloc(integrals,nbast,nbast,nbast,nbast)
     call II_get_4center_eri(LUPRI,LU_ERR,ls%setting,integrals,nbast,nbast,nbast,nbast)
     WRITE(lupri,*)'di_debug_4center_eri: some random 4 center integrals'
     !SAME ATOM ON ALL 4
     iA = 1; iB = 1; iC = 1; iD = 1
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 4; iD = 4
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 9; iD = 9
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 6; iD = 9
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON 1,2,3 and new on 4.
     iA = 1; iB = 1; iC = 1; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 4; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 9; iD = 12
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 6; iD = 12
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON LHS and same on RHS.
     iA = 1; iB = 1; iC = 11; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 11; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 11; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 11; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON LHS and different on RHS.
     iA = 1; iB = 1; iC = 13; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 2; iC = 13; iD = 11
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 3; iB = 3; iC = 10; iD = 12
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 2; iB = 3; iC = 10; iD = 12
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     !SAME ATOM ON 1 and same on 2,3,4
     iA = 10; iB = 10; iC = 10; iD = 10
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 11; iB = 10; iC = 10; iD = 10
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 12; iB = 10; iC = 10; iD = 10
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 13; iB = 10; iC = 10; iD = 10
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 10; iB = 1; iC = 1; iD = 1
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 11; iB = 1; iC = 1; iD = 1
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 12; iB = 1; iC = 1; iD = 1
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     iA = 13; iB = 1; iC = 1; iD = 1
     WRITE(lupri,'(4I3,F16.12)')iA,iB,iC,iD, integrals(ia,ib,ic,id)
     call mem_dealloc(integrals)

   END SUBROUTINE di_debug_4center_eri

   SUBROUTINE di_debug_4center(lupri,lu_err,ls,nbast,D)
      IMPLICIT NONE
      TYPE(Matrix),intent(in) :: D
      integer,intent(in)      :: lupri,lu_err,nbast
      type(lsitem),intent(inout) :: ls
      !
      TYPE(Matrix)  :: tempm1,tempm2,tempm3
      real(realk),pointer   :: integrals(:,:,:,:)
      real(realk) :: CoulombFactor,FAB,int
      integer,pointer :: orb2batch(:),batchdim(:)
      type(batchTOorb),pointer :: batch2orb(:)
      integer :: Abatch,Bbatch,IORB,JORB,dimA,dimB,J,nbatches
      integer :: iB,iA,iC,iD,offset,offset2,BORB,AORB,iprint

      CoulombFactor = 2E0_realk
      !build the coulomb matrix 
      call mat_init(tempm1,nbast,nbast)
      call mat_init(tempm2,nbast,nbast)
      call mat_init(tempm3,nbast,nbast)
      CALL II_get_coulomb_mat(lupri,lu_err,ls%setting,D,tempm1,1)

      allocate(orb2batch(nbast))
      CALL II_getBatchOrbitalInfo(ls%setting,nBast,orb2Batch,&
           & nBatches,lupri,lu_err)
      

      allocate(batchdim(nbatches))
      batchdim = 0
      DO IORB=1,nbast
         batchdim(orb2batch(IORB)) = batchdim(orb2batch(IORB))+1 
      ENDDO

      allocate(batch2orb(nbatches))
      DO Abatch=1,nbatches
         allocate(batch2orb(Abatch)%orbindex(batchdim(Abatch)))
         batch2orb(Abatch)%orbindex = 0
         batch2orb(Abatch)%norbindex = 0
      ENDDO
      DO IORB=1,nbast
         Abatch = orb2batch(IORB)
         batch2orb(Abatch)%norbindex = batch2orb(Abatch)%norbindex+1 
         J = batch2orb(Abatch)%norbindex
         batch2orb(Abatch)%orbindex(J) = IORB
      ENDDO
      DO Bbatch=1,nbatches
         dimB = batchdim(Bbatch)
         DO Abatch=1,nbatches
            dimA = batchdim(Abatch)
            nullify(integrals)
            allocate(integrals(dimA,dimB,nbast,nbast))
!            integrals=0E0_realk
            CALL II_get_ABres_4CenterEri(ls%setting,integrals,Abatch,&
                 & Bbatch,dimA,dimB,nbast,lupri,lu_err)            
            DO iB = 1,dimB
               !translate to orbitalindex: B orbital
               BORB = batch2orb(Bbatch)%orbindex(iB)
               offset2 = (BORB-1)*nbast
               DO iA = 1,dimA
                  !translate to orbitalindex: B orbital
                  AORB = batch2orb(Abatch)%orbindex(iA) 
                  FAB =0E0_realk
                  DO id=1,nbast
                     offset=(id-1)*nbast
                     DO ic=1,nbast
                        FAB=FAB+D%elms(ic+offset)*CoulombFactor*integrals(ia,ib,ic,id)
                     ENDDO
                  ENDDO
                  tempm2%elms(AORB+offset2) = FAB
               ENDDO
            ENDDO
            deallocate(integrals)
            nullify(integrals)
!            call mem_dealloc(integrals)
         ENDDO
      ENDDO
      DO Abatch=1,nbatches
         deallocate(batch2orb(Abatch)%orbindex)
      ENDDO
      deallocate(batch2orb)
      deallocate(batchdim)
      deallocate(orb2batch)
      call mat_add(1E0_realk,tempm1,-1E0_realk,tempm2,tempm3)
      write(lupri,*) 'QQQ DI_DEBUG_4CENTER STD ',mat_trab(tempm1,tempm1)
      write(lupri,*) 'QQQ DI_DEBUG_4CENTER CHOL',mat_trab(tempm2,tempm2)
      IF(ABS(mat_trab(tempm3,tempm3)).LE. 1E-15_realk)THEN
         write(lupri,*)'QQQ SUCCESFUL CHOLESKY TEST'
      ELSE
         CALL lsQUIT('CHOLESKY TEST TEST',lupri)
      ENDIF
      call mat_free(tempm1)
      call mat_free(tempm2)
      call mat_free(tempm3)

    END SUBROUTINE DI_DEBUG_4CENTER

    SUBROUTINE di_debug_magderiv_4center_eri(lupri,luerr,ls,nbast,D)
      IMPLICIT NONE
      integer,intent(in)      :: lupri,luerr,nbast
      type(lsitem),intent(inout) :: ls
      type(matrix) :: D
      !
      integer :: I
      type(matrix) :: Fx(3),Kx(3),Jx(3),tempm1,tempm2,Fx2(3),NONsym
      real(realk),parameter :: D4=4.0E0_realk

      call mat_init(NONsym,nbast,nbast)
      call mat_trans(D,NONsym)
      !symmetric Dmat

      !magderiv_4center_eri
      call mat_init(Fx(1),nbast,nbast)
      call mat_init(Fx(2),nbast,nbast)
      call mat_init(Fx(3),nbast,nbast)
      call II_get_magderiv_4center_eri(LUPRI,LUERR,ls%setting,nbast,D,Fx) 
      write(lupri,'(A,F16.8)') 'QQQ FxFx magderiv',D4*mat_trab(Fx(1),Fx(1))
      write(lupri,'(A,F16.8)') 'QQQ FxFy magderiv',D4*mat_trab(Fx(1),Fx(2))
      write(lupri,'(A,F16.8)') 'QQQ FxFz magderiv',D4*mat_trab(Fx(1),Fx(3))
      write(lupri,'(A,F16.8)') 'QQQ FyFy magderiv',D4*mat_trab(Fx(2),Fx(2))
      write(lupri,'(A,F16.8)') 'QQQ FyFz magderiv',D4*mat_trab(Fx(2),Fx(3))
      write(lupri,'(A,F16.8)') 'QQQ FzFz magderiv',D4*mat_trab(Fx(3),Fx(3))

      !build antisym dmat for next testing
      call mat_daxpy(+0.449232E0_realk, D, NONsym)
      call mat_daxpy(-0.349232E0_realk, Fx(1), NONsym)
      call mat_daxpy(+0.149232E0_realk, Fx(2), NONsym)
      call mat_daxpy(-0.249232E0_realk, Fx(3), NONsym)

      !magderivF
      call mat_init(Fx2(1),nbast,nbast)
      call mat_init(Fx2(2),nbast,nbast)
      call mat_init(Fx2(3),nbast,nbast)
      call II_get_magderivF(LUPRI,LUERR,ls%setting,nbast,(/D/),Fx2) 
      call mat_init(tempm2,nbast,nbast)
      DO I =1,3
         call mat_add(1E0_realk,Fx2(I),-1E0_realk,Fx(I),tempm2)
         WRITE(lupri,*)'DIFF F1=',mat_trab(tempm2,tempm2)
         IF(ABS(mat_trab(tempm2,tempm2)).GT.1.0E-10)CALL LSQUIT('Magderiv Error',-1)
         call mat_free(Fx2(I))
      ENDDO
      call mat_free(tempm2)

      !magderivK + magderivJ
      call mat_init(Kx(1),nbast,nbast)
      call mat_init(Kx(2),nbast,nbast)
      call mat_init(Kx(3),nbast,nbast)
      call II_get_magderivK(LUPRI,LUERR,ls%setting,nbast,(/D/),Kx) 
      write(lupri,'(A,F16.8)') 'QQQ KxKx magderiv',D4*mat_trab(Kx(1),Kx(1))
      write(lupri,'(A,F16.8)') 'QQQ KxKy magderiv',D4*mat_trab(Kx(1),Kx(2))
      write(lupri,'(A,F16.8)') 'QQQ kxKz magderiv',D4*mat_trab(Kx(1),Kx(3))
      write(lupri,'(A,F16.8)') 'QQQ KyKy magderiv',D4*mat_trab(Kx(2),Kx(2))
      write(lupri,'(A,F16.8)') 'QQQ KyKz magderiv',D4*mat_trab(Kx(2),Kx(3))
      write(lupri,'(A,F16.8)') 'QQQ KzKz magderiv',D4*mat_trab(Kx(3),Kx(3))
      call mat_init(Jx(1),nbast,nbast)
      call mat_init(Jx(2),nbast,nbast)
      call mat_init(Jx(3),nbast,nbast)
      call II_get_magderivJ(LUPRI,LUERR,ls%setting,nbast,(/D/),Jx) 
      write(lupri,'(A,F16.8)') 'QQQ JxJx magderiv',D4*mat_trab(Jx(1),Jx(1))
      write(lupri,'(A,F16.8)') 'QQQ JxJy magderiv',D4*mat_trab(Jx(1),Jx(2))
      write(lupri,'(A,F16.8)') 'QQQ JxJz magderiv',D4*mat_trab(Jx(1),Jx(3))
      write(lupri,'(A,F16.8)') 'QQQ JyJy magderiv',D4*mat_trab(Jx(2),Jx(2))
      write(lupri,'(A,F16.8)') 'QQQ JyJz magderiv',D4*mat_trab(Jx(2),Jx(3))
      write(lupri,'(A,F16.8)') 'QQQ JzJz magderiv',D4*mat_trab(Jx(3),Jx(3))
      call mat_init(tempm1,nbast,nbast)
      call mat_init(tempm2,nbast,nbast)
      DO I =1,3
         call mat_add(1E0_realk,Jx(I),1E0_realk,Kx(I),tempm1)
         call mat_add(1E0_realk,tempm1,-1E0_realk,Fx(I),tempm2)
         WRITE(lupri,*)'DIFF F2=',mat_trab(tempm2,tempm2)
         IF(ABS(mat_trab(tempm2,tempm2)).GT.1.0E-10)CALL LSQUIT('Magderiv Error',-1)
         call mat_free(Kx(I))
         call mat_free(Jx(I))
         call mat_free(Fx(I))
      ENDDO
      call mat_free(tempm1)
      call mat_free(tempm2)

      !Nonsymmetric Dmat

      !magderiv_4center_eri
      call mat_init(Fx(1),nbast,nbast)
      call mat_init(Fx(2),nbast,nbast)
      call mat_init(Fx(3),nbast,nbast)
      call II_get_magderiv_4center_eri(LUPRI,LUERR,ls%setting,nbast,NONsym,Fx) 
      write(lupri,'(A,F16.8)') 'QQQ NONsym FxFx magderiv',D4*mat_trab(Fx(1),Fx(1))
      write(lupri,'(A,F16.8)') 'QQQ NONsym FxFy magderiv',D4*mat_trab(Fx(1),Fx(2))
      write(lupri,'(A,F16.8)') 'QQQ NONsym FxFz magderiv',D4*mat_trab(Fx(1),Fx(3))
      write(lupri,'(A,F16.8)') 'QQQ NONsym FyFy magderiv',D4*mat_trab(Fx(2),Fx(2))
      write(lupri,'(A,F16.8)') 'QQQ NONsym FyFz magderiv',D4*mat_trab(Fx(2),Fx(3))
      write(lupri,'(A,F16.8)') 'QQQ NONsym FzFz magderiv',D4*mat_trab(Fx(3),Fx(3))

      !magderivK
!      call mat_init(Kx(1),nbast,nbast)
!      call mat_init(Kx(2),nbast,nbast)
!      call mat_init(Kx(3),nbast,nbast)
!      call II_get_magderivK_4center_eri(LUPRI,LUERR,ls%setting,nbast,NONsym,Kx) 
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym KxKx magderiv',D4*mat_trab(Kx(1),Kx(1))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym KxKy magderiv',D4*mat_trab(Kx(1),Kx(2))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym kxKz magderiv',D4*mat_trab(Kx(1),Kx(3))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym KyKy magderiv',D4*mat_trab(Kx(2),Kx(2))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym KyKz magderiv',D4*mat_trab(Kx(2),Kx(3))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym KzKz magderiv',D4*mat_trab(Kx(3),Kx(3))
!      write(lupri,*)'Kx_4center_eri'
!      call mat_print(Kx(1),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Ky_4center_eri'
!      call mat_print(Kx(2),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Kz_4center_eri'
!      call mat_print(Kx(3),1,nbast,1,nbast,lupri)

!      call mat_init(Jx(1),nbast,nbast)
!      call mat_init(Jx(2),nbast,nbast)
!      call mat_init(Jx(3),nbast,nbast)
!      call II_get_magderivJ_4center_eri(LUPRI,LUERR,ls%setting,nbast,NONsym,Jx) 
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym JxJx magderiv',D4*mat_trab(Jx(1),Jx(1))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym JxJy magderiv',D4*mat_trab(Jx(1),Jx(2))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym JxJz magderiv',D4*mat_trab(Jx(1),Jx(3))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym JyJy magderiv',D4*mat_trab(Jx(2),Jx(2))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym JyJz magderiv',D4*mat_trab(Jx(2),Jx(3))
!      write(lupri,'(A,F16.8)') 'QQQ1 NONsym JzJz magderiv',D4*mat_trab(Jx(3),Jx(3))
!      write(lupri,*)'Jx_4center_eri'
!      call mat_print(Jx(1),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Jy_4center_eri'
!      call mat_print(Jx(2),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Jz_4center_eri'
!      call mat_print(Jx(3),1,nbast,1,nbast,lupri)

!      DO I =1,3
!         call mat_free(Kx(I))
!         call mat_free(Jx(I))
!      ENDDO

      !magderivK + magderivJ
      call mat_init(Kx(1),nbast,nbast)
      call mat_init(Kx(2),nbast,nbast)
      call mat_init(Kx(3),nbast,nbast)
      call II_get_magderivK(LUPRI,LUERR,ls%setting,nbast,(/NONsym/),Kx) 
      write(lupri,'(A,F16.8)') 'QQQ2 NONsym KxKx magderiv',D4*mat_trab(Kx(1),Kx(1))
      write(lupri,'(A,F16.8)') 'QQQ2 NONsym KxKy magderiv',D4*mat_trab(Kx(1),Kx(2))
      write(lupri,'(A,F16.8)') 'QQQ2 NONsym kxKz magderiv',D4*mat_trab(Kx(1),Kx(3))
      write(lupri,'(A,F16.8)') 'QQQ2 NONsym KyKy magderiv',D4*mat_trab(Kx(2),Kx(2))
      write(lupri,'(A,F16.8)') 'QQQ2 NONsym KyKz magderiv',D4*mat_trab(Kx(2),Kx(3))
      write(lupri,'(A,F16.8)') 'QQQ2 NONsym KzKz magderiv',D4*mat_trab(Kx(3),Kx(3))

!      write(lupri,*)'Kx'
!      call mat_print(Kx(1),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Ky'
!      call mat_print(Kx(2),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Kz'
!      call mat_print(Kx(3),1,nbast,1,nbast,lupri)

!      call mat_init(Jx(1),nbast,nbast)
!      call mat_init(Jx(2),nbast,nbast)
!      call mat_init(Jx(3),nbast,nbast)
!      call II_get_magderivJ(LUPRI,LUERR,ls%setting,nbast,NONsym,Jx) 
!      write(lupri,'(A,F16.8)') 'QQQ2 NONsym JxJx magderiv',D4*mat_trab(Jx(1),Jx(1))
!      write(lupri,'(A,F16.8)') 'QQQ2 NONsym JxJy magderiv',D4*mat_trab(Jx(1),Jx(2))
!      write(lupri,'(A,F16.8)') 'QQQ2 NONsym JxJz magderiv',D4*mat_trab(Jx(1),Jx(3))
!      write(lupri,'(A,F16.8)') 'QQQ2 NONsym JyJy magderiv',D4*mat_trab(Jx(2),Jx(2))
!      write(lupri,'(A,F16.8)') 'QQQ2 NONsym JyJz magderiv',D4*mat_trab(Jx(2),Jx(3))
!      write(lupri,'(A,F16.8)') 'QQQ2 NONsym JzJz magderiv',D4*mat_trab(Jx(3),Jx(3))
!      write(lupri,*)'Jx'
!      call mat_print(Jx(1),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Jy'
!      call mat_print(Jx(2),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Jz'
!      call mat_print(Jx(3),1,nbast,1,nbast,lupri)

!      call mat_init(tempm1,nbast,nbast)
!      call mat_init(tempm2,nbast,nbast)
      DO I =1,3
!         call mat_add(1E0_realk,Jx(I),1E0_realk,Kx(I),tempm1)
!         call mat_add(1E0_realk,tempm1,-1E0_realk,Fx(I),tempm2)
!         WRITE(lupri,*)'DIFF NONsym F2=',mat_trab(tempm2,tempm2)
!         IF(ABS(mat_trab(tempm2,tempm2)).GT.1.0E-10)CALL LSQUIT('Magderiv Error',-1)
         call mat_free(Kx(I))
!         call mat_free(Jx(I))
      ENDDO
!      call mat_free(tempm2)
!      call mat_free(tempm1)

      !magderivF
!      call mat_init(Fx2(1),nbast,nbast)
!      call mat_init(Fx2(2),nbast,nbast)
!      call mat_init(Fx2(3),nbast,nbast)
!      call II_get_magderivF(LUPRI,LUERR,ls%setting,nbast,NONsym,Fx2) 
!      write(lupri,'(A,F16.8)') 'QQQ Fx2Fx2 magderiv',D4*mat_trab(Fx2(1),Fx2(1))
!      write(lupri,'(A,F16.8)') 'QQQ Fx2Fy2 magderiv',D4*mat_trab(Fx2(1),Fx2(2))
!      write(lupri,'(A,F16.8)') 'QQQ Fx2Fz2 magderiv',D4*mat_trab(Fx2(1),Fx2(3))
!      write(lupri,'(A,F16.8)') 'QQQ Fy2Fy2 magderiv',D4*mat_trab(Fx2(2),Fx2(2))
!      write(lupri,'(A,F16.8)') 'QQQ Fy2Fz2 magderiv',D4*mat_trab(Fx2(2),Fx2(3))
!      write(lupri,'(A,F16.8)') 'QQQ Fz2Fz2 magderiv',D4*mat_trab(Fx2(3),Fx2(3))

!      write(lupri,*)'Fx2'
!      call mat_print(Fx2(1),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Fy2'
!      call mat_print(Fx2(2),1,nbast,1,nbast,lupri)
!      write(lupri,*)'Fz2'
!      call mat_print(Fx2(3),1,nbast,1,nbast,lupri)

!      call mat_init(tempm2,nbast,nbast)
      DO I =1,3
!         call mat_add(1E0_realk,Fx2(I),-1E0_realk,Fx(I),tempm2)
!         WRITE(lupri,*)'DIFF NONsym F1=',mat_trab(tempm2,tempm2)
!         IF(ABS(mat_trab(tempm2,tempm2)).GT.1.0E-10)CALL LSQUIT('Magderiv Error NONsym',-1)
!         call mat_free(Fx2(I))
         call mat_free(Fx(I))
      ENDDO
!      call mat_free(tempm2)

      call mat_free(NONsym)

   END SUBROUTINE di_debug_magderiv_4center_eri

    subroutine di_ep_test(setting,D,lupri,luerr)
      implicit none
      type(lssetting)                :: setting
      integer                        :: lupri,luerr
      type(matrix)                   :: D
      !
      type(matrix)                   :: EPint
      integer                        :: I,N
      real(realk)                    :: R(3),output(1),t1,t2,output1
      logical                        :: TEST1,TEST2,TEST3
      real(realk),pointer            :: R2(:,:),output2(:)
      real(realk),parameter          :: D2=2E0_realk
      N=1
      call mat_init(EPint,D%nrow,D%ncol)
      R(1)=0.1233546d0
      R(2)=0.3465768d0
      R(3)=0.6768797d0
      call LSTIMER('START',t1,t2,LUPRI)
      call II_get_ep_integrals(LUPRI,LUERR,SETTING,EPint,R)
      call LSTIMER('EP1  ',t1,t2,LUPRI)
      output1 = D2*mat_trab(EPint,D)
      write(lupri,*)'output1',output1
      
      call LSTIMER('START',t1,t2,LUPRI)
      call II_get_ep_integrals2(LUPRI,LUERR,SETTING,R,D,output)
      call LSTIMER('EP2  ',t1,t2,LUPRI)
      write(lupri,*)'output(1)',output(1)
      
      call mem_alloc(R2,3,N)
      call mem_alloc(output2,N)
      R2(1,1)=0.1233546d0
      R2(2,1)=0.3465768d0
      R2(3,1)=0.6768797d0
      call II_get_ep_integrals3(LUPRI,LUERR,SETTING,R2,N,D,output2)
      write(lupri,*)'output2(1)',output2(1)

      TEST1 = ABS(output1-output(1)).GT.1E-12_realk
      TEST2 = ABS(output1-output2(1)).GT.1E-12_realk
      TEST3 = ABS(output(1)-output2(1)).GT.1E-12_realk

      IF(TEST1.OR.(TEST2.OR.TEST3))THEN
         print*,'EP 1 =',output1
         print*,'EP 2 =',output(1)
         print*,'EP 3 =',output2(1)
         CALL LSQUIT('EP error EP integrals are wrong',lupri)
      ELSE
         WRITE(lupri,*)'EP integrals calculated succesfully'
      ENDIF
      call mem_dealloc(output2)
      call mem_dealloc(R2)
      call mat_free(EPint)
      
    end subroutine di_ep_test

    subroutine di_screen_test(setting,nbast,lupri,luerr)
      implicit none
      integer,intent(in)      :: lupri,luerr,nbast
      type(lssetting)         :: setting
      !
      real(realk),pointer   :: integrals(:,:,:,:)
      real(realk),pointer   :: integrals2(:,:,:,:)
      integer :: n(7:15),test(7:15),iB,iA,iC,iD,I,J
      real(realk) :: CS(7:15),thresholdSAVE,int,THR
      real(realk),parameter :: CS7=1E-7_REALK,CS8=1E-8_REALK,CS9=1E-9_REALK  
      real(realk),parameter :: CS10=1E-10_REALK,CS11=1E-11_REALK,CS12=1E-12_REALK  
      real(realk),parameter :: CS13=1E-13_REALK,CS14=1E-14_REALK,CS15=1E-15_REALK  
      logical :: CS_screenSAVE, PS_screenSAVE, OE_screenSAVE, PARI_screenSAVE
      logical :: OD_screenSAVE, MBIE_screenSAVE,SAVEReCalcGab
      print*,'di_screen_test(setting,nbast,lupri,luerr)',nbast
      CS(7)=CS7; CS(8)=CS8; CS(9)=CS9; CS(10)=CS10; CS(11)=CS11; CS(12)=CS12;
      CS(13)=CS13; CS(14)=CS14; CS(15)=CS15;
      print*,'di_screen_test',CS(7:15)
      thresholdSAVE = setting%scheme%threshold
      CS_screenSAVE = setting%scheme%CS_screen
      PS_screenSAVE = setting%scheme%PS_screen
      OE_screenSAVE = setting%scheme%OE_screen
      PARI_screenSAVE = setting%scheme%PARI_screen
      OD_screenSAVE = setting%scheme%OD_screen
      MBIE_screenSAVE = setting%scheme%MBIE_screen

      setting%scheme%CS_screen = .FALSE.
      setting%scheme%PS_screen = .FALSE.
      setting%scheme%OE_screen = .FALSE.
      setting%scheme%PARI_screen = .FALSE.
      setting%scheme%OD_screen = .FALSE.
      setting%scheme%MBIE_screen = .FALSE.
      call mem_alloc(integrals,nbast,nbast,nbast,nbast)
      call II_get_4center_eri(LUPRI,LUERR,setting,integrals,nbast,nbast,nbast,nbast)
      n(7:15)=0
      DO iD=1,nbast
         DO iC=1,nbast
            DO iB=1,nbast
               DO iA=1,nbast
                  int = ABS(integrals(iA,iB,iC,iD))
                  DO J=7,15,3
                     IF(int.GT.CS(J)) n(J)=n(J)+1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO I=7,15,3
         WRITE(lupri,'(A,ES16.4,A,I10)')'The number of integrals greater than ',CS(I),' is ',n(I)
      ENDDO
      WRITE(lupri,*)'Testing CS screening'
      call mem_alloc(integrals2,nbast,nbast,nbast,nbast)
      setting%scheme%CS_screen = .TRUE.
      setting%scheme%PS_screen = .FALSE.
      setting%scheme%OE_screen = .FALSE.
      setting%scheme%PARI_screen = .FALSE.
      setting%scheme%OD_screen = .FALSE.
      setting%scheme%MBIE_screen = .FALSE.
      !WARNING the screening integrals have been calculated using OD screening
      !so we need to recalc the screening matrix
!      SAVEReCalcGab=setting%scheme%ReCalcGab
!      setting%scheme%ReCalcGab = .TRUE.
      DO I=7,15,3
         THR=CS(I)         
         WRITE(lupri,*)'Testing threshold = ',CS(I)
         setting%scheme%threshold = CS(I)/setting%scheme%J_THR
         WRITE(lupri,*)'Testing setting%scheme%threshold = ',setting%scheme%threshold
         call II_get_4center_eri(LUPRI,LUERR,setting,integrals2,nbast,nbast,nbast,nbast)
         test(7:15)=0;
         DO iD=1,nbast
            DO iC=1,nbast
               DO iB=1,nbast
                  DO iA=1,nbast
                     if(ABS(integrals(iA,iB,iC,iD)-integrals2(iA,iB,iC,iD)).GT.THR)then
                        WRITE(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8)')&
&'WITHOUT SCREENING Int(',ia,',',ib,',',ic,',',id,')=',integrals(iA,iB,iC,iD)
                        WRITE(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8)')&
&'WITH    SCREENING Int(',ia,',',ib,',',ic,',',id,')=',integrals2(iA,iB,iC,iD)
                        WRITE(lupri,*)'DIFF',ABS(integrals(iA,iB,iC,iD)-integrals2(iA,iB,iC,iD))
                        call ls_flshfo(lupri)
                        call lsquit('Error in di_screen_test CS elements wrong',-1)
                     endif
                     int = ABS(integrals(iA,iB,iC,iD))
                     DO J=7,15,3
                        IF(int.GT.CS(J)) test(J)=test(J)+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO J=7,I,3
            WRITE(lupri,'(A,ES16.4,A,I10)')'The number of integrals greater than ',CS(J),' is ',test(J)            
            if(test(J).NE.n(J)) call lsquit('Error in di_screen_test',lupri)
         ENDDO
      enddo      
      WRITE(lupri,*)'Testing PS screening'
      setting%scheme%CS_screen = .FALSE.
      setting%scheme%PS_screen = .TRUE.
      setting%scheme%OE_screen = .FALSE.
      setting%scheme%PARI_screen = .FALSE.
      setting%scheme%OD_screen = .FALSE.
      setting%scheme%MBIE_screen = .FALSE.
!      setting%scheme%ReCalcGab = .TRUE.
      DO I=7,15,3
         THR=CS(I)         
         WRITE(lupri,*)'Testing threshold = ',CS(I)
         setting%scheme%threshold = CS(I)/setting%scheme%J_THR
         WRITE(lupri,*)'Testing setting%scheme%threshold = ',setting%scheme%threshold
         call II_get_4center_eri(LUPRI,LUERR,setting,integrals2,nbast,nbast,nbast,nbast)
         test(7:15)=0;
         DO iD=1,nbast
            DO iC=1,nbast
               DO iB=1,nbast
                  DO iA=1,nbast
                     if(ABS(integrals(iA,iB,iC,iD)-integrals2(iA,iB,iC,iD)).GT.THR)then
                        WRITE(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8)')&
&'WITHOUT SCREENING Int(',ia,',',ib,',',ic,',',id,')=',integrals(iA,iB,iC,iD)
                        WRITE(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8)')&
&'WITH    SCREENING Int(',ia,',',ib,',',ic,',',id,')=',integrals2(iA,iB,iC,iD)
                        WRITE(lupri,*)'DIFF',ABS(integrals(iA,iB,iC,iD)-integrals2(iA,iB,iC,iD))
                        call lsquit('Error in di_screen_test PS elements wrong',-1)
                     endif
                     int = ABS(integrals(iA,iB,iC,iD))
                     DO J=7,15,3
                        IF(int.GT.CS(J)) test(J)=test(J)+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO J=7,I,3
            WRITE(lupri,'(A,ES16.4,A,I10)')'The number of integrals greater than ',CS(J),' is ',test(J)            
            if(test(J).NE.n(J)) call lsquit('Error in di_screen_test',lupri)
         ENDDO
      enddo      
      WRITE(lupri,*)'Testing ALL screening'
      setting%scheme%CS_screen = .TRUE.
      setting%scheme%PS_screen = .TRUE.
      setting%scheme%OE_screen = .TRUE.
      setting%scheme%PARI_screen = .TRUE.
      setting%scheme%OD_screen = .TRUE.
      setting%scheme%MBIE_screen = .FALSE.
!      setting%scheme%ReCalcGab = .TRUE.
      DO I=7,15,3
         THR=CS(I)         
         WRITE(lupri,*)'Testing threshold = ',CS(I)
         setting%scheme%threshold = CS(I)/setting%scheme%J_THR
         WRITE(lupri,*)'Testing setting%scheme%threshold = ',setting%scheme%threshold
         call II_get_4center_eri(LUPRI,LUERR,setting,integrals2,nbast,nbast,nbast,nbast)
         test(7:15)=0;
         DO iD=1,nbast
            DO iC=1,nbast
               DO iB=1,nbast
                  DO iA=1,nbast
                     if(ABS(integrals(iA,iB,iC,iD)-integrals2(iA,iB,iC,iD)).GT.THR)then
                        WRITE(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8)')&
&'WITHOUT SCREENING Int(',ia,',',ib,',',ic,',',id,')=',integrals(iA,iB,iC,iD)
                        WRITE(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8)')&
&'WITH    SCREENING Int(',ia,',',ib,',',ic,',',id,')=',integrals2(iA,iB,iC,iD)
                        call lsquit('Error in di_screen_test elements wrong',-1)
                     endif
                     int = ABS(integrals(iA,iB,iC,iD))
                     DO J=7,15,3
                        IF(int.GT.CS(J)) test(J)=test(J)+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO J=7,I,3
            WRITE(lupri,'(A,ES16.4,A,I10)')'The number of integrals greater than ',CS(J),' is ',test(J)            
            if(test(J).NE.n(J)) call lsquit('Error in di_screen_test',lupri)
         ENDDO
      enddo      
      setting%scheme%CS_screen = CS_screenSAVE 
      setting%scheme%PS_screen = PS_screenSAVE
      setting%scheme%OE_screen = OE_screenSAVE
      setting%scheme%PARI_screen = PARI_screenSAVE 
      setting%scheme%OD_screen = OD_screenSAVE  
      setting%scheme%MBIE_screen = MBIE_screenSAVE
      setting%scheme%threshold = thresholdSAVE
!      setting%scheme%ReCalcGab = SAVEReCalcGab
      call mem_dealloc(integrals)
      call mem_dealloc(integrals2)
      WRITE(lupri,*)'di_screen_test SUCCESSFULL'

    end subroutine di_screen_test

    subroutine di_geoderivOverlap_test(ls,D,lupri,luerr)
      implicit none
      type(lsitem)                   :: ls
      integer                        :: lupri,luerr
      type(matrix),target            :: D
      !
      type(matrix),pointer           :: Sx(:),SLx(:)
      type(matrix)                   :: tmp1,tmp2
      type(matrixp)                  :: Dp(1)
      integer                        :: I,NATOMS,nbast,J
      real(realk),pointer            :: reOrtho(:,:)
      nbast=D%nrow
      Dp(1)%p => D
      NATOMS = ls%input%MOLECULE%nAtoms
      call mem_alloc(reOrtho,3,natoms)
      call II_get_reorthoNormalization(reOrtho,Dp,1,ls%setting,lupri,luerr)

      allocate(Sx(3*Natoms))
      do I=1,3*Natoms
         CALL mat_init(Sx(I),nbast,nbast)
      ENDDO
      write(lupri,*)'II_get_geoderivOverlap'
      call II_get_geoderivOverlap(Sx,natoms,ls%setting,lupri,luerr)
      do I=1,Natoms
         IF(ABS(mat_trAB(Sx(3*(I-1)+1),D)-reOrtho(1,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(Sx(3*(I-1)+1),D)',mat_trAB(Sx(3*(I-1)+1),D)
            print*,'reOrtho(1,I)',reOrtho(1,I)
            print*,'DIFF:',ABS(mat_trAB(Sx(3*(I-1)+1),D)-reOrtho(1,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivOverlap_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(Sx(3*(I-1)+2),D)-reOrtho(2,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(Sx(3*(I-1)+1),D)',mat_trAB(Sx(3*(I-1)+2),D)
            print*,'reOrtho(1,I)',reOrtho(2,I)
            print*,'DIFF:',ABS(mat_trAB(Sx(3*(I-1)+2),D)-reOrtho(2,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivOverlap_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(Sx(3*(I-1)+3),D)-reOrtho(3,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(Sx(3*(I-1)+1),D)',mat_trAB(Sx(3*(I-1)+3),D)
            print*,'reOrtho(1,I)',reOrtho(3,I)
            print*,'DIFF:',ABS(mat_trAB(Sx(3*(I-1)+3),D)-reOrtho(3,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivOverlap_test',lupri)
         ENDIF
      enddo
      allocate(SLx(3*Natoms))
      do I=1,3*Natoms
         CALL mat_init(SLx(I),nbast,nbast)
      ENDDO
      call II_get_geoderivOverlapL(SLx,natoms,ls%setting,lupri,luerr)
      call II_get_reorthoNormalization2(reOrtho,Dp,1,ls%setting,lupri,luerr)
      do I=1,Natoms
         IF(ABS(mat_trAB(SLx(3*(I-1)+1),D)-reOrtho(1,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(SLx(3*(I-1)+1),D)',mat_trAB(SLx(3*(I-1)+1),D)
            print*,'reOrtho(1,I)',reOrtho(1,I)
            print*,'DIFF:',ABS(mat_trAB(SLx(3*(I-1)+1),D)-reOrtho(1,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivOverlapL_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(SLx(3*(I-1)+2),D)-reOrtho(2,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(SLx(3*(I-1)+1),D)',mat_trAB(SLx(3*(I-1)+2),D)
            print*,'reOrtho(1,I)',reOrtho(2,I)
            print*,'DIFF:',ABS(mat_trAB(SLx(3*(I-1)+2),D)-reOrtho(2,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivOverlapL_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(SLx(3*(I-1)+3),D)-reOrtho(3,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(SLx(3*(I-1)+1),D)',mat_trAB(SLx(3*(I-1)+3),D)
            print*,'reOrtho(1,I)',reOrtho(3,I)
            print*,'DIFF:',ABS(mat_trAB(SLx(3*(I-1)+3),D)-reOrtho(3,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivOverlapL_test',lupri)
         ENDIF
      enddo
      CALL mat_init(tmp1,nbast,nbast) !SLx^transposed
      CALL mat_init(tmp2,nbast,nbast)
      do I=1,3*Natoms
         call mat_trans(SLx(I), tmp1)     !tmp1 = SLx^T
         call mat_assign(tmp2,SLx(I))     !tmp2 = SLx
         call mat_daxpy(1.0E0_realk, tmp1, tmp2)      !tmp2 = SLx + SLx^T
         call mat_add(1E0_realk,tmp2,-1E0_realk,Sx(I),tmp1)
         WRITE(lupri,*)'DIFF1=',mat_trab(tmp1,tmp1)
         IF(ABS(mat_trab(tmp1,tmp1)).GT.1.0E-10)CALL LSQUIT('GeoderivL Error',-1)
      ENDDO
      CALL mat_free(tmp1)
      CALL mat_free(tmp2)
      do I=1,3*Natoms
         CALL mat_free(Sx(I))
      ENDDO
      deallocate(Sx)
      do I=1,3*Natoms
         CALL mat_free(SLx(I))
      ENDDO
      deallocate(SLx)
      call mem_dealloc(reOrtho)
      WRITE(lupri,*)'di_geoderivOverlap_test Successful'
    end subroutine di_geoderivOverlap_test

    subroutine di_geoderivKinetic_test(ls,D,lupri,luerr)
      implicit none
      type(lsitem)                   :: ls
      integer                        :: lupri,luerr
      type(matrix),target            :: D
      !
      type(matrix),pointer           :: hx(:)
      type(matrix)                   :: tmp1,tmp2
      type(matrixp)                  :: Dp(1)
      integer                        :: I,NATOMS,nbast,J
      real(realk),pointer            :: grad(:,:)
      nbast=D%nrow
      Dp(1)%p => D
      NATOMS = ls%input%MOLECULE%nAtoms
      call mem_alloc(grad,3,natoms)
      call mem_alloc(hx,3*Natoms)
      do I=1,3*Natoms
         CALL mat_init(hx(I),nbast,nbast)
      ENDDO
      CALL II_get_kinetic_gradient(Grad,Dp,1,ls%setting,lupri,luerr)
      call II_get_geoderivKinetic(hx,natoms,ls%setting,lupri,luerr)
      do I=1,Natoms
         IF(ABS(mat_trAB(hx(3*(I-1)+1),D)-grad(1,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(hx(3*(I-1)+1),D)',mat_trAB(hx(3*(I-1)+1),D)
            print*,'grad(1,I)',grad(1,I)
            print*,'DIFF:',ABS(mat_trAB(hx(3*(I-1)+1),D)-grad(1,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivKinetic_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(hx(3*(I-1)+2),D)-grad(2,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(hx(3*(I-1)+1),D)',mat_trAB(hx(3*(I-1)+2),D)
            print*,'grad(1,I)',grad(2,I)
            print*,'DIFF:',ABS(mat_trAB(hx(3*(I-1)+2),D)-grad(2,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivKinetic_test',lupri)
         ENDIF
         IF(ABS(mat_trAB(hx(3*(I-1)+3),D)-grad(3,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(hx(3*(I-1)+1),D)',mat_trAB(hx(3*(I-1)+3),D)
            print*,'grad(1,I)',grad(3,I)
            print*,'DIFF:',ABS(mat_trAB(hx(3*(I-1)+3),D)-grad(3,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivKinetic_test',lupri)
         ENDIF
      enddo
      WRITE(lupri,*)'di_geoderivKinetic_test Successful'
      do I=1,3*Natoms
         CALL mat_free(hx(I))
      ENDDO
      call mem_dealloc(grad)
      call mem_dealloc(hx)
      
    end subroutine di_geoderivKinetic_test

    subroutine di_geoderivNucel_test(ls,D,lupri,luerr)
      implicit none
      type(lsitem)                   :: ls
      integer                        :: lupri,luerr
      type(matrix),target            :: D
      !
      type(matrix),pointer           :: hx(:)
      type(matrix)                   :: tmp1,tmp2
      type(matrixp)                  :: Dp(1)
      integer                        :: I,NATOMS,nbast,J
      real(realk),pointer            :: grad(:,:)
      nbast=D%nrow
      Dp(1)%p => D
      NATOMS = ls%input%MOLECULE%nAtoms
      call mem_alloc(grad,3,natoms)
      call mem_alloc(hx,3*Natoms)
      do I=1,3*Natoms
         CALL mat_init(hx(I),nbast,nbast)
      ENDDO
      CALL II_get_ne_gradient(Grad,Dp,1,ls%setting,lupri,luerr)
      call II_get_geoderivnucel(hx,natoms,ls%setting,lupri,luerr)
      do I=1,Natoms
         IF(ABS(mat_trAB(hx(3*(I-1)+1),D)-grad(1,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(hx(3*(I-1)+1),D)',mat_trAB(hx(3*(I-1)+1),D)
            print*,'grad(1,I)',grad(1,I)
            print*,'DIFF:',ABS(mat_trAB(hx(3*(I-1)+1),D)-grad(1,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivNucel_testX',lupri)
         ENDIF
         IF(ABS(mat_trAB(hx(3*(I-1)+2),D)-grad(2,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(hx(3*(I-1)+1),D)',mat_trAB(hx(3*(I-1)+2),D)
            print*,'grad(1,I)',grad(2,I)
            print*,'DIFF:',ABS(mat_trAB(hx(3*(I-1)+2),D)-grad(2,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivNucel_testY',lupri)
         ENDIF
         IF(ABS(mat_trAB(hx(3*(I-1)+3),D)-grad(3,I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(hx(3*(I-1)+1),D)',mat_trAB(hx(3*(I-1)+3),D)
            print*,'grad(1,I)',grad(3,I)
            print*,'DIFF:',ABS(mat_trAB(hx(3*(I-1)+3),D)-grad(3,I)).GT.1E-12_realk
            CALL LSQUIT('ERROR IN di_geoderivNucel_testZ',lupri)
         ENDIF
      enddo
      WRITE(lupri,*)'di_geoderivNucel_test Successful'
      do I=1,3*Natoms
         CALL mat_free(hx(I))
      ENDDO
      call mem_dealloc(hx)
      call mem_dealloc(grad)
    end subroutine di_geoderivNucel_test

    subroutine di_geoderivExchange_test(ls,D,lupri,luerr)
      implicit none
      type(lsitem)                   :: ls
      integer                        :: lupri,luerr
      type(matrix),target            :: D
      !
      type(matrix),pointer           :: Kx(:)
      type(matrixp)                  :: Dp(1)
      integer                        :: I,NATOMS,nbast,J
      real(realk),pointer            :: grad(:,:)
      nbast=D%nrow
      Dp(1)%p => D
      NATOMS = ls%input%MOLECULE%nAtoms
      call mem_alloc(grad,3,natoms)
      call II_get_K_gradient(grad,Dp,Dp,1,1,ls%setting,lupri,luerr)

      allocate(Kx(3*Natoms))
      do I=1,3*Natoms
         CALL mat_init(Kx(I),nbast,nbast)
      ENDDO
      call II_get_geoderivExchange(Kx,D,natoms,ls%setting,lupri,luerr)
      do I=1,Natoms
         IF(ABS(-0.5E0_realk*mat_trAB(Kx(3*(I-1)+1),D)-grad(1,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(Kx(3*(I-1)+1),D)',-0.5E0_realk*mat_trAB(Kx(3*(I-1)+1),D)
            print*,'grad(1,I)',grad(1,I)
            print*,'DIFF:',ABS(mat_trAB(Kx(3*(I-1)+1),D)-grad(1,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN X di_geoderivExchange_test',lupri)
         ENDIF
         IF(ABS(-0.5E0_realk*mat_trAB(Kx(3*(I-1)+2),D)-grad(2,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(Kx(3*(I-1)+1),D)',-0.5E0_realk*mat_trAB(Kx(3*(I-1)+2),D)
            print*,'grad(1,I)',grad(2,I)
            print*,'DIFF:',ABS(mat_trAB(Kx(3*(I-1)+2),D)-grad(2,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN Y di_geoderivExchange_test',lupri)
         ENDIF
         IF(ABS(-0.5E0_realk*mat_trAB(Kx(3*(I-1)+3),D)-grad(3,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(Kx(3*(I-1)+1),D)',-0.5E0_realk*mat_trAB(Kx(3*(I-1)+3),D)
            print*,'grad(1,I)',grad(3,I)
            print*,'DIFF:',ABS(mat_trAB(Kx(3*(I-1)+3),D)-grad(3,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN Z di_geoderivExchange_test',lupri)
         ENDIF
      enddo
      do I=1,3*Natoms
         CALL mat_free(Kx(I))
      ENDDO
      call mem_dealloc(grad)
      WRITE(lupri,*)'di_geoderivExchange_test Successful'
    end subroutine di_geoderivExchange_test

    subroutine di_geoderivCoulomb_test(ls,D,lupri,luerr)
      implicit none
      type(lsitem)                   :: ls
      integer                        :: lupri,luerr
      type(matrix),target            :: D
      !
      type(matrix),pointer           :: Jx(:)
      type(matrixp)                  :: Dp(1)
      integer                        :: I,NATOMS,nbast,J
      real(realk),pointer            :: grad(:,:)
      nbast=D%nrow
      Dp(1)%p => D
      NATOMS = ls%input%MOLECULE%nAtoms
      call mem_alloc(grad,3,natoms)
      call II_get_J_gradient(grad,Dp,Dp,1,1,ls%setting,lupri,luerr)

      allocate(Jx(3*Natoms))
      do I=1,3*Natoms
         CALL mat_init(Jx(I),nbast,nbast)
      ENDDO
      call II_get_geoderivCoulomb(Jx,D,natoms,ls%setting,lupri,luerr)
      do I=1,Natoms
         IF(ABS(0.5E0_realk*mat_trAB(Jx(3*(I-1)+1),D)-grad(1,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(Jx(3*(I-1)+1),D)',-0.5E0_realk*mat_trAB(Jx(3*(I-1)+1),D)
            print*,'grad(1,I)',grad(1,I)
            print*,'DIFF:',ABS(mat_trAB(Jx(3*(I-1)+1),D)-grad(1,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN X di_geoderivCoulomb_test',lupri)
         ENDIF
         IF(ABS(0.5E0_realk*mat_trAB(Jx(3*(I-1)+2),D)-grad(2,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(Jx(3*(I-1)+1),D)',-0.5E0_realk*mat_trAB(Jx(3*(I-1)+2),D)
            print*,'grad(1,I)',grad(2,I)
            print*,'DIFF:',ABS(mat_trAB(Jx(3*(I-1)+2),D)-grad(2,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN Y di_geoderivCoulomb_test',lupri)
         ENDIF
         IF(ABS(0.5E0_realk*mat_trAB(Jx(3*(I-1)+3),D)-grad(3,I)).GT.1E-8_realk)THEN
            print*,'mat_trAB(Jx(3*(I-1)+1),D)',-0.5E0_realk*mat_trAB(Jx(3*(I-1)+3),D)
            print*,'grad(1,I)',grad(3,I)
            print*,'DIFF:',ABS(mat_trAB(Jx(3*(I-1)+3),D)-grad(3,I)).GT.1E-8_realk
            CALL LSQUIT('ERROR IN Z di_geoderivCoulomb_test',lupri)
         ENDIF
      enddo
      do I=1,3*Natoms
         CALL mat_free(Jx(I))
      ENDDO
      call mem_dealloc(grad)
      WRITE(lupri,*)'di_geoderivCoulomb_test Successful'
    end subroutine di_geoderivCoulomb_test

    subroutine di_magderivOverlap_test(ls,D,lupri,luerr)
      implicit none
      type(lsitem)                   :: ls
      integer,intent(in)             :: lupri,luerr
      type(matrix),intent(in)        :: D
      !
      type(matrix),pointer           :: Sx(:)
      integer                        :: I,nbast,J
      type(matrix)                   :: D2(1)
      real(realk)                    :: maggrad(3)

      nbast=D%nrow
      allocate(Sx(3))
      do I=1,3
         CALL mat_init(Sx(I),nbast,nbast)
      ENDDO
      call II_get_magderivOverlap(Sx,ls%setting,lupri,luerr)
      write(lupri,*) 'QQQ SxSx magderiv',mat_trab(Sx(1),Sx(1))
      write(lupri,*) 'QQQ SxSy magderiv',mat_trab(Sx(1),Sx(2))
      write(lupri,*) 'QQQ SxSz magderiv',mat_trab(Sx(1),Sx(3))
      write(lupri,*) 'QQQ SySy magderiv',mat_trab(Sx(2),Sx(2))
      write(lupri,*) 'QQQ SySz magderiv',mat_trab(Sx(2),Sx(3))
      write(lupri,*) 'QQQ SzSz magderiv',mat_trab(Sx(3),Sx(3))
      call II_get_maggradOverlap(maggrad,Sx(1:1),1,ls%setting,lupri,luerr)
      do I=1,3
         IF(ABS(mat_trAB(Sx(I),Sx(1))-maggrad(I)).GT.1E-12_realk)THEN
            print*,'mat_trAB(Sx(I),Sx(1))',mat_trAB(Sx(I),Sx(1))
            print*,'maggrad(I)',maggrad(I)
            print*,'DIFF:',ABS(mat_trAB(Sx(I),Sx(1))-maggrad(I))
            CALL LSQUIT('ERROR IN di_magderivOverlap_test',lupri)
         ENDIF
         print*,'mat_trAB(Sx(I),Sx(1))',mat_trAB(Sx(I),Sx(1))
         print*,'maggrad(I)',maggrad(I)
         print*,'DIFF:',ABS(mat_trAB(Sx(I),Sx(1))-maggrad(I))
      ENDDO
      do I=1,3
!         WRITE(lupri,*)'Sx(I) I=',I
!         call mat_print(Sx(I),1,nbast,1,nbast,lupri)
         CALL mat_free(Sx(I))         
      ENDDO
      deallocate(Sx)
      WRITE(lupri,*)'di_magderivOverlap_test Successful'
    end subroutine di_magderivOverlap_test

   SUBROUTINE di_debug_ccfragment(lupri,luerr,ls,nbast,D,CMO,nocc)
      IMPLICIT NONE
      TYPE(Matrix),intent(in) :: D,CMO
      integer,intent(in)      :: lupri,luerr,nbast,nocc
      type(lsitem),intent(inout) :: ls
      !
      type(lsitem),pointer :: fragment_ls(:)
      type(AOITEM) :: AO
      TYPE(Matrix)  :: tempm1,tempm2,tempm3,tempm4,D2,refF
      real(realk),pointer   :: densPtmp(:,:)
      integer,pointer :: orb2batch(:),orb2atom(:),ATOMS(:)
      logical,pointer :: LATOMS(:)
      integer :: iatom,nbast1,nbast2,I,J,dim1,dim2,dim3,dim4,igamma
      integer :: delta,iorb,natoms,nBatches,natom2,iprint
      real(realk) :: MAXCMO
      iprint=0
      WRITE(LUPRI,*)'DEBUGGING THE CCFRAGMENT'
      WRITE(*,*)'DEBUGGING THE CCFRAGMENT'

      !build the coulomb matrix 
      call mat_init(tempm1,nbast,nbast)
      call mat_init(tempm2,nbast,nbast)
      call mat_init(tempm3,nbast,nbast)
      call mat_init(refF,nbast,nbast)
      CALL II_get_coulomb_mat(LUPRI,LUERR,ls%setting,D,refF,1)

      natoms = ls%input%molecule%natoms
      nullify(LATOMS)
      nullify(ATOMS)
      allocate(LATOMS(natoms))
      ALLOCATE(ATOMS(natoms))

      nullify(orb2batch)
      allocate(orb2batch(nbast))
      CALL II_getBatchOrbitalInfo(ls%setting,&
           &nBast,orb2Batch,nBatches,lupri,luerr)

      call BUILD_AO(LUPRI,ls%setting%SCHEME,0,ls%setting%MOLECULE(1)%p,ls%setting%BASIS(1)%p%REGULAR,AO&
     &,.false.,.false.)

      call mem_alloc(orb2atom,nbast)
      do iorb = 1,nbast
         orb2atom(iorb) = AO%BATCH(orb2batch(iorb))%atom
      enddo

      dim1=nbast
      dim2=nbast

      nullify(fragment_ls)
      allocate(fragment_ls(nocc))
      call mat_zero(tempm3)
      MAXCMO = 0E0_realk
      DO iorb = 1,nocc
         do igamma=1,nbast
            MAXCMO = MAX(MAXCMO,CMO%elms(igamma+(iorb-1)*nbast))
         enddo
         LATOMS = .FALSE. 
         do igamma=1,nbast
            IF(ABS(MAXCMO*CMO%elms(igamma+(iorb-1)*nbast)).GT. 1E-13_realk)THEN
               LATOMS(orb2atom(igamma)) = .TRUE. 
            ENDIF
         enddo
!         print*,'LATOMS',LATOMS
         nbast1 = 0
         do igamma=1,nbast
            IF(LATOMS(orb2atom(igamma)))THEN
               nbast1 = nbast1+1
            ENDIF
         enddo
         dim3=nbast1
         dim4=nbast1
!         print*,'dim3',dim3,'dim4',dim4
         natom2=0
         ATOMS=0
         DO iatom = 1,natoms
            IF(LATOMS(iatom))THEN
               natom2=natom2+1
               ATOMS(natom2)=iatom
            ENDIF
         ENDDO
!         print*,'ATOMS',ATOMS
         nbast1 = 0
         nullify(DensPtmp)
         allocate(DensPtmp(dim3,dim4))
         DensPtmp = 0E0_realk
         do igamma=1,nbast
            if(LATOMS(orb2atom(igamma)))THEN
               nbast1 = nbast1+1
               nbast2 = 0
               do delta=1,nbast
                  if(LATOMS(orb2atom(delta)))THEN
                     nbast2 = nbast2+1
                     DensPtmp(nbast1,nbast2) = CMO%elms(igamma+(iorb-1)*nbast)*CMO%elms(delta+(iorb-1)*nbast)
                  endif
               enddo
            endif
         enddo
         call build_ccfragmentlsitem(LS,fragment_ls(iorb),ATOMS(1:NATOM2),NATOM2,LUPRI,IPRINT)
         ! here we set the ab in (ab|cd) to full
         fragment_ls(iorb)%SETTING%MOLECULE(1)%p => LS%SETTING%MOLECULE(1)%p
         fragment_ls(iorb)%SETTING%BASIS(1)%p => LS%SETTING%BASIS(1)%p
         fragment_ls(iorb)%SETTING%FRAGMENT(1)%p => LS%SETTING%FRAGMENT(1)%p
         fragment_ls(iorb)%SETTING%MOLECULE(2)%p => LS%SETTING%MOLECULE(2)%p
         fragment_ls(iorb)%SETTING%BASIS(2)%p => LS%SETTING%BASIS(2)%p
         fragment_ls(iorb)%SETTING%FRAGMENT(2)%p => LS%SETTING%FRAGMENT(2)%p
         fragment_ls(iorb)%SETTING%molID(1)=0
         fragment_ls(iorb)%SETTING%molID(2)=0
         fragment_ls(iorb)%SETTING%molID(3)=iorb
         fragment_ls(iorb)%SETTING%molID(4)=iorb
         fragment_ls(iorb)%SETTING%sameMOL=.TRUE.
         DO I=1,2
            DO J=3,4          
               fragment_ls(iorb)%SETTING%sameMOL(I,J) = .FALSE.
               fragment_ls(iorb)%SETTING%sameMOL(J,I) = .FALSE.
            ENDDO
         ENDDO
         call mat_init(D2,dim3,dim4)
         CALL mat_set_from_full(DensPtmp,1E0_realk, D2,'DenP')
         deallocate(DensPtmp)
         call mat_zero(tempm1)
         CALL II_get_coulomb_mat(LUPRI,LUERR,fragment_ls(iorb)%setting,D2,tempm1,1)
         call mat_free(D2)

!         allocate(DensPtmp(nbast,nbast))
!         DensPtmp=0E0_realk
!         do delta=1,nbast
!            do gamma=1,nbast
!               DensPtmp(gamma,delta) = CMO%elms(gamma+(iorb-1)*nbast)*CMO%elms(delta+(iorb-1)*nbast)
!            enddo
!         enddo
!         call mat_init(D2,nbast,nbast)
!         CALL mat_set_from_full(DensPtmp,1E0_realk, D2,'DenP')
!         call mat_zero(tempm2)
!         CALL II_get_coulomb_mat(LUPRI,LUERR,ls%setting,D2,tempm2)
!         call mat_init(tempm4,nbast,nbast)
!         call mat_add(1E0_realk,tempm1,-1E0_realk,tempm2,tempm4)
!         write(lupri,*) 'QQQ DI_DEBUG',iorb,'             :',mat_trab(tempm2,tempm2)
!         write(lupri,*) 'QQQ DI_DEBUG lsfragment',iorb,'  :',mat_trab(tempm1,tempm1)
!         write(lupri,*) 'QQQ DIFF                         :',mat_trab(tempm4,tempm4)
!         IF(ABS(mat_trab(tempm4,tempm4)).LE. 1E-8_realk)THEN
!            write(lupri,*)'QQQ SUCCESFUL'
!         ELSE
!            WRITE(lupri,*)'ZZ THE FULL DMAT IORB =',iorb
!            call mat_print(D2,1,nbast,1,nbast,lupri)
!            WRITE(lupri,*)'THE FRAGMENT FOCK'
!            call mat_print(tempm1,1,nbast,1,nbast,lupri)
!            WRITE(lupri,*)'THE LS FOCK'
!            call mat_print(tempm2,1,nbast,1,nbast,lupri)
!            WRITE(lupri,*)'THE DIFF'
!            call mat_print(tempm4,1,nbast,1,nbast,lupri)            
!            WRITE(LUPRI,*)'ZZZ LS SETTING'
!            call PRINT_LSSETTING(ls%setting,LUPRI)
!            WRITE(LUPRI,*)'QQQ FRAGMNET SETTING'
!            call PRINT_LSSETTING(fragment(iorb)%setting,LUPRI)
!            
!            CALL lsQUIT('TEST TEST IORB')
!         ENDIF
!         call mat_free(D2)
!         call mat_free(tempm4)

         call ls_free(fragment_ls(iorb))
         !** F = h + G
         call mat_daxpy(1E0_realk,tempm1,tempm3)
      ENDDO

      call mat_init(tempm4,nbast,nbast)
      call mat_add(1E0_realk,tempm3,-1E0_realk,refF,tempm4)
      write(lupri,*) 'QQQ DI_DEBUG     STD   ',mat_trab(refF,refF)
      write(lupri,*) 'QQQ DI_DEBUG lsfragment',mat_trab(tempm3,tempm3)
      write(lupri,*) 'QQQ DIFF               ',mat_trab(tempm4,tempm4)
      IF(ABS(mat_trab(tempm4,tempm4)).LE. 1E-8_realk)THEN
         write(lupri,*)'QQQ SUCCESFUL'
      ELSE
         WRITE(lupri,*)'THE DIFF'
         call mat_print(tempm4,1,nbast,1,nbast,lupri)
         CALL lsQUIT('TEST TEST',lupri)
      ENDIF
      call mem_dealloc(orb2atom)
      call mat_free(tempm1)
      call mat_free(tempm2)
      call mat_free(tempm3)
      call mat_free(tempm4)
      call mat_free(refF)
      call free_aoitem(lupri,AO)

    END SUBROUTINE di_debug_ccfragment

   SUBROUTINE di_get_overlap_and_H1(S, H1,LSint,newlupri,newluerr,ls)
      IMPLICIT NONE
      TYPE(Matrix), target :: S, H1
      TYPE(Matrix) :: tempm1,tempm2,Kinetic
      integer  :: NNBAST, n2basx,newlupri,newluerr
      type(lsitem) :: ls
      logical :: LFOUND,LSint
      real(realk), ALLOCATABLE :: wrk(:), wrk1(:)
      !local pointer to H1
      LFOUND = .FALSE.
      lH1 => H1
      call II_get_overlap(newlupri,newluerr,ls%setting,S)
      call II_get_h1(newlupri,newluerr,ls%setting,H1)

   END SUBROUTINE di_get_overlap_and_H1

   SUBROUTINE di_get_fock_LSDALTON(D,h1,F,Etotal,lupri,luerr,ls)
   ! ===================================================================
   ! di_get_fock obtains total fock matrix and corresponding energy.
   ! WE have to go through the interface to dalton before the fock
   ! evaluator learns how to handle arbitrary-type arrays.
   ! ===================================================================
      IMPLICIT NONE
      TYPE(Matrix), INTENT(IN)    :: H1, D
      TYPE(Matrix), INTENT(INOUT) :: F
      type(lsitem) :: ls
      real(realk), INTENT(INOUT) :: Etotal
      !
      real(realk)   :: edfty(1),fac,hfweight,EdXC
      integer nbast, lupri,luerr,ndmat
      logical :: Dsym,ADMMexchange
      TYPE(Matrix) :: K,dXC

      nbast = D%nrow
      fac = 2E0_realk
      IF(matrix_type .EQ. mtype_unres_dense)fac = 1E0_realk
      Dsym = .TRUE. !symmetric Density matrix
      ls%input%nfock = ls%input%nfock + 1

!     Turn of ADMM at level 2 for ADMM_GCBASIS option (we then use the level 2
!     basis as ADMM basis for level 3)
      ADMMexchange = ls%input%dalton%ADMM_EXCHANGE.AND.(.NOT.ls%optlevel.EQ.1)
      IF ((ls%optlevel.EQ.2).AND.ls%input%dalton%ADMM_GCBASIS) ADMMexchange = .FALSE.

! *********************************************************************************
! *                       Fock matrix with ADMM exchange
! *********************************************************************************
      IF (ADMMexchange) THEN

         !FixMe Should also work for incremental scheme
         IF(incremental_scheme)THEN
            call lsquit('Auxiliary Density Matrix Calculation requires NOINCREM',-1)
         ENDIF
         IF (ls%input%dalton%ADMM_GCBASIS.AND.(ls%input%basis%valence%nAtomtypes.EQ.0)) THEN
           call lsquit('Auxiliary Density Matrix GC-basis type Calculation requires TRILEVEL start guess',-1)
         ENDIF


         !We do the full Coulomb part with Full density
         call mat_zero(F)
         ndmat = 1
         call II_get_coulomb_mat(LUPRI,LUERR,ls%SETTING,D,F,ndmat)  
         WRITE(lupri,*)'The Coulomb energy contribution ',fac*0.5E0_realk*mat_dotproduct(D,F)


         !Then we build the ADMM exact exchange matrix K and the XC correction dXC
         call mat_init(K,nbast,nbast)
         call mat_init(dXC,nbast,nbast)

         call mat_zero(K)
         call mat_zero(dXC)

         CALL II_get_admm_exchange_mat(LUPRI,LUERR,ls%SETTING,ls%optlevel,D,K,dXC,ndmat,EdXC,dsym)
         call mat_daxpy(1.E0_realk,K,F)
         Etotal = fockenergy_f(F,D,H1,ls%input%dalton%unres,ls%input%potnuc,lupri)+EdXC
         call mat_daxpy(1.E0_realk,dXC,F)

         call mat_free(K)
         call mat_free(dXC)

! *********************************************************************************
! *                              Regular case          
! *********************************************************************************
      ELSE
         ndmat = 1
         call II_get_Fock_mat(lupri,luerr,ls%setting,D,Dsym,F,ndmat,incremental_scheme)      
         Etotal = fockenergy_f(F,D,H1,ls%input%dalton%unres,ls%input%potnuc,lupri)
      ENDIF
      IF(ls%setting%do_dft) THEN
         ndmat = 1
         nbast = D%nrow
         call II_get_xc_fock_mat(lupri,luerr,ls%setting,nbast,D,F,Edfty,ndmat)
         Etotal = Etotal + Edfty(1)
      ENDIF
      !** F = h + G
      call mat_daxpy(1E0_realk,H1,F)
   END SUBROUTINE di_get_fock_LSDALTON

   real(realk) function fockenergy_F(F,D,H1,unres,pot_nuc,lupri)
     !E = Tr(h + F)D + hnuc ! No factor since molecular orbitals
     implicit none
     TYPE(matrix), intent(in) :: F,D,H1
     real(realk) :: hD, FD, fac,pot_nuc
     integer :: ndim,lupri
     logical :: unres
     
     ndim=F%nrow
     fac=2E0_realk
     if(unres) fac=1E0_realk
     !Get the one-electron hamiltonian
     !Tr(FD)
     !Tr(hD)
     hD =       mat_dotproduct(D,H1)
     FD = 0.5E0_realk*mat_dotproduct(D,F)
     !E(HF) = Tr(h +F)D + hnuc
     fockenergy_F = (hd + FD)*fac + POT_NUC !Stinne: Remove call to boxed DF contribution (obsolete), 08-06-2010
     ! FIXME + potnuc
   end function fockenergy_F

   SUBROUTINE di_SCF_EnergyCont(D,h1,energy,ndmat,newlupri,newluerr,ls,modthresh)
      IMPLICIT NONE
      integer :: newlupri,newluerr,ndmat
      TYPE(Matrix), INTENT(IN)    :: H1, D(ndmat)
      real(realk), INTENT(INOUT) :: energy(ndmat)
      type(lsitem) :: ls
      real(realk), INTENT(IN) :: modthresh
      !
      integer :: nbast,idmat
      logical :: Dsym
      real(realk) :: hD,TwoElectronEcont(ndmat),Edfty(ndmat)
      real(realk) :: fac
      
      Dsym = .TRUE. !symmetric Density matrix
      fac=2E0_realk
      if(ls%input%dalton%unres) fac=1E0_realk

      call II_get_Econt(newlupri,newluerr,ls%setting,D,TwoElectronEcont,&
           & ndmat,modthresh)
      DO idmat = 1,ndmat
         energy(idmat) = 0.5E0_realk*TwoElectronEcont(idmat)*fac
      ENDDO
      IF(ls%setting%do_dft) THEN
         nbast = D(1)%nrow
         call II_get_xc_energy(newlupri,newluerr,ls%setting,nbast,D,Edfty,ndmat)
         DO idmat = 1,ndmat
            energy(idmat) = energy(idmat) + Edfty(idmat)
         enddo
      ENDIF
      DO idmat = 1,ndmat
         hD = mat_dotproduct(D(idmat),H1)
         energy(idmat) = energy(idmat) + hD*fac + ls%input%potnuc
      enddo
   END SUBROUTINE di_SCF_EnergyCont

!##########################################################
      subroutine di_GET_GbDsSingle(lupri,luerr,Dens,GbDs,setting)
        !*********************************************************
        ! Determine the G matrix for the 2-e contribution to sigma
        ! vector in RSP
        ! G([b,D]s) = 2-e part of Fock Matrix with a modified
        !             density [b,D]s (here called Dens)
        ! Thomas, sep 2010 
        !*********************************************************
        implicit none
        integer, intent(in) :: lupri,luerr
        type(Matrix), intent(in) :: Dens
        type(Matrix), intent(inout) :: GbDs  
        type(lssetting),optional :: setting !intent(inout)
        integer :: ndmat 
!       
        logical :: Dsym
        Dsym = .FALSE. !NONsymmetric Density matrix
        ndmat = 1
        IF(present(setting))THEN
           !This should be changed to a test like the MATSYM function
           ! for full matrices
           call II_get_Fock_mat(lupri,luerr,setting,Dens,Dsym,&
                                 & GbDs,ndmat,.FALSE.)
        ELSE
           call II_get_Fock_mat(lupri,luerr,lsint_fock_data%ls%setting,Dens,&
                                & Dsym,GbDs,ndmat,.FALSE.)
        ENDIF
      end subroutine di_GET_GbDsSingle

    !> \brief Determine the G matrix for the 2-e contribution to 
    !>        sigma vector in RSP,
    !>        and includes the XC contrib. to linear response
    !> \author Patrick Merlot
    !> \date 2013-03-26
    !> \param GbDs  The 2-e contribution to sigma vector in RSP
    !> \param Gxc   The xc cont to the linear response
    !> \param lupri Default print unit
    !> \param luerr Default error print unit
    !> \param setting Integral evalualtion settings
    !> \param Bmat  The b matrix G(b)
    !> \param nBmat The number of Bmat matrices
    !> \param nbast The number of basisfunctions
    !> \param Dmat  The Density matrix
     SUBROUTINE di_GET_GbDs_and_XC_linrsp_Array(GbDs,Gxc,lupri,luerr,&
                        & Bmat,nBmat,nbast,Dmat,do_dft,setting)
        IMPLICIT NONE
        INTEGER, intent(in)           :: lupri,luerr,nBmat,nbast
        TYPE(LSsetting),intent(inout), optional :: setting
        TYPE(Matrix), intent(in)      :: Bmat(nBmat)
        TYPE(Matrix), intent(in)      :: Dmat
        TYPE(Matrix), intent(inout)   :: GbDs(nBmat)
        TYPE(Matrix), intent(inout)   :: Gxc(nBmat)
        LOGICAL, intent(in)           :: do_dft
        !
        INTEGER                       :: iBmat
        LOGICAL                       :: ADMMexchange  , ADMMGCBASIS    
        ! -- ADMM modifications
        !     replace GdBs = J(B) + K(B)
        !    by      GdBs = J(B) + K(b) + X(B) - X(b) if ADMM
        IF(present(setting))THEN
            ADMMexchange = setting%scheme%ADMM_EXCHANGE
            ADMMGCBASIS  = setting%scheme%ADMM_GCBASIS
        ELSE
            ADMMexchange = lsint_fock_data%ls%setting%scheme%ADMM_EXCHANGE 
            ADMMGCBASIS  = lsint_fock_data%ls%setting%scheme%ADMM_GCBASIS
        ENDIF
        IF (ADMMGCBASIS) THEN
            ADMMexchange = .FALSE.
            ENDIF
        IF (ADMMexchange) THEN 
            ! GdBs = J(B) + K(b) + X(B) - X(b)
            call di_GET_GbDsArray_ADMM(lupri,luerr,Bmat,GbDs,nBmat,setting)
        ELSE 
            ! GdBs = J(B) + K(B)
            call di_GET_GbDsArray(lupri,luerr,Bmat,GbDs,nBmat,setting)
        ENDIF
        
        IF (do_dft) THEN  
            ! Add extra XC contributions to G 
            IF(present(setting))THEN
                call II_get_xc_linrsp(lupri,luerr,&
                      & setting,nbast,Bmat,Dmat,Gxc,nBmat) 
            ELSE
                call II_get_xc_linrsp(lupri,luerr,&
                      & lsint_fock_data%ls%setting,nbast,Bmat,Dmat,Gxc,nBmat) 
            ENDIF
            DO iBmat=1,nBmat
                call mat_daxpy(1E0_realk,Gxc(iBmat),GbDs(iBmat)) 
            ENDDO
        ENDIF
    END SUBROUTINE di_GET_GbDs_and_XC_linrsp_Array
    
    
    SUBROUTINE di_GET_GbDs_and_XC_linrsp_Single(GbDs,Gxc,lupri,luerr,&
                                    & Bmat,nbast,Dmat,do_dft, setting)
        IMPLICIT NONE
        INTEGER, intent(in)           :: lupri,luerr,nbast
        TYPE(LSsetting),intent(inout), optional :: setting
        TYPE(Matrix), intent(in)      :: Bmat
        TYPE(Matrix), intent(in)      :: Dmat
        TYPE(Matrix), intent(inout)   :: GbDs
        TYPE(Matrix), intent(inout)   :: Gxc
        logical, intent(in)           :: do_dft
        !
        type(matrix) :: GbDSArray(1), GxcArray(1), BmatArray(1)
        !
        call mat_init(GbDsArray(1),nbast,nbast)
        call mat_assign(GbDSArray(1),GbDs)
        call mat_init(BmatArray(1),nbast,nbast)
        call mat_assign(BmatArray(1),Bmat)
        IF (do_dft) THEN 
            call mat_init(GxcArray(1),nbast,nbast)
            call mat_assign(GxcArray(1),Gxc)
        ENDIF
        call di_GET_GbDs_and_XC_linrsp_Array(GbDsArray,GxcArray,lupri,&
                        &luerr,BmatArray,1,nbast,Dmat,do_dft,setting)
        call mat_assign(GbDs,GbDSArray(1))
        call mat_free(GbDSArray(1))
        call mat_free(BmatArray(1))
        IF (do_dft) THEN 
            call mat_assign(Gxc,GxcArray(1))
            call mat_free(GxcArray(1))
        ENDIF
    END SUBROUTINE di_GET_GbDs_and_XC_linrsp_Single
    
    
      subroutine di_GET_GbDsArray(lupri,luerr,Dens,GbDs,nDmat,setting)
        !*********************************************************
        ! Determine the G matrix for the 2-e contribution to sigma
        ! vector in RSP
        ! G([b,D]s) = 2-e part of Fock Matrix with a modified
        !             density [b,D]s (here called Dens)
        ! Sonia, October 2004
        ! Thomas, Feb 2010 (fixed unrestricted + added lsdalton lsint)
        !*********************************************************
        implicit none
        integer, intent(in) :: lupri,luerr,ndmat
        type(Matrix), intent(in) :: Dens(nDmat)
        type(Matrix), intent(inout) :: GbDs(nDmat)  !output
        type(lssetting),optional :: setting !intent(inout)
        !
        logical :: Dsym

        !This should be changed to a test like the MATSYM
        ! function for full matrices   
        Dsym = .FALSE. !NONsymmetric Density matrix
        IF(present(setting))THEN
           call II_get_Fock_mat(lupri,luerr,&
              & setting,Dens,Dsym,GbDs,ndmat,.FALSE.)
        ELSE
           call II_get_Fock_mat(lupri,luerr,&
              & lsint_fock_data%ls%setting,Dens,Dsym,GbDs,ndmat,.FALSE.)
        ENDIF
!        write (lupri,*) "FOCK mat in noADMM di_GET_GbDsArray()"
!        call mat_print(GbDs(1),1,GbDs(1)%nrow,1,GbDs(1)%ncol,lupri)

      end subroutine di_GET_GbDsArray


    !*********************************************************
    ! Determine the G matrix for the 2-e contribution to sigma
    ! vector in RSP using the ADMM exchange approximation
    ! GdBs = J(B) + K(b) + X(B) - X(b)
    ! G([b,D]s) = 2-e part of Fock Matrix with a modified
    !             density [b,D]s (here called Dens)
    ! Patrick, April 2013
    !*********************************************************
    SUBROUTINE di_GET_GbDsArray_ADMM(lupri,luerr,Dens,GbDs,nDmat,setting)
        implicit none
        integer, intent(in)                      :: lupri,luerr,ndmat
        type(Matrix), intent(in)                 :: Dens(nDmat)
        type(Matrix), intent(inout)              :: GbDs(nDmat)  !output
        type(lssetting), intent(inout), optional :: setting
        !
        logical                :: Dsym, copy_IntegralTransformGC
        type(Matrix)           :: Dens_AO(nDmat), K(nDmat)
        integer                :: i
        !
        IF(present(setting))THEN
            call di_GET_GbDsArray_ADMM_setting(lupri,luerr,Dens,GbDs,nDmat,&
                      & setting)
        ELSE
            call di_GET_GbDsArray_ADMM_setting(lupri,luerr,Dens,GbDs,nDmat,&
                      & lsint_fock_data%ls%setting)
        ENDIF
        !
        CONTAINS
          SUBROUTINE di_GET_GbDsArray_ADMM_setting(lupri,luerr,Dens,&
                                                    & GbDs,nDmat,setting)
            implicit none
            integer, intent(in)            :: lupri,luerr,ndmat
            type(Matrix), intent(in)       :: Dens(nDmat)
            type(Matrix), intent(inout)    :: GbDs(nDmat)
            type(lssetting), intent(inout) :: setting
            !
            logical                :: Dsym, copy_IntegralTransformGC
            logical                :: inc_scheme, do_inc
            type(Matrix)           :: Dens_AO(nDmat), K(nDmat)
            integer                :: i
            !
            IF(matrix_type .EQ. mtype_unres_dense) THEN
                call lsquit('di_GET_GbDsArray_ADMM does not support &
                             &unrestricted matrices',lupri)
            ENDIF
            call get_incremental_settings(inc_scheme,do_inc)
            IF(inc_scheme .OR. do_inc)THEN
                call lsquit('II_get_Fock_mat incremental scheme not &
                           & allowed in di_GET_GbDsArray_ADMM_setting()',lupri)
            ENDIF
            !This should be changed to a test like the MATSYM
            ! function for full matrices   
            Dsym = .FALSE. !NONsymmetric Density matrix

            ! Get rid of Grand canonical
            copy_IntegralTransformGC = setting%IntegralTransformGC
            IF(copy_IntegralTransformGC) THEN
                setting%IntegralTransformGC = .FALSE.
                ! we do this manually in order to get incremental correct
                ! change D to AO basis (currently in GCAO basis)
                DO i=1,nDmat
                    call mat_init(Dens_AO(i),Dens(1)%nrow,Dens(1)%ncol)
                    call GCAO2AO_transform_matrixD2(Dens(i),Dens_AO(i),&
                                                    & setting,lupri)
                ENDDO
            ELSE ! no GC transformation
                DO i=1,nDmat
                  call mat_init(Dens_AO(i),Dens(1)%nrow,Dens(1)%ncol)
                  call mat_assign(Dens_AO(i),Dens(i))
                ENDDO            
            ENDIF

            ! J(B)
            DO i=1,nDmat
                call mat_zero(GbDs(i))
            ENDDO
            call II_get_coulomb_mat(LUPRI,LUERR,setting,Dens_AO,GbDs,nDmat)   
            !write (lupri,*) "Coulomb mat in di_GET_GbDsArray_ADMM_setting()"
            !call mat_print(GbDs(1),1,GbDs(1)%nrow,1,GbDs(1)%ncol,lupri)

            ! K(B)
            DO i=1,nDmat
                call mat_init(K(i),GbDs(1)%nrow,GbDs(1)%ncol)
                call mat_zero(K(i))
            ENDDO
            IF(setting%scheme%daLinK) THEN
                call lsquit('di_GET_GbDsArray_ADMM not supported in daLink',lupri)
            ENDIF
            call II_get_exchange_mat(LUPRI,LUERR,setting,Dens_AO,nDmat,Dsym,K)
            !write (lupri,*) "Exchange mat in ADMM di_GET_GbDsArray_ADMM_setting()"
            !call mat_print(K(1),1,K(1)%nrow,1,K(1)%ncol,lupri)
            DO i=1,ndmat
                call mat_daxpy(1E0_realk,K(i),GbDs(i))
                call mat_free(K(i))
            ENDDO

            ! Transform back to GCAO basis 
            IF(copy_IntegralTransformGC) THEN
                ! Reset to the original value of IntegralTransformGC
                setting%IntegralTransformGC = copy_IntegralTransformGC 
                DO i=1,nDmat
                    call AO2GCAO_transform_matrixF(GbDs(i),setting,lupri)
                ENDDO
            ENDIF
            
            
            ! DEBUG PAT START -----------------------------    
            write (lupri,*) "J+K mat in ADMM di_GET_GbDsArray_ADMM_setting()"
            call mat_print(GbDs(1),1,GbDs(1)%nrow,1,GbDs(1)%ncol,lupri)

            ! FREE MEMORY
            DO I=1,ndmat
            call mat_free(Dens_AO(I))
            ENDDO
          END SUBROUTINE di_GET_GbDsArray_ADMM_setting
      !CONTAINS END
      END SUBROUTINE di_GET_GbDsArray_ADMM
      
      
      subroutine di_lowdin(S,S_sqrt,S_minus_sqrt,n,lupri)
        implicit none
        real(realk), target :: S(n*n), S_sqrt(n*n), S_minus_sqrt(n*n)
        type(Matrix)        :: tS, tS_sqrt, tS_minus_sqrt
        integer             :: n, lupri
        
        call lowdin_diag(n,S,S_sqrt, S_minus_sqrt, lupri)
#if 0 && defined(HAVE_BSM)
        select case (matrix_type)
        case(mtype_sparse_block)
           
           call mat_init(tS,n,n)
           call mat_init(tS_sqrt,n,n)
           call mat_init(tS_minus_sqrt,n,n)
           
           call mat_set_from_full(S,1.0E0_realk,tS)
           
           call lowdin_schulz(tS,tS_sqrt,tS_minus_sqrt,lupri)
           
           call mat_to_full(tS_sqrt,1.0E0_realk,S_sqrt)
           call mat_to_full(tS_minus_sqrt,1.0E0_realk,S_minus_sqrt)
           
           call mat_free(tS)
           call mat_free(tS_sqrt)
           call mat_free(tS_minus_sqrt)
        case default
           call lsquit('Only dense and BSM matrices are supported by di_lowdin()!',-1)
        end select
#endif
      end subroutine di_lowdin

!> \brief Get the three matrices of dipole integrals
!> \author C. Nygaard
!> \date Nov. 19 2012
!> \param V The array of dipole matrices (3 of them)
      subroutine di_get_dipole (V)

      implicit none

      type(matrix), intent(inout) :: V(:)

      if (size(V)/=3) call lsquit('Wrong dimension of dipole matrix &
                                 &array in di_get_dipole',lsint_fock_data%lupri)

      call II_get_prop (lsint_fock_data%lupri,lsint_fock_data%luerr,&
                       &lsint_fock_data%ls%setting,V(1:3),3,'DIPLEN ')

      end subroutine di_get_dipole

      subroutine di_get_GAOL_lsdalton(Lmat,Gmat)
        !*********************************************************
        ! Determine the the 2-e contribution (Coulomb and exchange) 
        ! contribution to the G(L) matrix. (ask Cecilie Nygaard)   
        ! G(L) = (g_{abcd}- 1/2 * g_{adcb}) Lmat_{cd}
        ! Thomas, July 2011
        !*********************************************************
        implicit none
        !Input/Output
        type(Matrix), intent(in) :: Lmat
        type(Matrix), intent(inout) :: Gmat  !output
        !local
        type(Matrix)             :: temp2
        integer :: nbast,LUINT,A,B,C,D
        logical :: Dsym,file_exists
        real(realk) :: JFAC,KFAC
        real(realk),pointer :: g(:,:,:,:),LmatFull(:,:),GmatFull(:,:)
        real(realk),pointer :: LmatFull2(:,:),GmatFull2(:,:)
        nbast = Lmat%nrow
        IF(.FALSE.)THEN
           call mat_init(TEMP2,Gmat%nrow,Gmat%ncol)
           Dsym = .false. !NONsymmetric Density matrix
           call mat_zero(Gmat)
           call mat_zero(temp2)
           call II_get_coulomb_mat(lsint_fock_data%lupri,lsint_fock_data%luerr,lsint_fock_data%ls%setting,Lmat,Gmat,1)
           call II_get_exchange_mat(lsint_fock_data%lupri,lsint_fock_data%luerr, &
                & lsint_fock_data%ls%setting,Lmat,1,Dsym,TEMP2)           
           call mat_daxpy (1.0E0_realk, temp2, Gmat)
           call mat_free(temp2)
           if (matrix_type /= mtype_unres_dense) then
              call mat_scal (0.5E0_realk, Gmat)
           endif
        ELSE
           call mem_alloc(g,nbast,nbast,nbast,nbast)
           INQUIRE(file='ERI.Integrals',EXIST=file_exists)
           IF(file_exists)THEN
              LUINT = -1
              call LSOPEN(LUINT,'ERI.Integrals','OLD','UNFORMATTED')
              do D=1,nbast
                 do C=1,nbast
                    do B=1,nbast
                       READ(LUINT)(g(A,B,C,D),A=1,nbast)
                    ENDDO
                 ENDDO
              ENDDO
           ELSE
              call II_get_4center_eri(lsint_fock_data%LUPRI,lsint_fock_data%LUERR,&
                   & lsint_fock_data%ls%SETTING,g,nbast,nbast,nbast,nbast)
              LUINT = -1
              call LSOPEN(LUINT,'ERI.Integrals','UNKNOWN','UNFORMATTED')
              do D=1,nbast
                 do C=1,nbast
                    do B=1,nbast
                       WRITE(LUINT)(g(A,B,C,D),A=1,nbast)
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
           call LSCLOSE (LUINT,"KEEP")
           JFAC = 1.0E0_realk
           KFAC = - lsint_fock_data%ls%setting%scheme%exchangeFactor * 0.5E0_realk
           call mat_init(TEMP2,Gmat%nrow,Gmat%ncol)
           IF(matrix_type.EQ.mtype_dense)THEN
              call mat_zero(TEMP2)
              call GAOLcont(TEMP2%elms,Lmat%elms,nbast,g,JFAC,KFAC)
           ELSEIF(matrix_type.EQ.mtype_unres_dense)THEN
              KFAC = KFAC * 2.E0_realk
              call mem_alloc(Lmatfull,nbast,nbast)
              call mem_alloc(Gmatfull,nbast,nbast)
              call mem_alloc(Lmatfull2,nbast,nbast)
              call mem_alloc(Gmatfull2,nbast,nbast)
              CALL DCOPY(nbast*nbast,Lmat%elms, 1,Lmatfull,1)
              CALL DCOPY(nbast*nbast,Lmat%elmsb,1,Lmatfull2,1)
              call ls_dzero(Gmatfull,nbast*nbast)
              call ls_dzero(Gmatfull2,nbast*nbast)
              call GAOLcont2(Gmatfull,Lmatfull,Gmatfull2,Lmatfull2,nbast,g,JFAC,KFAC)
              CALL DCOPY(nbast*nbast,Gmatfull, 1,TEMP2%elms,1)
              CALL DCOPY(nbast*nbast,Gmatfull2,1,TEMP2%elmsb,1)
              call mem_dealloc(Lmatfull)
              call mem_dealloc(Gmatfull)
              call mem_dealloc(Lmatfull2)
              call mem_dealloc(Gmatfull2)
           ELSE
              call mem_alloc(Lmatfull,nbast,nbast)
              call mem_alloc(Gmatfull,nbast,nbast)
              call mat_to_full(Lmat, 1.d0, LmatFull)
              call ls_dzero(Gmatfull,nbast*nbast)
              call GAOLcont(Gmatfull,Lmatfull,nbast,g,JFAC,KFAC)
              call mem_dealloc(Lmatfull)
              call mat_set_from_full(Gmatfull,1.d0, TEMP2)
              call mem_dealloc(Gmatfull)
           ENDIF
           call mem_dealloc(g)

           call mat_assign(Gmat,TEMP2)
           if (matrix_type /= mtype_unres_dense) then
              !call mat_scal (0.5E0_realk, Gmat)
           endif
           call mat_free(TEMP2)
        ENDIF        

    end subroutine di_get_GAOL_lsdalton

    subroutine GAOLcont(Gmat,Lmat,nbast,g,JFAC,KFAC)
      implicit none
      integer :: nbast
      real(realk) :: Gmat(nbast,nbast),Lmat(nbast,nbast)
      real(realk) :: g(nbast,nbast,nbast,nbast),JFAC,KFAC
      !
      integer ::A,B,C,D
      real(realk) :: Lmatdc

      do D=1,nbast
         do C=1,nbast
            Lmatdc=Lmat(d,c)
            do B=1,nbast
               do A=1,nbast
                  Gmat(a,b)=Gmat(a,b)+Lmatdc*& 
                       & (JFAC*g(a,b,c,d) + KFAC*g(a,d,c,b) ) 
               enddo
            enddo
         enddo
      enddo
    end subroutine GAOLcont

    subroutine GAOLcont2(Gmat,Lmat,Gmat2,Lmat2,nbast,g,JFAC,KFAC)
      implicit none
      integer :: nbast
      real(realk) :: Gmat(nbast,nbast),Lmat(nbast,nbast)
      real(realk) :: Gmat2(nbast,nbast),Lmat2(nbast,nbast)
      real(realk) :: g(nbast,nbast,nbast,nbast),JFAC,KFAC
      !
      integer ::A,B,C,D
      real(realk) :: Lmatdc,Lmatdc1,Lmatdc2
      do D=1,nbast
         do C=1,nbast
            Lmatdc1=Lmat(d,c)
            Lmatdc2=Lmat2(d,c)
            Lmatdc=Lmat(d,c)+Lmat2(d,c)
            do B=1,nbast
               do A=1,nbast
                  Gmat(a,b)=Gmat(a,b)+Lmatdc*JFAC*g(a,b,c,d) &
                  & + Lmatdc1*KFAC*g(a,d,c,b) 
                  Gmat2(a,b)=Gmat2(a,b)+Lmatdc*JFAC*g(a,b,c,d) &
                  & + Lmatdc2*KFAC*g(a,d,c,b) 
               enddo
            enddo
         enddo
      enddo
    end subroutine GAOLcont2

      subroutine di_get_sigma_xc_contSingle_lsd(b,D,sigma_xc_cont)   
        implicit none
        type(Matrix), intent(in) :: b,D
        type(Matrix), intent(inout) :: sigma_xc_cont  !output
        integer :: nbast,nbmat
        TYPE(MATRIX)  :: Sarray(1)
        nbast = D%nrow
        nbmat = 1
        !NOT OPTIMAL USE OF MEMORY OR ANYTHING
        call mat_init(Sarray(1),sigma_xc_cont%nrow,sigma_xc_cont%ncol)
        call II_get_xc_linrsp(lsint_fock_data%lupri,lsint_fock_data%luerr,&
             &lsint_fock_data%ls%setting,nbast,(/B/),D,Sarray,nbmat)
        call mat_assign(sigma_xc_cont,Sarray(1))
        call mat_free(Sarray(1))
      end subroutine di_get_sigma_xc_contSingle_lsd

    subroutine di_get_sigma_xc_contArray_lsd(b,D,sigma_xc_cont,nbmat)   
      implicit none
      integer, intent(in) :: nbmat
      type(Matrix), intent(in) :: b(nbmat),D
      type(Matrix), intent(inout) :: sigma_xc_cont(nbmat)  !output
      integer :: nbast
      nbast = D%nrow
      call II_get_xc_linrsp(lsint_fock_data%lupri,lsint_fock_data%luerr,&
           &lsint_fock_data%ls%setting,nbast,b,D,sigma_xc_cont,nbmat)
    end subroutine di_get_sigma_xc_contArray_lsd

    SUBROUTINE di_decpacked(lupri,luerr,ls,nbast,D)
      IMPLICIT NONE
      TYPE(Matrix),intent(in) :: D
      integer,intent(in)      :: lupri,luerr,nbast
      type(lsitem),intent(inout) :: ls
      !
      TYPE(Matrix)  :: K,Kdec,tempm3
      real(realk),pointer   :: integrals(:,:,:,:)
      real(realk) :: ExchangeFactor,Kcont
      real(realk),pointer :: Dfull(:,:),KdecLocFull(:,:)
      integer :: nbatches,iorb,JK,ao_iX,ao_iY,lu_pri, lu_err,thread_idx,nthreads,idx,nbatchesXY
      integer :: X,Y,dimX,dimY,batch_iX,batch_iY,i,j,MinAObatch,MaxAllowedDim,MaxActualDim,ib,id,ix,iy
      logical :: doscreen,fullrhs
      TYPE(DECscreenITEM)    :: DecScreen
      integer, pointer :: orb2batch(:), batchdim(:),batchsize(:), batchindex(:)
      type(batchtoorb), pointer :: batch2orb(:)
      character :: INTSPEC(5)
      doscreen = ls%setting%SCHEME%CS_SCREEN.OR.ls%setting%SCHEME%PS_SCREEN
      nullify(orb2batch)
      nullify(batchdim)
      nullify(batch2orb)
      nullify(batchsize)
      nullify(batchindex)
      
      call determine_maxBatchOrbitalsize(lupri,ls%setting,MinAObatch)
      print*,'MinAObatch',MinAObatch
      MaxAllowedDim = MinAObatch
      call mem_alloc(orb2batch,nbast)
      call build_batchesofAOS(lupri,ls%setting,MaxAllowedDim,&
           & nbast,MaxActualDim,batchsize,batchdim,batchindex,nbatchesXY,orb2Batch)
      print*,'nbatchesXY',nbatchesXY
      call mem_alloc(batch2orb,nbatchesXY)
      do idx=1,nbatchesXY
         call mem_alloc(batch2orb(idx)%orbindex,batchdim(idx) )
         batch2orb(idx)%orbindex = 0
         batch2orb(idx)%norbindex = 0
      end do
      do iorb=1,nbast
         idx = orb2batch(iorb)
         batch2orb(idx)%norbindex = batch2orb(idx)%norbindex+1
         j = batch2orb(idx)%norbindex
         batch2orb(idx)%orbindex(j) = iorb
      end do
      INTSPEC(1) = 'R' 
      INTSPEC(2) = 'R'
      INTSPEC(3) = 'R'
      INTSPEC(4) = 'R'
      INTSPEC(5) = 'C' !operator
      call II_precalc_DECScreenMat(DecScreen,lupri,luerr,ls%setting,nbatches,nbatchesXY,nbatchesXY,INTSPEC)
      
      IF(doscreen)then
         call II_getBatchOrbitalScreenK(DecScreen,ls%setting,&
              & nbast,nbatchesXY,nbatchesXY,batchsize,batchsize,batchindex,batchindex,&
              & batchdim,batchdim,lupri,luerr)
      endif
      FullRHS = .FALSE.

      !############################################################################
      !           Exchange Part
      !############################################################################
      call mem_alloc(Dfull,D%nrow,D%ncol)
      call mat_to_full(D,1E0_realk,DFULL)
      ExchangeFactor = -1E0_realk
      call mat_init(K,nbast,nbast)
      call mat_zero(K)
      CALL II_get_exchange_mat(lupri,luerr,ls%setting,D,1,.TRUE.,K)

      call mem_alloc(KdecLocFULL,nbast,nbast)
      KdecLocFULL = 0E0_realk

      BatchY: do Y = 1,nbatchesXY
         dimY = batchdim(Y)

         BatchX: do X = 1,nbatchesXY
            dimX = batchdim(X)

            nullify(integrals)
            allocate(integrals(dimX,nbast,dimY,nbast))
            integrals = 0E0_realk
            IF(doscreen)ls%setting%LST_GAB_LHS => DECSCREEN%batchGabKLHS(X)%p
            IF(doscreen)ls%setting%LST_GAB_RHS => DECSCREEN%batchGabKRHS(Y)%p

            call II_GET_DECPACKED4CENTER_K_ERI(LUPRI,LUERR,ls%SETTING,&
                 & integrals,batchindex(X),batchindex(Y),batchsize(X),batchsize(Y),&
                 & dimX,nbast,dimY,nbast,nbatches,INTSPEC,fullRHS)

            do batch_iY = 1,dimY
               iY = batch2orb(Y)%orbindex(batch_iY)
               do batch_iX = 1,dimX
                  iX = batch2orb(X)%orbindex(batch_iX)
                  Kcont = 0E0_realk
                  DO iD=1,nbast
                     DO iB=1,nbast
                        Kcont = Kcont + integrals(batch_iX,iB,batch_iY,iD)*Dfull(iB,iD)
                     ENDDO
                  ENDDO
                  KdecLocFULL(iX,iY)=KdecLocFULL(iX,IY)+Kcont
               ENDDO
            ENDDO
            deallocate(integrals)
            nullify(integrals)
         ENDDO BatchX
      ENDDO BatchY
      call mat_init(Kdec,nbast,nbast)
      call mat_set_from_full(KdecLocFULL,exchangeFactor,Kdec)
      call mem_dealloc(KdecLocFULL)

      nullify(ls%setting%LST_GAB_LHS)
      nullify(ls%setting%LST_GAB_RHS)
      call free_decscreen(DECSCREEN)
      call mem_dealloc(orb2batch)
      call mem_dealloc(batchdim)
      call mem_dealloc(batchsize)
      call mem_dealloc(batchindex)
      do X=1,nbatchesXY
         call mem_dealloc(batch2orb(X)%orbindex)
         nullify(batch2orb(X)%orbindex)
      end do
      call mem_dealloc(batch2orb)
      
      call mat_init(tempm3,nbast,nbast)
      call mat_add(1E0_realk,Kdec,-1E0_realk,K,tempm3)
      write(lupri,*) 'QQQ DI_DEBUG_DECPACK K STD    ',mat_trab(K,K)
      write(lupri,*) 'QQQ DI_DEBUG_DECPACK K DECPACK',mat_trab(Kdec,Kdec)
      write(lupri,*) 'QQQ DIFF',ABS(mat_trab(tempm3,tempm3))
      IF(ABS(mat_trab(tempm3,tempm3)).LE. 1E-15_realk)THEN
         write(lupri,*)'QQQ SUCCESFUL DECPACK K TEST'
      ELSE
         WRITE(lupri,*)'the Kref'
         call mat_print(K,1,nbast,1,nbast,lupri)
         WRITE(lupri,*)'the Kdec'
         call mat_print(Kdec,1,nbast,1,nbast,lupri)
         WRITE(lupri,*)'the Diff'
         call mat_print(tempm3,1,nbast,1,nbast,lupri)
         CALL lsQUIT('DECPACKED K TEST FAILED',lupri)
      ENDIF
      call mat_free(K)
      call mat_free(Kdec)
      call mat_free(tempm3)

    END SUBROUTINE DI_DECPACKED

    SUBROUTINE di_debug_4center_eri_interest(lupri,lu_err,ls,nbast)
      IMPLICIT NONE
      integer,intent(in)      :: lupri,lu_err,nbast
      type(lsitem),intent(inout) :: ls
      !
      real(realk),pointer   :: integrals(:,:,:,:)
      real(realk),pointer   :: integrals2(:,:,:,:)
      integer :: iB,iA,iC,iD
      call mem_alloc(integrals2,nbast,nbast,nbast,nbast)
      call II_get_4center_eri(LUPRI,LU_ERR,ls%setting,integrals2,nbast,nbast,nbast,nbast)

      ls%setting%scheme%interest = .FALSE.
      call mem_alloc(integrals,nbast,nbast,nbast,nbast)
      call II_get_4center_eri(LUPRI,LU_ERR,ls%setting,integrals,nbast,nbast,nbast,nbast)
      do iD=1,nbast
         do iC=1,nbast
            do iB=1,nbast
               do iA=1,nbast
                  IF (ABS(integrals(iA,iB,iC,iD) - integrals2(iA,iB,iC,iD)).GT. 1E-9_realk) THEN
                     WRITE(lupri,'(4I3,A,F16.12)')iA,iB,iC,iD,'Thermite',integrals(ia,ib,ic,id)
                     WRITE(lupri,'(4I3,A,F16.12)')iA,iB,iC,iD,'Interest',integrals2(ia,ib,ic,id)

                  ENDIF
               enddo
            enddo
         enddo
      enddo
      call mem_dealloc(integrals)
      call mem_dealloc(integrals2)
    END SUBROUTINE di_debug_4center_eri_interest

end MODULE dal_interface
