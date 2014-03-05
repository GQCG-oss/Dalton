MODULE IntegralInterfaceIchorMod
  use precision
  use TYPEDEFTYPE, only: LSSETTING, LSINTSCHEME, LSITEM, integralconfig,&
       & BASISSETLIBRARYITEM
  use basis_type, only: free_basissetinfo
  use basis_typetype,only: BASISSETINFO,BASISINFO
  use BuildBasisSet, only: Build_BASIS
  use Matrix_module, only: MATRIX, MATRIXP
  use Integralparameters
  use molecule_typetype, only: MOLECULEINFO, MOLECULE_PT, ATOMITEM
  use molecule_type, only: build_pointmolecule, DETERMINE_MAXCOOR, &
       & free_moleculeinfo, build_atomicmolecule, print_mol
  use integral_type, only: INTEGRALINPUT
  use TYPEDEF, only: getNbasis, retrieve_output, gcao2ao_transform_matrixd2, &
       & retrieve_screen_output, ao2gcao_transform_matrixf, &
       & gcao2ao_transform_fulld, ao2gcao_transform_fullf, &
       & ao2gcao_half_transform_matrix
  use molecule_module, only: DETERMINE_NBAST2
  use matrix_operations, only: mat_dotproduct, matrix_type, mtype_unres_dense,&
       & mat_daxpy, mat_init, mat_free, mat_write_to_disk, mat_print, mat_zero,&
       & mat_scal, mat_mul, mat_assign, mat_trans
  use matrix_util, only: mat_get_isym, util_get_symm_part,util_get_antisymm_part, matfull_get_isym
  use memory_handling, only: mem_alloc, mem_dealloc, debug_mem_stats
  use IntegralInterfaceMOD, only: II_get_4center_eri
  use IchorErimoduleHost!,only: MAIN_ICHORERI_DRIVER
  use lsmatrix_operations_dense
  use LSmatrix_type
  
  public::  II_unittest_Ichor
  private
CONTAINS
!> \brief Calculates overlap integral matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param S the overlap matrix
SUBROUTINE II_unittest_Ichor(LUPRI,LUERR,SETTING,DebugIchorOption)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,DebugIchorOption
!
#ifdef VAR_ICHOR
real(realk),pointer   :: integralsII(:,:,:,:),integralsIchor(:,:,:,:)
integer :: dim1,dim2,dim3,dim4,A,B,C,D,iprint,nbast(4),ibasiselm(4)
integer :: iBasis1,ibasis2,ibasis3,ibasis4,icharge,nbasis,nPass,ipass,itest
logical :: dirac,doprint,debug
TYPE(MOLECULEINFO),pointer :: Originalmolecule
TYPE(MOLECULEINFO),pointer :: atomicmolecule(:)
TYPE(BASISSETLIBRARYITEM) :: LIBRARY
CHARACTER(len=9)     :: BASISLABEL
TYPE(BASISINFO),pointer :: unittestBASIS(:)
TYPE(BASISINFO),pointer :: originalBASIS
CHARACTER(len=80)    :: BASISSETNAME
CHARACTER(len=20)    :: BASISTYPE(10)
real(realk)          :: Rxyz(3)
type(lsmatrix)       :: FINALVALUE(2)
logical      :: spherical,savedospherical,SpecialPass
logical      :: FAIL(10,10,10,10),ALLPASS,SameMOL
Character    :: intSpec(5)
integer :: iBasis1Q,iBasis2Q,iBasis3Q,iBasis4Q
integer :: nBasisA,nBasisB,nBasisC,nBasisD,iPassStart,iPassEnd
integer,pointer :: iBasisA(:),iBasisB(:),iBasisC(:),iBasisD(:)

intSpec(1) = 'R'
intSpec(2) = 'R'
intSpec(3) = 'R'
intSpec(4) = 'R'
intSpec(5) = 'C'
!WRITE(lupri,*)'Before IchorUnitTest'
!call debug_mem_stats(LUPRI)

!Special types
!   A          B          C          D          OVERALL
!Seg1Prim     Seg        Gen      Seg1Prim      SegP      (Tested as Option 1)
!Seg1Prim   Seg1Prim   Seg1Prim   Seg1Prim      Seg       (Tested as Option 2)
!Seg1Prim   Seg1Prim   Seg1Prim   Seg1Prim      Seg1Prim  (Tested as Option 3)  
!Seg1Prim     Gen        Seg      Seg1Prim      SegQ      (Tested as Option 4)
!  Gen        Gen        Gen        Gen         Gen       (Tested as Option 5)

SpecialPass = .FALSE.
SELECT CASE(DebugIchorOption)
CASE(0)
   !all types and all passes
   nbasisA = 9; nbasisB = 9; nbasisC = 9; nbasisD = 9
   IpassStart = 1; IpassEnd = 2
   call mem_alloc(iBasisA,nBasisA)
   call mem_alloc(iBasisB,nBasisB)
   call mem_alloc(iBasisC,nBasisC)
   call mem_alloc(iBasisD,nBasisD)
   do iBasis1Q=1,nbasisA
      iBasisA(iBasis1Q) = iBasis1Q
   enddo
   do iBasis2Q=1,nbasisB
      iBasisB(iBasis2Q) = iBasis2Q
   enddo
   do iBasis3Q=1,nbasisC
      iBasisC(iBasis3Q) = iBasis3Q
   enddo
   do iBasis4Q=1,nbasisD
      iBasisD(iBasis4Q) = iBasis4Q
   enddo
CASE(1)
   !Special types
   !   A          B          C          D          OVERALL
   !Seg1Prim     Seg        Gen      Seg1Prim      SegP
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   IpassStart = 1; IpassEnd = 2
   SpecialPass = .TRUE.
   call mem_alloc(iBasisA,nBasisA)
   call mem_alloc(iBasisB,nBasisB)
   call mem_alloc(iBasisC,nBasisC)
   call mem_alloc(iBasisD,nBasisD)
   iBasisA(1) = 1; iBasisA(2) = 2; iBasisA(3) = 3  !Seg1Prim
   iBasisB(1) = 4; iBasisB(2) = 5; iBasisB(3) = 6  !Seg
   iBasisC(1) = 7; iBasisC(2) = 8; iBasisC(3) = 9  !Gen
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3  !Seg1Prim
CASE(2)
   !Special types
   !   A          B          C           D          OVERALL
   !  Seg        Seg1Prim   Seg1Prim  Seg1Prim      Seg
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   IpassStart = 1; IpassEnd = 2
   call mem_alloc(iBasisA,nBasisA)
   call mem_alloc(iBasisB,nBasisB)
   call mem_alloc(iBasisC,nBasisC)
   call mem_alloc(iBasisD,nBasisD)
   !Pure Seg (S,P,D ; S,P,D | S,P,D ; S,P,D)
   iBasisA(1) = 4; iBasisA(2) = 5; iBasisA(3) = 6
   iBasisB(1) = 1; iBasisB(2) = 2; iBasisB(3) = 3
   iBasisC(1) = 1; iBasisC(2) = 2; iBasisC(3) = 3
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3
CASE(3)
   !Special types
   !   A          B          C          D          OVERALL
   !Seg1Prim   Seg1Prim   Seg1Prim   Seg1Prim      Seg1Prim 
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   IpassStart = 1; IpassEnd = 2
   call mem_alloc(iBasisA,nBasisA)
   call mem_alloc(iBasisB,nBasisB)
   call mem_alloc(iBasisC,nBasisC)
   call mem_alloc(iBasisD,nBasisD)
   !Pure Seg1Prim (S,P,D ; S,P,D | S,P,D ; S,P,D)
   iBasisA(1) = 1; iBasisA(2) = 2; iBasisA(3) = 3
   iBasisB(1) = 1; iBasisB(2) = 2; iBasisB(3) = 3
   iBasisC(1) = 1; iBasisC(2) = 2; iBasisC(3) = 3
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3
CASE(4)
   !Special types
   !   A          B          C          D          OVERALL
   !Seg1Prim     Gen        Seg      Seg1Prim      SegQ
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   IpassStart = 1; IpassEnd = 2
   SpecialPass = .TRUE.
   call mem_alloc(iBasisA,nBasisA)
   call mem_alloc(iBasisB,nBasisB)
   call mem_alloc(iBasisC,nBasisC)
   call mem_alloc(iBasisD,nBasisD)
   iBasisA(1) = 1; iBasisA(2) = 2; iBasisA(3) = 3  !Seg1Prim
   iBasisB(1) = 7; iBasisB(2) = 8; iBasisB(3) = 9  !Gen
   iBasisC(1) = 4; iBasisC(2) = 5; iBasisC(3) = 6  !Seg
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3  !Seg1Prim
CASE(5)
   !Gen1
   !Special types
   !   A          B          C          D          OVERALL
   !  Gen        Gen        Gen        Gen         Gen    
   !Pure Gen (S,P,D ; S,P,D | S,P,D ; S,P,D)
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   IpassStart = 1; IpassEnd = 2
   SpecialPass = .TRUE. !otherwise too expensive
   call mem_alloc(iBasisA,nBasisA)
   call mem_alloc(iBasisB,nBasisB)
   call mem_alloc(iBasisC,nBasisC)
   call mem_alloc(iBasisD,nBasisD)
   iBasisA(1) = 7; iBasisA(2) = 8; iBasisA(3) = 9
   iBasisB(1) = 7; iBasisB(2) = 8; iBasisB(3) = 9
   iBasisC(1) = 7; iBasisC(2) = 8; iBasisC(3) = 9
   iBasisD(1) = 7; iBasisD(2) = 8; iBasisD(3) = 9
CASE DEFAULT
   CALL LSQUIT('unknown option in Debug Ichor.',-1)
END SELECT

WRITE(lupri,*)'II_test_Ichor'
ALLPASS = .TRUE.
do A=1,80
   BASISSETNAME(A:A) = ' '
enddo
!
BASISTYPE(1) = 'UnitTest_segS1p     '
BASISTYPE(2) = 'UnitTest_segP1p     '
BASISTYPE(3) = 'UnitTest_segD1p     '
BASISTYPE(4) = 'UnitTest_segS       '
BASISTYPE(5) = 'UnitTest_segP       '
BASISTYPE(6) = 'UnitTest_segD       '
BASISTYPE(7) = 'UnitTest_genS       '
BASISTYPE(8) = 'UnitTest_genP       '
BASISTYPE(9) = 'UnitTest_genD       '
BASISTYPE(10) = 'UnitTest_segSP      '
!issues with SP in that 
!       IF(dim2.GT.dim1)CYCLE
!       IF(dim3.GT.dim1)CYCLE
!       IF(dim4.GT.dim3)CYCLE
!       IF(dim3+dim4.GT.dim1+dim2)CYCLE
! general do not work either
! do not work


Originalmolecule => SETTING%MOLECULE(1)%p
originalBasis => SETTING%BASIS(1)%p

dirac = .FALSE.
iprint=0
setting%scheme%intprint = 0
doprint = .FALSE.!.TRUE.
itest = 1

do Ipass = IpassStart,IpassEnd
   WRITE(lupri,*)'Number of Passes',Ipass
   !=========================================================================================================
   !                    Build Molecule
   !=========================================================================================================
   allocate(atomicmolecule(4))
   
   !   ICHARGE=6; Rxyz(1)=0.813381591740787E0_realk; Rxyz(2)=1.059191498062862E-002_realk;Rxyz(3)=0.689158554601339E0_realk; 
   !   call build_unittest_atomicmolecule(atomicmolecule(1),ICHARGE,Rxyz,1,lupri)
   !   ICHARGE=8; Rxyz(1)=0.762266389544351E0_realk; Rxyz(2)=2.983877565461657E-003_realk;Rxyz(3)=0.524979148086261E0_realk; 
   !   call build_unittest_atomicmolecule(atomicmolecule(2),ICHARGE,Rxyz,1,lupri)
   !   ICHARGE=9; Rxyz(1)=0.736938390171405E0_realk; Rxyz(2)=0.108186821166992E0_realk;Rxyz(3)=0.113699152299640E0_realk; 
   !   nPass = 2
   !   call build_unittest_atomicmolecule(atomicmolecule(3),ICHARGE,Rxyz,nPass,lupri)
   !   ICHARGE=17; Rxyz(1)=0.574178167982901E0_realk; Rxyz(2)=6.886728949849219E-2_realk;Rxyz(3)=9.213548738546455E-2_realk; 
   !   call build_unittest_atomicmolecule(atomicmolecule(4),ICHARGE,Rxyz,1,lupri)
   !nPassesP = 20 (OpenMP improved)
   nPass = 1
   IF(iPass.EQ.2)nPass = 5
   IF(SpecialPass)nPass = 3
   ICHARGE=6; Rxyz(1)=0.813381591740787E0_realk; Rxyz(2)=1.059191498062862E0_realk;Rxyz(3)=0.889158554601339E0_realk; 
   call build_unittest_atomicmolecule(atomicmolecule(1),ICHARGE,Rxyz,nPass,lupri)


   IF(iPass.EQ.2)nPass = 4
   IF(SpecialPass)nPass = 2
   ICHARGE=8; Rxyz(1)=0.762266389544351E0_realk; Rxyz(2)=0.983877565461657E0_realk;Rxyz(3)=0.624979148086261E0_realk; 
   call build_unittest_atomicmolecule(atomicmolecule(2),ICHARGE,Rxyz,nPass,lupri)
   IF(iPass.EQ.2)nPass = 3
   IF(SpecialPass)nPass = 1
   ICHARGE=9; Rxyz(1)=0.736938390171405E0_realk; Rxyz(2)=1.108186821166992E0_realk;Rxyz(3)=0.713699152299640E0_realk; 
   call build_unittest_atomicmolecule(atomicmolecule(3),ICHARGE,Rxyz,nPass,lupri)
   IF(iPass.EQ.2)nPass = 2
   IF(SpecialPass)nPass = 1
   ICHARGE=17; Rxyz(1)=0.574178167982901E0_realk; Rxyz(2)=1.086728949849219E0_realk;Rxyz(3)=0.913548738546455E0_realk; 
   call build_unittest_atomicmolecule(atomicmolecule(4),ICHARGE,Rxyz,nPass,lupri)
   
   SETTING%MOLECULE(1)%p => atomicmolecule(1)
   SETTING%MOLECULE(2)%p => atomicmolecule(2)
   SETTING%MOLECULE(3)%p => atomicmolecule(3)
   SETTING%MOLECULE(4)%p => atomicmolecule(4)
   !      Write(lupri,*)'Atomicmolecule 1'
   !      call print_mol(atomicmolecule(1),lupri)
   !      Write(lupri,*)'Atomicmolecule 2'
   !      call print_mol(atomicmolecule(2),lupri)
   !      Write(lupri,*)'Atomicmolecule 3'
   !      call print_mol(atomicmolecule(3),lupri)
   !      Write(lupri,*)'Atomicmolecule 4'
   !      call print_mol(atomicmolecule(4),lupri)
   
   allocate(UNITTESTBASIS(4))
   BASISLABEL='REGULAR  '
   
   FAIL = .FALSE.
   spherical = .TRUE.!.FALSE.
   nullify(integralsII)
   
   basisloop: do iBasis1Q = 1,nBasisA
    iBasis1 = iBasisA(iBasis1Q)
    do iBasis2Q = 1,nBasisB
     iBasis2 = iBasisB(iBasis2Q)
     do iBasis3Q = 1,nBasisC
      iBasis3 = iBasisC(iBasis3Q)
      do iBasis4Q = 1,nBasisD
      iBasis4 = iBasisD(iBasis4Q)
         !Skip forward 
!         IF(iBasis4+(iBasis3-1)*nBasis+(iBasis2-1)*nBasis*nBasis+&
!              & (iBasis1-1)*nBasis*nBasis*nBasis.LT.14) CYCLE

         print*,'iBasis:',iBasis4Q+(iBasis3Q-1)*nBasisD+(iBasis2Q-1)*nBasisD*nBasisC+&
              & (iBasis1Q-1)*nBasisD*nBasisC*nBasisB,'of',nBasisA*nBasisB*nBasisC*nbasisD
       ibasiselm(1) = iBasis1
       ibasiselm(2) = iBasis2
       ibasiselm(3) = iBasis3
       ibasiselm(4) = iBasis4
       do A = 1,4       
          BASISSETNAME(1:20) = BASISTYPE(iBasiselm(A))
          CALL Build_basis(LUPRI,IPRINT,&
               &SETTING%MOLECULE(A)%p,UNITTESTBASIS(A)%REGULAR,LIBRARY,&
               &BASISLABEL,.FALSE.,.FALSE.,doprint,spherical,BASISSETNAME)
          SETTING%BASIS(A)%p => UNITTESTBASIS(A)
          call determine_nbast2(SETTING%MOLECULE(A)%p,SETTING%BASIS(A)%p%REGULAR,spherical,.FALSE.,nbast(A))
       enddo
       dim1 = nbast(1); dim2 = nbast(2); dim3 = nbast(3); dim4 = nbast(4)
       !due to current code restrictions
!       IF(dim2.GT.dim1)CYCLE
!       IF(dim3/nPass.GT.dim1)CYCLE
!       IF(dim4.GT.dim3/nPass)CYCLE
!       IF(dim3/nPass+dim4.GT.dim1+dim2)CYCLE
       
       write(lupri,'(A,A,A,A,A,A,A,A,A)')'BASIS(',BASISTYPE(iBasis1),',',BASISTYPE(iBasis2),&
            & ',',BASISTYPE(iBasis3),',',BASISTYPE(iBasis4),')'
       
       write(lupri,*)'dim:',dim1,dim2,dim3,dim4
       
       Setting%sameMol = .FALSE.
       Setting%sameFrag = .FALSE.
       if(associated(integralsII))THEN
          call mem_dealloc(integralsII)
       ENDIF
       call mem_alloc(integralsII,dim1,dim2,dim3,dim4)
       integralsII = 0.0E0_realk
!       setting%scheme%intprint = 1000
       savedospherical = setting%scheme%dospherical
       setting%scheme%dospherical = spherical
       setting%scheme%OD_SCREEN = .FALSE.
       setting%scheme%CS_SCREEN = .FALSE.
       setting%scheme%PS_SCREEN = .FALSE.
!       WRITE(lupri,*)'ThermiteDriver'
       call II_get_4center_eri(LUPRI,LUERR,SETTING,integralsII,dim1,dim2,dim3,dim4,dirac)
       !   print*,'integralsII',integralsII
       setting%scheme%dospherical = savedospherical
       setting%scheme%OD_SCREEN = .TRUE.
       setting%scheme%CS_SCREEN = .TRUE.
       setting%scheme%PS_SCREEN = .TRUE.
       setting%scheme%intprint = 0
       Setting%sameFrag = .FALSE.
       Setting%sameMol = .FALSE.

!       print*,'dim1,dim2,dim3,dim4',dim1,dim2,dim3,dim4
       call mem_alloc(integralsIchor,dim1,dim2,dim3,dim4)
!       integralsIchor = 0.0E0_realk
!          setting%scheme%intprint = 1000
!       WRITE(lupri,*)'IchorDriver'
!       CALL FLUSH(LUPRI)
!       print*,'call ICHOR WITH BASIS(',BASISTYPE(iBasis1),',',BASISTYPE(iBasis2),',',&
!            & BASISTYPE(iBasis3),',',BASISTYPE(iBasis4),')'
       SameMOL = .FALSE.       
       call SCREEN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,INTSPEC,SameMOL)
       call MAIN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,integralsIchor,intspec,.TRUE.,1,1,1,1,1,1,1,1)!,spherical)
       call FREE_SCREEN_ICHORERI
!       print*,'DONE call ICHOR WITH BASIS'
       !setting%scheme%intprint = 0
       write(lupri,'(A,A,A,A,A,A,A,A,A)')'BASIS(',BASISTYPE(iBasis1),',',BASISTYPE(iBasis2),',',&
            & BASISTYPE(iBasis3),',',BASISTYPE(iBasis4),') TESTING'
       FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .FALSE.
       DO D=1,dim4
          DO C=1,dim3
             DO B=1,dim2
                DO A=1,dim1
                   IF(ABS(integralsII(A,B,C,D)).GT.1.0E-10_realk)THEN
                      IF(ABS(integralsII(A,B,C,D)-integralsIchor(A,B,C,D)).GT. &
                           & 1.0E-10_realk/(ABS(integralsII(A,B,C,D))))THEN
                         FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
                         write(lupri,'(A,ES16.8)')'THRESHOLD=',1.0E-10_realk/(ABS(integralsII(A,B,C,D)))
                         write(lupri,'(A,4I4)')'ELEMENTS: (A,B,C,D)=',A,B,C,D
                         write(lupri,'(A,ES16.8)')'integralsII(A,B,C,D)   ',integralsII(A,B,C,D)
                         write(lupri,'(A,ES16.8)')'integralsIchor(A,B,C,D)',integralsIchor(A,B,C,D)
                         write(lupri,'(A,ES16.8)')'DIFF                   ',&
                              & ABS(integralsII(A,B,C,D)-integralsIchor(A,B,C,D))
                         call lsquit('ERROR',-1)
                      ENDIF
                   ELSE
                      IF(ABS(integralsII(A,B,C,D)-integralsIchor(A,B,C,D)).GT. &
                           & 1.0E-10_realk)THEN
                         FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
                         write(lupri,'(A,ES16.8)')'THRESHOLD=',1.0E-10_realk/(ABS(integralsII(A,B,C,D)))
                         write(lupri,'(A,4I4)')'ELEMENTS: (A,B,C,D)=',A,B,C,D
                         write(lupri,'(A,ES16.8)')'integralsII(A,B,C,D)   ',integralsII(A,B,C,D)
                         write(lupri,'(A,ES16.8)')'integralsIchor(A,B,C,D)',integralsIchor(A,B,C,D)
                         write(lupri,'(A,ES16.8)')'DIFF                   ',&
                              & ABS(integralsII(A,B,C,D)-integralsIchor(A,B,C,D))
                         call lsquit('ERROR',-1)
                         !                   ELSE
                         !                      write(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8,A,ES16.8)')&
                         !                           & 'SUCCESS(',A,',',B,',',C,',',D,')=',integralsIchor(A,B,C,D),'  DIFF',&
                         !                           & ABS(integralsII(A,B,C,D)-integralsIchor(A,B,C,D))
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF(FAIL(iBasis1,ibasis2,ibasis3,ibasis4))THEN
          WRITE(lupri,'(A)')'CALC FAILED'
       ENDIF
       call mem_dealloc(integralsII)
       call mem_dealloc(integralsIchor)
       call free_basissetinfo(UNITTESTBASIS(1)%REGULAR)
       call free_basissetinfo(UNITTESTBASIS(2)%REGULAR)
       call free_basissetinfo(UNITTESTBASIS(3)%REGULAR)
       call free_basissetinfo(UNITTESTBASIS(4)%REGULAR)
       Setting%sameMol = .TRUE.
       Setting%sameFrag = .TRUE.
    ENDDO
   ENDDO
  ENDDO
 ENDDO basisloop

 write(lupri,'(A,I4)')'Summary of Unit Test for nPasses=',ipass
 do iBasis1Q = 1,nBasisA
  iBasis1 = iBasisA(iBasis1Q)
  do iBasis2Q = 1,nBasisB
   iBasis2 = iBasisB(iBasis2Q)
   do iBasis3Q = 1,nBasisC
    iBasis3 = iBasisC(iBasis3Q)
    do iBasis4Q = 1,nBasisD
     iBasis4 = iBasisD(iBasis4Q)
       !warning this only works if number is realted to angmom
!     IF(iBasis2.GT.iBasis1)CYCLE
!     IF(iBasis3.GT.iBasis1)CYCLE
!     IF(iBasis4.GT.iBasis3)CYCLE
!     IF(iBasis3+iBasis4.GT.iBasis1+iBasis2)CYCLE  

     IF(FAIL(iBasis1,ibasis2,ibasis3,ibasis4)) THEN
        write(lupri,'(A,A,A,A,A,A,A,A,A,I1,A)')'BASIS(',BASISTYPE(iBasis1)(10:15),',',&
             & BASISTYPE(iBasis2)(10:15),',',BASISTYPE(iBasis3)(10:15),',',&
             & BASISTYPE(iBasis4)(10:15),',',ipass,') FAILED'
        ALLPASS = .FALSE.
     ELSE
        write(lupri,'(A,A,A,A,A,A,A,A,A,I1,A)')'BASIS(',BASISTYPE(iBasis1)(10:15),',',&
             & BASISTYPE(iBasis2)(10:15),',',BASISTYPE(iBasis3)(10:15),',',&
             & BASISTYPE(iBasis4)(10:15),',',ipass,') SUCCESSFUL'
     ENDIF

!     write(lupri,'(A,A,A,A,A,A,A,A,A,I1,A)')'CRIT1=`$GREP "BASIS\(',BASISTYPE(iBasis1)(10:15),',',BASISTYPE(iBasis2)(10:15),',',BASISTYPE(iBasis3)(10:15),',',BASISTYPE(iBasis4)(10:15),',',ipass,'\) SUCCESSFUL" $log | wc -l`'
!     write(lupri,'(A,I3,A)')'TEST[',itest,']=`expr  $CRIT1`'
!     write(lupri,'(A,I3,A)')'CTRL[',itest,']=1'
!     write(lupri,'(A,I3,A,A,A,A,A,A,A,A,A,I1,A)')'ERROR[',itest,']="BASIS(',BASISTYPE(iBasis1)(10:15),',',BASISTYPE(iBasis2)(10:15),',',BASISTYPE(iBasis3)(10:15),',',BASISTYPE(iBasis4)(10:15),',',ipass,') NOT CORRECT -"'
!     itest = itest + 1
    enddo
   enddo
  enddo
 enddo
 call free_moleculeinfo(atomicmolecule(1))
 call free_moleculeinfo(atomicmolecule(2))
 call free_moleculeinfo(atomicmolecule(3))
 call free_moleculeinfo(atomicmolecule(4))
 deallocate(atomicmolecule)
 deallocate(UNITTESTBASIS)
enddo !ipass 
iprint=0
setting%scheme%intprint = 0

SETTING%MOLECULE(1)%p => Originalmolecule
SETTING%MOLECULE(2)%p => Originalmolecule
SETTING%MOLECULE(3)%p => Originalmolecule
SETTING%MOLECULE(4)%p => Originalmolecule

SETTING%BASIS(1)%p => originalBasis
SETTING%BASIS(2)%p => originalBasis
SETTING%BASIS(3)%p => originalBasis
SETTING%BASIS(4)%p => originalBasis

IF(ALLPASS)THEN
   WRITE(lupri,'(A)')'Ichor Integrals tested against Thermite: SUCCESSFUL'
   print*,'Ichor Integrals tested against Thermite: SUCCESSFUL'
ELSE
   WRITE(lupri,'(A)')'Ichor Integrals tested against Thermite: FAILED'
   print*,'Ichor Integrals tested against Thermite: FAILED'
ENDIF
WRITE(lupri,*)'done II_test_Ichor'
!call lsquit('II_test_Ichor done',-1)

!WRITE(lupri,*)'After IchorUnitTest'
!call debug_mem_stats(LUPRI)

#endif
END SUBROUTINE II_unittest_Ichor

subroutine build_unittest_atomicmolecule(atomicmolecule,ICHARGE,Rxyz,nAtoms,lupri)
implicit none
type(moleculeinfo) :: atomicmolecule
real(realk)        :: Rxyz(3)
integer,intent(in) :: ICHARGE,lupri,nAtoms
#ifdef VAR_ICHOR
character(len=22) :: label
character(len=4) :: Name
integer :: I
do I=1,22
   label(I:I)=' '
enddo
label(1:11)='UnittestMol'
do I=1,4
   name(I:I)=' '
enddo
name(1:1)='U'
IF(icharge.LT.10)THEN
   write(label(12:12),'(I1)') icharge
   write(name(2:2),'(I1)') icharge
ELSEIF(icharge.LT.100)THEN
   write(label(12:13),'(I2)') icharge
   write(name(2:3),'(I2)') icharge
ELSEIF(icharge.LT.1000)THEN
   write(name(2:4),'(I3)') icharge
   write(label(12:14),'(I3)') icharge
ENDIF
atomicmolecule%label = label
call mem_alloc(atomicmolecule%ATOM,nAtoms)
atomicmolecule%nAtoms = nAtoms
atomicmolecule%nAtomsNPC = nAtoms
atomicmolecule%nelectrons = ICHARGE*nAtoms
atomicmolecule%charge = 0.0E0_realk
atomicmolecule%nbastREG = 0
atomicmolecule%nbastAUX = 0
atomicmolecule%nbastCABS = 0
atomicmolecule%nbastJK = 0
atomicmolecule%nbastVAL = 0
atomicmolecule%nprimbastREG = 0
atomicmolecule%nprimbastAUX = 0
atomicmolecule%nprimbastCABS = 0
atomicmolecule%nprimbastJK = 0
atomicmolecule%nprimbastVAL = 0
atomicmolecule%pointMolecule = .FALSE.
do I=1,nAtoms
   atomicmolecule%ATOM(I)%Isotope = 1 
   atomicmolecule%ATOM(I)%Name = Name
   atomicmolecule%ATOM(I)%Mass = 0.0E0_realk    
   atomicmolecule%ATOM(I)%CovRad = 0.0E0_realk      
   atomicmolecule%ATOM(I)%Frag = 0.0E0_realk
   atomicmolecule%ATOM(I)%CENTER(1) = Rxyz(1)*I
   atomicmolecule%ATOM(I)%CENTER(2) = Rxyz(2)*I
   atomicmolecule%ATOM(I)%CENTER(3) = Rxyz(3)*I
   atomicmolecule%ATOM(I)%Atomic_number = 0 
   atomicmolecule%ATOM(I)%Charge = Icharge       
   atomicmolecule%ATOM(I)%nbasis=1
   
   atomicmolecule%ATOM(I)%basislabel(1) = 'None'
   atomicmolecule%ATOM(I)%basislabel(2) = 'None'
   atomicmolecule%ATOM(I)%basislabel(3) = 'None'
   atomicmolecule%ATOM(I)%basislabel(4) = 'None'
   atomicmolecule%ATOM(I)%basislabel(5) = 'None'
   atomicmolecule%ATOM(I)%Basisindex(1:5) = 0 
   atomicmolecule%ATOM(I)%IDtype(1:5) = 0 
   atomicmolecule%ATOM(I)%Phantom = .FALSE.
   atomicmolecule%ATOM(I)%Pointcharge = .FALSE. 
   atomicmolecule%ATOM(I)%nContOrbREG=0 
   atomicmolecule%ATOM(I)%nPrimOrbREG =0 
   atomicmolecule%ATOM(I)%nContOrbAUX =0 
   atomicmolecule%ATOM(I)%nPrimOrbAUX =0 
   atomicmolecule%ATOM(I)%nContOrbCABS =0 
   atomicmolecule%ATOM(I)%nPrimOrbCABS =0 
   atomicmolecule%ATOM(I)%nContOrbJK =0 
   atomicmolecule%ATOM(I)%nPrimOrbJK =0 
   atomicmolecule%ATOM(I)%nContOrbVAL =0 
   atomicmolecule%ATOM(I)%nPrimOrbVAL =0 
   atomicmolecule%ATOM(I)%molecularIndex =0 
ENDDO

#endif
end subroutine build_unittest_atomicmolecule

End MODULE IntegralInterfaceIchorMod

