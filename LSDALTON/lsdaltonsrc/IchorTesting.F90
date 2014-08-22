MODULE IntegralInterfaceIchorMod
  use precision
  use TYPEDEFTYPE, only: LSSETTING, LSINTSCHEME, LSITEM, integralconfig,&
       & BASISSETLIBRARYITEM
  use basis_type, only: free_basissetinfo
  use basis_typetype,only: BASISSETINFO,BASISINFO,RegBasParam,nBasisBasParam
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
  use IntegralInterfaceMOD, only: II_get_4center_eri, ii_get_2int_screenmat
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
SUBROUTINE II_unittest_Ichor(LUPRI,LUERR,SETTING,DebugIchorOption,dolink)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,DebugIchorOption
logical :: dolink
!
#ifdef VAR_ICHOR
real(realk),pointer   :: integralsII(:,:,:,:),integralsIchor(:,:,:,:)
real(realk),pointer   :: KmatII(:,:,:),KmatIchor(:,:,:)
integer :: dim1,dim2,dim3,dim4,A,B,C,D,iprint,nbast(4),ibasiselm(4)
integer :: iBasis1,ibasis2,ibasis3,ibasis4,icharge,nbasis,nPass,ipass,itest
logical :: dirac,doprint,debug
TYPE(MOLECULEINFO),pointer :: Originalmolecule
TYPE(MOLECULEINFO),pointer :: atomicmolecule(:)
TYPE(BASISSETLIBRARYITEM) :: LIBRARY(nBasisBasParam)
CHARACTER(len=9)     :: BASISLABEL
TYPE(BASISINFO),pointer :: unittestBASIS(:)
TYPE(BASISINFO),pointer :: originalBASIS
CHARACTER(len=80)    :: BASISSETNAME
CHARACTER(len=20)    :: BASISTYPE(10)
real(realk)          :: Rxyz(3)
type(lsmatrix)       :: FINALVALUE(2)
logical      :: spherical,savedospherical,SpecialPass,LHS
logical      :: FAIL(10,10,10,10),ALLPASS,SameMOL,TestScreening
logical      :: MoTrans,NoSymmetry
Character    :: intSpec(5)
integer :: iBasis1Q,iBasis2Q,iBasis3Q,iBasis4Q
integer :: nBasisA,nBasisB,nBasisC,nBasisD,iPassStart,iPassEnd
integer,pointer :: iBasisA(:),iBasisB(:),iBasisC(:),iBasisD(:)
real(realk),pointer :: IIBATCHGAB(:,:),BATCHGAB(:,:),Dmat(:,:,:)
integer :: nBatchA,nBatchB,iAO,iatom,jatom,i,j,idx,jdx,nDmat,idmat
type(matrix) :: GAB

NoSymmetry = .FALSE. !activate permutational symmetry
MoTrans=.FALSE.
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
TestScreening = .FALSE.
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
CASE(6)
   !all types and all passes - but only screening
   nbasisA = 9; nbasisB = 9; nbasisC = 1; nbasisD = 1
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
   iBasisC(1) = 1
   iBasisD(1) = 1
   TestScreening = .TRUE.
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
   DO iao=1,4
      SETTING%MOLECULE(iao)%p => atomicmolecule(iao)
      setting%fragment(iao)%p => setting%molecule(iao)%p
   ENDDO
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

!         print*,'iBasis:',iBasis4Q+(iBasis3Q-1)*nBasisD+(iBasis2Q-1)*nBasisD*nBasisC+&
!              & (iBasis1Q-1)*nBasisD*nBasisC*nBasisB,'of',nBasisA*nBasisB*nBasisC*nbasisD
       ibasiselm(1) = iBasis1
       ibasiselm(2) = iBasis2
       ibasiselm(3) = iBasis3
       ibasiselm(4) = iBasis4
       do A = 1,4       
          BASISSETNAME(1:20) = BASISTYPE(iBasiselm(A))
          CALL Build_basis(LUPRI,IPRINT,&
               &SETTING%MOLECULE(A)%p,UNITTESTBASIS(A)%BINFO(RegBasParam),LIBRARY,&
               &BASISLABEL,.FALSE.,.FALSE.,doprint,spherical,RegBasParam,BASISSETNAME)
          SETTING%BASIS(A)%p => UNITTESTBASIS(A)
          call determine_nbast2(SETTING%MOLECULE(A)%p,SETTING%BASIS(A)%p%BINFO(RegBasParam),spherical,.FALSE.,nbast(A))
       enddo
       dim1 = nbast(1); dim2 = nbast(2); dim3 = nbast(3); dim4 = nbast(4)
       !due to current code restrictions
!       IF(dim2.GT.dim1)CYCLE
!       IF(dim3/nPass.GT.dim1)CYCLE
!       IF(dim4.GT.dim3/nPass)CYCLE
!       IF(dim3/nPass+dim4.GT.dim1+dim2)CYCLE

       IF(.NOT.TestScreening)THEN
          write(lupri,'(A,A,A,A,A,A,A,A,A)')'BASIS(',BASISTYPE(iBasis1),',',BASISTYPE(iBasis2),&
               & ',',BASISTYPE(iBasis3),',',BASISTYPE(iBasis4),')'
          
          write(lupri,*)'dim:',dim1,dim2,dim3,dim4
          
          Setting%sameMol = .FALSE.
          Setting%sameFrag = .FALSE.

          IF(doLink)THEN
             nDmat = 2
             call mem_alloc(Dmat,dim2,dim4,nDmat)
             call mem_alloc(KmatII,dim1,dim3,nDmat)
             KmatII = 0.0E0_realk
             call mem_alloc(KmatIchor,dim1,dim3,nDmat)
             KmatIchor = 0.0E0_realk
             call MakeRandomDmat(Dmat,dim2,dim4,nDmat) !NON SYMMETRIC DMAT
!             WRITE(lupri,*)'The Density Matrix: '
!             call output(Dmat(:,:,2),1,dim2,1,dim4,dim2,dim4,1,lupri)

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
             call BuildExchange(integralsII,dim1,dim2,dim3,dim4,nDmat,Dmat,KmatII)
             call mem_dealloc(integralsII)

             SameMOL = .FALSE.       
             call SCREEN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,INTSPEC,SameMOL)
             call MAIN_LINK_ICHORERI_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,&
                  & nDmat,KmatIchor,Dmat,intspec,.TRUE.,1,1,1,1,1,1,1,1)
             call FREE_SCREEN_ICHORERI

!             WRITE(lupri,*)'The Exchange Matrix: '
!             call output(KmatII(:,:,2),1,dim1,1,dim3,dim1,dim3,1,lupri)
!             WRITE(lupri,*)'The Exchange Matrix: '
!             call output(KmatIchor(:,:,2),1,dim1,1,dim3,dim1,dim3,1,lupri)
!             WRITE(lupri,*)'The Density Matrix: '
!             call output(Dmat(:,:,2),1,dim2,1,dim4,dim2,dim4,1,lupri)

             write(lupri,'(A,A,A,A,A,A,A,A,A)')'LINK BASIS(',BASISTYPE(iBasis1),',',&
                  & BASISTYPE(iBasis2),',',BASISTYPE(iBasis3),',',BASISTYPE(iBasis4),') TESTING'
             FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .FALSE.
             DO idmat=1,nDmat
              DO C=1,dim3
               DO A=1,dim1
                IF(ABS(KmatII(A,C,idmat)).GT.1.0E-10_realk)THEN
                 IF(ABS(KmatII(A,C,idmat)-KmatIchor(A,C,idmat)).GT. &
                      & 1.0E-10_realk/(ABS(KmatII(A,C,idmat))))THEN
                    FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
                    write(lupri,'(A,ES16.8)')'THRESHOLD=',1.0E-10_realk/(ABS(KmatII(A,C,iDmat)))
                    write(lupri,'(A,4I4)')'ELEMENTS: (A,C,iDmat)= ',A,C,iDmat
                    write(lupri,'(A,ES16.8)')'Kmat(A,C,iDmat)     ',KmatII(A,C,iDmat)
                    write(lupri,'(A,ES16.8)')'KmatIchor(A,C,iDmat)',KmatIchor(A,C,iDmat)
                    write(lupri,'(A,ES16.8)')'DIFF                   ',&
                         & ABS(KmatII(A,C,iDmat)-KmatIchor(A,C,iDMat))
                    call lsquit('ERROR',-1)
                 ENDIF
                ELSE
                 IF(ABS(KmatII(A,C,iDmat)-KmatIchor(A,C,iDmat)).GT. &
                      & 1.0E-10_realk)THEN
                    FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
                    write(lupri,'(A,ES16.8)')'THRESHOLD=',1.0E-10_realk/(ABS(KmatII(A,C,iDmat)))
                    write(lupri,'(A,4I4)')'ELEMENTS: (A,C,iDmat)=   ',A,C,iDmat
                    write(lupri,'(A,ES16.8)')'KmatII(A,1,C,iDmat)   ',KmatII(A,C,iDmat)
                    write(lupri,'(A,ES16.8)')'KmatIchor(A,C,iDmat)  ',KmatIchor(A,C,iDmat)
                    write(lupri,'(A,ES16.8)')'DIFF                   ',&
                         & ABS(KmatII(A,C,iDmat)-KmatIchor(A,C,iDmat))
                    call lsquit('ERROR',-1)
                 ENDIF
                ENDIF
               ENDDO
              ENDDO
             ENDDO
             IF(FAIL(iBasis1,ibasis2,ibasis3,ibasis4))THEN
                WRITE(lupri,'(A)')'LINK CALC FAILED'
             ENDIF
          ELSE
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
             call MAIN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,&
                  & integralsIchor,intspec,.TRUE.,1,1,1,1,1,1,1,1,&
                  & MoTrans,dim1,dim2,dim3,dim4,NoSymmetry)
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

          ENDIF
       ELSE
          !test screening integrals
          call mat_init(GAB,dim1,dim2)
          Setting%sameMol = .FALSE.
          Setting%sameFrag = .FALSE.
          savedospherical = setting%scheme%dospherical
          setting%scheme%dospherical = spherical
          SETTING%batchindex=0
          call II_get_2int_ScreenMat(LUPRI,LUERR,SETTING,GAB)
          setting%scheme%dospherical = savedospherical
          Setting%sameFrag = .FALSE.
          Setting%sameMol = .FALSE.
!          WRITE(lupri,*)'OUTPUT'
!          call mat_print(GAB,1,dim1,1,dim2,lupri)

          SameMOL = .FALSE.       
          call SCREEN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,INTSPEC,SameMOL)
          LHS = .TRUE.
          call SCREEN_ICHORERI_RETRIEVE_GABDIM(LUPRI,IPRINT,setting,nBatchA,nBatchB,LHS)
          call mem_alloc(BATCHGAB,nBatchA,nBatchB)
          call SCREEN_ICHORERI_RETRIEVE_GAB(LUPRI,IPRINT,setting,nBatchA,nBatchB,LHS,BATCHGAB)
!          WRITE(lupri,*)'The LHS BatchGab  DIM:',nBatchA,nBatchB,'atoms=',&
!               & atomicmolecule(1)%nAtoms,atomicmolecule(2)%nAtoms
!          call OUTPUT(BATCHGAB,1,nBatchA,1,nBatchB,nBatchA,nBatchB,1,lupri)
          IF(atomicmolecule(1)%nAtoms.NE.nBatchA)call lsquit('Error dim1 in IchorTestingScreen',-1)
          IF(atomicmolecule(2)%nAtoms.NE.nBatchB)call lsquit('Error dim2 in IchorTestingScreen',-1)
          call mem_alloc(IIBATCHGAB,nBatchA,nBatchB)
          JDX = 0
          DO JAtom = 1,atomicmolecule(2)%nAtoms
             DO IAtom = 1,atomicmolecule(1)%nAtoms
                IIBATCHGAB(Iatom,Jatom) = 0.0E0_realk
             ENDDO
          ENDDO
          DO JAtom = 1,atomicmolecule(2)%nAtoms
             DO J=1,dim2/atomicmolecule(2)%nAtoms
                JDX = JDX + 1 
                IDX = 0
                DO IAtom = 1,atomicmolecule(1)%nAtoms
                   DO I=1,dim1/atomicmolecule(1)%nAtoms
                      IDX = IDX + 1 
                      IIBATCHGAB(Iatom,Jatom) = MAX(IIBATCHGAB(Iatom,Jatom),GAB%elms(IDX + (JDX-1)*dim1))
                   ENDDO
                ENDDO                     
             ENDDO
          ENDDO

          FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .FALSE.
          DO JAtom = 1,atomicmolecule(2)%nAtoms
             DO IAtom = 1,atomicmolecule(1)%nAtoms
                IF(ABS(IIBATCHGAB(Iatom,Jatom)-BATCHGAB(Iatom,Jatom)).GT.1.0E-10_realk)THEN
!                   print*,'ABS(IIBATCHGAB(Iatom,Jatom)-BATCHGAB(Iatom,Jatom))',ABS(IIBATCHGAB(Iatom,Jatom)-BATCHGAB(Iatom,Jatom))
                   FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.

                   WRITE(lupri,*)'ABS(IIBATCHGAB(Iatom,Jatom)-BATCHGAB(Iatom,Jatom))',&
                        & ABS(IIBATCHGAB(Iatom,Jatom)-BATCHGAB(Iatom,Jatom))
                   WRITE(lupri,*)'Thermite Matrix OUTPUT'
                   call mat_print(GAB,1,dim1,1,dim2,lupri)
                   WRITE(lupri,*)'The Ichor BatchGab'
                   call OUTPUT(BATCHGAB,1,nBatchA,1,nBatchB,nBatchA,nBatchB,1,lupri)
                   WRITE(lupri,*)'The Thermite BatchGab'
                   call OUTPUT(IIBATCHGAB,1,nBatchA,1,nBatchB,nBatchA,nBatchB,1,lupri)
                   call lsquit('Error in Ichor Screening Testing',-1)
!                ELSE
!                   print*,'ABS(IIBATCHGAB(Iatom,Jatom)-BATCHGAB(Iatom,Jatom))',ABS(IIBATCHGAB(Iatom,Jatom)-BATCHGAB(Iatom,Jatom))
                ENDIF
             ENDDO
          ENDDO
          call mat_free(GAB)
          call mem_dealloc(IIBATCHGAB)
          call mem_dealloc(BATCHGAB)
          call FREE_SCREEN_ICHORERI
       ENDIF
       call free_basissetinfo(UNITTESTBASIS(1)%BINFO(RegBasParam))
       call free_basissetinfo(UNITTESTBASIS(2)%BINFO(RegBasParam))
       call free_basissetinfo(UNITTESTBASIS(3)%BINFO(RegBasParam))
       call free_basissetinfo(UNITTESTBASIS(4)%BINFO(RegBasParam))
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

DO iao=1,4
   SETTING%MOLECULE(iao)%p => Originalmolecule
   setting%fragment(iao)%p => setting%molecule(iao)%p
ENDDO

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

subroutine MakeRandomDmat(Dmat,dim2,dim4,nDmat)
  implicit none
  integer :: dim2,dim4,nDmat
  real(realk) :: Dmat(dim2*dim4*nDmat)
  CALL RANDOM_SEED 
  CALL RANDOM_NUMBER(Dmat)
end subroutine MakeRandomDmat

subroutine BuildExchange(integralsII,dim1,dim2,dim3,dim4,nDmat,Dmat,Kmat)
  implicit none
  integer :: dim1,dim2,dim3,dim4,nDmat
  real(realk) :: integralsII(dim1,dim2,dim3,dim4)
  real(realk) :: Dmat(dim2,dim4,nDmat)
  real(realk) :: Kmat(dim1,dim3,nDmat)
  !
  integer :: A,B,C,D,iDmat
  do iDmat = 1,nDmat
     do D=1,dim4
        do C=1,dim3
           do B=1,dim2
              do A=1,dim1
                 Kmat(A,C,iDmat) = Kmat(A,C,iDmat) + integralsII(A,B,C,D)*Dmat(B,D,idmat)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine BuildExchange

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
nullify(atomicmolecule%SubsystemLabel)
atomicmolecule%nSubsystems = 0
atomicmolecule%nAtoms = nAtoms
atomicmolecule%nAtomsNPC = nAtoms
atomicmolecule%nelectrons = ICHARGE*nAtoms
atomicmolecule%charge = 0.0E0_realk
atomicmolecule%nbastREG = 0
atomicmolecule%nbastAUX = 0
atomicmolecule%nbastCABS = 0
atomicmolecule%nbastJK = 0
atomicmolecule%nbastADMM = 0
atomicmolecule%nbastVAL = 0
atomicmolecule%nprimbastREG = 0
atomicmolecule%nprimbastAUX = 0
atomicmolecule%nprimbastCABS = 0
atomicmolecule%nprimbastJK = 0
atomicmolecule%nprimbastADMM = 0
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
   atomicmolecule%ATOM(I)%nContOrbADMM =0 
   atomicmolecule%ATOM(I)%nPrimOrbADMM =0 
   atomicmolecule%ATOM(I)%nContOrbVAL =0 
   atomicmolecule%ATOM(I)%nPrimOrbVAL =0 
   atomicmolecule%ATOM(I)%molecularIndex =0 
   atomicmolecule%ATOM(I)%SubSystemIndex =0 
ENDDO

#endif
end subroutine build_unittest_atomicmolecule

End MODULE IntegralInterfaceIchorMod

