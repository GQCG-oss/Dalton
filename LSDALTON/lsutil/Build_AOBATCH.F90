!> @file
!> Module containing subroutines for building AOitem (information about AO-shells/batches)
!> \author T. Kjaergaard
!> \date 2008
MODULE BUILDAOBATCH
  use precision
  use TYPEDEF
  use TYPEDEFTYPE
  use lsmatrix_operations_dense
  use memory_handling
  use molecule_type
  use AO_type
  use AO_typetype

TYPE ORGANIZE
   INTEGER                    :: type,Segment(maxAOangmom),atom
   INTEGER                    :: exponentindex,angmoms,ANGMOM(maxAOangmom)
   INTEGER                    :: CC(maxAOangmom)
END TYPE ORGANIZE
  
TYPE AOorganizer
   INTEGER                    :: nbatches
   TYPE(ORGANIZE),pointer     :: ORG(:)
END TYPE AOORGANIZER

CONTAINS
!> \brief builds the AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> The most general routine which builds the AOitem from the BASISSET
!> in the BASISINFO and from the molecule 
!>
SUBROUTINE BUILD_AO(LUPRI,SCHEME,IPRINT,MOLECULE,BASISINFO,AO&
     &,UNCONTRACTED,INTNRM,EXTEND,NORMA)
implicit none
!> the logical unit number for the output file
INTEGER,intent(in)        :: LUPRI
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME),intent(in) :: SCHEME
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
!> the AOitem to be build
TYPE(AOITEM),intent(inout):: AO
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)        :: IPRINT
!> if the AOitem should be uncontracted
LOGICAL,intent(in)        :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL,intent(in)        :: INTNRM 
!> extend (add to) the input AO 
LOGICAL,intent(in)        :: EXTEND
!> if AOS should be the normalized orbitals (default true)
LOGICAL,optional          :: NORMA
!
TYPE(AOITEM),pointer      :: AOmodel(:)
TYPE(AOITEM)              :: AOTMP
TYPE(AOorganizer)         :: AOorganize 
INTEGER                   :: nsegments,nAngmom,nbatches,L,A,B
INTEGER                   :: UNIQUEMATRIXES,UNIQUEEXPONENTS,SUM
INTEGER                   :: I,J,K,type,nPrimitives,nContracted
INTEGER                   :: R,aobatches,UNIQUE2
INTEGER                   :: orbitalIndex,tmporbitalindex
INTEGER                   :: GHOSTFUNCS,AOmodelbat,nrow,ncol
INTEGER,pointer           :: MODELTYPES(:)
INTEGER                   :: tmpprimorbitalindex,primorbitalindex,icharge
INTEGER                   :: nbat,SUM1,nMODELEXP,itype,offsetAT
INTEGER,pointer           :: MODELEXP(:)
LOGICAL                   :: NOFAMILY,OLDATOM,NOSEGMENT,NORM
IF(BASISINFO%natomtypes.EQ. 0)THEN
   print*,'Error BUILD_AO called with empty basis'
   print*,'BASISINFO%label        ',BASISINFO%label
   print*,'BASISINFO%nChargeindex ',BASISINFO%nChargeindex
   print*,'BASISINFO%nbast        ',BASISINFO%nbast
   print*,'BASISINFO%nprimbast    ',BASISINFO%nprimbast
   CALL LSQUIT('Error BUILD_AO called with empty basis',lupri)
ENDIF
IF(EXTEND)THEN
   IF(IPRINT .GT. 20) THEN
      WRITE(LUPRI,*)'The Original AO '
      CALL PRINT_AO(LUPRI,AO)
   ENDIF

   call nullifyAOITEM(AOTMP)
   CALL COPY_AOITEM(AO,AOTMP)
   call free_AOitem(lupri,AO)

   IF(IPRINT .GT. 20) THEN
      WRITE(LUPRI,*)'The Copyied AO '
      CALL PRINT_AO(LUPRI,AOTMP)      
   ENDIF
ENDIF
call nullifyAOITEM(AO)
IF(PRESENT(NORMA))THEN
   NORM=NORMA
ELSE
   NORM=.TRUE.
ENDIF
J=0
DO I=1,MOLECULE%natoms  
   IF(MOLECULE%ATOM(I)%pointcharge)CYCLE 
   J=J+1
ENDDO
AO%natoms = J
IF(EXTEND)THEN
   AO%natoms = J + AOTMP%natoms
ENDIF
CALL MEM_ALLOC(AO%ATOMICnORB,AO%natoms)
CALL MEM_ALLOC(AO%ATOMICnBATCH,AO%natoms)
J=0
offsetAT = 0
IF(EXTEND)THEN
   DO I=1,AOTMP%natoms   
      J=J+1
      AO%ATOMICnORB(J)=AOTMP%ATOMICnORB(J)
      AO%ATOMICnBATCH(J)=AOTMP%ATOMICnBATCH(J)
   ENDDO
   offsetAT = AOTMP%natoms
ENDIF
DO I=1,MOLECULE%natoms   
   IF(MOLECULE%ATOM(I)%pointcharge)CYCLE 
   J=J+1
   AO%ATOMICnORB(J)=0
   AO%ATOMICnBATCH(J)=0
ENDDO
AO%empty=.FALSE.
NOFAMILY=SCHEME%NOFAMILY
NOSEGMENT=SCHEME%NOSEGMENT

call mem_alloc(MODELTYPES,BASISINFO%natomtypes)
SUM=0
L=0
IF(UNCONTRACTED)THEN
   DO I=1,BASISINFO%natomtypes    
      SUM1=0
      L=L+1
      MODELTYPES(I)=L
      DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
         SUM1=SUM1+BASISINFO%ATOMTYPE(I)%SHELL(K)%nprim
      ENDDO
      SUM=MAX(SUM1,SUM)         
   ENDDO
ELSE !DEFAULT
   DO I=1,BASISINFO%natomtypes
      SUM1=0
      L=L+1
      MODELTYPES(I)=L
      IF(NOSEGMENT)THEN
         SUM1=SUM1+BASISINFO%ATOMTYPE(I)%nAngmom
      ELSE
         DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
            SUM1=SUM1+BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments
         ENDDO
      ENDIF
      SUM=MAX(SUM1,SUM)
   ENDDO
ENDIF
IF(SUM .EQ. 0) THEN
   WRITE(lupri,*)'BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file'
   CALL FLUSH(LUPRI)
   print*,'BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file'
   CALL LSQUIT('BUILD_AO CALLED BUT NO BATCHES - SOMETHING WRONG&
        & This may be because you want to use densityfitting but&
        & have not specified a auxillary basis set in the molecule input file',lupri)
ENDIF
!SUM IS NOW THE MAXIMUM NUMBER OF BATCHES A SINGLE ATOM CAN HAVE 
AOmodelbat=L
GHOSTFUNCS=0 !SHOULD BE CHANGED TO ACCOUNT FOR GHOST FUNCTIONS
aobatches=(MOLECULE%natoms-GHOSTFUNCS)*SUM
IF(EXTEND)THEN
   aobatches=(MOLECULE%natoms-GHOSTFUNCS)*SUM + AOTMP%nbatches
ENDIF
CALL MEM_ALLOC(AO%BATCH,aobatches)
call nullifyAOBATCH(AO)
IF(IPRINT .GT. 15)WRITE(lupri,*)aobatches,' aobatches should be more than sufficient'
call init_AOmodel(AOmodel,AOmodelbat,SUM)
IF(SUM .EQ. 0) CALL LSQUIT('SUM EQ ZERO SOMETHINGS WRONG',lupri)
nMODELEXP=SUM
call mem_alloc(MODELEXP,SUM)

NULLIFY(AOorganize%ORG)
ALLOCATE(AOorganize%ORG(aobatches))
DO I=1,aobatches
   AOorganize%ORG(I)%exponentindex=0
   DO J=1,maxAOangmom
      AOorganize%ORG(I)%CC(J)=0
   ENDDO
   AOorganize%ORG(I)%atom=0
ENDDO
IF(EXTEND)THEN
   CALL MEM_ALLOC(AO%CC,AOTMP%nCC+SUM*AOmodelbat)
   CALL MEM_ALLOC(AO%Exponents,AOTMP%nEXP+SUM*AOmodelbat)
   CALL MEM_ALLOC(AO%angmom,AOTMP%nCC+SUM*AOmodelbat)
   !COPY PARTIAL
   nbatches=0
   DO I=1,SIZE(AOTMP%BATCH)
      CALL COPY_AOBATCH2(AOTMP%BATCH(I),AO%BATCH(I))
   ENDDO
   DO I=1,AOTMP%nCC
      nrow = AOTMP%CC(I)%nrow
      ncol = AOTMP%CC(I)%ncol
      CALL lsmat_dense_init(AO%CC(I),AOTMP%CC(I)%nrow,AOTMP%CC(I)%ncol)
      call dcopy(nrow*ncol,AOTMP%CC(I)%elms,1,AO%CC(I)%elms,1)
      AO%angmom(I) = AOTMP%angmom(I)
   ENDDO
   DO I=1,AOTMP%nEXP
      nrow = AOTMP%Exponents(I)%nrow
      CALL lsmat_dense_init(AO%Exponents(I),AOTMP%Exponents(I)%nrow,1)
      call dcopy(nrow,AOTMP%Exponents(I)%elms,1,AO%Exponents(I)%elms,1)
   ENDDO
   DO I=1,AOTMP%nbatches
      do J = 1, AOTMP%BATCH(I)%nAngmom
         nullify(AO%BATCH(I)%pCC(J)%p)
         AO%BATCH(I)%pCC(J)%p => AO%CC(AO%BATCH(I)%CCindex(J))
      enddo
      nullify(AO%BATCH(I)%pExponents)
      AO%BATCH(I)%pExponents => AO%Exponents(AO%BATCH(I)%ExpIndex)
   ENDDO
   UNIQUEEXPONENTS=AOTMP%nEXP
   UNIQUEMATRIXES=AOTMP%nCC
   ITYPE = AOTMP%ntype
   AO%nbatches = AOTMP%nbatches
   orbitalIndex = AOTMP%nbast + 1
   primorbitalIndex = AOTMP%nprimbast + 1
ELSE
   CALL MEM_ALLOC(AO%CC,SUM*AOmodelbat)
   CALL MEM_ALLOC(AO%Exponents,SUM*AOmodelbat)
   CALL MEM_ALLOC(AO%angmom,SUM*AOmodelbat)
   nbatches=0
   UNIQUEMATRIXES=0
   UNIQUEEXPONENTS=0
   ITYPE = 0 
   AO%nbatches = 0 
   orbitalIndex = 1
   primorbitalIndex = 1
ENDIF
AOorganize%nbatches =  0 
R = BASISINFO%Labelindex
DO I=1,MOLECULE%natoms  
   IF(MOLECULE%ATOM(I)%pointcharge)CYCLE !no basis functions on this point charge
   IF(IPRINT .GT. 2) THEN
      WRITE(LUPRI,'(2X,A6,I4,2X,A8,F8.5)')'Atom: ',I,'Charge: ',&
           &MOLECULE%ATOM(I)%Charge
   ENDIF
   IF(R.EQ. 0)THEN
      ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
      type = BASISINFO%CHARGEINDEX(ICHARGE)
   ELSE
      type = MOLECULE%ATOM(I)%IDtype(R)
   ENDIF
   IF(IPRINT .GT. 2) THEN
      WRITE(LUPRI,'(2X,A10,A20,2X,A6,I4)')'Basisset: ',BASISINFO%ATOMTYPE(type)%NAME(1:20),'Type: ',type
      WRITE(LUPRI,'(2X,A26)')'-------------------------------------'
   ENDIF
   
   !TEST IF THIS TYPE OF ATOM HAS ALREADY BEEN PROCESSED -> OLDATOM=TRUE
   OLDATOM=.FALSE.
   IF(MODELTYPES(type) .NE. 0)THEN
      IF(AOmodel(MODELTYPES(type))%nbatches .NE. 0) OLDATOM=.TRUE.
   ENDIF
   IF(BASISINFO%ATOMTYPE(type)%nAngmom.EQ.0) CYCLE
   IF(OLDATOM) THEN
      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'AS THIS IS AN OLD ATOM WE COPY FROM MODELAOBATCH'
      L=MODELTYPES(type)
      CALL COPY_FROM_MODEL_AO(AOmodel(L),AO,I,MOLECULE,orbitalindex,primorbitalindex,lupri,offsetAT)
   ELSE
      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'THIS IS A NEW ATOM SO WE BUILD A MODELAO BATCH'
      L=MODELTYPES(type)

      AOmodel(L)%nbatches = 0
      TMPorbitalIndex = 0
      TMPprimorbitalIndex = 0
      
      nAngmom=BASISINFO%ATOMTYPE(type)%nAngmom
      IF(NOSEGMENT)THEN
         nsegments=nAngmom
      ELSE
         nsegments=0
         DO B=1,nAngmom
            nsegments=nsegments+BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
         ENDDO
      ENDIF

      DO B=1,nAngmom
         IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A,I3,A,I3)')'Angmom:',B,&
              &' of ',nAngmom
         IF(BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments .EQ. 0)THEN
            !NON OF THESE ORBITALS - DO NOT ADD BATCH            
         ELSE
            IF(NOSEGMENT)THEN
             IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A)')'This is a new segment so we &
                  &increase the number of batches'
             J=1
             CALL ADD_BATCH(SCHEME,MOLECULE,BASISINFO,AOmodel(L),AOorganize,I,&
                  &type,B,J,lupri,iprint,nPrimitives,nContracted,&
                  &nbatches,UNCONTRACTED,itype)
             nPrimitives = BASISINFO%ATOMTYPE(type)%SHELL(B)%nprim
             nContracted = BASISINFO%ATOMTYPE(type)%SHELL(B)%norb
             CALL ADD_SINGLESEGMENT_AND_CC(AO,AOmodel(L),AOorganize,BASISINFO,nbatches,nPrimitives,&
                  &nContracted,type,B,lupri,iprint,UNIQUEMATRIXES,UNIQUEEXPONENTS,UNCONTRACTED,&
                  &INTNRM,MODELEXP,nMODELEXP,NORM)
             CALL DIRECT_POINTERS(AO,AOmodel(L),AOorganize,nContracted,nPrimitives,nbatches,B,&
                  &lupri,iprint,TMPorbitalIndex,TMPprimorbitalIndex,UNCONTRACTED,INTNRM)
            ELSE
             DO J=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
               IF(IPRINT .GT. 2)  WRITE(LUPRI,'(2X,A,I3,A,I3)')&
                    &'Segment:',J,' of ',&
                    &BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
               
               !TEST IF EXPONENTS is already in AO%Exponents(:) 
               !OR PUT IN ANOTHER WAY 
               !TEST IF IT SHARES EXPONENTS WITH SOME OTHER SEGMENT
               IF(AOmodel(L)%nExp .EQ. 0 )THEN
                  UNIQUE2=0
               ELSE
                  IF(BASISINFO%ATOMTYPE(type)%FAMILY)THEN
                     CALL SHARES_EXPONENTS(LUPRI,IPRINT,AOmodel(L),AOmodel(L)%nExp,BASISINFO,type,B,J,&
                          &UNIQUE2,UNCONTRACTED,MODELEXP,nMODELEXP)
                  ELSE
                     UNIQUE2 = 0
                  ENDIF
               ENDIF
               IF(UNIQUE2 .EQ. 0 .OR. NOFAMILY) THEN
                  !***********************************************************************
                  !*
                  !*   NEW SEGMENT SO NEW BATCH IS ADDED 
                  !*
                  !***********************************************************************
                  
                  IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A)')'This is a new segment so we &
                       &increase the number of batches'
                  
                  CALL ADD_BATCH(SCHEME,MOLECULE,BASISINFO,AOmodel(L),AOorganize,I,&
                       &type,B,J,lupri,iprint,nPrimitives,nContracted,&
                       &nbatches,UNCONTRACTED,itype)
                  CALL ADD_SEGMENT_AND_CC(AO,AOmodel(L),AOorganize,BASISINFO,nbatches,nPrimitives,&
                       &nContracted,type,B,J,lupri,iprint,UNIQUEMATRIXES,UNIQUEEXPONENTS,UNCONTRACTED,&
                       &INTNRM,MODELEXP,nMODELEXP,NORM)
                  
                  CALL DIRECT_POINTERS(AO,AOmodel(L),AOorganize,nContracted,nPrimitives,nbatches,B,&
                       &lupri,iprint,TMPorbitalIndex,TMPprimorbitalIndex,UNCONTRACTED,INTNRM)
                  
               ELSE
                  !***********************************************************************
                  !*
                  !*   THIS SEGMENT SHARES EXPONENTS WITH AN EARLY SEGMENT
                  !*
                  !*   SO A.FAMILY BASISSETIS IS EMPLOYED
                  !*
                  !***********************************************************************
                  IF(IPRINT .GT. 2) WRITE(LUPRI,'(2X,A,I3)')&
                       &'SHARES EXPONENTS WITH ',UNIQUE2
                  CALL DETERMINE_nbathces(AOorganize,I,UNIQUE2,nbat)
                  CALL EXPAND_BATCH(AOmodel(L),AOorganize,basisinfo,type,B,J,nbat,&
                       &A,lupri,UNCONTRACTED)
                  CALL ADD_CC(AO,AOorganize,nbat,nPrimitives,nContracted,lupri,iprint,&
                       &BASISINFO,type,B,J,A,UNIQUEMATRIXES,UNCONTRACTED,INTNRM,NORM)
                  
                  CALL DIRECT_CC_POINTERS(AO,AOmodel(L),AOorganize,nContracted,&
                       &nPrimitives,A,nbat,B,lupri,iprint,TMPorbitalIndex,TMPprimorbitalindex,UNCONTRACTED,intnrm)
               ENDIF
             ENDDO
            ENDIF
         ENDIF
      ENDDO      
      CALL SET_AO_EXTENT(AOmodel(L),SCHEME)
      AOmodel(L)%ATOMICnBATCH(1)=AOmodel(L)%nbatches
      IF(IPRINT .GT. 40)THEN
         WRITE(LUPRI,'(2X,A)')'The Model AO for this atom of this type and basisset'
         CALL PRINT_AO(LUPRI,AOmodel(L))
      ENDIF
      CALL COPY_FROM_MODEL_AO(AOmodel(L),AO,I,MOLECULE,orbitalindex,primorbitalindex,lupri,offsetAT)
   ENDIF
ENDDO

if (.not.ASSOCIATED(AOorganize%ORG)) then
  print*,'memory previously released!!'
  STOP 'Error in BUILD_AO - memory previously released'
endif

DO I =1,AOmodelbat
   call free_AOitem(lupri,AOmodel(I))
ENDDO
DEALLOCATE(AOmodel)
DEALLOCATE(AOorganize%ORG)
NULLIFY(AOorganize%ORG)
call mem_dealloc(MODELTYPES)
call mem_dealloc(MODELEXP)

!AO%nbatches=nbatches
AO%ntype = itype
AO%nCC=UNIQUEMATRIXES
AO%nExp=UNIQUEEXPONENTS
AO%nbast = orbitalindex-1
AO%nprimbast = primorbitalIndex-1
call set_redtype(AO,lupri)
call set_maxJ(AO,lupri)
!CALL DETERMINE_AOBATCH_MEM(AO)
IF(IPRINT .GT. 20) THEN
   WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH.  UNCONTRACTED',UNCONTRACTED
   IF(EXTEND)WRITE(LUPRI,*)'EXTENDED From previous AO '
   CALL PRINT_AO(LUPRI,AO)
   IF(EXTEND)THEN
      WRITE(LUPRI,*)'The Original AO '
      CALL PRINT_AO(LUPRI,AOTMP)      
   ENDIF
ENDIF
IF(EXTEND)THEN
   call free_aoitem(lupri,AOTMP)
ENDIF
END SUBROUTINE BUILD_AO

subroutine init_AOmodel(AOmodel,AOmodelbat,SUM)
implicit none
integer,intent(in) :: AOmodelbat,SUM
TYPE(AOITEM),pointer      :: AOmodel(:)
!
integer :: I

ALLOCATE(AOmodel(AOmodelbat))
DO I=1,AOmodelbat
   call nullifyAOITEM(AOmodel(I))
   AOmodel(I)%nbast = 0
   AOmodel(I)%nprimbast = 0
   AOmodel(I)%natoms=1
   CALL MEM_ALLOC(AOmodel(I)%ATOMICnORB,1)
   AOmodel(I)%ATOMICnORB(1)=1
   CALL MEM_ALLOC(AOmodel(I)%ATOMICnBATCH,1)
   AOmodel(I)%ATOMICnBATCH(1)=0
   CALL MEM_ALLOC(AOmodel(I)%BATCH,SUM)
   call nullifyAOBATCH(AOmodel(I))
   AOmodel(I)%nbatches = 0
   CALL MEM_ALLOC(AOmodel(I)%CC,1) !not used
   CALL lsmat_dense_init(AOmodel(I)%CC(1),1,1)
   AOmodel(I)%CC(1)%elms=0.0E0_realk
   CALL MEM_ALLOC(AOmodel(I)%angmom,1) !not used
   AOmodel(I)%angmom(1) = 0
   AOmodel(I)%nCC = 1 !not used
   AOmodel(I)%nExp = 0   
   CALL MEM_ALLOC(AOmodel(I)%Exponents,SUM)
ENDDO
end subroutine init_AOmodel

subroutine set_redtype(AO,lupri)
implicit none
TYPE(AOITEM)        :: AO
integer :: lupri
!
integer,pointer :: typetoAO(:),typetoredtype(:),redtypetoAO(:)
integer :: itype,I,J,iredtype,nredtype,iang
logical :: newredtype,tmp

call mem_alloc(typetoAO,AO%ntype)
TYPELOOP: do itype = 1,AO%ntype 
   BATCHLOOP:   DO I = 1,AO%nbatches
      IF(AO%BATCH(I)%ITYPE.EQ.itype)then
         typetoAO(itype) = I
         EXIT BATCHLOOP
      ENDIF
   ENDDO BATCHLOOP
enddo TYPELOOP

call mem_alloc(typetoredtype,AO%ntype)
call mem_alloc(redtypetoAO,AO%ntype)
typetoredtype(1) = 1
redtypetoAO(1) = typetoAO(1)
nredtype=1

do itype = 2,AO%ntype 
   I = typetoAO(itype)
   newredtype = .TRUE.
   DO iredtype = 1,nredtype
      !determine if this type is a new redtype or old
      J = redtypetoAO(iredtype)
      tmp = .true.
      IF(AO%BATCH(I)%nAngmom.NE.AO%BATCH(J)%nAngmom) tmp = .false.
      IF(AO%BATCH(I)%nPrimitives.NE.AO%BATCH(J)%nPrimitives) tmp = .false.
      if(tmp)THEN
         DO iang = 1,AO%BATCH(I)%nAngmom
            IF(AO%BATCH(I)%nContracted(iang).NE.AO%BATCH(J)%nContracted(iang)) tmp = .false.
            IF(AO%BATCH(I)%angmom(iang).NE.AO%BATCH(J)%angmom(iang)) tmp = .false.
         ENDDO
      ENDIF
      IF(tmp)THEN
!         WRITE(lupri,*)'This is not a new type  I=',I,'itype',itype
!         call PRINT_AOBATCH(AO%BATCH(I),lupri)         
!         WRITE(lupri,*)'It has the same type as J=',J,'iretypde',iredtype
!         call PRINT_AOBATCH(AO%BATCH(J),lupri)         
         typetoredtype(itype) = iredtype
         newredtype = .FALSE.
         EXIT
      ENDIF
   ENDDO
   IF(newredtype)THEN
!      WRITE(lupri,*)'This is a new type  I=',I
!      call PRINT_AOBATCH(AO%BATCH(I),lupri)         
      nredtype=nredtype+1
      typetoredtype(itype) = nredtype
      redtypetoAO(nredtype) = I
   ENDIF
enddo
AO%nredtype = nredtype
DO I = 1,AO%nbatches
   AO%BATCH(I)%redtype = typetoredtype(AO%BATCH(I)%itype)
ENDDO
call mem_dealloc(typetoAO)
call mem_dealloc(typetoredtype)
call mem_dealloc(redtypetoAO)
end subroutine set_redtype


!!$SUBROUTINE DETERMINE_AOBATCH_MEM(AO)
!!$IMPLICIT NONE
!!$TYPE(AOITEM)              :: AO
!!$INTEGER(KIND=long)        :: alloc_memory
!!$INTEGER                    :: I
!!$alloc_memory = 0
!!$alloc_memory = alloc_memory + sizeof(AO%EMPTY) 
!!$alloc_memory = alloc_memory + sizeof(AO%nbatches) 
!!$alloc_memory = alloc_memory + sizeof(AO%nCC) 
!!$alloc_memory = alloc_memory + sizeof(AO%nExp) 
!!$DO I = 1,AO%nbatches
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%type)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%spherical)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%CENTER)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nPrimitives)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%maxContracted)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%maxAngmom)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%pExponents)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nAngmom)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%extent)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%ANGMOM)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nContracted)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%startOrbital)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%startprimOrbital)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nOrbComp)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nPrimOrbComp)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%nOrbitals)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%pCC)
!!$   alloc_memory = alloc_memory + sizeof(AO%BATCH(I)%CCindex)
!!$ENDDO
!!$DO I = 1,AO%nCC
!!$   alloc_memory = alloc_memory+sizeof(AO%CC(I)%elms)+sizeof(AO%CC(I)%nrow)+sizeof(AO%CC(I)%ncol)+sizeof(AO%CC(I)%complex)
!!$ENDDO
!!$DO I = 1,AO%nEXP
!!$   alloc_memory = alloc_memory+sizeof(AO%Exponents(I)%elms)+sizeof(AO%Exponents(I)%nrow)+sizeof(AO%Exponents(I)%ncol)+sizeof(AO%Exponents(I)%complex)
!!$ENDDO
!!$END SUBROUTINE DETERMINE_AOBATCH_MEM

!> \brief builds an empty AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> build an empty AOitem used, for when calculating non 4center integrals
!>
SUBROUTINE BUILD_EMPTY_AO(AO,LUPRI)
IMPLICIT NONE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the logical unit number for the output file
INTEGER                   :: LUPRI
AO%natoms = 1                       
CALL MEM_ALLOC(AO%ATOMICnORB,AO%natoms) 
CALL MEM_ALLOC(AO%ATOMICnBATCH,AO%natoms)
AO%ATOMICnORB(1)=1                  
AO%ATOMICnBATCH(1)=1                
AO%nbast = 1                
AO%nprimbast = 1

AO%empty=.TRUE.
AO%nbatches=1
AO%nCC=1
AO%nExp=1
CALL MEM_ALLOC(AO%BATCH,1)   
call nullifyAOBATCH(AO)
CALL MEM_ALLOC(AO%CC,1)
CALL MEM_ALLOC(AO%angmom,1)
CALL MEM_ALLOC(AO%Exponents,1)
CALL lsmat_dense_init(AO%CC(1),1,1)
CALL lsmat_dense_init(AO%Exponents(1),1,1)
AO%angmom(1)=0
AO%CC(1)%elms(1)=1
AO%Exponents(1)%elms(1)=0
AO%BATCH(1)%itype = 1
AO%BATCH(1)%redtype = 1
AO%BATCH(1)%type_Empty = .TRUE.
AO%BATCH(1)%type_Nucleus = .FALSE.
AO%BATCH(1)%type_elField = .FALSE.
AO%BATCH(1)%type_pCharge = .FALSE.
AO%BATCH(1)%spherical=.false.
AO%BATCH(1)%atom=1
AO%BATCH(1)%molecularIndex=1
AO%BATCH(1)%batch=1
AO%BATCH(1)%CENTER(1)=0E0_realk
AO%BATCH(1)%CENTER(2)=0E0_realk
AO%BATCH(1)%CENTER(3)=0E0_realk
AO%BATCH(1)%nPrimitives=1
AO%BATCH(1)%maxContracted=1
AO%BATCH(1)%maxAngmom=0
AO%BATCH(1)%pExponents => AO%Exponents(1)
AO%BATCH(1)%nAngmom=1
AO%BATCH(1)%extent=0E0_realk
AO%BATCH(1)%ANGMOM(1)=0
AO%BATCH(1)%nContracted(1)=1
AO%BATCH(1)%startOrbital(1)=1
AO%BATCH(1)%startprimOrbital(1)=1
AO%BATCH(1)%nOrbComp(1)=1
AO%BATCH(1)%nPrimOrbComp(1)=1
AO%BATCH(1)%nOrbitals(1)=1
AO%BATCH(1)%pCC(1)%p => AO%CC(1)
AO%BATCH(1)%CCindex(1) = 0
AO%BATCH(1)%Expindex = 0

!WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH'
!CALL PRINT_AO(LUPRI,AO)
!CALL DETERMINE_AOBATCH_MEM(AO)
AO%ntype = 1
AO%nredtype = 1
AO%maxJ = 0
END SUBROUTINE BUILD_EMPTY_AO

!> \brief builds an empty nuclear AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> build an empty nuclear AOitem, used for nuclear attraction integrals
!>
SUBROUTINE BUILD_EMPTY_NUCLEAR_AO(AO,MOLECULE,LUPRI)
use molecule_type
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!
INTEGER                   :: aobatches,I,ncharges,J,NEWCHARGE,nAtoms,K
REAL(REALK),pointer       :: CHARGES(:)

J=0
DO I=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(I)%phantom)CYCLE !no actul charge
   J = J +1
ENDDO
nAtoms = J

AO%natoms = 1                            
CALL MEM_ALLOC(AO%ATOMICnORB,AO%natoms)       
CALL MEM_ALLOC(AO%ATOMICnBATCH,AO%natoms)     
AO%ATOMICnORB(1)=0                       
AO%ATOMICnBATCH(1)=1                     
AO%nbast = 1                     
AO%nprimbast = 1                     

AO%empty=.TRUE.
aobatches=nAtoms    !MOLECULE%natoms
AO%nbatches=aobatches
!AO%nCC=1
AO%nExp=1

CALL MEM_ALLOC(AO%BATCH,aobatches)
call nullifyAOBATCH(AO)
CALL MEM_ALLOC(AO%Exponents,1)
CALL lsmat_dense_init(AO%Exponents(1),1,1)
call mem_alloc(CHARGES,MOLECULE%natoms)
ncharges=0
I=0
DO K=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(K)%phantom)CYCLE !no actul charge
   I=I+1
   NEWCHARGE=0
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(K)%Charge-charges(J) ).LE. 1E-10_realk ) NEWCHARGE=J
   ENDDO
   IF(NEWCHARGE .EQ. 0)THEN
      ncharges=ncharges+1
      CHARGES(ncharges)=MOLECULE%ATOM(K)%Charge
   ENDIF
ENDDO

AO%nCC=nCharges

CALL MEM_ALLOC(AO%CC,nCharges)
CALL MEM_ALLOC(AO%angmom,nCharges)

DO I=1,nCharges
   CALL lsmat_dense_init(AO%CC(I),1,1)
   AO%CC(I)%elms(1)=-CHARGES(I)
   AO%angmom(I) = 0
ENDDO

I=0
DO K=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(K)%phantom)CYCLE !no actul charge
   I = I + 1
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(K)%Charge-charges(J)).LE. 1E-10_realk) NEWCHARGE=J
   ENDDO
   AO%BATCH(I)%pCC(1)%p => AO%CC(NEWCHARGE)   
   AO%BATCH(I)%CCindex(1) = NEWCHARGE   
ENDDO
AO%BATCH(I)%Expindex = 1

call mem_dealloc(CHARGES)

AO%Exponents(1)%elms(1)=0E0_realk

I=0
DO J=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(J)%phantom)CYCLE !no actul charge
   I=I+1
   AO%BATCH(I)%itype = I
   AO%BATCH(I)%redtype = 1
   AO%BATCH(I)%type_Nucleus = .TRUE.
   AO%BATCH(I)%type_elField = .FALSE.
   AO%BATCH(I)%type_pCharge = .FALSE.
   AO%BATCH(I)%TYPE_Empty = .FALSE.
   AO%BATCH(I)%spherical=.false.
   AO%BATCH(I)%atom=J
   AO%BATCH(I)%molecularIndex=MOLECULE%ATOM(J)%molecularIndex
   AO%BATCH(I)%batch=1
   AO%BATCH(I)%CENTER(1)=MOLECULE%ATOM(J)%CENTER(1)
   AO%BATCH(I)%CENTER(2)=MOLECULE%ATOM(J)%CENTER(2)
   AO%BATCH(I)%CENTER(3)=MOLECULE%ATOM(J)%CENTER(3)
   AO%BATCH(I)%nPrimitives=1
   AO%BATCH(I)%maxContracted=1
   AO%BATCH(I)%maxAngmom=0
   AO%BATCH(I)%pExponents => AO%Exponents(1)
   AO%BATCH(I)%nAngmom=1
   AO%BATCH(I)%extent=0E0_realk
   AO%BATCH(I)%ANGMOM(1)=0
   AO%BATCH(I)%nContracted(1)=1
   AO%BATCH(I)%startOrbital(1)=1
   AO%BATCH(I)%startprimOrbital(1)=I
   AO%BATCH(I)%nOrbComp(1)=1
   AO%BATCH(I)%nPrimOrbComp(1)=1
   AO%BATCH(I)%nOrbitals(1)=1
ENDDO
AO%ntype = natoms!MOLECULE%natoms
AO%nredtype = 1
AO%maxJ = 0
!WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH'
!CALL PRINT_AO(LUPRI,AO)
!CALL DETERMINE_AOBATCH_MEM(AO)

END SUBROUTINE BUILD_EMPTY_NUCLEAR_AO

!> \brief builds an empty single nuclear AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> build an empty nuclear AOitem, used for nuclear attraction integrals
!>
SUBROUTINE BUILD_EMPTY_SINGLE_NUCLEAR_AO(AO,MOLECULE,LUPRI,IATOM)
use molecule_type
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> which atom to choose
INTEGER                   :: IATOM
!
INTEGER                   :: aobatches,I,ncharges,J,NEWCHARGE,nAtoms,K
REAL(REALK),pointer       :: CHARGES(:)

IF(MOLECULE%ATOM(IATOM)%phantom)THEN
   CALL LSQUIT('BUILD_EMPTY_SINGLE_NUCLEAR_AO: Error Phantom atom',-1)
ENDIF
nAtoms = 1
AO%natoms = 1                            
CALL MEM_ALLOC(AO%ATOMICnORB,AO%natoms)       
CALL MEM_ALLOC(AO%ATOMICnBATCH,AO%natoms)     
AO%ATOMICnORB(1)=0                       
AO%ATOMICnBATCH(1)=1                     
AO%nbast = 1                     
AO%nprimbast = 1                     

AO%empty=.TRUE.
aobatches=nAtoms  ! = 1
AO%nbatches=aobatches
!AO%nCC=1
AO%nExp=1

CALL MEM_ALLOC(AO%BATCH,aobatches)
call nullifyAOBATCH(AO)
CALL MEM_ALLOC(AO%Exponents,1)
CALL lsmat_dense_init(AO%Exponents(1),1,1)
AO%nCC=1
CALL MEM_ALLOC(AO%CC,1)
CALL MEM_ALLOC(AO%angmom,1)
CALL lsmat_dense_init(AO%CC(1),1,1)
AO%CC(1)%elms(1)=-MOLECULE%ATOM(IATOM)%Charge
AO%angmom(1) = 0

AO%BATCH(1)%pCC(1)%p => AO%CC(1)
AO%BATCH(1)%CCindex(1) = 1
AO%Exponents(1)%elms(1)=0E0_realk
AO%BATCH(1)%Expindex = 1

I=1
AO%BATCH(I)%itype = I
AO%BATCH(I)%redtype = 1
AO%BATCH(I)%type_Nucleus = .TRUE.
AO%BATCH(I)%type_elField = .FALSE.
AO%BATCH(I)%type_pCharge = .FALSE.
AO%BATCH(I)%TYPE_Empty = .FALSE.
AO%BATCH(I)%spherical=.false.
AO%BATCH(I)%atom=IATOM
AO%BATCH(I)%molecularIndex=MOLECULE%ATOM(IATOM)%molecularIndex
AO%BATCH(I)%batch=1
AO%BATCH(I)%CENTER(1)=MOLECULE%ATOM(IATOM)%CENTER(1)
AO%BATCH(I)%CENTER(2)=MOLECULE%ATOM(IATOM)%CENTER(2)
AO%BATCH(I)%CENTER(3)=MOLECULE%ATOM(IATOM)%CENTER(3)
AO%BATCH(I)%nPrimitives=1
AO%BATCH(I)%maxContracted=1
AO%BATCH(I)%maxAngmom=0
AO%BATCH(I)%pExponents => AO%Exponents(1)
AO%BATCH(I)%nAngmom=1
AO%BATCH(I)%extent=0E0_realk
AO%BATCH(I)%ANGMOM(1)=0
AO%BATCH(I)%nContracted(1)=1
AO%BATCH(I)%startOrbital(1)=1
AO%BATCH(I)%startprimOrbital(1)=I
AO%BATCH(I)%nOrbComp(1)=1
AO%BATCH(I)%nPrimOrbComp(1)=1
AO%BATCH(I)%nOrbitals(1)=1
AO%ntype = 1
AO%nredtype = 1
AO%maxJ = 0

END SUBROUTINE BUILD_EMPTY_SINGLE_NUCLEAR_AO

!> \brief builds an empty electric-field AOitem
!> \author S. Reine
!> \date 2014
!>
!> build an empty electric-field  AOitem, used for electric-field integrals
!>
SUBROUTINE BUILD_EMPTY_ELFIELD_AO(AO,MOLECULE,LUPRI)
use molecule_type
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!
INTEGER                   :: aobatches,I,ncharges,J,NEWCHARGE,nAtoms,K
REAL(REALK),pointer       :: CHARGES(:)

J=0
DO I=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(I)%phantom)CYCLE !no actul charge
   J = J +1
ENDDO
nAtoms = J

AO%natoms = 1                            
CALL MEM_ALLOC(AO%ATOMICnORB,AO%natoms)       
CALL MEM_ALLOC(AO%ATOMICnBATCH,AO%natoms)     
AO%ATOMICnORB(1)=0                       
AO%ATOMICnBATCH(1)=1                     
AO%nbast = 1
AO%nprimbast = 1

AO%empty=.TRUE.
aobatches=nAtoms    !MOLECULE%natoms
AO%nbatches=aobatches
!AO%nCC=1
AO%nExp=1

CALL MEM_ALLOC(AO%BATCH,aobatches)
call nullifyAOBATCH(AO)
CALL MEM_ALLOC(AO%Exponents,1)
CALL lsmat_dense_init(AO%Exponents(1),1,1)
call mem_alloc(CHARGES,MOLECULE%natoms)
ncharges=0
I=0
DO K=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(K)%phantom)CYCLE !no actul charge
   I=I+1
   NEWCHARGE=0
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(K)%Charge-charges(J) ).LE. 1E-10_realk ) NEWCHARGE=J
   ENDDO
   IF(NEWCHARGE .EQ. 0)THEN
      ncharges=ncharges+1
      CHARGES(ncharges)=MOLECULE%ATOM(K)%Charge
   ENDIF
ENDDO

AO%nCC=nCharges

CALL MEM_ALLOC(AO%CC,nCharges)
CALL MEM_ALLOC(AO%angmom,nCharges)

DO I=1,nCharges
   CALL lsmat_dense_init(AO%CC(I),1,1)
   AO%CC(I)%elms(1)=-CHARGES(I)
   AO%angmom(I) = 0
ENDDO

I=0
DO K=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(K)%phantom)CYCLE !no actul charge
   I = I + 1
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(K)%Charge-charges(J)).LE. 1E-10_realk) NEWCHARGE=J
   ENDDO
   AO%BATCH(I)%pCC(1)%p => AO%CC(NEWCHARGE)   
   AO%BATCH(I)%CCindex(1) = NEWCHARGE   
   AO%BATCH(I)%Expindex = 1
ENDDO

call mem_dealloc(CHARGES)

AO%Exponents(1)%elms(1)=0E0_realk

I=0
DO J=1,MOLECULE%natoms
   IF(MOLECULE%ATOM(J)%phantom)CYCLE !no actul charge
   I=I+1
   AO%BATCH(I)%itype = I
   AO%BATCH(I)%redtype = 1
   AO%BATCH(I)%type_Nucleus = .FALSE.
   AO%BATCH(I)%type_elField = .TRUE.
   AO%BATCH(I)%type_pCharge = .FALSE.
   AO%BATCH(I)%TYPE_Empty = .FALSE.
   AO%BATCH(I)%spherical=.false.
   AO%BATCH(I)%atom=J
   AO%BATCH(I)%molecularIndex=MOLECULE%ATOM(J)%molecularIndex
   AO%BATCH(I)%batch=1
   AO%BATCH(I)%CENTER(1)=MOLECULE%ATOM(J)%CENTER(1)
   AO%BATCH(I)%CENTER(2)=MOLECULE%ATOM(J)%CENTER(2)
   AO%BATCH(I)%CENTER(3)=MOLECULE%ATOM(J)%CENTER(3)
   AO%BATCH(I)%nPrimitives=1
   AO%BATCH(I)%maxContracted=1
   AO%BATCH(I)%maxAngmom=0
   AO%BATCH(I)%pExponents => AO%Exponents(1)
   AO%BATCH(I)%nAngmom=1
   AO%BATCH(I)%extent=0E0_realk
   AO%BATCH(I)%ANGMOM(1)=0
   AO%BATCH(I)%nContracted(1)=1
   AO%BATCH(I)%startOrbital(1)=1
   AO%BATCH(I)%startprimOrbital(1)=I
   AO%BATCH(I)%nOrbComp(1)=1
   AO%BATCH(I)%nPrimOrbComp(1)=1
   AO%BATCH(I)%nOrbitals(1)=1
ENDDO
AO%ntype = natoms!MOLECULE%natoms
AO%nredtype = 1
AO%maxJ = 0
!WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH'
!CALL PRINT_AO(LUPRI,AO)
!CALL DETERMINE_AOBATCH_MEM(AO)

END SUBROUTINE BUILD_EMPTY_ELFIELD_AO

!> \brief builds an empty nuclear AOitem
!> \author T. Kjaergaard
!> \date 2008
!>
!> build an empty nuclear AOitem, used for nuclear attraction integrals
!>
SUBROUTINE BUILD_EMPTY_PCHARGE_AO(AO,MOLECULE,LUPRI)
use molecule_type
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!
INTEGER                   :: aobatches,I,ncharges,J,NEWCHARGE
REAL(REALK),pointer       :: CHARGES(:)

AO%natoms = 1!MOLECULE%NATOMS
CALL MEM_ALLOC(AO%ATOMICnORB,AO%natoms)       
CALL MEM_ALLOC(AO%ATOMICnBATCH,AO%natoms)
DO I=1,AO%natoms     
   AO%ATOMICnORB(I)=0                       
   AO%ATOMICnBATCH(I)=MOLECULE%NATOMS
ENDDO
AO%nbast = MOLECULE%NATOMS    
AO%nprimbast = MOLECULE%NATOMS    

AO%empty=.TRUE.
aobatches=MOLECULE%natoms
AO%nbatches=aobatches
!AO%nCC=1
AO%nExp=1

CALL MEM_ALLOC(AO%BATCH,aobatches)
call nullifyAOBATCH(AO)
CALL MEM_ALLOC(AO%Exponents,1)

CALL lsmat_dense_init(AO%Exponents(1),1,1)

call mem_alloc(CHARGES,MOLECULE%natoms)
CHARGES(1)=MOLECULE%ATOM(1)%Charge
ncharges=1
DO I=2,MOLECULE%natoms
   NEWCHARGE=0
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(I)%Charge-charges(J) ).LE.1.d-10 ) NEWCHARGE=J
   ENDDO
   IF(NEWCHARGE .EQ. 0)THEN
      ncharges=ncharges+1
      CHARGES(ncharges)=MOLECULE%ATOM(I)%Charge
   ENDIF
ENDDO

AO%nCC=nCharges

CALL MEM_ALLOC(AO%CC,nCharges)
CALL MEM_ALLOC(AO%angmom,nCharges)
DO I=1,nCharges
   CALL lsmat_dense_init(AO%CC(I),1,1)
   AO%CC(I)%elms(1)=-CHARGES(I)
   AO%angmom(I) = 0
ENDDO

DO I=1,MOLECULE%natoms
   DO J=1,ncharges
      IF(ABS(MOLECULE%ATOM(I)%Charge-charges(J)).LE.1.d-10) NEWCHARGE=J
   ENDDO
   AO%BATCH(I)%pCC(1)%p => AO%CC(NEWCHARGE)   
   AO%BATCH(I)%CCindex(1) = NEWCHARGE   
   AO%BATCH(I)%Expindex = 1
ENDDO

call mem_dealloc(CHARGES)

AO%Exponents(1)%elms(1)=0.d0

DO I=1,MOLECULE%natoms
   AO%BATCH(I)%itype = I
   AO%BATCH(I)%redtype = 1
   AO%BATCH(I)%type_Nucleus = .TRUE.
   AO%BATCH(I)%type_elField = .FALSE.
   AO%BATCH(I)%type_pCharge = .TRUE.
   AO%BATCH(I)%TYPE_Empty = .FALSE.
   AO%BATCH(I)%spherical=.false.
   AO%BATCH(I)%atom=I
   AO%BATCH(I)%molecularIndex = MOLECULE%ATOM(I)%molecularIndex
   AO%BATCH(I)%batch=1
   AO%BATCH(I)%CENTER(1)=MOLECULE%ATOM(I)%CENTER(1)
   AO%BATCH(I)%CENTER(2)=MOLECULE%ATOM(I)%CENTER(2)
   AO%BATCH(I)%CENTER(3)=MOLECULE%ATOM(I)%CENTER(3)
   AO%BATCH(I)%nPrimitives=1
   AO%BATCH(I)%maxContracted=1
   AO%BATCH(I)%maxAngmom=0
   AO%BATCH(I)%pExponents => AO%Exponents(1)
   AO%BATCH(I)%nAngmom=1
   AO%BATCH(I)%extent=0.0E0_realk
   AO%BATCH(I)%ANGMOM(1)=0
   AO%BATCH(I)%nContracted(1)=1
   AO%BATCH(I)%startOrbital(1)=1
   AO%BATCH(I)%startprimOrbital(1)=I
   AO%BATCH(I)%nOrbComp(1)=1
   AO%BATCH(I)%nPrimOrbComp(1)=1
   AO%BATCH(I)%nOrbitals(1)=1
ENDDO
AO%ntype = MOLECULE%natoms
AO%nredtype = 1
AO%maxJ = 0
!WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL AO BATCH'
!CALL PRINT_AO(LUPRI,AO)
!CALL DETERMINE_AOBATCH_MEM(AO)

END SUBROUTINE BUILD_EMPTY_PCHARGE_AO

SUBROUTINE determinenAObatches(nAObatches,LUPRI,SCHEME,&
              & IPRINT,molecule,BASISINFO,UNCONTRACTED,INTNRM)
use molecule_type
implicit none
!> number of AO batches
INTEGER,intent(inout)    :: nAObatches
!> the logical unit number for the output file
INTEGER,intent(in)           :: LUPRI
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME),intent(in) :: SCHEME
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)       :: IPRINT
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in):: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in) :: BASISINFO
!> if the AOitem should be uncontracted
LOGICAL,intent(in)       :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL,intent(in)       :: INTNRM 
!
TYPE(AOITEM)              :: AOfull
CALL BUILD_AO(LUPRI,SCHEME,IPRINT,MOLECULE,BASISINFO,AOfull&
     &,UNCONTRACTED,INTNRM,.FALSE.)
nAObatches = AOfull%nbatches
call free_aoitem(lupri,AOfull)
end SUBROUTINE determinenAObatches
!> \brief builds an AOitem with a single shellbatch (an S, or PxPyPz, etc.), or several
!> \author T. Kjaergaard
!> \date 2008
!>
!> builds an AOitem with a single Orbitalbatch used for coupled-cluster
!> calculations, where the 4 center 2 electron integrals are calculated 
!> in batches to avoid memory problems 
!>
SUBROUTINE BUILD_SHELLBATCH_AO(LUPRI,SCHEME,IPRINT,MOLECULE,BASISINFO,&
     &AO_output,UNCONTRACTED,INTNRM,RequestedBatchIndex,dim,SizeRequestedBatch)
use molecule_type
implicit none
!> the logical unit number for the output file
INTEGER,intent(in)           :: LUPRI
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME),intent(in) :: SCHEME
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)       :: IPRINT
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in):: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in) :: BASISINFO
!> the AOitem to be build
TYPE(AOITEM),intent(inout) :: AO_output
!> if the AOitem should be uncontracted
LOGICAL,intent(in)       :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL,intent(in)       :: INTNRM 
!> the requested batch index to build the AOitem from
INTEGER,intent(in)       :: RequestedBatchIndex
!> dimension of the orbital batch (1 for s, 3 for p, 5 for d,etc)
INTEGER,intent(inout)    :: dim
!> the requested size of the batch
INTEGER,intent(in)       :: SizeRequestedBatch
!
TYPE(AOITEM)              :: AOfull
INTEGER :: I,A,iexp,nrow,ncol,Ifull,IATOMfullold,IATOM,norbitals,iATOMfull
INTEGER :: nPrimOrbitals,nprim,ntype,unique,itest,itype,batch
INTEGER :: K
logical :: extend2
integer,pointer :: uniquetypes(:)
real(realk)               :: SUM
REAL(REALK),PARAMETER     :: NULL=0.0E0_realk
extend2 = .FALSE.
CALL BUILD_AO(LUPRI,SCHEME,IPRINT,MOLECULE,BASISINFO,AOfull&
     &,UNCONTRACTED,INTNRM,extend2)

IF(AOfull%nbatches.LT.RequestedBatchIndex+SizeRequestedBatch-1)THEN
   CALL LSQUIT('Error in BUILD_SHELLBATCH_AO',-1)
ENDIF
call nullifyAOITEM(AO_output)
AO_output%EMPTY = AOfull%EMPTY

call mem_alloc(AO_output%BATCH,SizeRequestedBatch)
call nullifyAOBATCH(AO_output)

AO_output%nCC = AOfull%nCC
call mem_alloc(AO_output%CC,AO_output%nCC)
call mem_alloc(AO_output%angmom,AO_output%nCC)
do I=1,AO_output%nCC
   nrow = AOfull%CC(I)%nrow
   ncol = AOfull%CC(I)%ncol
   CALL lsmat_dense_init(AO_output%CC(I),nrow,ncol)
   call dcopy(nrow*nCol,AOfull%CC(I)%elms,1,AO_output%CC(I)%elms,1)   
   AO_output%angmom(I) = AOfull%angmom(I)
enddo

AO_output%nEXP = AOfull%nEXP
call mem_alloc(AO_output%Exponents,AO_output%nEXP)
do I=1,AO_output%nEXP
   nrow = AOfull%Exponents(I)%nrow
   CALL lsmat_dense_init(AO_output%Exponents(I),nrow,1)
   call dcopy(nrow,AOfull%Exponents(I)%elms,1,AO_output%Exponents(I)%elms,1)
enddo

IATOMfullold = 0
IATOM=0
DO I=1,SizeRequestedBatch
   Ifull = RequestedBatchIndex+I-1
   IATOMfull = AOfull%BATCH(Ifull)%atom 
   IF(IATOMfull.NE.IATOMfullold)THEN
      IATOM=IATOM+1
   ENDIF
   IATOMfullold = IATOMfull
ENDDO
AO_output%nAtoms = IATOM
call mem_alloc(AO_output%ATOMICnORB,AO_output%nAtoms)
call mem_alloc(AO_output%ATOMICnBATCH,AO_output%nAtoms)
IATOMfullold = 0
IATOM=0
nOrbitals = 0
nPrimOrbitals = 0
AO_output%nbatches = SizeRequestedBatch
BATCH = 0
DO I=1,SizeRequestedBatch
   Ifull = RequestedBatchIndex+I-1
   IATOMfull = AOfull%BATCH(Ifull)%atom 
   IF(IATOMfull.NE.IATOMfullold)THEN
      IATOM=IATOM+1
      IATOMfullold = IATOMfull
      BATCH = 0
   ENDIF
   call copy_aobatch2(AOfull%BATCH(Ifull),AO_output%BATCH(I))
   BATCH = BATCH + 1 
   AO_output%BATCH(I)%batch = BATCH
   AO_output%BATCH(I)%atom = IATOM
   nprim = AO_output%BATCH(I)%nPrimitives
   DO A=1,AO_output%BATCH(I)%nAngmom
      AO_output%BATCH(I)%pCC(A)%p=> AO_output%CC(AO_output%BATCH(I)%CCindex(A))
      AO_output%BATCH(I)%startOrbital(A) =  1 + nOrbitals
      nOrbitals = nOrbitals + AO_output%BATCH(I)%nOrbitals(A)
      AO_output%BATCH(I)%startPrimOrbital(A) = 1 + nPrimOrbitals
      nPrimOrbitals = nPrimOrbitals + nPrim*AO_output%BATCH(I)%nOrbComp(A)
   ENDDO
   AO_output%BATCH(I)%pExponents => AO_output%Exponents(AO_output%BATCH(I)%Expindex)
ENDDO
DO I=1,AO_output%natoms
   AO_output%ATOMICnOrb(I) = 0
   AO_output%ATOMICnBatch(I) = 0
ENDDO
DO I=1,AO_output%nbatches
   IATOM = AO_output%BATCH(I)%atom 
   DO A=1,AO_output%BATCH(I)%nAngmom
      AO_output%ATOMICnOrb(IATOM) = AO_output%ATOMICnOrb(IATOM)+AO_output%BATCH(I)%nOrbitals(A)
   ENDDO
   AO_output%ATOMICnBatch(IATOM) = AO_output%ATOMICnBatch(IATOM)+1
ENDDO


call mem_alloc(uniquetypes,AOfull%ntype)
ntype = 0
DO I=1,AO_output%nbatches
   Ifull = RequestedBatchIndex+I-1
   itest = AOfull%BATCH(Ifull)%itype
   unique = 0
   do itype = 1,ntype
      if(itest.EQ.uniquetypes(itype)) unique = itype !not unique
   enddo
   IF(unique.NE.0)THEN
      !not unique
      AO_output%BATCH(I)%itype = unique
   else
      !unique
      ntype = ntype + 1
      AO_output%BATCH(I)%itype = ntype
      uniquetypes(ntype) = itest
   endif
ENDDO
call mem_dealloc(uniquetypes)
AO_output%ntype = ntype
call set_redtype(AO_output,lupri)
call set_maxJ(AO_output,lupri)
AO_output%nbast = nOrbitals
AO_output%nprimbast = nPrimOrbitals
dim = nOrbitals
DO I=1,AO_output%nbatches
   IATOM = AO_output%BATCH(I)%atom 
   BATCH = AO_output%BATCH(I)%batch
   IF(BATCH.GT.AO_output%ATOMICnBatch(IATOM))THEN
      call lsquit('error in BUILD_SHELLBATCH_AO',-1)
   ENDIF
ENDDO
IF(IPRINT .GT. 20) THEN
   WRITE(LUPRI,*)'BUILD_AOBATCH:  PRINTING FINAL SHELLBATCH AO'
   CALL PRINT_AO(LUPRI,AO_output)
ENDIF

call free_aoitem(lupri,AOfull)
END SUBROUTINE BUILD_SHELLBATCH_AO

subroutine build_batchesOfAOs(lupri,setting,maxallowedorbitals,nbast,&
     &MaxOrbitals,batchsize,batchdim,batchindex,nbatches,orbTobatch,AOspec)
implicit none
integer,intent(in)         :: lupri,maxallowedorbitals,nbast
type(lssetting) :: setting
integer,pointer :: batchsize(:),batchdim(:),batchindex(:)
integer :: orbtoBatch(nbast)
integer,intent(inout) :: nbatches,MaxOrbitals
character(len=1),intent(in) :: AOspec
!
integer :: I,A,norbitals,nbatLoc,iOrb,tmporb,allocnbatches,extend2
logical :: uncont,intnrm,extend
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
uncont=.FALSE.
intnrm = .false.; extend=.false.
extend = .false.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'D')THEN
   !   The Aux AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(AuxBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)
   extend = .true.
ELSEIF(AOspec.EQ.'O')THEN
   ! only The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(CABBasParam)
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
     & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
call build_batchesOfAOs1(lupri,maxallowedorbitals,nbast,&
     &MaxOrbitals,batchsize,batchdim,batchindex,nbatches,orbTobatch,AO,0)
IF(extend)THEN
   extend2 = AO%nbatches
   call free_aoitem(lupri,AO)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,.FALSE.)
   call build_batchesOfAOs1(lupri,maxallowedorbitals,nbast,&
        &MaxOrbitals,batchsize,batchdim,batchindex,nbatches,orbTobatch,AO,extend2)
ENDIF
call free_aoitem(lupri,AO)
end subroutine build_batchesOfAOs

subroutine build_batchesOfAOs1(lupri,maxallowedorbitals,nbast,&
     &MaxOrbitals,batchsize,batchdim,batchindex,nbatches,orbTobatch,AO,extend)
implicit none
integer,intent(in)         :: lupri,maxallowedorbitals,nbast
integer,pointer :: batchsize(:),batchdim(:),batchindex(:)
integer,intent(inout) :: orbtoBatch(nbast)
integer,intent(inout) :: nbatches,MaxOrbitals
type(AOITEM),intent(in) :: AO
integer,intent(in) :: extend
!
integer,pointer :: tmp(:)
integer :: I,A,norbitals,nbatLoc,iOrb,tmporb,allocnbatches,nbatchesOrig
logical :: uncont,intnrm
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2

nbatchesOrig = 0
IF(extend.GT.0)nbatchesOrig = nbatches 
nbatches = 0
norbitals = 0
do I=1,AO%nbatches
   tmporb = 0
   DO A=1,AO%BATCH(I)%nAngmom
      tmporb = tmporb + AO%BATCH(I)%norbitals(A)
   ENDDO
   IF(norbitals + tmporb.GT.maxallowedorbitals)THEN
      nbatches = nbatches + 1
      norbitals = tmporb
      IF(tmporb.GT.maxallowedorbitals)THEN
         CALL LSQUIT('nOrbitals is greater then maxallowedorbitals in typedef_buildbatchesOfAOs',-1)
      ENDIF
   ELSE
      norbitals = norbitals + tmporb
   ENDIF
enddo
IF(norbitals.GT.0)THEN
   nbatches = nbatches + 1
ENDIF

IF(extend.GT.0)THEN
   nbatches = nbatchesOrig+nbatches
   call mem_alloc(tmp,nbatchesOrig)
   do I=1,nbatchesOrig
      tmp(I)=batchsize(I)
   enddo
   call mem_dealloc(batchsize)   
   call mem_alloc(batchsize,nbatches)
   do I=1,nbatchesOrig
      batchsize(I)=tmp(I)
      tmp(I)=batchindex(I)
   enddo
   call mem_dealloc(batchindex)
   call mem_alloc(batchindex,nbatches)
   do I=1,nbatchesOrig
      batchindex(I)=tmp(I)
      tmp(I)=batchdim(I)
   enddo
   call mem_dealloc(batchdim)
   call mem_alloc(batchdim,nbatches)
   do I=1,nbatchesOrig
      batchdim(I)=tmp(I)
   enddo
   call mem_dealloc(tmp)
   do I=nbatchesOrig+1,nbatches
      batchsize(I) = 0
      batchindex(I) = 0
      batchdim(I) = 0
   enddo
ELSE
   call mem_alloc(batchsize,nbatches)
   call mem_alloc(batchindex,nbatches)
   call mem_alloc(batchdim,nbatches)
   do I=1,nbatches
      batchsize(I) = 0
      batchindex(I) = 0
      batchdim(I) = 0
   enddo
ENDIF
allocnbatches = nbatches
nbatLoc = 0
nbatches = nbatchesOrig
norbitals = 0
batchindex(nbatches+1) = extend + 1 
do I=1,AO%nbatches
   tmporb = 0
   DO A=1,AO%BATCH(I)%nAngmom
      tmporb = tmporb + AO%BATCH(I)%norbitals(A)
   ENDDO
   IF(norbitals + tmporb.GT.maxallowedorbitals)THEN
      nbatches = nbatches + 1
      batchsize(nbatches) = nbatLoc
      batchdim(nbatches) = norbitals
      IF(nbatches+1.LE.allocnbatches)THEN
         !extend = AO%nbatches of previous AO that need extending
         batchindex(nbatches+1) = extend + I 
      ENDIF
      nbatLoc = 1
      norbitals = tmporb
   ELSE
      nbatLoc = nbatLoc + 1
      norbitals = norbitals + tmporb
   ENDIF
enddo
IF(norbitals.GT.0)THEN
   nbatches = nbatches + 1
   batchsize(nbatches) = nbatLoc
   batchdim(nbatches) = norbitals
ENDIF

norbitals = 0         
nbatLoc = 0
do I=1,nbatches
   norbitals=norbitals + batchdim(I)
enddo
do I=nbatchesOrig+1,nbatches
   nbatLoc = nbatLoc + batchsize(I) 
enddo
IF(extend.GT.0)THEN
   IF(norbitals.NE.nbast)call lsquit('basfunc mismatch1 in build_batchesOfAOs',-1)
   IF(nbatLoc.NE.AO%nbatches)call lsquit('basfunc mismatch2 in build_batchesOfAOs',-1)
ELSE
   IF(norbitals.NE.AO%nbast)call lsquit('basfunc mismatch3 in build_batchesOfAOs',-1)
   IF(nbatLoc.NE.AO%nbatches)call lsquit('basfunc mismatch4 in build_batchesOfAOs',-1)
ENDIF
norbitals = 0
do I=1,nbatches
   do iOrb = 1,batchdim(I)
      norbitals=norbitals+1
      orbtoBatch(norbitals)=I      
   enddo   
enddo
MaxOrbitals = 0
do I=1,nbatches
   MaxOrbitals = MAX(MaxOrbitals,batchdim(I))
enddo

end subroutine build_batchesOfAOs1

subroutine determine_ActualDim(lupri,setting,maxallowedorbitals,nbast,&
     &MaxOrbitals,AOspec)
implicit none
integer,intent(in)         :: lupri,maxallowedorbitals,nbast
type(lssetting) :: setting
integer,intent(inout) :: MaxOrbitals
character(len=1),intent(in) :: AOspec
!
integer :: I,A,norbitals,nbatLoc,iOrb,tmporb,allocnbatches
logical :: uncont,intnrm,extend
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
uncont=.FALSE.
intnrm = .false.; extend=.false.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'D')THEN
   !   The Aux AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(AuxBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)
   extend = .true.
ELSEIF(AOspec.EQ.'O')THEN
   ! only The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(CABBasParam)
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
IF(extend)THEN
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,extend)
ELSE
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,extend)
ENDIF

norbitals = 0
MaxOrbitals = 0
do I=1,AO%nbatches
   tmporb = 0
   DO A=1,AO%BATCH(I)%nAngmom
      tmporb = tmporb + AO%BATCH(I)%norbitals(A)
   ENDDO
   IF(norbitals + tmporb.GT.maxallowedorbitals)THEN      
      MaxOrbitals = MAX(MaxOrbitals,norbitals)
      norbitals = tmporb
   ELSE
      norbitals = norbitals + tmporb
   ENDIF
enddo
IF(norbitals.GT.0)THEN
   MaxOrbitals = MAX(MaxOrbitals,norbitals)
ENDIF
call free_aoitem(lupri,AO)
end subroutine determine_ActualDim

subroutine DetermineBatchIndexAndSize(lupri,setting,startA,startB,startC,startD,&
     & ndimA,ndimB,ndimC,ndimD,batchsizeA,batchsizeB,batchsizeC,batchsizeD,&
     & batchindexA,batchindexB,batchindexC,batchindexD,&
     & offsetA,offsetB,offsetC,offsetD,ndimAs,ndimBs,ndimCs,ndimDs,AOspec)
implicit none
integer,intent(in)         :: lupri,startA,startB,startC,startD
integer,intent(in)         :: ndimA,ndimB,ndimC,ndimD
type(lssetting) :: setting
integer,intent(inout) :: batchsizeA,batchsizeB,batchsizeC,batchsizeD
integer,intent(inout) :: batchindexA,batchindexB,batchindexC,batchindexD
integer,intent(inout) :: ndimAs,ndimBs,ndimCs,ndimDs
integer,intent(inout) :: offsetA,offsetB,offsetC,offsetD
character(len=1),intent(in) :: AOspec
!
integer :: batchsize(4),batchindex(4),ndims(4),offset(4),start(4),ndim(4)
integer :: I,A,norbitals,iOrb,tmporb,J
logical :: uncont,intnrm,extend,OutsideBatch,Family
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2

start(1) = startA 
start(2) = startB 
start(3) = startC 
start(4) = startD 
ndim(1) = ndimA
ndim(2) = ndimB
ndim(3) = ndimC
ndim(4) = ndimD

uncont=.FALSE.
intnrm = .false.; extend=.false.
IF(AOspec.EQ.'R')THEN   
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)!   The regular AO-basis
ELSEIF(AOspec.EQ.'D')THEN   
   AObasis1 => setting%basis(1)%p%BINFO(AuxBasParam)!   The Aux AO-type basis
ELSEIF(AOspec.EQ.'C')THEN   
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)!   The CABS AO-type basis
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)!   The CABS AO-type basis
   extend=.true.
ELSEIF(AOspec.EQ.'O')THEN   
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)!   The CABS AO-type basis
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF

IF(extend)THEN
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,extend)
ELSE
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,extend)
ENDIF

Family = .FALSE.
do I=1,AO%nbatches
   IF(AO%BATCH(I)%nAngmom.GT.1) Family = .TRUE.
enddo
IF(Family)THEN
   call lsquit('Family basis set not allowed in build_minimalbatchesOfAOs',-1)
ENDIF

do J=1,4
   norbitals = 0
   OutsideBatch = .TRUE.
   do I=1,AO%nbatches
      tmporb = AO%BATCH(I)%norbitals(1)
      IF(OutsideBatch)THEN
         !determine if we go inside batch
         IF(norbitals + tmporb.GE.start(J))THEN !
            OutsideBatch = .FALSE.
            !first batch
            batchindex(J) = I
            batchsize(J) = 1
            ndims(J) = tmporb            
            offset(J) = start(J) - norbitals - 1
            norbitals = norbitals + tmporb
         ELSE
            !still outside batch
            norbitals = norbitals + tmporb
         ENDIF
      ELSE
         !inside batch          
         batchsize(J) = batchsize(J) + 1
         ndims(J) = ndims(J) + tmporb
         norbitals = norbitals + tmporb
         !determine if this is the last batch
         IF(norbitals.GE.start(J)+ndim(J)-1)THEN
            EXIT
         ENDIF
      ENDIF
   enddo
enddo
call free_aoitem(lupri,AO)
batchsizeA = batchsize(1)
batchsizeB = batchsize(2)
batchsizeC = batchsize(3)
batchsizeD = batchsize(4)
batchindexA = batchindex(1)
batchindexB = batchindex(2)
batchindexC = batchindex(3)
batchindexD = batchindex(4)
ndimAs = ndims(1)
ndimBs = ndims(2)
ndimCs = ndims(3)
ndimDs = ndims(4)
offsetA = offset(1)
offsetB = offset(2)
offsetC = offset(3)
offsetD = offset(4)
end subroutine DetermineBatchIndexAndSize

subroutine build_minimalbatchesOfAOs(lupri,setting,nbast,&
     & batchsize,batchdim,batchindex,nbatches,orbTobatch,AOspec)
implicit none
integer,intent(in)         :: lupri,nbast
type(lssetting) :: setting
integer,pointer :: batchdim(:),batchsize(:),batchindex(:)
integer :: orbtoBatch(nbast)
integer,intent(inout) :: nbatches
character(len=1),intent(in) :: AOspec
!
integer :: I,A,norbitals,iOrb,tmporb,allocnbatches
logical :: uncont,intnrm,extend,Family
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
uncont=.FALSE.
intnrm = .false.; extend=.false.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'D')THEN
   !   The Aux AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(AuxBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)
   extend=.true.
ELSEIF(AOspec.EQ.'O')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(CABBasParam)
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
Family = setting%SCHEME%NOFAMILY
setting%SCHEME%NOFAMILY = .FALSE.
IF(extend)THEN
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,extend)
ELSE
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,extend)
ENDIF
setting%SCHEME%NOFAMILY = Family

Family = .FALSE.
do I=1,AO%nbatches
   IF(AO%BATCH(I)%nAngmom.GT.1) Family = .TRUE.
enddo
IF(Family)THEN
   call lsquit('Family basis set not allowed in build_minimalbatchesOfAOs',-1)
ENDIF

nbatches = AO%nbatches
call mem_alloc(batchdim,nbatches)
call mem_alloc(batchsize,nbatches)
call mem_alloc(batchindex,nbatches)
do I=1,nbatches
   batchdim(I) = 0
   batchsize(I) = 1
   batchindex(I) = I
enddo
do I=1,AO%nbatches
   norbitals = 0
   DO A=1,AO%BATCH(I)%nAngmom
      norbitals = norbitals + AO%BATCH(I)%norbitals(A)
   ENDDO
   batchdim(I) = norbitals
enddo

norbitals = 0         
do I=1,nbatches
   norbitals=norbitals + batchdim(I)
enddo
IF(norbitals.NE.nbast)call lsquit('basfunc mismatch in build_minimalbatchesOfAOs',-1)

call free_aoitem(lupri,AO)
norbitals = 0
do I=1,nbatches
   do iOrb = 1,batchdim(I)
      norbitals=norbitals+1
      orbtoBatch(norbitals)=I      
   enddo   
enddo
end subroutine build_minimalbatchesOfAOs

subroutine build_minimalbatchesOfAOs2(lupri,setting,batchdim,nbatches,AOspec)
implicit none
integer,intent(in) :: lupri
type(lssetting) :: setting
integer,pointer :: batchdim(:)
integer,intent(inout) :: nbatches
character(len=1),intent(in) :: AOspec
!
integer :: I,A,norbitals,iOrb,tmporb,allocnbatches
logical :: uncont,intnrm,extend,Family
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
uncont=.FALSE.
intnrm = .false.; extend=.false.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'D')THEN
   !   The Aux AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(AuxBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)
   extend=.true.
ELSEIF(AOspec.EQ.'O')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(CABBasParam)
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
Family = setting%SCHEME%NOFAMILY
setting%SCHEME%NOFAMILY = .FALSE.
IF(extend)THEN
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,extend)
ELSE
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,extend)
ENDIF
setting%SCHEME%NOFAMILY = Family

Family = .FALSE.
do I=1,AO%nbatches
   IF(AO%BATCH(I)%nAngmom.GT.1) Family = .TRUE.
enddo
IF(Family)THEN
   call lsquit('Family basis set not allowed in build_minimalbatchesOfAOs',-1)
ENDIF

nbatches = AO%nbatches
call mem_alloc(batchdim,nbatches)
do I=1,nbatches
   batchdim(I) = 0
enddo
do I=1,AO%nbatches
   norbitals = 0
   DO A=1,AO%BATCH(I)%nAngmom
      norbitals = norbitals + AO%BATCH(I)%norbitals(A)
   ENDDO
   batchdim(I) = norbitals
enddo
call free_aoitem(lupri,AO)
end subroutine build_minimalbatchesOfAOs2

subroutine MPIdistributeAOs(setting,AOspec,nAuxMPI,numnodes,IndexToGlobal,&
     & MaxnAuxMPI,GindexToLocal,nAuxBasis)
implicit none
integer,intent(in)    :: numnodes,nAuxBasis
integer,pointer :: IndexToGlobal(:,:)
integer,intent(inout) :: nAuxMPI(numnodes),MaxnAuxMPI,GIndexToLocal(nAuxBasis)
type(lssetting) :: setting
character(len=1),intent(in) :: AOspec
!
integer :: I,A,idx(1),lupri,J,K,NA
integer,pointer :: load(:)
logical :: uncont,intnrm,extend
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
uncont=.FALSE.
intnrm = .false.; extend=.false.
IF(AOspec.EQ.'R')THEN      !    The regular AO-basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'D')THEN  !    The Aux AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(AuxBasParam)
ELSEIF(AOspec.EQ.'C')THEN  !    The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)
   extend=.true.
ELSEIF(AOspec.EQ.'O')THEN  !    The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(CABBasParam)
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
lupri = 6
IF(extend)THEN
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,extend)
ELSE
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,extend)
ENDIF

call mem_alloc(load,numnodes)
load = 0
nAuxMPI = 0 
do I=1,AO%nbatches
   idx = MINLOC(load)
   DO A=1,AO%BATCH(I)%nAngmom
      nAuxMPI(idx(1)) = nAuxMPI(idx(1)) + AO%BATCH(I)%norbitals(A)
      load(idx(1)) = load(idx(1)) + AO%BATCH(I)%norbitals(A)
   ENDDO
enddo
MaxnAuxMPI = MAXVAL(nAuxMPI)
call mem_alloc(IndexToGlobal,MaxnAuxMPI,numnodes)
load = 0
nAuxMPI = 0 
NA = 0
do I=1,AO%nbatches
   idx = MINLOC(load)
   load(idx(1)) = load(idx(1)) + (2*AO%BATCH(I)%maxAngmom+1)
   K=0
   DO A=1,AO%BATCH(I)%nAngmom
      DO J=1,AO%BATCH(I)%nOrbitals(A)
         nA = NA + 1
         nAuxMPI(idx(1)) = nAuxMPI(idx(1)) + 1
         GIndexToLocal(NA) = nAuxMPI(idx(1))
         IndexToGlobal(nAuxMPI(idx(1)),idx(1)) = AO%BATCH(I)%startOrbital(A)+K
         K=K+1
      ENDDO
   ENDDO
enddo
call mem_dealloc(load)
call free_aoitem(lupri,AO)
end subroutine MPIdistributeAOs

subroutine determine_MaxOrbitals(lupri,setting,maxallowedorbitals,MaxOrbitals,AOspec)
implicit none
integer,intent(in)    :: lupri,maxallowedorbitals
type(lssetting)       :: setting
integer,intent(inout) :: MaxOrbitals
character(len=1),intent(in) :: AOspec
!
integer :: I,A,norbitals,tmporb
logical :: uncont,intnrm,extend
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
uncont=.FALSE.
intnrm = .false.; extend=.false.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)
   extend=.true.
ELSEIF(AOspec.EQ.'O')THEN  !    The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(CABBasParam)
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
MaxOrbitals=0
IF(extend)THEN
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,extend)
ELSE
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,extend)
ENDIF

norbitals = 0
do I=1,AO%nbatches
   tmporb = 0
   DO A=1,AO%BATCH(I)%nAngmom
      tmporb = tmporb + AO%BATCH(I)%norbitals(A)
   ENDDO
   IF(norbitals + tmporb.GT.maxallowedorbitals)THEN
      norbitals = tmporb
      MaxOrbitals = MAX(MaxOrbitals,norbitals)
   ELSE
      norbitals = norbitals + tmporb
      MaxOrbitals = MAX(MaxOrbitals,norbitals)
   ENDIF
enddo
MaxOrbitals = MAX(MaxOrbitals,norbitals)
call free_aoitem(lupri,AO)
end subroutine determine_MaxOrbitals

subroutine determine_maxBatchOrbitalsize(lupri,setting,maxBatchOrbitalsize,AOspec)
implicit none
integer,intent(in)         :: lupri
integer,intent(inout)      :: maxBatchOrbitalsize
type(lssetting) :: setting
character(len=1),intent(in) :: AOspec
!
integer :: I,A,tmporb
logical :: uncont,intnrm,extend
type(AOITEM) :: AO
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
uncont=.FALSE.
intnrm = .false.; extend=.false.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'D')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(AuxBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(RegBasParam)
   AObasis2 => setting%basis(1)%p%BINFO(CABBasParam)
   extend=.true.
ELSEIF(AOspec.EQ.'O')THEN
   !   The CABS AO-type basis
   AObasis1 => setting%basis(1)%p%BINFO(CABBasParam)
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF

IF(extend)THEN
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,.FALSE.)
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis2,AO,uncont,intnrm,extend)
ELSE
   call build_AO(lupri,setting%scheme,setting%scheme%AOprint,&
        & setting%molecule(1)%p,AObasis1,AO,uncont,intnrm,extend)
ENDIF

maxBatchOrbitalsize = 0
do I=1,AO%nbatches
   tmporb = 0
   IF(AO%BATCH(I)%nAngmom.GT.1)THEN
      print*,'The Framework around determine_maxBatchOrbitalsize, '
      print*,'build_batchesofAOS and II_GET_DECPACKED4CENTER_J_ERI'
      print*,'requires .NOFAMILY keyword'
      CALL LSQUIT('NOFAMILY required in determine_maxBatchOrbitalsize',-1)
   ENDIF
   DO A=1,AO%BATCH(I)%nAngmom
      tmporb = tmporb + AO%BATCH(I)%norbitals(A)
   ENDDO
   maxBatchOrbitalsize = MAX(maxBatchOrbitalsize,tmporb) 
enddo
call free_aoitem(lupri,AO)
end subroutine determine_maxBatchOrbitalsize

!> \brief set the AO extent
!> \author S. Reine
!> \date 2010
SUBROUTINE SET_AO_EXTENT(AO,SCHEME)
implicit none
!> the AOitem to be build
TYPE(AOITEM)      :: AO
!> contains all info about the integralscheme requested (thresholds)
TYPE(LSINTSCHEME) :: SCHEME
!
Integer :: iBatch

DO iBatch=1,AO%nBatches
   AO%BATCH(iBatch)%extent = getAObatchExtent(AO%BATCH(iBatch),SCHEME)
ENDDO

END SUBROUTINE SET_AO_EXTENT

!> \brief calculate the AObatch extent
!> \author S. Reine
!> \date 2010
REAL(REALK) FUNCTION getAObatchExtent(AOB,SCHEME)
implicit none
!> the AObatch from which to calc extent
TYPE(AOBATCH)     :: AOB
!> contains all info about the integralscheme requested (thresholds)
TYPE(LSINTSCHEME) :: SCHEME
!
Integer     :: iAngmom,angmom,nPrim,nCont
Real(realk) :: extent

extent = 0E0_realk
nPrim = AOB%nPrimitives

DO iAngmom=1,AOB%nAngmom
  angmom = AOB%Angmom(iAngmom)
  nCont  = AOB%pCC(iAngmom)%p%ncol
  extent = max(extent,getExtent(AOB%pExponents%elms,AOB%pCC(iAngmom)%p%elms,&
     &           nPrim,nCont,angmom,SCHEME%OD_SCREEN,SCHEME%THRESHOLD*SCHEME%OD_THRESHOLD))
ENDDO

getAObatchExtent = extent

END FUNCTION getAObatchExtent

!> \brief calculate the extent of a single angular momentum for a single AObatch
!> \author S. Reine documented by T. Kjaergaard
!> \date 2010
REAL(REALK) FUNCTION getExtent(Exponents,CC,nPrim,nCont,angmom,Screen,ThresholdI)
implicit none
!> the number of primitive orbitals
Integer     :: nPrim
!> the number of contractedorbitals
Integer     :: nCont
!> the exponents 
Real(realk) :: Exponents(nPrim)
!> the contractioncoefficients 
Real(realk) :: CC(nPrim,nCont)
!> the angular momentum
Integer     :: angmom
!> if Overlap distribution screening ODscreening is used in this calc
Logical     :: Screen
!> the threshold for the Overlap distribution screening
Real(realk) :: ThresholdI
!
real(realk) :: Threshold
Integer :: iPrim,iCont
real(realk) :: extent,maxContraction,extent2,r2,fun,funD,rold,rnew
Threshold=ThresholdI
extent = 0E0_realk
IF (Screen) THEN
  extent2 = 0E0_realk
  DO iPrim=1,nPrim
    maxContraction = tiny(1E0_realk)
    DO iCont=1,nCont
      maxContraction = max(maxContraction,abs(CC(iPrim,iCont)))
    ENDDO
    IF ((Threshold.LT.tiny(1E0_realk)).OR.(nCont.LT. 1).OR.(Exponents(iPrim).LT.tiny(1E0_realk))) THEN
      WRITE(*,*) 'Error in FUNCTION getExtent', Threshold, nCont,Exponents(iPrim)
      CALL LSQUIT('Error in FUNCTION getExtent',-1)
    ENDIF
!    r2 = (-log(Threshold)+log(maxContraction)+log(real(nCont)))/Exponents(iPrim)
    r2 = (-log(Threshold)+log(maxContraction))/Exponents(iPrim)
!   Take the above expression to be a constant A. We should then solve 
!      r^2 = A + l/a*ln(r)       (i)
!   to get a proper extent. We instead take A to be an estimate of r^2
!   and make the correction r^2 = A + l/a*ln(sqrt(A)). The equation (i)
!   can of course instead be solved iteratively.
    IF (r2.GT. 0E0_realk) THEN
       IF(angmom .GT. 0)THEN
          Rold = sqrt(r2 + angmom*log(sqrt(r2))/Exponents(iPrim))
          fun = ABS(maxContraction*(rold**angmom)*exp(-Exponents(iPrim)*Rold**2))
          IF(threshold .LE. fun)THEN
             DO
                funD = ABS( maxContraction*angmom*(rold**(angmom-1))*&
                     &exp(-Exponents(iPrim)*Rold**2)+maxContraction*(rold**angmom)&
                     &*(-2*Exponents(iPrim)*Rold)*exp(-Exponents(iPrim)*Rold**2))
                Rnew = ROLD - (Threshold-fun)/funD
                IF(ABS(Rnew-ROLD) .LE. 1.0E-24_8)THEN
                   write(*,*)'calculation stalled in getExtent'
                   EXIT
                ENDIF
                fun = ABS(maxContraction*(Rnew)**angmom*exp(-Exponents(iPrim)*Rnew**2))
                Rold=Rnew
                IF(ABS(fun-Threshold).LE. Threshold*1.0E-11_8)EXIT
                IF(Rold.NE.Rold)THEN
                   write(*,*)'found NaN in aobatch iteratvie extent'
                   EXIT
                ENDIF
             ENDDO
          ENDIF
          r2=Rold*Rold
       ENDIF
    ENDIF
    extent2 = max(extent2,r2)
  ENDDO
  IF (extent2.LT. 0E0_realk) THEN
      WRITE(*,*) 'Negative squared distance in FUNCTION getExtent',extent2
      CALL LSQUIT('Negative squared distance in FUNCTION getExtent',-1)
  ENDIF
  extent = sqrt(extent2)
ENDIF
getExtent = extent

END FUNCTION getExtent

!> \brief determines if this segment of exponents are identical to earlier stored exponentes - used for family basis set
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE SHARES_EXPONENTS(LUPRI,IPRINT,AO,UNIQUEEXPONENTS,BASISINFO,type,B,J,UNIQUE,&
                                &UNCONTRACTED,MODELEXP,nMODELEXP)
implicit none
      !> the logical unit number for the output file
      INTEGER                   :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)           :: IPRINT
!> The AOitem that contains already processed exponents
TYPE(AOITEM)              :: AO
!> the number of unique exponents saved in the AO%Exponents 
INTEGER                   :: UNIQUEEXPONENTS
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO
!> Atomtype for the current atom being processed
INTEGER                   :: type
!> Angular moment for the current atom being processed
INTEGER                   :: B
!> The segment for the current atom being processed
INTEGER                   :: J
!> Output: either 0 if the exponents are unique or the index of the place in AO%Exponents where these exponents are already stored.
INTEGER :: UNIQUE
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> the number of differents segments of exponents saved in MODELEXP
INTEGER :: nMODELEXP
!> contains the index of where the differents segments of exponents are placed in AO%Exponents(I)
INTEGER :: MODELEXP(nMODELEXP)
!
INTEGER                   :: I,K,nPrim
real(realk)               :: SUM
REAL(REALK),PARAMETER     :: NULL=0E0_realk
IF(UNIQUEEXPONENTS .GT. 0) THEN
   IF(UNCONTRACTED)THEN
      UNIQUE=0
      DO I=1,UNIQUEEXPONENTS
         IF(ABS(BASISINFO%ATOMTYPE(type)%SHELL(B)%segment(J)%Exponents(1)-AO%Exponents(I)%elms(1)) < ExpThr)THEN
            UNIQUE=MODELEXP(I)       
         ENDIF
      ENDDO
   ELSE !DEFAULT
      UNIQUE=0
      nPrim=BASISINFO%ATOMTYPE(type)%SHELL(B)%segment(J)%nrow
      DO I=1,UNIQUEEXPONENTS
         IF(AO%Exponents(I)%nrow .EQ. nPrim)THEN
            SUM=NULL
            DO K=1,nPrim
               SUM=SUM+ABS(BASISINFO%ATOMTYPE(type)%SHELL(B)%segment(J)%Exponents(K)&
                    &-AO%Exponents(I)%elms(K))
            ENDDO
            IF (SUM<ExpThr) THEN
               UNIQUE=MODELEXP(I)       
            ENDIF
         ENDIF
      ENDDO
   ENDIF
ELSE
UNIQUE=0
ENDIF
IF(IPRINT .GT. 2)THEN
   IF(UNIQUE/=0) WRITE(LUPRI,*)'IT WAS DETERMINED THAT THIS SEGMENT SHARES'
   IF(UNIQUE/=0) WRITE(LUPRI,*)'EXPONENTS WITH SEGMENT:',UNIQUE
ENDIF

END SUBROUTINE SHARES_EXPONENTS

!> \brief add batch to aoitem
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE ADD_BATCH(SCHEME,MOLECULE,BASISINFO,AO,AOorganize,iATOM,type,B,J,&
     &lupri,iprint,nPrim,nCont,nbatches,UNCONTRACTED,itype)
use molecule_type
implicit none
!> type index
INTEGER                   :: itype
!> contains all info about the integralscheme requested (thresholds,use cartesian,..)
TYPE(LSINTSCHEME)         :: SCHEME
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO 
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> contains info about the full AOorganize
TYPE(AOorganizer)         :: AOorganize
!> the atomtype for the current atom
INTEGER                   :: type
!> the angular momentum for the current atom
INTEGER                   :: B
!> the segment for the current atom
INTEGER                   :: J
!> the logical unit number for the output file
INTEGER                   :: lupri
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: iprint
!> number of contracted functions
INTEGER,intent(out)       :: nCont
!> number of primitive functions
INTEGER,intent(out)       :: nPrim
!> nbatches = old number of batches + new number of batches
INTEGER                   :: nbatches
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!
INTEGER                   :: iATOM,I

IF(UNCONTRACTED)THEN
   !WE ADD ONE BATCH FOR EACH PRIMITIVES
   nPrim=BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
   nCont=nPrim
   DO I=1,nPrim
      nbatches = AO%nbatches + I
      itype = itype + 1
      AOorganize%ORG(nbatches)%type=itype
      AOorganize%ORG(nbatches)%angmoms=1
      AOorganize%ORG(nbatches)%ANGMOM(1)=B-1
      AOorganize%ORG(nbatches)%Segment(1)=J
      AOorganize%ORG(nbatches)%atom=iATOM
      AO%BATCH(nbatches)%itype = itype
      AO%BATCH(nbatches)%nAngmom = 1
      AO%BATCH(nbatches)%molecularIndex = MOLECULE%ATOM(iATOM)%molecularIndex
      AO%BATCH(nbatches)%type_Empty  = .FALSE. 
      AO%BATCH(nbatches)%type_Nucleus  = .FALSE.
      AO%BATCH(nbatches)%type_elField = .FALSE.
      AO%BATCH(nbatches)%type_pCharge = .FALSE.
      AO%BATCH(nbatches)%spherical= SCHEME%DoSpherical
      AO%BATCH(nbatches)%Angmom(1)= B - 1
      AO%BATCH(nbatches)%maxAngmom= B - 1
      AO%BATCH(nbatches)%CENTER(1)=MOLECULE%ATOM(iATOM)%CENTER(1)
      AO%BATCH(nbatches)%CENTER(2)=MOLECULE%ATOM(iATOM)%CENTER(2)
      AO%BATCH(nbatches)%CENTER(3)=MOLECULE%ATOM(iATOM)%CENTER(3)
   ENDDO
   AO%nbatches = nbatches
   AOorganize%nbatches = nbatches
ELSE !DEFAULT
   nbatches = AO%nbatches + 1
   itype = itype + 1
   AO%nbatches = nbatches 
   AOorganize%nbatches = nbatches
   AOorganize%ORG(nbatches)%type=itype
   AOorganize%ORG(nbatches)%angmoms=1
   AOorganize%ORG(nbatches)%ANGMOM(1)=B-1
   AOorganize%ORG(nbatches)%Segment(1)=J
   AOorganize%ORG(nbatches)%atom=iATOM
   AO%BATCH(nbatches)%itype = itype
   AO%BATCH(nbatches)%nAngmom = 1
   AO%BATCH(nbatches)%molecularIndex = MOLECULE%ATOM(iATOM)%molecularIndex
   AO%BATCH(nbatches)%type_Empty  = .FALSE. 
   AO%BATCH(nbatches)%type_Nucleus  = .FALSE.
   AO%BATCH(nbatches)%type_elField = .FALSE.
   AO%BATCH(nbatches)%type_pCharge  = .FALSE.
   AO%BATCH(nbatches)%spherical= SCHEME%DoSpherical
   AO%BATCH(nbatches)%Angmom(1)= B - 1
   AO%BATCH(nbatches)%maxAngmom= B - 1
   AO%BATCH(nbatches)%CENTER(1)=MOLECULE%ATOM(iATOM)%CENTER(1)
   AO%BATCH(nbatches)%CENTER(2)=MOLECULE%ATOM(iATOM)%CENTER(2)
   AO%BATCH(nbatches)%CENTER(3)=MOLECULE%ATOM(iATOM)%CENTER(3)
   nPrim=BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
   nCont=BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
ENDIF

END SUBROUTINE ADD_BATCH

!> \brief Does different initializations before get_density is called.
!> \author T. Kjaergaard
!> \date 2008
!>
!> Add a segment of exponents and the corresponding contraction coefficients  
!> to the unique exponents and contraction coefficients
!>
SUBROUTINE ADD_SEGMENT_AND_CC(AO,AOmodel,AOorganize,BASISINFO,nbatches,nPrim,nCont,&
     & type,B,J,lupri,iprint,UNIQUEMATRIXES,UNIQUEEXPONENTS,UNCONTRACTED,&
     & INTNRM,MODELEXP,nMODELEXP,NORM)
IMPLICIT NONE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the model AOitem for the atom 
TYPE(AOITEM)              :: AOmodel
!> the AOorganizer that contain the unique exponents
TYPE(AOorganizer)         :: AOorganize 
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO 
!> the number of batches already in the full AO
INTEGER                   :: nbatches
!> the number of primitive functions
INTEGER                   :: nPrim
!> the number of contracted functions
INTEGER                   :: nCont
!> the atomtype for the current atom
INTEGER                   :: type
!> the angular momentum for the current atom
INTEGER                   :: B
!> the segment for the current atom
INTEGER                   :: J
!> the logical unit number for the output file
INTEGER                   :: lupri
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: iprint
!> the number of unique contractioncoefficientmatrices saved
INTEGER                   :: UNIQUEMATRIXES
!> the number of unique exponents saved
INTEGER                   :: UNIQUEEXPONENTS
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL                   :: INTNRM 
!> the number of differents segments of exponents saved in MODELEXP
INTEGER :: nMODELEXP
!> contains the index of where the differents segments of exponents are placed in AO%Exponents(I)
! use normalized contractioncoefficients
LOGICAL,intent(in)       :: NORM
INTEGER :: MODELEXP(nMODELEXP)
INTEGER                   :: nbat,I,icont,iprim
REAL(REALK)               :: PI,EXPO,PIPPI,maxcont,element
PI=3.14159265358979323846E0_realk

IF(UNCONTRACTED)THEN
 PIPPI = (0.5E0_realk/PI)**(0.75E0_realk)
 nbat = nbatches-nPrim
 IF(INTNRM)THEN
   !DETERMINE MAX CONTRACTION ELEMENT
   MAXcont=0E0_realk
   DO iprim=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
     DO iCont=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
       element=ABS(BASISINFO%ATOMTYPE(type)&
       &%SHELL(B)%SEGMENT(J)%elms(iprim+(iCont-1)*nPrim))
       MAXcont=MAX(MAXcont,element)
     ENDDO
   ENDDO
 ENDIF

 DO I=1,nPrim
    nbatches = nbat+I
    UNIQUEMATRIXES=UNIQUEMATRIXES+1
    UNIQUEEXPONENTS=UNIQUEEXPONENTS+1
    AOmodel%nExp=AOmodel%nExp+1
    MODELEXP(AOmodel%nExp)=UNIQUEEXPONENTS
    AOorganize%ORG(nbatches)%exponentindex=UNIQUEEXPONENTS
    AOorganize%ORG(nbatches)%CC(1)=UNIQUEMATRIXES
    CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),1,1)
    CALL lsmat_dense_init(AO%Exponents(UNIQUEEXPONENTS),1,1)
    CALL lsmat_dense_init(AOmodel%Exponents(AOmodel%nExp),1,1)
    
    Expo = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%Exponents(I)

    AO%Exponents(UNIQUEEXPONENTS)%elms(1)= Expo
    AOmodel%Exponents(AOmodel%nExp)%elms(1)= Expo

    IF(IPRINT .GT. 2)THEN
       WRITE(LUPRI,*)'THE NEW EXPONENTS',UNIQUEEXPONENTS
       call lsmat_dense_print(AO%Exponents(UNIQUEEXPONENTS),&
            &1,1,1,1,lupri)
    ENDIF
    
    IF(INTNRM)THEN
       AO%CC(UNIQUEMATRIXES)%elms(1)= MAXcont*&
            & (4*Expo)**(0.5E0_realk*B+0.25E0_realk)*PIPPI 
    ELSE
       AO%CC(UNIQUEMATRIXES)%elms(1)= &
            & (4*Expo)**(0.5E0_realk*B+0.25E0_realk)*PIPPI 
    ENDIF
    AO%angmom(UNIQUEMATRIXES) = B-1
    
    IF(IPRINT .GT. 2)THEN
       WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
       call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
            &1,1,1,1,lupri)
    ENDIF
 ENDDO
ELSE !DEFAULT
   IF(INTNRM)nCont = nPrim
   UNIQUEMATRIXES=UNIQUEMATRIXES+1
   UNIQUEEXPONENTS=UNIQUEEXPONENTS+1
   AOmodel%nExp=AOmodel%nExp+1
   MODELEXP(AOmodel%nExp)=UNIQUEEXPONENTS
   AOorganize%ORG(nbatches)%exponentindex=UNIQUEEXPONENTS
   AOorganize%ORG(nbatches)%CC(1)=UNIQUEMATRIXES
   CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),nPrim,nCont)
   CALL lsmat_dense_init(AO%Exponents(UNIQUEEXPONENTS),nPrim,1)
   CALL lsmat_dense_init(AOmodel%Exponents(AOmodel%nExp),nPrim,1)

   IF(INTNRM)THEN
      PIPPI = (0.5E0_realk/PI)**(0.75E0_realk)
      !DETERMINE MAX CONTRACTION ELEMENT
      MAXcont=0E0_realk
      DO iCont=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
         DO iprim=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
            element=ABS(BASISINFO%ATOMTYPE(type)&
                 &%SHELL(B)%SEGMENT(J)%elms(iprim+(iCont-1)*nPrim))
            MAXcont=MAX(MAXcont,element)
         ENDDO
      ENDDO
      CALL lsmat_dense_zero(AO%CC(UNIQUEMATRIXES))
      DO I=1,nPrim
         Expo = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%Exponents(I)
         AO%CC(UNIQUEMATRIXES)%elms(I+(I-1)*nprim)= MAXcont*&
              & (4*Expo)**(0.5E0_realk*B+0.25E0_realk)*PIPPI 
      ENDDO
   ELSE
      IF(NORM)THEN
         call dcopy(nPrim*nCont,BASISINFO%ATOMTYPE(type)%SHELL(B)%&
              &SEGMENT(J)%elms,1,AO%CC(UNIQUEMATRIXES)%elms,1)
      ELSE
         call dcopy(nPrim*nCont,BASISINFO%ATOMTYPE(type)%SHELL(B)%&
              &SEGMENT(J)%UCCelms,1,AO%CC(UNIQUEMATRIXES)%elms,1)
      ENDIF
   ENDIF
   AO%angmom(UNIQUEMATRIXES) = B-1
   call dcopy(nPrim,BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%Exponents,1,&
        &AO%Exponents(UNIQUEEXPONENTS)%elms,1)
   call dcopy(nPrim,BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%Exponents,1,&
        &AOmodel%Exponents(AOmodel%nExp)%elms,1)
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
      call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
           &1,nPrim,1,nCont,lupri)
   ENDIF
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW EXPONENTS',UNIQUEEXPONENTS
      call lsmat_dense_print(AO%Exponents(UNIQUEEXPONENTS),&
           &1,nPrim,1,1,lupri)
   ENDIF
ENDIF

END SUBROuTINE ADD_SEGMENT_AND_CC

!> \brief Does different initializations before get_density is called.
!> \author T. Kjaergaard
!> \date 2008
!>
!> Add a segment of exponents and the corresponding contraction coefficients  
!> to the unique exponents and contraction coefficients
!>
SUBROUTINE ADD_SINGLESEGMENT_AND_CC(AO,AOmodel,AOorganize,BASISINFO,nbatches,nPrim,nCont,&
     & type,B,lupri,iprint,UNIQUEMATRIXES,UNIQUEEXPONENTS,UNCONTRACTED,&
     & INTNRM,MODELEXP,nMODELEXP,NORM)
IMPLICIT NONE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the model AOitem for the atom 
TYPE(AOITEM)              :: AOmodel
!> the AOorganizer that contain the unique exponents
TYPE(AOorganizer)         :: AOorganize 
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASISINFO 
!> the number of batches already in the full AO
INTEGER                   :: nbatches
!> the number of primitive functions
INTEGER                   :: nPrim
!> the number of contracted functions
INTEGER                   :: nCont
!> the atomtype for the current atom
INTEGER                   :: type
!> the angular momentum for the current atom
INTEGER                   :: B
!> the segment for the current atom
INTEGER                   :: J
!> the logical unit number for the output file
INTEGER                   :: lupri
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: iprint
!> the number of unique contractioncoefficientmatrices saved
INTEGER                   :: UNIQUEMATRIXES
!> the number of unique exponents saved
INTEGER                   :: UNIQUEEXPONENTS
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> if INTERMIDIATE NORMALIZATION shoule be used, USED FOR SCREENING
LOGICAL                   :: INTNRM 
!> the number of differents segments of exponents saved in MODELEXP
INTEGER :: nMODELEXP
!> contains the index of where the differents segments of exponents are placed in AO%Exponents(I)
! use normalized contractioncoefficients
LOGICAL,intent(in)       :: NORM
INTEGER :: MODELEXP(nMODELEXP)
INTEGER                   :: nbat,I,icont,iprim,iprimLoc,icontloc,nprimloc,ielm
integer                   :: ncontloc
REAL(REALK)               :: PI,EXPO,PIPPI,maxcont,element
PI=3.14159265358979323846E0_realk

IF(UNCONTRACTED)THEN
   call lsquit('Uncontracted and nosegment - does not work at the same time - what did you expect',lupri)
ELSE !DEFAULT
   IF(INTNRM)nCont = nPrim
   UNIQUEMATRIXES=UNIQUEMATRIXES+1
   UNIQUEEXPONENTS=UNIQUEEXPONENTS+1
   AOmodel%nExp=AOmodel%nExp+1
   MODELEXP(AOmodel%nExp)=UNIQUEEXPONENTS
   AOorganize%ORG(nbatches)%exponentindex=UNIQUEEXPONENTS
   AOorganize%ORG(nbatches)%CC(1)=UNIQUEMATRIXES
!   write(lupri,*)'nPrim,nCont',nPrim,nCont
   CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),nPrim,nCont)
   CALL lsmat_dense_init(AO%Exponents(UNIQUEEXPONENTS),nPrim,1)
   CALL lsmat_dense_init(AOmodel%Exponents(AOmodel%nExp),nPrim,1)

   IF(INTNRM)THEN
      PIPPI = (0.5E0_realk/PI)**(0.75E0_realk)
      !DETERMINE MAX CONTRACTION ELEMENT
      MAXcont=0E0_realk
      DO J = 1,BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
       DO iCont=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
        DO iprim=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
           element=ABS(BASISINFO%ATOMTYPE(type)&
                &%SHELL(B)%SEGMENT(J)%elms(iprim+(iCont-1)*nPrim))
           MAXcont=MAX(MAXcont,element)
        ENDDO
       ENDDO
      ENDDO
      CALL lsmat_dense_zero(AO%CC(UNIQUEMATRIXES))
      iprim = 0
      do J = 1,BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
         do iprimLoc = 1, BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
            Expo = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%Exponents(iprimLoc)
            AO%CC(UNIQUEMATRIXES)%elms(iprimLoc+iprim+(iprimLoc+iprim-1)*nprim)= MAXcont*&
                 & (4*Expo)**(0.5E0_realk*B+0.25E0_realk)*PIPPI 
         ENDDO
         iprim = iprim + BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
      ENDDO
   ELSE
      CALL lsmat_dense_zero(AO%CC(UNIQUEMATRIXES))
      iprim = 0
      icont = 0
      DO J = 1,BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
         nprimLoc = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
         ncontloc = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
!         WRITE(LUPRI,*)'THE J=',J,'CCelms',nprimloc,ncontloc
!         call ls_output(BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%elms,&
!              &1,nPrimloc,1,nContloc,nPrimloc,nContloc,1,lupri)
         IF(NORM)THEN
            DO iprimLoc=1,nPrimLoc
               DO iContLoc=1,ncontloc
                  ielm = iprim+iprimloc+(icont+icontloc-1)*nprim
                  !               write(lupri,*)'ielm',ielm,'nprim*ncont',nprim*ncont,'iprimLoc+(icontLoc-1)*nPrimLoc',&
                  !                    &nPrimLoc*BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
                  AO%CC(UNIQUEMATRIXES)%elms(ielm) = &
                       & BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%elms(iprimLoc+(icontLoc-1)*nPrimLoc)
               ENDDO
            ENDDO
         ELSE
            DO iprimLoc=1,nPrimLoc
               DO iContLoc=1,ncontloc
                  ielm = iprim+iprimloc+(icont+icontloc-1)*nprim
                  !               write(lupri,*)'ielm',ielm,'nprim*ncont',nprim*ncont,'iprimLoc+(icontLoc-1)*nPrimLoc',&
                  !                    &nPrimLoc*BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
                  AO%CC(UNIQUEMATRIXES)%elms(ielm) = &
                       & BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%UCCelms(iprimLoc+(icontLoc-1)*nPrimLoc)
               ENDDO
            ENDDO
         ENDIF
         iprim = iprim + nPrimLoc
         icont = icont + BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
      ENDDO
!      WRITE(LUPRI,*)'THE new CCelms',nprim,ncont
!      call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
!           &1,nPrim,1,nCont,lupri)

!      call dcopy(nPrim*nCont,BASISINFO%ATOMTYPE(type)%SHELL(B)%&
!           &SEGMENT(J)%elms,1,AO%CC(UNIQUEMATRIXES)%elms,1)
   ENDIF
   AO%angmom(UNIQUEMATRIXES) = B-1

   iprim = 0
   do J = 1,BASISINFO%ATOMTYPE(type)%SHELL(B)%nsegments
      nprimLoc = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
      do iprimLoc = 1,nprimLoc
         AO%Exponents(UNIQUEEXPONENTS)%elms(iprim+iPrimLoc)=&
              &BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%Exponents(iprimLoc)
         AOmodel%Exponents(AOmodel%nExp)%elms(iprim+iPrimLoc)=&
              &BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%Exponents(iprimLoc)
      ENDDO
      iprim = iprim + nprimLoc
   ENDDO
!   call dcopy(nPrim,BASISINFO%ATOMTYPE(type)&
!        &%SHELL(B)%SEGMENT(J)%Exponents,1,&
!        &AO%Exponents(UNIQUEEXPONENTS)%elms,1)
!   call dcopy(nPrim,BASISINFO%ATOMTYPE(type)&
!        &%SHELL(B)%SEGMENT(J)%Exponents,1,&
!        &AOmodel%Exponents(AOmodel%nExp)%elms,1)
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
      call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
           &1,nPrim,1,nCont,lupri)
   ENDIF
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW EXPONENTS',UNIQUEEXPONENTS
      call lsmat_dense_print(AO%Exponents(UNIQUEEXPONENTS),&
           &1,nPrim,1,1,lupri)
   ENDIF
ENDIF

END SUBROuTINE ADD_SINGLESEGMENT_AND_CC

!> \brief Direct pointers to point to the unique exponents and contraction coefficients. 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DIRECT_POINTERS(AO,AOmodel,AOorganize,nCont,nPrim,nbatches,B,lupri,iprint,orbitalIndex,&
     &primorbitalindex,UNCONTRACTED,INTNRM)
IMPLICIT NONE
!> the AOitem to be build
TYPE(AOITEM)              :: AO
!> the AOitem for the current atom
TYPE(AOITEM)              :: AOmodel
!> the AOorganizer
TYPE(AOorganizer)         :: AOorganize 
!> the number of contracted functions
INTEGER                   :: nCont
!> the number of primitive functions
INTEGER                   :: nPrim
!> the number of batches already in the full AO
INTEGER                   :: nbatches
!> the angular momentum for the current atom
INTEGER                   :: B
!> the logical unit number for the output file
INTEGER                   :: lupri
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: iprint
!> the orbital index
INTEGER                   :: orbitalIndex
!> the primitive orbital index
INTEGER                   :: primorbitalindex
!> if the AOitem should be uncontracted
LOGICAL                   :: UNCONTRACTED
!> if the AOitem should be intermediate normalized
LOGICAL                   :: iNTNRM
INTEGER                   :: nOrbitals,I,nbat,A

IF(UNCONTRACTED)THEN
   nbat = nbatches-nPrim
   DO I=1,nPrim
      nbatches = nbat+I
      A=AOmodel%BATCH(nbatches)%nangmom
      IF(A .GT. maxAOangmom)THEN
         WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
              & - more than just the normal s, p, d, f, g, h, i, j, k and l&
              & orbitals so you need to increase the maxAOangmom parameter&
              & in TYPE-DEF.f90 file and recompile. It should be even to or&
              & greater than the number of different orbitals you need.&
              & Thomas Kjaergaard'
         CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
      ENDIF
      AOmodel%BATCH(nbatches)%nContracted(A) = 1
      AOmodel%BATCH(nbatches)%maxContracted  = 1
      AOmodel%BATCH(nbatches)%nPrimitives=1
   
      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to ',AOorganize%ORG(nbatches)%CC(1)
   
      AOmodel%BATCH(nbatches)%pCC(A)%p=>&
        &AO%CC(AOorganize%ORG(nbatches)%CC(1))
      AOmodel%BATCH(nbatches)%CCindex(A)=AOorganize%ORG(nbatches)%CC(1)
      IF(IPRINT .GT. 30) THEN
         WRITE(LUPRI,*)'Which is'
         CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pCC(A)%p,1,1,1,1,lupri)
      ENDIF
   
      IF(IPRINT .GT. 2)WRITE(LUPRI,*)'Exp pointer to ',AOorganize%ORG(nbatches)%exponentindex

      AOmodel%BATCH(nbatches)%Expindex=AOorganize%ORG(nbatches)%exponentindex
      AOmodel%BATCH(nbatches)%pExponents=>&
           &AO%Exponents(AOorganize%ORG(nbatches)%exponentindex)

      IF(IPRINT .GT. 30) THEN
         WRITE(LUPRI,*)'Which is'
!         CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pExponents,1,1,1,1,lupri)
         CALL LSMAT_DENSE_PRINT(AO%Exponents(AOorganize%ORG(nbatches)%exponentindex),1,1,1,1,lupri)
      ENDIF
   
      ! Orbital index (used in for example density-matrix, Fock/KS-matrix,
      ! fitting coefficient etc.)
      AOmodel%BATCH(nbatches)%startOrbital(A) = orbitalIndex
      AOmodel%BATCH(nbatches)%startprimOrbital(A) = primorbitalIndex
      AOmodel%BATCH(nbatches)%nPrimOrbComp(A) = B*(B+1)/2 
      IF (AOmodel%BATCH(nbatches)%spherical) THEN
         AOmodel%BATCH(nbatches)%nOrbComp(A) = 2*B-1
      ELSE
         AOmodel%BATCH(nbatches)%nOrbComp(A) = B*(B+1)/2
      ENDIF
      nOrbitals = AOmodel%BATCH(nbatches)%nOrbComp(A)
      orbitalIndex = orbitalIndex + nOrbitals
      primorbitalIndex = primorbitalIndex + nOrbitals
      AOmodel%BATCH(nbatches)%nOrbitals(A)    = nOrbitals
   ENDDO
ELSE
   A=AOmodel%BATCH(nbatches)%nangmom
   IF(A .GT. maxAOangmom)THEN
      WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
           & - more than just the normal s, p, d, f, g, h, i, j, k and l&
           & orbitals so you need to increase the maxAOangmom parameter&
           & in TYPE-DEF.f90 file and recompile. It should be even to or&
           & greater than the number of different orbitals you need.&
           & Thomas Kjaergaard'
      CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
   ENDIF
   IF(INTNRM)THEN
      AOmodel%BATCH(nbatches)%nContracted(A) = nPrim
      AOmodel%BATCH(nbatches)%maxContracted  = nPrim
   ELSE
      AOmodel%BATCH(nbatches)%nContracted(A) = nCont
      AOmodel%BATCH(nbatches)%maxContracted  = nCont
   ENDIF
   AOmodel%BATCH(nbatches)%nPrimitives=nPrim
   
   IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to ',AOorganize%ORG(nbatches)%CC(1)
   
   AOmodel%BATCH(nbatches)%pCC(A)%p=>&
        &AO%CC(AOorganize%ORG(nbatches)%CC(1))
   AOmodel%BATCH(nbatches)%CCindex(A)=AOorganize%ORG(nbatches)%CC(1)
   AOmodel%BATCH(nbatches)%Expindex=AOorganize%ORG(nbatches)%exponentindex
   
   IF(IPRINT .GT. 30) THEN
      WRITE(LUPRI,*)'Which is'
      CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pCC(A)%p,1,nPrim,1,nCont,lupri)
   ENDIF
   
   
   ! Orbital index (used in for example density-matrix, Fock/KS-matrix,
   ! fitting coefficient etc.)
   AOmodel%BATCH(nbatches)%startOrbital(A) = orbitalIndex
   AOmodel%BATCH(nbatches)%startprimOrbital(A) = primorbitalIndex
   AOmodel%BATCH(nbatches)%nPrimOrbComp(A) = B*(B+1)/2
   IF (AOmodel%BATCH(nbatches)%spherical) THEN
      AOmodel%BATCH(nbatches)%nOrbComp(A) = 2*B-1
   ELSE
      AOmodel%BATCH(nbatches)%nOrbComp(A) = B*(B+1)/2
   ENDIF
   nOrbitals = AOmodel%BATCH(nbatches)%nOrbComp(A)*nCont
   orbitalIndex = orbitalIndex + nOrbitals
   primorbitalIndex = primorbitalIndex &
        &+ AOmodel%BATCH(nbatches)%nOrbComp(A)*nprim
   AOmodel%BATCH(nbatches)%nOrbitals(A)    = nOrbitals
   IF(IPRINT .GT. 2)WRITE(LUPRI,*)'Exp pointer to ',AOorganize%ORG(nbatches)%exponentindex
   AOmodel%BATCH(nbatches)%pExponents=>&
        &AO%Exponents(AOorganize%ORG(nbatches)%exponentindex)
   
   IF(IPRINT .GT. 30) THEN
      WRITE(LUPRI,*)'Which is'
      CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbatches)%pExponents,1,nPrim,1,1,lupri)
   ENDIF
ENDIF
END SUBROUTINE DIRECT_POINTERS

!> \brief Determine the number of batches for the atom index atom 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DETERMINE_nbathces(AOorganize,atom,UNIQUE,nbatches)
implicit none
TYPE(AOorganizer) :: AOorganize
INTEGER           :: atom,nbatches,I,UNIQUE
nbatches=0
DO I=AOorganize%nbatches,1,-1
   IF(AOorganize%ORG(I)%atom==atom)THEN
      IF(AOorganize%ORG(I)%exponentindex==UNIQUE)THEN
         nbatches=I
         RETURN
      ENDIF
   ENDIF
ENDDO

END SUBROUTINE DETERMINE_nbathces

!> \brief adding extra angular moments to the already existing BATCH(family basisset)
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE EXPAND_BATCH(AO,AOorganize,BASISINFO,type,B,J,nbat,A,lupri,&
     &UNCONTRACTED)
IMPLICIT NONE
INTEGER                   :: A,lupri,I,nprim
INTEGER,intent(in)        :: B,nbat,J,type
TYPE(BASISSETINFO)        :: BASISINFO 
TYPE(AOITEM)              :: AO
TYPE(AOorganizer)         :: AOorganize 
LOGICAL                   :: UNCONTRACTED

IF(UNCONTRACTED)THEN
   nPrim=BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%nrow
   DO I=0,nPrim-1
      A=AO%BATCH(nbat+I)%nAngmom+1
      IF(A .GT. maxAOangmom)THEN
        WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
             & - more than just the normal s, p, d, f, g, h, i, j, k and l&
             & orbitals so you need to increase the maxAOangmom parameter&
             & in TYPE-DEF.f90 file and recompile. It should be even to or&
             & greater than the number of different orbitals you need.&
             & Thomas Kjaergaard'
        CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
      ENDIF
      AO%BATCH(nbat+I)%nAngmom = A
      AO%BATCH(nbat+I)%Angmom(A)   = B - 1
      AO%BATCH(nbat+I)%maxAngmom = B - 1!MAX(AO%BATCH(UNIQUE)%maxAngmom,B)
      AOorganize%ORG(nbat+I)%angmoms= A
      AOorganize%ORG(nbat+I)%ANGMOM(A)=B-1
   ENDDO
ELSE
   A=AO%BATCH(nbat)%nAngmom+1
   IF(A .GT. maxAOangmom)THEN
      WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
           & - more than just the normal s, p, d, f, g, h, i, j, k and l&
           & orbitals so you need to increase the maxAOangmom parameter&
           & in TYPE-DEF.f90 file and recompile. It should be even to or&
           & greater than the number of different orbitals you need.&
           & Thomas Kjaergaard'
      CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
   ENDIF
   AO%BATCH(nbat)%nAngmom = A
   AO%BATCH(nbat)%Angmom(A)   = B - 1
   AO%BATCH(nbat)%maxAngmom = B - 1!MAX(AO%BATCH(UNIQUE)%maxAngmom,B)
   AOorganize%ORG(nbat)%angmoms= A
   AOorganize%ORG(nbat)%ANGMOM(A)=B-1
ENDIF

END SUBROUTINE EXPAND_BATCH

!> \brief add a contraction coefficient matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ADD_CC(AO,AOorganize,nbatches,nPrim,nCont,lupri,iprint,&
     &BASISINFO,type,B,J,A,UNIQUEMATRIXES,UNCONTRACTED,INTNRM,NORM)
IMPLICIT NONE
! use normalized contractioncoefficients
LOGICAL,intent(in)       :: NORM
TYPE(AOITEM)              :: AO
TYPE(AOorganizer)         :: AOorganize
TYPE(BASISSETINFO)        :: BASISINFO 
INTEGER                   :: UNIQUEMATRIXES,nbatches,nPrim
INTEGER                   :: lupri,nCont,iprint
INTEGER                   :: type,B,J,A,I,icont,iprim
REAL(REALK)               :: PI,EXP,PIPPI,maxcont,element,Expo
LOGICAL                   :: UNCONTRACTED,INTNRM
PI=3.14159265358979323846E0_realk

   nPrim=BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%nrow
   nCont=BASISINFO%ATOMTYPE(type)%SHELL(B)&
        &%SEGMENT(J)%ncol

IF(UNCONTRACTED)THEN  
   PIPPI = (0.5E0_realk/PI)**(0.75E0_realk)
   IF(INTNRM)THEN
      !DETERMINE MAX CONTRACTION ELEMENT
      MAXcont=0E0_realk
      DO iprim=1,nPrim
         DO iCont=1,nCont
            MAXcont=MAX(MAXcont,ABS(BASISINFO%ATOMTYPE(type)&
                 &%SHELL(B)%SEGMENT(J)%elms(iprim+(iCont-1)*nPrim)))
         ENDDO
      ENDDO
   ENDIF
   DO I=0,nPrim-1
      UNIQUEMATRIXES=UNIQUEMATRIXES+1
      AOorganize%ORG(nbatches+I)%CC(A)=UNIQUEMATRIXES
      CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),1,1)
      Exp = BASISINFO%ATOMTYPE(type)&
        &%SHELL(B)%SEGMENT(J)%Exponents(I+1)
      IF(INTNRM)THEN
         AO%CC(UNIQUEMATRIXES)%elms(1)= MAXcont* &
              & (4*Exp)**(0.5E0_realk*B+0.25E0_realk)*PIPPI 
      ELSE
         AO%CC(UNIQUEMATRIXES)%elms(1)=&
              & (4*Exp)**(0.5E0_realk*B+0.25E0_realk)*PIPPI 
      ENDIF
      AO%angmom(UNIQUEMATRIXES) = B-1
      IF(IPRINT .GT. 2)THEN
         WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
         call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
              &1,1,1,1,lupri)
      ENDIF
   ENDDO
ELSE
   if(intnrm)ncont = nprim
   UNIQUEMATRIXES=UNIQUEMATRIXES+1
   AOorganize%ORG(nbatches)%CC(A)=UNIQUEMATRIXES
   CALL lsmat_dense_init(AO%CC(UNIQUEMATRIXES),nPrim,nCont)
   IF(INTNRM)THEN
      PIPPI = (0.5E0_realk/PI)**(0.75E0_realk)
      !DETERMINE MAX CONTRACTION ELEMENT
      MAXcont=0E0_realk
      DO iCont=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%ncol
         DO iprim=1,BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%nrow
            element=ABS(BASISINFO%ATOMTYPE(type)&
                 &%SHELL(B)%SEGMENT(J)%elms(iprim+(iCont-1)*nPrim))
            MAXcont=MAX(MAXcont,element)
         ENDDO
      ENDDO
      CALL lsmat_dense_zero(AO%CC(UNIQUEMATRIXES))
      DO I=1,nPrim
         Expo = BASISINFO%ATOMTYPE(type)%SHELL(B)%SEGMENT(J)%Exponents(I)
         AO%CC(UNIQUEMATRIXES)%elms(I+(I-1)*nprim)= MAXcont*&
              & (4*Expo)**(0.5E0_realk*B+0.25E0_realk)*PIPPI 
      ENDDO
   ELSE
      IF(NORM)THEN
         call dcopy(nPrim*nCont,BASISINFO%ATOMTYPE(type)&
              &%SHELL(B)%SEGMENT(J)%elms,1,&
              &AO%CC(UNIQUEMATRIXES)%elms,1)
      ELSE
         call dcopy(nPrim*nCont,BASISINFO%ATOMTYPE(type)&
              &%SHELL(B)%SEGMENT(J)%UCCelms,1,&
              &AO%CC(UNIQUEMATRIXES)%elms,1)
      ENDIF
   ENDIF
   AO%angmom(UNIQUEMATRIXES) = B-1
   IF(IPRINT .GT. 2)THEN
      WRITE(LUPRI,*)'THE NEW CC',UNIQUEMATRIXES
      call lsmat_dense_print(AO%CC(UNIQUEMATRIXES),&
           &1,nPrim,1,nCont,lupri)
   ENDIF
ENDIF

END SUBROuTINE ADD_CC

!> \brief direct contraction coefficient matrix pointers to the contraction coefficients in the AOorganize
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DIRECT_CC_POINTERS(AO,AOmodel,AOorganize,nCont,nPrim,A,nbat,B,&
     &lupri,iprint,orbitalIndex,primorbitalIndex,UNCONTRACTED,INTNRM)
IMPLICIT NONE
INTEGER,intent(in)        :: A
TYPE(AOITEM)              :: AO,AOmodel
TYPE(AOorganizer)         :: AOorganize 
INTEGER                   :: nOrbitals,orbitalIndex,lupri,iprint,B
INTEGER                   :: nCont,nPrim,nbat,I,primorbitalindex
LOGICAL                   :: UNCONTRACTED,INTNRM
IF(UNCONTRACTED)THEN
   DO I=0,nPrim-1
      AOmodel%BATCH(nbat+I)%pCC(A)%p=>&
           &AO%CC(AOorganize%ORG(nbat+I)%CC(A))
      AOmodel%BATCH(nbat+I)%CCindex(A)=AOorganize%ORG(nbat+I)%CC(A)
      AOmodel%BATCH(nbat+I)%Expindex=AOorganize%ORG(nbat+I)%exponentindex

      IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to        ',AOorganize%ORG(nbat)%CC(A)
      IF(IPRINT .GT. 5)THEN
         WRITE(LUPRI,*)'Exp already points to '
         CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbat+I)%pExponents,1,1,1,1,lupri)
      ENDIF
      AOmodel%BATCH(nbat+I)%nContracted(A) = 1
      AOmodel%BATCH(nbat+I)%maxContracted= 1
      AOmodel%BATCH(nbat+I)%startOrbital(A) = orbitalIndex
      AOmodel%BATCH(nbat+I)%startprimOrbital(A) = primorbitalIndex
      AOmodel%BATCH(nbat+I)%nPrimOrbComp(A) = B*(B+1)/2
      IF (AOmodel%BATCH(nbat+I)%spherical) THEN
         AOmodel%BATCH(nbat+I)%nOrbComp(A) = 2*B-1
      ELSE
         AOmodel%BATCH(nbat+I)%nOrbComp(A) = B*(B+1)/2
      ENDIF
      nOrbitals = AOmodel%BATCH(nbat+I)%nOrbComp(A)
      orbitalIndex = orbitalIndex + nOrbitals
      primorbitalIndex = primorbitalIndex + nOrbitals
   ENDDO
ELSE !DEFAULT
   AOmodel%BATCH(nbat)%pCC(A)%p=>&
        &AO%CC(AOorganize%ORG(nbat)%CC(A))
   AOmodel%BATCH(nbat)%CCindex(A)=AOorganize%ORG(nbat)%CC(A)
   AOmodel%BATCH(nbat)%Expindex=AOorganize%ORG(nbat)%exponentindex

   IF(IPRINT .GT. 2) WRITE(LUPRI,*)'CC pointer to        ',AOorganize%ORG(nbat)%CC(A)
   IF(IPRINT .GT. 5)THEN
      WRITE(LUPRI,*)'Exp already points to '
      CALL LSMAT_DENSE_PRINT(AOmodel%BATCH(nbat)%pExponents,1,&
           &AOmodel%BATCH(nbat)%pExponents%nrow,1,&
           &AOmodel%BATCH(nbat)%pExponents%ncol,lupri)
   ENDIF
   IF(INTNRM)THEN 
      AOmodel%BATCH(nbat)%nContracted(A) = nprim
      AOmodel%BATCH(nbat)%maxContracted=MAX(AOmodel%BATCH(nbat)%maxContracted,nprim)
   ELSE
      AOmodel%BATCH(nbat)%nContracted(A) = nCont
      AOmodel%BATCH(nbat)%maxContracted=MAX(AOmodel%BATCH(nbat)%maxContracted,nCont)
   ENDIF
   AOmodel%BATCH(nbat)%startOrbital(A) = orbitalIndex
   AOmodel%BATCH(nbat)%startprimOrbital(A) = primorbitalIndex
   AOmodel%BATCH(nbat)%nPrimOrbComp(A) = B*(B+1)/2
   IF (AOmodel%BATCH(nbat)%spherical) THEN
      AOmodel%BATCH(nbat)%nOrbComp(A) = 2*B-1
   ELSE
      AOmodel%BATCH(nbat)%nOrbComp(A) = B*(B+1)/2
   ENDIF
   nOrbitals = AOmodel%BATCH(nbat)%nOrbComp(A)*nCont
   orbitalIndex = orbitalIndex + nOrbitals
   primorbitalIndex = primorbitalIndex&
        & + AOmodel%BATCH(nbat)%nOrbComp(A)*nprim
   AOmodel%BATCH(nbat)%nOrbitals(A)    = nOrbitals      
ENDIF

END SUBROUTINE DIRECT_CC_POINTERS

!> \brief in case of old atoms we copy AObatches from a modelAO to the AO
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE COPY_FROM_MODEL_AO(AOmodel,AO,iATOM,MOLECULE,orbitalIndex,primorbitalindex,lupri,offsetAT)
use molecule_type
IMPLICIT NONE
TYPE(MOLECULEINFO)        :: MOLECULE
TYPE(AOITEM)              :: AO,AOmodel
INTEGER                   :: I,J,bat,iATOM,orbitalIndex,lupri,SUM
INTEGER                   :: primorbitalindex,pSUM,norbitals,offsetAT

SUM=0
pSUM=0
AO%ATOMICnBatch(offsetAT+IATOM) = AOmodel%nbatches 
DO I=1,AOmodel%nbatches
   bat=AO%nbatches+1
   AO%nbatches = bat
   AO%BATCH(bat)%atom = offsetAT+IATOM
   AO%BATCH(bat)%molecularIndex = MOLECULE%ATOM(iATOM)%molecularIndex
   AO%BATCH(bat)%batch = I
   AO%BATCH(bat)%itype = AOmodel%BATCH(I)%itype
   AO%BATCH(bat)%type_empty = AOmodel%BATCH(I)%type_empty
   AO%BATCH(bat)%type_nucleus = AOmodel%BATCH(I)%type_nucleus
   AO%BATCH(bat)%type_elField = AOmodel%BATCH(I)%type_elField
   AO%BATCH(bat)%type_pCharge = AOmodel%BATCH(I)%type_pCharge
   AO%BATCH(bat)%spherical = AOmodel%BATCH(I)%spherical
   AO%BATCH(bat)%CENTER(1) = MOLECULE%ATOM(iATOM)%CENTER(1)
   AO%BATCH(bat)%CENTER(2) = MOLECULE%ATOM(iATOM)%CENTER(2)
   AO%BATCH(bat)%CENTER(3) = MOLECULE%ATOM(iATOM)%CENTER(3)
   AO%BATCH(bat)%nPrimitives = AOmodel%BATCH(I)%nPrimitives
   AO%BATCH(bat)%maxContracted = AOmodel%BATCH(I)%maxContracted
   AO%BATCH(bat)%maxAngmom = AOmodel%BATCH(I)%maxAngmom
   AO%BATCH(bat)%pExponents => AOmodel%BATCH(I)%pExponents
   AO%BATCH(bat)%nAngmom = AOmodel%BATCH(I)%nAngmom
   AO%BATCH(bat)%extent  = AOmodel%BATCH(I)%extent
   AO%BATCH(bat)%Expindex= AOmodel%BATCH(I)%ExpIndex
   DO J=1,AOmodel%BATCH(I)%nAngmom
      AO%BATCH(bat)%Angmom(J) = AOmodel%BATCH(I)%Angmom(J)
      AO%BATCH(bat)%nContracted(J) = AOmodel%BATCH(I)%nContracted(J)
      AO%BATCH(bat)%nOrbComp(J) = AOmodel%BATCH(I)%nOrbComp(J)
      AO%BATCH(bat)%nPrimOrbComp(J) = AOmodel%BATCH(I)%nPrimOrbComp(J)
      AO%BATCH(bat)%pCC(J)%p =>  AOmodel%BATCH(I)%pCC(J)%p
      AO%BATCH(bat)%CCindex(J) =  AOmodel%BATCH(I)%CCindex(J)

      AO%BATCH(bat)%startOrbital(J) = orbitalIndex &
           &+ AOmodel%BATCH(I)%startOrbital(J)  
      AO%BATCH(bat)%startprimOrbital(J) = primorbitalIndex &
           &+ AOmodel%BATCH(I)%startprimOrbital(J)  

      nOrbitals = AOmodel%BATCH(I)%nOrbComp(J)*AOmodel%BATCH(I)%nContracted(J)
      AO%BATCH(bat)%nOrbitals(J) = nOrbitals
      SUM=SUM+nOrbitals
      pSUM=pSUM+AOmodel%BATCH(I)%nOrbComp(J)*AOmodel%BATCH(I)%nPrimitives
   ENDDO
ENDDO
orbitalIndex = orbitalIndex + SUM
AO%ATOMICnORB(offsetAT+IATOM) = SUM
primorbitalIndex = primorbitalIndex + pSUM

END SUBROUTINE COPY_FROM_MODEL_AO

!> \brief build basinf structure used in exchange-correlation calculations
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,GRDONE,dopriexp)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER           :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER           :: IPRINT
!> the basinf structure to be build
TYPE(BASINF)      :: BAS
!> logical for grid. 1 if grid have been build and 0 if not 
INTEGER           :: GRDONE 
!> the setting structure containing info about molecule,basisset,...
TYPE(LSSETTING)   :: SETTING
!> should the priexpM be constructed
logical        :: dopriexp 
!
INTEGER           :: MAXNSHELL !number of shells see type def of BASINF
INTEGER           :: MAXANGMOM
INTEGER           :: nprim,TOTPRIM,nrow,ncol,N,OLDnrow
INTEGER           :: TYP,NR,USHELL,nshells,ushell2
INTEGER           :: natoms,type,SHELL,orbitalindex,L,I,K,norb,R,M,MXPRIM
INTEGER,pointer :: ATOMtype(:)
INTEGER,pointer :: MODELTYPES(:)
LOGICAL           :: NEWATOM
REAL(REALK)       :: FACL(10),R2,THLOG,EXP,LOGVAL
REAL(REALK),PARAMETER     :: DMIN=1E-13_realk
INTEGER           :: irow,nrow2,nrow3,NRSIZE,ICHARGE,J,nOrbComp
LOGICAL           :: END
INTEGER,pointer   :: NSHELLINDEX(:),USHELLINDEX(:)
INTEGER,pointer   :: UATOM(:)
TYPE(BASISSETINFO),pointer :: BASIS

IF(AORdefault.EQ.AOregular)THEN
   BASIS => SETTING%BASIS(1)%p%BINFO(RegBasParam)
ELSEIF(AORdefault.EQ.AOVAL)THEN
   BASIS => SETTING%BASIS(1)%p%BINFO(VALBasParam)
ELSEIF(AORdefault.EQ.AOdfAux)THEN
   BASIS => SETTING%BASIS(1)%p%BINFO(AuxBasParam)
ELSEIF(AORdefault.EQ.AOdfJK)THEN
   BASIS => SETTING%BASIS(1)%p%BINFO(JKBasParam)
ELSEIF(AORdefault.EQ.AOdfCABS)THEN
   BASIS => SETTING%BASIS(1)%p%BINFO(CABBasParam)
ELSEIF(AORdefault.EQ.AOadmm)THEN
   BASIS => SETTING%BASIS(1)%p%BINFO(ADMBasParam)
ELSE
   CALL LSQUIT('ERROR in BASINF unknown AO',-1)
ENDIF


!natoms = SETTING%MOLECULE(1)%p%nAtoms
natoms = 0
DO I=1,SETTING%MOLECULE(1)%p%nAtoms
   IF(.NOT.SETTING%MOLECULE(1)%p%ATOM(I)%pointcharge)THEN
      natoms = natoms+1
   ENDIF
ENDDO

call mem_alloc(ATOMtype,natoms)
BAS%nAtoms = nAtoms
CALL MEM_ALLOC(BAS%X,nAtoms)
CALL MEM_ALLOC(BAS%Y,nAtoms)
CALL MEM_ALLOC(BAS%Z,nAtoms)
CALL MEM_ALLOC(BAS%CHARGE,nAtoms)

I=0
DO J=1,SETTING%MOLECULE(1)%p%nAtoms
   IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
   I=I+1
   BAS%X(I) = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(1)
   BAS%Y(I) = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(2)
   BAS%Z(I) = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(3)
   BAS%CHARGE(I) = INT(SETTING%MOLECULE(1)%p%ATOM(J)%CHARGE)
   IF(SETTING%MOLECULE(1)%p%ATOM(J)%phantom)BAS%CHARGE(I) = 0
ENDDO

MAXNSHELL=0
MXPRIM=0
BAS%spherical = SETTING%BASIS(1)%p%BINFO(REGBASPARAM)%spherical
R = SETTING%BASIS(1)%p%BINFO(REGBASPARAM)%Labelindex
IF(R.EQ. 0)THEN
   I=0
   DO J=1,SETTING%MOLECULE(1)%p%nAtoms
      IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
      I=I+1
      ICHARGE = INT(SETTING%MOLECULE(1)%p%ATOM(J)%Charge)
      type = SETTING%BASIS(1)%p%BINFO(REGBASPARAM)%Chargeindex(ICHARGE)
      ATOMtype(I)=type
      !determine MAXNSHELL = number of shells 
      DO K=1,BASIS%ATOMTYPE(type)%nAngmom
         MAXNSHELL=MAXNSHELL+BASIS%ATOMTYPE(type)%SHELL(K)%norb
         MXPRIM=MXPRIM+BASIS%ATOMTYPE(type)%SHELL(K)%nprim
      ENDDO
   ENDDO
ELSE
   I=0
   DO J=1,SETTING%MOLECULE(1)%p%nAtoms
      IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
      I=I+1
      ICHARGE = INT(SETTING%MOLECULE(1)%p%ATOM(J)%Charge)
      type = SETTING%MOLECULE(1)%p%ATOM(J)%IDtype(R)
      ATOMtype(I)=type
      !determine MAXNSHELL = number of shells 
      DO K=1,BASIS%ATOMTYPE(type)%nAngmom
         MAXNSHELL=MAXNSHELL+BASIS%ATOMTYPE(type)%SHELL(K)%norb
         MXPRIM=MXPRIM+BASIS%ATOMTYPE(type)%SHELL(K)%nprim
      ENDDO
   ENDDO
ENDIF
BAS%MAXNSHELL=MAXNSHELL

CALL MEM_ALLOC(BAS%SHELL2ATOM,MAXNSHELL)
CALL MEM_ALLOC(BAS%NSTART,MAXNSHELL)
CALL MEM_ALLOC(BAS%SHELLANGMOM,MAXNSHELL)
CALL MEM_ALLOC(BAS%SHELLNPRIM,MAXNSHELL)
CALL MEM_ALLOC(BAS%PRIEXPSTART,MAXNSHELL)
CALL MEM_ALLOC(BAS%CENT,3,MAXNSHELL)

orbitalindex = 0
SHELL = 1
MAXANGMOM=0
TOTPRIM=0
I=0
DO J=1,SETTING%MOLECULE(1)%p%nAtoms
   IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
   I=I+1
   type=ATOMtype(I)
   DO K=1,BASIS%ATOMTYPE(type)%nAngmom
      norb=BASIS%ATOMTYPE(type)%SHELL(K)%norb
      nprim=BASIS%ATOMTYPE(type)%SHELL(K)%nprim
      DO L=1,norb
         BAS%SHELL2ATOM(SHELL) = I
         BAS%CENT(1,SHELL) = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(1) 
         BAS%CENT(2,SHELL) = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(2) 
         BAS%CENT(3,SHELL) = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(3) 
         BAS%NSTART(SHELL) = orbitalindex
         BAS%PRIEXPSTART(SHELL) = TOTPRIM !the accumulated number of primitives
         BAS%SHELLANGMOM(SHELL) = K
         BAS%SHELLNPRIM(SHELL) = nprim
         IF(BAS%spherical)THEN
            nOrbComp = 2*K-1
         ELSE
            nOrbComp = K*(K+1)/2
         ENDIF
         orbitalindex=orbitalindex+nOrbComp
         SHELL = SHELL+1
      ENDDO
      TOTPRIM=TOTPRIM+nprim
      MAXANGMOM = MAX(K,MAXANGMOM)
   ENDDO
ENDDO

CALL MEM_ALLOC(BAS%PRIEXP,MXPRIM)
MXPRIM = 1
I=0
DO J=1,SETTING%MOLECULE(1)%p%nAtoms
   IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
   I=I+1
   type=ATOMtype(I)
   DO K=1,BASIS%ATOMTYPE(type)%nAngmom
      DO L=1,BASIS%ATOMTYPE(type)%SHELL(K)%nsegments
         DO M= 1,BASIS%ATOMTYPE(type)%SHELL(K)%segment(L)%nrow
            BAS%PRIEXP(MXPRIM) = BASIS%ATOMTYPE(type)%SHELL(K)%segment(L)%Exponents(M)
            MXPRIM = MXPRIM+1
         ENDDO
      ENDDO
   ENDDO
ENDDO
BAS%GRDONE = 0
BAS%MAXANGMOM=MAXANGMOM
BAS%MXPRIM=MXPRIM-1
MXPRIM = MXPRIM-1

!-------------------------------------------------------
!THIS SHOULD NOT BE CONSTUCTED EVERY TIME - JUST SAVE CC,CCINDEX,CCSTART in DALTON
!------------------------------------------------------------
TYP=BASIS%nAtomtypes
call mem_alloc(MODELTYPES,TYP)
NR=0
DO I=1,BASIS%natomtypes    
   NR=NR+1
   MODELTYPES(I)=NR
ENDDO
NRSIZE=NR
call mem_alloc(UATOM,NR)
call mem_alloc(NSHELLINDEX,NR)
call mem_alloc(USHELLINDEX,NR)
DO I=1,NR
   UATOM(I) = 0 
ENDDO

USHELL=1
I=0
DO J=1,SETTING%MOLECULE(1)%p%nAtoms
   IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
   I=I+1
   type=ATOMtype(I)  
   ! TEST IF THIS TYPE OF ATOM HAS ALREADY BEEN PROCESSED -> OLDATOM=TRUE
   NEWATOM=.FALSE.
   NR=MODELTYPES(type)
   IF(UATOM(NR) .EQ. 0) NEWATOM=.TRUE.
   IF(NEWATOM) THEN
     UATOM(NR)=1
     DO K=1,BASIS%ATOMTYPE(type)%nAngmom
      DO L=1,BASIS%ATOMTYPE(type)%SHELL(K)%nsegments
         ncol=BASIS%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%ncol
         DO N=1,ncol
            USHELL = USHELL+1
         ENDDO
      ENDDO
     ENDDO
   ENDIF
ENDDO

DO I=1,NRSIZE
   UATOM(I) = 0 
ENDDO

BAS%Ushells=USHELL-1
CALL MEM_ALLOC(BAS%CC,USHELL)
CALL MEM_ALLOC(BAS%CCSTART,USHELL)
CALL MEM_ALLOC(BAS%CCINDEX,MAXNSHELL)
IF(dopriexp)THEN
   CALL MEM_ALLOC(BAS%priexpM,USHELL)
ELSE
   nullify(BAS%priexpM)
ENDIF

SHELL=1
USHELL=1
I=0
DO J=1,SETTING%MOLECULE(1)%p%nAtoms
   IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
   I=I+1
   type=ATOMtype(I)  
   ! TEST IF THIS TYPE OF ATOM HAS ALREADY BEEN PROCESSED -> OLDATOM=TRUE
   NEWATOM=.FALSE.
   NR=MODELTYPES(type)
   IF(UATOM(NR) .EQ. 0) NEWATOM=.TRUE.
   IF(NEWATOM) THEN
      UATOM(NR)=I
      USHELLINDEX(NR)=USHELL
      NSHELLINDEX(NR)=0
     DO K=1,BASIS%ATOMTYPE(type)%nAngmom
         OLDnrow=1
      DO L=1,BASIS%ATOMTYPE(type)%SHELL(K)%nsegments
         nrow=BASIS%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%nrow
         ncol=BASIS%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%ncol
         DO N=1,ncol
            nrow2=0
            DO irow = 1,nrow
               IF(ABS(BASIS%ATOMTYPE(type)%&
                    &SHELL(K)%SEGMENT(L)%elms(irow+(N-1)*nrow)) > DMIN)THEN
                  nrow2=nrow2+1
               ENDIF
            ENDDO
            if(nrow2 .EQ. nrow)then !all element nonzero!      
               call lsmat_dense_init(BAS%CC(USHELL),nrow,1)
               call dcopy(nrow,BASIS%ATOMTYPE(type)%&
                    &SHELL(K)%SEGMENT(L)%elms(1+(N-1)*nrow:N*nrow),&
                    &1,BAS%CC(USHELL)%elms,1)
               IF(dopriexp)then
                  call lsmat_dense_init(BAS%priexpM(USHELL),nrow,1)
                  call dcopy(nrow,BASIS%ATOMTYPE(type)%&
                       &SHELL(K)%SEGMENT(L)%Exponents(1:nrow),&
                       &1,BAS%priexpM(USHELL)%elms,1)
               endif
               BAS%CCSTART(USHELL)=OLDnrow
               BAS%CCINDEX(SHELL)=USHELL
               SHELL = SHELL+1
               USHELL = USHELL+1
               NSHELLINDEX(NR)=NSHELLINDEX(NR)+1
            ELSE !some zeros
               nrow3=nrow
               END=.FALSE.
               DO irow = 1,nrow
                  IF(END)EXIT
                  IF(ABS(BASIS%ATOMTYPE(type)%&
                       &SHELL(K)%SEGMENT(L)%elms(irow+(N-1)*nrow)) < DMIN)THEN
                     nrow3=nrow3-1
                  ELSE
                     call lsmat_dense_init(BAS%CC(USHELL),nrow3,1)
                     call dcopy(nrow3,BASIS%ATOMTYPE(type)%&
                          &SHELL(K)%SEGMENT(L)%elms(irow+(N-1)*nrow:N*nrow),&
                          &1,BAS%CC(USHELL)%elms,1)
                     IF(dopriexp)then
                        call lsmat_dense_init(BAS%priexpM(USHELL),nrow3,1)
                        call dcopy(nrow3,BASIS%ATOMTYPE(type)%&
                             &SHELL(K)%SEGMENT(L)%Exponents(irow:nrow),&
                             &1,BAS%priexpM(USHELL)%elms,1)
                     endif
                     BAS%CCSTART(USHELL)=OLDnrow+irow-1
                     END=.TRUE.
                  ENDIF
               ENDDO
               BAS%CCINDEX(SHELL)=USHELL
               SHELL = SHELL+1
               USHELL = USHELL+1
               NSHELLINDEX(NR)=NSHELLINDEX(NR)+1
            ENDIF
         ENDDO
         OLDnrow=nrow+OLDnrow
      ENDDO
     ENDDO
   ELSE ! OLD ATOM TYPE -> COPY FROM Uniqueatom
      NR=MODELTYPES(type)
      nshells = NSHELLINDEX(NR)
      USHELL2 = USHELLINDEX(NR) 
      DO L=1,nSHELLS
         BAS%CCINDEX(SHELL)=USHELL2
         SHELL = SHELL+1
         USHELL2 = USHELL2+1
      ENDDO
   ENDIF
ENDDO

call mem_dealloc(NSHELLINDEX)
call mem_dealloc(USHELLINDEX)
call mem_dealloc(UATOM)

!---------------------------------------------------------------------------

IF(SHELL-1 .NE. MAXNSHELL) CALL LSQUIT('number of shell is wrong in Build_basinf',lupri) 


IF(GRDONE .EQ. 0)THEN
   ! get radii of all shells as defined by specified threshold DFTHRI.
   FACL(1)= 1E0_realk
   FACL(2)= 1.3333E0_realk
   FACL(3)= 1.6E0_realk
   FACL(4)= 1.83E0_realk
   FACL(5)= 2.03E0_realk
   FACL(6)= 2.22E0_realk
   FACL(7)= 2.39E0_realk
   FACL(8)= 2.55E0_realk
   FACL(9)= 2.70E0_realk
   FACL(10)= 2.84E0_realk

   CALL MEM_ALLOC(BAS%RSHEL,MAXNSHELL)
   CALL LS_DZERO(BAS%RSHEL,MAXNSHELL)
   THLOG = LOG(SETTING%SCHEME%DFT%DFTHRI)
   SHELL = 1
   I=0
   DO J=1,SETTING%MOLECULE(1)%p%nAtoms
      IF(SETTING%MOLECULE(1)%p%ATOM(J)%pointcharge)CYCLE
      I=I+1
      type=ATOMtype(I)
      DO K=1,BASIS%ATOMTYPE(type)%nAngmom
         DO L=1,BASIS%ATOMTYPE(type)%SHELL(K)%nsegments
            nrow=BASIS%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%nrow
            ncol=BASIS%ATOMTYPE(type)%SHELL(K)%SEGMENT(L)%ncol
            DO N=1,ncol
               DO M=1,nrow
                  EXP=BASIS%ATOMTYPE(type)%SHELL(K)%&
                       &SEGMENT(L)%Exponents(M)
                  LOGVAL=ABS(BASIS%ATOMTYPE(type)%&
                       &SHELL(K)%SEGMENT(L)%elms(M+(N-1)*nrow) )  
                  IF(LOGVAL.GT.0E0_realk)THEN
                     R2 = (LOG(LOGVAL) -THLOG)/EXP
                     IF(BAS%RSHEL(SHELL).LT.R2) BAS%RSHEL(SHELL) = R2
                  ENDIF
               ENDDO
               SHELL = SHELL+1
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   DO SHELL=1,MAXNSHELL
      BAS%RSHEL(SHELL) = SQRT(BAS%RSHEL(SHELL))*FACL(BAS%SHELLANGMOM(SHELL))
   END DO
ELSE
   !not used but call mem_allocd to satisfy some compilers
   CALL MEM_ALLOC(BAS%RSHEL,1)
ENDIF

call mem_dealloc(ATOMtype)
call mem_dealloc(MODELTYPES)

END SUBROUTINE BUILD_BASINF

!> \brief free the basinf structure used in exchange-correlation calculations
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE FREE_BASINF(BAS)
TYPE(BASINF)      :: BAS
INTEGER           :: ushells,I
ushells = BAS%ushells

CALL MEM_DEALLOC(BAS%X)
CALL MEM_DEALLOC(BAS%Y)
CALL MEM_DEALLOC(BAS%Z)
CALL MEM_DEALLOC(BAS%CHARGE)
CALL MEM_DEALLOC(BAS%SHELL2ATOM)
CALL MEM_DEALLOC(BAS%NSTART)
CALL MEM_DEALLOC(BAS%SHELLANGMOM)
CALL MEM_DEALLOC(BAS%SHELLNPRIM)
CALL MEM_DEALLOC(BAS%PRIEXPSTART)
CALL MEM_DEALLOC(BAS%CENT)
CALL MEM_DEALLOC(BAS%PRIEXP)
DO I=1,ushells
   call lsmat_dense_free(BAS%CC(I))
ENDDO   
CALL MEM_DEALLOC(BAS%CC)
CALL MEM_DEALLOC(BAS%CCSTART)
CALL MEM_DEALLOC(BAS%CCINDEX)
CALL MEM_DEALLOC(BAS%RSHEL)
IF(associated(BAS%priexpM))THEN
   CALL MEM_DEALLOC(BAS%priexpM)
ENDIF

END SUBROUTINE FREE_BASINF

SUBROUTINE ADD_TO_CC(CC,STARTPRIM,NROW,TMP,N,MXPRIM,MAXNSHELL)
IMPLICIT NONE
INTEGER     :: M,STARTPRIM,NROW,N,MAXNSHELL,MXPRIM
REAL(REALK) :: CC(MXPRIM*MAXNSHELL),TMP(:)

DO M=1,nrow
   CC(STARTPRIM+M) = TMP(M+(N-1)*nrow)  
ENDDO

END SUBROUTINE ADD_TO_CC

SUBROUTINE SHELL_NUMBER(LUPRI,IPRINT,DALTON,BASISINFO,MAXNSHELL)
IMPLICIT NONE
INTEGER           :: LUPRI,IPRINT
TYPE(DALTONINPUT) :: DALTON
TYPE(BASISSETINFO)        :: BASISINFO
INTEGER           :: MAXNSHELL !number of shells see type def of BASINF
!
INTEGER           :: I,R,type,K,icharge

MAXNSHELL=0
R = BASISINFO%Labelindex

IF(R.EQ. 0)THEN
   DO I=1,DALTON%MOLECULE%nAtoms
      ICHARGE = INT(DALTON%MOLECULE%ATOM(I)%Charge)
      type = BASISINFO%Chargeindex(ICHARGE)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         MAXNSHELL=MAXNSHELL+BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
      ENDDO
   ENDDO
ELSE
   DO I=1,DALTON%MOLECULE%nAtoms
      type = DALTON%MOLECULE%ATOM(I)%IDtype(R)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         MAXNSHELL=MAXNSHELL+BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE SHELL_NUMBER

SUBROUTINE ORBITAL_DEGENERACY(LUPRI,IPRINT,DALTON,BASISINFO,MAXNSHELL,DEGENERACY)
IMPLICIT NONE
INTEGER           :: LUPRI,IPRINT,MAXNSHELL!number of shells
TYPE(DALTONINPUT) :: DALTON
TYPE(BASISSETINFO)        :: BASISINFO
INTEGER           :: DEGENERACY(MAXNSHELL) 
!
INTEGER           :: I,R,type,K,shell,norb,L,ICHARGE
R = BASISINFO%Labelindex

IF(R.EQ. 0)THEN
   SHELL = 1
   DO I=1,DALTON%MOLECULE%nAtoms
      ICHARGE = INT(DALTON%MOLECULE%ATOM(I)%Charge)
      type = BASISINFO%Chargeindex(ICHARGE)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         DO L=1,norb
            DEGENERACY(SHELL)=(2*(K-1)+1)
            SHELL = SHELL+1
         ENDDO
      ENDDO
   ENDDO
ELSE
   SHELL = 1
   DO I=1,DALTON%MOLECULE%nAtoms
      type = DALTON%MOLECULE%ATOM(I)%IDtype(R)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         DO L=1,norb
            DEGENERACY(SHELL)=(2*(K-1)+1)
            SHELL = SHELL+1
         ENDDO
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE ORBITAL_DEGENERACY

SUBROUTINE ATOMIC_ORBITALS(LUPRI,IPRINT,DALTON,BASISINFO,natoms,orbitals)
IMPLICIT NONE
INTEGER           :: LUPRI,IPRINT,natoms!number of shells
TYPE(DALTONINPUT) :: DALTON
TYPE(BASISSETINFO)        :: BASISINFO
INTEGER           :: orbitals(natoms) 
!
INTEGER           :: I,R,type,K,norb,ICHARGE

CALL LS_IZERO(ORBITALS,natoms)

R = BASISINFO%Labelindex
IF(R.EQ. 0)THEN
   DO I=1,DALTON%MOLECULE%nAtoms
      ICHARGE = INT(DALTON%MOLECULE%ATOM(I)%Charge)
      type = BASISINFO%Chargeindex(ICHARGE)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         ORBITALS(I)=ORBITALS(I)+(2*(K-1)+1)*norb
      ENDDO
   ENDDO
ELSE
   DO I=1,DALTON%MOLECULE%nAtoms
      type = DALTON%MOLECULE%ATOM(I)%IDtype(R)
      DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
         norb=BASISINFO%ATOMTYPE(type)%SHELL(K)%norb
         ORBITALS(I)=ORBITALS(I)+(2*(K-1)+1)*norb
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE ATOMIC_ORBITALS

END MODULE BUILDAOBATCH


