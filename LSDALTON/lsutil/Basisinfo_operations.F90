!> @file
!> basisinfo type module, contains also standard operation subroutines for this type
!> \brief basisinfo type module and associated subroutines for this type 
!> \author T. Kjaergaard
!> \date 2010
MODULE basis_type
 use lsmatrix_operations_dense
 use precision
 use memory_handling
 use AO_type
 use basis_typetype
 private
 public :: write_atomtypeitem_to_disk,& 
      & write_basissetinfo_to_disk, &
      & transform_basis_to_GCAO, &
      & write_basisinfo_to_disk, &
      & read_atomtypeitem_from_disk, &
      & read_basissetinfo_from_disk, &
      & read_basisinfo_from_disk, &
      & init_basissetinfo_types, &
      & init_basissetinfo_ContractionM, &
      & init_basissetinfo_elms, &
      & alloc_and_take_subbasissetinfo, &
      & copy_basissetinfo, &
      & free_basissetinfo, &
      & lsmpi_alloc_basissetinfo, &
      & print_basissetinfo, &
      & freeBrakebasinfo,initBrakebasinfo,&
      & copybrakebasinfo,buildbasisfrombrakebasinf,&
      & print_brakebas, add_basissetinfo

CONTAINS
!#################################################################
!#
!# SUBROUTINES THAT ONLY AFFECT THE TYPES DEFINED IN THIS FILE
!# transform_basis_to_GCAO
!# copy_basissetinfo
!# write_atomtypeitem_to_disk
!# write_basissetinfo_to_disk
!# write_basisinfo_to_disk
!# read_atomtypeitem_from_disk
!# read_basissetinfo_from_disk
!# read_basisinfo_from_disk
!# INIT_BASISSETINFO_Types
!# INIT_BASISSETINFO_ContractionMatrix
!# INIT_BASISSETINFO_elms
!# ALLOC_AND_TAKE_SUBBASISSETINFO
!# FREE_BASISSETINFO
!# LSMPI_ALLOC_BASISSETINFO
!# PRINT_BASISSETINFO
!################################################################

!> \brief write the atomtype structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param Atomtype Contains the basis information for the atomtype, and type specifiers
!> \param lun logic unit number of file to write to
SUBROUTINE write_atomtypeitem_to_disk(lun,ATOMTYPE)
implicit none
TYPE(ATOMTYPEITEM),intent(in)    :: ATOMTYPE
INTEGER,intent(in)               :: lun
!
Integer :: I,J,K,nrow,ncol

write(lun)ATOMTYPE%nAngmom
write(lun)ATOMTYPE%FAMILY
write(lun)ATOMTYPE%ToTnorb
write(lun)ATOMTYPE%ToTnprim
write(lun)ATOMTYPE%Charge
DO I=1,ATOMTYPE%nAngmom
   WRITE(lun)ATOMTYPE%SHELL(I)%nprim
   WRITE(lun)ATOMTYPE%SHELL(I)%norb
   WRITE(lun)ATOMTYPE%SHELL(I)%nsegments
   DO J=1,ATOMTYPE%SHELL(I)%nsegments
      nrow = ATOMTYPE%SHELL(I)%segment(J)%nrow
      ncol = ATOMTYPE%SHELL(I)%segment(J)%ncol
      WRITE(lun)nrow
      WRITE(lun)ncol
      WRITE(lun)(ATOMTYPE%SHELL(I)%segment(J)%elms(K),K=1,nrow*ncol)
      WRITE(lun)(ATOMTYPE%SHELL(I)%segment(J)%UCCelms(K),K=1,nrow*ncol)
      WRITE(lun)(ATOMTYPE%SHELL(I)%segment(J)%Exponents(K),K=1,nrow)
   ENDDO
ENDDO
write(lun)ATOMTYPE%Name

END SUBROUTINE write_atomtypeitem_to_disk

!> \brief write the basissetinfo structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basissetinfo Contains the basis information for given set
!> \param lun logic unit number of file to write to
SUBROUTINE write_basissetinfo_to_disk(lun,BASISSET)
implicit none
TYPE(BASISSETINFO),intent(in)    :: BASISSET
INTEGER,intent(in)               :: lun
!
INTEGER :: I

write(lun)BASISSET%natomtypes
write(lun)BASISSET%DunningsBasis
write(lun)BASISSET%GeminalScalingFactor
write(lun)BASISSET%GCbasis
write(lun)BASISSET%Spherical
write(lun)BASISSET%GCont
DO I=1,BASISSET%natomtypes
   call write_atomtypeitem_to_disk(lun,BASISSET%ATOMTYPE(I))
ENDDO
write(lun)BASISSET%labelindex
write(lun)BASISSET%nChargeindex
DO I=0,BASISSET%nChargeindex
   write(lun)BASISSET%Chargeindex(I)
ENDDO
write(lun)BASISSET%nbast
write(lun)BASISSET%nprimbast
write(lun)BASISSET%label

END SUBROUTINE write_basissetinfo_to_disk

subroutine transform_basis_to_GCAO(BASIS)
implicit none
TYPE(BASISINFO),intent(in)    :: BASIS
real(realk),pointer :: CCtmp(:,:),bCMO(:,:)
integer :: J,K,iprim,nrow2,ncol2
integer :: iCont,iseg,iContLoc,ic1,ip1,nprim,norb,ielm,iprimLoc
type(SHELL),pointer :: shell2
DO J=1,BASIS%BINFO(RegBasParam)%nAtomtypes
 DO K=1,BASIS%BINFO(RegBasParam)%ATOMTYPE(J)%nAngmom
  shell2 => BASIS%BINFO(RegBasParam)%ATOMTYPE(J)%SHELL(K)

  nrow2 = BASIS%BINFO(GCTBasParam)%ATOMTYPE(J)%SHELL(K)%nprim
  if (nrow2.eq. 0) cycle
  ncol2 = BASIS%BINFO(GCTBasParam)%ATOMTYPE(J)%SHELL(K)%norb
  call mem_alloc(bCMO,nrow2,ncol2)
  bCMO = reshape(BASIS%BINFO(GCTBasParam)%ATOMTYPE(J)%SHELL(K)%segment(1)%elms,(/ nrow2,ncol2 /))

  nprim     = shell2%nprim
  norb      = shell2%norb

  IF(norb.NE.nrow2)call lsquit('dim mismatch in transform_basis_to_GCAO',-1)

  call mem_alloc(CCtmp,nprim,norb)
  do iCont=1,norb
     do iprim=1,nprim
        CCtmp(iprim,iCont) = 0E0_realk
     enddo
  enddo
  iPrim=0
  iCont=0      
  DO iseg = 1,shell2%nsegments
     ielm = 0
     iContLoc = iCont
     do ic1=1,shell2%segment(iseg)%ncol
        iContLoc = iContLoc+1 
        iprimLoc = iPrim            
        do ip1=1,shell2%segment(iseg)%nrow
           iprimLoc = iprimLoc+1 
           ielm = ielm +1
           CCtmp(iPrimLoc,iContLoc) = shell2%segment(iseg)%elms(ielm)
        enddo
     enddo
     iPrim = iPrim + shell2%segment(iseg)%nrow
     iCont = iCont + shell2%segment(iseg)%ncol
  ENDDO  
  CCtmp = matmul(CCtmp,bCMO)   
  call mem_dealloc(shell2%segment(1)%elms)
  call mem_alloc(shell2%segment(1)%elms,nprim*norb)
  shell2%segment(1)%elms = reshape(CCtmp, (/ nprim*norb /))
  
  iPrim=0
  DO iseg = 1,shell2%nsegments
     do ip1=1,shell2%segment(iseg)%nrow
        iprim = iprim+1 
        CCtmp(iPrim,1) = shell2%segment(iseg)%exponents(iP1)
     enddo
  enddo
  call mem_dealloc(shell2%segment(1)%exponents)
  call mem_alloc(shell2%segment(1)%exponents,nprim)
  DO iPrim=1,nPrim
     shell2%segment(1)%exponents(iPrim) = CCtmp(iPrim,1)
  ENDDO
  call mem_dealloc(shell2%segment(1)%UCCelms)
  call mem_alloc(shell2%segment(1)%UCCelms,nprim*norb)
  shell2%segment(1)%UCCelms = 0E0_realk
  DO iseg = 2,shell2%nsegments
     call mem_dealloc(shell2%segment(iseg)%elms)
     call mem_dealloc(shell2%segment(iseg)%UCCelms)
     call mem_dealloc(shell2%segment(iseg)%exponents)
  enddo
  shell2%nsegments=1      
  shell2%segment(1)%ncol = norb
  shell2%segment(1)%nrow = nprim
  call mem_dealloc(CCtmp)
  call mem_dealloc(bCMO)
 ENDDO
ENDDO
end subroutine transform_basis_to_GCAO
!> \brief write the basisinfo structure to disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information
!> \param lun logic unit number of file to write to
SUBROUTINE write_basisinfo_to_disk(lun,BASIS)
implicit none
TYPE(BASISINFO),intent(in)    :: BASIS
INTEGER,intent(in)            :: lun
!
integer :: i
DO I=1,nBasisBasParam
   write(lun) BASIS%WBASIS(I)
   IF(BASIS%WBASIS(I))THEN
      CALL WRITE_BASISSETINFO_TO_DISK(lun,BASIS%BINFO(I))
   ENDIF
ENDDO
END SUBROUTINE write_basisinfo_to_disk

!> \brief read the atomtype structure from disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param Atomtype Contains the basis information for the atomtype, and type specifiers
!> \param lun logic unit number of file to write to
SUBROUTINE read_atomtypeitem_from_disk(lun,ATOMTYPE)
implicit none
TYPE(ATOMTYPEITEM),intent(inout)    :: ATOMTYPE
INTEGER,intent(in)               :: lun
!
Integer :: I,J,K,nrow,ncol

read(lun)ATOMTYPE%nAngmom
read(lun)ATOMTYPE%FAMILY
read(lun)ATOMTYPE%ToTnorb
read(lun)ATOMTYPE%ToTnprim
read(lun)ATOMTYPE%Charge
DO I=1,ATOMTYPE%nAngmom
   READ(lun)ATOMTYPE%SHELL(I)%nprim
   READ(lun)ATOMTYPE%SHELL(I)%norb
   READ(lun)ATOMTYPE%SHELL(I)%nsegments
   DO J=1,ATOMTYPE%SHELL(I)%nsegments
      READ(lun)nrow
      READ(lun)ncol
      ATOMTYPE%SHELL(I)%segment(J)%nrow = nrow
      ATOMTYPE%SHELL(I)%segment(J)%ncol = ncol
      CALL MEM_ALLOC(ATOMTYPE%SHELL(I)%segment(J)%elms,nrow*ncol)
      READ(lun)(ATOMTYPE%SHELL(I)%segment(J)%elms(K),K=1,nrow*ncol)
      CALL MEM_ALLOC(ATOMTYPE%SHELL(I)%segment(J)%UCCelms,nrow*ncol)
      READ(lun)(ATOMTYPE%SHELL(I)%segment(J)%UCCelms(K),K=1,nrow*ncol)
      CALL MEM_ALLOC(ATOMTYPE%SHELL(I)%segment(J)%Exponents,nrow)
      READ(lun)(ATOMTYPE%SHELL(I)%segment(J)%Exponents(K),K=1,nrow)
   ENDDO
ENDDO
read(lun)ATOMTYPE%Name

END SUBROUTINE read_atomtypeitem_from_disk

!> \brief read the basissetinfo structure from disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basissetinfo Contains the basis information for given set
!> \param lun logic unit number of file to read from
SUBROUTINE read_basissetinfo_from_disk(lun,BASISSET)
implicit none
TYPE(BASISSETINFO),intent(inout)    :: BASISSET
INTEGER,intent(in)               :: lun
!
INTEGER :: I
call nullifybasisset(BASISSET)
read(lun)BASISSET%natomtypes
IF(BASISSET%natomtypes.GT. 0)THEN
   read(lun)BASISSET%DunningsBasis
   read(lun)BASISSET%GeminalScalingFactor
   read(lun)BASISSET%GCbasis
   read(lun)BASISSET%Spherical
   read(lun)BASISSET%GCont
   CALL MEM_ALLOC(BASISSET%ATOMTYPE,BASISSET%natomtypes)
   DO I=1,BASISSET%natomtypes
      call read_atomtypeitem_from_disk(lun,BASISSET%ATOMTYPE(I))
   ENDDO
   read(lun)BASISSET%labelindex
   read(lun)BASISSET%nChargeindex
   CALL MEM_ALLOC(BASISSET%Chargeindex,BASISSET%nChargeindex,.TRUE.)
   DO I=0,BASISSET%nChargeindex
      read(lun)BASISSET%Chargeindex(I)
   ENDDO
   read(lun)BASISSET%nbast
   read(lun)BASISSET%nprimbast
   read(lun)BASISSET%label
ENDIF

END SUBROUTINE read_basissetinfo_from_disk

!> \brief read the basisinfo structure from disk
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information
!> \param lun logic unit number of file to read from
SUBROUTINE read_basisinfo_from_disk(lun,BASIS)
implicit none
TYPE(BASISINFO),intent(inout)    :: BASIS
INTEGER,intent(in)            :: lun
!
integer :: I
call nullifyMainbasis(BASIS)
DO I=1,nBasisBasParam
   read(lun)BASIS%WBASIS(I)
   IF(BASIS%WBASIS(I))THEN
      CALL READ_BASISSETINFO_FROM_DISK(lun,BASIS%BINFO(I))
   ENDIF
ENDDO

END SUBROUTINE read_basisinfo_from_disk

!> \brief initialise basissetinfo structure
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information for given set
!> \param natomtypes is number of atomtypes
SUBROUTINE init_basissetinfo_types(BasisInfo,natomtypes)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
INTEGER,intent(in)               :: natomtypes
!
integer :: i

BasisInfo%natomtypes=natomtypes
CALL MEM_ALLOC(BasisInfo%AtomType,natomtypes)
do I=1,natomtypes
   call nullifyAtomtype(BasisInfo%AtomType(I))
enddo
END SUBROUTINE init_basissetinfo_types

!> \brief initialise shell in basissetinfo structure
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information for given set
!> \param natomtype is number of atomtypes
!> \param nAngmom is the current angularmoment
!> \param segments is number of segments in the shell for this angular moment
SUBROUTINE init_basissetinfo_ContractionM(BasisInfo,natomtype,&
                                                         &nAngmom,segments)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
INTEGER,intent(in)            :: natomtype,nAngmom,segments

BasisInfo%AtomType(natomtype)%nAngmom=nAngmom
BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%nsegments=segments
IF(segments .GT. maxBASISsegment) THEN
   print*,'You have alot of segments you probably use a very&
        & large basisset which is not generally contracted. This&
        & is okay, but you need to increase the parameter&
        & maxBASISsegment in the lsutil/BasisinfoType.f90 file and recompile.&
        & Currently the parameter is set to ',maxBASISsegment
   print*,'You must increase the number in TYPE-DEF.f90 to at least ',segments
   print*,'Thomas Kjaergaard'
   CALL LSQUIT('Increase maxBASISsegment in TYPE-DEF.f90 file',-1)
ENDIF

END SUBROUTINE init_basissetinfo_ContractionM

!> \brief initialise exponents and contractionmatrix in basissetinfo structure
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param basisinfo Contains the all basis information for given set
!> \param natomtype is number of atomtypes
!> \param nAngmom is the current angularmoment
!> \param nsegments is number of segments in the shell for this angular moment
!> \param nrow is the number of exponents and the number of rows for the contraction matrix
!> \param ncol is the number of collums for the contraction matrix
SUBROUTINE init_basissetinfo_elms(BasisInfo,natomtype,nAngmom,&
                                                       &segments,nrow,ncol)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
INTEGER,intent(in)            :: natomtype,nAngmom,segments,nrow,ncol
!
INTEGER            :: nsize,i
nsize=nrow*ncol
BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                     &segment(segments)%nrow=nrow
BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                      &segment(segments)%ncol=ncol
CALL MEM_ALLOC(BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                    &segment(segments)%elms,nsize)
CALL MEM_ALLOC(BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                                    &segment(segments)%UCCelms,nsize)
do i = 1,nsize
  BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                   &segment(segments)%elms(i) = 0.0E0_realk
enddo
do i = 1,nsize
  BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                   &segment(segments)%UCCelms(i) = 0.0E0_realk
enddo
CALL MEM_ALLOC(BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                       &segment(segments)%Exponents,nrow)
do i = 1,nrow
   BasisInfo%AtomType(natomtype)%SHELL(nAngmom)%&
                   &segment(segments)%Exponents(i) = 0.0E0_realk
enddo

END SUBROUTINE init_basissetinfo_elms

!> \brief call mem_alloc and build a basissetinfo for a given type of a full basissetinfo  
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param oldbas is the old full BASISSETINFO
!> \param itype is the type in the full basissetinfo that is requested 
!> \param newbas is the new small BASISSETINFO which only contain 1 type, the itype in OLDBAS
SUBROUTINE alloc_and_take_subbasissetinfo(OLDBAS,itype,NEWBAS)
  implicit none
  TYPE(BASISSETINFO),intent(in)    :: OLDBAS
  TYPE(BASISSETINFO),intent(inout) :: NEWBAS
  INTEGER,intent(in)  :: itype
!
  INTEGER            :: nsegments,isegment,nbast,nprimbast
  INTEGER            :: iang,nrow,ncol,icharge
  
  nbast = 0
  nprimbast = 0
  NEWBAS%Gcont = OLDBAS%Gcont
  NEWBAS%DunningsBasis = OLDBAS%DunningsBasis
  NEWBAS%GeminalScalingFactor = OLDBAS%GeminalScalingFactor
  NEWBAS%GCbasis = OLDBAS%GCbasis
  NEWBAS%Spherical = OLDBAS%Spherical
  NEWBAS%labelindex = OLDBAS%labelindex
  CALL MEM_ALLOC(NEWBAS%ATOMTYPE,1)
  NEWBAS%ATOMTYPE(1)%name = OLDBAS%ATOMTYPE(itype)%name
  NEWBAS%nATOMTYPES = 1
  NEWBAS%ATOMTYPE(1) = OLDBAS%ATOMTYPE(itype)
  DO iang = 1,OLDBAS%ATOMTYPE(itype)%nangmom
     NEWBAS%ATOMTYPE(1)%SHELL(iang) = OLDBAS%ATOMTYPE(itype)%SHELL(iang)
     nsegments = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%nsegments
     DO isegment = 1,nsegments
        nrow = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%nrow
        ncol = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%ncol
        nbast = nbast+ncol
        nprimbast = nprimbast+nrow
        NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment) = OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)
        CALL MEM_ALLOC(NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%elms,nrow*ncol)
        CALL MEM_ALLOC(NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%UCCelms,nrow*ncol)
        CALL MEM_ALLOC(NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%Exponents,nrow)
        NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%elms = &
             &OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%elms
        NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%UCCelms = &
             &OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%UCCelms
        NEWBAS%ATOMTYPE(1)%SHELL(iang)%segment(isegment)%Exponents = &
             &OLDBAS%ATOMTYPE(itype)%SHELL(iang)%segment(isegment)%Exponents
     ENDDO
  ENDDO
  NEWBAS%nbast = nbast
  NEWBAS%nprimbast = nprimbast
  IF(OLDBAS%labelindex .EQ. 0)THEN
     IF(OLDBAS%nChargeIndex.NE. 0)THEN
        icharge = OLDBAS%ATOMTYPE(itype)%Charge
        CALL MEM_ALLOC(NEWBAS%Chargeindex,icharge,.TRUE.)
        NEWBAS%Chargeindex = 0
        NEWBAS%Chargeindex(icharge) = 1
        NEWBAS%nChargeindex = icharge
     ELSE
        NEWBAS%nChargeindex = 0
     ENDIF
  ENDIF
 
END SUBROUTINE alloc_and_take_subbasissetinfo

!> \brief copy basissetinfo to another basissetinfo  
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param oldbas is the original BASISSETINFO
!> \param newbas is the copied BASISSETINFO
SUBROUTINE copy_basissetinfo(OLDBAS,NEWBAS)
  implicit none
  TYPE(BASISSETINFO),intent(in)    :: OLDBAS
  TYPE(BASISSETINFO),intent(inout) :: NEWBAS
!
  INTEGER            :: I,J,K,nrow,ncol

  call nullifyBasisset(NEWBAS)
  NEWBAS%natomtypes = OLDBAS%natomtypes
  NEWBAS%Gcont = OLDBAS%Gcont
  NEWBAS%DunningsBasis = OLDBAS%DunningsBasis
  NEWBAS%GeminalScalingFactor = OLDBAS%GeminalScalingFactor
  NEWBAS%Gcbasis = OLDBAS%GCbasis
  NEWBAS%Spherical = OLDBAS%Spherical
  NEWBAS%labelindex = OLDBAS%labelindex
  NEWBAS%nChargeindex = OLDBAS%nChargeindex
  NEWBAS%nbast = OLDBAS%nbast
  NEWBAS%nprimbast = OLDBAS%nprimbast
  NEWBAS%label = OLDBAS%label
  NULLIFY(NEWBAS%ATOMTYPE)
  IF(NEWBAS%natomtypes.GT. 0)THEN
     CALL MEM_ALLOC(NEWBAS%ATOMTYPE,NEWBAS%natomtypes)
     DO I = 1,NEWBAS%natomtypes
        NEWBAS%ATOMTYPE(I)%nAngmom = OLDBAS%ATOMTYPE(I)%nAngmom 
        NEWBAS%ATOMTYPE(I)%family = OLDBAS%ATOMTYPE(I)%family 
        NEWBAS%ATOMTYPE(I)%ToTnorb = OLDBAS%ATOMTYPE(I)%ToTnorb
        NEWBAS%ATOMTYPE(I)%ToTnprim = OLDBAS%ATOMTYPE(I)%ToTnprim
        NEWBAS%ATOMTYPE(I)%Charge = OLDBAS%ATOMTYPE(I)%Charge
        NEWBAS%ATOMTYPE(I)%NAME  = OLDBAS%ATOMTYPE(I)%NAME 
        DO J=1,NEWBAS%ATOMTYPE(I)%nAngmom
           NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = OLDBAS%ATOMTYPE(I)%SHELL(J)%nprim
           NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = OLDBAS%ATOMTYPE(I)%SHELL(J)%norb
           NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments  = OLDBAS%ATOMTYPE(I)%SHELL(J)%nsegments
           DO K=1,NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments
              nrow = OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow
              ncol = OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow  = nrow
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol  = ncol
              CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms,nrow*ncol)
              CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms,nrow*ncol)
              CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents,nrow)
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(1:nrow*ncol) = &
                   &OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(1:nrow*ncol)
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol) = &
                   &OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol)
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(1:nrow) = &
                   &OLDBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(1:nrow)
           ENDDO
        ENDDO
        DO J=NEWBAS%ATOMTYPE(I)%nAngmom+1,maxAOangmom
           NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim = -1
           NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = -1
           NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments  = -1
           DO K=1,maxBASISsegment
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow  = -1
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol  = -1
              NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms)
              NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms)
              NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  IF(NEWBAS%nChargeindex .NE. 0)THEN
     CALL MEM_ALLOC(NEWBAS%Chargeindex,NEWBAS%nChargeindex,.TRUE.)
     DO I = 0,NEWBAS%nChargeindex
        NEWBAS%Chargeindex(I) = OLDBAS%Chargeindex(I)  
     ENDDO
  ELSE
     NULLIFY(NEWBAS%Chargeindex)
  ENDIF
 
END SUBROUTINE copy_basissetinfo

!> \brief add two basissetinfos together to generat a third basissetinfo  
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param lupri output logical unit number. 
!> \param iprint print level integer
!> \param bas1 is the original BASISSETINFO1
!> \param bas2 is the original BASISSETINFO2
!> \param newbas is the new output BASISSETINFO
!> \param Newlabelindex  labelindex of new basis
SUBROUTINE add_basissetinfo(lupri,iprint,BAS1,BAS2,NEWBAS,NewLabelindex)
  implicit none
  integer,intent(in)               :: lupri,iprint,NewLabelindex
  TYPE(BASISSETINFO),intent(in)    :: BAS1,BAS2
  TYPE(BASISSETINFO),intent(inout) :: NEWBAS
!
  INTEGER            :: I,J,K,nrow,ncol,nK,MAXANGMOM

  call nullifyBasisset(NEWBAS)
  IF(BAS1%natomtypes.NE.BAS2%natomtypes)THEN
     call lsquit('natomtypes mismatch in add_basissetinfo',-1)
  ENDIF
  NEWBAS%natomtypes = BAS1%natomtypes
  NEWBAS%Gcont = BAS1%Gcont .OR. BAS2%Gcont
  NEWBAS%DunningsBasis = BAS1%DunningsBasis .OR. BAS2%DunningsBasis
  NEWBAS%GeminalScalingFactor = BAS1%GeminalScalingFactor
  NEWBAS%Gcbasis = BAS1%GCbasis
  NEWBAS%Spherical = BAS1%Spherical
  IF(BAS1%labelindex.EQ.0.AND.BAS2%labelindex.EQ.0)THEN
     NEWBAS%labelindex = 0
     NEWBAS%nChargeindex = BAS1%nChargeindex
  ELSEIF(BAS1%labelindex.EQ.0.AND.BAS2%labelindex.NE.0)THEN
     call lsquit('labelindex mismatch in add_basissetinfo',-1)
  ELSE
     NEWBAS%labelindex = NewLabelindex
     NEWBAS%nChargeindex = 0
  ENDIF
  NEWBAS%nbast = BAS1%nbast + BAS2%nbast
  NEWBAS%nprimbast = BAS1%nprimbast + BAS2%nprimbast
  NEWBAS%label = BasParamLABEL(NewLabelindex)
  NULLIFY(NEWBAS%ATOMTYPE)
  IF(NEWBAS%natomtypes.GT. 0)THEN
   CALL MEM_ALLOC(NEWBAS%ATOMTYPE,NEWBAS%natomtypes)
   DO I = 1,NEWBAS%natomtypes
!    IF(BAS1%ATOMTYPE(I)%nAngmom .NE.BAS2%ATOMTYPE(I)%nAngmom)THEN
!       call lsquit('nAngmom mismatch in add_basissetinfo, not implemented contact TK',-1)
!    ENDIF
    MAXANGMOM=MAX(BAS1%ATOMTYPE(I)%nAngmom,BAS2%ATOMTYPE(I)%nAngmom)
    NEWBAS%ATOMTYPE(I)%nAngmom = MAXANGMOM
    NEWBAS%ATOMTYPE(I)%family = BAS1%ATOMTYPE(I)%family 
    NEWBAS%ATOMTYPE(I)%ToTnorb = BAS1%ATOMTYPE(I)%ToTnorb + BAS2%ATOMTYPE(I)%ToTnorb
    NEWBAS%ATOMTYPE(I)%ToTnprim = BAS1%ATOMTYPE(I)%ToTnprim + BAS2%ATOMTYPE(I)%ToTnprim
    NEWBAS%ATOMTYPE(I)%Charge = BAS1%ATOMTYPE(I)%Charge
    NEWBAS%ATOMTYPE(I)%NAME  = TRIM(BAS1%ATOMTYPE(I)%NAME)//'+'//TRIM(BAS2%ATOMTYPE(I)%NAME)
    DO J=1,MAXANGMOM
     NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = 0
     NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = 0
     NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments =0
     IF(BAS1%ATOMTYPE(I)%nAngmom.GE.J)THEN
        NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = BAS1%ATOMTYPE(I)%SHELL(J)%nprim 
        NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = BAS1%ATOMTYPE(I)%SHELL(J)%norb 
        NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments = BAS1%ATOMTYPE(I)%SHELL(J)%nsegments
     ENDIF
     IF(BAS2%ATOMTYPE(I)%nAngmom.GE.J)THEN
        NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  + BAS2%ATOMTYPE(I)%SHELL(J)%nprim
        NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = NEWBAS%ATOMTYPE(I)%SHELL(J)%norb + BAS2%ATOMTYPE(I)%SHELL(J)%norb
        NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments =NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments + BAS2%ATOMTYPE(I)%SHELL(J)%nsegments
     ENDIF
     IF(BAS1%ATOMTYPE(I)%nAngmom.GE.J)THEN
      DO K=1,BAS1%ATOMTYPE(I)%SHELL(J)%nsegments
       nrow = BAS1%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow
       ncol = BAS1%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow  = nrow
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol  = ncol
       CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms,nrow*ncol)
       CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms,nrow*ncol)
       CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents,nrow)
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(1:nrow*ncol) = &
            &BAS1%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(1:nrow*ncol)
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol) = &
            &BAS1%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol)
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(1:nrow) = &
            &BAS1%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(1:nrow)
      ENDDO
      nK = BAS1%ATOMTYPE(I)%SHELL(J)%nsegments
     ELSE
      nK = 0
     ENDIF
     IF(BAS2%ATOMTYPE(I)%nAngmom.GE.J)THEN
      DO K=1,BAS2%ATOMTYPE(I)%SHELL(J)%nsegments
       nrow = BAS2%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow
       ncol = BAS2%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%nrow  = nrow
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%ncol  = ncol
       CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%elms,nrow*ncol)
       CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%UCCelms,nrow*ncol)
       CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%Exponents,nrow)
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%elms(1:nrow*ncol) = &
            &BAS2%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(1:nrow*ncol)
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%UCCelms(1:nrow*ncol) = &
            &BAS2%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol)
       NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(nK+K)%Exponents(1:nrow) = &
            &BAS2%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(1:nrow)
      ENDDO
     ENDIF
    ENDDO !nangmom
    DO J=NEWBAS%ATOMTYPE(I)%nAngmom+1,maxAOangmom
     NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim = -1
     NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = -1
     NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments  = -1
     DO K=1,maxBASISsegment
      NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow  = -1
      NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol  = -1
      NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms)
      NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms)
      NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents)
     ENDDO
    ENDDO
   ENDDO
  ENDIF
  IF(NEWBAS%nChargeindex .NE. 0)THEN
     CALL MEM_ALLOC(NEWBAS%Chargeindex,NEWBAS%nChargeindex,.TRUE.)
     DO I = 0,NEWBAS%nChargeindex
        NEWBAS%Chargeindex(I) = BAS1%Chargeindex(I)  
     ENDDO
  ELSE
     NULLIFY(NEWBAS%Chargeindex)
  ENDIF
 
END SUBROUTINE add_basissetinfo

!> \brief call mem_dealloc BASISSETINFO
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param BasisInfo the BasissetInfo to be call mem_deallocd
SUBROUTINE free_basissetinfo(BasisInfo)
implicit none
TYPE(BASISSETINFO),intent(inout) :: BasisInfo
!
INTEGER  :: J,K,L,natomtypes,nAngmom,nsegments!,nsize,nrow,ncol
natomtypes=BASISINFO%natomtypes
IF(natomtypes.GT. 0)THEN
 DO J=1,natomtypes
   nangmom=BASISINFO%ATOMTYPE(J)%nAngmom
   DO K=1,nAngmom
      nsegments=BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
      IF(nsegments .NE. 0)THEN
         DO L=1,nsegments
            if (.not.ASSOCIATED(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)) then
               print*,'memory previously released!!'
               call lsquit('Error in FREE_BASISSETINFO1 - memory previously released',-1)
            endif
            CALL MEM_DEALLOC(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)
            if (.not.ASSOCIATED(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%UCCelms)) then
               print*,'memory previously released!!'
               call lsquit('Error in FREE_BASISSETINFO1 - memory previously released',-1)
            endif
            CALL MEM_DEALLOC(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%UCCelms)
            if (.not.ASSOCIATED(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)) then
               print*,'memory previously released!!'
               call lsquit('Error in FREE_BASISSETINFO2 - memory previously released',-1)
            endif
            CALL MEM_DEALLOC(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)
         ENDDO
      ENDIF
   ENDDO
 ENDDO
 if ((natomtypes.NE. 0).AND.(.not.ASSOCIATED(BASISINFO%ATOMTYPE))) then
   print*,'memory previously released!!'
   call lsquit('Error in FREE_BASISSETINFO4 - memory previously released',-1)
 endif
 CALL MEM_DEALLOC(BASISINFO%ATOMTYPE)
 IF(BASISINFO%Labelindex.EQ. 0)THEN
   IF(.not.ASSOCIATED(BASISINFO%Chargeindex))THEN
     print*,'memory previously released!!'
     call lsquit('Error in FREE_BASISSETINFO5  - memory previously released',-1)
   ENDIF
   CALL MEM_DEALLOC(BASISINFO%Chargeindex)   
   BASISINFO%nChargeindex=0
 ENDIF
 BASISINFO%natomtypes=0
ENDIF
END SUBROUTINE free_basissetinfo

!> \brief call mem_alloc BASISSETINFO, from already set values in BASISSETINFO, used in MPI parallelization
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param BasisInfo the BasissetInfo to be call mem_allocd
SUBROUTINE lsmpi_alloc_basissetinfo(BASISINFO)
IMPLICIT NONE
TYPE(BASISSETINFO),intent(inout) :: BASISINFO
!
INTEGER            :: J,K,L,nrow,nsize

CALL MEM_ALLOC(BASISINFO%ATOMTYPE,BASISINFO%natomtypes)
DO J=1,BASISINFO%nAtomtypes
   !NO need to call mem_alloc SHELL
   DO K=1,BASISINFO%ATOMTYPE(J)%nAngmom
      !NO need to call mem_alloc segments
      DO L=1,BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
         nrow=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%nrow
         nsize=nrow*BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%ncol
         CALL MEM_ALLOC(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms,nSIZE)
         CALL MEM_ALLOC(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%UCCelms,nSIZE)
         CALL MEM_ALLOC(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents,nrow)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE lsmpi_alloc_basissetinfo

SUBROUTINE print_brakebas(BB,LUPRI)
implicit none
type(BRAKEBASINFO),intent(in) :: BB
INTEGER,intent(in)      :: LUPRI

WRITE(Lupri,*)'OriginalType:',BB%OriginalType
WRITE(Lupri,*)'OriginalnAngmom:',BB%OriginalnAngmom
WRITE(Lupri,*)'nAngmom:',BB%nAngmom
WRITE(Lupri,*)'BB%OriginalnOrb:',BB%OriginalnOrb(1:BB%OriginalnAngmom)
WRITE(Lupri,*)'BB%nOrb:',BB%nOrb(1:BB%nAngmom)
WRITE(Lupri,*)'BB%Orb:',BB%Orb
WRITE(Lupri,*)'BB%Angmom:',BB%Angmom(1:BB%nAngmom)
END SUBROUTINE PRINT_BRAKEBAS

!> \brief print BASISSETINFO routine
!> \author T. Kjaergaard
!> \date  2010-02-24
!> \param lupri the logical unit number of the output file
!> \param BasisInfo the BasissetInfo to be printed
SUBROUTINE print_basissetinfo(LUPRI,BASISINFO)
implicit none
TYPE(BASISSETINFO),intent(in)  :: BASISINFO
INTEGER,intent(in)      :: LUPRI
!
INTEGER             :: nAngmom,natomtypes,J,K,L,Charge,segments
INTEGER             :: nPrimitives,ncol
CHARACTER(len=1)    :: SPDFGH(10)=(/'S','P','D','F','G','H','I','J','K','L'/) 

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'BASISSETINFORMATION'
WRITE(LUPRI,*)'nbast',BASISINFO%nbast
natomtypes=BASISINFO%natomtypes
WRITE(LUPRI,*)'Number of different types of Atoms with this basisset:',natomtypes
DO J=1,natomtypes
   nAngmom=BASISINFO%ATOMTYPE(J)%nAngmom
   Charge=BASISINFO%ATOMTYPE(J)%Charge
   WRITE(LUPRI,*)'------------------------------------------------------------'
   WRITE(LUPRI,'(A8,I4,A14,I4,A,A)')' TYPE   :',J,' has charge ',Charge,' Name:',TRIM(BASISINFO%ATOMTYPE(J)%NAME)
   IF(BASISINFO%ATOMTYPE(J)%FAMILY) WRITE(LUPRI,'(A)')' This atomtype is a Family basis set type'
   WRITE(LUPRI,*)'------------------------------------------------------------'
   DO K=1,nAngmom
      segments=BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
      WRITE(LUPRI,*)'Number of ContractionCoefficient blocks =',segments
      IF(segments .NE. 0)THEN
         DO L=1,segments
            nPrimitives=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%nrow
            ncol=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%ncol
            IF (K.LE. 10) THEN
               WRITE(LUPRI,'(2X,A1,A)')SPDFGH(K),'-TYPE FUNCTIONS'
            ELSE 
               WRITE(LUPRI,'(2X,I2,A)')K,'-TYPE FUNCTIONS'
            ENDIF
            WRITE(LUPRI,'(A22,I3,A22,I3)')'Exponents BLOCK=',L
            call LS_OUTPUT(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%&
                 &Exponents, 1, nPrimitives, 1, 1, nPrimitives, 1, 1, LUPRI)
            WRITE(LUPRI,*)'segment BLOCK=',L
            CALL LS_OUTPUT(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms,1,nPrimitives,1,ncol,nPrimitives,ncol,1,LUPRI)
            WRITE(LUPRI,*)'The unmodified and non normalized segment read from basis set file BLOCK=',L
            CALL LS_OUTPUT(BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%UCCelms,1,nPrimitives,1,ncol,nPrimitives,ncol,1,LUPRI)
         ENDDO
      ELSE
         IF (K.LE. 10) THEN
            WRITE(LUPRI,'(2X,A3,A1,A)')'NO ',SPDFGH(K),'-TYPE FUNCTIONS'
         ELSE 
            WRITE(LUPRI,'(2X,A3,I2,A)')'NO ',K,'-TYPE FUNCTIONS'
         ENDIF
      ENDIF
   ENDDO
ENDDO

WRITE(LUPRI,*)' '

END SUBROUTINE print_basissetinfo

subroutine freeBrakebasinfo(BASINFO)
implicit none
type(BRAKEBASINFO) :: BASINFO
call mem_dealloc(BASINFO%OriginalnOrb)
call mem_dealloc(BASINFO%nOrb)
call mem_dealloc(BASINFO%Orb) 
call mem_dealloc(BASINFO%angmom)
end subroutine freeBrakebasinfo

subroutine initBrakebasinfo(BASINFO,itype,nangmom,MaxnOrb)
implicit none
type(BRAKEBASINFO) :: BASINFO
integer,intent(in) :: itype,nangmom,MaxnOrb
BASINFO%OriginalType = itype
BASINFO%OriginalnAngmom = nAngmom
call mem_alloc(BASINFO%OriginalnOrb,nangmom)
BASINFO%OriginalnOrb = 0
call mem_alloc(BASINFO%nOrb,nangmom)
BASINFO%nOrb = 0
call mem_alloc(BASINFO%Orb,MaxnOrb,nangmom)
BASINFO%Orb = 0
call mem_alloc(BASINFO%Angmom,nangmom)
BASINFO%Angmom = 0
end subroutine initBrakebasinfo

subroutine copybrakebasinfo(OLDBASINFO,NEWBASINFO)
  implicit none
  type(BRAKEBASINFO),intent(in) :: OLDBASINFO
  type(BRAKEBASINFO),intent(inout) :: NEWBASINFO
  NEWBASINFO%OriginalType = OLDBASINFO%OriginalType
  NEWBASINFO%OriginalnAngmom = OLDBASINFO%OriginalnAngmom
  NEWBASINFO%nAngmom = OLDBASINFO%nAngmom
  call mem_alloc(NEWBASINFO%OriginalnOrb,SIZE(OLDBASINFO%OriginalnOrb))
  call mem_alloc(NEWBASINFO%nOrb,SIZE(OLDBASINFO%nOrb))
  call mem_alloc(NEWBASINFO%Orb,SIZE(OLDBASINFO%Orb,1),SIZE(OLDBASINFO%Orb,2))
  call mem_alloc(NEWBASINFO%Angmom,SIZE(OLDBASINFO%Angmom))
  NEWBASINFO%OriginalnOrb = OLDBASINFO%OriginalnOrb
  NEWBASINFO%nOrb = OLDBASINFO%nOrb
  NEWBASINFO%Orb = OLDBASINFO%Orb
  NEWBASINFO%Angmom = OLDBASINFO%Angmom
end subroutine copybrakebasinfo

subroutine buildbasisfrombrakebasinf(BB,NATOMS,nBB,oldBAS,newBAS,labelindex)
  implicit none
  integer,intent(in) :: nBB,natoms,labelindex
  type(BRAKEBASINFO),intent(in) :: BB(nBB)
  type(BASISSETINFO) :: oldBAS,newBAS

  INTEGER :: I,J,K,nrow,ncol,OrigI,OrigJ,iOrbIndex,iOrbIndex2,nOrbComp
  INTEGER :: newK,nSegmentOrb,nSegmentOrb2,L,LL,nprim,ToTnprim,ToTnOrb,JJ
  INTEGER :: FoundAngmom
  call nullifyBasisset(NEWBAS)
  NEWBAS%natomtypes = nBB
  NEWBAS%Gcont = OLDBAS%Gcont
  NEWBAS%DunningsBasis = OLDBAS%DunningsBasis
  NEWBAS%GeminalScalingFactor = OLDBAS%GeminalScalingFactor
  NEWBAS%Gcbasis = OLDBAS%GCbasis
  NEWBAS%Spherical = OLDBAS%Spherical
  !The labelindex cannot be 0 and nChargeindex must be zero
  IF(labelindex.EQ.0)call lsquit('labelindex cannot be 0 in buildbasisfrombrakebasinf',-1)
  NEWBAS%labelindex = labelindex
  NEWBAS%nChargeindex = 0 !OLDBAS%nChargeindex
  NEWBAS%label = OLDBAS%label
  NULLIFY(NEWBAS%ATOMTYPE)
  IF(NEWBAS%natomtypes.GT. 0)THEN
     CALL MEM_ALLOC(NEWBAS%ATOMTYPE,NEWBAS%natomtypes)
     DO I = 1,NEWBAS%natomtypes
        OrigI = BB(I)%OriginalType
        JJ = 1
        DO J=1,BB(I)%nAngmom
           JJ = MAX(JJ,BB(I)%Angmom(J)+1)
        ENDDO
        NEWBAS%ATOMTYPE(I)%nAngmom = JJ
        NEWBAS%ATOMTYPE(I)%family = OLDBAS%ATOMTYPE(OrigI)%family 
        NEWBAS%ATOMTYPE(I)%Charge = OLDBAS%ATOMTYPE(OrigI)%Charge
!WRONG CHANGE AT SOME POINT
        NEWBAS%ATOMTYPE(I)%NAME  = OLDBAS%ATOMTYPE(OrigI)%NAME 

        IF(BB(I)%nAngmom.NE.OLDBAS%ATOMTYPE(OrigI)%nAngmom)THEN
           NEWBAS%ATOMTYPE(I)%NAME  = TRIM(OLDBAS%ATOMTYPE(OrigI)%NAME)//'Reduced'
        ENDIF
        DO J=1,NEWBAS%ATOMTYPE(I)%nAngmom !the number in the the Original Basis, and new
         FoundAngmom = -1
         DO JJ=1,BB(I)%nAngmom
          IF(BB(I)%Angmom(JJ)+1.EQ.J)THEN
             FoundAngmom = JJ !the number in BB
          ENDIF
         ENDDO
         IF(FoundAngmom.NE.-1)THEN
           JJ = FoundAngmom !the number in BB
           OrigJ = J!BB(I)%Angmom(J)+1
           IF(BB(I)%nOrb(JJ).EQ.OLDBAS%ATOMTYPE(OrigI)%SHELL(OrigJ)%nOrb)THEN
              !all Orbitals included 
              NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = OLDBAS%ATOMTYPE(OrigI)%SHELL(OrigJ)%nprim
              NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = OLDBAS%ATOMTYPE(OrigI)%SHELL(OrigJ)%norb
              NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments  = OLDBAS%ATOMTYPE(OrigI)%SHELL(OrigJ)%nsegments
              DO K=1,NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments
                 nrow = OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%nrow
                 ncol = OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%ncol
                 NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow  = nrow
                 NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol  = ncol
                 CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms,nrow*ncol)
                 CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms,nrow*ncol)
                 CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents,nrow)
                 NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms(1:nrow*ncol) = &
                      &OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%elms(1:nrow*ncol)
                 NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol) = &
                      &OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol)
                 NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents(1:nrow) = &
                      &OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%Exponents(1:nrow)
              ENDDO
           ELSE
              OrigJ = J
              NEWBAS%ATOMTYPE(I)%NAME  = TRIM(OLDBAS%ATOMTYPE(OrigI)%NAME)//'Reduced'
              !not all orbitals included 
              NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = BB(I)%nOrb(JJ)
              
              iOrbIndex = 0
              newK = 0 
              nprim = 0
              DO K=1,OLDBAS%ATOMTYPE(OrigI)%SHELL(OrigJ)%nsegments
                 ncol = OLDBAS%ATOMTYPE(OrigI)%SHELL(OrigJ)%segment(K)%ncol
                 nrow = OLDBAS%ATOMTYPE(OrigI)%SHELL(OrigJ)%segment(K)%nrow
                 nSegmentOrb = 0
                 iOrbIndex2 = iOrbIndex
                 DO L=1,ncol
                    iOrbIndex2 = iOrbIndex2 + 1
                    DO LL = 1, BB(I)%nOrb(JJ)
                       IF(iOrbIndex2.EQ.BB(I)%Orb(LL,JJ))THEN
                          nSegmentOrb = nSegmentOrb + 1
                       ENDIF
                    ENDDO
                 ENDDO
                 IF(nSegmentOrb.EQ.ncol)THEN
                    newK = newK + 1
                    nprim = nprim + nrow
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%nrow  = nrow
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%ncol  = ncol
                    CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%elms,nrow*ncol)
                    CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%UCCelms,nrow*ncol)
                    CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%Exponents,nrow)
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%Exponents(1:nrow) = &
                         &OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%Exponents(1:nrow)
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%elms(1:nrow*ncol) = &
                         &OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%elms(1:nrow*ncol)
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%UCCelms(1:nrow*ncol) = &
                         &OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%UCCelms(1:nrow*ncol)
                 ELSEIF(nSegmentOrb.NE.0)THEN
                    newK = newK + 1
                    nprim = nprim + nrow
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%nrow  = nrow
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%ncol  = nSegmentOrb
                    CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%elms,nrow*nSegmentOrb)
                    CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%UCCelms,nrow*nSegmentOrb)
                    CALL MEM_ALLOC(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%Exponents,nrow)
                    NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%Exponents(1:nrow) = &
                         &OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%Exponents(1:nrow)
                    
                    nSegmentOrb2 = 0
                    iOrbIndex2 = iOrbIndex
                    DO L=1,ncol
                       iOrbIndex2 = iOrbIndex2 + 1
                       DO LL = 1, BB(I)%nOrb(JJ)
                          IF(iOrbIndex2.EQ.BB(I)%Orb(LL,JJ))THEN
                             nSegmentOrb2 = nSegmentOrb2 + 1
                             CALL COPYPARTSEG(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%elms,&
                                  & OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%elms,&
                                  & NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(newK)%UCCelms,&
                                  & OLDBAS%ATOMTYPE(OrigI)%SHELL(J)%segment(K)%UCCelms,nrow,ncol,&
                                  & nSegmentOrb,L,nSegmentOrb2)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDIF
                 iOrbIndex = iOrbIndex + ncol
              ENDDO
              NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = nprim
              NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments = newK 
           ENDIF
         ELSE
            !angmom not included 
            NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim  = 0
            NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = 0
            NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments = 0             
         ENDIF
        ENDDO
        DO J=NEWBAS%ATOMTYPE(I)%nAngmom+1,maxAOangmom
           NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim = -1
           NEWBAS%ATOMTYPE(I)%SHELL(J)%norb  = -1
           NEWBAS%ATOMTYPE(I)%SHELL(J)%nsegments  = -1
           DO K=1,maxBASISsegment
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%nrow  = -1
              NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%ncol  = -1
              NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%elms)
              NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%UCCelms)
              NULLIFY(NEWBAS%ATOMTYPE(I)%SHELL(J)%segment(K)%Exponents)
           ENDDO
        ENDDO
        ToTnorb = 0
        ToTnprim = 0 
        DO J=1,NEWBAS%ATOMTYPE(I)%nAngmom
           IF(NEWBAS%Spherical)THEN
              nOrbComp = 2*J-1
           ELSE
              nOrbComp = J*(J+1)/2
           ENDIF
           ToTnprim = ToTnprim + NEWBAS%ATOMTYPE(I)%SHELL(J)%nprim*nOrbComp
           ToTnorb = ToTnorb + NEWBAS%ATOMTYPE(I)%SHELL(J)%norb*nOrbComp
        ENDDO
        NEWBAS%ATOMTYPE(I)%ToTnorb = ToTnorb
        NEWBAS%ATOMTYPE(I)%ToTnprim = ToTnprim
     ENDDO
  ENDIF
  IF(NEWBAS%nChargeindex .NE. 0)THEN
     CALL MEM_ALLOC(NEWBAS%Chargeindex,NEWBAS%nChargeindex,.TRUE.)
     DO I = 0,NEWBAS%nChargeindex
        NEWBAS%Chargeindex(I) = OLDBAS%Chargeindex(I)  
     ENDDO
  ELSE
     NULLIFY(NEWBAS%Chargeindex)
  ENDIF

! will be determined using DETERMINE_NBAST
!  NEWBAS%nbast = OLDBAS%nbast
!  NEWBAS%nprimbast = OLDBAS%nprimbast
 
end subroutine buildbasisfrombrakebasinf

subroutine COPYPARTSEG(newelms,oldelms,newUCCelms,oldUCCelms,nrow,ncol,&
     & nSegmentOrb,L,nSegmentOrb2)
  implicit none
  integer :: nrow,ncol,nSegmentOrb,L,nSegmentOrb2
  real(realk) :: newelms(nrow,nSegmentOrb),newUCCelms(nrow,nSegmentOrb)
  real(realk) :: oldelms(nrow,ncol),oldUCCelms(nrow,ncol)
  newelms(1:nrow,nSegmentOrb2) = oldelms(1:nrow,L)
  newUCCelms(1:nrow,nSegmentOrb2) = oldUCCelms(1:nrow,L)
end subroutine COPYPARTSEG

END MODULE basis_type

