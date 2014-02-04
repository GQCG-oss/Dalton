!> @file
!> Module containing information about the OD-batches (combination of zero, one or two AO-batches)
MODULE ODbatches
 use precision
 use typedeftype
 use typedef
 use memory_handling
 use ls_util
 use OD_type
 !> An OD-batch is a set of:
 !>  i)   product overlaps between two batches of basis-functions -
 !>       where each AO-batch belong to a given center, and
 !>       share a common set of primitive basis-functions
 !>  ii)  AO basis-functions (beloning to one
 !>       center and with a set of primitive functions).
 !>  iii) an empty batch (used for debugging purposes only -
 !>       useful when testing contruction of E-coefficients).
 !>
 private
 public :: Create_ODbatches
CONTAINS
!> \brief Creates OD-batches from two sets of AO-batches
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param OD ODbatch to be built
!> \param INPUT integral input which contain all input to the integral eval.
!> \param SIDE label LHS or RHS 
!> \param IUNIT the logical unit number for the output file
SUBROUTINE Create_ODbatches(OD,INPUT,SIDE,IUNIT)
IMPLICIT NONE
 TYPE(ODITEM)      :: OD
 TYPE(INTEGRALINPUT),target :: INPUT
 LOGICAL           :: sameAOs
! Character(len=80) :: IDENTIFIER
 Character*(*)     :: SIDE
 Integer           :: IOD,IA,IB,IBSTART
 Integer           :: SA
 Integer           :: AOA,AOB,IUNIT
 Logical           :: screen,OPT1,OPT2,doMagscreen
 integer(kind=short) :: maxGab,maxpGab,maxGabElm,maxpGabElm
 Integer           :: i2,atomA,atomB
 Integer           :: batchA,batchB,ntypeA,nredtypeA
 real(realk)       :: d2,distance
! real(realk),parameter    :: R2PI52 = 5.91496717279561287782E0_realk
 type(lstensor),pointer :: GAB
 real(realk),parameter :: TEN=1E1_realk
!                     R2PI52 = sqrt(2 * sqrt(PI^5) )
 OPT1=(INPUT%MBIE_SCREEN.OR.INPUT%MBIE_INT)
 OPT2=(INPUT%DO_MULMOM.OR.INPUT%NonClassical_SCREEN)
 IF(INPUT%OE_Screen.AND.OPT1)THEN
    CALL LSQUIT('Input error to Create_ODbatches: multiple use of ODcenter, ODextent',iunit)
 ENDIF
 IF(INPUT%OE_Screen.AND.OPT2)THEN
    CALL LSQUIT('Input error to Create_ODbatches: multiple use of ODcenter, ODextent',iunit)
 ENDIF
 IF(OPT1.AND.OPT2)THEN
    CALL LSQUIT('Input error to Create_ODbatches: multiple use of ODcenter, ODextent',iunit)
 ENDIF

 SELECT CASE(SIDE)
 CASE('LHS')
    sameAOs=INPUT%sameLHSaos
    AOA = 1
    AOB = 2
    IF(INPUT%CS_SCREEN.OR.INPUT%PS_SCREEN)GAB => Input%LST_GAB_LHS
    IF(INPUT%CS_SCREEN)maxGabElm = INPUT%CS_MAXELM_RHS
    IF(INPUT%PS_SCREEN)maxpGabElm = INPUT%PS_MAXELM_RHS
    doMagscreen = INPUT%doMagScreen.AND.(INPUT%magderOrderP.EQ.1.OR.INPUT%LinComCarmomType.GT.1)
    CALL INIT_OD(OD,INPUT%AO(AOA)%p,INPUT%AO(AOB)%p,sameAOs,GAB,&
         & Input%CS_Screen,Input%PS_Screen,Input%OD_Screen,Input%CS_thrlog,&
         & Input%PS_thrlog,Input%CS_MAXELM_RHS,Input%PS_MAXELM_RHS,doMagscreen)
 CASE('RHS')
    sameAOs=INPUT%sameRHSaos
    AOA = 3
    AOB = 4
    IF(INPUT%CS_SCREEN.OR.INPUT%PS_SCREEN)GAB => Input%LST_GAB_RHS
    IF(INPUT%CS_SCREEN)maxGabElm = INPUT%CS_MAXELM_LHS
    IF(INPUT%PS_SCREEN)maxpGabElm = INPUT%PS_MAXELM_LHS
    doMagscreen = INPUT%doMagScreen.AND.INPUT%magderOrderQ.EQ.1
    CALL INIT_OD(OD,INPUT%AO(AOA)%p,INPUT%AO(AOB)%p,sameAOs,GAB,&
         & Input%CS_Screen,Input%PS_Screen,Input%OD_Screen,Input%CS_thrlog,&
         & Input%PS_thrlog,Input%CS_MAXELM_LHS,Input%PS_MAXELM_LHS,doMagscreen)
 CASE DEFAULT
    WRITE(IUNIT,'(1X,2A)') 'Wrong case in Create_ODbatches =',SIDE
    CALL LSQUIT('Wrong case in Create_ODbatches',-1)
 END SELECT

 OD%maxJ = INPUT%AO(AOA)%p%maxJ+INPUT%AO(AOB)%p%maxJ

! OD%identifier = IDENTIFIER  -  not used
 OD%sameAOs    = sameAOs
 ntypeA = INPUT%AO(AOA)%p%ntype
 nredtypeA = INPUT%AO(AOA)%p%nredtype
 IOD = 0
 DO IA=1,INPUT%AO(AOA)%p%nbatches
   IBSTART = 1
   IF (sameAOs) IBSTART = IA
   DO IB=IBSTART,INPUT%AO(AOB)%p%nbatches
!    Determine screening criteria
     screen = .FALSE.
!    Extent screening
     IF (INPUT%OD_SCREEN) THEN
       CALL getODscreening(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB),screen) 
     ENDIF
     IF(doMagscreen)THEN
         !Magnetic differentiation consist of two terms
         ! \f$ Q_{MN} (\mu \nu| frac{r_{1}}{r_{12}} | \rho \sigma ) \f$ 
         ! and 
         ! \f$ - Q_{RS} (\mu \nu| frac{r_{2}}{r_{12}} | \rho \sigma ) \f$ 
         ! where M,N,R,S is the center of atomic orbital \mu,\nu,\rho,\sigma
         ! and where Q_{MN} is a factor which depend on the distance 
         !between the two AOs (\mu and \nu) on the overlap to be differentiated
         !if the two AOs belong to the same atom, the magnetic differentiation 
         !will not contribute. See (JCP 95,2595)
         distance = INPUT%AO(AOA)%p%BATCH(IA)%center(1)-INPUT%AO(AOB)%p%BATCH(IB)%center(1)
         d2 = distance*distance
         distance = INPUT%AO(AOA)%p%BATCH(IA)%center(2)-INPUT%AO(AOB)%p%BATCH(IB)%center(2)
         d2 = d2 + distance*distance
         distance = INPUT%AO(AOA)%p%BATCH(IA)%center(3)-INPUT%AO(AOB)%p%BATCH(IB)%center(3)
         d2 = d2 + distance*distance
         IF(d2.LT.1E-5_realk) screen = .TRUE.     
     ENDIF
!    OD screenin must come before CS screening because Gab matrices are 
!    allocd according to OD screening. 
!    Cauchy-Scwartz screening |(a|b)| <= sqrt{(a|a)} sqrt{(b|b)}
     IF(.NOT.screen.AND.INPUT%CS_SCREEN)THEN
        maxGab = GAB%maxgab(IA,IB)
        IF (maxGab.LE.Input%CS_Thrlog-maxGabElm)screen = .TRUE.
     ENDIF
     IF(.NOT.screen.AND.INPUT%PS_SCREEN)THEN
        maxpGab = GAB%maxprimgab(IA,IB)
        IF (maxpGab.LE.Input%PS_Thrlog-maxpGabElm)screen = .TRUE.
     ENDIF
     IF (.NOT.screen) THEN
!      Create new OD-batch
       IOD         = IOD + 1
       OD%BATCH(IOD)%IA            =  IA
       OD%BATCH(IOD)%IB            =  IB
       OD%BATCH(IOD)%AO(1)%p       => INPUT%AO(AOA)%p%BATCH(IA)
       OD%BATCH(IOD)%AO(2)%p       => INPUT%AO(AOB)%p%BATCH(IB)
       OD%BATCH(IOD)%nAngmom       =  INPUT%AO(AOA)%p%BATCH(IA)%nAngmom*INPUT%AO(AOB)%p%BATCH(IB)%nAngmom
       OD%BATCH(IOD)%nPrimitives   =  INPUT%AO(AOA)%p%BATCH(IA)%nPrimitives*INPUT%AO(AOB)%p%BATCH(IB)%nPrimitives
       OD%BATCH(IOD)%maxContracted =  INPUT%AO(AOA)%p%BATCH(IA)%maxContracted*INPUT%AO(AOB)%p%BATCH(IB)%maxContracted
       OD%BATCH(IOD)%maxGAB = shortzero
       OD%BATCH(IOD)%ITYPE         = INPUT%AO(AOA)%p%BATCH(IA)%ITYPE+(INPUT%AO(AOB)%p%BATCH(IB)%ITYPE-1)*ntypeA
       OD%BATCH(IOD)%redTYPE       = INPUT%AO(AOA)%p%BATCH(IA)%redTYPE+(INPUT%AO(AOB)%p%BATCH(IB)%redTYPE-1)*nredtypeA
       OD%BATCH(IOD)%sameAO = (sameAOs .AND. IA .EQ. IB)
       OD%BATCH(IOD)%ODcenter(1) = 0.0E0_realk
       OD%BATCH(IOD)%ODcenter(2) = 0.0E0_realk
       OD%BATCH(IOD)%ODcenter(3) = 0.0E0_realk
       OD%BATCH(IOD)%ODextent = 0.0E0_realk
       
       IF (OD%BATCH(IOD)%sameAO) OD%BATCH(IOD)%nAngmom &
     &       = INPUT%AO(AOA)%p%BATCH(IA)%nAngmom*(INPUT%AO(AOA)%p%BATCH(IA)%nAngmom+1)/2
!     
       IF(INPUT%CS_SCREEN) OD%BATCH(IOD)%maxGAB  = maxGab
       IF (INPUT%MBIE_SCREEN .OR. INPUT%MBIE_INT) THEN
          CALL getMBIEcenterAndExtent(INPUT%AO(AOA)%p%BATCH(IA),&
               & INPUT%AO(AOB)%p%BATCH(IB),OD%BATCH(IOD),IUNIT)
       ENDIF
!      Non-classical screening screen away all contributions that can be calculated 
!      classically to within some threshold. For this screning we therefore set up
!      an (effective) center as well as the non-classical extent of the OD-batches
       IF(INPUT%DO_MULMOM.OR.INPUT%NonClassical_SCREEN)THEN
         CALL getODcenter(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB), &
     &                    OD%BATCH(IOD),INPUT%MM_SCREENTHR,IUNIT)
         CALL getODextentNonClassical(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB),&
     &                                OD%BATCH(IOD),INPUT%MM_SCREENTHR,IUNIT)
!      For overlap-integrals we screen if the distance between the LHS and RHS OB-batches exceeds
!      the sum of their extents. The OD-batch centers are the same weigthed center as for the 
!      non-classical screening, the primitive extents are determined by when the (gaussian product)
!      OD-function becomes smaller than a given threshold. As for non-classical screening, the
!      maximum primitive extent, augmented by the distance to the weigthed OD-center, is 
!      taken as the extent.
       ELSEIF (INPUT%OE_Screen) THEN
         CALL getODcenter(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB), &
     &                    OD%BATCH(IOD),INPUT%OE_THRESHOLD,IUNIT)
         CALL getODextentOverlap(INPUT%AO(AOA)%p%BATCH(IA),INPUT%AO(AOB)%p%BATCH(IB),&
     &                           OD%BATCH(IOD),INPUT%OE_THRESHOLD,IUNIT)
       ENDIF
    ENDIF
 ENDDO
ENDDO
 
 IF (IOD.NE.OD%nbatches) THEN
   WRITE(IUNIT,'(1X,1A,I5)') 'Number of batches from Create_ODbatches', IOD
   WRITE(IUNIT,'(1X,1A,I5)') 'Number of batches from Get_Nodbatches  ', OD%nbatches
   WRITE(IUNIT,'(1X,1A)')    'Programming error: Mismatch between number of OD-batches'
   CALL LSQUIT('Mismatch between number of OD-batches in Create_ODbatches',-1)
 ENDIF
! CALL DETERMINE_ODITEM_MEM(OD)

END SUBROUTINE Create_ODbatches

!!$SUBROUTINE DETERMINE_ODITEM_MEM(OD)
!!$IMPLICIT NONE
!!$TYPE(ODitem)        :: OD
!!$INTEGER(KIND=long)  :: mem_alloc
!!$INTEGER             :: I
!!$
!!$mem_alloc = 0
!!$mem_alloc = mem_alloc+sizeof(OD%sameAOs)
!!$DO I = 1,OD%nbatches
!!$!mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%identifier)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%AO)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%nAngmom)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%nPrimitives)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%maxContracted)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%maxGAB)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%sameAO)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%IA)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%IB)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%ODcenter)
!!$mem_alloc = mem_alloc+SIZEOF(OD%BATCH(I)%ODextent)
!!$ENDDO
!!$
!!$END SUBROUTINE DETERMINE_ODITEM_MEM

!> \brief Determine the center of the overlap distributions
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param ODB the overlap distribution 
!> \param THRESH threshold for zero contribution   
!> \param IUNIT the logical unit number for the output file
SUBROUTINE getODcenter(AOA,AOB,ODB,THRESH,IUNIT)
IMPLICIT NONE
 TYPE(AOBATCH)       :: AOA,AOB
 TYPE(ODBATCH)       :: ODB
 Integer           :: IOD
 Integer           :: IUNIT
 real(realk)       :: X,Y,Z,DIST12,THRESH,SM,AM,BM
 !real(realk),allocatable  :: EXTPMAX(:),primcenter(:,:) 
 real(realk)       :: e1,e2,p,mu,MAXCC1SUM,CCSUM,MAXCC2SUM,pm1,WT
 real(realk)       :: FAC,DIST
 real(realk),pointer :: CC1SUM(:),CC2SUM(:)
 Integer           :: nP1,nP2,i2,i1,i12,iA1,iA2,iC1,iC2,n12
 real(realk)       :: R2PI52 = 5.91496717279561287782E0_realk
 logical           :: nonzero,nucleus

nucleus = AOA%type_Nucleus .OR. AOB%type_Nucleus
IF (nucleus) THEN
  ODB%ODextent = 0E0_realk
!Simen     Set up an artificial extent of 1 because the error in the electron-nuclear attraction
!Simen     is larger than the electronic repulsion contribution. Should be analyzed!
!Simen     Matched by an exual extent set up in mm_read_in_raw_data (mm_interface.f90)
  ODB%ODextent = 1E0_realk
  IF (AOA%type_Nucleus) THEN
    ODB%ODcenter(1) = AOA%center(1)
    ODB%ODcenter(2) = AOA%center(2)
    ODB%ODcenter(3) = AOA%center(3)
  ELSE IF (AOB%type_Nucleus) THEN
    ODB%ODcenter(1) = AOB%center(1)
    ODB%ODcenter(2) = AOB%center(2)
    ODB%ODcenter(3) = AOB%center(3)
  ELSE
    CALL LSQUIT('Error in getODcenter, both AO-batches are nuclei!',-1)
  ENDIF
ELSE
   IF(AOA%type_empty.AND.AOB%type_empty)THEN
      ODB%ODcenter(1) = 0.0E0_realk
      ODB%ODcenter(2) = 0.0E0_realk
      ODB%ODcenter(3) = 0.0E0_realk
   ELSE
      X = AOA%center(1)-AOB%center(1)
      Y = AOA%center(2)-AOB%center(2)
      Z = AOA%center(3)-AOB%center(3)
      DIST12 = X*X+Y*Y+Z*Z
      !NOT YET IMPLEMENTED HARDCODED FOR NOW
      SM=0.0E0_realk
      AM=0.0E0_realk
      BM=0.0E0_realk
      NONZERO = .FALSE.
      call mem_alloc(CC1SUM,AOA%nAngmom)
      call mem_alloc(CC2SUM,AOB%nAngmom)
      nP2=AOB%nPrimitives
      nP1=AOA%nPrimitives
      i12=0
      DO i2=1,nP2
         e2  = ODB%AO(2)%p%pexponents%elms(i2)
         DO i1=1,nP1
            e1  = ODB%AO(1)%p%pexponents%elms(i1)
            p = e1+e2
            pm1 = 1/p
            mu = e1*e2*pm1
            CC1SUM = 0E0_realk
            CC2SUM = 0E0_realk
            DO iA1 = 1,AOA%nAngmom
               DO iC1=1,AOA%nContracted(iA1)
                  CC1SUM(iA1) = CC1SUM(iA1) + ABS(AOA%pCC(iA1)%p%elms(i1+(iC1-1)*nP1))
               ENDDO
            ENDDO
            MAXCC1SUM = 0E0_realk
            DO iA1 = 1,AOA%nAngmom
               MAXCC1SUM = MAX(MAXCC1SUM,CC1SUM(iA1)) 
            ENDDO
            
            DO iA2 = 1,AOB%nAngmom
               DO iC2=1,AOB%nContracted(iA2)
                  CC2SUM(iA2) = CC2SUM(iA2) + ABS(AOB%pCC(iA2)%p%elms(i2+(iC2-1)*nP2))
               ENDDO
            ENDDO
            MAXCC2SUM = 0E0_realk
            DO iA2 = 1,AOB%nAngmom
               MAXCC2SUM = MAX(MAXCC2SUM,CC2SUM(iA2)) 
            ENDDO
            CCSUM = MAXCC1SUM*MAXCC2SUM
            
            FAC = R2PI52*EXP(-MU*DIST12)*pm1*CCSUM
            IF(FAC .GT. THRESH)THEN
               NONZERO = .TRUE.
               MAXCC1SUM = 0E0_realk
               DO iA1 = 1,AOA%nAngmom
                  MAXCC1SUM = MAXCC1SUM+CC1SUM(iA1) 
               ENDDO
               MAXCC2SUM = 0E0_realk
               DO iA2 = 1,AOB%nAngmom
                  MAXCC2SUM = MAXCC2SUM+CC2SUM(iA2) 
               ENDDO
               
               WT = MAXCC1SUM*MAXCC2SUM
               
               i12 = i12 + 1
               
               WT = WT * EXP(-MU*DIST12)
               SM = SM + WT
               AM = AM + WT * e1 * pm1
               BM = BM + WT * e2 * pm1
            ENDIF
         ENDDO
      ENDDO
      n12=i12
      call mem_dealloc(CC1SUM)
      call mem_dealloc(CC2SUM)
      
      IF(NONZERO)THEN
         AM = AM/SM
         BM = BM/SM
         ODB%ODcenter(1) = AM*AOA%center(1) + BM*AOB%center(1)
         ODB%ODcenter(2) = AM*AOA%center(2) + BM*AOB%center(2)
         ODB%ODcenter(3) = AM*AOA%center(3) + BM*AOB%center(3)
      ELSE
         ODB%ODcenter(1) = (AOA%center(1) + AOB%center(1))/2E0_realk
         ODB%ODcenter(2) = (AOA%center(2) + AOB%center(2))/2E0_realk
         ODB%ODcenter(3) = (AOA%center(3) + AOB%center(3))/2E0_realk
      ENDIF
   ENDIF
ENDIF
END SUBROUTINE getODcenter

!> \brief Determine the non classical overlap distribution extent
!> \author S. Reine
!> \date 2010
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param ODB the overlap distribution 
!> \param THRESH threshold for extent (when the OD falls below THRESH)
!> \param IUNIT the logical unit number for the output file
SUBROUTINE getODextentNonClassical(AOA,AOB,ODB,THRESH,IUNIT)
IMPLICIT NONE
 TYPE(AOBATCH)       :: AOA,AOB
 TYPE(ODBATCH)       :: ODB
 Integer           :: IOD
 Integer           :: IUNIT
 real(realk)       :: LS_ERFCIV, X,Y,Z,DIST12,THRESH,EXTPMAX2,FACCLS
 real(realk),pointer  :: EXTPMAX(:),primcenter(:,:) 
 real(realk)       :: e1,e2,p,mu,MAXCC1SUM,CCSUM,MAXCC2SUM,pm1,EXTP
 real(realk)       :: FAC,dist
 real(realk),pointer :: CC1SUM(:),CC2SUM(:)
 Integer           :: nP1,nP2,i2,i1,i12,iA1,iA2,iC1,iC2,n12
 real(realk),parameter       :: R2PI52 = 5.91496717279561287782E0_realk
 logical           :: nucleus

nucleus = AOA%type_Nucleus .OR. AOB%type_Nucleus
IF (nucleus) THEN
  ODB%ODextent = 0E0_realk
!Simen     Set up an artificial extent of 1 because the error in the electron-nuclear attraction
!Simen     is larger than the electronic repulsion contribution. Should be analyzed!
!Simen     Matched by an exual extent set up in mm_read_in_raw_data (mm_interface.f90)
  ODB%ODextent = 1E0_realk
ELSE
   IF(AOA%type_empty.AND.AOB%type_empty)THEN
      ODB%ODextent = 0E0_realk
   ELSE
      X = AOA%center(1)-AOB%center(1)
      Y = AOA%center(2)-AOB%center(2)
      Z = AOA%center(3)-AOB%center(3)
      DIST12 = X*X+Y*Y+Z*Z
      !NOT YET IMPLEMENTED HARDCODED FOR NOW
      FACCLS = LS_ERFCIV(THRESH)
      EXTPMAX2=0.0E0_realk
      call mem_alloc(CC1SUM,AOA%nAngmom)
      call mem_alloc(CC2SUM,AOB%nAngmom)
      nP2=AOB%nPrimitives
      nP1=AOA%nPrimitives
      call mem_alloc(EXTPMAX,nP1*nP2)
      call mem_alloc(Primcenter,3,nP1*nP2)
      i12=0
      DO i2=1,nP2
         e2  = ODB%AO(2)%p%pexponents%elms(i2)
         DO i1=1,nP1
            e1  = ODB%AO(1)%p%pexponents%elms(i1)
            p = e1+e2
            pm1 = 1/p
            mu = e1*e2*pm1
            CC1SUM = 0E0_realk
            CC2SUM = 0E0_realk
            DO iA1 = 1,AOA%nAngmom
               DO iC1=1,AOA%nContracted(iA1)
                  CC1SUM(iA1) = CC1SUM(iA1) + ABS(AOA%pCC(iA1)%p%elms(i1+(iC1-1)*nP1))
               ENDDO
            ENDDO
            MAXCC1SUM = 0E0_realk
            DO iA1 = 1,AOA%nAngmom
               MAXCC1SUM = MAX(MAXCC1SUM,CC1SUM(iA1)) 
            ENDDO
            
            DO iA2 = 1,AOB%nAngmom
               DO iC2=1,AOB%nContracted(iA2)
                  CC2SUM(iA2) = CC2SUM(iA2) + ABS(AOB%pCC(iA2)%p%elms(i2+(iC2-1)*nP2))
               ENDDO
            ENDDO
            MAXCC2SUM = 0E0_realk
            DO iA2 = 1,AOB%nAngmom
               MAXCC2SUM = MAX(MAXCC2SUM,CC2SUM(iA2)) 
            ENDDO
            CCSUM = MAXCC1SUM*MAXCC2SUM
            
            FAC = R2PI52*EXP(-MU*DIST12)*pm1*CCSUM
            IF(FAC .GT. THRESH)THEN
               i12 = i12 + 1
               
               EXTP = SQRT(pm1) * FACCLS
               EXTPMAX(i12) = EXTP
               
               Primcenter(1,i12) = (e1*AOA%center(1) + e2*AOB%center(1))/(e1+e2)
               Primcenter(2,i12) = (e1*AOA%center(2) + e2*AOB%center(2))/(e1+e2)
               Primcenter(3,i12) = (e1*AOA%center(3) + e2*AOB%center(3))/(e1+e2)
            ENDIF
         ENDDO
      ENDDO
      n12=i12
      call mem_dealloc(CC1SUM)
      call mem_dealloc(CC2SUM)
      
      EXTPMAX2 = 0E0_realk
      DO I12=1,n12
         X = ODB%ODcenter(1) - Primcenter(1,i12)
         Y = ODB%ODcenter(2) - Primcenter(2,i12)
         Z = ODB%ODcenter(3) - Primcenter(3,i12)
         DIST = SQRT(X*X+Y*Y+Z*Z)
         EXTPMAX2 = MAX(EXTPMAX2,EXTPMAX(i12)+DIST)
      ENDDO
      ODB%ODextent = EXTPMAX2
      call mem_dealloc(EXTPMAX)
      call mem_dealloc(Primcenter)
   ENDIF
ENDIF
END SUBROUTINE getODextentNonClassical

!> \brief Determine the overlap distribution extent
!> \author S. Reine
!> \date 2010
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param ODB the overlap distribution 
!> \param THRESH threshold for extent (when the OD falls below THRESH)
!> \param IUNIT the logical unit number for the output file
SUBROUTINE getODextentOverlap(AOA,AOB,ODB,THRESH,IUNIT)
IMPLICIT NONE
 TYPE(AOBATCH)       :: AOA,AOB
 TYPE(ODBATCH)       :: ODB
 Integer           :: IOD
 Integer           :: IUNIT
 real(realk)       :: X,Y,Z,DIST12,THRESH,EXTPMAX2
 real(realk),pointer  :: EXTPMAX(:),primcenter(:,:) 
 real(realk)       :: e1,e2,p,mu,MAXCC1SUM,CCSUM,MAXCC2SUM,pm1,EXTP
 real(realk)       :: FACSCREEN,FAC,DIST
 real(realk),pointer :: CC1SUM(:),CC2SUM(:)
 Integer           :: nP1,nP2,i2,i1,i12,iA1,iA2,iC1,iC2,n12
 real(realk)       :: R2PI52 = 5.91496717279561287782E0_realk
 logical           :: nucleus

nucleus = AOA%type_Nucleus .OR. AOB%type_Nucleus
IF (nucleus) THEN
  CALL LSQUIT('Not an option getODextentOverlap with nuclei',-1)
ELSE
  X = AOA%center(1)-AOB%center(1)
  Y = AOA%center(2)-AOB%center(2)
  Z = AOA%center(3)-AOB%center(3)
  DIST12 = X*X+Y*Y+Z*Z
  EXTPMAX2=0.0E0_realk
  call mem_alloc(CC1SUM,AOA%nAngmom)
  call mem_alloc(CC2SUM,AOB%nAngmom)
  nP2=AOB%nPrimitives
  nP1=AOA%nPrimitives
  call mem_alloc(EXTPMAX,nP1*nP2)
  call mem_alloc(Primcenter,3,nP1*nP2)
  i12=0
  DO i2=1,nP2
     e2  = ODB%AO(2)%p%pexponents%elms(i2)
     DO i1=1,nP1
        e1  = ODB%AO(1)%p%pexponents%elms(i1)
        p = e1+e2
        pm1 = 1/p
        mu = e1*e2*pm1
        CC1SUM = 0E0_realk
        CC2SUM = 0E0_realk
        DO iA1 = 1,AOA%nAngmom
           DO iC1=1,AOA%nContracted(iA1)
              CC1SUM(iA1) = CC1SUM(iA1) + ABS(AOA%pCC(iA1)%p%elms(i1+(iC1-1)*nP1))
           ENDDO
        ENDDO
        MAXCC1SUM = 0E0_realk
        DO iA1 = 1,AOA%nAngmom
           MAXCC1SUM = MAX(MAXCC1SUM,CC1SUM(iA1)) 
        ENDDO
        
        DO iA2 = 1,AOB%nAngmom
           DO iC2=1,AOB%nContracted(iA2)
              CC2SUM(iA2) = CC2SUM(iA2) + ABS(AOB%pCC(iA2)%p%elms(i2+(iC2-1)*nP2))
           ENDDO
        ENDDO
        MAXCC2SUM = 0E0_realk
        DO iA2 = 1,AOB%nAngmom
           MAXCC2SUM = MAX(MAXCC2SUM,CC2SUM(iA2)) 
        ENDDO
        CCSUM = MAXCC1SUM*MAXCC2SUM
        
!Simen  Not proper for overlap integrals, refine!
        FAC = R2PI52*EXP(-MU*DIST12)*pm1*CCSUM
        IF(FAC .GT. THRESH)THEN
           i12 = i12 + 1
           
           EXTP = getODbatchOverlapExtent(p,AOA%maxAngmom+AOB%maxAngmom,&
     &                                    EXP(-MU*DIST12)*CCSUM,THRESH)
           
           EXTPMAX(i12) = EXTP

           Primcenter(1,i12) = (e1*AOA%center(1) + e2*AOB%center(1))/(e1+e2)
           Primcenter(2,i12) = (e1*AOA%center(2) + e2*AOB%center(2))/(e1+e2)
           Primcenter(3,i12) = (e1*AOA%center(3) + e2*AOB%center(3))/(e1+e2)
        ENDIF
     ENDDO
  ENDDO
  n12=i12
  call mem_dealloc(CC1SUM)
  call mem_dealloc(CC2SUM)
  
  EXTPMAX2 = 0E0_realk
  DO I12=1,n12
     X = ODB%ODcenter(1) - Primcenter(1,i12)
     Y = ODB%ODcenter(2) - Primcenter(2,i12)
     Z = ODB%ODcenter(3) - Primcenter(3,i12)
     DIST = SQRT(X*X+Y*Y+Z*Z)
     EXTPMAX2 = MAX(EXTPMAX2,EXTPMAX(i12)+DIST)
  ENDDO
  ODB%ODextent = EXTPMAX2
  call mem_dealloc(EXTPMAX)
  call mem_dealloc(Primcenter)
ENDIF
END SUBROUTINE getODextentOverlap

!> \brief Determine the overlap distribution extent
!> \author S. Reine
!> \date 2010
!> \param exponent the p=p1+p2 exponent of the OD 
!> \param angmax the maximum angular momentum
!> \param prefactor the prefactor exp(-mu *Rab), mu is reduced exponent
!> \param Threshold the threshold
!>
!>Function to return the extent of gaussian solid harmonic function;
!>   prefactor*S^{angmax,m}(r_P)*exp(-p r_R^2) 
!>      .LE. prefactor*r_P^l * exp(-p r_R^2) .LE. threshold
!>
REAL(REALK) FUNCTION getODbatchOverlapExtent(exponent,angmax,prefactor,threshold)
implicit none
REAL(REALK) :: exponent,prefactor,threshold
INTEGER     :: angmax
!
REAL(REALK) :: r2

r2 = -(1/exponent)*log(threshold/prefactor)
!Simen: Not exact for l>0 refine!
IF (r2.GT. 1.0E0_realk) r2=r2+angmax*log(sqrt(r2))
getODbatchOverlapExtent = sqrt(r2)

END FUNCTION getODbatchOverlapExtent

!> \brief Determine the number of OD-batches
!> \author S. Reine
!> \date 2010
!> \param nbatches the number of batches to be determined
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param sameAOs if the AOs on AOA and AOB is the same
!> \param GAB the screening lstensor
!> \param CS_screen if Cauchy-Schwarz screening is used
!> \param OD_screen if the overlap distribution screening is used
!> \param CS_thrlog the Cauchy-Schwarz screening threshold
!> \param maxgabelm the maximum element of the gab matrix
SUBROUTINE GET_NODBATCHES(nbatches,AOA,AOB,sameAOs,GAB,&
     CS_Screen,PS_screen,OD_Screen,CS_thrlog,PS_thrlog,maxGabElm,maxpGabElm,doMagscreen)
IMPLICIT NONE
 Integer      :: nbatches,IA,IB,IBSTART
 TYPE(AOITEM) :: AOA, AOB
 LOGICAL      :: sameAOs,doMagscreen
! Real(realk)  :: GAB(:,:)
type(lstensor),pointer :: GAB
 Logical      :: CS_Screen, OD_Screen,PS_Screen
 integer(kind=short) :: CS_thrlog,maxGabElm,PS_thrlog,maxpGabElm
!
 Integer      :: atomA,atomB,batcha,BatchB
 integer(kind=short)  :: maxGab,maxpGAB
 real(realk)  :: d2,distance
 Logical      :: screen

 nbatches=0
 DO IA=1,AOA%nbatches
   IBSTART = 1
   IF (sameAOs) IBSTART = IA
   DO IB=IBSTART,AOB%nbatches
     screen = .FALSE.
     IF (OD_SCREEN) THEN
        CALL getODscreening(AOA%BATCH(IA),AOB%BATCH(IB),screen) 
     ENDIF
     IF(doMagscreen)THEN
         !Magnetic differentiation consist of two terms
         ! \f$ Q_{MN} (\mu \nu| frac{r_{1}}{r_{12}} | \rho \sigma ) \f$ 
         ! and 
         ! \f$ - Q_{RS} (\mu \nu| frac{r_{2}}{r_{12}} | \rho \sigma ) \f$ 
         ! where M,N,R,S is the center of atomic orbital \mu,\nu,\rho,\sigma
         ! and where Q_{MN} is a factor which depend on the distance 
         !between the two AOs (\mu and \nu) on the overlap to be differentiated
         !if the two AOs belong to the same atom, the magnetic differentiation 
         !will not contribute. See (JCP 95,2595)
         distance = AOA%BATCH(IA)%center(1)-AOB%BATCH(IB)%center(1)
         d2 = distance*distance
         distance = AOA%BATCH(IA)%center(2)-AOB%BATCH(IB)%center(2)
         d2 = d2 + distance*distance
         distance = AOA%BATCH(IA)%center(3)-AOB%BATCH(IB)%center(3)
         d2 = d2 + distance*distance
         IF(d2.LT.1E-5_realk) screen = .TRUE.     
     ENDIF
!    OD screening must come before CS screening because Gab matrices are 
!    allocd according to OD screening. 
     IF(.NOT.screen.AND.CS_SCREEN)THEN
        maxGab = GAB%maxgab(IA,IB)
        IF (maxGab.LE.CS_Thrlog-maxGabElm)THEN
           screen = .TRUE.
        ENDIF
     ENDIF
     IF(.NOT.screen.AND.PS_SCREEN)THEN
        maxpGab = GAB%maxprimgab(IA,IB)
        IF (maxpGab.LE.PS_Thrlog-maxpGabElm) screen = .TRUE.
     ENDIF
     IF (.NOT.screen) THEN
       nbatches = nbatches + 1
     ENDIF
   ENDDO
 ENDDO
END SUBROUTINE GET_NODBATCHES

!> \brief Initialize OD-item
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param OD the Overlap distribution 
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param sameAOs if the AOs on AOA and AOB is the same
!> \param GAB the screening lstensor
!> \param CS_screen if Cauchy-Schwarz screening is used
!> \param OD_screen if the overlap distribution screening is used
!> \param CS_thrlog the Cauchy-Schwarz screening threshold
!> \param maxgabelm the maximum element of the gab matrix
SUBROUTINE Init_OD(OD,AOA,AOB,sameAOs,GAB,CS_Screen,PS_Screen,&
     &OD_Screen,CS_thrlog,PS_thrlog,maxGabElm,maxpGabElm,doMagscreen)
  use memory_handling
IMPLICIT NONE
 TYPE(ODITEM)   :: OD
 TYPE(AOITEM) :: AOA, AOB
 Integer      :: nbatches,SIZE
 LOGICAL      :: sameAOs,doMagscreen
 type(lstensor),pointer :: GAB
! Real(realk)  :: GAB(:,:) 
 Logical      :: OD_Screen,CS_Screen,PS_Screen
 integer (kind=short):: CS_thrlog,maxGabElm,PS_thrlog,maxpGabElm
 integer (kind=long) :: nsize
!
 CALL GET_NODBATCHES(nbatches,AOA,AOB,sameAOs,GAB,&
 &CS_Screen,PS_Screen,OD_Screen,CS_thrlog,PS_thrlog,maxGabElm,maxpGabElm,doMagscreen)
 OD%nbatches = nbatches
 IF(nbatches.GT. 0)THEN
    call mem_Alloc(OD%BATCH,nbatches)
 ELSE
    Nullify(OD%BATCH) 
 ENDIF
 !SIZE = SIZE OF ODBATCH
!$OMP CRITICAL (memory)
 SIZE = 5*mem_realsize+5*mem_intsize+1*mem_logicalsize
 nsize = SIZE*nbatches 
 call mem_allocated_mem_ODitem(nsize)
!$OMP END CRITICAL (memory)
END SUBROUTINE Init_OD

!> \brief calculate the maximum gab element
!> \author T. Kjaergaard
!> \date 2010
!> \param maxgab the maximum gab element
!> \param GAB an array/block/part of the screening tensor 
!> \param ndim the dimension of the small array 
SUBROUTINE CALC_MAXGAB(maxGab,GAB,ndim)
IMPLICIT NONE
INTEGER       :: I,ndim
real(realk)   :: maxGab
real(realk)   :: GAB(ndim)

maxGab = 0.0E0_realk
DO I=1,ndim
   maxGAB = MAX(maxGAB,GAB(I)) 
ENDDO

END SUBROUTINE CALC_MAXGAB

!> \brief Determine the center of the overlap distributions, and extent used in MBIE
!> \author T. Kjaergaard
!> \date 2010
!> \param AOA the first AO batch
!> \param AOB the second AO batch
!> \param ODB the overlap distribution 
!> \param IUNIT the logical unit number for the output file
SUBROUTINE getMBIEcenterAndExtent(AOA,AOB,ODB,IUNIT)
IMPLICIT NONE
 TYPE(AOBATCH)     :: AOA,AOB
 TYPE(ODBATCH)     :: ODB
 Integer           :: IOD
 Integer           :: IUNIT
 real(realk)       :: X,Y,Z,DIST12,THRESH,SM,AM,BM,DISTANCE
 real(realk)       :: e1,e2,p,mu,MAXCC1SUM,CCSUM,MAXCC2SUM,pm1,WT
 real(realk)       :: FAC,DIST,Sx,Sy,Sz,SCC
 real(realk),pointer :: CC(:),primcenter(:,:,:),center(:,:),primextent(:,:)
 Integer           :: nP1,nP2,i2,i1,i12,iA1,iA2,iC1,iC2,n12,iC12,la,lb,nC12
 real(realk),parameter       :: R2PI52 = 5.91496717279561287782E0_realk
 logical           :: nonzero,nucleus
 real(realk)       :: eta,thr,extent,primextent2,primextentTMP
 
 eta = -4.5E0_realk
 !in the paper (JCP 123, 184101) they suggest
 ! eta = -LOG10(THRESH)/3
 !but they also say that they use eta=4, which works well for THRESH=10**(-6) to THRESH=10**(-8)
 !(which according to their suggestions will give eta=2 and eta=2.666)
 !we like THRESH=10**(-10) as a default threshold so we need to increase num val of eta.  
 thr = 10.0E0_realk**eta
 !1. CALCULATE CONTRACTED CENTER OF THE CHARGE DISTRIBUTION          
 nP2=AOB%nPrimitives
 nP1=AOA%nPrimitives
 nC12=0
 DO iA1 = 1,AOA%nAngmom
    DO iA2 = 1,AOB%nAngmom
       nC12 = nC12+AOA%nContracted(iA1)*AOB%nContracted(iA2)
    ENDDO
 ENDDO
 call mem_alloc(CC,nP1*nP2)
 call mem_alloc(Primcenter,3,nP1*nP2,nC12)
 call mem_alloc(Primextent,nP1*nP2,nC12)
 call mem_alloc(Center,3,nC12)
 iC12=0
 DO iA1 = 1,AOA%nAngmom
    la=AOA%ANGMOM(iA1)
    DO iA2 = 1,AOB%nAngmom
       lb=AOB%ANGMOM(iA2)
       DO iC2=1,AOB%nContracted(iA2)
          DO iC1=1,AOA%nContracted(iA1)
             iC12 = iC12+1
             i12=0
             DO i2=1,nP2
                DO i1=1,nP1
                   i12=i12+1
                   e1  = ODB%AO(1)%p%pexponents%elms(i1)
                   e2  = ODB%AO(2)%p%pexponents%elms(i2)
                   Primcenter(1,i12,iC12) = (e1*AOA%center(1) + e2*AOB%center(1))/(e1+e2)
                   Primcenter(2,i12,iC12) = (e1*AOA%center(2) + e2*AOB%center(2))/(e1+e2)
                   Primcenter(3,i12,iC12) = (e1*AOA%center(3) + e2*AOB%center(3))/(e1+e2)
                   X = AOA%center(1)-AOB%center(1)
                   Y = AOA%center(2)-AOB%center(2)
                   Z = AOA%center(3)-AOB%center(3)
                   CC(i12)=AOA%pCC(iA1)%p%elms(i1+(iC1-1)*nP1)*AOB%pCC(iA2)%p%elms(i2+(iC2-1)*nP2)
                   CALL CALCULATE_ITERATIVE_EXTENT(primextentTMP,thr,la,lb,la+lb,e1,e2,e1+e2,e1*e2/(e1+e2),X,Y,Z,iunit)
                   primextent(i12,iC12) = CC(i12)*primextentTMP
                ENDDO
             ENDDO
             n12=i12
             Sx=0
             Sy=0
             Sz=0
             SCC=0
             DO i12=1,n12
                Sx = Sx+CC(i12)*primcenter(1,i12,iC12)
                Sy = Sy+CC(i12)*primcenter(2,i12,iC12)
                Sz = Sz+CC(i12)*primcenter(3,i12,iC12)
                SCC = SCC+CC(i12)
             ENDDO
             Center(1,iC12) = Sx/SCC
             Center(2,iC12) = Sy/SCC
             Center(3,iC12) = Sz/SCC
          ENDDO
       ENDDO
    ENDDO
 ENDDO
 nC12 = iC12
 Sx=0
 Sy=0
 Sz=0
 DO iC12=1,nC12
    Sx = Sx+Center(1,iC12)
    Sy = Sy+Center(2,iC12)
    Sz = Sz+Center(3,iC12)
 ENDDO
 ODB%ODCENTER(1) = Sx/nC12
 ODB%ODCENTER(2) = Sy/nC12
 ODB%ODCENTER(3) = Sz/nC12
 !2. EXTENT IS NOW THE DISTANCE FROM THE CONTRACTED CENTER TO TO PRIMITIV 
 !   CENTER FARTHEST AWAY PLUS THE EXTENT OF THE MAX PRIMITIVES
 iC12=0
 DISTANCE = 0.0E0_realk
 DO iA1 = 1,AOA%nAngmom
    DO iA2 = 1,AOB%nAngmom
       DO iC2=1,AOB%nContracted(iA2)
          DO iC1=1,AOA%nContracted(iA1)
             iC12 = iC12+1
             i12=0
             DO i2=1,nP2
                DO i1=1,nP1
                   i12=i12+1
                   X = Primcenter(1,i12,iC12)-ODB%ODCENTER(1) 
                   Y = Primcenter(2,i12,iC12)-ODB%ODCENTER(2) 
                   Z = Primcenter(3,i12,iC12)-ODB%ODCENTER(3)
                   DIST = SQRT(X*X+Y*Y+Z*Z)+primextent(i12,iC12) 
                   DISTANCE = MAX(DISTANCE,DIST)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO
 ODB%ODEXTENT = DISTANCE
 call mem_dealloc(CC)
 call mem_dealloc(Primcenter)
 call mem_dealloc(Primextent)
 call mem_dealloc(Center)

END SUBROUTINE getMBIEcenterAndExtent

SUBROUTINE CALCULATE_ITERATIVE_EXTENT(extent,thr,la,lb,angmom,a,b,p,mu,X,Y,Z,lupri)
  implicit none
  INTEGER,intent(in)    :: la,lb,angmom,lupri
  REAL(REALK),intent(out):: extent
  REAL(REALK),intent(in):: thr,a,b,p,mu,X,Y,Z
  !
  REAL(realk):: XAB,ROLD,fun0,fun1,Rnew,XAB2,n,l,K
  REAL(realk),parameter :: e=2.71828182845904523536028747135266249775724709369995_realk
  REAL(realk),parameter :: D1=1E0_realk,D2=2E0_realk

  XAB2 = X*X+Y*Y+Z*Z
  XAB=SQRT(XAB2)
  IF(ANGMOM .EQ. 0)THEN
     IF(thr .GT. exp(-(a*b)/p*(XAB2)))THEN
        extent = 0.0E0_realk
     ELSE
        extent = sqrt(-(a*b/p*(XAB2) + LOG(thr))/p)
     ENDIF
  ELSE
     IF(thr .GT. exp(-(a*b)/p*(XAB2)))THEN
        extent = 0.0E0_realk
     ELSE
        ROLD = sqrt(-( mu*XAB2 + LOG(thr))/(p))     
        fun0 = ABS((ROLD+b/p*XAB)**lb*(ROLD+a/p*XAB)**la*exp(-mu*XAB2)*exp(-p*ROLD**2))
        IF(thr .LE. fun0)THEN
           DO
              fun1 = ABS(lb*(ROLD+b/p*XAB)**(lb-1)*(ROLD+a/p*XAB)**la*exp(-mu*XAB2)*exp(-p*ROLD**2)+&
                   & (ROLD+b/p*XAB)**lb*la*(ROLD+a/p*XAB)**(la-1)*exp(-mu*XAB2)*exp(-p*ROLD**2)+&
                   & (ROLD+b/p*XAB)**lb*(ROLD+a/p*XAB)**la*exp(-mu*XAB2)*exp(-p*ROLD**2)*(-2*p*ROLD))
              Rnew = ROLD - (thr-fun0)/fun1
              IF(Rnew-ROLD .LE. 1.0E-24_realk)THEN
                 write(lupri,*)'calculation stalled in CALCULATE_ITERATIVE_EXTENT'
                 CALL LSQUIT('calculation stalled in CALCULATE_ITERATIVE_EXTENT',-1)
              ENDIF
              fun0 = ABS((Rnew-b/p*XAB)**lb*(Rnew+a/p*XAB)**la*exp(-mu*XAB2)*exp(-p*Rnew**2))   
              ROLD = Rnew
              IF((fun0-thr).LE. thr*1.0E-10_realk)EXIT
              IF(Rold.NE.Rold)THEN
                 write(lupri,*)'found NaN in CALCULATE_ITERATIVE_EXTENT'
                 CALL LSQUIT('found NaN in CALCULATE_ITERATIVE_EXTENT',lupri)
              ENDIF
           ENDDO
           extent = ROLD
        ELSE
           extent = ROLD
        ENDIF
     ENDIF
  ENDIF

!  n=100
!  l=la+lb
!  K = (((n+D1+l)**((n+D1+l)/(D2*n)))/l**(l/(D2*n)))*(D1/e)**((n+D1)/(D2*n))
!  WRITE(lupri,*)'EQ. 48=',K*(1/(2*p))**((n+D1)/(D2*n))
!  extent = K*(1/(2*p))**((n+D1)/(D2*n))

END SUBROUTINE CALCULATE_ITERATIVE_EXTENT

END MODULE ODbatches

!> \brief 
!> \author 
!> \date 
!> \param THR the threshold 
FUNCTION LS_ERFCIV(THR)
 use precision
implicit none
real(realk),PARAMETER :: D0 = 0.0E0_realk,THRMAX = 1.0E-15_realk
real(realk)           :: ERFCI(0:15),LS_ERFCIV,THR
REAL(realk)           :: ARG

  ERFCI = (/0.00000E0_realk, 1.16309E0_realk, 1.82139E0_realk, 2.32675E0_realk, &
     &                       2.75106E0_realk, 3.12341E0_realk, 3.45891E0_realk, 3.76656E0_realk,& 
     &                       4.05224E0_realk, 4.32001E0_realk, 4.57282E0_realk, 4.81292E0_realk,&
     &                       5.04203E0_realk, 5.26151E0_realk, 5.47248E0_realk, 5.67585E0_realk/)
  ARG = MAX(MAX(ABS(THR),THRMAX),D0)
  LS_ERFCIV = ERFCI(-INT(LOG10(ARG)))
END FUNCTION
