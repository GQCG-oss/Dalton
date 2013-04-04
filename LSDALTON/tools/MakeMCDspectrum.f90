program makeMCDspectra
IMPLICIT NONE
INTEGER, PARAMETER :: realk = KIND(1D0)
integer, parameter :: long = selected_int_kind(10)
Real(realk),parameter :: PI=3.14159265358979323846D0
Real(realk),parameter :: BtermConstant=5.88764*1.d-05
Real(realk),parameter :: AtermConstant=12.920947806163007
Real(realk),parameter :: VerdetConstant=5.88764*1.d-05*33.53!0.152123007*1.0d-7
character(len=4) :: STRING4
real(realk) :: eigenvalue(3),eigenvalueL(3),sign,signL,X1,Y1,Z1,X2,Y2,Z2
real(realk),allocatable :: twophoton(:,:,:),Atermmoment(:,:),MCDBterm(:),GAMMAATERM(:)
real(realk),allocatable :: twophotonL(:,:,:),AtermmomentL(:,:),MCDBtermL(:),GAMMABTERM(:)
real(realk),allocatable :: Bmoment(:,:),Amoment(:,:,:),MCDAterm(:),MCDAtermL(:),freqA(:),freqB(:)
integer,allocatable :: BtermRspFunc(:),AtermRspFunc(:)
integer,allocatable :: BtermRspFuncLondon(:),AtermRspFuncLondon(:)
!
INTEGER :: nBterms,LUG,Nstep,lengthX,maxA,maxB,lumcd,lumcd2,nXcoor2,nXcoor3,lumcd3,lumcd4,i,lupri,nAterms
real(realk),allocatable  :: Y(:),YL(:)
real(realk)    :: Range,step,Estart,Eend,FUN,INT,x,DISPLACEMENT,FA
real(realk) :: DTHR
real(realk),allocatable :: PUREB(:),PUREBL(:),PUREA(:),PUREAL(:)
real(realk),allocatable :: Xcoor(:),COMBINED(:),COMBINEDL(:),Deltavec(:),Xcoor2(:),Xcoor3(:)
complex(8),allocatable :: verdet(:),verdetL(:) 
character(len=30)  :: INPUTFILENAME,OUTPUTFILENAME
logical            :: input_exists,LORENTZ,unique
logical            :: degenerateStates,ATOMICUNITS
real(realk),parameter :: TWOPIM1 = 0.15915494309189535d0
real(realk),parameter :: hartree=27.21138386D0
Real(realk), parameter ::eVTOcm= 8065.54D0 !eV to cm^-1
real(realk), parameter :: autocm=hartree*eVTOcm 
character(len=30)  :: FILENAME
character(len=7)   :: TYPE

!=====================================================================
!  Input Reading 
!=====================================================================

WRITE(*,*)' WHICH INPUT FILE SHOUD BE USED (MAX 30 characters)              '
WRITE(*,*)' if you do not know the input format write non existing file name'
READ(*,'(A30)') FILENAME 
INQUIRE(file=FILENAME,EXIST=input_exists)
IF(.NOT.input_exists)THEN
   WRITE(*,*)'NO INPUT FILE WITH THE NAME'
   WRITE(*,*)'-------------------------------'
   WRITE(*,*)FILENAME
   WRITE(*,*)'-------------------------------'
   WRITE(*,*)'THE INPUT SHOULD LOOK LIKE THE FOLLOWING EXAMPLE '
   WRITE(*,'(A)')'  '
   WRITE(*,'(A)')'  '
   WRITE(*,'(A)')'UNITS  '
   WRITE(*,'(A)')'nBTERMS (number of B terms)'
   WRITE(*,'(A)')'freq1         BTERMS1        LAOBTERMS1        LSPARAM1'
   WRITE(*,'(A)')'freq2         BTERMS2        LAOBTERMS2        LSPARAM2'
   WRITE(*,'(A)')'...           ...            ...               ...'        
   WRITE(*,'(A)')'freqnBTERMS   BTERMSnBTERMS  LAOBTERMSnBTERMS  LSPARAMnBTERMS'
   WRITE(*,'(A)')'nATERMS'
   WRITE(*,'(A)')'freq1         ATERMS1        LAOATERMS1        LSPARAM1'
   WRITE(*,'(A)')'freq2         ATERMS2        LAOATERMS2        LSPARAM2'
   WRITE(*,'(A)')'...           ...            ...               ...'
   WRITE(*,'(A)')'freqnATERMS   ATERMSnATERMS  LAOATERMSnATERMS  LSPARAMnBTERMS'
   WRITE(*,'(A)')'nSTEPS'
   WRITE(*,'(A)')'TYPE  '
   WRITE(*,'(A)')'DISPLACEMENT'
   WRITE(*,'(A)')'  '
   WRITE(*,'(A)')'  '
   WRITE(*,'(A)')'UNITS should be either STD (for standard units) or AU (atomic units)  '
   WRITE(*,'(A)')'  '
   WRITE(*,'(A)')'In the case of STD:'
   WRITE(*,'(A)')'All frequences should be written in units of inverse cm.'
   WRITE(*,'(A)')'LAOBTERMS means values for the B terms where London atomic orbitals have been used'
   WRITE(*,'(A)')'for the calculation.'
   WRITE(*,'(A)')'All B terms should be written in standard units - which means molar ellipticity.'
   WRITE(*,'(A)')'All A terms should be written in standard units - which means molar ellipticity.'
   WRITE(*,'(A)')'The individuel lineshape parameters (LSPARAM) should be written in inverse cm.'
   WRITE(*,'(A)')'nSTEPS is the number of points in the simulated specta (suggestion 5000)'
   WRITE(*,'(A)')'TYPE can be LORENTZ or GAUSS.'
   WRITE(*,'(A)')'DISPLACEMENT should be written in inverse cm and is a option to parallel displace'
   WRITE(*,'(A)')'the entire spectra in order to more easily compare to experiments (suggestion 0.0)'   
   WRITE(*,'(A)')'  '
   WRITE(*,'(A)')'In the case of AU'
   WRITE(*,'(A)')'All frequences should be written in atomic units (hartree)'
   WRITE(*,'(A)')'LAOBTERMS means values for the B terms where London atomic orbitals have been used'
   WRITE(*,'(A)')'for the calculation.'
   WRITE(*,'(A)')'All B terms should be written in atomic units'
   WRITE(*,'(A)')'All A terms should be written in atomic units'
   WRITE(*,'(A)')'The individuel lineshape parameters (LSPARAM) should be written in atomic units'
   WRITE(*,'(A)')'nSTEPS is the number of points in the simulated specta (suggestion 5000)'
   WRITE(*,'(A)')'TYPE can be LORENTZ or GAUSS.'
   WRITE(*,'(A)')'DISPLACEMENT should be written in atomic units and is a option to parallel displace'
   WRITE(*,'(A)')'the entire spectra in order to more easily compare to experiments (suggestion 0.0)'   
ELSE
   LUG=1
   LUPRI=2
   open( unit=LUG, file=FILENAME,status="old")
   open( unit=LUPRI, file="MCDspectra.dat",status="unknown")
   READ(LUG,*) TYPE
   IF(TYPE(1:3).EQ.'STD')THEN
      ATOMICUNITS=.FALSE.
   ELSEIF(TYPE(1:2).EQ.'AU')THEN
      ATOMICUNITS=.TRUE.
   ELSE
      WRITE(*,*)'NO UNITS CHOSEN - DEFAULT IS AU'
      ATOMICUNITS=.TRUE.
   ENDIF
   READ(LUG,'(I4)') nBterms
   IF(nBterms .NE. 0)THEN
      ALLOCATE(freqB(nBterms))
      ALLOCATE(MCDBterm(nBterms))
      ALLOCATE(MCDBtermL(nBterms))
      ALLOCATE(GAMMABTERM(nBterms))
      DO i=1,nBterms
         READ(LUG,*)freqB(i),MCDBterm(i),MCDBtermL(i),GAMMABTERM(i)
      ENDDO
      maxB=nBterms
      WRITE(*,*)'THE FOLLOWING ',nBterms,'B terms wil be used to make a MCD spectra'
      DO i=1,nBterms
         WRITE(*,*)'BTERM(',i,')=',MCDBterm(i) ,MCDBtermL(i)
      ENDDO
   ELSE
      WRITE(*,*)'NO B terms'
   ENDIF

   READ(LUG,'(I4)') nAterms
   IF(nAterms .NE. 0)THEN
      ALLOCATE(freqA(nAterms))
      ALLOCATE(MCDAterm(nAterms))
      ALLOCATE(MCDAtermL(nAterms))
      ALLOCATE(GAMMAATERM(nAterms))
      DO i=1,nAterms
         READ(LUG,*)freqA(i),MCDAterm(i),MCDAtermL(i),GAMMAATERM(i)
      ENDDO
      WRITE(*,*)'THE FOLLOWING ',nBterms,'A terms wil be used to make a MCD spectra'
      DO i=1,nAterms
         WRITE(*,*)'ATERM(',i,')=',MCDAterm(i),MCDAtermL(i)
      ENDDO
   ELSE
      WRITE(*,*)'NO A terms'
   ENDIF

   READ(LUG,*) Nstep   ! = 5000
!   GAMMA = 0.005D0
   READ(LUG,*) TYPE
   IF(TYPE(1:7).EQ.'LORENTZ')THEN
      LORENTZ=.TRUE.
   ELSEIF(TYPE(1:5).EQ.'GAUSS')THEN
      LORENTZ=.FALSE.
   ELSE
      WRITE(*,*)'NO TYPE CHOSEN - DEFAULT IS LORENTZ'
      LORENTZ=.TRUE.
   ENDIF
   READ(LUG,*)DISPLACEMENT
!=====================================================================
! Change input to AU
!=====================================================================
   IF(.NOT.ATOMICUNITS)THEN

      IF(nBterms .NE. 0)THEN
         DO i=1,nBterms
            freqB(i) = freqB(i)/auTOcm
         ENDDO
         DO i=1,nBterms
            GammaBterm(i) = GammaBterm(i)/auTOcm
         ENDDO
         DO i=1,nBterms
            MCDBterm(i) = MCDBterm(i)/BtermConstant
         ENDDO
         DO i=1,nBterms
            MCDBtermL(i) = MCDBtermL(i)/BtermConstant
         ENDDO
      ENDIF
      IF(nAterms .NE. 0)THEN
         DO i=1,nAterms
            freqA(i) = freqA(i)/auTOcm
         ENDDO
         DO i=1,nAterms
            GammaAterm(i) = GammaAterm(i)/auTOcm
         ENDDO
         DO i=1,nAterms
            MCDAterm(i) = MCDAterm(i)/AtermConstant
         ENDDO
         DO i=1,nAterms
            MCDAtermL(i) = MCDAtermL(i)/AtermConstant
         ENDDO
      ENDIF
   ENDIF
!=====================================================================
! Done Input Reading freqA(i) in AU 
!=====================================================================

   call simulateSpectra(6,nBterms,nAterms,MCDBterm,MCDBtermL,&
        & MCDAterm,MCDAtermL,freqA,freqB,LORENTZ,Nstep,&
        & GAMMAATERM,GAMMABTERM)

ENDIF

CONTAINS
subroutine simulateSpectra(lupri,nBterms,nAterms,MCDBterm,MCDBtermL,MCDAterm,&
     &MCDAtermL,freqA,freqB,LORENTZ,Nsteps,GAMMAATERM,GAMMABTERM)
  implicit none           
  !> logical unit number of output 
  integer,intent(in)     :: lupri
  integer,intent(in)     :: Nsteps,nAterms,nBterms
  !> MCD B terms in atomic units
  real(realk),intent(in) :: MCDBterm(nBterms),MCDBtermL(nBterms)
  !> MCD A terms in atomic units
  real(realk),intent(in) :: MCDAterm(nAterms),MCDAtermL(nAterms)
  !> MCD A and B term freqs
  real(realk),intent(in) :: freqA(nAterms),freqB(nBterms)
  real(realk),intent(in) :: GAMMAATERM(nAterms)
  real(realk),intent(in) :: GAMMABTERM(nBterms)
  logical,intent(in)     :: LORENTZ
!
  integer   :: i,j,k,l,LUMCD,LUMCD2,Nstep,lengthX,ios
  real(realk),allocatable  :: Y(:),YL(:)
  real(realk),allocatable :: PUREB(:),PUREBL(:),PUREA(:),PUREAL(:)
  real(realk),allocatable :: Xcoor(:),COMBINED(:),COMBINEDL(:),Deltavec(:),Xcoor2(:),Xcoor3(:)
  complex(realk),allocatable :: verdet(:),verdetL(:) 
  real(realk) :: Range,step,Estart,Eend,FUN,INT,x,DISPLACEMENT,FA
!
  Real(realk),parameter :: PI=3.14159265358979323846E0_realk
  !the Bterm is in standard units given in D**2*muB/(cm^-1) (au : e**3*a0**4/hbar)
  ! D is the unit debye        :   1 D = 0.393430307 (au.) (au : e*a0)
  ! muB is the Bohr magneton   : 1 muB = 0.5 (au.)
  ! and inverse centimeter     : 1cm^⁻1= 219475 (au.)
  ! so the standard units are : D**2*muB/(cm^-1) = 16985.974069713367 (au.)
  ! so to convert from au to standard we inverse this to obtain
  Real(realk),parameter :: BtermConstant=5.88764*1E-05_realk
  !the Aterm is in standard units given in D**2*muB (au : e**3*hbar*a0**2/me)
  !so the standard units are : D**2*muB = 0.07739370323305711 (au.)
  !so to convert from au to standard we inverse this to obtain
  Real(realk),parameter :: AtermConstant=12.920947806163007
  Real(realk),parameter :: VerdetConstant=5.88764*1E-05_realk*33.53!0.152123007*1.0E-7_realk
  !hc = 1239.8419300923944 eV*nm  
  !so  lambda (in nm) = hc/E (E in eV)
  !and E (in eV) = hc/lambda (lambda in nm)
  Real(realk),parameter :: hcConstant_eVnm=1239.8419300923944
  real(realk),parameter :: TWOPIM1 = 0.15915494309189535E0_realk
  real(realk),parameter :: hartree=27.21138386E0_realk
  Real(realk), parameter ::eVTOcm= 8065.54E0_realk !eV to cm^-1
  real(realk), parameter :: autocm=hartree*eVTOcm 

  !#################################################################
  !#   Atomic units UNITS
  !#################################################################
  print*,'freqA',freqA
  print*,'freqB',freqB
  print*,'GAMMAATERM',GAMMAATERM
  print*,'GAMMABTERM',GAMMABTERM
  print*,'MCDBterm',MCDBterm
  print*,'MCDAterm',MCDAterm
  print*,'MCDBtermL',MCDBtermL
  print*,'MCDAtermL',MCDAtermL

  IF(nAterms.NE. 0 .AND. nBterms.NE. 0)then
     Eend=MAX(freqA(nAterms),freqB(nBterms))+0.1*MAX(freqA(nAterms),freqB(nBterms))
  ELSE
     IF(nAterms.NE. 0)THEN
        Eend=freqA(nAterms)+0.1*freqA(nAterms)
     ELSEIF(nBterms.NE. 0)THEN
        Eend = freqB(nBterms)+0.1*freqB(nBterms)
     ELSE
        Eend = -1
     ENDIF
  ENDIF
  print*,'Eend',eEnd
  IF(nAterms.NE. 0)THEN
     IF(nBterms.NE. 0)THEN
        Estart = MAX(MIN(freqA(1),freqB(1))-0.1*MAX(freqA(nAterms),freqB(nBterms)),0E0_realk)
     ELSE
        Estart = MAX(freqA(1)-0.1*MAX(freqA(nAterms),freqB(nBterms)),0E0_realk)
     ENDIF
  ELSE
     IF(nBterms.NE. 0)THEN
        Estart = MAX(freqB(1)-0.1*freqB(nBterms),0E0_realk)
     ELSE
        Estart = 0
     ENDIF
  ENDIF
  print*,'Estart',Estart
  IF(Eend.NE.-1)THEN
     Nstep = Nsteps
     Range = Eend-Estart
     step = Range/Nstep
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
     ENDDO
     Allocate(PUREB(lengthX))
     Allocate(PUREBL(lengthX))
     PUREB = 0E0_realk
     PUREBL = 0E0_realk
     Allocate(PUREA(lengthX))
     Allocate(PUREAL(lengthX))
     PUREA = 0E0_realk
     PUREAL = 0E0_realk
     Allocate(Xcoor(lengthX))
     Xcoor = 0E0_realk
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
        XCOOR(lengthX)=x
     ENDDO
     Allocate(COMBINED(lengthX))
     Allocate(COMBINEDL(lengthX))
     COMBINED = 0E0_realk
     COMBINEDL = 0E0_realk
     DISPLACEMENT = 0E0_realk          
     !#################################################################
     !#    THE Atomic unit B TERMS   
     !#################################################################
     IF(nBterms .NE. 0)THEN
        Allocate(Y(nBterms))
        Allocate(YL(nBterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=GAMMABTERM(i)/((x-freqB(i))**2+GAMMABTERM(i)**2)
                 INT=INT+(FUN/x)*step
              ENDDO
              Y(i)=0.0
              YL(i)=0.0
              IF(ABS(MCDBterm(i)).GT.1.0d-16) Y(i)=-MCDBterm(i)/INT 
              IF(ABS(MCDBtermL(i)).GT.1.0d-16) YL(i)=-MCDBtermL(i)/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 IF(ABS(MCDBterm(i)).GT.1.0d-16) PUREB(lengthX)=Y(i)*GAMMABTERM(i)/((x-freqB(i))**2+GAMMABTERM(i)**2) 
                 IF(ABS(MCDBtermL(i)).GT.1.0d-16) PUREBL(lengthX)=YL(i)*GAMMABTERM(i)/((x-freqB(i))**2+GAMMABTERM(i)**2) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX)  
              ENDDO
           ENDDO
        ELSE !GAUSSIAN
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=1/(gammaBTERM(i)*sqrt(2*PI))*EXP(-((x-freqB(i))**2)/(2*gammaBTERM(i)*gammaBTERM(i))) 
                 INT=INT+(FUN/x)*step
              ENDDO
              !PLOT IN ATOMIC UNIT THEN ELSE
              IF(ABS(MCDBterm(i)).GT.1.0d-16) Y(i)=-MCDBterm(i)/INT 
              IF(ABS(MCDBtermL(i)).GT.1.0d-16) YL(i)=-MCDBtermL(i)/INT
!              Y(i)=-MCDBterm(i)/INT 
!              YL(i)=-MCDBtermL(i)/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 IF(ABS(MCDBterm(i)).GT.1.0d-16) &
                 & PUREB(lengthX)=Y(i)/(gammaBTERM(i)*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i))**2)/(2*gammaBTERM(i)*gammaBTERM(i))) 
                 IF(ABS(MCDBtermL(i)).GT.1.0d-16) &
                      & PUREBL(lengthX)=YL(i)/(gammaBTERM(i)*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i))**2)/(2*gammaBTERM(i)*gammaBTERM(i))) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    THE Atomic unit A TERMS   
     !#################################################################
     IF(nAterms .NE. 0)THEN
        Allocate(Y(nAterms))
        Allocate(YL(nAterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=(2E0_realk*GAMMAATERM(i)*(x-freqA(i)))/((x-freqA(i))**2+GAMMAATERM(i)**2)**2
                 INT=INT+(x-freqA(i))*(FUN/x)*step
              ENDDO
              IF(ABS(MCDAterm(i)).GT.1.0d-16) Y(i)=MCDAterm(i)/INT 
              IF(ABS(MCDAtermL(i)).GT.1.0d-16) YL(i)=MCDAtermL(i)/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 IF(ABS(MCDAtermL(i)).GT.1.0d-16) &
                      & PUREA(lengthX)  =Y(i)*(2E0_realk*GAMMAATERM(i)*(x-freqA(i)))/((x-freqA(i))**2+GAMMAATERM(i)**2)**2
                 IF(ABS(MCDAtermL(i)).GT.1.0d-16) &
                      & PUREAL(lengthX)=YL(i)*(2E0_realk*GAMMAATERM(i)*(x-freqA(i)))/((x-freqA(i))**2+GAMMAATERM(i)**2)**2
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX) + PUREAL(lengthX)  
              ENDDO
           ENDDO
        ELSE !GAUSSIAN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=((x-freqA(i))/(gammaATERM(i)*gammaATERM(i)))/(gammaATERM(i)*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i))**2)/(2*gammaATERM(i)*gammaATERM(i))) 
                 INT=INT+(x-freqA(i))*(FUN/x)*step
              ENDDO
              IF(ABS(MCDAterm(i)).GT.1.0d-16) Y(i)=MCDAterm(i)/INT 
              IF(ABS(MCDAtermL(i)).GT.1.0d-16) YL(i)=MCDAtermL(i)/INT 
!              Y(i)=MCDAterm(i)/INT
!              YL(i)=MCDAtermL(i)/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 IF(ABS(MCDAtermL(i)).GT.1.0d-16) &
                      & PUREA(lengthX)=Y(i)*((x-freqA(i))/(gammaATERM(i)*gammaATERM(i)))/(gammaATERM(i)*sqrt(2*PI))*&
                      & EXP(-((x-freqA(i))**2)/(2*gammaATERM(i)*gammaATERM(i))) 
                 IF(ABS(MCDAtermL(i)).GT.1.0d-16) &
                      & PUREAL(lengthX)=YL(i)*((x-freqA(i))/(gammaATERM(i)*gammaATERM(i)))/(gammaATERM(i)*sqrt(2*PI))*&
                      & EXP(-((x-freqA(i))**2)/(2*gammaATERM(i)*gammaATERM(i))) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX)+ PUREAL(lengthX)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    PRINT THE Atomic unit MCD SPECTRA   
     !#################################################################
     LUMCD=23
     OPEN(UNIT=LUMCD,FILE='MCDspectraAU.dat',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
     DO I=1,lengthX
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREB(1,I),PUREB(2,I),PUREB(3,I),COMBINED(I)
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREA(1,I),COMBINED(I)
        WRITE(LUMCD,*)Xcoor(I),COMBINEDL(I),COMBINED(I)
     ENDDO
     close(LUMCD,STATUS='KEEP',IOSTAT=IOS)

     LUMCD=24
     OPEN(UNIT=LUMCD,FILE='BtermAU.dat',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
     DO I=1,lengthX
        WRITE(LUMCD,*)Xcoor(I),PUREBL(I),PUREB(I)
     ENDDO
     close(LUMCD,STATUS='KEEP',IOSTAT=IOS)

     LUMCD=25
     OPEN(UNIT=LUMCD,FILE='AtermAU.dat',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
     DO I=1,lengthX
        WRITE(LUMCD,*)Xcoor(I),PUREAL(I),PUREA(I)
     ENDDO
     close(LUMCD,STATUS='KEEP',IOSTAT=IOS)

     LUMCD2=26
     OPEN(UNIT=LUMCD2,FILE='gnuplot_MCDAU.gnu',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
     WRITE(LUMCD2,'(A)')'reset'
     WRITE(LUMCD2,'(A)')'set title '//Achar(39)//'MCD Spectra with london orbitals'//Achar(39)
     WRITE(LUMCD2,'(A)')'set terminal postscript eps enhanced color'
     WRITE(LUMCD2,'(A)')'set output ''simulatedMCDspectraAU.ps'' '
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A)')'#set linestyle  1' 
     WRITE(LUMCD2,'(A)')'#set linestyle  2'
     WRITE(LUMCD2,'(A)')'#set linestyle  7'
     WRITE(LUMCD2,'(A)')'set size 2.0,1.0'
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A,F8.4,A,F8.4,A)')'set xrange [',Xcoor(1),':',Xcoor(lengthX),']'
     WRITE(LUMCD2,'(A)')'P(x)=0.000'
     WRITE(LUMCD2,'(A)')'plot ''MCDspectraAU.dat'' using 1:2 title ''london'' w l ls 1,'//Achar(92)
     WRITE(LUMCD2,'(A)')'     ''MCDspectraAU.dat'' using 1:3 title ''nolondon'' w l ls 2,'//Achar(92)
     WRITE(LUMCD2,'(TL11,A)')'     P(x) w l ls 7'
     CLOSE(UNIT=LUMCD2,STATUS='KEEP')

     deAllocate(PUREB)
     deAllocate(PUREBL)
     deAllocate(PUREA)
     deAllocate(PUREAL)
     deAllocate(Xcoor)
     deAllocate(COMBINED)
     deAllocate(COMBINEDL)

  endif

  !#################################################################
  !#   Standard UNITS which means CM^-1 for the freq and 
  !#   molar ellipticity for the ellipticity
  !#################################################################
  IF(nAterms.NE. 0 .AND. nBterms.NE. 0)then
     Eend=MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm)+0.1*MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm)
  ELSE
     IF(nAterms.NE. 0)THEN
        Eend=freqA(nAterms)*auTOcm+0.1*freqA(nAterms)*auTOcm
     ELSEIF(nBterms.NE. 0)THEN
        Eend = freqB(nBterms)*auTOcm+0.1*freqB(nBterms)*auTOcm
     ELSE
        Eend = -1
     ENDIF
  ENDIF
  IF(nAterms.NE. 0)THEN
     IF(nBterms.NE. 0)THEN
        Estart = MAX(MIN(freqA(1)*auTOcm,freqB(1)*auTOcm)-0.1*MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm),0E0_realk)
     ELSE
        Estart = MAX(freqA(1)*auTOcm-0.1*MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm),0E0_realk)
     ENDIF
  ELSE
     IF(nBterms.NE. 0)THEN
        Estart = MAX(freqB(1)*auTOcm-0.1*freqB(nBterms)*auTOcm,0E0_realk)
     ELSE
        Estart = 0
     ENDIF
  ENDIF
  IF(Eend.NE.-1)THEN
     Range = Eend-Estart
     step = Range/Nstep
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
     ENDDO
     Allocate(PUREB(lengthX))
     Allocate(PUREBL(lengthX))
     PUREB = 0E0_realk
     PUREBL = 0E0_realk
     Allocate(PUREA(lengthX))
     Allocate(PUREAL(lengthX))
     PUREA = 0E0_realk
     PUREAL = 0E0_realk
     Allocate(Xcoor(lengthX))
     Xcoor = 0E0_realk
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
        XCOOR(lengthX)=x
     ENDDO
     Allocate(COMBINED(lengthX))
     Allocate(COMBINEDL(lengthX))
     COMBINED = 0E0_realk
     COMBINEDL = 0E0_realk
     DISPLACEMENT = 0E0_realk
     !#################################################################
     !#    THE B TERMS   
     !#################################################################
     IF(nBterms .NE. 0)THEN
        Allocate(Y(nBterms))
        Allocate(YL(nBterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=GAMMABTERM(i)*auTOcm/((x-freqB(i)*auTOcm)**2+(GAMMABTERM(i)*auTOcm)**2)
                 INT=INT+(FUN/x)*step
              ENDDO
              Y(i)=-BtermConstant*MCDBterm(i)*33.53/INT 
              YL(i)=-BtermConstant*MCDBtermL(i)*33.53/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREB(lengthX)=Y(i)*GAMMABTERM(i)*auTOcm/((x-freqB(i)*auTOcm)**2+(GAMMABTERM(i)*auTOcm)**2) 
                 PUREBL(lengthX)=YL(i)*GAMMABTERM(i)*auTOcm/((x-freqB(i)*auTOcm)**2+(GAMMABTERM(i)*auTOcm)**2) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX)  
              ENDDO
           ENDDO
        ELSE !gaussian
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=1/(gammaBTERM(i)*auTOcm*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i)*auTOcm)**2)/(2*(gammaBTERM(i)*gammaBTERM(i)*auTOcm*auTOcm))) 
                 INT=INT+(FUN/x)*step
              ENDDO
              !PLOT IN ATOMIC UNIT THEN ELSE
              Y(i)=-BtermConstant*MCDBterm(i)*33.53/INT 
              YL(i)=-BtermConstant*MCDBtermL(i)*33.53/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREB(lengthX)=Y(i)/(gammaBTERM(i)*auTOcm*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i)*auTOcm)**2)/(2*gammaBTERM(i)*auTOcm*gammaBTERM(i)*auTOcm)) 
                 PUREBL(lengthX)=YL(i)/(gammaBTERM(i)*auTOcm*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i)*auTOcm)**2)/(2*gammaBTERM(i)*auTOcm*gammaBTERM(i)*auTOcm)) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    THE A TERMS   
     !#################################################################
     IF(nAterms .NE. 0)THEN
        Allocate(Y(nAterms))
        Allocate(YL(nAterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=(2E0_realk*GAMMAATERM(i)*auTOcm*(x-freqA(i)*auTOcm))/((x-freqA(i)*auTOcm)**2+(GAMMAATERM(i)*auTOcm)**2)**2
                 INT=INT+(x-freqA(i)*auTOcm)*(FUN/x)*step
              ENDDO
              Y(i)=AtermConstant*MCDAterm(i)*33.53/INT 
              YL(i)=AtermConstant*MCDAtermL(i)*33.53/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREA(lengthX)  = Y(i)*(2E0_realk*GAMMAATERM(i)*auTOcm*(x-freqA(i)*auTOcm))/((x-freqA(i)*auTOcm)**2+(GAMMAATERM(i)*auTOcm)**2)**2
                 PUREAL(lengthX) =YL(i)*(2E0_realk*GAMMAATERM(i)*auTOcm*(x-freqA(i)*auTOcm))/((x-freqA(i)*auTOcm)**2+(GAMMAATERM(i)*auTOcm)**2)**2
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX) + PUREAL(lengthX)  
              ENDDO
           ENDDO
        ELSE !GAUSSIAN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=((x-freqA(i)*auTOcm)/(gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm))/(gammaATERM(i)*auTOcm*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i)*auTOcm)**2)/(2*gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm)) 
                 INT=INT+(x-freqA(i)*auTOcm)*(FUN/x)*step
              ENDDO
              Y(i)=AtermConstant*MCDAterm(i)*33.53/INT
              YL(i)=AtermConstant*MCDAtermL(i)*33.53/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREA(lengthX) = Y(i)*((x-freqA(i)*auTOcm)/(gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm))/(gammaATERM(i)*auTOcm*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i)*auTOcm)**2)/(2*gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm)) 
                 PUREAL(lengthX)=YL(i)*((x-freqA(i)*auTOcm)/(gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm))/(gammaATERM(i)*auTOcm*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i)*auTOcm)**2)/(2*gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm)) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX)+ PUREAL(lengthX)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    PRINT THE MCD SPECTRA   
     !#################################################################
     LUMCD=123
     OPEN(UNIT=LUMCD,FILE='MCDspectra.dat',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
!     CALL LSOPEN(LUMCD,'MCDspectra.dat','REPLACE','FORMATTED')
     DO I=1,lengthX
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREB(1,I),PUREB(2,I),PUREB(3,I),COMBINED(I)
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREA(1,I),COMBINED(I)
        WRITE(LUMCD,*)Xcoor(I),COMBINEDL(I),COMBINED(I)
     ENDDO
     close(LUMCD,STATUS='KEEP',IOSTAT=IOS)
!     CALL LSCLOSE(LUMCD,'KEEP')

     LUMCD=124
     OPEN(UNIT=LUMCD,FILE='Bterm.dat',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
!     CALL LSOPEN(LUMCD,'Bterm.dat','REPLACE','FORMATTED')
     DO I=1,lengthX
        WRITE(LUMCD,*)Xcoor(I),PUREBL(I),PUREB(I)
     ENDDO
     close(LUMCD,STATUS='KEEP',IOSTAT=IOS)
!     CALL LSCLOSE(LUMCD,'KEEP')

     LUMCD=125
     OPEN(UNIT=LUMCD,FILE='Aterm.dat',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
!     CALL LSOPEN(LUMCD,'Aterm.dat','REPLACE','FORMATTED')
     DO I=1,lengthX
        WRITE(LUMCD,*)Xcoor(I),PUREAL(I),PUREA(I)
     ENDDO
     close(LUMCD,STATUS='KEEP',IOSTAT=IOS)
!     CALL LSCLOSE(LUMCD,'KEEP')

     LUMCD2=126
     OPEN(UNIT=LUMCD2,FILE='gnuplot_MCD.gnu',STATUS='REPLACE',FORM='FORMATTED',IOSTAT=IOS)
!     CALL LSOPEN(LUMCD2,'gnuplot_MCD.gnu','REPLACE','FORMATTED')
     WRITE(LUMCD2,'(A)')'reset'
     WRITE(LUMCD2,'(A)')'set title '//Achar(39)//'MCD Spectra with london orbitals'//Achar(39)
     WRITE(LUMCD2,'(A)')'set terminal postscript eps enhanced color'
     WRITE(LUMCD2,'(A)')'set output ''simulatedMCDspectra.ps'' '
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A)')'#set linestyle  1' 
     WRITE(LUMCD2,'(A)')'#set linestyle  2'
     WRITE(LUMCD2,'(A)')'#set linestyle  7'
     WRITE(LUMCD2,'(A)')'set size 2.0,1.0'
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A,F9.2,A,F9.2,A)')'set xrange [',Xcoor(1),':',Xcoor(lengthX),']'
     WRITE(LUMCD2,'(A)')'P(x)=0.000'
     WRITE(LUMCD2,'(A)')'plot ''MCDspectra.dat'' using 1:2 title ''london'' w l ls 1,'//Achar(92)
     WRITE(LUMCD2,'(A)')'     ''MCDspectra.dat'' using 1:3 title ''nolondon'' w l ls 2,'//Achar(92)
     WRITE(LUMCD2,'(TL11,A)')'     P(x) w l ls 7'
     close(LUMCD2,STATUS='KEEP',IOSTAT=IOS)
!     CALL LSCLOSE(LUMCD2,'KEEP')
     deAllocate(PUREB)
     deAllocate(PUREBL)
     deAllocate(PUREA)
     deAllocate(PUREAL)
     deAllocate(Xcoor)
     deAllocate(COMBINED)
     deAllocate(COMBINEDL)
  endif
end subroutine simulateSpectra

END program makeMCDspectra
