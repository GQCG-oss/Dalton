!> Module containing soubroutine for calculation of the exchange-correlation contribution to KS-matrix
MODULE IIDFTKSMWORK
use precision
!use TYPEDEF
use dft_type
use LS_UTIL,only: DGEMM_TS
use dft_memory_handling
!WARNING you must not add memory_handling, all memory goes through 
!grid_memory_handling  module so as to determine the memory used in this module.
#ifdef VAR_XCFUN
use xcfun_host
#endif

logical,save :: XCintNoOMP

CONTAINS
SUBROUTINE SetNoOMP(InputNoOMP)
implicit none
logical,intent(in) :: InputNoOMP
XCintNoOMP = InputNoOMP
END SUBROUTINE SetNoOMP

SUBROUTINE DFT_DOGGA_DOMETA(DOGGA,DOMETA)
implicit none
LOGICAL,intent(inout) :: DOGGA,DOMETA
LOGICAL,EXTERNAL :: DFT_ISGGA
#ifdef VAR_XCFUN
IF(USEXCFUN)THEN
   call xcfun_host_type(DOGGA,DOMETA)
ELSE
   DOGGA  = DFT_ISGGA()
   DOMETA = .FALSE.
ENDIF
#else
   DOGGA  = DFT_ISGGA()
   DOMETA = .FALSE.
#endif
END SUBROUTINE DFT_DOGGA_DOMETA

!==============================================================================
!
! Worker routines called from II_DFTINT (the CB function name)
!
!==============================================================================

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Worker routine that for a batch of gridpoints build the LDA kohn-sham matrix
!>
SUBROUTINE II_DFT_KSMLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                   RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                   GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Max Number of active basis functions
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK) :: VX(5),DFTENE
INTEGER     :: IPNT,I,J,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk
REAL(REALK),pointer :: VXC(:,:)
#ifdef VAR_XCFUN
REAL(REALK),pointer :: XCFUNINPUT(:,:),XCFUNOUTPUT(:,:)
INTEGER,pointer     :: IFULL(:)
INTEGER             :: NPNT
#endif
EXTERNAL DFTENE
call mem_dft_alloc(VXC,NBLEN,NDMAT)
! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IDMAT = 1, NDMAT
#ifdef VAR_XCFUN
IF(.NOT.USEXCFUN)THEN
#endif
 DO IPNT = 1, NBLEN
   IF(RHO(IPNT,IDMAT) .GT. RHOTHR)THEN
         !get the functional derivatives 
         !vx(1) = drvs.df1000   = \frac{\partial  f}{\partial \rho_{\alpha}}
         CALL dft_funcderiv1(RHO(IPNT,IDMAT),DUMMY,WGHT(IPNT),VX)
         DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + DFTENE(RHO(IPNT,IDMAT),DUMMY)*WGHT(IPNT)
         !WARNING the factor 2 is due to F=F_alpha + F_beta = 2 F_alpha 
         VXC(IPNT,IDMAT) = D2*VX(1) 
   ELSE
      VXC(IPNT,IDMAT) = 0.0E0_realk
   ENDIF
 END DO
#ifdef VAR_XCFUN
ELSE
  CALL xcfun_lda_init(RHO,VXC,MXBLLEN,NBLEN,NPNT,IDMAT,NDMAT,RHOTHR,IFULL,XCFUNINPUT,XCFUNOUTPUT)
  call xcfun_lda_xc_eval(XCFUNINPUT,XCFUNOUTPUT,NPNT)
  DO IPNT=1,NPNT
    DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + XCFUNOUTPUT(1,IPNT)*WGHT(IFULL(IPNT))
    !XCFUNOUTPUT(2,1) = d Exc/d rho
    VXC(IFULL(IPNT),IDMAT) = D2*XCFUNOUTPUT(2,IPNT)*WGHT(IFULL(IPNT))
  ENDDO
  CALL xcfun_lda_free(XCFUNINPUT,XCFUNOUTPUT,IFULL)
ENDIF
#endif
ENDDO
W1 = 1
W2 = NBLEN*Nactbast                        !W1 - 1 + NBLEN*Nactbast    -> GAORED
W3 = NBLEN*Nactbast+1                      !W2 + 1
W4 = (NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast -> EXCRED
W5 = (NBLEN+Nactbast)*Nactbast + 1         !W4 + 1 
W6 = (NBLEN+Nactbast + 1) * Nactbast       !W5 - 1 + Nactbast          -> GAOGMX
W7 = (NBLEN+Nactbast + 1) * Nactbast + 1   !W6 + 1
W8 = (2*NBLEN+Nactbast + 1) * Nactbast     !W7 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
IF(W8.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_KSMLDA',lupri)
#endif
! LDA Exchange-correlation contribution to Kohn-Sham matrix
DO IDMAT = 1, NDMAT
   CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
        & VXC(:,IDMAT),GAO(:,:,1),SHAREDDFTDATA%FKSM(:,:,IDMAT),DFTHRI,WORK(W1:W2),&
        & WORK(W3:W4),WORK(W5:W6),WORK(W7:W8))
ENDDO
call mem_dft_dealloc(VXC)

END SUBROUTINE II_DFT_KSMLDA

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Worker routine that for a batch of gridpoints build the unrestricted LDA kohn-sham matrix
!>
SUBROUTINE II_DFT_KSMLDAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                        RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                        GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk
INTEGER     :: IPNT,I,J,IDMAT,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT1,IDMAT2
REAL(REALK) :: VX(5),DFTENEUNRES
REAL(REALK),pointer :: VXC(:,:)
REAL(REALK) :: XCFUNINPUT(2,1),XCFUNOUTPUT(3,1)
EXTERNAL DFTENEUNRES
call mem_dft_alloc(VXC,NBLEN,NDMAT)

! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IDMAT = 1, NDMAT/2
 IDMAT1 = 1 + (IDMAT-1)*2
 IDMAT2 = 2 + (IDMAT-1)*2
 DO IPNT = 1, NBLEN
  IF(RHO(IPNT,IDMAT1) .GT. RHOTHR .OR. RHO(IPNT,IDMAT2) .GT. RHOTHR)THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000   = \frac{\partial  f}{\partial \rho_{\alpha}} 
      !and same for beta
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         CALL dft_funcderiv1unres(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),DUMMY,DUMMY,WGHT(IPNT),VX)
         DFTDATA%ENERGY(IDMAT1) = DFTDATA%ENERGY(IDMAT1) + DFTENEUNRES(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),DUMMY,DUMMY)*WGHT(IPNT)
         VXC(IPNT,IDMAT1) = VX(1)
         VXC(IPNT,IDMAT2) = VX(2)
#ifdef VAR_XCFUN
      ELSE
         ! rho_alpha = XCFUNINPUT(1,1)
         ! rho_beta = XCFUNINPUT(2,1)
         XCFUNINPUT(1,1) = RHO(IPNT,IDMAT1)
         XCFUNINPUT(2,1) = RHO(IPNT,IDMAT2)
         call xcfun_lda_unres_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
         DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + XCFUNOUTPUT(1,1)*WGHT(IPNT)
         ! XCFUNOUTPUT(1,1) - Exc
         ! XCFUNOUTPUT(2,1) - d Exc/d rho_alpha
         ! XCFUNOUTPUT(3,1) - d Exc/d rho_beta
         VXC(IPNT,IDMAT1) = XCFUNOUTPUT(2,1)*WGHT(IPNT)
         VXC(IPNT,IDMAT2) = XCFUNOUTPUT(3,1)*WGHT(IPNT)
      ENDIF
#endif
   ELSE
      VXC(IPNT,IDMAT1) = 0.0E0_realk
      VXC(IPNT,IDMAT2) = 0.0E0_realk
   ENDIF
 END DO
ENDDO
! LDA Exchange-correlation contribution to Kohn-Sham matrix
W1 = 1
W2 = NBLEN*Nactbast                        !W1 - 1 + NBLEN*Nactbast    -> GAORED 
W3 = NBLEN*Nactbast+1                      !W2 + 1
W4 = (NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast -> EXCRED
W5 = (NBLEN+Nactbast)*Nactbast + 1         !W4 + 1 
W6 = (NBLEN+Nactbast + 1) * Nactbast       !W5 - 1 + Nactbast          -> GAOGMX
W7 = (NBLEN+Nactbast + 1) * Nactbast + 1   !W6 + 1
W8 = (2*NBLEN+Nactbast + 1) * Nactbast     !W7 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
IF(W8.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_KSMLDAUNRES',lupri)
#endif
DO IDMAT = 1, NDMAT 
   CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
        & VXC(:,IDMAT),GAO(:,:,1),SHAREDDFTDATA%FKSM(:,:,IDMAT),DFTHRI,WORK(W1:W2),WORK(W3:W4),&
        & WORK(W5:W6),WORK(W7:W8))
!   CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
!        & VXC(:,2),GAO(:,:,1),SHAREDDFTDATA%FKSM(:,:,2),DFTHRI,WORK(W1:W2),WORK(W3:W4),&
!        & WORK(W5:W6),WORK(W7:W8))
ENDDO
call mem_dft_dealloc(VXC)
END SUBROUTINE II_DFT_KSMLDAUNRES

!> \brief main closed shell GGA kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DFT_KSMGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                   RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                   GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D8 = 8.0E0_realk,D4 = 4.0E0_realk 
INTEGER     :: IPNT,I,J,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT
REAL(REALK) :: VX(3),DFTENE,GRD,GRDA,A
REAL(REALK),pointer :: VXC(:,:,:)
REAL(REALK) :: XCFUNINPUT(4,1),XCFUNOUTPUT(5,1)
EXTERNAL DFTENE
call mem_dft_alloc(VXC,4,NBLEN,NDMAT)
!     GGA Exchange-correlation contribution to Kohn-Sham matrix
DO IDMAT=1,NDMAT
 DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,IDMAT)*GRAD(1,IPNT,IDMAT)+GRAD(2,IPNT,IDMAT)*GRAD(2,IPNT,IDMAT)&
        &+GRAD(3,IPNT,IDMAT)*GRAD(3,IPNT,IDMAT))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,IDMAT).GT.RHOTHR) THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000 = \frac{\partial f}{\partial \rho_{\alpha}}        
      !vx(2) = drvs.df0010 = \frac{\partial f}{\partial |\nabla \rho_{\alpha}|^{2}}
      !vx(3) = drvs.df00001= \frac{\partial f}{\partial \nabla \rho_{\alpha} \nabla \rho_{\beta}}  
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         CALL dft_funcderiv1(RHO(IPNT,IDMAT),GRD,WGHT(IPNT),VX)
         IF(DFTDATA%LB94)THEN
            CALL LB94correction(rho(IPNT,IDMAT),GRD,DFTDATA%HFexchangeFac,&
                 & WGHT(IPNT),VX(1))
         ELSEIF(DFTDATA%CS00)THEN
            CALL CS00correction(rho(IPNT,IDMAT),GRD,DFTDATA%HFexchangeFac,&
                 & WGHT(IPNT),VX(1),DFTDATA%CS00SHIFT,DFTDATA%CS00eHOMO,DFTDATA%CS00ZND1,DFTDATA%CS00ZND2)
         ENDIF
         DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + DFTENE(RHO(IPNT,IDMAT),GRD)*WGHT(IPNT)
         VXC(1,IPNT,IDMAT) = D2*VX(1) !this is 2 times greater than the UNRES VXC BECAUSE we need alpha + beta
         IF(GRD.GT. 1E-40_realk) THEN
            GRDA = D05*GRD
            !         VXC(2,IPNT) = D2*(VX(2)/GRDA + VX(3))  
            A = D2*(VX(2)/GRDA + VX(3))
            VXC(2,IPNT,IDMAT) = A*GRAD(1,IPNT,IDMAT)
            VXC(3,IPNT,IDMAT) = A*GRAD(2,IPNT,IDMAT)
            VXC(4,IPNT,IDMAT) = A*GRAD(3,IPNT,IDMAT)
            !WARNING: this is the same as for the unrestricted case, but since we use \nabla rho 
            !and not \nabla rho_{\alpha} in II_DISTGGA it is the same as having a factor 2 on this 
            !meaning that we take alpha + beta.
         ELSE
            VXC(2,IPNT,IDMAT) = 0E0_realk
            VXC(3,IPNT,IDMAT) = 0E0_realk
            VXC(4,IPNT,IDMAT) = 0E0_realk
         ENDIF
#ifdef VAR_XCFUN
      ELSE
         XCFUNINPUT(1,1) = RHO(IPNT,IDMAT)
         XCFUNINPUT(2,1) = GRAD(1,IPNT,IDMAT)
         XCFUNINPUT(3,1) = GRAD(2,IPNT,IDMAT)
         XCFUNINPUT(4,1) = GRAD(3,IPNT,IDMAT)
         ! Input:
         !rho   = XCFUNINPUT(1,1)
         !grad_x = XCFUNINPUT(2,1)
         !grad_y = XCFUNINPUT(3,1)
         !grad_z = XCFUNINPUT(4,1)
         call xcfun_gga_components_xc_single_eval(XCFUNINPUT,5,XCFUNOUTPUT,1)
         ! Output
         ! Order 0
         ! out(1,1) Exc
         ! Order 1
         ! out(2,1) d^1 Exc / d rho
         ! out(3,1) d^1 Exc / d grad_x
         ! out(4,1) d^1 Exc / d grad_y
         ! out(5,1) d^1 Exc / d grad_z

         DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + XCFUNOUTPUT(1,1)*WGHT(IPNT)
         IF(DFTDATA%LB94)THEN
            call lsquit('error lb94 not implemented for xcfun',-1)
         ELSEIF(DFTDATA%CS00)THEN
            call lsquit('error cs00 not implemented for xcfun',-1)
         ENDIF
         !the \Omega_{\mu \nu} part
         VXC(1,IPNT,IDMAT) = D2*XCFUNOUTPUT(2,1)*WGHT(IPNT)
         VXC(2,IPNT,IDMAT) = D4*XCFUNOUTPUT(3,1)*WGHT(IPNT)
         VXC(3,IPNT,IDMAT) = D4*XCFUNOUTPUT(4,1)*WGHT(IPNT)
         VXC(4,IPNT,IDMAT) = D4*XCFUNOUTPUT(5,1)*WGHT(IPNT)
      ENDIF
#endif
   ELSE
      VXC(1,IPNT,IDMAT) = 0E0_realk
      VXC(2,IPNT,IDMAT) = 0E0_realk
      VXC(3,IPNT,IDMAT) = 0E0_realk
      VXC(4,IPNT,IDMAT) = 0E0_realk
   END IF
 END DO
ENDDO
W1 = 1
W2 = NBLEN*Nactbast*4                        !W1 - 1 + NBLEN*Nactbast*4 -> GAORED 
W3 = 4*NBLEN*Nactbast+1                      !W2 + 1
W4 = (4*NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast  -> EXCRED
W5 = (4*NBLEN+Nactbast)*Nactbast + 1         !W4 + 1 
W6 = (5*NBLEN+Nactbast)*Nactbast             !W5 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
IF(W6.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_KSMLDA',lupri)
#endif

CALL II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,NactBast,NBAST,NTYPSO,NDMAT,&
     & VXC,GAO,SHAREDDFTDATA%FKSM,DFTHRI,WORK(W1:W2),WORK(W3:W4),GAOGMX,GAOMAX,WORK(W5:W6),MaxNactbast)
call mem_dft_dealloc(VXC)

END SUBROUTINE II_DFT_KSMGGA

!> \brief main closed shell meat-GGA kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DFT_KSMMETA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                    RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                    GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D8 = 8.0E0_realk
INTEGER     :: IPNT,I,J,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT
REAL(REALK) :: VX(3),DFTENE,GRD,GRDA,A
REAL(REALK),pointer :: VXC(:,:,:)
REAL(REALK) :: XCFUNINPUT(3,1),XCFUNOUTPUT(4,1)
EXTERNAL DFTENE
call mem_dft_alloc(VXC,5,NBLEN,NDMAT)
!     GGA Exchange-correlation contribution to Kohn-Sham matrix
DO IDMAT=1,NDMAT
 DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,IDMAT)*GRAD(1,IPNT,IDMAT)+GRAD(2,IPNT,IDMAT)*GRAD(2,IPNT,IDMAT)&
        &+GRAD(3,IPNT,IDMAT)*GRAD(3,IPNT,IDMAT))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,IDMAT).GT.RHOTHR .OR. TAU(IPNT,IDMAT).GT.RHOTHR) THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000 = \frac{\partial f}{\partial \rho_{\alpha}}        
      !vx(2) = drvs.df0010 = \frac{\partial f}{\partial |\nabla \rho_{\alpha}|^{2}}
      !vx(3) = drvs.df00001= \frac{\partial f}{\partial \nabla \rho_{\alpha} \nabla \rho_{\beta}}  
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
!         CALL LSQUIT('meta-GGA only implemented using XC-fun',-1)
         CALL LSQUIT('meta-GGA not implemented',-1)
#ifdef VAR_XCFUN
      ELSE
         call lsquit('xcfun version of II_DFT_KSMMETA not implemented',-1)
      ENDIF
#endif
   ELSE
      VXC(1,IPNT,IDMAT) = 0E0_realk
      VXC(2,IPNT,IDMAT) = 0E0_realk
      VXC(3,IPNT,IDMAT) = 0E0_realk
      VXC(4,IPNT,IDMAT) = 0E0_realk
      VXC(5,IPNT,IDMAT) = 0E0_realk
   END IF
 END DO
ENDDO
W1 = 1
W2 = NBLEN*Nactbast*4                        !W1 - 1 + NBLEN*Nactbast*4 -> GAORED 
W3 = 4*NBLEN*Nactbast+1                      !W2 + 1
W4 = (4*NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast  -> EXCRED
W5 = (4*NBLEN+Nactbast)*Nactbast + 1         !W4 + 1 
W6 = (4*NBLEN+Nactbast + 1) * Nactbast       !W5 - 1 + Nactbast          -> GAOGMX
W7 = (4*NBLEN+Nactbast + 1) * Nactbast + 1   !W6 + 1
W8 = (5*NBLEN+Nactbast + 1) * Nactbast       !W7 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
IF(W8.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_KSMLDA',lupri)
#endif

CALL II_DISTMETA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,NactBast,NBAST,NTYPSO,NDMAT,&
     & VXC,GAO,SHAREDDFTDATA%FKSM,DFTHRI,WORK(W1:W2),WORK(W3:W4),WORK(W5:W6),WORK(W7:W8))
call mem_dft_dealloc(VXC)

END SUBROUTINE II_DFT_KSMMETA

!> \brief main unrestricted GGA kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DFT_KSMGGAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                        RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                        GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D4 = 4.0E0_realk,D2 = 2.0E0_realk,DUMMY = 0E0_realk
INTEGER     :: IPNT,I,J,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT,IDMAT1,IDMAT2
REAL(REALK) :: VXC(NBLEN,2,2),VXM(NBLEN),VX(5),DFTENEUNRES,GRDA,GRDB
REAL(REALK) :: XCFUNINPUT(5,1),XCFUNOUTPUT(6,1)
EXTERNAL DFTENEUNRES
!     GGA Exchange-correlation contribution to Kohn-Sham matrix
DO IDMAT = 1,NDMAT/2
 IDMAT1 = 1 + (IDMAT-1)*2
 IDMAT2 = 2 + (IDMAT-1)*2
 DO IPNT = 1, NBLEN
   GRDA = SQRT(GRAD(1,IPNT,IDMAT1)*GRAD(1,IPNT,IDMAT1)+GRAD(2,IPNT,IDMAT1)*GRAD(2,IPNT,IDMAT1)+&
        &GRAD(3,IPNT,IDMAT1)*GRAD(3,IPNT,IDMAT1))
   GRDB = SQRT(GRAD(1,IPNT,IDMAT2)*GRAD(1,IPNT,IDMAT2)+GRAD(2,IPNT,IDMAT2)*GRAD(2,IPNT,IDMAT2)+&
        &GRAD(3,IPNT,IDMAT2)*GRAD(3,IPNT,IDMAT2))
   IF((GRDA .GT. RHOTHR .OR. RHO(IPNT,IDMAT1).GT.RHOTHR).OR.&
        &(GRDB .GT. RHOTHR .OR. RHO(IPNT,IDMAT2).GT.RHOTHR))then 
      !get the functional derivatives 
      !vx(1) = drvs.df1000 = \frac{\partial f}{\partial \rho_{\alpha}}        
      !vx(2) = drvs.df0100 = \frac{\partial f}{\partial \rho_{\beta }}        
      !vx(3) = drvs.df0010 = \frac{\partial f}{\partial |\nabla \rho_{\alpha}|}
      !vx(4) = drvs.df0001 = \frac{\partial f}{\partial |\nabla \rho_{\alpha}|}
      !vx(5) = drvs.df00001= \frac{\partial f}{\partial \nabla \rho_{\alpha} \nabla \rho_{\beta}}  
      IF(GRDA.LT. 1E-40_realk) GRDA = 1E-40_realk
      IF(GRDB.LT. 1E-40_realk) GRDB = 1E-40_realk
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         CALL dft_funcderiv1unres(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),GRDA,GRDB,WGHT(IPNT),VX)
         DFTDATA%ENERGY(IDMAT1) = DFTDATA%ENERGY(IDMAT1) + DFTENEUNRES(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),GRDA,GRDB)*WGHT(IPNT)
         VXC(IPNT,1,1) = VX(1)
         VXC(IPNT,1,2) = VX(2)
         !WARNING: The factor of 2 in these next coeffficients comes from the fact 
         !that we only build half of these terms and use a symmetrization at the end to get the full contribution
         VXC(IPNT,2,1) = D2*VX(3)/GRDA  
         VXC(IPNT,2,2) = D2*VX(4)/GRDB
         VXM(IPNT) = D2*VX(5) !mixed derivate
#ifdef VAR_XCFUN
      ELSE
         XCFUNINPUT(1,1) = RHO(IPNT,IDMAT1)
         XCFUNINPUT(2,1) = RHO(IPNT,IDMAT2)
         XCFUNINPUT(3,1) = GRAD(1,IPNT,IDMAT1)*GRAD(1,IPNT,IDMAT1)&
              &+GRAD(2,IPNT,IDMAT1)*GRAD(2,IPNT,IDMAT1)&
              &+GRAD(3,IPNT,IDMAT1)*GRAD(3,IPNT,IDMAT1)
         XCFUNINPUT(4,1) = GRAD(1,IPNT,IDMAT1)*GRAD(1,IPNT,IDMAT2)&
              &+GRAD(2,IPNT,IDMAT1)*GRAD(2,IPNT,IDMAT2)&
              &+GRAD(3,IPNT,IDMAT1)*GRAD(3,IPNT,IDMAT2)
         XCFUNINPUT(5,1) = GRAD(1,IPNT,IDMAT2)*GRAD(1,IPNT,IDMAT2)&
              &+GRAD(2,IPNT,IDMAT2)*GRAD(2,IPNT,IDMAT2)&
              &+GRAD(3,IPNT,IDMAT2)*GRAD(3,IPNT,IDMAT2)
         call xcfun_gga_unres_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
         DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + XCFUNOUTPUT(1,1)*WGHT(IPNT)

         IF(DFTDATA%LB94)THEN
            call lsquit('error lb94 xcfun',-1)
         ELSEIF(DFTDATA%CS00)THEN
            call lsquit('error cs00 xcfun',-1)
         ENDIF
         VXC(IPNT,1,1) =    XCFUNOUTPUT(2,1)*WGHT(IPNT)
         VXC(IPNT,1,2) =    XCFUNOUTPUT(3,1)*WGHT(IPNT)
         VXC(IPNT,2,1) = D4*XCFUNOUTPUT(4,1)*WGHT(IPNT)
         VXC(IPNT,2,2) = D4*XCFUNOUTPUT(6,1)*WGHT(IPNT)
         VXM(IPNT)     = D2*XCFUNOUTPUT(5,1)*WGHT(IPNT)
      ENDIF
#endif
   ELSE
      VXC(IPNT,1,1) = 0E0_realk
      VXC(IPNT,2,1) = 0E0_realk
      VXC(IPNT,1,2) = 0E0_realk
      VXC(IPNT,2,2) = 0E0_realk
      VXM(IPNT) = 0E0_realk
   END IF
 END DO

 W1 = 1
 W2 = NBLEN*Nactbast*4                        !W1 - 1 + NBLEN*Nactbast*4 -> GAORED 
 W3 = 4*NBLEN*Nactbast+1                      !W2 + 1
 W4 = (4*NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast  -> EXCRED
 W5 = (4*NBLEN+Nactbast)*Nactbast + 1         !W4 + 1 
 W6 = (4*NBLEN+Nactbast + 1) * Nactbast       !W5 - 1 + Nactbast          -> GAOGMX
 W7 = (4*NBLEN+Nactbast + 1) * Nactbast + 1   !W6 + 1
 W8 = (5*NBLEN+Nactbast + 1) * Nactbast       !W7 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
 IF(W8.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_KSMLDAUNRES',lupri)
#endif

 ! call with drho_alpha dgradrho_alpha dmixed and gradA, gradB
 CALL II_DISTGGABUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,NactBast,NBAST,&
      &  VXC(:,1,1),VXC(:,2,1),VXM,GAO(:,:,1:4),GRAD(:,:,1),GRAD(:,:,2),SHAREDDFTDATA%FKSM(:,:,IDMAT1),&
      &DFTHRI,WORK(W1:W2),WORK(W3:W4),WORK(W5:W6),WORK(W7:W8))
 ! call with drho_beta dmixed dgradrho_beta and gradA, gradB
 CALL II_DISTGGABUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,NactBast,NBAST,&
      &  VXC(:,1,2),VXM,VXC(:,2,2),GAO(:,:,1:4),GRAD(:,:,1),GRAD(:,:,2),SHAREDDFTDATA%FKSM(:,:,IDMAT2),&
      &DFTHRI,WORK(W1:W2),WORK(W3:W4),WORK(W5:W6),WORK(W7:W8))
ENDDO

END SUBROUTINE II_DFT_KSMGGAUNRES

!> \brief GGA Exchange-correlation contribution to molecular gradient
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_geoderiv_molgrad_worker_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     & NDMAT,DMAT,NTYPSO,GAOS,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,&
     & WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
  IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D8 = 8.0E0_realk, D4 = 4E0_realk
REAL(REALK) :: VXC(NBLEN,1),VX(5),GRDA,GRD,DMAX
INTEGER  :: IPNT,I,J,JBL,IBL,K
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL :: DOGGA
INTEGER     :: KVALS(3,3)
REAL(REALK),pointer :: GAORED(:,:,:),GDRED(:,:,:)
REAL(REALK),pointer :: DRED(:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,NRED,orb2atom(nbast)
INTEGER     :: atom(NACTBAST),iatom,IX,K1,K2,K3,KA,ik,jk
REAL(REALK) :: FRC,GA,GA2,GFS
#ifdef VAR_XCFUN
REAL(REALK),pointer :: XCFUNINPUT(:,:),XCFUNOUTPUT(:,:)
INTEGER,pointer     :: IFULL(:)
INTEGER             :: NPNT
#endif



orb2atom = DFTDATA%orb2atom
KVALS(1:3,1) = (/1, 2, 3/)
KVALS(1:3,2) = (/2, 4, 5/)
KVALS(1:3,3) = (/3, 5, 6/)
#ifdef VAR_XCFUN
IF(.NOT.USEXCFUN)THEN
#endif
GRD = 0.D0
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1).GT.RHOTHR) THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000 = \frac{\partial f}{\partial \rho_{\alpha}}        
         CALL dft_funcderiv1(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
         VXC(IPNT,1) = D2*VX(1) 
   ELSE
      VXC(IPNT,1) = 0E0_realk
   END IF
END DO
#ifdef VAR_XCFUN
ELSE
  CALL xcfun_lda_init(RHO,VXC,MXBLLEN,NBLEN,NPNT,1,1,RHOTHR,IFULL,XCFUNINPUT,XCFUNOUTPUT)
  call xcfun_lda_xc_eval(XCFUNINPUT,XCFUNOUTPUT,NPNT)
  IF(DFTDATA%LB94)THEN
     call lsquit('error lb94 not implemented for xcfun',-1)
  ELSEIF(DFTDATA%CS00)THEN
     call lsquit('error cs00 not implemented for xcfun',-1)
  ENDIF
  DO IPNT=1,NPNT
    !XCFUNOUTPUT(2,1) = d Exc/d rho
    VXC(IFULL(IPNT),1) = D2*XCFUNOUTPUT(2,IPNT)*WGHT(IFULL(IPNT))
  ENDDO
  CALL xcfun_lda_free(XCFUNINPUT,XCFUNOUTPUT,IFULL)
ENDIF
#endif

! Set up maximum density-matrix elements
DMAX = 0.0E0_realk
DO JBL=1, NBLOCKS
   DO IBL=1, NBLOCKS
      DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)        !J is active index
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)  !I is active index
            DMAX = MAX(DMAX,ABS(DMAT(I,J,1)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
! Count reduced number of AO's
NRED = 0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
         NRED = NRED + 1
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT. 0) THEN
   call mem_dft_alloc(DRED,NRED,NRED)
   call mem_dft_alloc(GAORED,NBLEN,NRED,NTYPSO)
   call mem_dft_alloc(GDRED,NBLEN,NRED,NTYPSO)
   ! Set up reduced Gaussian AO's
   IRED = 0
   DO IBL=1, NBLOCKS
      DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
         IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
            IRED = IRED + 1
            INXRED(IRED) = I
            DO J=1,NTYPSO
               DO K = 1, NBLEN
                  GAORED(K,IRED,J)  = GAOS(K,I,J)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ! Set up reduced density-matrix
   DO JRED=1,NRED            !Jred is reduced index
      J = INXRED(JRED)       !J is active index
      DO IRED=1,NRED         !Ired is reduced index
         I = INXRED(IRED)    !I is active index
         DRED(IRED,JRED) = DMAT(I,J,1)
      ENDDO
   ENDDO
   ! Set up reduced coordinate index in gradient
   DO IRED=1,NRED
      I = INXACT(INXRED(IRED)) !I is orbital index
      atom(IRED) = orb2atom(I)
   ENDDO
   
   ! Density-matrix contraction
   IF(XCintNoOMP)THEN
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED,&
           &                 NBLEN,DRED,NRED,0.0E0_realk,GDRED,NBLEN    )
   ELSE !Use Thread Safe version 
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED,&
           &                 NBLEN,DRED,NRED,0.0E0_realk,GDRED,NBLEN    )
   ENDIF
   DO IRED=1,NRED
      iatom = atom(IRED)
      KA = INXRED(IRED)  !KA is active index
      DO IX=1,3
         FRC = 0E0_realk
         DO I = 1, NBLEN
            FRC = FRC + VXC(I,1)*GDRED(I,IRED,1)*GAOS(I,KA,IX+1)
         END DO
         DFTDATA%GRAD(IX,iatom) = DFTDATA%GRAD(IX,iatom) - FRC
      ENDDO ! IX
   ENDDO ! IA
   call mem_dft_dealloc(DRED)
   call mem_dft_dealloc(GAORED)
   call mem_dft_dealloc(GDRED)
ENDIF !NRED GT 0

END SUBROUTINE II_GEODERIV_MOLGRAD_WORKER_LDA

!> \brief GGA Exchange-correlation contribution to molecular gradient
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_geoderiv_molgrad_worker_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     & NDMAT,DMAT,NTYPSO,GAOS,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,&
     & WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
  IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D8 = 8.0E0_realk, D4 = 4E0_realk
REAL(REALK) :: VXC(4,NBLEN),VX(5),GRDA,GRD,DMAX
INTEGER  :: IPNT,I,J,JBL,IBL,K
LOGICAL,EXTERNAL :: DFT_ISGGA
LOGICAL :: DOGGA
INTEGER     :: KVALS(3,3)
REAL(REALK),pointer :: GAORED(:,:,:),GDRED(:,:,:)
REAL(REALK),pointer :: DRED(:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,NRED,orb2atom(nbast)
INTEGER     :: atom(NACTBAST),iatom,IX,K1,K2,K3,KA,ik,jk
REAL(REALK) :: FRC,GA,GA2,GFS,XCFUNINPUT(4,1),XCFUNOUTPUT(5,1)

orb2atom = DFTDATA%orb2atom
KVALS(1:3,1) = (/1, 2, 3/)
KVALS(1:3,2) = (/2, 4, 5/)
KVALS(1:3,3) = (/3, 5, 6/)
DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
      !get the functional derivatives 
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         CALL dft_funcderiv1(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
         IF(DFTDATA%LB94)THEN
            CALL LB94correction(rho(IPNT,1),GRD,DFTDATA%HFexchangeFac,&
                 & WGHT(IPNT),VX(1))            
         ELSEIF(DFTDATA%CS00)THEN
            CALL CS00correction(rho(IPNT,1),GRD,DFTDATA%HFexchangeFac,&
                 & WGHT(IPNT),VX(1),DFTDATA%CS00SHIFT,DFTDATA%CS00eHOMO,DFTDATA%CS00ZND1,DFTDATA%CS00ZND2)
         ENDIF
         VXC(1,IPNT) = D2*VX(1) 
         IF(GRD.GT. 1E-40_realk) THEN
            GRDA = D05*GRD
            VXC(2,IPNT) = (VX(2)/GRDA + VX(3))*GRAD(1,IPNT,1)
            VXC(3,IPNT) = (VX(2)/GRDA + VX(3))*GRAD(2,IPNT,1)
            VXC(4,IPNT) = (VX(2)/GRDA + VX(3))*GRAD(3,IPNT,1)
         ELSE
            VXC(2,IPNT) = 0E0_realk
            VXC(3,IPNT) = 0E0_realk
            VXC(4,IPNT) = 0E0_realk
         ENDIF
#ifdef VAR_XCFUN
      ELSE
         XCFUNINPUT(1,1) = RHO(IPNT,1)
         XCFUNINPUT(2,1) = GRAD(1,IPNT,1)
         XCFUNINPUT(3,1) = GRAD(2,IPNT,1)
         XCFUNINPUT(4,1) = GRAD(3,IPNT,1)
         ! Input:
         !rho   = XCFUNINPUT(1,1)
         !grad_x = XCFUNINPUT(2,1)
         !grad_y = XCFUNINPUT(3,1)
         !grad_z = XCFUNINPUT(4,1)
         call xcfun_gga_components_xc_single_eval(XCFUNINPUT,5,XCFUNOUTPUT,1)
         ! Output
         ! Order 0
         ! out(1,1) Exc
         ! Order 1
         ! out(2,1) d^1 Exc / d rho
         ! out(3,1) d^1 Exc / d grad_x
         ! out(4,1) d^1 Exc / d grad_y
         ! out(5,1) d^1 Exc / d grad_z
         IF(DFTDATA%LB94)THEN
            call lsquit('error lb94 not implemented for xcfun',-1)
         ELSEIF(DFTDATA%CS00)THEN
            call lsquit('error cs00 not implemented for xcfun',-1)
         ENDIF
         VXC(1,IPNT) = D2*XCFUNOUTPUT(2,1)*WGHT(IPNT)
         VXC(2,IPNT) = D2*XCFUNOUTPUT(3,1)*WGHT(IPNT)
         VXC(3,IPNT) = D2*XCFUNOUTPUT(4,1)*WGHT(IPNT)
         VXC(4,IPNT) = D2*XCFUNOUTPUT(5,1)*WGHT(IPNT)
      ENDIF
#endif
   ELSE
      VXC(1,IPNT) = 0E0_realk
      VXC(2,IPNT) = 0E0_realk
      VXC(3,IPNT) = 0E0_realk
      VXC(4,IPNT) = 0E0_realk
   END IF
END DO

! Set up maximum density-matrix elements
DMAX = 0.0E0_realk
DO JBL=1, NBLOCKS
   DO IBL=1, NBLOCKS
      DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)        !J is active index
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)  !I is active index
            DMAX = MAX(DMAX,ABS(DMAT(I,J,1)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
! Count reduced number of AO's
NRED = 0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
         NRED = NRED + 1
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT. 0) THEN
   call mem_dft_alloc(DRED,NRED,NRED)
   call mem_dft_alloc(GAORED,NBLEN,NRED,NTYPSO)
   call mem_dft_alloc(GDRED,NBLEN,NRED,NTYPSO)
   ! Set up reduced Gaussian AO's
   IRED = 0
   DO IBL=1, NBLOCKS
      DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
         IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
            IRED = IRED + 1
            INXRED(IRED) = I
            DO J=1,NTYPSO
               DO K = 1, NBLEN
                  GAORED(K,IRED,J)  = GAOS(K,I,J)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ! Set up reduced density-matrix
   DO JRED=1,NRED            !Jred is reduced index
      J = INXRED(JRED)       !J is active index
      DO IRED=1,NRED         !Ired is reduced index
         I = INXRED(IRED)    !I is active index
         DRED(IRED,JRED) = DMAT(I,J,1)
      ENDDO
   ENDDO
   ! Set up reduced coordinate index in gradient
   DO IRED=1,NRED
      I = INXACT(INXRED(IRED)) !I is orbital index
      atom(IRED) = orb2atom(I)
   ENDDO
   
   ! Density-matrix contraction  
   ! \chi_{\mu}_{A,B} D_{\mu \nu}
   ! \frac{\partial \chi_{\mu}}{\frac \partial x} D_{\mu \nu}
   ! \frac{\partial \chi_{\mu}}{\frac \partial y} D_{\mu \nu}
   ! \frac{\partial \chi_{\mu}}{\frac \partial z} D_{\mu \nu}
   DO J=1,4
      IF(XCintNoOMP)THEN
         CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(1,1,J),&
           &        NBLEN,DRED,NRED,0.0E0_realk,GDRED(1,1,J),NBLEN )
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,J),&
              &        NBLEN,DRED,NRED,0.0E0_realk,GDRED(:,:,J),NBLEN )
      ENDIF
   ENDDO
   DO IRED=1,NRED
      iatom = atom(IRED)
      KA = INXRED(IRED) !KA is active index
      DO IX=1,3
         K1 = KVALS(1,IX) + 4
         K2 = KVALS(2,IX) + 4
         K3 = KVALS(3,IX) + 4
         FRC = 0E0_realk
         ! Assuming E_{\xc}=\int f[\rho ,\nabla \rho] d\textbf{r}
         ! \frac{ \partial E^{xc}[\rho]}{\partial R} =  
         ! \int \frac{\partial f }{\partial \rho} \frac{\partial \rho(\textbf{r})}{\partial R} d\textbf{r} + 
         ! \int \frac{\partial f }{\partial |\nabla \rho_{\alpha}|} \frac{\nabla \rho(\textbf{r})}{|\nabla \rho_{\alpha}|}  
         ! \frac{\partial \nabla \rho(\textbf{r})}{\partial R} d\textbf{r}
         ! VXC(1,I) = \frac{\partial f }{\partial \rho}
         ! VXC(2,I) = \int \frac{\partial f }{\partial |\nabla \rho_{\alpha}|}
         DO I = 1, NBLEN
            !\frac{\partial \chi_{\mu}}{\partial R_{\gamma}}
!            GA  = GAOS(I,KA,IX+1)
            !\nabla \rho \nabla \chi_{\mu} D_{\mu \nu}
!            GFS = VXC(2,I)*GDRED(I,IRED,2)+VXC(3,I)*GDRED(I,IRED,3)+VXC(4,I)*GDRED(I,IRED,4)
            !\nabla \rho \frac{\partial \nabla \chi_{\mu}}{\partial R_{\gamma}}
!            GA2 = VXC(2,I)*GAOS(I,KA,K1)  +VXC(3,I)*GAOS(I,KA,K2)  +VXC(4,I)*GAOS(I,KA,K3)
            !\int VXC(1)\frac{\partial \chi_{\mu}}{\partial R_{\gamma}}\chi_{\nu}D_{\mu \nu}
            ! + VXC(2) \nabla \rho \frac{\partial \nabla \chi_{\mu}}{\partial R_{\gamma}} \chi_{\nu}D_{\mu \nu}
            ! + VXC(2) \nabla \rho \nabla \chi_{\mu} \frac{\partial \chi_{\mu}}{\partial R_{\gamma}} D_{\mu \nu}
!            FRC = FRC + VXC(1,I)*GDRED(I,IRED,1)*GA + GDRED(I,IRED,1)*GA2 + GFS*GA

            FRC = FRC + VXC(1,I)*GDRED(I,IRED,1)*GAOS(I,KA,IX+1) &
            &+ VXC(2,I)*(GDRED(I,IRED,1)*GAOS(I,KA,K1) + GDRED(I,IRED,2)*GAOS(I,KA,IX+1)) &
            &+ VXC(3,I)*(GDRED(I,IRED,1)*GAOS(I,KA,K2) + GDRED(I,IRED,3)*GAOS(I,KA,IX+1)) &
            &+ VXC(4,I)*(GDRED(I,IRED,1)*GAOS(I,KA,K3) + GDRED(I,IRED,4)*GAOS(I,KA,IX+1))

         END DO
         DFTDATA%GRAD(IX,iatom) = DFTDATA%GRAD(IX,iatom) - FRC
      ENDDO ! IX
   ENDDO ! IRED
   call mem_dft_dealloc(DRED)
   call mem_dft_dealloc(GAORED)
   call mem_dft_dealloc(GDRED)
ENDIF !NRED GT 0

END SUBROUTINE II_GEODERIV_MOLGRAD_WORKER_GGA

!> \brief Main LDA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_linrsplda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                     RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                     GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK),pointer :: EXPVAL(:,:),VXC(:,:)
!EXPVAL(NBLEN,DFTDATA%NBMAT),VXC(NBLEN,DFTDATA%NBMAT),VX(9),DFTENE
REAL(REALK) :: VX(14),DFTENE
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred,W1,W2,W3,W4,W5,W6,W7,W8
LOGICAL     :: DOCALC
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk
Real(realk), parameter :: D4 = 4.0E0_realk
REAL(REALK) :: fRR
REAL(REALK) :: XCFUNINPUT(1,1),XCFUNOUTPUT(4,1)

W1 = 1
W2 = NBLEN*Nactbast                        !W1 - 1 + NBLEN*Nactbast    -> GAORED 
W3 = NBLEN*Nactbast+1                      !W2 + 1
W4 = (NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast -> EXCRED
W5 = (NBLEN+Nactbast)*Nactbast + 1         !W4 + 1 
W6 = (NBLEN+Nactbast + 1) * Nactbast       !W5 - 1 + Nactbast          -> GAOGMX
W7 = (NBLEN+Nactbast + 1) * Nactbast + 1   !W6 + 1
W8 = (2*NBLEN+Nactbast + 1) * Nactbast     !W7 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
IF(W8.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in ii_dft_linrsplda',lupri)
#endif

call mem_dft_alloc(EXPVAL,NBLEN,DFTDATA%NBMAT)
call mem_dft_alloc(VXC,NBLEN,DFTDATA%NBMAT)
NBMAT = DFTDATA%NBMAT
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
   !get expectation value of BMAT = \sum_{\mu \nu} \chi_{\mu} \chi_{\nu} BMAT_{\mu \nu}
 call II_get_expval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
        &Nactbast,NBAST,GAO,EXPVAL,SHAREDDFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
    !get the functional derivatives 
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         CALL dft_funcderiv2(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
         !fRR = 0.5*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
         fRR = VX(6) + VX(7)     
         DO IBMAT = 1,NBMAT
            !THE factor 0.5 IN fRR cancels a factor 2.0 in VXC so  
            VXC(IPNT,IBMAT) = D2*fRR*EXPVAL(IPNT,IBMAT) ! D4*fRR*EXPVAL(IPNT,IBMAT)
            !WARNING the factor 4 is due to 
            ! G_{\rho \sigma} &=& \frac{\delta (F^{\alpha} + F^{\beta})}{\delta D^{\alpha}_{\nu \mu}}\kappa^{\alpha}_{\nu \mu}\\
            ! &+& \frac{\delta (F^{\alpha} + F^{\beta})}{\delta D^{\beta}_{\nu \mu}}\kappa^{\beta}_{\nu \mu}\\
            ! G_{\rho \sigma} &=& 4 \frac{\delta (F^{\alpha}{\delta D^{\alpha}_{\nu \mu}} \kappa^{\alpha}_{\nu \mu}\\
            ! G_{\rho \sigma} &=& 4 \int \frac{\partial^{2} f }{\partial \rho^{2}_{\alpha}} \Omega_{\rho \sigma}
            ! \Omega_{\mu \nu} \kappa^{\alpha}_{\mu \nu} d\textbf{r}
         ENDDO
#ifdef VAR_XCFUN
      ELSE
         XCFUNINPUT(1,1) = RHO(IPNT,1)
         call xcfun2_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
         !1 = E
         !2 = fR
         !3 = fRR
         DO IBMAT = 1,NBMAT
            VXC(IPNT,IBMAT) =D4*XCFUNOUTPUT(3,1)*WGHT(IPNT)*EXPVAL(IPNT,IBMAT)
         ENDDO
      ENDIF
#endif
   ELSE
         VXC(IPNT,:) = 0.0E0_realk
   ENDIF
  END DO
  DO IBMAT = 1,NBMAT
   CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
        & VXC(:,IBMAT),GAO(:,:,1),SHAREDDFTDATA%FKSM(:,:,IBMAT),DFTHRI,&
        & WORK(W1:W2),WORK(W3:W4),WORK(W5:W6),WORK(W7:W8))
  ENDDO
 ENDIF
ENDIF
call mem_dft_dealloc(EXPVAL)
call mem_dft_dealloc(VXC)

END SUBROUTINE II_DFT_LINRSPLDA

!> \brief main unrestricted LDA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_LINRSPLDAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,&
     & NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     & GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
CALL LSQUIT('II_DFT_LINRSPLDAUNRES not implemented yet',lupri)

END SUBROUTINE II_DFT_LINRSPLDAUNRES

!> \brief main GGA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_LINRSPGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                      RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                      GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK),pointer :: EXPVAL(:,:),VXC(:,:,:),EXPGRAD(:,:,:)
REAL(REALK) :: VX(14),DFTENE,GRD,GRDA,MIXEDGRDA
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred,W1,W2,W3,W4,W5,W6,W7,W8
LOGICAL     :: DOCALC
Real(realk), parameter :: D4 = 4.0E0_realk, DUMMY = 0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D2 = 2.0E0_realk, D8 = 8.0E0_realk,D025 = 0.25E0_realk
Real(realk), parameter :: D16 = 16.0E0_realk,D32 = 32.0E0_realk
REAL(REALK) :: fR,fZ,fRR,fRZ,fZZ,fRG,fZG,fGG,fG,A,B
REAL(REALK) :: XCFUNINPUT(4,1),XCFUNOUTPUT(15,1)
NBMAT = DFTDATA%NBMAT
call mem_dft_alloc(EXPVAL,NBLEN,NBMAT)
call mem_dft_alloc(EXPGRAD,3,NBLEN,NBMAT)
call mem_dft_alloc(VXC,4,NBLEN,NBMAT)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
  !CALC
  !EXPVAL : the expectation value of BMAT = \sum_{\mu \nu} \chi_{\mu} \chi_{\nu} BMAT_{\mu \nu}
  !EXPGRAD: the gradient components of BMAT = \sum_{\mu \nu} \nabla (\chi_{\mu} \chi_{\nu}) BMAT_{\mu \nu}
  W1 = 1
  W2 = NBLEN*Nactbast             !W1 - 1 + NBLEN*Nactbast -> GAORED1 
  W3 = NBLEN*Nactbast+1
  W4 = NBLEN*Nactbast*4           !W1 - 1 + NBLEN*Nactbast*3 -> GAORED2 
  W5 = (4*NBLEN)*Nactbast + 1     !W4 + 1
  W6 = (5*NBLEN)*Nactbast + 1     !W5 - 1 + NBLEN*Nactbast -> TMP
#ifdef VAR_LSDEBUGINT
  IF(W6.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_LINRSPGGA',lupri)
#endif
  call II_get_expval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
       & Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,SHAREDDFTDATA%BMAT,&
       & DFTDATA%nBMAT,DFTHRI,NRED,WORK(W1:W2),WORK(W3:W4),GAOGMX,GAOMAX,&
       & WORK(W5:W6),MaxNactBast)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD.LT. 1E-40_realk) GRD = 1E-40_realk
   GRDA = D05*GRD
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
    !get the functional derivatives 
#ifdef VAR_XCFUN
    IF(.NOT.USEXCFUN)THEN
#endif
       CALL dft_funcderiv2(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
       fR  = D05*(VX(1) + VX(2))   !0.5*(drvs.df1000 + drvs.df0100);
       fZ  = VX(3)                    !drvs.df0010;
       fRR = D05*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
       fRZ = D05*(VX(8) + VX(9))   !0.5*(drvs.df1010 + drvs.df1001);
       fZZ = D05*(VX(11) + VX(12)) !0.5*(drvs.df0020 + drvs.df0011);
       fRG = D05*VX(10)             !0.5*drvs.df10001;   
       fZG = D05*VX(13)             !0.5*drvs.df00101; 
       fGG = D025*VX(14)            !0.25*drvs.df00002; 
       fG  = D05*VX(5)               !0.5*drvs.df00001;  
       DO IBMAT = 1,NBMAT
          MIXEDGRDA = (EXPGRAD(1,IPNT,IBMAT)*GRAD(1,IPNT,1)&
               &+EXPGRAD(2,IPNT,IBMAT)*GRAD(2,IPNT,1)&
               &+EXPGRAD(3,IPNT,IBMAT)*GRAD(3,IPNT,1))
          !the LDA part
          VXC(1,IPNT,IBMAT) =D4*fRR*EXPVAL(IPNT,IBMAT)+D4*(fRZ/GRD+fRG)*MIXEDGRDA
          !the non LDA parts
          A = D8*((fRZ/GRD + fRG)*EXPVAL(IPNT,IBMAT)&
               & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDA)
          B= D8*(fZ/GRD + fG)
          VXC(2,IPNT,IBMAT) = A*GRAD(1,IPNT,1)+B*EXPGRAD(1,IPNT,IBMAT)
          VXC(3,IPNT,IBMAT) = A*GRAD(2,IPNT,1)+B*EXPGRAD(2,IPNT,IBMAT)
          VXC(4,IPNT,IBMAT) = A*GRAD(3,IPNT,1)+B*EXPGRAD(3,IPNT,IBMAT)
       ENDDO
#ifdef VAR_XCFUN
    ELSE
       XCFUNINPUT(1,1) = RHO(IPNT,1)
       XCFUNINPUT(2,1) = GRAD(1,IPNT,1)
       XCFUNINPUT(3,1) = GRAD(2,IPNT,1)
       XCFUNINPUT(4,1) = GRAD(3,IPNT,1)
       ! Input:
       !rho   = XCFUNINPUT(1,1)
       !grad_x = XCFUNINPUT(2,1)
       !grad_y = XCFUNINPUT(3,1)
       !grad_z = XCFUNINPUT(4,1)
       call xcfun_gga_components_xc_single_eval(XCFUNINPUT,15,XCFUNOUTPUT,2)
       ! Output
       ! Order 0
       ! out(1,1) Exc
       ! Order 1
       ! out(2,1) d^1 Exc / d n
       ! out(3,1) d^1 Exc / d nx
       ! out(4,1) d^1 Exc / d ny
       ! out(5,1) d^1 Exc / d nz
       ! Order 2
       ! out(6,1) d^2 Exc / d n n
       ! out(7,1) d^2 Exc / d n nx
       ! out(8,1) d^2 Exc / d n ny
       ! out(9,1) d^2 Exc / d n nz
       ! out(10,1) d^2 Exc / d nx nx
       ! out(11,1) d^2 Exc / d nx ny
       ! out(12,1) d^2 Exc / d nx nz
       ! out(13,1) d^2 Exc / d ny ny
       ! out(14,1) d^2 Exc / d ny nz
       ! out(15,1) d^2 Exc / d nz nz
       DO IBMAT = 1,NBMAT
          !the \Omega_{\mu \nu} part
          VXC(1,IPNT,IBMAT) = &
               &   D4*XCFUNOUTPUT(6,1)*WGHT(IPNT)*EXPVAL(IPNT,IBMAT) &
               & + D4*XCFUNOUTPUT(7,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,IBMAT)&
               & + D4*XCFUNOUTPUT(8,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,IBMAT)&
               & + D4*XCFUNOUTPUT(9,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,IBMAT)
          !the \frac{\partial \Omega_{\mu \nu}}{\partial x} part
          VXC(2,IPNT,IBMAT) = &
               &   D8*XCFUNOUTPUT(7,1)*WGHT(IPNT)*EXPVAL(IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(10,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(11,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(12,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,IBMAT)
          !the \frac{\partial \Omega_{\mu \nu}}{\partial y} part             
          VXC(3,IPNT,IBMAT) = &
               &   D8*XCFUNOUTPUT(8,1)*WGHT(IPNT)*EXPVAL(IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(11,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(13,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(14,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,IBMAT)
          !the \frac{\partial \Omega_{\mu \nu}}{\partial z} part =
          !     (d^2 Exc/d rho d gz)kappa + (d^2 Exc/d gx d gz)dkappa/dx
          !+ (d^2 Exc/d gy d gz)dkappa/dy + (d^2 Exc/d gz d gz)dkappa/dz 
          VXC(4,IPNT,IBMAT) = &
               &   D8*XCFUNOUTPUT(9,1)*WGHT(IPNT)*EXPVAL(IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(12,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(14,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,IBMAT)&
               & + D8*XCFUNOUTPUT(15,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,IBMAT)
       ENDDO

    ENDIF
#endif
   ELSE
    VXC(:,IPNT,:) = 0.0E0_realk
   ENDIF
  END DO
  W1 = 1
  W2 = NBLEN*Nactbast*4                        !W1 - 1 + NBLEN*Nactbast*4 -> GAORED 
  W3 = 4*NBLEN*Nactbast+1                      !W2 + 1
  W4 = (4*NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast  -> EXCRED
  W5 = (4*NBLEN+Nactbast)*Nactbast + 1         !W4 + 1
  W6 = (5*NBLEN+Nactbast)*Nactbast             !W5 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
  IF(W6.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_KSMLDAUNRES',lupri)
#endif
  CALL II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NTYPSO,NBMAT,&
       & VXC,GAO,SHAREDDFTDATA%FKSM,DFTHRI,WORK(W1:W2),WORK(W3:W4),GAOGMX,GAOMAX,WORK(W5:W6),MaxNactBast)
 ENDIF
ENDIF
call mem_dft_dealloc(EXPVAL)
call mem_dft_dealloc(EXPGRAD)
call mem_dft_dealloc(VXC)

END SUBROUTINE II_DFT_LINRSPGGA

!> \brief main unrestricted GGA linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_LINRSPGGAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                           RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                           GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
CALL LSQUIT('II_DFT_LINRSPGGAUNRES not implemented yet',lupri)

END SUBROUTINE II_DFT_LINRSPGGAUNRES

!> \brief LDA Exchange-correlation contribution to quadratic response 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_quadrsplda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                       RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                       GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK),pointer :: EXPVAL(:,:),VXC(:)
!EXPVAL(NBLEN,DFTDATA%NBMAT),VXC(NBLEN,DFTDATA%NBMAT),VX(9),DFTENE
REAL(REALK) :: VX(27),DFTENE
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred
INTEGER     :: W1,W2,W3,W4,W5,W6,W7,W8
LOGICAL     :: DOCALC
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk,D3 = 3.0E0_realk, D8=8.0E0_realk
REAL(REALK) :: fRRR
REAL(REALK) :: XCFUNINPUT(1,1),XCFUNOUTPUT(4,1)

W1 = 1
W2 = NBLEN*Nactbast                        !W1 - 1 + NBLEN*Nactbast    -> GAORED 
W3 = NBLEN*Nactbast+1                      !W2 + 1
W4 = (NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast -> EXCRED
W5 = (NBLEN+Nactbast)*Nactbast + 1         !W4 + 1 
W6 = (NBLEN+Nactbast + 1) * Nactbast       !W5 - 1 + Nactbast          -> GAOGMX
W7 = (NBLEN+Nactbast + 1) * Nactbast + 1   !W6 + 1
W8 = (2*NBLEN+Nactbast + 1) * Nactbast     !W7 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
IF(W8.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_QUADRSPLDA',lupri)
#endif

call mem_dft_alloc(EXPVAL,NBLEN,DFTDATA%NBMAT)
call mem_dft_alloc(VXC,NBLEN)
NBMAT = DFTDATA%NBMAT
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
!get expectation value of BMAT=\sum_{\mu \nu}\chi_{\mu}\chi_{\nu}BMAT_{\mu \nu}
 call II_get_expval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
        &Nactbast,NBAST,GAO,EXPVAL,SHAREDDFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
    !get the functional derivatives 
#ifdef VAR_XCFUN
    IF(.NOT.USEXCFUN)THEN
#endif
       CALL dft_funcderiv3(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
       !    fR  = VX(1)              !drvs.df1000;
       !radovan: factor 0.5 "missing" here compared to "true" derivative wrt RR i'm 
       !sure it's compensated elsewhere but can be confusing when comparing with 
       !other codes
       !     fRR = (VX(6) + VX(7))   !(drvs.df2000 + drvs.df1100);
       !radovan: factor 0.25 "missing" here compared to "true" derivative wrt RRR
       !        (i'm sure it's compensated elsewhere)
       !         but can be confusing when comparing with other codes
       fRRR = (VX(15)+D3*VX(16))     !(drvs.df3000 + 3*drvs.df2100);
       !    VXCB(IPNT) = fRR*EXPVALB(IPNT) 
       !    VXCC(IPNT) = fRR*EXPVALC(IPNT) 
       VXC(IPNT) = D2*fRRR*EXPVAL(IPNT,1)*EXPVAL(IPNT,2)  
#ifdef VAR_XCFUN
    ELSE
        XCFUNINPUT(1,1) = RHO(IPNT,1)
        call xcfun3_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
        !1 = E
        !2 = fR
        !3 = fRR
        !4 = fRRR
        DO IBMAT = 1,NBMAT
           fRRR = XCFUNOUTPUT(4,1)*WGHT(IPNT)
           VXC(IPNT) =D8*fRRR*EXPVAL(IPNT,1)*EXPVAL(IPNT,2)
        ENDDO
    ENDIF
#endif
   ELSE
    VXC(IPNT) = 0.0E0_realk
   ENDIF
  END DO
   CALL II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
        & VXC(:),GAO(:,:,1),SHAREDDFTDATA%FKSM(:,:,1),DFTHRI,WORK(W1:W2),&
        & WORK(W3:W4),WORK(W5:W6),WORK(W7:W8))
 ENDIF
ENDIF
call mem_dft_dealloc(EXPVAL)
call mem_dft_dealloc(VXC)

END SUBROUTINE II_DFT_QUADRSPLDA

!> \brief GGA Exchange-correlation contribution to quadratic response 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_QUADRSPGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                       RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                       GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK),pointer :: EXPVAL(:,:),VXC(:,:),EXPGRAD(:,:,:)
REAL(REALK) :: VX(27),DFTENE,GRD,GRDA,MIXEDGRDA
INTEGER     :: I,J,NBMAT,IPNT,IBMAT,nred,W1,W2,W3,W4,W5,W6,W7,W8,X,Y,Z,XYZ
LOGICAL     :: DOCALC
Real(realk), parameter :: D4 = 4.0E0_realk, DUMMY = 0E0_realk,D05 = 0.5E0_realk,D3 = 3E0_realk
Real(realk), parameter :: D2 = 2.0E0_realk, D8 = 8.0E0_realk,D025 = 0.25E0_realk
REAL(REALK) :: GRD2,GRDA2,GRDA3
REAL(REALK) :: fRZ,fRG,fZZ,fRRR,fRRZ,fRRG,fRRGX,fRZZ,fZZZ,gradY,gradZ,gradYZ,A,B,C
REAL(REALK) :: XCFUNINPUT(4,1),XCFUNOUTPUT(35,1)
REAL(realk) :: FDERIV(4,4,4),E1(4),E2(4),TMP
NBMAT = DFTDATA%NBMAT
IF(NBMAT.NE. 2)call lsquit('QRSP XC error',lupri)
call mem_dft_alloc(EXPVAL,NBLEN,NBMAT)
call mem_dft_alloc(EXPGRAD,3,NBLEN,NBMAT)
call mem_dft_alloc(VXC,4,NBLEN)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 W1 = 1
 W2 = NBLEN*Nactbast             !W1 - 1 + NBLEN*Nactbast -> GAORED1 
 W3 = NBLEN*Nactbast+1
 W4 = NBLEN*Nactbast*4           !W1 - 1 + NBLEN*Nactbast*3 -> GAORED2 
 W5 = 4*NBLEN*Nactbast + 1   !W6 + 1
 W6 = 5*NBLEN*Nactbast + 1   !W7 - 1 + NBLEN*Nactbast -> TMP
#ifdef VAR_LSDEBUGINT
 IF(W6.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_QUADRSPGGA',lupri)
#endif
 call II_get_expval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
      & Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,SHAREDDFTDATA%BMAT,&
      & DFTDATA%nBMAT,DFTHRI,NRED,WORK(W1:W2),WORK(W3:W4),&
      & GAOGMX,GAOMAX,WORK(W5:W6),MaxNactBast)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD.LT. 1E-40_realk) GRD = 1E-40_realk
   GRDA = D05*GRD
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
#ifdef VAR_XCFUN
    IF(.NOT.USEXCFUN)THEN
#endif
       !get the functional derivatives 
       CALL dft_funcderiv3(RHO(IPNT,1),GRD,WGHT(IPNT),VX)

       GRD2 = GRD*GRD
       GRDA2 = GRDA*GRDA
       GRDA3 = GRDA2*GRDA
       gradY = D05*(EXPGRAD(1,IPNT,1)*GRAD(1,IPNT,1) &
            &+EXPGRAD(2,IPNT,1)*GRAD(2,IPNT,1) &
            &+EXPGRAD(3,IPNT,1)*GRAD(3,IPNT,1))
       gradZ = D05*(EXPGRAD(1,IPNT,2)*GRAD(1,IPNT,1) &
            &+EXPGRAD(2,IPNT,2)*GRAD(2,IPNT,1) &
            &+EXPGRAD(3,IPNT,2)*GRAD(3,IPNT,1))
       gradYZ = (EXPGRAD(1,IPNT,2)*EXPGRAD(1,IPNT,1) &
            &+EXPGRAD(2,IPNT,2)*EXPGRAD(2,IPNT,1) &
            &+EXPGRAD(3,IPNT,2)*EXPGRAD(3,IPNT,1))

       fRZ = (VX(8)+VX(9))/GRD              !(drvs.df1010 + drvs.df1001)/(2*grada)
       fRG = D2*VX(10)                      !2*drvs.df10001   
       fZZ = (VX(11)+VX(12))/GRD2-VX(3)/(GRD2*GRDA) !(drvs.df0020 + drvs.df0011)/(4*grada2)-drvs.df0010/(4*grada3)
       fRRR = VX(15)+D3*VX(16)              !(drvs.df3000 + 3*drvs.df2100) 
       fRRZ = (VX(17)+VX(18)+D2*VX(20))/GRD !(drvs.df2010+drvs.df2001+2*drvs.df1110)/(2*grada)
       fRRG = VX(19)+VX(21)                 !drvs.df20001+drvs.df11001
       fRRGX = D2*(VX(19)+VX(21))           !2*(drvs.df20001+drvs.df11001)
       fRZZ = (VX(22)+VX(24)+D2*VX(23))/GRD2 - (VX(8)+VX(9))/(GRD2*GRDA)  !(drvs.df1020+drvs.df0120+2*drvs.df1011)/(4*grada2)-(drvs.df1010+drvs.df1001)/(4*grada3)
       fZZZ = ((VX(25)+D3*VX(26))/(GRDA3)-D3*(VX(11)+VX(12))/(GRDA2*GRDA2)+D3*VX(3)/(GRDA3*GRDA2))/D8
       !((drvs.df0030 + 3*drvs.df0021)/grada3& 
       !         &-3*(drvs.df0020 + drvs.df0011)/(grada2*grada2)&
       !         &+3*drvs.df0010/(grada3*grada2))/8.0
       VXC(1,IPNT) = D2*fRRR*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) &!OK
            &+D4*(fRRZ*EXPVAL(IPNT,1)*gradZ+fRRZ*EXPVAL(IPNT,2)*gradY) & !OK
            &+D8*gradZ*gradY*fRZZ &! OK
            &+D4*gradYZ*fRZ &! OK  
            &+D4*fRRG*EXPVAL(IPNT,1)*gradZ &! OK
            &+D4*fRRG*EXPVAL(IPNT,2)*gradY+D2*fRG*gradYZ !OK
       
       A = D8*fZZZ*gradY*gradZ &
            & + D4*(fRZZ*EXPVAL(IPNT,1)*gradZ + fRZZ*EXPVAL(IPNT,2)*gradY) &
            & + D2*fRRZ*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) &
            & + fRRGX*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) + D4*fZZ*gradYZ
       B = D8*fZZ*gradY + D4*fRZ*EXPVAL(IPNT,1) + D2*fRG*EXPVAL(IPNT,1)
       C = D8*fZZ*gradZ + D4*fRZ*EXPVAL(IPNT,2) + D2*fRG*EXPVAL(IPNT,2)

       VXC(2,IPNT) = D2*A*GRAD(1,IPNT,1) + D2*B*EXPGRAD(1,IPNT,2) + D2*C*EXPGRAD(1,IPNT,1)
       VXC(3,IPNT) = D2*A*GRAD(2,IPNT,1) + D2*B*EXPGRAD(2,IPNT,2) + D2*C*EXPGRAD(2,IPNT,1)
       VXC(4,IPNT) = D2*A*GRAD(3,IPNT,1) + D2*B*EXPGRAD(3,IPNT,2) + D2*C*EXPGRAD(3,IPNT,1)
       
#ifdef VAR_XCFUN
    ELSE
       DO IBMAT = 1,NBMAT
       ! Input:
       !rho    = XCFUNINPUT(1,1)
       !grad_x = XCFUNINPUT(2,1)
       !grad_y = XCFUNINPUT(3,1)
       !grad_z = XCFUNINPUT(4,1)
       XCFUNINPUT(1,1) = RHO(IPNT,1)
       XCFUNINPUT(2,1) = GRAD(1,IPNT,1)
       XCFUNINPUT(3,1) = GRAD(2,IPNT,1)
       XCFUNINPUT(4,1) = GRAD(3,IPNT,1)
       call xcfun_gga_components_xc_single_eval(XCFUNINPUT,35,XCFUNOUTPUT,3)
       ! Output
       ! Order 0
       ! out(1,1) Exc
       ! Order 1
       ! out(2,1) d^1 Exc / d n
       ! out(3,1) d^1 Exc / d nx
       ! out(4,1) d^1 Exc / d ny
       ! out(5,1) d^1 Exc / d nz
       ! Order 2
       ! out(6,1) d^2 Exc / d n n
       ! out(7,1) d^2 Exc / d n nx
       ! out(8,1) d^2 Exc / d n ny
       ! out(9,1) d^2 Exc / d n nz
       ! out(10,1) d^2 Exc / d nx nx
       ! out(11,1) d^2 Exc / d nx ny
       ! out(12,1) d^2 Exc / d nx nz
       ! out(13,1) d^2 Exc / d ny ny
       ! out(14,1) d^2 Exc / d ny nz
       ! out(15,1) d^2 Exc / d nz nz
       ! Order 3
       ! out(16,1) d^3 Exc / d n n n
       ! out(17,1) d^3 Exc / d n n nx
       ! out(18,1) d^3 Exc / d n n ny
       ! out(19,1) d^3 Exc / d n n nz
       ! out(20,1) d^3 Exc / d n nx nx
       ! out(21,1) d^3 Exc / d n nx ny
       ! out(22,1) d^3 Exc / d n nx nz
       ! out(23,1) d^3 Exc / d n ny ny
       ! out(24,1) d^3 Exc / d n ny nz
       ! out(25,1) d^3 Exc / d n nz nz
       ! out(26,1) d^3 Exc / d nx nx nx
       ! out(27,1) d^3 Exc / d nx nx ny
       ! out(28,1) d^3 Exc / d nx nx nz
       ! out(29,1) d^3 Exc / d nx ny ny
       ! out(30,1) d^3 Exc / d nx ny nz
       ! out(31,1) d^3 Exc / d nx nz nz
       ! out(32,1) d^3 Exc / d ny ny ny
       ! out(33,1) d^3 Exc / d ny ny nz
       ! out(34,1) d^3 Exc / d ny nz nz
       ! out(35,1) d^3 Exc / d nz nz nz
       
       ! Constructing intermediates for the VXC automated loop
       XYZ = 15
       DO X=1,4
         DO Y=X,4
           DO Z=Y,4
             XYZ = XYZ + 1
             FDERIV(X,Y,Z) = XCFUNOUTPUT(XYZ,1)
             FDERIV(X,Z,Y) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Y,X,Z) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Y,Z,X) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Z,X,Y) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Z,Y,X) = XCFUNOUTPUT(XYZ,1)
           ENDDO
         ENDDO
       ENDDO
       E1(1) = EXPVAL(IPNT,1)
       E1(2) = EXPGRAD(1,IPNT,1)
       E1(3) = EXPGRAD(2,IPNT,1)
       E1(4) = EXPGRAD(3,IPNT,1)
       E2(1) = EXPVAL(IPNT,2)
       E2(2) = EXPGRAD(1,IPNT,2)
       E2(3) = EXPGRAD(2,IPNT,2)
       E2(4) = EXPGRAD(3,IPNT,2)

       ! VXC loop
       DO X=1,4
         VXC(X,IPNT) = 0E0_realk
         DO Y=1,4
           TMP = D8*E1(Y)*WGHT(IPNT)
           IF (X.GE.2) TMP = D2*TMP
           DO Z=1,4
             VXC(X,IPNT) = VXC(X,IPNT) + FDERIV(Z,Y,X)*E2(Z)*TMP
           ENDDO
         ENDDO
       ENDDO

       ENDDO !IPNT
     ENDIF !XCFUN
#endif
   ELSE
    VXC(:,IPNT) = 0.0E0_realk
   ENDIF
  END DO
  W1 = 1
  W2 = NBLEN*Nactbast*4                        !W1 - 1 + NBLEN*Nactbast*4 -> GAORED 
  W3 = 4*NBLEN*Nactbast+1                      !W2 + 1
  W4 = (4*NBLEN+Nactbast)*Nactbast             !W3 - 1 + Nactbast*Nactbast  -> EXCRED
  W5 = (4*NBLEN+Nactbast)*Nactbast + 1         !W4 + 1
  W6 = (5*NBLEN+Nactbast)*Nactbast             !W5 - 1 + NBLEN*Nactbast    -> TMP
#ifdef VAR_LSDEBUGINT
  IF(W6.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_QUADRSPGGA',lupri)
#endif
  CALL II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NTYPSO,1,&
       & VXC,GAO,SHAREDDFTDATA%FKSM,DFTHRI,WORK(W1:W2),WORK(W3:W4),GAOGMX,GAOMAX,WORK(W5:W6),MaxNactBast)
 ENDIF
ENDIF
call mem_dft_dealloc(EXPVAL)
call mem_dft_dealloc(EXPGRAD)
call mem_dft_dealloc(VXC)

END SUBROUTINE II_DFT_QUADRSPGGA

!> \brief magnetic derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_magderiv_kohnshamLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D4 = 4.0E0_realk,DUMMY = 0E0_realk
REAL(REALK) :: VXC(NBLEN),VX(5),DFTENE
INTEGER     :: I,J,IPNT
REAL(REALK) :: XCFUNINPUT(1,1),XCFUNOUTPUT(2,1)

! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IPNT = 1, NBLEN
!   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000   = \frac{\partial  f}{\partial \rho_{\alpha}}
#ifdef VAR_XCFUN
   IF(.NOT.USEXCFUN)THEN
#endif
      CALL dft_funcderiv1(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
      !WARNING the factor 2 is due to F=F_alpha + F_beta = 2 F_alpha 
      VXC(IPNT) = VX(1)*D4 
#ifdef VAR_XCFUN
   ELSE
      call lsquit('xcfun version of II_DFT_MAGDERIV_KOHNSHAMLDA not implemented',-1)
      XCFUNINPUT(1,1) = RHO(IPNT,1)
      call xcfun_lda_xc_eval(XCFUNINPUT,XCFUNOUTPUT,1)
      !WARNING the factor 2 is due to F=F_alpha + F_beta = 2 F_alpha 
      VXC(IPNT) = D4*XCFUNOUTPUT(2,1)*WGHT(IPNT)
      call lsquit('xcfun inconsistency',-1)
   ENDIF
#endif
!   ELSE
!      VXC(IPNT) = 0.0E0_realk
!   ENDIF
END DO
!ntypso should be 4
CALL II_DFT_distmagderiv_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     &VXC,GAO(:,:,:),SHAREDDFTDATA%FKSM(:,:,1:3),COORD,DFTHRI,NTYPSO)

END SUBROUTINE II_DFT_MAGDERIV_KOHNSHAMLDA

!> \brief  magnetic derivative Kohn-sham matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_magderiv_kohnshamGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
implicit none
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk,D3 = 3.0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D4 = 4.0E0_realk
REAL(REALK) :: VX(5),DFTENE,GRD,GRDA,A,XCFUNINPUT(4,1),XCFUNOUTPUT(5,1)
REAL(REALK),pointer :: VXC(:,:)
INTEGER     :: I,J,W1,W2,W3,W4,W5,W6,IDMAT,IPNT
IDMAT = 1
IF(NDMAT.GT.1)CALL LSQUIT('II_DFT_MAGDERIV_KOHNSHAMGGA ndmat',-1)
call mem_dft_alloc(VXC,4,NBLEN)
DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         CALL dft_funcderiv1(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
         IF(DFTDATA%LB94)THEN
            CALL LB94correction(rho(IPNT,IDMAT),GRD,DFTDATA%HFexchangeFac,&
                 & WGHT(IPNT),VX(1))            
         ELSEIF(DFTDATA%CS00)THEN
            CALL CS00correction(rho(IPNT,IDMAT),GRD,DFTDATA%HFexchangeFac,&
                 & WGHT(IPNT),VX(1),DFTDATA%CS00SHIFT,DFTDATA%CS00eHOMO,DFTDATA%CS00ZND1,DFTDATA%CS00ZND2)
         ENDIF
         VXC(1,IPNT) = D2*VX(1) 
         IF(GRD.GT. 1E-40_realk) THEN
            GRDA = D05*GRD
            A = D2*(VX(2)/GRDA + VX(3))
            VXC(2,IPNT) = A*GRAD(1,IPNT,1)
            VXC(3,IPNT) = A*GRAD(2,IPNT,1)
            VXC(4,IPNT) = A*GRAD(3,IPNT,1)
         ELSE
            VXC(2,IPNT) = 0E0_realk
            VXC(3,IPNT) = 0E0_realk
            VXC(4,IPNT) = 0E0_realk
         ENDIF
#ifdef VAR_XCFUN
      ELSE
         XCFUNINPUT(1,1) = RHO(IPNT,IDMAT)
         XCFUNINPUT(2,1) = GRAD(1,IPNT,1)
         XCFUNINPUT(3,1) = GRAD(2,IPNT,1)
         XCFUNINPUT(4,1) = GRAD(3,IPNT,1)
         call xcfun_gga_components_xc_single_eval(XCFUNINPUT,5,XCFUNOUTPUT,1)
         ! Output
         ! Order 0
         ! out(1,1) Exc
         ! Order 1
         ! out(2,1) d^1 Exc / d rho
         ! out(3,1) d^1 Exc / d grad_x
         ! out(4,1) d^1 Exc / d grad_y
         ! out(5,1) d^1 Exc / d grad_z

         IF(DFTDATA%LB94)THEN
            call lsquit('error lb94 xcfun',-1)
         ELSEIF(DFTDATA%CS00)THEN
            call lsquit('error cs00 xcfun',-1)
         ENDIF

         VXC(1,IPNT) = D2*XCFUNOUTPUT(2,1)*WGHT(IPNT)
         VXC(2,IPNT) = D4*XCFUNOUTPUT(3,1)*WGHT(IPNT)
         VXC(3,IPNT) = D4*XCFUNOUTPUT(4,1)*WGHT(IPNT)
         VXC(4,IPNT) = D4*XCFUNOUTPUT(5,1)*WGHT(IPNT)
      ENDIF
#endif
   ELSE
      VXC(:,IPNT) = 0E0_realk
   END IF
END DO


W1 = 1
W2 = NBLEN*Nactbast*NTYPSO        !W1 - 1 + NBLEN*Nactbast*3 -> GAORED 
W3 = W2+1
W4 = W3 - 1 + NBLEN*Nactbast*3    !W3 - 1 + NBLEN*Nactbast*3 -> TMP
W5 = W4+1                         
W6 = W5 - 1 + Nactbast            !W5 - 1 + Nactbast -> GAOGMX
CALL II_DFT_distmagderiv_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     &VXC,GAO,SHAREDDFTDATA%FKSM,COORD,DFTHRI,NTYPSO,WORK(W1:W2),WORK(W3:W4),WORK(W5:W6))
call mem_dft_dealloc(VXC)
END SUBROUTINE II_DFT_MAGDERIV_KOHNSHAMGGA

!> \brief magnetic derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_MAGDERIV_LINRSPLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D4 = 4.0E0_realk,DUMMY = 0E0_realk
REAL(REALK) :: fRR,VX(14)
INTEGER     :: I,J,NBMAT,ibmat,N,NRED,IPNT
LOGICAL     :: DOCALC,dosympart
REAL(REALK),parameter :: D2=2E0_realk
REAL(REALK),pointer :: EXPVAL(:,:),VXC(:,:)
REAL(REALK),pointer :: EXPGRAD(:,:,:),VXC2(:,:,:)

NBMAT = DFTDATA%NBMAT
dosympart = DFTDATA%dosympart
call mem_dft_alloc(EXPVAL,NBLEN,NBMAT)
call mem_dft_alloc(EXPGRAD,NBLEN,3,NBMAT)!magnetic gradients of EXPVAL=BMAT*chi(r)*chi(r)
call mem_dft_alloc(VXC,NBLEN,NBMAT)
call mem_dft_alloc(VXC2,NBLEN,NBMAT,3)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
!get expectation value of BMAT
 call II_get_magderiv_expval_lda(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,&
        & INXACT,Nactbast,NBAST,GAO,COORD,EXPVAL,EXPGRAD,SHAREDDFTDATA%BMAT,nBMAT,&
        & NRED,DFTHRI,dosympart)
 IF(NRED.GT. 0)THEN
  VXC2 = 0.0E0_realk
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
    !get the functional derivatives 
#ifdef VAR_XCFUN
     IF(.NOT.USEXCFUN)THEN
#endif
        CALL dft_funcderiv2(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
        !fRR = 0.5*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
        fRR = VX(6) + VX(7)     
        !THE factor 0.5 IN fRR cancels a factor 2.0 in VXC so  
        DO IBMAT = 1,NBMAT
           VXC(IPNT,IBMAT) = D2*fRR*EXPVAL(IPNT,IBMAT)!D4*fRR*EXPVAL(IPNT,IBMAT)
        ENDDO
        DO N=1,3
           DO IBMAT = 1,NBMAT
              VXC2(IPNT,IBMAT,N)=D2*fRR*EXPGRAD(IPNT,N,IBMAT)
           ENDDO
        ENDDO
#ifdef VAR_XCFUN
     ELSE
        call lsquit('xcfun version of II_DFT_MAGDERIV_LINRSPLDA not implemented',-1)
     ENDIF
#endif
   ELSE
    VXC(IPNT,:) = 0.0E0_realk
   ENDIF
  END DO
  CALL II_DFT_DISTMAGDERIV_linrsp_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
    &Nactbast,NBAST,VXC,VXC2,NBMAT,GAO,SHAREDDFTDATA%FKSM,SHAREDDFTDATA%FKSMS,COORD,DFTHRI,NTYPSO,dosympart)
 ENDIF
ENDIF
call mem_dft_dealloc(EXPVAL)
call mem_dft_dealloc(EXPGRAD)
call mem_dft_dealloc(VXC)
call mem_dft_dealloc(VXC2)

END SUBROUTINE II_DFT_MAGDERIV_LINRSPLDA

!> \brief magnetic derivative linrsp matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_MAGDERIV_LINRSPGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
INTEGER     :: I,J,NBMAT,ibmat,N,NRED,IPNT,IJ
REAL(REALK) :: VX(14),GRD,GRDA,MIXEDGRDA,MIXEDGRDAMAG(3)
LOGICAL     :: DOCALC,dosympart
REAL(REALK),pointer :: VXC1(:,:),VXC1MAG(:,:,:)
REAL(REALK),pointer :: EXPGRAD(:,:,:,:),VXC2(:,:,:),VXC2MAG(:,:,:,:)
Real(realk), parameter :: D0 = 0.0E0_realk
Real(realk), parameter :: D4 = 4.0E0_realk, DUMMY = 0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D2 = 2.0E0_realk, D8 = 8.0E0_realk,D025 = 0.25E0_realk
REAL(REALK) :: fR,fZ,fRR,fRZ,fZZ,fRG,fZG,fGG,fG,A,B,AMAG(3)
REAL(REALK) :: XCFUNINPUT(4,1),XCFUNOUTPUT(15,1)
REAL(REALK) :: fxc(4,4),vxc(4,4),FAC
NBMAT = DFTDATA%NBMAT
dosympart = DFTDATA%dosympart
call mem_dft_alloc(EXPGRAD,NBLEN,4,4,NBMAT)!mixed geo,magn gradients of EXPVAL
call mem_dft_alloc(VXC1,NBLEN,NBMAT)
call mem_dft_alloc(VXC1MAG,NBLEN,NBMAT,3)
call mem_dft_alloc(VXC2,NBLEN,NBMAT,3)
call mem_dft_alloc(VXC2MAG,NBLEN,NBMAT,3,3)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
!get expectation value of BMAT
 call II_get_magderiv_expval_gga(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,&
        & INXACT,Nactbast,NBAST,GAO,COORD,EXPGRAD,SHAREDDFTDATA%BMAT,nBMAT,&
        & NRED,DFTHRI,dosympart)
 IF(NRED.GT. 0)THEN
  VXC2 = 0.0E0_realk
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
    IF(GRD.LT. 1E-40_realk) GRD = 1E-40_realk
    GRDA = D05*GRD
    !get the functional derivatives 
#ifdef VAR_XCFUN
    IF(.NOT.USEXCFUN)THEN
#endif
       CALL dft_funcderiv2(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
       fR  = D05*(VX(1) + VX(2))   !0.5*(drvs.df1000 + drvs.df0100);
       fZ  = VX(3)                    !drvs.df0010;
       fRR = D05*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
       fRZ = D05*(VX(8) + VX(9))   !0.5*(drvs.df1010 + drvs.df1001);
       fZZ = D05*(VX(11) + VX(12)) !0.5*(drvs.df0020 + drvs.df0011);
       fRG = D05*VX(10)             !0.5*drvs.df10001;   
       fZG = D05*VX(13)             !0.5*drvs.df00101; 
       fGG = D025*VX(14)            !0.25*drvs.df00002; 
       fG  = D05*VX(5)               !0.5*drvs.df00001;  
       DO IBMAT = 1,NBMAT
          MIXEDGRDA = (EXPGRAD(IPNT,2,1,IBMAT)*GRAD(1,IPNT,1)&
               &+EXPGRAD(IPNT,3,1,IBMAT)*GRAD(2,IPNT,1)&
               &+EXPGRAD(IPNT,4,1,IBMAT)*GRAD(3,IPNT,1))
          VXC1(IPNT,IBMAT) =D4*fRR*EXPGRAD(IPNT,1,1,IBMAT)+D4*(fRZ/GRD+fRG)*MIXEDGRDA
          !magnetic differentiation of VXC1
          DO N=1,3
             MIXEDGRDAMAG(N) = (EXPGRAD(IPNT,2,N+1,IBMAT)*GRAD(1,IPNT,1)&
                  &+EXPGRAD(IPNT,3,N+1,IBMAT)*GRAD(2,IPNT,1)&
                  &+EXPGRAD(IPNT,4,N+1,IBMAT)*GRAD(3,IPNT,1))
          ENDDO
          DO N=1,3
             VXC1MAG(IPNT,IBMAT,N)=D4*(fRZ/GRD+fRG)*MIXEDGRDAMAG(N)+D4*fRR*EXPGRAD(IPNT,1,1+N,IBMAT)
          ENDDO
          A = D8*((fRZ/GRD + fRG)*EXPGRAD(IPNT,1,1,IBMAT)&
               & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDA)
          B= D8*(fZ/GRD + fG)
          VXC2(IPNT,IBMAT,1) = A*GRAD(1,IPNT,1)+B*EXPGRAD(IPNT,2,1,IBMAT)
          VXC2(IPNT,IBMAT,2) = A*GRAD(2,IPNT,1)+B*EXPGRAD(IPNT,3,1,IBMAT)
          VXC2(IPNT,IBMAT,3) = A*GRAD(3,IPNT,1)+B*EXPGRAD(IPNT,4,1,IBMAT)
          !magnetic differentiation of VXC2
          DO N=1,3
             AMAG(N) = D8*((fRZ/GRD + fRG)*EXPGRAD(IPNT,1,N+1,IBMAT)&
                  & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDAMAG(N))
          ENDDO
          DO N=1,3
             VXC2MAG(IPNT,IBMAT,1,N)=AMAG(N)*GRAD(1,IPNT,1)+B*EXPGRAD(IPNT,2,1+N,IBMAT)
             VXC2MAG(IPNT,IBMAT,2,N)=AMAG(N)*GRAD(2,IPNT,1)+B*EXPGRAD(IPNT,3,1+N,IBMAT)
             VXC2MAG(IPNT,IBMAT,3,N)=AMAG(N)*GRAD(3,IPNT,1)+B*EXPGRAD(IPNT,4,1+N,IBMAT)
          ENDDO
       ENDDO
#ifdef VAR_XCFUN
    ELSE
       XCFUNINPUT(1,1) = RHO(IPNT,1)
       XCFUNINPUT(2,1) = GRAD(1,IPNT,1)
       XCFUNINPUT(3,1) = GRAD(2,IPNT,1)
       XCFUNINPUT(4,1) = GRAD(3,IPNT,1)
       ! Input:
       !rho   = XCFUNINPUT(1,1)
       !grad_x = XCFUNINPUT(2,1)
       !grad_y = XCFUNINPUT(3,1)
       !grad_z = XCFUNINPUT(4,1)
       call xcfun_gga_components_xc_single_eval(XCFUNINPUT,15,XCFUNOUTPUT,2)
       ! Output
       ! Order 0
       ! out(1,1) Exc
       ! Order 1
       ! out(2,1) d^1 Exc / d n
       ! out(3,1) d^1 Exc / d nx
       ! out(4,1) d^1 Exc / d ny
       ! out(5,1) d^1 Exc / d nz
       ! Order 2
       ! out(6,1) d^2 Exc / d n n
       ! out(7,1) d^2 Exc / d n nx
       ! out(8,1) d^2 Exc / d n ny
       ! out(9,1) d^2 Exc / d n nz
       ! out(10,1) d^2 Exc / d nx nx
       ! out(11,1) d^2 Exc / d nx ny
       ! out(12,1) d^2 Exc / d nx nz
       ! out(13,1) d^2 Exc / d ny ny
       ! out(14,1) d^2 Exc / d ny nz
       ! out(15,1) d^2 Exc / d nz nz
       ij=5
       DO i=1,4
         DO j=i,4
           ij=ij+1
           fxc(i,j) = XCFUNOUTPUT(ij,1)*WGHT(IPNT)
           fxc(j,i) = fxc(i,j)
           vxc(i,j) = D0
           vxc(j,i) = D0
         ENDDO
       ENDDO
       DO IBMAT = 1,NBMAT
          VXC1(IPNT,IBMAT) = D0
          DO i=1,4
            FAC = D4
            IF (i.GE.2) FAC = D8
            DO j=1,4
              DO N=1,4
                vxc(i,j) = vxc(i,j) + FAC*fxc(i,N)*EXPGRAD(IPNT,N,j,IBMAT)
              ENDDO
            ENDDO
          ENDDO
          VXC1(IPNT,IBMAT)        = vxc(1,1)
          VXC2(IPNT,IBMAT,1:3)    = vxc(2:4,1)
          VXC1MAG(IPNT,IBMAT,1:3) = vxc(1,2:4)
          VXC2MAG(IPNT,IBMAT,1:3,1:3) = vxc(2:4,2:4)
       ENDDO
    ENDIF
#endif
   ELSE
    VXC1(IPNT,:) = 0.0E0_realk
    VXC1MAG(IPNT,:,:) = 0.0E0_realk
    VXC2(IPNT,:,:) = 0.0E0_realk
    VXC2MAG(IPNT,:,:,:) = 0.0E0_realk
   ENDIF
  END DO
  CALL II_DFT_DISTMAGDERIV_linrsp_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
    &Nactbast,NBAST,VXC1,VXC1MAG,VXC2,VXC2MAG,NBMAT,GAO,&
    &SHAREDDFTDATA%FKSM,SHAREDDFTDATA%FKSMS,COORD,DFTHRI,NTYPSO,dosympart)
 ENDIF
ENDIF
call mem_dft_dealloc(EXPGRAD)
call mem_dft_dealloc(VXC1)
call mem_dft_dealloc(VXC2)
call mem_dft_dealloc(VXC1MAG)
call mem_dft_dealloc(VXC2MAG)

END SUBROUTINE II_DFT_MAGDERIV_LINRSPGGA

!> \brief geometrical derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_kohnshamLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK) :: VX(14),DFTENE,EXPVAL(NBLEN,1)
REAL(REALK) :: fR(NBLEN),fRR(NBLEN)
REAL(REALK),parameter :: D2=2E0_realk,D4=4E0_realk,DUMMY = 0E0_realk 
INTEGER     :: I,J,nred,nbmat,ipnt
logical :: DOCALC
REAL(REALK) :: XCFUNINPUT(1,1),XCFUNOUTPUT(4,1)
nbmat=1
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 call II_get_expval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
!      &Nactbast,NBAST,GAO,EXPVAL,SHAREDDFTDATA%BMAT,DFTDATA%nBMAT,DFTHRI,NRED)
      &Nactbast,NBAST,GAO,EXPVAL,SHAREDDFTDATA%BMAT,1,DFTHRI,NRED)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
#ifdef VAR_XCFUN
    IF(.NOT.USEXCFUN)THEN
#endif
      CALL dft_funcderiv2(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
      fR(IPNT) = VX(1)*D4 
      fRR(IPNT) = D2*(VX(6) + VX(7))*EXPVAL(IPNT,1) 
#ifdef VAR_XCFUN
    ELSE
       XCFUNINPUT(1,1) = RHO(IPNT,1)
       call xcfun2_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
       !1 = E
       !2 = fR
       !3 = fRR
       fR(IPNT)  = D4*XCFUNOUTPUT(2,1)*WGHT(IPNT)
       fRR(IPNT) = D4*XCFUNOUTPUT(3,1)*WGHT(IPNT)*EXPVAL(IPNT,1)
    ENDIF
#endif
   ELSE
      fR(IPNT) = 0.0E0_realk
      fRR(IPNT) = 0.0E0_realk
   ENDIF
  END DO
   CALL II_DFT_distgeoderiv_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
        &NBAST,fR,fRR,GAO(:,:,:),DFTDATA%GRAD,DFTDATA%orb2atom,&
        &SHAREDDFTDATA%BMAT,nbmat,DMAT,ndmat,DFTDATA%natoms,COORD,DFTHRI,NTYPSO)
 ENDIF
ENDIF
END SUBROUTINE II_DFT_GEODERIV_KOHNSHAMLDA

!> \brief geometrical derivative Kohn-sham matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_kohnshamGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK) :: VX(14),DFTENE,EXPVAL(NBLEN,1)
REAL(REALK) :: EXPGRAD(3,NBLEN,1),GRD,GRDA,MIXEDGRDA
REAL(REALK) :: fR,fRR,fZ,fRZ,fZZ,fRG,fZG,fGG,fG,A,B
REAL(REALK) :: VXC1(NBLEN),VXC2(3,NBLEN),VXC3(NBLEN),VXC4(3,NBLEN)
REAL(REALK),parameter :: D2=2E0_realk, D4=4E0_realk, DUMMY = 0E0_realk, D05=0.5E0_realk, D025=0.25E0_realk 
REAL(REALK),parameter :: D8=8E0_realk
REAL(REALK) :: XCFUNINPUT(4,1),XCFUNOUTPUT(15,1)
INTEGER     :: I,J,nred,nbmat,IPNT,W1,W2,W3,W4,W7,W8,W5,W6,W9,W10
logical :: DOCALC
nbmat=1
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 W1 = 1
 W2 = NBLEN*Nactbast             !W1 - 1 + NBLEN*Nactbast -> GAORED1 
 W3 = NBLEN*Nactbast+1
 W4 = NBLEN*Nactbast*4           !W1 - 1 + NBLEN*Nactbast*3 -> GAORED2 
 W5 = 4*NBLEN*Nactbast + 1       !W4 + 1
 W6 = 5*NBLEN*Nactbast + 1       !W5 - 1 + NBLEN*Nactbast -> TMP
#ifdef VAR_LSDEBUGINT
 IF(W6.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_geoderiv_kohnshamGGA',lupri)
#endif
 call II_get_expval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
      & Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,SHAREDDFTDATA%BMAT,&
      & DFTDATA%nBMAT,DFTHRI,NRED,WORK(W1:W2),WORK(W3:W4),&
      & GAOGMX,GAOMAX,WORK(W5:W6),MaxNactBast)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
     GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
          &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
     GRDA = D05*GRD
     IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
       CALL dft_funcderiv2(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
       fR  = D05*(VX(1) + VX(2))   !0.5*(drvs.df1000 + drvs.df0100);
       fZ  = VX(3)                    !drvs.df0010;
       fRR = D05*(VX(6) + VX(7))   !0.5*(drvs.df2000 + drvs.df1100);
       fRZ = D05*(VX(8) + VX(9))   !0.5*(drvs.df1010 + drvs.df1001);
       fZZ = D05*(VX(11) + VX(12)) !0.5*(drvs.df0020 + drvs.df0011);
       fRG = D05*VX(10)             !0.5*drvs.df10001;   
       fZG = D05*VX(13)             !0.5*drvs.df00101; 
       fGG = D025*VX(14)            !0.25*drvs.df00002; 
       fG  = D05*VX(5)               !0.5*drvs.df00001;  
       MIXEDGRDA = (EXPGRAD(1,IPNT,1)*GRAD(1,IPNT,1)&
            &+EXPGRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
            &+EXPGRAD(3,IPNT,1)*GRAD(3,IPNT,1))

       VXC1(IPNT) = D4*VX(1)
       A = D2*(VX(3)/GRDA + VX(5))
       VXC2(1,IPNT) = A*GRAD(1,IPNT,1)
       VXC2(2,IPNT) = A*GRAD(2,IPNT,1)
       VXC2(3,IPNT) = A*GRAD(3,IPNT,1)
       !the LDA part
       VXC3(IPNT) =D4*fRR*EXPVAL(IPNT,1)+D4*(fRZ/GRD+fRG)*MIXEDGRDA
       !the non LDA parts
       A = D4*((fRZ/GRD + fRG)*EXPVAL(IPNT,1)&
            & + (((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)*MIXEDGRDA)
       B= D4*(fZ/GRD + fG)
       VXC4(1,IPNT) = A*GRAD(1,IPNT,1)+B*EXPGRAD(1,IPNT,1)
       VXC4(2,IPNT) = A*GRAD(2,IPNT,1)+B*EXPGRAD(2,IPNT,1)
       VXC4(3,IPNT) = A*GRAD(3,IPNT,1)+B*EXPGRAD(3,IPNT,1)
#ifdef VAR_XCFUN
      ELSE
       XCFUNINPUT(1,1) = RHO(IPNT,1)
       XCFUNINPUT(2,1) = GRAD(1,IPNT,1)
       XCFUNINPUT(3,1) = GRAD(2,IPNT,1)
       XCFUNINPUT(4,1) = GRAD(3,IPNT,1)
       ! Input:
       !rho    = XCFUNINPUT(1,1)
       !grad_x = XCFUNINPUT(2,1)
       !grad_y = XCFUNINPUT(3,1)
       !grad_z = XCFUNINPUT(4,1)
       call xcfun_gga_components_xc_single_eval(XCFUNINPUT,15,XCFUNOUTPUT,2)
       ! Output
       ! Order 0
       ! out(1,1) Exc
       ! Order 1
       ! out(2,1) d^1 Exc / d n
       ! out(3,1) d^1 Exc / d nx
       ! out(4,1) d^1 Exc / d ny
       ! out(5,1) d^1 Exc / d nz
       ! Order 2
       ! out(6,1) d^2 Exc / d n n
       ! out(7,1) d^2 Exc / d n nx
       ! out(8,1) d^2 Exc / d n ny
       ! out(9,1) d^2 Exc / d n nz
       ! out(10,1) d^2 Exc / d nx nx
       ! out(11,1) d^2 Exc / d nx ny
       ! out(12,1) d^2 Exc / d nx nz
       ! out(13,1) d^2 Exc / d ny ny
       ! out(14,1) d^2 Exc / d ny nz
       ! out(15,1) d^2 Exc / d nz nz
       !the \Omega_{\mu \nu} part
       VXC1(IPNT)   = D4*XCFUNOUTPUT(2,1)*WGHT(IPNT)
       VXC2(1,IPNT) = D4*XCFUNOUTPUT(3,1)*WGHT(IPNT)
       VXC2(2,IPNT) = D4*XCFUNOUTPUT(4,1)*WGHT(IPNT)
       VXC2(3,IPNT) = D4*XCFUNOUTPUT(5,1)*WGHT(IPNT)

       VXC3(IPNT) = &
               &   D4*XCFUNOUTPUT(6,1)*WGHT(IPNT)*EXPVAL(IPNT,1) &
               & + D4*XCFUNOUTPUT(7,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,1)&
               & + D4*XCFUNOUTPUT(8,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,1)&
               & + D4*XCFUNOUTPUT(9,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,1)
       VXC4(1,IPNT) = &
            &   D4*XCFUNOUTPUT(7,1)*WGHT(IPNT)*EXPVAL(IPNT,1)&
            & + D4*XCFUNOUTPUT(10,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,1)&
            & + D4*XCFUNOUTPUT(11,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,1)&
            & + D4*XCFUNOUTPUT(12,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,1)
       VXC4(2,IPNT) = &
            &   D4*XCFUNOUTPUT(8,1)*WGHT(IPNT)*EXPVAL(IPNT,1)&
            & + D4*XCFUNOUTPUT(11,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,1)&
            & + D4*XCFUNOUTPUT(13,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,1)&
            & + D4*XCFUNOUTPUT(14,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,1)
       VXC4(3,IPNT) = &
            &   D4*XCFUNOUTPUT(9,1)*WGHT(IPNT)*EXPVAL(IPNT,1)&
            & + D4*XCFUNOUTPUT(12,1)*WGHT(IPNT)*EXPGRAD(1,IPNT,1)&
            & + D4*XCFUNOUTPUT(14,1)*WGHT(IPNT)*EXPGRAD(2,IPNT,1)&
            & + D4*XCFUNOUTPUT(15,1)*WGHT(IPNT)*EXPGRAD(3,IPNT,1)
      ENDIF
#endif
    ELSE
       VXC1(IPNT) = 0.0E0_realk
       VXC2(:,IPNT) = 0.0E0_realk
       VXC3(IPNT) = 0.0E0_realk
       VXC4(:,IPNT) = 0.0E0_realk
    ENDIF
  END DO
  W1  = 1
  W2  = NBLEN*Nactbast*NTYPSO      !W1 - 1 + NBLEN*Nactbast*NTYPSO -> GAORED
  W3  = W2+1                       
  W4  = W2+NBLEN*NactBast          !TMP(NBLEN,NactBast)
  W5  = W4+1                       
  W6  = W4+NBLEN*NactBast          !TMPX(NBLEN,NactBast)
  W7  = W6+1                       
  W8  = W6+NBLEN*NactBast          !TMPY(NBLEN,NactBast)
  W9  = W8+1                       
  W10 = W8+NBLEN*NactBast          !TMPZ(NBLEN,NactBast)
#ifdef VAR_LSDEBUGINT
  IF(W10.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_geoderiv_kohnshamGGA',lupri)
#endif
  CALL II_DFT_distgeoderiv_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
      &NBAST,VXC1,VXC2,VXC3,VXC4,GAO(:,:,:),DFTDATA%GRAD,&
      &DFTDATA%orb2atom,SHAREDDFTDATA%BMAT,DMAT,DFTDATA%natoms,&
      &COORD,DFTHRI,NTYPSO,WORK(W1:W2),WORK(W3:W4),WORK(W5:W6),&
      &WORK(W7:W8),WORK(W9:W10))
 ENDIF
ENDIF
END SUBROUTINE II_DFT_GEODERIV_KOHNSHAMGGA

!> \brief LDA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_linrspLDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK) :: VX(27),DFTENE,EXPVAL(NBLEN,2)
REAL(REALK) :: fRRR(NBLEN),fRR(2,NBLEN),A
REAL(REALK),parameter :: D2=2E0_realk,D4=4E0_realk,D3=3E0_realk,DUMMY = 0E0_realk 
REAL(REALK),parameter :: D8=8E0_realk
INTEGER     :: I,J,nred,nbmat,IPNT
logical :: DOCALC
REAL(REALK) :: XCFUNINPUT(1,1),XCFUNOUTPUT(4,1)

nbmat=DFTDATA%nBMAT
IF(nbmat.NE. 2)call LSQUIT('II_DFT_geoderiv_linrspLDA requires 2 matrices',lupri)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 call II_get_expval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
      &Nactbast,NBAST,GAO,EXPVAL,SHAREDDFTDATA%BMAT,nbmat,DFTHRI,NRED)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
#ifdef VAR_XCFUN
     IF(.NOT.USEXCFUN)THEN
#endif
      CALL dft_funcderiv3(RHO(IPNT,1),DUMMY,WGHT(IPNT),VX)
      A = D4*(VX(6) + VX(7))
      fRR(1,IPNT) = A*EXPVAL(IPNT,1) 
      fRR(2,IPNT) = A*EXPVAL(IPNT,2) 
      A = D2*(VX(15) + D3*VX(16))
      fRRR(IPNT) = A*EXPVAL(IPNT,1)*EXPVAL(IPNT,2)
#ifdef VAR_XCFUN
     ELSE
       XCFUNINPUT(1,1) = RHO(IPNT,1)
       call xcfun3_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
       !1 = E
       !2 = fR
       !3 = fRR
       !4 = fRRR
       fRR(1,IPNT) =  D8*XCFUNOUTPUT(3,1)*WGHT(IPNT)*EXPVAL(IPNT,1)
       fRR(2,IPNT) =  D8*XCFUNOUTPUT(3,1)*WGHT(IPNT)*EXPVAL(IPNT,2)
       fRRR(IPNT)  =  D8*XCFUNOUTPUT(4,1)*WGHT(IPNT)*EXPVAL(IPNT,1)*EXPVAL(IPNT,2)
     ENDIF
#endif
   ELSE
      fRR(1,IPNT) = 0.0E0_realk
      fRR(2,IPNT) = 0.0E0_realk
      fRRR(IPNT) = 0.0E0_realk
   ENDIF
  END DO
   CALL II_DFT_distgeoderiv2_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
        &NBAST,fRR,fRRR,GAO(:,:,:),DFTDATA%GRAD,DFTDATA%orb2atom,&
        &SHAREDDFTDATA%BMAT,nbmat,DMAT,ndmat,DFTDATA%natoms,COORD,DFTHRI,NTYPSO)
 ENDIF
ENDIF
END SUBROUTINE II_DFT_GEODERIV_LINRSPLDA

!> \brief the GGA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_linrspGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     & Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,&
     & DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER,intent(in)     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK) :: VX(27),DFTENE,EXPVAL(NBLEN,2)
REAL(REALK) :: EXPGRAD(3,NBLEN,2)
REAL(REALK) :: VXC1(2,NBLEN),VXC2(3,NBLEN,2),VXC3(4,NBLEN)
REAL(REALK) :: ARHOGRAD,BRHOGRAD,ABRHOGRAD,GRDA,GRD2,GRD,GRDA2,GRDA3
REAL(REALK),parameter :: D2=2E0_realk,D4=4E0_realk,DUMMY = 0E0_realk,D05=0.5E0_realk,D025=0.25E0_realk 
REAL(REALK),parameter :: D8=8E0_realk,D3=3E0_realk,D16=16E0_realk
REAL(REALK) :: fR,fZ,fRZ,fZZ,fRG,fZG,fGG,fG,B,facW,factorRZ,A,fRR,fRRR,fRRZ,fRRG
REAL(REALK) :: fRRGX,fRZZ,fZZZ,gradA,gradB,gradAB,C
INTEGER     :: I,J,nred,nbmat,IPNT,W1,W2,W3,W4,W7,W8,W5,W6,X,Y,Z,XY,XYZ
logical :: DOCALC
REAL(REALK) :: XCFUNINPUT(4,1),XCFUNOUTPUT(35,1)
REAL(realk) :: FDERIV(4,4,4),E1(4),E2(4),TMP,F2(4,4),VXCTMP(4,2)
nbmat=DFTDATA%nBMAT
IF(nbmat.NE. 2)call LSQUIT('II_DFT_geoderiv_linrspGGA requires 2 matrices',lupri)
DOCALC = .FALSE.
DO IPNT = 1, NBLEN
   IF(RHO(IPNT,1) .GT. RHOTHR)THEN
      DOCALC = .TRUE.
      EXIT
   ENDIF
ENDDO
IF(DOCALC)THEN
 W1 = 1
 W2 = NBLEN*Nactbast             !W1 - 1 + NBLEN*Nactbast -> GAORED1 
 W3 = NBLEN*Nactbast+1
 W4 = NBLEN*Nactbast*4           !W1 - 1 + NBLEN*Nactbast*3 -> GAORED2 
 W5 = 4*NBLEN*Nactbast + 1       !W4 + 1
 W6 = 5*NBLEN*Nactbast           !W5 - 1 + NBLEN*Nactbast -> TMP
#ifdef VAR_LSDEBUGINT
 IF(W6.GT.WORKLENGTH) CALL LSQUIT('WORKLENGTH error in II_DFT_geoderiv_linrspGGA',lupri)
#endif
 call II_get_expval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
      & Nactbast,NBAST,GAO,EXPVAL,EXPGRAD,SHAREDDFTDATA%BMAT,&
      & DFTDATA%nBMAT,DFTHRI,NRED,WORK(W1:W2),WORK(W3:W4),&
      & GAOGMX,GAOMAX,WORK(W5:W6),MaxNactBast)
 IF(NRED.GT. 0)THEN
  DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,1)*GRAD(1,IPNT,1)+GRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
        &+GRAD(3,IPNT,1)*GRAD(3,IPNT,1))
   GRDA = D05*GRD
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,1).GT.RHOTHR) THEN
    ARHOGRAD = (EXPGRAD(1,IPNT,1)*GRAD(1,IPNT,1)&
         &+EXPGRAD(2,IPNT,1)*GRAD(2,IPNT,1)&
         &+EXPGRAD(3,IPNT,1)*GRAD(3,IPNT,1))
    BRHOGRAD = (EXPGRAD(1,IPNT,2)*GRAD(1,IPNT,1)&
         &+EXPGRAD(2,IPNT,2)*GRAD(2,IPNT,1)&
         &+EXPGRAD(3,IPNT,2)*GRAD(3,IPNT,1))
    ABRHOGRAD = (EXPGRAD(1,IPNT,1)*EXPGRAD(1,IPNT,2)&
         &+EXPGRAD(2,IPNT,1)*EXPGRAD(2,IPNT,2)&
         &+EXPGRAD(3,IPNT,1)*EXPGRAD(3,IPNT,2))     
    GRD2 = GRD*GRD
    GRDA2 = GRDA*GRDA
    GRDA3 = GRDA2*GRDA
#ifdef VAR_XCFUN
    IF(.NOT.USEXCFUN)THEN
#endif
       CALL dft_funcderiv3(RHO(IPNT,1),GRD,WGHT(IPNT),VX)
       !THE NON DIFFERENTIATED PART
       fR  = D05*(VX(1) + VX(2))    !0.5*(drvs.df1000 + drvs.df0100);
       fZ  = VX(3)                  !drvs.df0010;
       fRR = D05*(VX(6) + VX(7))    !0.5*(drvs.df2000 + drvs.df1100);
       fRZ = D05*(VX(8) + VX(9))    !0.5*(drvs.df1010 + drvs.df1001);
       fZZ = D05*(VX(11) + VX(12))  !0.5*(drvs.df0020 + drvs.df0011);
       fRG = D05*VX(10)             !0.5*drvs.df10001;   
       fZG = D05*VX(13)             !0.5*drvs.df00101; 
       fGG = D025*VX(14)            !0.25*drvs.df00002; 
       fG  = D05*VX(5)              !0.5*drvs.df00001;  
       VXC1(1,IPNT) = D8*fRR*EXPVAL(IPNT,1) + D8*(fRZ/GRD+fRG)*ARHOGRAD
       VXC1(2,IPNT) = D8*fRR*EXPVAL(IPNT,2) + D8*(fRZ/GRD+fRG)*BRHOGRAD
       B= D8*(fZ/GRD + fG)
       FacW = D8*(((-fZ/GRD+fZZ)/GRD + D2*fZG)/GRD + fGG)
       FactorRZ = D8*(fRZ/GRD + fRG)
       A = FactorRZ*EXPVAL(IPNT,1)+FacW*ARHOGRAD
       VXC2(1,IPNT,1) = B*EXPGRAD(1,IPNT,1)+A*GRAD(1,IPNT,1)
       VXC2(2,IPNT,1) = B*EXPGRAD(2,IPNT,1)+A*GRAD(2,IPNT,1)
       VXC2(3,IPNT,1) = B*EXPGRAD(3,IPNT,1)+A*GRAD(3,IPNT,1)
       A = FactorRZ*EXPVAL(IPNT,2)+FacW*BRHOGRAD
       VXC2(1,IPNT,2) = B*EXPGRAD(1,IPNT,2)+A*GRAD(1,IPNT,1)
       VXC2(2,IPNT,2) = B*EXPGRAD(2,IPNT,2)+A*GRAD(2,IPNT,1)
       VXC2(3,IPNT,2) = B*EXPGRAD(3,IPNT,2)+A*GRAD(3,IPNT,1)
       !THE DIFFERENTIATED PART
       fRRR = (VX(15) + D3*VX(16))
       fRZ = (VX(8)+VX(9))/GRD  !warning change definition of fRZ
       fRG = D4*fRG                 !D2*VX(10)                      
       fZZ = (VX(11)+VX(12))/GRD2-VX(3)/(GRD2*GRDA)
       fRRZ = (VX(17)+VX(18)+D2*VX(20))/GRD 
       fRRG = VX(19)+VX(21)                 
       fRRGX = D2*(VX(19)+VX(21))           
       fRZZ = (VX(22)+VX(24)+D2*VX(23))/GRD2 - (VX(8)+VX(9))/(GRD2*GRDA)  
       fZZZ = ((VX(25)+D3*VX(26))/(GRDA3)-D3*(VX(11)+VX(12))/(GRDA2*GRDA2)+D3*VX(3)/(GRDA3*GRDA2))/D8
       
       VXC3(1,IPNT) = D2*fRRR*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) & 
            &+D2*fRRZ*(EXPVAL(IPNT,1)*BRHOGRAD+EXPVAL(IPNT,2)*ARHOGRAD) & 
            &+D2*BRHOGRAD*ARHOGRAD*fRZZ& 
            &+D4*ABRHOGRAD*fRZ &
            &+D2*fRRG*(EXPVAL(IPNT,1)*BRHOGRAD+EXPVAL(IPNT,2)*ARHOGRAD)&
            &+D2*fRG*ABRHOGRAD 
       
       A = D2*fZZZ*ARHOGRAD*BRHOGRAD &
            & + D2*(fRZZ*EXPVAL(IPNT,1)*BRHOGRAD + fRZZ*EXPVAL(IPNT,2)*ARHOGRAD) &
            & + D2*fRRZ*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) &
            & + fRRGX*EXPVAL(IPNT,1)*EXPVAL(IPNT,2) + D4*fZZ*ABRHOGRAD
       B = D4*fZZ*ARHOGRAD + D4*fRZ*EXPVAL(IPNT,1) + D2*fRG*EXPVAL(IPNT,1)
       C = D4*fZZ*BRHOGRAD + D4*fRZ*EXPVAL(IPNT,2) + D2*fRG*EXPVAL(IPNT,2)
       
       VXC3(2,IPNT) = A*GRAD(1,IPNT,1) + B*EXPGRAD(1,IPNT,2) + C*EXPGRAD(1,IPNT,1)
       VXC3(3,IPNT) = A*GRAD(2,IPNT,1) + B*EXPGRAD(2,IPNT,2) + C*EXPGRAD(2,IPNT,1)
       VXC3(4,IPNT) = A*GRAD(3,IPNT,1) + B*EXPGRAD(3,IPNT,2) + C*EXPGRAD(3,IPNT,1)
#ifdef VAR_XCFUN
    ELSE
       ! Input:
       !rho    = XCFUNINPUT(1,1)
       !grad_x = XCFUNINPUT(2,1)
       !grad_y = XCFUNINPUT(3,1)
       !grad_z = XCFUNINPUT(4,1)
       XCFUNINPUT(1,1) = RHO(IPNT,1)
       XCFUNINPUT(2,1) = GRAD(1,IPNT,1)
       XCFUNINPUT(3,1) = GRAD(2,IPNT,1)
       XCFUNINPUT(4,1) = GRAD(3,IPNT,1)
       call xcfun_gga_components_xc_single_eval(XCFUNINPUT,35,XCFUNOUTPUT,3)
       ! Output
       ! Order 0
       ! out(1,1) Exc
       ! Order 1
       ! out(2,1) d^1 Exc / d n
       ! out(3,1) d^1 Exc / d nx
       ! out(4,1) d^1 Exc / d ny
       ! out(5,1) d^1 Exc / d nz
       ! Order 2
       ! out(6,1) d^2 Exc / d n n
       ! out(7,1) d^2 Exc / d n nx
       ! out(8,1) d^2 Exc / d n ny
       ! out(9,1) d^2 Exc / d n nz
       ! out(10,1) d^2 Exc / d nx nx
       ! out(11,1) d^2 Exc / d nx ny
       ! out(12,1) d^2 Exc / d nx nz
       ! out(13,1) d^2 Exc / d ny ny
       ! out(14,1) d^2 Exc / d ny nz
       ! out(15,1) d^2 Exc / d nz nz
       ! Order 3
       ! out(16,1) d^3 Exc / d n n n
       ! out(17,1) d^3 Exc / d n n nx
       ! out(18,1) d^3 Exc / d n n ny
       ! out(19,1) d^3 Exc / d n n nz
       ! out(20,1) d^3 Exc / d n nx nx
       ! out(21,1) d^3 Exc / d n nx ny
       ! out(22,1) d^3 Exc / d n nx nz
       ! out(23,1) d^3 Exc / d n ny ny
       ! out(24,1) d^3 Exc / d n ny nz
       ! out(25,1) d^3 Exc / d n nz nz
       ! out(26,1) d^3 Exc / d nx nx nx
       ! out(27,1) d^3 Exc / d nx nx ny
       ! out(28,1) d^3 Exc / d nx nx nz
       ! out(29,1) d^3 Exc / d nx ny ny
       ! out(30,1) d^3 Exc / d nx ny nz
       ! out(31,1) d^3 Exc / d nx nz nz
       ! out(32,1) d^3 Exc / d ny ny ny
       ! out(33,1) d^3 Exc / d ny ny nz
       ! out(34,1) d^3 Exc / d ny nz nz
       ! out(35,1) d^3 Exc / d nz nz nz
       
       ! Constructing intermediates for the VXC automated loop
       XYZ = 15
       XY  = 5
       DO X=1,4
         DO Y=X,4
           XY = XY + 1
           F2(X,Y) = XCFUNOUTPUT(XY,1)
           F2(Y,X) = XCFUNOUTPUT(XY,1)
           DO Z=Y,4
             XYZ = XYZ + 1
             FDERIV(X,Y,Z) = XCFUNOUTPUT(XYZ,1)
             FDERIV(X,Z,Y) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Y,X,Z) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Y,Z,X) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Z,X,Y) = XCFUNOUTPUT(XYZ,1)
             FDERIV(Z,Y,X) = XCFUNOUTPUT(XYZ,1)
           ENDDO
         ENDDO
       ENDDO
       E1(1) = EXPVAL(IPNT,1)
       E1(2) = EXPGRAD(1,IPNT,1)
       E1(3) = EXPGRAD(2,IPNT,1)
       E1(4) = EXPGRAD(3,IPNT,1)
       E2(1) = EXPVAL(IPNT,2)
       E2(2) = EXPGRAD(1,IPNT,2)
       E2(3) = EXPGRAD(2,IPNT,2)
       E2(4) = EXPGRAD(3,IPNT,2)

       ! VXC3 loop
       DO X=1,4
         VXC3(X,IPNT) = 0E0_realk
         DO Y=1,4
           TMP = D8*E1(Y)*WGHT(IPNT)
           DO Z=1,4
             VXC3(X,IPNT) = VXC3(X,IPNT) + FDERIV(Z,Y,X)*E2(Z)*TMP
           ENDDO
         ENDDO
       ENDDO

       ! VXC2 loop
       DO X=1,4
         VXCTMP(X,1) = 0E0_realk
         VXCTMP(X,2) = 0E0_realk
         DO Y=1,4
           TMP = -D8*WGHT(IPNT)*F2(Y,X)
           VXCTMP(X,1) = VXCTMP(X,1) - E1(Y)*TMP
           VXCTMP(X,2) = VXCTMP(X,2) - E2(Y)*TMP
         ENDDO
       ENDDO

       VXC1(1,IPNT)   = VXCTMP(1,1)
       VXC2(1,IPNT,1) = VXCTMP(2,1)
       VXC2(2,IPNT,1) = VXCTMP(3,1)
       VXC2(3,IPNT,1) = VXCTMP(4,1)
       VXC1(2,IPNT)   = VXCTMP(1,2)
       VXC2(1,IPNT,2) = VXCTMP(2,2)
       VXC2(2,IPNT,2) = VXCTMP(3,2)
       VXC2(3,IPNT,2) = VXCTMP(4,2)

     ENDIF !XCFUN
#endif
   ELSE
      VXC1(:,IPNT) = 0.0E0_realk
      VXC2(:,IPNT,:) = 0.0E0_realk
      VXC3(:,IPNT) = 0.0E0_realk
   ENDIF !RHO/GRAD THR
  END DO !IPNT
  CALL II_DFT_distgeoderiv2_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,&
       &NBAST,VXC1,VXC2,VXC3,GAO(:,:,:),DFTDATA%GRAD,DFTDATA%orb2atom,&
       &SHAREDDFTDATA%BMAT,nbmat,DMAT,ndmat,DFTDATA%natoms,COORD,DFTHRI,NTYPSO)
 ENDIF
ENDIF

END SUBROUTINE II_DFT_GEODERIV_LINRSPGGA

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Worker routine that for a batch of gridpoints build the LDA kohn-sham matrix
!>
SUBROUTINE II_DFT_KSMELDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                    RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                    GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
REAL(REALK) :: VX(5),DFTENE, Econt,XCFUNINPUT(1,1),XCFUNOUTPUT(2,1)
INTEGER     :: IPNT,I,J,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk
EXTERNAL DFTENE
! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IDMAT = 1, NDMAT
 Econt = 0.0E0_realk
 DO IPNT = 1, NBLEN
   IF(RHO(IPNT,IDMAT) .GT. RHOTHR)THEN
#ifdef VAR_XCFUN        
     IF(.NOT.USEXCFUN)THEN
#endif
        Econt = Econt + DFTENE(RHO(IPNT,IDMAT),DUMMY)*WGHT(IPNT)
#ifdef VAR_XCFUN        
     ELSE
        XCFUNINPUT(1,1) = RHO(IPNT,IDMAT)
        call xcfun_lda_xc_eval(XCFUNINPUT,XCFUNOUTPUT,1)
        DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + XCFUNOUTPUT(1,1)*WGHT(IPNT)
     ENDIF
#endif
   ENDIF
 END DO
 DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + Econt
ENDDO

END SUBROUTINE II_DFT_KSMELDA

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Worker routine that for a batch of gridpoints build the unrestricted LDA kohn-sham matrix
!>
SUBROUTINE II_DFT_KSMELDAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,&
     &                         GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,&
     &                         WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk
INTEGER     :: IPNT,I,J,IDMAT,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT1,IDMAT2
REAL(REALK) :: VX(5),DFTENEUNRES
EXTERNAL DFTENEUNRES
! LDA Exchange-correlation contribution to Kohn-Sham energy
DO IDMAT = 1, NDMAT/2
 IDMAT1 = 1 + (IDMAT-1)*2
 IDMAT2 = 2 + (IDMAT-1)*2
 DO IPNT = 1, NBLEN
  IF(RHO(IPNT,IDMAT1) .GT. RHOTHR .OR. RHO(IPNT,IDMAT2) .GT. RHOTHR)THEN
      !get the functional derivatives 
      !vx(1) = drvs.df1000   = \frac{\partial  f}{\partial \rho_{\alpha}} 
      !and same for beta
#ifdef VAR_XCFUN
     IF(.NOT.USEXCFUN)THEN
#endif
        !      CALL dft_funcderiv1unres(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),DUMMY,DUMMY,WGHT(IPNT),VX)
        DFTDATA%ENERGY(IDMAT1) = DFTDATA%ENERGY(IDMAT1) + DFTENEUNRES(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),DUMMY,DUMMY)*WGHT(IPNT)
#ifdef VAR_XCFUN
     ELSE
        call lsquit('xcfun version of II_DFT_KSMELDAUNRES not implemented',-1)
     ENDIF
#endif
   ENDIF
 END DO
ENDDO

END SUBROUTINE II_DFT_KSMELDAUNRES

!> \brief main closed shell GGA kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DFT_KSMEGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,GAO,&
     &                    RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,&
     &                    GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk,D3 = 3.0E0_realk,D05 = 0.5E0_realk
INTEGER     :: IPNT,I,J,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT
REAL(REALK) :: VX(5),DFTENE,GRD,GRDA,A,Econt
REAL(REALK) :: XCFUNINPUT(4,1),XCFUNOUTPUT(1,1)
REAL(REALK),pointer :: VXC(:,:,:)
EXTERNAL DFTENE
!     GGA Exchange-correlation contribution to Kohn-Sham matrix
DO IDMAT=1,NDMAT
 Econt = 0.0E0_realk
 DO IPNT = 1, NBLEN
   GRD = SQRT(GRAD(1,IPNT,IDMAT)*GRAD(1,IPNT,IDMAT)+GRAD(2,IPNT,IDMAT)*GRAD(2,IPNT,IDMAT)&
        &+GRAD(3,IPNT,IDMAT)*GRAD(3,IPNT,IDMAT))
   IF(GRD .GT. RHOTHR .OR. RHO(IPNT,IDMAT).GT.RHOTHR) THEN
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         Econt = Econt + DFTENE(RHO(IPNT,IDMAT),GRD)*WGHT(IPNT)
#ifdef VAR_XCFUN
      ELSE
         XCFUNINPUT(1,1) = RHO(IPNT,IDMAT)
!         XCFUNINPUT(2,1) = GRAD(1,IPNT,IDMAT)*GRAD(1,IPNT,IDMAT)&
!              &+GRAD(2,IPNT,IDMAT)*GRAD(2,IPNT,IDMAT)&
!              &+GRAD(3,IPNT,IDMAT)*GRAD(3,IPNT,IDMAT)
         XCFUNINPUT(2,1) = GRAD(1,IPNT,IDMAT)
         XCFUNINPUT(3,1) = GRAD(2,IPNT,IDMAT)
         XCFUNINPUT(4,1) = GRAD(3,IPNT,IDMAT)
!         call xcfun_gga_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
         call xcfun_gga_components_xc_single_eval(XCFUNINPUT,1,XCFUNOUTPUT,0)
         DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + XCFUNOUTPUT(1,1)*WGHT(IPNT)
      ENDIF
#endif
   END IF
 END DO
 DFTDATA%ENERGY(IDMAT) = DFTDATA%ENERGY(IDMAT) + Econt
ENDDO
END SUBROUTINE II_DFT_KSMEGGA

!> \brief main unrestricted GGA kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DFT_KSMEGGAUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,NDMAT,DMAT,NTYPSO,&
     & GAO,RHO,GRAD,TAU,MXBLLEN,COORD,WGHT,DFTDATA,sharedDFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,MaxNactbast)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Max Number of active basis functions
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NACTBAST,NACTBAST,NDMAT)
!> The number of gaussian atomic orbitals ( and geometrical deriv, london derivatives)
INTEGER     :: NTYPSO
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAO(NBLEN,NACTBAST,NTYPSO)
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> the gradient of electron density for all gridpoints
REAL(REALK),intent(in) :: GRAD(3,MXBLLEN,NDMAT)
!> the grad dot grad tau
REAL(REALK),intent(in) :: TAU(MXBLLEN,NDMAT)
!> max number of gridpoints
INTEGER     :: MXBLLEN
!> grippoint coordinates
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> grippoint weights
REAL(REALK),intent(in) :: WGHT(NBLEN)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> contains all info required which is not directly related to the integration (OpenMP shared)
TYPE(DFTDATATYPE),intent(inout) :: sharedDFTDATA
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> length of tmp array
integer,intent(in)        :: WORKLENGTH
!> tmp array to avoid allocation and deallocation of mem
REAL(REALK),intent(inout) :: WORK(WORKLENGTH)
!> Maximum gaussian atomic orbital values for each basis function
REAL(REALK),intent(in) :: GAOGMX(MAXNACTBAST)
!> Maximum gaussian atomic orbital value (non differentiated)
REAL(REALK),intent(in) :: GAOMAX
!
Real(realk), parameter :: D2 = 2.0E0_realk,DUMMY = 0E0_realk
INTEGER     :: IPNT,I,J,W1,W2,W3,W4,W5,W6,W7,W8,IDMAT,IDMAT1,IDMAT2
REAL(REALK) :: VXC(NBLEN,2,2),VXM(NBLEN),VX(5),DFTENEUNRES,GRDA,GRDB
EXTERNAL DFTENEUNRES

!     GGA Exchange-correlation contribution to Kohn-Sham matrix
DO IDMAT = 1,NDMAT
 IDMAT1 = 1 + (IDMAT-1)*2
 IDMAT2 = 2 + (IDMAT-1)*2
 DO IPNT = 1, NBLEN
   GRDA = SQRT(GRAD(1,IPNT,IDMAT1)*GRAD(1,IPNT,IDMAT1)+GRAD(2,IPNT,IDMAT1)*GRAD(2,IPNT,IDMAT1)+&
        &GRAD(3,IPNT,IDMAT1)*GRAD(3,IPNT,IDMAT1))
   GRDB = SQRT(GRAD(1,IPNT,IDMAT2)*GRAD(1,IPNT,IDMAT2)+GRAD(2,IPNT,IDMAT2)*GRAD(2,IPNT,IDMAT2)+&
        &GRAD(3,IPNT,IDMAT2)*GRAD(3,IPNT,IDMAT2))
   IF((GRDA .GT. RHOTHR .OR. RHO(IPNT,IDMAT1).GT.RHOTHR).OR.&
        &(GRDB .GT. RHOTHR .OR. RHO(IPNT,IDMAT2).GT.RHOTHR))then 
      IF(GRDA.LT. 1E-40_realk) GRDA = 1E-40_realk
      IF(GRDB.LT. 1E-40_realk) GRDB = 1E-40_realk
#ifdef VAR_XCFUN
      IF(.NOT.USEXCFUN)THEN
#endif
         !      CALL dft_funcderiv1unres(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),GRDA,GRDB,WGHT(IPNT),VX)
         DFTDATA%ENERGY(IDMAT1) = DFTDATA%ENERGY(IDMAT1) + DFTENEUNRES(RHO(IPNT,IDMAT1),RHO(IPNT,IDMAT2),GRDA,GRDB)*WGHT(IPNT)
#ifdef VAR_XCFUN
      ELSE
         call lsquit('xcfun version of II_DFT_KSMEGGAUNRES not implemented',-1)
      ENDIF
#endif
   END IF
 END DO
ENDDO
END SUBROUTINE II_DFT_KSMEGGAUNRES

!=====================================================================================
!
! Routines that distribute contributions to Result like Fock Matrix or gradient
! And Routines that calculated expectation values. 
!
!=====================================================================================

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> routine that does the actual work see II_dft_ksm.tex
!>
SUBROUTINE II_DFT_DIST_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     &                    COEF,GAOS,EXCMAT,DFTHRI,GAORED,EXCRED,GAOGMX,TMP)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
! TEMPORARY MEM FROM WORK
REAL(REALK),intent(inout) :: TMP(NBLEN,NACTBAST) 
REAL(REALK),intent(inout) :: GAOGMX(NACTBAST)
REAL(REALK),intent(inout) :: GAORED(NBLEN,NactBAST)
REAL(REALK),intent(inout) :: EXCRED(Nactbast*NactBAST)
!REAL(REALK) :: GAOGMX(NACTBAST)
!REAL(REALK) :: GAORED(NBLEN,NactBAST)
!REAL(REALK),pointer :: TMP(:,:)!NBLEN,NACTBAST) 
!REAL(REALK),pointer :: EXCRED(:,:)!Nactbast,NactBAST)
!
REAL(REALK) :: GAOMAX
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
INTEGER     :: INXRED(NACTBAST),IRED,JRED,offset

NRED = 0 
GAOMAX = 0.0E0_realk
!        Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0E0_realk
      DO K = 1, NBLEN
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J)))
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J)))
      ENDDO
   ENDDO
ENDDO
!        Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J) 
         DO K = 1, NBLEN
            GAORED(K,NRED) = GAOS(K,J)
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT. 0) THEN
   !  First half-contraction of GAO's with potential
!   call mem_dft_alloc(TMP,NBLEN,NRED)
    DO J=1,NRED
     DO K=1, NBLEN
      TMP(K,J) =  coef(K)* GAORED(K,J)
     ENDDO
    ENDDO
!   call mem_dft_alloc(EXCRED,NRED,NRED)
   !  Second half-contraction of GAO's with potential
!   call mem_dft_alloc(EXCRED,NRED,NRED)
    IF(XCintNoOMP)THEN
       CALL DGEMM('T','N',NRED,NRED,NBLEN,1E0_realk,&
            &                GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,&
            &                EXCRED(1:NRED*NRED),NRED)
    ELSE !Use Thread Safe version 
       CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1E0_realk,&
            &                GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,&
            &                EXCRED(1:NRED*NRED),NRED)
    ENDIF
   !  Distribute contributions to KS-matrix

   DO JRED=1,NRED         !Jred is reduced index
      J = INXRED(JRED)    !J is orbitalindex
      offset = (JRED-1)*NRED
      DO IRED=1,NRED      !Ired is reduced index
         I = INXRED(IRED) !I is orbitalindex
!$OMP ATOMIC
         EXCMAT(I,J) = EXCMAT(I,J) + EXCRED(IRED+offset)
      ENDDO
   ENDDO
!   call mem_dft_dealloc(TMP)
!   call mem_dft_dealloc(EXCRED)
ENDIF

END SUBROUTINE II_DFT_DIST_LDA

!> \brief a distribution routine 
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DISTGGABUNRES(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     &                    COEF1,COEFA,COEFB,GAOS,GRADA,GRADB,EXCMAT,DFTHRI,GAORED,EXCRED,GAOGMX,TMP)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEFA(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEFB(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,4)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
!> tmp 
REAL(REALK),intent(inout) :: EXCRED(NactBast*NactBast)
REAL(REALK),intent(in) :: GRADA(3,NBLEN),GRADB(3,NBLEN)
REAL(REALK),intent(inout) :: TMP(NBLEN,NACTBAST)
REAL(REALK),intent(inout) :: GAORED(NBLEN,NACTBAST,4),GAOGMX(NACTBAST)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,NRED,offset
REAL(REALK) :: GAOMAX
!REAL(REALK),pointer :: EXCRED(:,:)
!REAL(REALK) :: EXCRED(NactBast*NactBast)
NRED = 0
GAOMAX = 0.0E0_realk
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0E0_realk
      DO K = 1, NBLEN
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
         DO I=1,4
            GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,I)))
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO K = 1, NBLEN
            GAORED(K,NRED,1) = GAOS(K,J,1)
            GAORED(K,NRED,2) = GAOS(K,J,2)
            GAORED(K,NRED,3) = GAOS(K,J,3)
            GAORED(K,NRED,4) = GAOS(K,J,4)
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT. 0) THEN
!  First half-contraction of GAO's with potential
   DO J=1,NRED !J is reduced index
      DO K=1, NBLEN
         TMP(K,J) =  coef1(K)* GAORED(K,J,1)&
              &    + coefA(K)*(GAORED(K,J,2)*GRADA(1,K)+&
              &                GAORED(K,J,3)*GRADA(2,K)+&
              &                GAORED(K,J,4)*GRADA(3,K))&
              &    + coefB(K)*(GAORED(K,J,2)*GRADB(1,K)+&
              &                GAORED(K,J,3)*GRADB(2,K)+&
              &                GAORED(K,J,4)*GRADB(3,K))
      ENDDO
   ENDDO
!  Second half-contraction of GAO's with potential
!   CALL MEM_DFT_ALLOC(EXCRED,NRED,NRED)
   IF(XCintNoOMP)THEN
      CALL DGEMM('T','N',NRED,NRED,NBLEN,1.0E0_realk,&
           &                GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,&
           &                EXCRED(1:NRED*NRED),NRED)
   ELSE !Use Thread Safe version 
      CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1.0E0_realk,&
           &                GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,&
           &                EXCRED(1:NRED*NRED),NRED)
   ENDIF
!  Distribute contributions to KS-matrix
   DO JRED=1,NRED         !Jred is reduced index
      J = INXRED(JRED)    !J is orbital index
      offset = (JRED-1)*NRED
      DO IRED=1,NRED      !Ired is reduced index
         I = INXRED(IRED) !I is orbital index
!$OMP ATOMIC
         EXCMAT(I,J) = EXCMAT(I,J) + EXCRED(IRED+offset)
      ENDDO
   ENDDO
!   CALL MEM_DFT_DEALLOC(EXCRED)
ENDIF

END SUBROUTINE II_DISTGGABUNRES

!> \brief the main GGA distribution routine 
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DISTGGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     & NTYPSO,NDMAT,COEF,GAOS,EXCMAT,DFTHRI,GAORED,EXCRED,GAOGMX,GAOMAX,&
     & TMP,MAXNACTBAST)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Maximum number of active basis functions
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of GAOs types
INTEGER,intent(in) :: NTYPSO
!> Number of matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(4,NBLEN,NDMAT)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,NDMAT)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK),intent(inout) :: TMP(NBLEN,NACTBAST)
REAL(REALK),intent(in)    :: GAOGMX(MAXNACTBAST),GAOMAX
REAL(REALK),intent(inout) :: GAORED(NBLEN,NACTBAST,4)
REAL(REALK),intent(inout) :: EXCRED(NactBAST*NactBAST)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN
INTEGER     :: IRED,JRED,NRED,offset,IDMAT
INTEGER,pointer :: INXRED(:)
call mem_dft_alloc(INXRED,NACTBAST)
NRED = 0
! TODO: change the NRED from being individual to be BLOCKS => more efficient EXCRED
 
! Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO I=1,4
            DO K = 1, NBLEN
               GAORED(K,NRED,I) = GAOS(K,J,I)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT. 0) THEN
   DO IDMAT = 1,NDMAT
      !  First half-contraction of GAO's with potential
      DO J=1,NRED
         DO K=1, NBLEN
            TMP(K,J) =  coef(1,K,IDMAT)*GAORED(K,J,1)&
                 &   + coef(2,K,IDMAT)*GAORED(K,J,2)&
                 &   + coef(3,K,IDMAT)*GAORED(K,J,3)&
                 &   + coef(4,K,IDMAT)*GAORED(K,J,4)
         ENDDO
      ENDDO
      !  Second half-contraction of GAO's with potential
      !   CALL MEM_DFT_ALLOC(EXCRED,NRED,NRED)
      IF(XCintNoOMP)THEN
         CALL DGEMM('T','N',NRED,NRED,NBLEN,1.0E0_realk,&
              & GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,EXCRED(1:NRED*NRED),NRED)
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1.0E0_realk,&
              & GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,EXCRED(1:NRED*NRED),NRED)
      ENDIF
      !  Distribute contributions to KS-matrix
      DO JRED=1,NRED          !Jred is reduced index
         J = INXRED(JRED)     !J is orbital index
         offset = (JRED-1)*NRED
         DO IRED=1,NRED       !Ired is reduced index
            I = INXRED(IRED)  !I is orbital index
!$OMP ATOMIC
            EXCMAT(I,J,IDMAT) = EXCMAT(I,J,IDMAT) + EXCRED(IRED+offset)
         ENDDO
      ENDDO
      !   CALL MEM_DFT_DEALLOC(EXCRED)
   ENDDO
ENDIF
call mem_dft_dealloc(INXRED)

END SUBROUTINE II_DISTGGA

!> \brief the main GGA distribution routine 
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_DISTMETA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,Nactbast,NBAST,&
     & NTYPSO,NDMAT,COEF,GAOS,EXCMAT,DFTHRI,GAORED,EXCRED,GAOGMX,TMP)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of GAOs types
INTEGER,intent(in) :: NTYPSO
!> Number of matrices
INTEGER,intent(in) :: NDMAT
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(5,NBLEN,NDMAT)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,NDMAT)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!
REAL(REALK),intent(inout) :: TMP(NBLEN,NACTBAST)
REAL(REALK),intent(inout) :: GAOGMX(NACTBAST)
REAL(REALK),intent(inout) :: GAORED(NBLEN,NACTBAST,4)
REAL(REALK),intent(inout) :: EXCRED(NactBAST*NactBAST)
REAL(REALK) :: GAOMAX,GAOMAXTMP
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN
INTEGER     :: IRED,JRED,NRED,offset,IDMAT
INTEGER,pointer :: INXRED(:)
call mem_dft_alloc(INXRED,NACTBAST)
DO IDMAT = 1,NDMAT
 NRED = 0
 GAOMAX = 0.0E0_realk
 ! Set up maximum Gaussian AO elements
 DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOMAXTMP = 0.0E0_realk
      DO K = 1, NBLEN
         GAOMAXTMP = MAX(GAOMAXTMP,ABS(GAOS(K,J,1)))
      ENDDO
      GAOMAX = MAX(GAOMAX,GAOMAXTMP)
      DO I=2,4
         DO K = 1, NBLEN
            GAOMAXTMP = MAX(GAOMAXTMP,ABS(GAOS(K,J,I)))
         ENDDO
      ENDDO
      GAOGMX(J) = GAOMAXTMP      
   ENDDO
 ENDDO

 ! Set up reduced Gaussian AO's
 DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO I=1,4
            DO K = 1, NBLEN
               GAORED(K,NRED,I) = GAOS(K,J,I)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
 ENDDO
 IF (NRED.GT. 0) THEN
 !  First half-contraction of GAO's with potential
   DO J=1,NRED
      DO K=1, NBLEN
         TMP(K,J) =  coef(1,K,IDMAT)*GAORED(K,J,1)&
              &   + coef(2,K,IDMAT)*GAORED(K,J,2)&
              &   + coef(3,K,IDMAT)*GAORED(K,J,3)&
              &   + coef(4,K,IDMAT)*GAORED(K,J,4)
      ENDDO
   ENDDO
!  Second half-contraction of GAO's with potential
!   CALL MEM_DFT_ALLOC(EXCRED,NRED,NRED)
   IF(XCintNoOMP)THEN
      CALL DGEMM('T','N',NRED,NRED,NBLEN,1.0E0_realk,&
           &                GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,&
           &                EXCRED(1:NRED*NRED),NRED)
   ELSE !Use Thread Safe version 
      CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1.0E0_realk,&
           &                GAORED,NBLEN,TMP,NBLEN,0.0E0_realk,&
           &                EXCRED(1:NRED*NRED),NRED)
   ENDIF
!  Distribute contributions to KS-matrix
   DO JRED=1,NRED          !Jred is reduced index
      J = INXRED(JRED)     !J is orbital index
      offset = (JRED-1)*NRED
      DO IRED=1,NRED       !Ired is reduced index
         I = INXRED(IRED)  !I is orbital index
!$OMP ATOMIC
         EXCMAT(I,J,IDMAT) = EXCMAT(I,J,IDMAT) + EXCRED(IRED+offset)
      ENDDO
   ENDDO
!   CALL MEM_DFT_DEALLOC(EXCRED)
 ENDIF
ENDDO
call mem_dft_dealloc(INXRED)

END SUBROUTINE II_DISTMETA

SUBROUTINE LB94correction(rho,GRD,HFexchangeFac,WGHT,VX)
implicit none
REAL(REALK),intent(in)    :: rho,GRD,HFexchangeFac,WGHT
REAL(REALK),intent(inout) :: VX
Real(realk), parameter :: D3 = 3.0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D005 = 0.05E0_realk,D1=1.0E0_realk,D08=0.8E0_realk
REAL(REALK) :: scaled_grad,rho13,rho43,sg2,gdenom,g

!LB94 correction
rho13 = rho**(1.d0/3.d0)
rho43 = rho13*rho
IF(rho43.GT.1.0E-13_realk)THEN
   scaled_grad = GRD/rho43
ELSE
   scaled_grad = GRD/1.0E-13_realk
ENDIF
!BETA=D005
sg2 = scaled_grad*scaled_grad
gdenom = D1 + D3*D005*scaled_grad*log(scaled_grad+dsqrt(1d0+scaled_grad*scaled_grad))
g = -D005*sg2/gdenom
VX = VX + (D1-HFexchangeFac)*rho13*g*WGHT
end SUBROUTINE LB94correction

SUBROUTINE CS00correction(rho,GRD,HFexchangeFac,WGHT,VX,CS00SHIFT,CS00eHOMO,CS00ZND1,CS00ZND2)
implicit none
REAL(REALK),intent(in)    :: rho,GRD,HFexchangeFac,WGHT,CS00SHIFT,CS00eHOMO,CS00ZND1,CS00ZND2
REAL(REALK),intent(inout) :: VX
Real(realk), parameter :: D3 = 3.0E0_realk,D05 = 0.5E0_realk
Real(realk), parameter :: D005 = 0.05E0_realk,D1=1.0E0_realk,D08=0.8E0_realk
REAL(REALK) :: scaled_grad,rho13,rho43,sg2,gdenom,g,delta,P1,P2

!LB94 correction
IF(ABS(CS00SHIFT).LT.1.0E-12_realk)THEN
   delta = -CS00ZND1*CS00eHOMO + CS00ZND2
ELSE
   delta = CS00SHIFT
ENDIF

rho13 = rho**(1.d0/3.d0)
rho43 = rho13*rho
IF(rho43.GT.1.0E-13_realk)THEN
   scaled_grad = GRD/rho43
ELSE
   scaled_grad = GRD/1.0E-13_realk
ENDIF
!BETA=D005
sg2 = scaled_grad*scaled_grad
gdenom = D1 + D3*D005*scaled_grad*log(scaled_grad+dsqrt(1d0+scaled_grad*scaled_grad))
g = -D005*sg2/gdenom
P1 = VX + (D1-HFexchangeFac)*rho13*g*WGHT
P2 = VX - D05*delta*WGHT
VX = MAX(P1,P2)
end SUBROUTINE CS00correction

!> \brief computes the expectation value of a matrix BMAT
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_expval_lda(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,EXPVAL,BMAT,nBmat,DFTHRI,NRED)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,1)
!> The expectation value of the B matrix
REAL(REALK),intent(inout) :: EXPVAL(NBLEN,NBMAT)
!> The B matrix (some perturbed density matrix) 
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in) :: NBMAT
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of reduced orbitals
INTEGER,intent(inout) :: NRED
!
REAL(REALK) :: GAOMAX,BMAX
REAL(REALK),pointer :: GAORED(:,:),GAOGMX(:),TMP(:,:)
INTEGER,pointer     :: INXRED(:)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN
INTEGER     :: IRED,JRED,IORB,JORB,IBMAT
real(realk),parameter :: D05 = 0.5E0_realk
REAL(REALK),pointer :: BRED(:,:)

call mem_dft_alloc(GAORED,NBLEN,NACTBAST)
call mem_dft_alloc(GAOGMX,NACTBAST)
call mem_dft_alloc(TMP,NBLEN,NACTBAST)
call mem_dft_alloc(INXRED,NACTBAST)


!IF(NBMAT .GT. 1)CALL LSQUIT('II_get_expval_lda does not work correctly for more than 1 nbmat',lupri)
NRED = 0
GAOMAX = 0.0E0_realk
! Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0E0_realk
      DO K = 1, NBLEN
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,1)))
      ENDDO
      GAOMAX = MAX(GAOMAX,GAOGMX(J))
   ENDDO
ENDDO

! Set up maximum response-vector elements
BMAX = 0.0E0_realk
DO IBMAT=1,NBMAT
DO JBL=1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL) !J is active index
      JORB = INXACT(J)                 !JORB is orbital index 
      DO IBL=1, NBLOCKS
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)  !I is active index
            IORB = INXACT(I)                  !IORB is orbital index 
            BMAX = MAX(BMAX,ABS(BMAT(IORB,JORB,IBMAT)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
ENDDO

! Set up reduced Gaussian AO's
NRED=0
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL) !J is active index
      IF (GAOGMX(J)*GAOMAX*BMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO K = 1, NBLEN 
            GAORED(K,NRED) = GAOS(K,J,1)
         ENDDO
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT. 0) THEN
   ! Set up reduced density-matrix
   call mem_dft_alloc(BRED,NRED,NRED)
   DO IBMAT=1,NBMAT
    DO JRED=1,NRED         !Jred is reduced index 
      J = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED       !Ired is reduced index 
         I = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(I,J,IBMAT) + BMAT(J,I,IBMAT)
      ENDDO
    ENDDO
    ! First half-contraction of Gaussian AO with density-matrix
    IF(XCintNoOMP)THEN
       CALL DGEMM("N","N",NBLEN,NRED,NRED,1.0E0_realk,GAORED,NBLEN,&
            &     BRED,NRED,0.0E0_realk,TMP,NBLEN)
    ELSE !Use Thread Safe version 
       CALL DGEMM_TS("N","N",NBLEN,NRED,NRED,1.0E0_realk,GAORED,NBLEN,&
            &     BRED,NRED,0.0E0_realk,TMP,NBLEN)
    ENDIF
    ! Second half-contraction 
    DO K = 1, NBLEN
       EXPVAL(K,IBMAT) = GAORED(K,1)*TMP(K,1)
    END DO
    DO I = 2, NRED  !I is reduced index
       DO K = 1, NBLEN
          EXPVAL(K,IBMAT) = EXPVAL(K,IBMAT) + GAORED(K,I)*TMP(K,I)
       END DO
    END DO
    DO K = 1, NBLEN
       EXPVAL(K,IBMAT) = D05*EXPVAL(K,IBMAT)
    END DO
   ENDDO
   call mem_dft_dealloc(BRED)
ENDIF

call mem_dft_dealloc(GAORED)
call mem_dft_dealloc(GAOGMX)
call mem_dft_dealloc(TMP)
call mem_dft_dealloc(INXRED)

END SUBROUTINE II_GET_EXPVAL_LDA

!> \brief computes the expectation value and and gradient of it
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_expval_gga(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,EXPVAL,EXPGRAD,BMAT,nBMAT,DFTHRI,NRED,&
     &GAORED1,GAORED2,GAOGMX,GAOMAX,TMP,MAXNACTBAST)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Maximum number of active basis functions
INTEGER,intent(in) :: MaxNactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,4)
!> The expectation value of the B MATRIX
REAL(REALK),intent(inout) :: EXPVAL(NBLEN,NBMAT)
!> The gradient of the expectation value of the B MATRIX
REAL(REALK),intent(inout) :: EXPGRAD(3,NBLEN,NBMAT)
!> The B matrix (some perturbed density matrix) 
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in) :: NBMAT
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of reduced orbitals
INTEGER,intent(inout) :: NRED
REAL(REALK),intent(in) :: GAOMAX,GAOGMX(MAXNACTBAST)
!
REAL(REALK) :: BREDIJ,BREDJI,BREDCOMBI,BMAX,TMPGAOMAX
REAL(REALK),intent(inout) :: GAORED1(NBLEN,NACTBAST),GAORED2(NBLEN,NACTBAST,3),TMP(NBLEN,NACTBAST)
integer     :: INXRED(NACTBAST)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN
INTEGER     :: IRED,JRED,IORB,JORB,IBMAT
REAL(REALK),parameter :: D05 = 0.5E0_realk,D2 = 2.0E0_realk
REAL(REALK),pointer :: BRED(:,:)
INTEGER     :: ORBBLOCKS(2,NBLOCKS)

! Set up maximum Gaussian AO elements
DO IBL=1, NBLOCKS
   ORBBLOCKS(1,IBL) = INXACT(BLOCKS(1,IBL))
   ORBBLOCKS(2,IBL) = INXACT(BLOCKS(2,IBL))
ENDDO

! Set up maximum response-vector elements
BMAX = 0.0E0_realk
DO IBMAT=1,NBMAT
 DO JBL=1, NBLOCKS
  DO IBL=1, NBLOCKS
   DO JORB = ORBBLOCKS(1,JBL), ORBBLOCKS(2,JBL)  
    DO IORB = ORBBLOCKS(1,IBL), ORBBLOCKS(2,IBL) 
     BMAX = MAX(BMAX,ABS(BMAT(IORB,JORB,IBMAT)))
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
! Set up reduced Gaussian AO's
NRED=0
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX*BMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO K = 1, NBLEN
            GAORED1(K,NRED) = GAOS(K,J,1)
            GAORED2(K,NRED,1) = GAOS(K,J,2)
            GAORED2(K,NRED,2) = GAOS(K,J,3)
            GAORED2(K,NRED,3) = GAOS(K,J,4)
         ENDDO
      ENDIF
   ENDDO
ENDDO
   
IF (NRED.GT. 0) THEN
   ! Set up reduced density-matrix
   call mem_dft_alloc(BRED,NRED,NRED)
   DO IBMAT=1,NBMAT
    DO JRED=1,NRED         !Jred is reduced index 
      J = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED       !Ired is reduced index 
         I = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(I,J,IBMAT) + BMAT(J,I,IBMAT)
      ENDDO
    ENDDO
    ! First half-contraction of Gaussian AO with density-matrix
    IF(XCintNoOMP)THEN
       CALL DGEMM("N","N",NBLEN,NRED,NRED,1.0E0_realk,GAORED1,NBLEN,&
            &     BRED,NRED,0.0E0_realk,TMP,NBLEN)
    ELSE !Use Thread Safe version 
       CALL DGEMM_TS("N","N",NBLEN,NRED,NRED,1.0E0_realk,GAORED1,NBLEN,&
            &     BRED,NRED,0.0E0_realk,TMP,NBLEN)
    ENDIF
    ! Second half-contraction 
    DO K = 1, NBLEN
       EXPVAL(K,IBMAT) = GAORED1(K,1)*TMP(K,1)
    END DO
    DO I = 2, NRED  !I is reduced index
       DO K = 1, NBLEN
          EXPVAL(K,IBMAT) = EXPVAL(K,IBMAT) + GAORED1(K,I)*TMP(K,I)
       END DO
    END DO
    DO K = 1, NBLEN
       EXPVAL(K,IBMAT) = D05*EXPVAL(K,IBMAT)
    END DO
    DO K = 1, NBLEN
       EXPGRAD(1,K,IBMAT) = GAORED2(K,1,1)*TMP(K,1)
       EXPGRAD(2,K,IBMAT) = GAORED2(K,1,2)*TMP(K,1)
       EXPGRAD(3,K,IBMAT) = GAORED2(K,1,3)*TMP(K,1)
    ENDDO
    DO I = 2, NRED
       DO K = 1, NBLEN
          EXPGRAD(1,K,IBMAT) = EXPGRAD(1,K,IBMAT) + GAORED2(K,I,1)*TMP(K,I)
          EXPGRAD(2,K,IBMAT) = EXPGRAD(2,K,IBMAT) + GAORED2(K,I,2)*TMP(K,I)
          EXPGRAD(3,K,IBMAT) = EXPGRAD(3,K,IBMAT) + GAORED2(K,I,3)*TMP(K,I)
       ENDDO
    ENDDO
   ENDDO
   call mem_dft_dealloc(BRED)
ENDIF

END SUBROUTINE II_GET_EXPVAL_GGA

!> \brief distribution routine for magnetic derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF,GAOS,EXCMAT,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
REAL(REALK),pointer :: TMP(:,:,:)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST),GAOMAX
REAL(REALK),pointer :: GAORED(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),pointer :: EXCRED(:,:,:)
REAL(REALK) :: Rbeta,Rgamma
real(realk),parameter :: D2=2E0_realk

call mem_dft_alloc(GAORED,NBLEN,NACTBAST,4)
call mem_dft_alloc(TMP,NBLEN,NACTBAST,3)
NRED = 0 
GAOMAX = 0.0E0_realk
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0E0_realk
      DO K = 1, NBLEN
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,N)))
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
      ENDDO
   ENDDO
 ENDDO
ENDDO
!        Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = INXACT(J) 
   DO N=1,NTYPSO
    DO K = 1, NBLEN
       GAORED(K,NRED,N) = GAOS(K,J,N)
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT. 0) THEN
   !  First half-contraction of GAO's with potential
 DO ALPHA = 1,3
  beta = betaList(ALPHA)
  gamma= gammaList(ALPHA)
  DO J=1,NRED
   DO K=1, NBLEN
    Rgamma = COORD(gamma,K) !- Origin(gamma)       
    Rbeta = COORD(beta,K) !- Origin(beta)
    TMP(K,J,ALPHA) = coef(K)*GAORED(K,J,1+beta)*Rgamma & 
         & - coef(K)*GAORED(K,J,1+gamma)*Rbeta
   ENDDO
  ENDDO
 ENDDO
 !  Second half-contraction of GAO's
 call mem_dft_alloc(EXCRED,NRED,NRED,3)
 DO ALPHA = 1,3
    IF(XCintNoOMP)THEN
       CALL DGEMM('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
            &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
    ELSE !Use Thread Safe version 
       CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
            &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
    ENDIF
 ENDDO
 DO ALPHA=1,3
  DO JRED=1,NRED    !Jred is reduced index
   J = INXRED(JRED) !J is orbital index
   DO IRED=1,NRED   !Ired is reduced index
    I = INXRED(IRED)!I is orbital index 
!$OMP ATOMIC
    EXCMAT(I,J,ALPHA) = EXCMAT(I,J,ALPHA) + EXCRED(IRED,JRED,ALPHA)
   ENDDO
  ENDDO
 ENDDO
 call mem_dft_dealloc(EXCRED)
ENDIF
call mem_dft_dealloc(GAORED)
call mem_dft_dealloc(TMP)

END SUBROUTINE II_DFT_DISTMAGDERIV_LDA

!> \brief distribution routine for magnetic derivative Kohn-sham matrix GGA
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF,GAOS,EXCMAT,COORD,DFTHRI,NTYPSO,GAORED,TMP,GAOGMX)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(4,NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> TMP 
REAL(REALK),intent(inout) :: GAORED(NBLEN,NACTBAST,NTYPSO)
REAL(REALK),intent(inout) :: TMP(NBLEN,NACTBAST,3)
REAL(REALK),intent(inout) :: GAOGMX(NACTBAST)
!
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
INTEGER     :: beta1,gamma1,beta2X,beta2Y,beta2Z,gamma2X,gamma2Y,gamma2Z
INTEGER     :: IRED,JRED,alpha,beta,gamma,N,COORDINATE
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK) :: Rbeta,Rgamma,GAOMAX
real(realk),parameter :: D2=2E0_realk,D05=0.5E0_realk
INTEGER :: INXRED(NACTBAST)
REAL(REALK),pointer :: EXCRED(:,:,:)
NRED = 0 
GAOMAX = 0.0E0_realk
!        Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOGMX(J) = 0.0E0_realk
  DO N=1,NTYPSO
   DO K = 1, NBLEN
    GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,N)))
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
!        Set up reduced Gaussian AO's
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = INXACT(J) 
   DO N=1,NTYPSO
    DO K = 1, NBLEN
       GAORED(K,NRED,N) = GAOS(K,J,N)
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT. 0) THEN
 DO ALPHA = 1,3
  beta = betaList(ALPHA)
  gamma= gammaList(ALPHA)
  beta1 = 4+betaList(ALPHA)
  gamma1= 4+gammaList(ALPHA)
  beta2X = 7+(betaList(ALPHA)-1)*3+1
  beta2Y = 7+(betaList(ALPHA)-1)*3+2
  beta2Z = 7+(betaList(ALPHA)-1)*3+3
  gamma2X= 7+(gammaList(ALPHA)-1)*3+1
  gamma2Y= 7+(gammaList(ALPHA)-1)*3+2
  gamma2Z= 7+(gammaList(ALPHA)-1)*3+3
  DO J=1,NRED
   DO K=1, NBLEN
    Rgamma = COORD(gamma,K) !- Origin(gamma)       
    Rbeta = COORD(beta,K) !- Origin(beta)
    TMP(K,J,ALPHA) = &
         &+(D2*coef(1,K)*Rgamma+coef(1+gamma,K))*GAORED(K,J,beta1)&
         &-(D2*coef(1,K)*Rbeta+coef(1+beta,K))*GAORED(K,J,gamma1)&
         &-(GAORED(K,J,gamma2X)*Rbeta-GAORED(K,J,beta2X)*Rgamma)*coef(2,K)&
         &-(GAORED(K,J,gamma2Y)*Rbeta-GAORED(K,J,beta2Y)*Rgamma)*coef(3,K)&
         &-(GAORED(K,J,gamma2Z)*Rbeta-GAORED(K,J,beta2Z)*Rgamma)*coef(4,K)
   ENDDO
  ENDDO
 ENDDO
 call mem_dft_alloc(EXCRED,NRED,NRED,3)
 DO ALPHA = 1,3
    IF(XCintNoOMP)THEN
       CALL DGEMM('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
            &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
    ELSE !Use Thread Safe version 
       CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
            &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
    ENDIF
 ENDDO
 DO COORDINATE=2,4
    DO ALPHA = 1,3
     beta = betaList(ALPHA)
     gamma= gammaList(ALPHA)
     beta1 = 4+betaList(ALPHA)
     gamma1= 4+gammaList(ALPHA)
     DO J=1,NRED
      DO K=1, NBLEN
       Rgamma = COORD(gamma,K) !- Origin(gamma)       
       Rbeta = COORD(beta,K) !- Origin(beta)
       TMP(K,J,ALPHA) = (-GAORED(K,J,gamma1)*Rbeta+GAORED(K,J,beta1)*Rgamma)&
            &*coef(COORDINATE,K)
      ENDDO
     ENDDO
    ENDDO
    DO ALPHA = 1,3
       IF(XCintNoOMP)THEN
          CALL DGEMM('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,COORDINATE),NBLEN,&
               &     TMP(:,:,ALPHA),NBLEN,1.0E0_realk,EXCRED(:,:,ALPHA),NRED)
       ELSE !Use Thread Safe version 
          CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,COORDINATE),NBLEN,&
               &     TMP(:,:,ALPHA),NBLEN,1.0E0_realk,EXCRED(:,:,ALPHA),NRED)
       ENDIF
    ENDDO
 ENDDO
 DO ALPHA=1,3
  DO JRED=1,NRED     !Jred is reduced index
   J = INXRED(JRED)  !J is orbital index
   DO IRED=1,NRED    !Ired is reduced index
    I = INXRED(IRED) !I is orbital index
!$OMP ATOMIC
    EXCMAT(I,J,ALPHA) = EXCMAT(I,J,ALPHA) + EXCRED(IRED,JRED,ALPHA)
   ENDDO
  ENDDO
 ENDDO
 call mem_dft_dealloc(EXCRED)
ENDIF

END SUBROUTINE II_DFT_DISTMAGDERIV_GGA

!> \brief computes the expectation value and the magnetic derivative of it
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_magderiv_expval_lda(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,COORD,EXPVAL,EXPGRAD,BMAT,nBmat,NRED,DFTHRI,DOSYMPART)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in)     :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in)     :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> The expvalue
REAL(REALK),intent(inout) :: EXPVAL(NBLEN,NBMAT)
!> The magnetic gradient of expvalue
REAL(REALK),intent(inout) :: EXPGRAD(NBLEN,3,NBMAT)
!> B matrix some perturbed density matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> Number of reduced orbitals
INTEGER,intent(inout)  :: NRED
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
REAL(REALK),pointer :: TMP(:,:,:),GAORED(:,:,:)
REAL(REALK) :: GAOGMX(NACTBAST),X,Y,Z,GAOMAX,BMAX
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,INXRED(NACTBAST)
INTEGER     :: IRED,JRED,IORB,JORB,gabX,gabY,gabZ,N,IBMAT
REAL(REALK),pointer :: BRED(:,:,:)
call mem_dft_alloc(TMP,NBLEN,NACTBAST,4)
call mem_dft_alloc(GAORED,NBLEN,NACTBAST,NTYPSO)
!IF(NBMAT .GT. 1)CALL LSQUIT('II_get_expval_lda does not work correctly for more than 1 nbmat',lupri)
NRED = 0
GAOMAX = 0.0E0_realk
! Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      GAOGMX(J) = 0.0E0_realk
      DO K = 1,NBLEN
         GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,N)))
         GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
      ENDDO
   ENDDO
 ENDDO
ENDDO

! Set up maximum response-vector elements
BMAX = 0.0E0_realk
DO IBMAT=1,NBMAT
DO JBL=1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL) !J is active index 
      JORB = INXACT(J) !JORB is orbital index 
      DO IBL=1, NBLOCKS  
         DO I = BLOCKS(1,IBL), BLOCKS(2,IBL) !I is active index 
            IORB = INXACT(I) !IORB is orbital index 
            BMAX = MAX(BMAX,ABS(BMAT(IORB,JORB,IBMAT)))
         ENDDO
      ENDDO
   ENDDO
ENDDO
ENDDO

! Set up reduced Gaussian AO's
NRED=0
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      IF (GAOGMX(J)*GAOMAX*BMAX.GT.DFTHRI) THEN
         NRED = NRED + 1
         INXRED(NRED) = INXACT(J)
         DO N=1,NTYPSO
          DO K = 1, NBLEN
             GAORED(K,NRED,N) = GAOS(K,J,N)
          ENDDO 
         ENDDO
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT. 0) THEN
   ! Set up reduced density-matrix
 call mem_dft_alloc(BRED,NRED,NRED,NBMAT)
 DO IBMAT=1,NBMAT
  DO JRED=1,NRED     !Jred is reduced index 
   J = INXRED(JRED)  !J is orbital index 
   DO IRED=1,NRED    !Ired is reduced index 
    I = INXRED(IRED) !I is orbital index 
    BRED(IRED,JRED,IBMAT) = BMAT(I,J,IBMAT)
   ENDDO
  ENDDO
 ENDDO

   ! First half-contraction of Gaussian AO with density-matrix
 DO IBMAT=1,NBMAT
    IF(XCintNoOMP)THEN
       CALL DGEMM("N","N",NBLEN,NRED,NRED,1E0_realk,GAORED(:,:,1),NBLEN,&
            &     BRED(:,:,IBMAT),NRED,0.0E0_realk,TMP(:,:,1),NBLEN)
    ELSE !Use Thread Safe version 
       CALL DGEMM_TS("N","N",NBLEN,NRED,NRED,1E0_realk,GAORED(:,:,1),NBLEN,&
            &     BRED(:,:,IBMAT),NRED,0.0E0_realk,TMP(:,:,1),NBLEN)
    ENDIF
    ! Second half-contraction 
    DO K = 1, NBLEN
       EXPVAL(K,IBMAT) = GAORED(K,1,1)*TMP(K,1,1)
    END DO
    DO I = 2, NRED
       DO K = 1, NBLEN
          EXPVAL(K,IBMAT) = EXPVAL(K,IBMAT) + GAORED(K,I,1)*TMP(K,I,1)
       END DO
    END DO
 ENDDO

 IF(DOSYMPART)THEN 
    !IF MATRIX B IS SYMMETRIC THIS WILL BE ZERO
    !SO WE ONLY CALCULATE THIS IF THE B MATRIX IS NONSYMMETRIC
    !RESULTING IN A SYMMETRIC FINAL CONTRIBUTION 
   gabX = 2
   gabY = 3
   gabZ = 4
   DO IBMAT=1,NBMAT
   ! First half-contraction of Gaussian AO with density-matrix
    !TMP(K,J,N) = GAORED(K,I,N)*BRED(I,J) - TMP(K,J,1) is already built
    DO N=2,4
       IF(XCintNoOMP)THEN
          CALL DGEMM("N","N",NBLEN,NRED,NRED,1E0_realk,GAORED(:,:,N),NBLEN,&
               &     BRED(:,:,IBMAT),NRED,0.0E0_realk,TMP(:,:,N),NBLEN)
       ELSE !Use Thread Safe version 
          CALL DGEMM_TS("N","N",NBLEN,NRED,NRED,1E0_realk,GAORED(:,:,N),NBLEN,&
               &     BRED(:,:,IBMAT),NRED,0.0E0_realk,TMP(:,:,N),NBLEN)
       ENDIF
    ENDDO
   ! Second half-contraction of Gaussian AOs
    J=1
    DO K = 1, NBLEN
     Y = COORD(2,K) !- Origin(2)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,1,IBMAT) = &
          & +GAORED(K,J,1)*(Z*TMP(K,J,gabY)-Y*TMP(K,J,gabZ))&
          & +TMP(K,J,1)*(Y*GAORED(K,J,gabZ)-Z*GAORED(K,J,gabY))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        Y = COORD(2,K) !- Origin(2)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,1,IBMAT) = EXPGRAD(K,1,IBMAT)&
             & +GAORED(K,J,1)*(Z*TMP(K,J,gabY)-Y*TMP(K,J,gabZ))&
             & +TMP(K,J,1)*(Y*GAORED(K,J,gabZ)-Z*GAORED(K,J,gabY))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,2,IBMAT) = &
          & +GAORED(K,J,1)*(X*TMP(K,J,gabZ)-Z*TMP(K,J,gabX))&
          & +TMP(K,J,1)*(Z*GAORED(K,J,gabX)-X*GAORED(K,J,gabZ))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,2,IBMAT) = EXPGRAD(K,2,IBMAT)&
             & +GAORED(K,J,1)*(X*TMP(K,J,gabZ)-Z*TMP(K,J,gabX))&
             & +TMP(K,J,1)*(Z*GAORED(K,J,gabX)-X*GAORED(K,J,gabZ))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Y = COORD(2,K) !- Origin(2)      
     EXPGRAD(K,3,IBMAT) = &
          & +GAORED(K,J,1)*(Y*TMP(K,J,gabX)-X*TMP(K,J,gabY))&
          & +TMP(K,J,1)*(X*GAORED(K,J,gabY)-Y*GAORED(K,J,gabX))
    END DO
    DO J = 2, NRED
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Y = COORD(2,K) !- Origin(2)      
        EXPGRAD(K,3,IBMAT) = EXPGRAD(K,3,IBMAT)&
             & +GAORED(K,J,1)*(Y*TMP(K,J,gabX)-X*TMP(K,J,gabY))&
             & +TMP(K,J,1)*(X*GAORED(K,J,gabY)-Y*GAORED(K,J,gabX))
     END DO
    END DO
   ENDDO
 ENDIF
 call mem_dft_dealloc(BRED)
ENDIF
call mem_dft_dealloc(TMP)
call mem_dft_dealloc(GAORED)

END SUBROUTINE II_GET_MAGDERIV_EXPVAL_LDA

!> \brief A distribution routine for the magnetic derivative of the linear response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF,COEF2,NBMAT,GAOS,EXCMAT,EXCMATSYM,COORD,DFTHRI,NTYPSO,dosympart)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: nBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF(NBLEN,NBMAT)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN,NBMAT,3)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3*NBMAT)
!> The sym part magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMATSYM(NBAST,NBAST,3*NBMAT)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST),GAOMAX
REAL(REALK),pointer :: GAORED(:,:,:),TMP(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,IBMAT
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),pointer :: EXCRED(:,:,:)
REAL(REALK) :: Rbeta,Rgamma
real(realk),parameter :: D2=2E0_realk,D05=0.5E0_realk
call mem_dft_alloc(TMP,NBLEN,NACTBAST,3)
call mem_dft_alloc(GAORED,NBLEN,NACTBAST,4)
GAOMAX = 0.0E0_realk
!        Set up maximum Gaussian AO elements
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOGMX(J) = 0.0E0_realk
  DO K = 1, NBLEN
     GAOMAX = MAX(GAOMAX,ABS(GAOS(K,J,1)))
     DO N=1,NTYPSO
        GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
     ENDDO
  ENDDO
 ENDDO
ENDDO
!        Set up reduced Gaussian AO's
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = INXACT(J) 
   DO N=1,NTYPSO
    DO K = 1, NBLEN
       GAORED(K,NRED,N) = GAOS(K,J,N)
    ENDDO
   ENDDO
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT. 0) THEN
   !  First half-contraction of GAO's with potential
 call mem_dft_alloc(EXCRED,NRED,NRED,3)
 DO IBMAT = 1,NBMAT
  DO ALPHA = 1,3
   beta = betaList(ALPHA)
   gamma= gammaList(ALPHA)
   DO J=1,NRED
    DO K=1, NBLEN
     Rgamma = COORD(gamma,K) !- Origin(gamma)       
     Rbeta = COORD(beta,K) !- Origin(beta)
     TMP(K,J,ALPHA) = &
          &-coef(K,IBMAT)*GAORED(K,J,1+beta)*Rgamma &
          &+coef(K,IBMAT)*GAORED(K,J,1+gamma)*Rbeta 
    ENDDO
   ENDDO
  ENDDO
  !  Second half-contraction of GAO's
  DO ALPHA = 1,3
     IF(XCintNoOMP)THEN
        CALL DGEMM('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
             &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
     ELSE !Use Thread Safe version 
        CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
             &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
     ENDIF
  ENDDO
  DO ALPHA=1,3
   N = ALPHA+(IBMAT-1)*NBMAT
   DO JRED=1,NRED      !Jred is reduced index 
    J = INXRED(JRED)   !J is orbital index 
    DO IRED=1,NRED     !Ired is reduced index 
     I = INXRED(IRED)  !I is orbital index 
!$OMP ATOMIC
     EXCMAT(I,J,N) = EXCMAT(I,J,N) + EXCRED(IRED,JRED,ALPHA)
    ENDDO
   ENDDO
  ENDDO
 ! second part magnetic differentiation on Bmat part give a sym mat and will only contribute if BMAT is non symmetric
  IF(DOSYMPART)THEN
   DO ALPHA = 1,3
    DO J=1,NRED
     DO K=1, NBLEN
      TMP(K,J,ALPHA) = -GAORED(K,J,1)*coef2(K,IBMAT,ALPHA)
     ENDDO
    ENDDO
   ENDDO
   !  Second half-contraction of GAO's
   DO ALPHA = 1,3
      IF(XCintNoOMP)THEN
         CALL DGEMM('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
              &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS('T','N',NRED,NRED,NBLEN,1E0_realk,GAORED(:,:,1),NBLEN,&
              &     TMP(:,:,ALPHA),NBLEN,0.0E0_realk,EXCRED(:,:,ALPHA),NRED)
      ENDIF
   ENDDO
   DO ALPHA=1,3
    N = ALPHA+(IBMAT-1)*NBMAT
    DO JRED=1,NRED     !Jred is reduced index 
     J = INXRED(JRED)  !J is orbital index  
     DO IRED=1,NRED    !Ired is reduced index  
      I = INXRED(IRED) !I is orbital index 
!$OMP ATOMIC
      EXCMATSYM(I,J,N) = EXCMATSYM(I,J,N) + EXCRED(IRED,JRED,ALPHA)
     ENDDO
    ENDDO
   ENDDO
  ENDIF
 ENDDO
 call mem_dft_dealloc(EXCRED)
ENDIF
call mem_dft_dealloc(GAORED)
call mem_dft_dealloc(TMP)

END SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_LDA

!> \brief distribution routine for magnetic derivative linrsp matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF1,COEF1MAG,COEF2,COEF2MAG,NBMAT,GAOS,&
     &EXCMAT,EXCMATSYM,COORD,DFTHRI,NTYPSO,dosympart)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN,NBMAT)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN,NBMAT,3)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1MAG(NBLEN,NBMAT,3)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2MAG(NBLEN,NBMAT,3,3)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> The magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMAT(NBAST,NBAST,3*NBMAT)
!> The sym part magnetic derivative LDA Kohn-Sham matrix
REAL(REALK),intent(inout) :: EXCMATSYM(NBAST,NBAST,3*NBMAT)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
REAL(REALK),pointer :: TMP(:,:)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED,beta1,gamma1
INTEGER     :: gab1beta,gab1gamma,gab2Xbeta,gab2Xgamma,gab2Ybeta,gab2Ygamma
INTEGER     :: gab2Zbeta,gab2Zgamma,COORDINATE
INTEGER     :: IRED,JRED,alpha,beta,gamma,N,IBMAT
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
REAL(REALK),pointer :: EXCRED(:,:)
REAL(REALK) :: Rbeta,Rgamma
real(realk),parameter :: D2=2E0_realk,D05=0.5E0_realk
call mem_dft_alloc(TMP,NBLEN,NACTBAST)
call mem_dft_alloc(EXCRED,NACTBAST,NACTBAST)
!--------------------------------------------------
!   First term containing differentiation on Omega
!-------------------------------------------------
!  First half-contraction of GAO's with potential
DO IBMAT = 1,NBMAT
   DO ALPHA = 1,3
      beta = betaList(ALPHA)
      gamma= gammaList(ALPHA)
      beta1 = 4+beta
      gamma1= 4+gamma
      DO J=1,NACTBAST
         DO K=1, NBLEN
            Rgamma = COORD(gamma,K) !- Origin(gamma)       
            Rbeta = COORD(beta,K) !- Origin(beta)
            TMP(K,J) = &
                 &-coef1(K,IBMAT)*GAOS(K,J,beta1)*Rgamma &
                 &+coef1(K,IBMAT)*GAOS(K,J,gamma1)*Rbeta 
         ENDDO
      ENDDO
      !  Second half-contraction of GAO's
      IF(XCintNoOMP)THEN
         CALL DGEMM('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,GAOS,NBLEN,&
              &     TMP,NBLEN,0.0E0_realk,EXCRED,NACTBAST)
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,GAOS,NBLEN,&
           &     TMP,NBLEN,0.0E0_realk,EXCRED,NACTBAST)
      ENDIF
      N = ALPHA+(IBMAT-1)*3
      DO JRED=1,NACTBAST       !Jred is active index 
         J = INXACT(JRED)     !J is orbital index 
         DO IRED=1,NACTBAST   !Ired is active index 
            I = INXACT(IRED)  !I is orbital index 
!$OMP ATOMIC
            EXCMAT(I,J,N) = EXCMAT(I,J,N) + EXCRED(IRED,JRED)
         ENDDO
      ENDDO
   ENDDO
   !------------------------------------------------------------
   !   Second term containing differentiation on Bmat(LDA part) 
   !------------------------------------------------------------
   ! second part magnetic differentiation on Bmat part give a sym mat and will only contribute if BMAT is non symmetric
   IF(DOSYMPART)THEN
      DO ALPHA = 1,3
         DO J=1,NACTBAST
            DO K=1, NBLEN
               TMP(K,J) = -GAOS(K,J,1)*coef1MAG(K,IBMAT,ALPHA)
            ENDDO
         ENDDO
      !  Second half-contraction of GAO's
         IF(XCintNoOMP)THEN
            CALL DGEMM('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,GAOS,NBLEN,&
                 &     TMP,NBLEN,0.0E0_realk,EXCRED,NACTBAST)
         ELSE !Use Thread Safe version 
            CALL DGEMM_TS('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,GAOS,NBLEN,&
                 &     TMP,NBLEN,0.0E0_realk,EXCRED,NACTBAST)
         ENDIF
         N = ALPHA+(IBMAT-1)*NBMAT
         DO JRED=1,NACTBAST       !Jred is active index 
            J = INXACT(JRED)      !J is orbital index 
            DO IRED=1,NACTBAST    !Ired is active index 
               I = INXACT(IRED)   !I is orbital index 
!$OMP ATOMIC
               EXCMATSYM(I,J,N) = EXCMATSYM(I,J,N) + EXCRED(IRED,JRED)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   !-------------------------------------------------------
   !   Third term containing differentiation on Nabla*Omega
   !-------------------------------------------------------
   DO ALPHA = 1,3
      beta = betaList(ALPHA)
      gamma= gammaList(ALPHA)
      gab1beta = 4+beta
      gab1gamma= 4+gamma
      gab2Xbeta = 5+3*beta
      gab2Xgamma = 5+3*gamma
      gab2Ybeta = 6+3*beta
      gab2Ygamma = 6+3*gamma
      gab2Zbeta = 7+3*beta
      gab2Zgamma = 7+3*gamma
      DO J=1,NACTBAST
         DO K=1, NBLEN
            Rgamma = COORD(gamma,K) !- Origin(gamma)       
            Rbeta = COORD(beta,K) !- Origin(beta)
            TMP(K,J) = &
                 &-D05*coef2(K,IBMAT,gamma)*GAOS(K,J,gab1beta)+D05*coef2(K,IBMAT,beta)*GAOS(K,J,gab1gamma)&
                 &+D05*(GAOS(K,J,gab2Xgamma)*Rbeta-GAOS(K,J,gab2Xbeta)*Rgamma)*coef2(K,IBMAT,1)&
                 &+D05*(GAOS(K,J,gab2Ygamma)*Rbeta-GAOS(K,J,gab2Ybeta)*Rgamma)*coef2(K,IBMAT,2)&
                 &+D05*(GAOS(K,J,gab2Zgamma)*Rbeta-GAOS(K,J,gab2Zbeta)*Rgamma)*coef2(K,IBMAT,3)
         ENDDO
      ENDDO
      IF(XCintNoOMP)THEN
         CALL DGEMM('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,GAOS,NBLEN,&
              &     TMP,NBLEN,0.0E0_realk,EXCRED,NACTBAST)
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,GAOS,NBLEN,&
              &     TMP,NBLEN,0.0E0_realk,EXCRED,NACTBAST)
      ENDIF
      DO COORDINATE=1,3
         beta = betaList(ALPHA)
         gamma= gammaList(ALPHA)
         gab1beta = 4+beta
         gab1gamma= 4+gamma
         DO J=1,NACTBAST
            DO K=1, NBLEN
               Rgamma = COORD(gamma,K) !- Origin(gamma)       
               Rbeta = COORD(beta,K) !- Origin(beta)
               TMP(K,J) = -D05*(-GAOS(K,J,gab1gamma)*Rbeta+GAOS(K,J,gab1beta)*Rgamma)&
                    &*coef2(K,IBMAT,COORDINATE)
            ENDDO
         ENDDO
         IF(XCintNoOMP)THEN
            CALL DGEMM('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,&
                 & GAOS(:,:,1+COORDINATE),NBLEN,TMP,NBLEN,1.0E0_realk,&
                 & EXCRED,NACTBAST)
         ELSE !Use Thread Safe version 
            CALL DGEMM_TS('T','N',NACTBAST,NACTBAST,NBLEN,1E0_realk,&
                 & GAOS(:,:,1+COORDINATE),NBLEN,TMP,NBLEN,1.0E0_realk,&
                 & EXCRED,NACTBAST)
         ENDIF
      ENDDO
      DO JRED=1,NACTBAST      !Jred is active index 
         J = INXACT(JRED)     !J is orbital index    
         DO IRED=1,NACTBAST   !Ired is active index 
            I = INXACT(IRED)  !I is orbital index   
!$OMP ATOMIC
            EXCMAT(I,J,ALPHA) = EXCMAT(I,J,ALPHA) + EXCRED(IRED,JRED)
         ENDDO
      ENDDO
   ENDDO
   !-------------------------------------------------------
   !   Fourth term containing differentiation on Bmat(GGA)
   !-------------------------------------------------------
   !  First half-contraction of GAO's with potential
   DO ALPHA=1,3
      DO J=1,NACTBAST
         DO K=1, NBLEN
            TMP(K,J) = -coef2MAG(K,IBMAT,1,ALPHA)*GAOS(K,J,2)&
                 &   - coef2MAG(K,IBMAT,2,ALPHA)*GAOS(K,J,3)&
                 &   - coef2MAG(K,IBMAT,3,ALPHA)*GAOS(K,J,4)
         ENDDO
      ENDDO
      !  Second half-contraction of GAO's with potential
      IF(XCintNoOMP)THEN
         CALL DGEMM('T','N',NACTBAST,NACTBAST,NBLEN,1.0E0_realk,&
              &                GAOS,NBLEN,TMP,NBLEN,0.0E0_realk,&
              &                EXCRED,NACTBAST)
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS('T','N',NACTBAST,NACTBAST,NBLEN,1.0E0_realk,&
              &                GAOS,NBLEN,TMP,NBLEN,0.0E0_realk,&
              &                EXCRED,NACTBAST)
      ENDIF
      !  Distribute contributions to KS-matrix
      DO JRED=1,NACTBAST      !Jred is active index   
         J = INXACT(JRED)     !J is orbital index   
         DO IRED=1,NACTBAST   !Ired is active index   
            I = INXACT(IRED)  !I is orbital index   
!$OMP ATOMIC
            EXCMATSYM(I,J,ALPHA) = EXCMATSYM(I,J,ALPHA) + EXCRED(IRED,JRED)
         ENDDO
      ENDDO
   ENDDO
ENDDO
call mem_dft_dealloc(EXCRED)
call mem_dft_dealloc(TMP)

END SUBROUTINE II_DFT_DISTMAGDERIV_linrsp_GGA

!> \brief computes the expectation value(expval), geometrical gradient of expval and magnetic derivative of both.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_magderiv_expval_gga(LUPRI,NTYPSO,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,GAOS,COORD,EXPGRAD,BMAT,nBmat,NRED,DFTHRI,DOSYMPART)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> expectation value and geometrical and magnetic gradient
REAL(REALK),intent(inout) :: EXPGRAD(NBLEN,4,4,NBMAT)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> Number of reduced orbitals
INTEGER,intent(inout)  :: NRED
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> do the symmetric part 
LOGICAL,intent(in) :: DOSYMPART
!
REAL(REALK),pointer :: TMP(:,:,:)!,GAORED(:,:,:)
REAL(REALK) :: BMAX,BMATIJ,BMATJI,X,Y,Z,Rgamma,Rbeta
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN
INTEGER     :: IRED,JRED,IORB,JORB,gabX,gabY,gabZ,N,IBMAT,gabx1,gaby1,gabz1,M
INTEGER     :: alpha,beta,gamma,gab1beta,gab1gamma,gab2Xbeta,gab2Xgamma
INTEGER     :: gab2Ybeta,gab2Ygamma,gab2Zbeta,gab2Zgamma
integer,parameter     :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
!REAL(REALK),pointer :: BRED(:,:,:)
REAL(REALK),parameter :: zero=0E0_realk,D05=0.5E0_realk
call mem_dft_alloc(TMP,NBLEN,NACTBAST,4)
!IF(NBMAT .GT. 1)CALL LSQUIT('II_get_expval_lda does not work correctly for more than 1 nbmat',lupri)
!-------------------------------------------------------------
!  The expectation value
!-------------------------------------------------------------
   ! First half-contraction of Gaussian AO with density-matrix
 DO IBMAT=1,NBMAT
    IF(XCintNoOMP)THEN
       CALL DGEMM("N","N",NBLEN,NACTBAST,NACTBAST,1E0_realk,GAOS,NBLEN,&
            &     BMAT(:,:,IBMAT),NACTBAST,0.0E0_realk,TMP,NBLEN)
    ELSE !Use Thread Safe version 
       CALL DGEMM_TS("N","N",NBLEN,NACTBAST,NACTBAST,1E0_realk,GAOS,NBLEN,&
            &     BMAT(:,:,IBMAT),NACTBAST,0.0E0_realk,TMP,NBLEN)
    ENDIF
    ! Second half-contraction 
    DO K = 1, NBLEN
       EXPGRAD(K,1,1,IBMAT) = GAOS(K,1,1)*TMP(K,1,1)
    END DO
    DO I = 2, NACTBAST
       DO K = 1, NBLEN
          EXPGRAD(K,1,1,IBMAT) = EXPGRAD(K,1,1,IBMAT) + GAOS(K,I,1)*TMP(K,I,1)
       END DO
    END DO
 ENDDO

!-------------------------------------------------------------
!  The magnetic derivative of the expectation value
!-------------------------------------------------------------
 IF(DOSYMPART)THEN 
    !IF MATRIX B IS SYMMETRIC THIS WILL BE ZERO
    !SO WE ONLY CALCULATE THIS IF THE B MATRIX IS NONSYMMETRIC
    !RESULTING IN A SYMMETRIC FINAL CONTRIBUTION 
   gabX = 5
   gabY = 6
   gabZ = 7
   gabX1 = 2
   gabY1 = 3
   gabZ1 = 4
   DO IBMAT=1,NBMAT
   ! First half-contraction of Gaussian AO with density-matrix
    !TMP(K,J,N) = GAOS(K,I,N)*BMAT(I,J) - TMP(K,J,1) is already built
    DO N=2,4
       M=N+3
       IF(XCintNoOMP)THEN
          CALL DGEMM("N","N",NBLEN,NACTBAST,NACTBAST,1E0_realk,GAOS(:,:,M),NBLEN,&
               &     BMAT(:,:,IBMAT),NACTBAST,0.0E0_realk,TMP(:,:,N),NBLEN)
       ELSE !Use Thread Safe version 
          CALL DGEMM_TS("N","N",NBLEN,NACTBAST,NACTBAST,1E0_realk,GAOS(:,:,M),NBLEN,&
               &     BMAT(:,:,IBMAT),NACTBAST,0.0E0_realk,TMP(:,:,N),NBLEN)
       ENDIF
    ENDDO
   ! Second half-contraction of Gaussian AOs
    J=1
    DO K = 1, NBLEN
     Y = COORD(2,K) !- Origin(2)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,1,2,IBMAT) = &
          & +GAOS(K,J,1)*(Z*TMP(K,J,gabY1)-Y*TMP(K,J,gabZ1))&
          & +TMP(K,J,1)*(Y*GAOS(K,J,gabZ)-Z*GAOS(K,J,gabY))
    END DO
    DO J = 2, NACTBAST
     DO K = 1, NBLEN
        Y = COORD(2,K) !- Origin(2)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,1,2,IBMAT) = EXPGRAD(K,1,2,IBMAT)&
             & +GAOS(K,J,1)*(Z*TMP(K,J,gabY1)-Y*TMP(K,J,gabZ1))&
             & +TMP(K,J,1)*(Y*GAOS(K,J,gabZ)-Z*GAOS(K,J,gabY))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Z = COORD(3,K) !- Origin(3)      
     EXPGRAD(K,1,3,IBMAT) = &
          & +GAOS(K,J,1)*(X*TMP(K,J,gabZ1)-Z*TMP(K,J,gabX1))&
          & +TMP(K,J,1)*(Z*GAOS(K,J,gabX)-X*GAOS(K,J,gabZ))
    END DO
    DO J = 2, NACTBAST
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Z = COORD(3,K) !- Origin(3)      
        EXPGRAD(K,1,3,IBMAT) = EXPGRAD(K,1,3,IBMAT)&
             & +GAOS(K,J,1)*(X*TMP(K,J,gabZ1)-Z*TMP(K,J,gabX1))&
             & +TMP(K,J,1)*(Z*GAOS(K,J,gabX)-X*GAOS(K,J,gabZ))
     END DO
    ENDDO
    J=1
    DO K = 1, NBLEN
     X = COORD(1,K) !- Origin(1)      
     Y = COORD(2,K) !- Origin(2)      
     EXPGRAD(K,1,4,IBMAT) = &
          & +GAOS(K,J,1)*(Y*TMP(K,J,gabX1)-X*TMP(K,J,gabY1))&
          & +TMP(K,J,1)*(X*GAOS(K,J,gabY)-Y*GAOS(K,J,gabX))
    END DO
    DO J = 2, NACTBAST
     DO K = 1, NBLEN
        X = COORD(1,K) !- Origin(1)      
        Y = COORD(2,K) !- Origin(2)      
        EXPGRAD(K,1,4,IBMAT) = EXPGRAD(K,1,4,IBMAT)&
             & +GAOS(K,J,1)*(Y*TMP(K,J,gabX1)-X*TMP(K,J,gabY1))&
             & +TMP(K,J,1)*(X*GAOS(K,J,gabY)-Y*GAOS(K,J,gabX))
     END DO
    END DO
   ENDDO
 ELSE
   DO IBMAT=1,NBMAT
    DO N=1,3
     DO K = 1, NBLEN
      EXPGRAD(K,1,N,IBMAT)=zero
     ENDDO
    ENDDO
   ENDDO
 ENDIF
!-------------------------------------------------------------
!  The geometrical derivative of the expectation value
!-------------------------------------------------------------
! I should make these as subroutines also called from DISTGGA
! WARNING THIS COULD BE IMPROVED MAYBE SOME DGEMM
 DO IBMAT=1,NBMAT
    DO J = 1, NACTBAST
     BMATIJ = BMAT(1,J,IBMAT)
     BMATJI = BMAT(J,1,IBMAT)
     DO K = 1, NBLEN
      TMP(K,J,1) = GAOS(K,1,1)*(BMATIJ + BMATJI)
     END DO
     DO I = 2, NACTBAST
      BMATIJ = BMAT(I,J,IBMAT)
      BMATJI = BMAT(J,I,IBMAT)
      DO K = 1, NBLEN
       TMP(K,J,1) = TMP(K,J,1) + GAOS(K,I,1)*(BMATIJ + BMATJI)
      END DO
     ENDDO
    ENDDO
      ! Second half-contraction
    DO K = 1, NBLEN
      EXPGRAD(K,2,1,IBMAT) = GAOS(K,1,2)*TMP(K,1,1)
      EXPGRAD(K,3,1,IBMAT) = GAOS(K,1,3)*TMP(K,1,1)
      EXPGRAD(K,4,1,IBMAT) = GAOS(K,1,4)*TMP(K,1,1)
    ENDDO
    DO I = 2, NACTBAST
      DO K = 1, NBLEN
         EXPGRAD(K,2,1,IBMAT) = EXPGRAD(K,2,1,IBMAT) + GAOS(K,I,2)*TMP(K,I,1)
         EXPGRAD(K,3,1,IBMAT) = EXPGRAD(K,3,1,IBMAT) + GAOS(K,I,3)*TMP(K,I,1)
         EXPGRAD(K,4,1,IBMAT) = EXPGRAD(K,4,1,IBMAT) + GAOS(K,I,4)*TMP(K,I,1)
      ENDDO
    ENDDO
 ENDDO
!-------------------------------------------------------------
!  The mixed geometrical and magnetic derivative of the expectation value
!-------------------------------------------------------------
 do IBMAT=1,NBMAT
  do I=2,4
   do J=2,4
    do K=1,NBLEN
      EXPGRAD(K,J,I,IBMAT) = 0E0_realk
    enddo
   enddo
  enddo
 enddo
 DO IBMAT=1,NBMAT
  DO J = 1, NACTBAST
   DO I = 1, NACTBAST
    DO ALPHA =1,3
     beta = betalist(ALPHA)
     gamma = gammalist(ALPHA)
     gab1beta = 4+beta
     gab1gamma = 4+gamma
     gab2Xbeta = 5+3*beta
     gab2Xgamma = 5+3*gamma
     gab2Ybeta = 6+3*beta
     gab2Ygamma = 6+3*gamma
     gab2Zbeta = 7+3*beta
     gab2Zgamma = 7+3*gamma
     DO K = 1, NBLEN
        Rgamma = COORD(gamma,K) !- Origin(gamma)       
        Rbeta = COORD(beta,K) !- Origin(beta)
        EXPGRAD(K,gamma+1,alpha+1,IBMAT) = EXPGRAD(K,gamma+1,alpha+1,IBMAT)&
             &+GAOS(K,I,gab1beta)*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &-GAOS(K,I,1)*BMAT(I,J,IBMAT)*GAOS(K,J,gab1beta)

        EXPGRAD(K,beta+1,alpha+1,IBMAT) = EXPGRAD(K,beta+1,alpha+1,IBMAT)&
             &-GAOS(K,I,gab1gamma)*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &+GAOS(K,I,1)*BMAT(I,J,IBMAT)*GAOS(K,J,gab1gamma)

        EXPGRAD(K,2,alpha+1,IBMAT) = EXPGRAD(K,2,alpha+1,IBMAT)&
             &+GAOS(K,I,gab2Xbeta)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &-GAOS(K,I,gab2Xgamma)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &-GAOS(K,I,1)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,gab2Xbeta)&
             &+GAOS(K,I,1)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,gab2Xgamma)&
             &+GAOS(K,I,gab1beta)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,2)&
             &-GAOS(K,I,gab1gamma)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,2)&
             &-GAOS(K,I,2)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,gab1beta)&
             &+GAOS(K,I,2)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,gab1gamma)
        
        EXPGRAD(K,3,alpha+1,IBMAT) = EXPGRAD(K,3,alpha+1,IBMAT)&
             &+GAOS(K,I,gab2Ybeta)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &-GAOS(K,I,gab2Ygamma)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &-GAOS(K,I,1)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,gab2Ybeta)&
             &+GAOS(K,I,1)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,gab2Ygamma)&
             &+GAOS(K,I,gab1beta)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,3)&
             &-GAOS(K,I,gab1gamma)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,3)&
             &-GAOS(K,I,3)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,gab1beta)&
             &+GAOS(K,I,3)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,gab1gamma)
        
        EXPGRAD(K,4,alpha+1,IBMAT) = EXPGRAD(K,4,alpha+1,IBMAT)&
             &+GAOS(K,I,gab2Zbeta)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &-GAOS(K,I,gab2Zgamma)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,1)&
             &-GAOS(K,I,1)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,gab2Zbeta)&
             &+GAOS(K,I,1)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,gab2Zgamma)&
             &+GAOS(K,I,gab1beta)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,4)&
             &-GAOS(K,I,gab1gamma)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,4)&
             &-GAOS(K,I,4)*Rgamma*BMAT(I,J,IBMAT)*GAOS(K,J,gab1beta)&
             &+GAOS(K,I,4)*Rbeta*BMAT(I,J,IBMAT)*GAOS(K,J,gab1gamma)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
call mem_dft_dealloc(TMP)
NRED = NACTBAST
END SUBROUTINE II_GET_MAGDERIV_EXPVAL_GGA

!> \brief distribution routine for geometrical derivative Kohn-sham matrix LDA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF1,COEF2,GAOS,GRAD,orb2atom,&
     &BMAT,nbmat,DMAT,ndmat,natoms,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(nbast)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST,NDMAT)
!> Number of D matrices
INTEGER,intent(in)  :: NDMAT
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
REAL(REALK),pointer :: TMPB(:,:),TMPD(:,:)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST),TMP2,BT,GAOMAX
REAL(REALK),pointer :: GAORED(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB
real(realk),parameter :: D2=2E0_realk,D05=0.5E0_realk
REAL(REALK),pointer :: DRED(:,:),BRED(:,:)
call mem_dft_alloc(TMPB,NBLEN,NACTBAST)
call mem_dft_alloc(TMPD,NBLEN,NACTBAST)
call mem_dft_alloc(GAORED,NBLEN,NACTBAST,4)
GAOMAX = 0.0E0_realk
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
   GAOGMX(J) = 0.0E0_realk
   DO K = 1, NBLEN
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOMAX = MAX(GAOMAX,GAOGMX(J))
 ENDDO
ENDDO
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT. 0) THEN
! Set up reduced Gaussian AO's
   call mem_dft_alloc(DRED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !J is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !I is active index
         DRED(IRED,JRED) = DMAT(IACT,JACT,1)
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = 1, NBLEN
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   call mem_dft_alloc(BRED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(IORB,JORB,1)
      ENDDO
   ENDDO

   IF(XCintNoOMP)THEN
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED,&
           &          NBLEN,DRED,NRED,0.0E0_realk,TMPD,NBLEN)
   ELSE !Use Thread Safe version 
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED,&
           &          NBLEN,DRED,NRED,0.0E0_realk,TMPD,NBLEN)
   ENDIF
   DO JRED = 1,NRED
      IRED =1
      BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
      DO K=1, NBLEN
         TMPB(K,JRED) = GAORED(K,IRED,1)*BT
      ENDDO
      DO IRED =2,NRED
         BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
         DO K=1, NBLEN
            TMPB(K,JRED) = TMPB(K,JRED)+GAORED(K,IRED,1)*BT
         ENDDO
      ENDDO
   ENDDO
 DO IRED = 1,NRED
   iatom = atom(IRED)
   DO COORDINATE = 1,3
      DO K=1, NBLEN
         GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
              & - coef1(K)*GAORED(K,IRED,1+COORDINATE)*TMPB(K,IRED)&
              & - coef2(K)*GAORED(K,IRED,1+COORDINATE)*TMPD(K,IRED)
      ENDDO
   ENDDO
 ENDDO
 call mem_dft_dealloc(DRED)
 call mem_dft_dealloc(BRED)
ENDIF
call mem_dft_dealloc(TMPB)
call mem_dft_dealloc(TMPD)
call mem_dft_dealloc(GAORED)

END SUBROUTINE II_DFT_DISTGEODERIV_LDA

!> \brief distribution routine for geometrical derivative Kohn-sham matrix GGA closed shell worker routine
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF1,COEF2,COEF3,COEF4,GAOS,GRAD,orb2atom,&
     &BMAT,DMAT,natoms,COORD,DFTHRI,NTYPSO,GAORED,&
     &TMP,TMPX,TMPY,TMPZ)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(3,NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF3(NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF4(3,NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(NBAST)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST)
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST)
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
REAL(REALK) :: GAOMAX
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST),TMP2,BT,DT,TMPGAOMAX
!REAL(REALK),pointer :: GAORED(:,:,:)
!REAL(REALK),pointer :: TMP(:,:),TMPX(:,:)
!REAL(REALK),pointer :: TMPY(:,:),TMPZ(:,:)
REAL(REALK),intent(inout) :: GAORED(NBLEN,NACTBAST,NTYPSO)
REAL(REALK),intent(inout) :: TMP(NBLEN,NactBast),TMPX(NBLEN,NactBast)
REAL(REALK),intent(inout) :: TMPY(NBLEN,NactBast),TMPZ(NBLEN,NactBast)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB,R,RX,RY,RZ
integer,parameter     :: xlist(3)=(/5,6,7/),ylist(3)=(/6,8,9/)
integer,parameter     :: zlist(3)=(/7,9,10/)
real(realk),parameter :: D2=2E0_realk,D05=0.5E0_realk,D8=8E0_realk
REAL(REALK),pointer :: TMPRED(:,:)

! Set up maximum Gaussian AO elements
GAOMAX = 0.0E0_realk
DO JBL = 1, NBLOCKS
   DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
      TMPGAOMAX = 0.0E0_realk
      DO K = 1,NBLEN
         TMPGAOMAX = MAX(TMPGAOMAX,ABS(GAOS(K,J,1)))
      ENDDO
      GAOMAX = MAX(GAOMAX,TMPGAOMAX)
      DO N=2,NTYPSO
         DO K = 1,NBLEN
            TMPGAOMAX = MAX(TMPGAOMAX,ABS(GAOS(K,J,N)))
         ENDDO
      ENDDO
      GAOGMX(J) = TMPGAOMAX
   ENDDO
ENDDO
!!$GAOMAX = 0.0E0_realk
!!$!        Set up maximum Gaussian AO elements
!!$DO N=1,NTYPSO
!!$ DO JBL = 1, NBLOCKS
!!$  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
!!$   GAOGMX(J) = 0.0E0_realk
!!$   DO K = 1, NBLEN
!!$    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
!!$   ENDDO
!!$  ENDDO
!!$ ENDDO
!!$ENDDO
!!$DO JBL = 1, NBLOCKS
!!$ DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
!!$  GAOMAX = MAX(GAOMAX,GAOGMX(J))
!!$ ENDDO
!!$ENDDO

NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT. 0) THEN
! Set up reduced Gaussian AO's
   call mem_dft_alloc(TMPRED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !JACT is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !IACT is active index
         TMPRED(IRED,JRED) = D05*(DMAT(IACT,JACT) + DMAT(JACT,IACT))
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = 1, NBLEN
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   IF(XCintNoOMP)THEN
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,1),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMP,NBLEN)
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,2),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPX,NBLEN)
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,3),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPY,NBLEN)
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,4),&
        &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPZ,NBLEN)
   ELSE !Use Thread Safe version 
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,1),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMP,NBLEN)
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,2),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPX,NBLEN)
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,3),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPY,NBLEN)
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,4),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPZ,NBLEN)
   ENDIF
   !TMPRED contain DRED
   DO COORDINATE = 1,3
    R=1+COORDINATE
    RX=xlist(COORDINATE) 
    RY=ylist(COORDINATE) 
    RZ=zlist(COORDINATE) 
    DO IRED = 1,NRED
     iatom = atom(IRED)
     DO K=1, NBLEN
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef3(K)*GAORED(K,IRED,R)*TMP(K,IRED)
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef4(1,K)*(GAORED(K,IRED,R)*TMPX(K,IRED)+GAORED(K,IRED,RX)*TMP(K,IRED))&
        &-coef4(2,K)*(GAORED(K,IRED,R)*TMPY(K,IRED)+GAORED(K,IRED,RY)*TMP(K,IRED))&
        &-coef4(3,K)*(GAORED(K,IRED,R)*TMPZ(K,IRED)+GAORED(K,IRED,RZ)*TMP(K,IRED))
     ENDDO
    ENDDO
   ENDDO

   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         TMPRED(IRED,JRED) = D05*(BMAT(IORB,JORB) + BMAT(JORB,IORB))
      ENDDO
   ENDDO
   IF(XCintNoOMP)THEN
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,1),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMP,NBLEN)
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,2),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPX,NBLEN)
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,3),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPY,NBLEN)
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,4),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPZ,NBLEN)
   ELSE !Use Thread Safe version 
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,1),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMP,NBLEN)
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,2),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPX,NBLEN)
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,3),&
           &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPY,NBLEN)
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,4),&
        &          NBLEN,TMPRED,NRED,0.0E0_realk,TMPZ,NBLEN)
   ENDIF
   !TMPRED CONTAINS BRED
   DO COORDINATE = 1,3
    R=1+COORDINATE
    RX=xlist(COORDINATE) 
    RY=ylist(COORDINATE) 
    RZ=zlist(COORDINATE) 
    DO IRED = 1,NRED
     iatom = atom(IRED)
     DO K=1, NBLEN
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef1(K)*GAORED(K,IRED,R)*TMP(K,IRED)
      GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
        &-coef2(1,K)*(GAORED(K,IRED,RX)*TMP(K,IRED)+GAORED(K,IRED,R)*TMPX(K,IRED))&
        &-coef2(2,K)*(GAORED(K,IRED,RY)*TMP(K,IRED)+GAORED(K,IRED,R)*TMPY(K,IRED))&
        &-coef2(3,K)*(GAORED(K,IRED,RZ)*TMP(K,IRED)+GAORED(K,IRED,R)*TMPZ(K,IRED))
     ENDDO
    ENDDO
   ENDDO
   call mem_dft_dealloc(TMPRED)
ENDIF

END SUBROUTINE II_DFT_DISTGEODERIV_GGA

!> \brief distribution routine LDA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV2_LDA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF1,COEF2,GAOS,GRAD,orb2atom,&
     &BMAT,nbmat,DMAT,ndmat,natoms,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(2,NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(NBAST)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST,NDMAT)
!> Number of D matrices
INTEGER,intent(in)  :: NDMAT
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
REAL(REALK),pointer :: TMPB(:,:),TMPD(:,:),TMPA(:,:)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOGMX(NACTBAST),TMP2,BT,AT,GAOMAX
REAL(REALK),pointer :: GAORED(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB
real(realk),parameter :: D2=2E0_realk,D05=0.5E0_realk
REAL(REALK),pointer :: DRED(:,:),BRED(:,:),ARED(:,:)
call mem_dft_alloc(GAORED,NBLEN,NACTBAST,4)
call mem_dft_alloc(TMPB,NBLEN,NACTBAST)
call mem_dft_alloc(TMPD,NBLEN,NACTBAST)
call mem_dft_alloc(TMPA,NBLEN,NACTBAST)


GAOMAX = 0.0E0_realk
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
   GAOGMX(J) = 0.0E0_realk
   DO K = 1, NBLEN
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOMAX = MAX(GAOMAX,GAOGMX(J))
 ENDDO
ENDDO
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT. 0) THEN
! Set up reduced Gaussian AO's
   call mem_dft_alloc(DRED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !J is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !I is active index
         DRED(IRED,JRED) = DMAT(IACT,JACT,1)
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = 1, NBLEN
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   call mem_dft_alloc(ARED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         ARED(IRED,JRED) = BMAT(IORB,JORB,1)
      ENDDO
   ENDDO
   call mem_dft_alloc(BRED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(IORB,JORB,2)
      ENDDO
   ENDDO

   IF(XCintNoOMP)THEN
      CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED,&
           &          NBLEN,DRED,NRED,0.0E0_realk,TMPD,NBLEN)
   ELSE !Use Thread Safe version 
      CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED,&
           &          NBLEN,DRED,NRED,0.0E0_realk,TMPD,NBLEN)
   ENDIF
   DO JRED = 1,NRED
      IRED =1
      AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
      BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
      DO K=1, NBLEN
         TMPB(K,JRED) = GAORED(K,IRED,1)*BT
      ENDDO
      DO IRED =2,NRED
         BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
         DO K=1, NBLEN
            TMPB(K,JRED) = TMPB(K,JRED)+GAORED(K,IRED,1)*BT
         ENDDO
      ENDDO
      IRED =1
      DO K=1, NBLEN
         TMPA(K,JRED) = GAORED(K,IRED,1)*AT
      ENDDO
      DO IRED =2,NRED
         AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
         DO K=1, NBLEN
            TMPA(K,JRED) = TMPA(K,JRED)+GAORED(K,IRED,1)*AT
         ENDDO
      ENDDO
   ENDDO
 DO IRED = 1,NRED
   iatom = atom(IRED)
   DO COORDINATE = 1,3
      DO K=1, NBLEN
         GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
              & - coef1(1,K)*GAORED(K,IRED,1+COORDINATE)*TMPB(K,IRED)&
              & - coef1(2,K)*GAORED(K,IRED,1+COORDINATE)*TMPA(K,IRED)&
              & - coef2(K)*GAORED(K,IRED,1+COORDINATE)*TMPD(K,IRED)
      ENDDO
   ENDDO
 ENDDO
 call mem_dft_dealloc(DRED)
 call mem_dft_dealloc(BRED)
 call mem_dft_dealloc(ARED)
ENDIF
call mem_dft_dealloc(GAORED)
call mem_dft_dealloc(TMPB)
call mem_dft_dealloc(TMPD)
call mem_dft_dealloc(TMPA)

END SUBROUTINE II_DFT_DISTGEODERIV2_LDA

!> \brief distribution routine GGA geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_DISTGEODERIV2_GGA(LUPRI,NBLEN,NBLOCKS,BLOCKS,INXACT,&
     &Nactbast,NBAST,COEF1,COEF2,COEF3,GAOS,GRAD,&
     &orb2atom,BMAT,nbmat,DMAT,ndmat,natoms,COORD,DFTHRI,NTYPSO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> number of blocks
INTEGER,intent(in) :: NBLOCKS
!> contains the start and end index of active orbitals
INTEGER,intent(in) :: BLOCKS(2,NBLOCKS)
!> for a given active index INXACT provide the orbital index 
INTEGER,intent(in) :: INXACT(NACTBAST)
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF1(2,NBLEN)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF2(3,NBLEN,2)
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: COEF3(4,NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST,NTYPSO)
!> molecular gradient
REAL(REALK),intent(inout) :: GRAD(3,NATOMS)
!> for given orbital, the atom it is attached to
INTEGER,intent(in) :: orb2atom(NBAST)
!> The B matrix
REAL(REALK),intent(in) :: BMAT(NBAST,NBAST,NBMAT)
!> Number of B matrices
INTEGER,intent(in)  :: NBMAT
!> The density matrix
REAL(REALK),intent(in) :: DMAT(NactBAST,NactBAST,NDMAT)
!> Number of D matrices
INTEGER,intent(in)  :: NDMAT
!> Number of atoms
INTEGER,intent(in)  :: NATOMS
!> coordinates of the gridpoints
REAL(REALK),intent(in) :: COORD(3,NBLEN)
!> threshold on the value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> number of gaos
INTEGER,intent(in) :: NTYPSO
!
INTEGER     :: atom(Nactbast)
INTEGER     :: ISTART,IBL,I,ILEN,JBL,J,K,JSTART,JLEN,NRED
REAL(REALK) :: GAOMAX,TMP2,BT,AT,GAOGMX(NACTBAST)
REAL(REALK),pointer :: TMPB(:,:),TMPD(:,:,:)
REAL(REALK),pointer :: TMPA(:,:)
REAL(REALK),pointer :: TMPAX(:,:),TMPAY(:,:),TMPAZ(:,:)
REAL(REALK),pointer :: TMPBX(:,:),TMPBY(:,:),TMPBZ(:,:)
REAL(REALK),pointer :: GAORED(:,:,:)
INTEGER     :: INXRED(NACTBAST),IRED,JRED,alpha,beta,gamma,N,Jact,Iact,iatom
INTEGER     :: COORDINATE,IORB,JORB,R,RX,RY,RZ
real(realk),parameter :: D2=2E0_realk,D05=0.5E0_realk
REAL(REALK),pointer :: DRED(:,:),BRED(:,:),ARED(:,:)
integer,parameter     :: xlist(3)=(/5,6,7/),ylist(3)=(/6,8,9/)
integer,parameter     :: zlist(3)=(/7,9,10/)

call mem_dft_alloc(TMPB,NBLEN,NACTBAST)
call mem_dft_alloc(TMPD,NBLEN,NACTBAST,4)
call mem_dft_alloc(TMPA,NBLEN,NACTBAST)
call mem_dft_alloc(TMPAX,NBLEN,NACTBAST)
call mem_dft_alloc(TMPAY,NBLEN,NACTBAST)
call mem_dft_alloc(TMPAZ,NBLEN,NACTBAST)
call mem_dft_alloc(TMPBX,NBLEN,NACTBAST)
call mem_dft_alloc(TMPBY,NBLEN,NACTBAST)
call mem_dft_alloc(TMPBZ,NBLEN,NACTBAST)
call mem_dft_alloc(GAORED,NBLEN,NACTBAST,NTYPSO)

GAOMAX = 0.0E0_realk
!        Set up maximum Gaussian AO elements
DO N=1,NTYPSO
 DO JBL = 1, NBLOCKS
  DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
   GAOGMX(J) = 0.0E0_realk
   DO K = 1, NBLEN
    GAOGMX(J) = MAX(GAOGMX(J),ABS(GAOS(K,J,N)))
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  GAOMAX = MAX(GAOMAX,GAOGMX(J))
 ENDDO
ENDDO
NRED = 0 
DO JBL = 1, NBLOCKS
 DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
  IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
   NRED = NRED + 1
   INXRED(NRED) = J 
  ENDIF
 ENDDO
ENDDO

IF (NRED.GT. 0) THEN
! Set up reduced Gaussian AO's
   call mem_dft_alloc(DRED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JACT = INXRED(JRED)     !J is active index
      DO IRED=1,NRED          !IRED is reduced index
         IACT = INXRED(IRED)  !I is active index
         DRED(IRED,JRED) = DMAT(IACT,JACT,1)
      ENDDO
   ENDDO

   NRED = 0 
   DO JBL = 1, NBLOCKS
    DO J = BLOCKS(1,JBL), BLOCKS(2,JBL)
     IF (GAOGMX(J)*GAOMAX.GT.DFTHRI) THEN
      NRED = NRED + 1
      !WARNING WE CHANGE THE USE OF INXRED
      INXRED(NRED) = INXACT(J)
      !inxred for a given redindex the corresponding orbitalindex
      DO N=1,NTYPSO
       DO K = 1, NBLEN
        GAORED(K,NRED,N) = GAOS(K,J,N)
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DO IRED=1,NRED
      I = INXRED(IRED) !oribtal index
      atom(IRED) = orb2atom(I)
   ENDDO

   call mem_dft_alloc(ARED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         ARED(IRED,JRED) = BMAT(IORB,JORB,1)
      ENDDO
   ENDDO
   call mem_dft_alloc(BRED,NRED,NRED)
   DO JRED=1,NRED             !JRED is reduced index
      JORB = INXRED(JRED)     !J is orbital index
      DO IRED=1,NRED          !IRED is reduced index
         IORB = INXRED(IRED)  !I is orbital index
         BRED(IRED,JRED) = BMAT(IORB,JORB,2)
      ENDDO
   ENDDO
   DO I=1,4
      IF(XCintNoOMP)THEN
         CALL DGEMM('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,I),&
              &          NBLEN,DRED,NRED,0.0E0_realk,TMPD(:,:,I),NBLEN)
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS('N','N',NBLEN,NRED,NRED,1.0E0_realk,GAORED(:,:,I),&
              &          NBLEN,DRED,NRED,0.0E0_realk,TMPD(:,:,I),NBLEN)
      ENDIF
   ENDDO
   DO JRED = 1,NRED
      IRED =1
      AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
      BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
      DO K=1, NBLEN
         TMPB(K,JRED) = GAORED(K,IRED,1)*BT
         TMPBX(K,JRED) = GAORED(K,IRED,2)*BT
         TMPBY(K,JRED) = GAORED(K,IRED,3)*BT
         TMPBZ(K,JRED) = GAORED(K,IRED,4)*BT
      ENDDO
      DO IRED =2,NRED
         BT = D05*(BRED(IRED,JRED)+BRED(JRED,IRED))
         DO K=1, NBLEN
            TMPB(K,JRED) = TMPB(K,JRED)+GAORED(K,IRED,1)*BT
            TMPBX(K,JRED) = TMPBX(K,JRED)+GAORED(K,IRED,2)*BT
            TMPBY(K,JRED) = TMPBY(K,JRED)+GAORED(K,IRED,3)*BT
            TMPBZ(K,JRED) = TMPBZ(K,JRED)+GAORED(K,IRED,4)*BT
         ENDDO
      ENDDO
      IRED =1
      DO K=1, NBLEN
         TMPA(K,JRED) = GAORED(K,IRED,1)*AT
         TMPAX(K,JRED) = GAORED(K,IRED,2)*AT
         TMPAY(K,JRED) = GAORED(K,IRED,3)*AT
         TMPAZ(K,JRED) = GAORED(K,IRED,4)*AT
      ENDDO
      DO IRED =2,NRED
         AT = D05*(ARED(IRED,JRED)+ARED(JRED,IRED))
         DO K=1, NBLEN
            TMPA(K,JRED) = TMPA(K,JRED)+GAORED(K,IRED,1)*AT
            TMPAX(K,JRED) = TMPAX(K,JRED)+GAORED(K,IRED,2)*AT
            TMPAY(K,JRED) = TMPAY(K,JRED)+GAORED(K,IRED,3)*AT
            TMPAZ(K,JRED) = TMPAZ(K,JRED)+GAORED(K,IRED,4)*AT
         ENDDO
      ENDDO
   ENDDO
 DO IRED = 1,NRED
   iatom = atom(IRED)
   DO COORDINATE = 1,3
      R=1+COORDINATE
      RX=xlist(COORDINATE) 
      RY=ylist(COORDINATE) 
      RZ=zlist(COORDINATE) 
      DO K=1, NBLEN
         GRAD(COORDINATE,iatom)=GRAD(COORDINATE,iatom)&
              & - coef1(1,K)*GAORED(K,IRED,R)*TMPB(K,IRED)&
              & - coef1(2,K)*GAORED(K,IRED,R)*TMPA(K,IRED)&
              & - coef2(1,K,1)*(TMPB(K,IRED)*GAORED(K,IRED,RX)+TMPBX(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(2,K,1)*(TMPB(K,IRED)*GAORED(K,IRED,RY)+TMPBY(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(3,K,1)*(TMPB(K,IRED)*GAORED(K,IRED,RZ)+TMPBZ(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(1,K,2)*(TMPA(K,IRED)*GAORED(K,IRED,RX)+TMPAX(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(2,K,2)*(TMPA(K,IRED)*GAORED(K,IRED,RY)+TMPAY(K,IRED)*GAORED(K,IRED,R))&
              & - coef2(3,K,2)*(TMPA(K,IRED)*GAORED(K,IRED,RZ)+TMPAZ(K,IRED)*GAORED(K,IRED,R))&
              & - coef3(1,K)*GAORED(K,IRED,R)*TMPD(K,IRED,1)&
              & - coef3(2,K)*(TMPD(K,IRED,1)*GAORED(K,IRED,RX)+TMPD(K,IRED,2)*GAORED(K,IRED,R))&
              & - coef3(3,K)*(TMPD(K,IRED,1)*GAORED(K,IRED,RY)+TMPD(K,IRED,3)*GAORED(K,IRED,R))&
              & - coef3(4,K)*(TMPD(K,IRED,1)*GAORED(K,IRED,RZ)+TMPD(K,IRED,4)*GAORED(K,IRED,R))

      ENDDO
   ENDDO
 ENDDO
 call mem_dft_dealloc(DRED)
 call mem_dft_dealloc(BRED)
 call mem_dft_dealloc(ARED)
ENDIF
call mem_dft_dealloc(TMPB)
call mem_dft_dealloc(TMPD)
call mem_dft_dealloc(TMPA)
call mem_dft_dealloc(TMPAX)
call mem_dft_dealloc(TMPAY)
call mem_dft_dealloc(TMPAZ)
call mem_dft_dealloc(TMPBX)
call mem_dft_dealloc(TMPBY)
call mem_dft_dealloc(TMPBZ)
call mem_dft_dealloc(GAORED)

END SUBROUTINE II_DFT_DISTGEODERIV2_GGA

SUBROUTINE xcfun_lda_init(RHO,VXC,MXBLLEN,NBLEN,NPNT,IDMAT,NDMAT,RHOTHR,IFULL,&
     &                    XCFUNINPUT,XCFUNOUTPUT)
implicit none
!> the number of gridpoints
INTEGER,intent(in)  :: NBLEN
!> max number of gridpoints
INTEGER,intent(in)  :: MXBLLEN
!> numer of points with rho greater then threshold§
INTEGER,intent(OUT) :: NPNT
!> index of current density matrix
INTEGER,intent(IN)  :: IDMAT
!> numer of density matrices
INTEGER,intent(IN)  :: NDMAT
!> the electron density for all gridpoints
REAL(REALK),intent(in) :: RHO(MXBLLEN,NDMAT)
!> threshold on the electron density
REAL(REALK),intent(in) :: RHOTHR
REAL(REALK),intent(inout) :: VXC(NBLEN,NDMAT)
REAL(REALK),pointer       :: XCFUNINPUT(:,:),XCFUNOUTPUT(:,:)
INTEGER,pointer           :: IFULL(:)
!
INTEGER     :: IPNT
  CALL mem_dft_alloc(IFULL,NBLEN)
  NPNT = 0
  DO IPNT = 1, NBLEN
    IF(RHO(IPNT,IDMAT) .GT. RHOTHR) THEN
      NPNT = NPNT + 1
      IFULL(NPNT) = IPNT
    ELSE
      VXC(IPNT,IDMAT) = 0.0E0_realk
    ENDIF
  ENDDO
  call mem_dft_alloc(XCFUNINPUT,1,NPNT)
  call mem_dft_alloc(XCFUNOUTPUT,2,NPNT)
  DO IPNT=1,NPNT
    XCFUNINPUT(1,IPNT) = RHO(IFULL(IPNT),IDMAT)
  ENDDO
END SUBROUTINE xcfun_lda_init

SUBROUTINE xcfun_lda_free(XCFUNINPUT,XCFUNOUTPUT,IFULL)
implicit none
Real(realk),pointer :: XCFUNINPUT(:,:),XCFUNOUTPUT(:,:)
Integer,pointer     :: IFULL(:)
  call mem_dft_dealloc(IFULL)
  call mem_dft_dealloc(XCFUNINPUT)
  call mem_dft_dealloc(XCFUNOUTPUT)
END SUBROUTINE xcfun_lda_free

END MODULE IIDFTKSMWORK

