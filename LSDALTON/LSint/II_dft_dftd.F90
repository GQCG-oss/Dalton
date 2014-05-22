      MODULE IIDFTD
      use gridgenerationmodule
      use dft_memory_handling
      use Integralparameters
      !use memory_handling
      !WARNING you must not add memory_handling, all memory goes through 
      !grid_memory_handling  module so as to determine the memory used in this
      !module.
      use precision
      use TYPEDEF
      use dft_type
      use dft_typetype
      use IIDFTKSMWORK
      private
      public :: DFT_D_LSDAL_IFC 
      CONTAINS

!      INTEGER FUNCTION GET_INDX(I,J)
!      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: I,J
!      INTEGER :: IMAX,IMIN
!      IMAX=max(I,J)
!      IMIN=min(I,J)
!      GET_INDX=IMIN+IMAX*(IMAX-1)/2
!      RETURN
!      END FUNCTION GET_INDX


      SUBROUTINE GET_INDX(I,J,INDX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I,J
      INTEGER, INTENT(OUT):: INDX
      INTEGER :: IMAX,IMIN
      IMAX=max(I,J)
      IMIN=min(I,J)
      INDX=IMIN+IMAX*(IMAX-1)/2
      RETURN
      END SUBROUTINE GET_INDX





!AMT--------------------------------------------------------------------
!AMT  DFT_D Calculates the Grimme Dispersion Correction for the 
!AMT  energy and molecular gradient. 
!AMT  
!AMT  The code is tested against DFT-D3 V3.0 Rev 0 by S. Grimme
!AMT
!AMT  DFT-D2
!AMT  ------
!AMT  S. Grimme, J. Comput. Chem., 27 (2006), 1787-1799
!AMT
!AMT  DFT-D3
!AMT  ------
!AMT  S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys. 132 (2010), 154104
!AMT
!AMT  DFT-D3 + BJ Damping
!AMT  -------------------
!AMT  S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem. 32 (2011), 1456-1465
!AMT 
!AMT  The DFT-D3 empirical PARAMETERs are taken from the DFT-D3 code v3.0 rev 0.
!AMT
!AMT----------------------------------------------------------------------
      SUBROUTINE DFT_D_LSDAL_IFC(SETTING,GRAD,DIM1,DIM2,NDERIV,LUPRI)
      use ls_util
      IMPLICIT NONE
      ! external real
      INTEGER              :: DIM1,DIM2
      REAL(REALK), INTENT(INOUT) :: GRAD(DIM1,DIM2)
      INTEGER              :: LUPRI
      INTEGER              :: N_ELEM 
      INTEGER              :: N_ELEM_D2
      INTEGER              :: MX_CORD 
      INTEGER              :: NDERIV 
      INTEGER              :: IATOM,ISCOOR,NATOMS 
      REAL(REALK)          :: EDISP  
      REAL(REALK), pointer :: C6AB(:)
      REAL(REALK), pointer :: CC6AB(:)
      REAL(REALK), pointer :: R0AB(:)
      REAL(REALK), pointer :: FCOORD(:)
      REAL(REALK), pointer :: TEMPG(:)
      INTEGER, pointer     :: MXC(:)
      REAL(REALK), pointer :: R0(:)
      REAL(REALK), pointer :: C6(:)
      REAL(REALK), pointer :: RCOV(:)
      REAL(REALK), pointer :: R2R4(:)
      INTEGER, pointer     :: ICOMP(:)
      REAL(REALK), pointer :: R2AB(:)
      REAL(REALK), pointer :: DMP(:)
      REAL(REALK), pointer :: SKIPFD(:)
      REAL(REALK), pointer :: DRIJ(:)
      REAL(REALK), pointer :: DC6_SAV(:)
      REAL(REALK), pointer :: DC6IJ(:)
      REAL(REALK), pointer :: DCN(:)
      REAL(REALK), pointer :: D2(:)
      REAL(REALK), pointer :: GRADFTD(:)
      REAL(REALK), pointer :: CN(:)
      LOGICAL, pointer :: SKIPDFD(:)
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
!AMT----------------------------------------------------------------------
!AMT  Variable info:
!AMT  
!AMT  EDISP     - The dispersion energy correction (output)
!AMT  NDERIV    - The number of derivatives to calculate (presently 
!AMT              max 1 -- add hessian later?) (input)
!AMT  N_ELEM    - the element Z value upto which the D3 PARAMETERs are
!AMT              defined
!AMT  N_ELEM_D2 - the element Z value upto which the D2 PARAMETERs are
!AMT              defined 
!AMT  MX_CORD   - Maximum number of coordination number references per
!AMT              element 
!AMT
!AMT  GRADFTD   - The gradient if requested.
!AMT----------------------------------------------------------------------

      ! only for DFT 
      IF (.NOT. SETTING%DO_DFT) RETURN

      ! we do not want dispersion correction -> return
      IF (.NOT. SETTING%SCHEME%DFT%DODISP) RETURN

      NATOMS = SETTING%MOLECULE(1)%p%NATOMS
      !WRITE(LUPRI,*)'NATOMS:',NATOMS
      !WRITE(LUPRI,*)'DFT CALC?:',SETTING%DO_DFT
      !WRITE(LUPRI,*)'DODISP?:',SETTING%SCHEME%DFT%DODISP
      !WRITE(LUPRI,*)'DODISP2?:',SETTING%SCHEME%DFT%DODISP2
      !WRITE(LUPRI,*)'DODISP3?:',SETTING%SCHEME%DFT%DODISP3
      !WRITE(LUPRI,*)'CALCULATE DERIVATIVE ORD:',NDERIV
      !SETTING%EDISP
      !SETTING%SCHEME%DFT%DISPDONE = .FALSE.

      N_ELEM=94
      N_ELEM_D2=86
      MX_CORD=5

      call mem_dft_alloc(C6AB,N_ELEM*N_ELEM*MX_CORD*MX_CORD*3)
      call mem_dft_alloc(R0AB,N_ELEM*N_ELEM)
      call mem_dft_alloc(FCOORD,NATOMS*3)
      call mem_dft_alloc(TEMPG,NATOMS*3)
      call mem_dft_alloc(MXC,N_ELEM)
      call mem_dft_alloc(R0,N_ELEM_D2)
      call mem_dft_alloc(C6,N_ELEM_D2)
      call mem_dft_alloc(RCOV,N_ELEM)
      call mem_dft_alloc(R2R4,N_ELEM)
      call mem_dft_alloc(ICOMP,NATOMS*NATOMS)
      call mem_dft_alloc(R2AB,NATOMS*NATOMS)
      call mem_dft_alloc(CC6AB,NATOMS*NATOMS)
      call mem_dft_alloc(DMP,NATOMS*NATOMS)
      call mem_dft_alloc(SKIPDFD,(NATOMS*(NATOMS+1)/2))
      call mem_dft_alloc(DRIJ,(NATOMS*(NATOMS+1)/2))
      call mem_dft_alloc(DC6_SAV,(NATOMS*(NATOMS+1)/2))
      call mem_dft_alloc(DC6IJ,NATOMS*NATOMS)
      call mem_dft_alloc(DCN,(NATOMS*(NATOMS+1)/2))
      call mem_dft_alloc(D2,3)
      call mem_dft_alloc(GRADFTD,NATOMS*3)
      call mem_dft_alloc(CN,NATOMS)

      CALL DFT_D(SETTING,EDISP,NDERIV,C6AB,R0AB,FCOORD,TEMPG,&
     &           MXC,R0,C6,RCOV,R2R4,&
     &           ICOMP,R2AB,CC6AB,DMP,CN,SKIPDFD,&
     &           DRIJ,DC6_SAV,DC6IJ,DCN, D2,GRADFTD,&
     &           N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,LUPRI) 

      SETTING%EDISP = EDISP

!     Add derivative
      IF (NDERIV.EQ. 1) THEN
         IF( (DIM1.NE. 3) .OR. (DIM2.NE.NATOMS)) THEN
            CALL LSQUIT('ERROR IN DFTDISP WITH GRADIENT',-1)
         ENDIF
         DO IATOM = 1, NATOMS
            ISCOOR = (IATOM-1)*3
            GRAD(1,IATOM) = GRAD(1,IATOM) + GRADFTD(ISCOOR+1)
            GRAD(2,IATOM) = GRAD(2,IATOM) + GRADFTD(ISCOOR+2)
            GRAD(3,IATOM) = GRAD(3,IATOM) + GRADFTD(ISCOOR+3)
         END DO
        
         !WRITE(LUPRI,*)'DISPERSION CONTRIB TO MOL. GRAD.' 
         !DO IATOM = 1, NATOMS
         !  ISCOOR = (IATOM-1)*3
         !  WRITE(LUPRI,'(I6,1X,3(F20.12))') IATOM,GRADFTD(ISCOOR+1),   &
     !&                          GRADFTD(ISCOOR+2), GRADFTD(ISCOOR+3)
         !ENDDO
      END IF

      call mem_dft_dealloc(C6AB)
      call mem_dft_dealloc(R0AB)
      call mem_dft_dealloc(FCOORD)
      call mem_dft_dealloc(TEMPG)
      call mem_dft_dealloc(MXC)
      call mem_dft_dealloc(R0)
      call mem_dft_dealloc(C6)
      call mem_dft_dealloc(RCOV)
      call mem_dft_dealloc(R2R4)
      call mem_dft_dealloc(ICOMP)
      call mem_dft_dealloc(R2AB)
      call mem_dft_dealloc(CC6AB)
      call mem_dft_dealloc(DMP)
      call mem_dft_dealloc(SKIPDFD)
      call mem_dft_dealloc(DRIJ)
      call mem_dft_dealloc(DC6_SAV)
      call mem_dft_dealloc(DC6IJ)
      call mem_dft_dealloc(DCN)
      call mem_dft_dealloc(D2)
      call mem_dft_dealloc(GRADFTD)
      call mem_dft_dealloc(CN)

      RETURN
      END  SUBROUTINE DFT_D_LSDAL_IFC

      SUBROUTINE DFT_D(SETTING,EDISP,NDERIV,C6AB, R0AB, FCOORD,&
     &                 TEMPG, MXC, &
     &                 R0, C6, RCOV, R2R4, ICOMP, R2AB,&
     &                 CC6AB, DMP, CN, SKIPDFD, DRIJ, DC6_SAV, DC6IJ,&
     &                 DCN,D2, GRADFTD,N_ELEM, N_ELEM_D2, MX_CORD,NATOMS,LUPRI) 
      use ls_util
      IMPLICIT NONE
      INTEGER :: N_ELEM, N_ELEM_D2, MX_CORD, NATOMS
      REAL(REALK) :: EDISP
      REAL(REALK) :: GRADFTD(NATOMS*3)
      REAL(REALK) :: CN(NATOMS)
!C6AB has REAL(REALK)s of No. elements x No. elements x Coord No. Refs x Coord No. Refs. x 3
      REAL(REALK) :: C6AB(N_ELEM,N_ELEM,MX_CORD,MX_CORD,3)      
      REAL(REALK) :: R0AB(N_ELEM,N_ELEM)            
      REAL(REALK) :: FCOORD(NATOMS*3)       
      REAL(REALK) :: TEMPG(NATOMS*3)        
      INTEGER     :: MXC(N_ELEM)                  
      REAL(REALK) :: R0(N_ELEM_D2),C6(N_ELEM_D2)         
      REAL(REALK) :: RCOV(N_ELEM)    
      REAL(REALK) :: R2R4(N_ELEM)     
      INTEGER     :: ICOMP(NATOMS*NATOMS)  
      REAL(REALK) :: R2AB(NATOMS*NATOMS), CC6AB(NATOMS*NATOMS)
      REAL(REALK) :: DMP(NATOMS*NATOMS) 
      REAL(REALK) :: D2(3) 
      LOGICAL     :: SKIPDFD(NATOMS*(NATOMS+1)/2) 
      REAL(REALK) :: DRIJ(NATOMS*(NATOMS+1)/2)       
      REAL(REALK) :: DC6_SAV(NATOMS*(NATOMS+1)/2)    
      REAL(REALK) :: DCN(NATOMS*(NATOMS+1)/2)        
      REAL(REALK) :: DC6IJ(NATOMS,NATOMS)            
      REAL(REALK) :: S6, RS6, S18, RS18, ALP, ALP6, ALP8, E6, E8
      logical     :: L_GRAD
      INTEGER     :: LUPRI, NVER, NDERIV, IPRLOC, IOFF, IATOM
      INTEGER     :: J
!AMT External c-code routines for functional dependent PARAMETER setup
      !DFT-D2
      REAL(REALK), EXTERNAL  ::     DFT_D2_S6, DFT_D2_RS6, DFT_D2_ALP 
      !DFT-D3
      REAL(REALK), EXTERNAL  ::     DFT_D3_S6, DFT_D3_RS6, DFT_D3_ALP, DFT_D3_RS18 
      REAL(REALK), EXTERNAL  ::     DFT_D3_S18
      !DFT-D3-BJ
      REAL(REALK), EXTERNAL  ::     DFT_D3BJ_S6, DFT_D3BJ_RS6, DFT_D3BJ_ALP, DFT_D3BJ_RS18 
      REAL(REALK), EXTERNAL  ::     DFT_D3BJ_S18
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING

!AMT Parameters k1, k2 and k3
      REAL(REALK) :: k1, k2, k3
      PARAMETER (k1=16.0E0_realk)
      PARAMETER (k2=4.0E0_realk/3.0E0_realk)
      PARAMETER (k3=-4.0E0_realk)
      NVER=0
!AMT  Determine which flavour of DFT-D we are to run (2,3,3+BJ)
      IF (SETTING%SCHEME%DFT%DO_DFTD2) THEN
        NVER=2
      ELSE IF (SETTING%SCHEME%DFT%DO_DFTD3) THEN
        NVER=3
      ENDIF

      IF (NVER.EQ.0) THEN
          WRITE(LUPRI,*)'NO DFT-D VERSION SPECIFIED, DEFAULTS TO DFT-D3'
          NVER=3
          SETTING%SCHEME%DFT%DO_DFTD3=.TRUE.
      ENDIF

      CALL LSHEADER(LUPRI,'DFT-D Empirical Dispersion Correction')
      WRITE(LUPRI,'(A31,I5)')'Running DFT-D Ver.   :',NVER
      IF (SETTING%SCHEME%DFT%DO_DFTD2) THEN
       WRITE(LUPRI,'(A31,A49)')'Ref. :',&
     & 'S. Grimme, J. Comput. Chem., 27 (2006), 1787-1799'
      ELSE IF (SETTING%SCHEME%DFT%DO_DFTD3) THEN
       WRITE(LUPRI,'(A31,A47)')'Ref. :',&
     & 'S. Grimme, J. Antony, S. Ehrlich and H. Krieg, '
       WRITE(LUPRI,'(A31,A33)')' ',&
     & 'J. Chem. Phys. 132 (2010), 154104'
        IF (SETTING%SCHEME%DFT%DO_BJDAMP) THEN
        WRITE(LUPRI,'(A31,A38)')'Ref. :',&
     & 'S. Grimme, S. Ehrlich and L. Goerigk, '
        WRITE(LUPRI,'(A31,A37)')' ',&
     & 'J. Comput. Chem. 32 (2011), 1456-1465'
        ENDIF
      ENDIF

      IF (SETTING%SCHEME%DFT%DO_DFTD3) THEN
        WRITE(LUPRI,'(A31,3X,L1)')'Becke-Johnson Damping   :',&
     &                SETTING%SCHEME%DFT%DO_BJDAMP
        WRITE(LUPRI,'(A31,3X,L1)')'3-body terms   :',&
     &                SETTING%SCHEME%DFT%DO_3BODY
      ENDIF

!AMT  Check not Becke-Johnson damping / 3-body terms when running with DFT-D2
      IF ((SETTING%SCHEME%DFT%DO_DFTD2 .AND. SETTING%SCHEME%DFT%DO_BJDAMP)&
     &     .OR.(SETTING%SCHEME%DFT%DO_DFTD2 .AND.SETTING%SCHEME%DFT%DO_3BODY)) THEN
       CALL LSQUIT('ERROR: INPUT PROBLEM IN DFT_D')
      ENDIF

!AMT  Determine what functional PARAMETERs are relevant and their values
      IF (SETTING%SCHEME%DFT%DO_DFTD2) THEN
         IF (SETTING%SCHEME%DFT%L_INP_D2PAR) THEN
           WRITE(LUPRI,*)' '
           WRITE(LUPRI,'(A55)')'** User Input DFT-D2 Parameters **'
           WRITE(LUPRI,*)' '
           S6  = SETTING%SCHEME%DFT%D2_s6_inp
           RS6 = SETTING%SCHEME%DFT%D2_rs6_inp
           ALP = SETTING%SCHEME%DFT%D2_alp_inp
           S18 = 0.0E0_realk
         ELSE
           S6  = DFT_D2_S6()
           RS6 = DFT_D2_RS6()
           ALP = DFT_D2_ALP()
           S18 = 0.0E0_realk
         ENDIF
         CALL LSHEADER(LUPRI,'DFT-D2 Functional Dependent Parameters')
         WRITE(LUPRI,'(2(3X,A8,F20.12))')'s_6 = ',S6,'sr_6 = ',RS6
         WRITE(LUPRI,'(2(3X,A8,F20.12))')'s_8 = ',S18,'alp = ',alp
      ENDIF

      IF (SETTING%SCHEME%DFT%DO_DFTD3) THEN
        IF (SETTING%SCHEME%DFT%L_INP_D3PAR) THEN
         WRITE(LUPRI,*)' '
         WRITE(LUPRI,'(A55)')'** User Input DFT-D3 Parameters **'
         WRITE(LUPRI,*)' '
         S6   = SETTING%SCHEME%DFT%D3_s6_inp
         ALP  = SETTING%SCHEME%DFT%D3_alp_inp
         RS6  = SETTING%SCHEME%DFT%D3_rs6_inp
         S18  = SETTING%SCHEME%DFT%D3_s18_inp
         RS18 = SETTING%SCHEME%DFT%D3_rs18_inp
        ELSE
         IF (SETTING%SCHEME%DFT%DO_BJDAMP) THEN
           S6   = DFT_D3BJ_S6() 
           ALP  = DFT_D3BJ_ALP()
           RS6  = DFT_D3BJ_RS6()
           S18  = DFT_D3BJ_S18()
           RS18 = DFT_D3BJ_RS18()
         ELSE
           S6   = DFT_D3_S6()
           ALP  = DFT_D3_ALP()
           RS6  = DFT_D3_RS6()
           S18  = DFT_D3_S18()
           RS18 = DFT_D3_RS18()
         ENDIF
        ENDIF
        CALL LSHEADER(LUPRI,'DFT-D3 Functional Dependent Parameters')
        WRITE(LUPRI,'(2(3X,A9,F20.12))')'s_6 = ',S6,'sr_6 = ',RS6
        WRITE(LUPRI,'(2(3X,A9,F20.12))')'s_8 = ',S18,'sr_8 = ',RS18
        WRITE(LUPRI,'(3X,A9,F20.12)')'alp = ',alp
      ENDIF

!AMT  Set the run type (currently only energy / gradients)
      IF (NDERIV.EQ.0) THEN
        L_GRAD=.FALSE.
      ELSE IF (NDERIV.EQ.1) THEN
        L_GRAD=.TRUE.
        IF (L_GRAD.AND.SETTING%SCHEME%DFT%DO_3BODY) THEN
           CALL LSQUIT('ERROR: Grad. for 3-body Corrections NYI, E Only')
        ENDIF
      ELSE
        CALL LSQUIT('ERROR: DFT_D - E and G only, use 1st Ord. Opts.')
      ENDIF

!AMT  Calculate the Energy Correction

      IF (SETTING%SCHEME%DFT%DO_DFTD2) THEN
        CALL DFTD2_ENERGY(SETTING,EDISP,S6,ALP,RS6,&
     &                 N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,R0,C6,LUPRI)
      ENDIF

      IF (SETTING%SCHEME%DFT%DO_DFTD3) THEN
        ALP6 = ALP
        ALP8 = ALP+2.0E0_realk
        CALL DFTD3_ENERGY(SETTING,S6,S18,EDISP,RS6,&
     &                    RS18,ALP6,ALP8,E6,E8,C6AB,&
     &                    R0AB,MXC,ICOMP,R2AB,CC6AB,DMP,CN,D2,&
     &                    N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,RCOV,R2R4,LUPRI)
      ENDIF
 
      CALL LSHEADER(LUPRI,'Dispersion Energy Correction')
      WRITE(LUPRI,'(A31,3X,F20.12)')'E_disp: ',EDISP
      
!AMT  Calculate the Gradient Contributions
      IF (L_GRAD) THEN
        IF (SETTING%SCHEME%DFT%DO_DFTD2) THEN
          CALL DFTD2_GRAD(SETTING,S6,ALP,RS6,&
     &                    GRADFTD,N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,R0,C6,LUPRI)
        ENDIF

        IF (SETTING%SCHEME%DFT%DO_DFTD3) THEN
          IF (SETTING%SCHEME%DFT%DO_BJDAMP) THEN
            CALL DFTD3_GRAD_BJ(SETTING,S6,S18,RS6,RS18,&
     &                         ALP6,ALP8,E6,E8,C6AB,&
     &                         R0AB,MXC,SKIPDFD,&
     &                         DRIJ,DC6_SAV,DCN,DC6IJ,GRADFTD,CN,&
     &                         N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,&
     &                         RCOV,R2R4,LUPRI)   
          ELSE
            CALL DFTD3_GRAD(SETTING,S6,S18,RS6,RS18,&
     &                      ALP6,ALP8,E6,E8,C6AB,&
     &                      R0AB,MXC,SKIPDFD,&
     &                      DRIJ,DC6_SAV,DCN,DC6IJ,GRADFTD,CN,&
     &                      N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,&
     &                      RCOV,R2R4,LUPRI)
          ENDIF
        ENDIF

      ENDIF
      RETURN
      END SUBROUTINE DFT_D


!AMT ------------------------------ DFT-D2 ENERGY -------------------------------------
      SUBROUTINE DFTD2_ENERGY(SETTING,EDISP,S6,ALP,&
     &                        RS6,N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,&
     &                        R0,C6,LUPRI)
      use ls_util
      IMPLICIT NONE
      REAL(REALK) :: EDISP
      REAL(REALK) :: S6, ALP, RS6, ALP6, ALP8
      REAL(REALK) :: E, CORDAX, CORDAY, CORDAZ
      REAL(REALK) :: CORDBX, CORDBY, CORDBZ
      REAL(REALK) :: C6A, RvdWA, RvdwB, C6FAC, RX, RY, RZ, R2
      REAL(REALK) :: C6B, R6FAC, R, RR, ALPHA, EXPOA, FDMP, EADD
      REAL(REALK) :: R2_THR
      INTEGER :: N_ELEM,N_ELEM_D2,MX_CORD, NATOMS
      INTEGER :: NCENTA, NCENTB, CHARGA, CHARGB,LUPRI
      REAL(REALK) :: R0(N_ELEM_D2),C6(N_ELEM_D2)         
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
!     convert angstrom to bohr 
      REAL(REALK) :: ang_to_au, c6_conv
      PARAMETER (ang_to_au=1.88972612E0_realk)  
!     convert Joule*nm^6/mol to Bohr^6*hartree 
!        1 Bohr = 52.9177 * 10^-3 nm 
!        1 Hartree = 2.6255*10^6 Joule/mol 
      PARAMETER (c6_conv=17.3452771E0_realk) 
!     Initialize
      E = 0.0E0_realk 
!
!     Setup PARAMETERs -- for now we assume runs are for E or (E+G) so these pass
!                         on to the gradient run. i.e. never just G
      CALL SET_R0(N_ELEM_D2,R0)
      CALL SET_C6(N_ELEM_D2,C6)
!
!     Loop over pairs of atoms A-B
! 
!     Run over nuclei A 
! 
      DO NCENTA = 1, NATOMS-1 
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM_D2 .AND. CHARGA .GT. 0) THEN 
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)
            C6A    = C6(INT(CHARGA))*c6_conv 
            RvdWA  = R0(INT(CHARGA))*ang_to_au 

! 
!           Run over nuclei B 
! 
            DO NCENTB =  NCENTA+1, NATOMS 
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM_D2 .AND. CHARGB .GT. 0) THEN 
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)
                  C6B    = C6(INT(CHARGB))*c6_conv 
                  RvdWB  = R0(INT(CHARGB))*ang_to_au 
 
!                 C6 factor for atom pair A-B 
                  C6FAC = SQRT(C6A*C6B) 

!                 distance R between atoms A and B 
!                 and R^6  
                  RX    = CORDAX-CORDBX 
                  RY    = CORDAY-CORDBY 
                  RZ    = CORDAZ-CORDBZ 
                  R2    = RX*RX + RY*RY + RZ*RZ 
                  R6FAC = R2*R2*R2 
                  R     = SQRT(R2) 
 
!                 sum of the van der Waals radii RR   
                  RR   = RS6*(RvdWA + RvdWB) 
 
!                 damping function FDMP 
                  ALPHA = -ALP*R/RR+ALP 
                  EXPOA = EXP(ALPHA) 
                  FDMP  = 1.0E0_realk/(1.0E0_realk + EXPOA) 
 
!                 dispersion correction contribution 
                  EADD  = C6FAC/R6FAC * FDMP 
                  E = E + EADD 
               END IF 
            END DO   
         END IF 
      END DO 
!     final dispersion correction 
!     S6 factor is functional dependent and passed in 
      EDISP = -S6 * E 
      RETURN 
      END SUBROUTINE DFTD2_ENERGY                                       

!AMT ------------------------------ DFT-D2 GRAD -------------------------------------
      SUBROUTINE DFTD2_GRAD(SETTING,S6,ALP,RS6,&
     &                      GRADFTD,N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,&
     &                      R0,C6,LUPRI)
      use ls_util
      IMPLICIT NONE
      REAL(REALK) :: S6, ALP, RS6
      INTEGER :: N_ELEM, N_ELEM_D2, MX_CORD,NATOMS
      REAL(REALK) :: R0(N_ELEM_D2),C6(N_ELEM_D2)       
      REAL(REALK) :: GRADFTD(NATOMS*3)   
      REAL(REALK) :: R2_THR, R2_THR2, E
      REAL(REALK) :: C6A, RvdwA, C6B
      INTEGER :: NCENTA, CHARGA, ISCOOA,LUPRI
      REAL(REALK) :: CORDAX, CORDAY, CORDAZ
      REAL(REALK) :: CORDBX, CORDBY, CORDBZ
      REAL(REALK) :: RX, RY, RZ, R2, R6FAC, R, R7FAC
      INTEGER :: NCENTB, CHARGB, ISCOOB
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
      REAL(REALK) :: ang_to_au, c6_conv
      REAL(REALK) :: RvdwB, C6FAC, RR, ALPHA, EXPOA, FDMP, T1, T2, GRDX1
      REAL(REALK) :: GRDY1, GRDZ1, GRDX2, GRDY2, GRDZ2
      REAL(REALK) :: CN_THRESH
      INTEGER     :: IPRLOC
!     convert angstrom to bohr 
      PARAMETER (ang_to_au=1.88972612E0_realk)
!     convert Joule*nm^6/mol to Bohr^6*hartree 
!        1 Bohr = 52.9177 * 10^-3 nm 
!        1 Hartree = 2.6255*10^6 Joule/mol 
      PARAMETER (c6_conv=17.3452771E0_realk)

      R2_THR=9000.0E0_realk   
      R2_THR2=1600.0E0_realk

!     Initialize
      E = 0.0E0_realk
      CALL LS_DZERO(GRADFTD,3*NATOMS)
!
!     Loop over pairs of atoms A-B
! 
!     Run over nuclei A 
! 
      DO NCENTA = 1, NATOMS-1
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM_D2 .AND. CHARGA .GT. 0) THEN
            ISCOOA = (NCENTA-1)*3 
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)
            C6A    = C6(INT(CHARGA))*c6_conv
            RvdWA  = R0(INT(CHARGA))*ang_to_au

! 
!           Run over nuclei B 
! 
            DO NCENTB =  NCENTA+1, NATOMS 
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM_D2 .AND. CHARGB .GT. 0) THEN
                  ISCOOB = (NCENTB-1)*3 
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)
                  C6B    = C6(INT(CHARGB))*c6_conv
                  RvdWB  = R0(INT(CHARGB))*ang_to_au

!                 C6 factor for atom pair A-B 
                  C6FAC = SQRT(C6A*C6B)*S6

!                 distance R between atoms A and B 
!                 and R^6  
                  RX    = CORDAX-CORDBX
                  RY    = CORDAY-CORDBY
                  RZ    = CORDAZ-CORDBZ
                  R2    = RX*RX + RY*RY + RZ*RZ
                  R6FAC = R2*R2*R2
                  R     = SQRT(R2)
                  R7FAC = R6FAC*R

!                 sum of the van der Waals radii RR   
                  RR   = RS6*(RvdWA + RvdWB)

!                 damping function FDMP 
                  ALPHA = -ALP*R/RR+ALP
                  EXPOA = EXP(ALPHA) !EXPOA=DAMP6

                  FDMP  = (1.0E0_realk + EXPOA) !FDMP=DAMP1

                  T1 = EXPOA/(FDMP*FDMP*R7FAC*RR)
                  T2 = 6.0E0_realk/(FDMP*R*R7FAC)
!  Gradient Contribs from first atom in pair
                  GRDX1 = (ALP*RX*T1-T2*RX)
                  GRDY1 = (ALP*RY*T1-T2*RY)
                  GRDZ1 = (ALP*RZ*T1-T2*RZ)
!  Gradient Contribs from second atom in pair
                  GRDX2 = (ALP*(-RX)*T1+T2*RX)
                  GRDY2 = (ALP*(-RY)*T1+T2*RY)
                  GRDZ2 = (ALP*(-RZ)*T1+T2*RZ)
                  GRADFTD(ISCOOA+1) = GRADFTD(ISCOOA+1) - GRDX1*C6FAC 
                  GRADFTD(ISCOOA+2) = GRADFTD(ISCOOA+2) - GRDY1*C6FAC 
                  GRADFTD(ISCOOA+3) = GRADFTD(ISCOOA+3) - GRDZ1*C6FAC 
                  GRADFTD(ISCOOB+1) = GRADFTD(ISCOOB+1) - GRDX2*C6FAC 
                  GRADFTD(ISCOOB+2) = GRADFTD(ISCOOB+2) - GRDY2*C6FAC 
                  GRADFTD(ISCOOB+3) = GRADFTD(ISCOOB+3) - GRDZ2*C6FAC 

               END IF
            END DO  
         END IF
      END DO
      RETURN
      END SUBROUTINE DFTD2_GRAD

!AMT ------------------------------ DFT-D3 ENERGY -------------------------------------
      SUBROUTINE DFTD3_ENERGY(SETTING,S6,S8,EDISP,RS6,RS8,ALP6,&
     &                        ALP8,E6,E8,C6AB,R0AB,MXC,&
     &                        ICOMP,R2AB,CC6AB,DMP,CN,D2,&
     &                        N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,&
     &                        RCOV,R2R4,LUPRI)
      use ls_util
      IMPLICIT NONE
      REAL(REALK) :: EDISP
      REAL(REALK) :: CN(NATOMS) 
      REAL(REALK) :: S6, S8, RS6, RS8, ALP6, ALP8, E6, E8
      INTEGER :: N_ELEM, N_ELEM_D2, MX_CORD, NATOMS
      REAL(REALK) :: RCOV(N_ELEM)     
      REAL(REALK) :: R2R4(N_ELEM)     
      REAL(REALK) :: C6AB(N_ELEM,N_ELEM,MX_CORD,MX_CORD,3) 
      REAL(REALK) :: R0AB(N_ELEM,N_ELEM)       

      INTEGER     :: ICOMP(NATOMS*NATOMS)  
      REAL(REALK) :: R2AB(NATOMS*NATOMS), CC6AB(NATOMS*NATOMS)  
      REAL(REALK) :: DMP(NATOMS*NATOMS) 
      REAL(REALK) :: D2(3) 
      REAL(REALK) :: R2_THR, R2_THR2, CN_THRESH
      REAL(REALK) :: E63, bj_par1, bj_par2
      REAL(REALK) :: CORDAX, CORDAY, CORDAZ 
      REAL(REALK) :: CORDBX, CORDBY, CORDBZ 
      REAL(REALK) :: RX, Ry, RZ, R2, R, R6, R8, RAV 
      REAL(REALK) :: RR, FDMP6, FDMP8,C6,C8,SRC8C6
      REAL(REALK) :: FDMP, C9, t1, t2, t3, ang
      INTEGER     :: MXC(N_ELEM) 
      !INTEGER     :: GET_INDX
      INTEGER     :: LUPRI,IPRLOC
      INTEGER     :: NCENTA, CHARGA, NCENTB, CHARGB, NCENTC
      INTEGER     :: IJ, IK, JK
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
      REAL(REALK) :: ang_to_au, c6_conv
!     convert angstrom to bohr 
      PARAMETER (ang_to_au=1.88972612E0_realk)  
!     convert Joule*nm^6/mol to Bohr^6*hartree 
!        1 Bohr = 52.9177 * 10^-3 nm 
!        1 Hartree = 2.6255*10^6 Joule/mol 
      PARAMETER (c6_conv=17.3452771E0_realk) 

      R2_THR=9000.0E0_realk   
      R2_THR2=1600.0E0_realk
      CN_THRESH = 1600.0E0_realk

      IPRLOC=0
      IF (IPRLOC.GT.1) THEN
        CALL LSHEADER(LUPRI,'Cutoff Parameters')
        WRITE(LUPRI,'(A30,F20.12)')'Cut-off:',DSQRT(R2_THR)
        WRITE(LUPRI,'(A30,F20.12)')'Coord No. Cutoff:',DSQRT(CN_THRESH)
      ENDIF

!     Initialize
      E6 = 0.0E0_realk 
      E8 = 0.0E0_realk 
      E63 = 0.0E0_realk
!
! Setup the PARAMETERs  -- Here we assume runs are for E or E+G, 
!                          so paras are set here and passed onto
!                          the gradient run 
      CALL SET_RCOV(N_ELEM,RCOV)
      CALL SET_R2R4(N_ELEM,R2R4)
      CALL SET_R0AB(N_ELEM,ang_to_au,R0AB)
      CALL LOAD_C6AB(MX_CORD,N_ELEM,c6ab,mxc)

!     Set PARAMETERs for Becke-Johnson Damping
      bj_par1 = RS6
      bj_par2 = RS8

!     Calculate the coordinaton numbers (set in CN in CBLOCK)
      CALL SET_COORD_NOS(SETTING,RCOV,CN,N_ELEM,NATOMS)

!     Loop over pairs of atoms A-B
! 
!     Run over nuclei A 
! 
      DO NCENTA = 1, NATOMS-1 
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM .AND. CHARGA .GT. 0) THEN 
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)
! 
!           Run over nuclei B 
! 
            DO NCENTB =  NCENTA+1, NATOMS 
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM .AND. CHARGB .GT. 0) THEN 
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)

!                 distance R between atoms A and B 
!                 and R^6  
                  RX    = CORDAX-CORDBX 
                  RY    = CORDAY-CORDBY 
                  RZ    = CORDAZ-CORDBZ 
                  R2    = RX*RX + RY*RY + RZ*RZ 

                  IF (R2.LT.R2_THR) THEN
                    R      = SQRT(R2) 
                    R6     = R2**3 
                    R8     = R6*R2 
 
!                   sum of the van der Waals radii RR   
 
                    RR = R0AB(INT(CHARGB),INT(CHARGA)) / R
 
!                   damping function FDMP 
!                   c6 part
                    FDMP6  = 1.0E0_realk/(1.0E0_realk + 6.0E0_realk*(RS6*RR)**ALP6) 
!                   c8 part
                    FDMP8  = 1.0E0_realk/(1.0E0_realk + 6.0E0_realk*(RS8*RR)**ALP8) 

!                   get c6 coeffient
                    CALL SET_C6_COEFF(SETTING,MXC,C6AB,CN,NCENTA,NCENTB,C6,&
     &                                N_ELEM,MX_CORD,NATOMS)
!                   calculate c8 coefficient
                    C8 = 3.0E0_realk*C6*R2R4(INT(CHARGA))*R2R4(INT(CHARGB))

!                   if requested to add 3-body terms store the required parts
                    IF (SETTING%SCHEME%DFT%DO_3BODY) THEN
                      !IJ = GET_INDX(NCENTB,NCENTA)
                      CALL GET_INDX(NCENTB,NCENTA,IJ)
                      ICOMP(IJ) = 1
                      CC6AB(IJ) = SQRT(C6)
                      R2AB(IJ)  = R2
                      DMP(IJ) = (1.0E0_realk/RR)**(1.0E0_realk/3.0E0_realk)
                    END IF                   

!                   accumulate the E6 and E8 dispersion contributions 
                    IF (SETTING%SCHEME%DFT%DO_BJDAMP) THEN
                      SRC8C6 = SQRT(C8/C6)
                      E6  = E6 + C6 / (R6+(bj_par1*SRC8C6+bj_par2)**6) 
                      E8  = E8 + C8 / (R8+(bj_par1*SRC8C6+bj_par2)**8) 
                    ELSE
                      E6  = E6 + C6*FDMP6/R6
                      E8  = E8 + C8*FDMP8/R8 
                    ENDIF
                  END IF 
               END IF 
            END DO   
         END IF 
      END DO

!     Construct 3-body corrections and add
      IF (SETTING%SCHEME%DFT%DO_3BODY) THEN
! 
!     Run over nuclei A 
! 
      DO NCENTA = 1, NATOMS-1
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM .AND. CHARGA .GT. 0) THEN
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)
! 
!           Run over nuclei B 
! 
            DO NCENTB =  NCENTA+1, NATOMS 
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM .AND. CHARGB .GT. 0) THEN
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)

                  !IJ = GET_INDX(NCENTB,NCENTA)
                  CALL GET_INDX(NCENTB,NCENTA,IJ)

                  IF (ICOMP(IJ).EQ.1) THEN
                    DO NCENTC = NCENTB+1, NATOMS 
                      !IK=GET_INDX(NCENTC,NCENTA)
                      CALL GET_INDX(NCENTC,NCENTA,IK)
                      !JK=GET_INDX(NCENTC,NCENTB)
                      CALL GET_INDX(NCENTC,NCENTB,JK)
                      IF ((ICOMP(IK).EQ.1).AND.(ICOMP(JK).EQ.1)) THEN
                        !damping function product
                        RAV = (4.0E0_realk/3.0E0_realk)/(DMP(IK)*DMP(JK)*DMP(IJ))
                        FDMP = 1.E0_realk/( 1.E0_realk+6.E0_realk*RAV**ALP8 )
                        !triple C6 coefficient
                        C9 = CC6AB(IJ)*CC6AB(IK)*CC6AB(JK)
                        !angular terms
                        d2(1) = r2ab(ij)
                        d2(2) = r2ab(jk)
                        d2(3) = r2ab(ik)
                        t1 = (d2(1)+d2(2)-d2(3))/sqrt(d2(1)*d2(2))
                        t2 = (d2(1)+d2(3)-d2(2))/sqrt(d2(1)*d2(3))
                        t3 = (d2(3)+d2(2)-d2(1))/sqrt(d2(2)*d2(3))
                        ang=0.375E0_realk*t1*t2*t3+1.0E0_realk
                        !compute E63
                        E63 = E63 - FDMP*C9*ANG/(d2(1)*d2(2)*d2(3))**1.5E0_realk
                      ENDIF
                    ENDDO
                  END IF
               END IF
            END DO
         END IF
      END DO

      ENDIF

      E6 = -E6*S6
      E8 = -E8*S8

      E63 = E63*S6

      EDISP = E6+E8+E63
 
      RETURN
      END SUBROUTINE DFTD3_ENERGY                                       

!AMT ------------------------------ DFT-D3 GRAD -------------------------------------
      SUBROUTINE DFTD3_GRAD(SETTING,S6,S18,RS6,RS8,ALP6,&
     &                      ALP8,E6,E8,C6AB,R0AB,MXC,&
     &                      SKIPDFD,&
     &                      DRIJ,DC6_SAV,DCN,DC6IJ,GRADFTD,CN,&
     &                      N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,&
     &                      RCOV,R2R4,LUPRI)
      use ls_util
      IMPLICIT NONE
      REAL(REALK) :: S6, S18, RS6, RS8, ALP6, ALP8, E6, E8
      INTEGER     :: N_ELEM, N_ELEM_D2, MX_CORD, NATOMS
      LOGICAL     :: SKIPDFD(NATOMS*(NATOMS+1)/2) 
      REAL(REALK) :: C6AB(N_ELEM,N_ELEM,MX_CORD,MX_CORD,3)   
      REAL(REALK) :: R0AB(N_ELEM,N_ELEM)         
      REAL(REALK) :: CN(NATOMS)         
      REAL(REALK) :: RCOV(N_ELEM),R2R4(N_ELEM)    
      INTEGER     :: MXC(N_ELEM)               

      REAL(REALK) :: GRADFTD(NATOMS*3)    

      REAL(REALK) :: DRIJ(NATOMS*(NATOMS+1)/2)       
      REAL(REALK) :: DC6_SAV(NATOMS*(NATOMS+1)/2)    
      REAL(REALK) :: DCN(NATOMS*(NATOMS+1)/2)        
      REAL(REALK) :: DC6IJ(NATOMS,NATOMS)            
      REAL(REALK) :: R2_THR, R2_THR2, CN_THRESH, S8, S10, DISP
      INTEGER     :: I, NCENTA, CHARGA,CHARGB,LUPRI
      !INTEGER     :: GET_INDX
      REAL(REALK) :: CORDAX, CORDAY, CORDAZ
      REAL(REALK) :: CORDBX, CORDBY, CORDBZ
      REAL(REALK) :: RX, RY, RZ, R2, R, R6, R7, R8, R9, RR, R42, RCOVIJ
      REAL(REALK) :: C6, T6, FDMP6, T8, FDMP8, EXPT, DC6
      REAL(REALK) :: GRDX1, GRDX2, GRDY1, GRDY2, GRDZ1, GRDZ2
      INTEGER     :: NCENTB, IJ, NCENTC, IK, JK, ISCOOA, ISCOOB
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
      REAL(REALK) :: ang_to_au

      PARAMETER (ang_to_au=1.88972612E0_realk)  

!  Parameters k1, k2 and k3
      REAL(REALK) :: k1, k2, k3
      PARAMETER (k1=16.0E0_realk)
      PARAMETER (k2=4.0E0_realk/3.0E0_realk)
      PARAMETER (k3=-4.0E0_realk)

      R2_THR=9000.0E0_realk   
      R2_THR2=1600.0E0_realk
      CN_THRESH = 1600.0E0_realk

!  Initialize
      CALL LS_DZERO(GRADFTD,3*NATOMS)
      S8  = S18
      S10 = S18

      DISP=0
      DO I = 1, NATOMS*(NATOMS+1)/2
        DRIJ(I)=0.0E0_realk
        DC6_SAV(I)=0.0E0_realk
        DCN(I)=0.0E0_realk
      ENDDO

      DO I = 1, NATOMS*(NATOMS+1)/2
        SKIPDFD(I) = .TRUE.
      ENDDO
! 
!     Run over nuclei A 
! 
      DO NCENTA = 1, NATOMS 
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM .AND. CHARGA .GT. 0) THEN
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)
! 
!           Run over nuclei B 
! 
            DO NCENTB =  1, NCENTA-1
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM .AND. CHARGB .GT. 0) THEN
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)
!
!                 distance R between atoms A and B 
!                 and R^6  
                  RX    = CORDAX-CORDBX
                  RY    = CORDAY-CORDBY
                  RZ    = CORDAZ-CORDBZ
                  R2    = RX*RX + RY*RY + RZ*RZ
                  R     = DSQRT(R2)
                  R6    = R2*R2*R2
                  R7    = R6*R
                  R8    = R6*R2
                  R9    = R8*R

!                 sum of the van der Waals radii RR   

                  RR     = R0AB(INT(CHARGB),INT(CHARGA)) 
                  R42    = R2R4(INT(CHARGA))*R2R4(INT(CHARGB))
                  RCOVIJ = RCOV(INT(CHARGA)) + RCOV(INT(CHARGB)) 

                  IF (R2.LT.R2_THR) THEN
                  !IJ = GET_INDX(NCENTB,NCENTA)
                  CALL GET_INDX(NCENTB,NCENTA,IJ)
                  SKIPDFD(IJ) = .FALSE.

! Calculate dC6(i,j)/dCN(i) and dC6(i,j)/dCN(j)
! results in dC6ij
                  CALL GET_dC6_dCNij(MX_CORD,N_ELEM,C6AB,&
     &                               MXC(INT(CHARGA)),&
     &                               MXC(INT(CHARGB)), CN(NCENTA),&
     &                               CN(NCENTB), INT(CHARGA), &
     &                               INT(CHARGB), NCENTA, NCENTB,&
     &                               C6, DC6IJ(NCENTA,NCENTB),&
     &                               DC6IJ(NCENTB,NCENTA) )
              
! Calculate the damping function derivatives
                  T6 = (R/(RS6*RR))**(-ALP6)
                  FDMP6 = 1.0E0_realk / (1.0E0_realk+6.0E0_realk*T6)
                  T8 = (R/(RS8*RR))**(-ALP8)
                  FDMP8 = 1.0E0_realk / (1.0E0_realk+6.0E0_realk*T8)

                  ! d(r^(-6)) / d(r_ij) ...
                  DRIJ(IJ) = DRIJ(IJ) - S6*(6.0E0_realk/R7*C6*FDMP6)&
     &                                - S8*(24.0E0_realk/R9*C6*R42*FDMP8)
                  ! d(f_dmp)/d(r_ij) ...
                  DRIJ(IJ) = DRIJ(IJ)&
     &                     + S6*C6/R7*6.0E0_realk*ALP6*T6*FDMP6*FDMP6&
     &                     + S8*C6*R42/R9*18.0E0_realk*ALP8*T8*FDMP8*FDMP8

! Store part of the dC6 term for later
                  DC6_SAV(IJ) = S6/R6*FDMP6+3.0E0_realk*S8*R42/R8*FDMP8
                  
! Calculate dCN(i)/dr_ij = dCN(j)/dr_ij
                  IF (R2.LT.CN_THRESH) THEN
                    EXPT = exp(-k1*(RCOVIJ/R-1.0E0_realk))
                    DCN(IJ) = -k1*RCOVIJ*EXPT/&
     &                         (R*R*(EXPT+1.0E0_realk)*(EXPT+1.0E0_realk)) 

! Calculate dC6/dCN * dCN/dr_ij = dC6/dr_ij
                    DC6 = (DC6IJ(NCENTA,NCENTB)+DC6IJ(NCENTB,NCENTA))&
     &                  * DCN(IJ)
                    ! d(C6(IJ))/d(r_ij) 
                    DRIJ(IJ) = DRIJ(IJ) + DC6_SAV(IJ)*DC6
                  ELSE
                    DC6 = 0.0E0_realk
                    SKIPDFD(IJ) = .TRUE.
                  ENDIF


!  Loop to calculate dC6(i,k)/dr_ij, dC6(j,k)/dr_ij, dC6(ij)/dr_ik, dC6(ij)/dr_jk
                  DO NCENTC = 1, NCENTB-1

                    !IK = GET_INDX(NCENTA,NCENTC)
                    CALL GET_INDX(NCENTA,NCENTC,IK)
                    !JK = GET_INDX(NCENTB,NCENTC)
                    CALL GET_INDX(NCENTB,NCENTC,JK)

                    IF (.NOT. SKIPDFD(IJ)) THEN
                      DRIJ(IJ) =DRIJ(IJ) +&
     &                          DC6_SAV(IK)*DC6IJ(NCENTA,NCENTC)*DCN(IJ)
                      DRIJ(IJ) =DRIJ(IJ) +&
     &                          DC6_SAV(JK)*DC6IJ(NCENTB,NCENTC)*DCN(IJ)
                    ENDIF

                    IF (.NOT. SKIPDFD(JK)) THEN
                      DRIJ(JK) =DRIJ(JK) +&
     &                          DC6_SAV(IK)*DC6IJ(NCENTC,NCENTA)*DCN(JK)
                      DRIJ(JK) =DRIJ(JK) +&
     &                          DC6_SAV(IJ)*DC6IJ(NCENTB,NCENTA)*DCN(JK)
                    ENDIF

                    IF (.NOT. SKIPDFD(IK)) THEN
                      DRIJ(IK) =DRIJ(IK) +&
     &                          DC6_SAV(JK)*DC6IJ(NCENTC,NCENTB)*DCN(IK)
                      DRIJ(IK) =DRIJ(IK) +&
     &                          DC6_SAV(IJ)*DC6IJ(NCENTA,NCENTB)*DCN(IK)
                    ENDIF


                  ENDDO  !NCENTC

                  ENDIF

               ENDIF
            ENDDO !NCENTB
         ENDIF
      ENDDO !NCENTA


! All dE/dr_ij calculated. Now calculate dE/dr_ij * dr_ij/dxyz_i

      DO NCENTA = 2, NATOMS 
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM .AND. CHARGA .GT. 0) THEN
            ISCOOA = (NCENTA-1)*3
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)

            DO NCENTB = 1, NCENTA-1
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM .AND. CHARGB .GT. 0) THEN
                  ISCOOB = (NCENTB-1)*3 
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)

!                 distance R between atoms A and B 
!                 and R^6  
                  RX    = CORDAX-CORDBX
                  RY    = CORDAY-CORDBY
                  RZ    = CORDAZ-CORDBZ
                  R2    = RX*RX + RY*RY + RZ*RZ
                  R     = DSQRT(R2)

                  !IJ = GET_INDX(NCENTA,NCENTB)
                  CALL GET_INDX(NCENTA,NCENTB,IJ)
                  GRDX1 = -DRIJ(IJ)*RX/R
                  GRDY1 = -DRIJ(IJ)*RY/R
                  GRDZ1 = -DRIJ(IJ)*RZ/R
                  GRDX2 =  DRIJ(IJ)*RX/R
                  GRDY2 =  DRIJ(IJ)*RY/R
                  GRDZ2 =  DRIJ(IJ)*RZ/R
                  GRADFTD(ISCOOA+1) = GRADFTD(ISCOOA+1) + GRDX1
                  GRADFTD(ISCOOA+2) = GRADFTD(ISCOOA+2) + GRDY1
                  GRADFTD(ISCOOA+3) = GRADFTD(ISCOOA+3) + GRDZ1
                  GRADFTD(ISCOOB+1) = GRADFTD(ISCOOB+1) + GRDX2
                  GRADFTD(ISCOOB+2) = GRADFTD(ISCOOB+2) + GRDY2
                  GRADFTD(ISCOOB+3) = GRADFTD(ISCOOB+3) + GRDZ2
               ENDIF          
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DFTD3_GRAD


!AMT ------------------------------ DFT-D3 GRAD with Becke-Johnson damping-------------------------------------
      SUBROUTINE DFTD3_GRAD_BJ(SETTING,S6,S18,RS6,RS8,ALP6,&
     &                         ALP8,E6,E8,C6AB,R0AB,MXC,&
     &                         SKIPDFD,&
     &                         DRIJ,DC6_SAV,DCN,DC6IJ,GRADFTD,CN,&
     &                         N_ELEM,N_ELEM_D2,MX_CORD,NATOMS,&
     &                         RCOV,R2R4,LUPRI)
      use ls_util
      IMPLICIT NONE
      REAL(REALK) :: S6, S18, RS6, RS8, ALP6, ALP8, E6, E8
      INTEGER     :: N_ELEM, N_ELEM_D2, MX_CORD, NATOMS
      LOGICAL     :: SKIPDFD(NATOMS*(NATOMS+1)/2)    
      REAL(REALK) :: C6AB(N_ELEM,N_ELEM,MX_CORD,MX_CORD,3)             
      REAL(REALK) :: R0AB(N_ELEM,N_ELEM)                   
      REAL(REALK) :: RCOV(N_ELEM),R2R4(N_ELEM)             
      INTEGER     :: MXC(N_ELEM)                         

      REAL(REALK) :: GRADFTD(NATOMS*3)             
      REAL(REALK) :: CN(NATOMS)             

      REAL(REALK) :: DRIJ(NATOMS*(NATOMS+1)/2)      
      REAL(REALK) :: DC6_SAV(NATOMS*(NATOMS+1)/2)   
      REAL(REALK) :: DCN(NATOMS*(NATOMS+1)/2)       
      REAL(REALK) :: DC6IJ(NATOMS,NATOMS)           
      REAL(REALK) :: R2_THR, R2_THR2, CN_THRESH, A1, A2, S8, DISP
      INTEGER     :: I, NCENTA, CHARGA, NCENTB, CHARGB
      !INTEGER     :: GET_INDX
      INTEGER     :: J, IJ, NCENTC, IK, JK, ISCOOA, ISCOOB, LUPRI
      REAL(REALK) :: CORDAX,CORDAY,CORDAZ
      REAL(REALK) :: CORDBX,CORDBY,CORDBZ
      REAL(REALK) :: RX, RY, RZ, R2, R, R4, R6, R7, R8, R9, RR, R42
      REAL(REALK) :: RCOVIJ, R0BJ, T6, T8, EXPT, C6, DC6
      REAL(REALK) :: GRDX1, GRDY1, GRDZ1, GRDX2, GRDY2, GRDZ2
      REAL(REALK) :: CN_I
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
      REAL(REALK) :: ang_to_au
      PARAMETER (ang_to_au=1.88972612E0_realk)

! Parameters k1, k2 and k3
      REAL(REALK) :: k1, k2, k3
      PARAMETER (k1=16.0E0_realk)
      PARAMETER (k2=4.0E0_realk/3.0E0_realk)
      PARAMETER (k3=-4.0E0_realk)

      R2_THR=9000.0E0_realk   
      R2_THR2=1600.0E0_realk
      CN_THRESH = 1600.0E0_realk

! Initialize
      CALL LS_DZERO(GRADFTD,3*NATOMS)
      A1 = RS6
      A2 = RS8
      S8  = S18

      DISP=0
      DO I = 1, NATOMS*(NATOMS+1)/2
        DRIJ(I)=0.0E0_realk
        DC6_SAV(I)=0.0E0_realk
        DCN(I)=0.0E0_realk
      ENDDO

      DO I = 1, NATOMS*(NATOMS+1)/2
        SKIPDFD(I) = .TRUE.
      ENDDO
!
!     Run over nuclei A 
! 
      DO NCENTA = 1, NATOMS 
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM .AND. CHARGA .GT. 0) THEN
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)
! 
!           Run over nuclei B 
! 
            DO NCENTB =  1, NCENTA-1
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM .AND. CHARGB .GT. 0) THEN
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)

!                 distance R between atoms A and B 
!                 and R^6  
                  RX    = CORDAX-CORDBX
                  RY    = CORDAY-CORDBY
                  RZ    = CORDAZ-CORDBZ
                  R2    = RX*RX + RY*RY + RZ*RZ
                  R     = DSQRT(R2)
                  R4    = R2*R2
                  R6    = R2*R2*R2
                  R7    = R6*R
                  R8    = R6*R2
                  R9    = R8*R

!                 sum of the van der Waals radii RR   

                  RR     = R0AB(INT(CHARGB),INT(CHARGA))
                  R42    = R2R4(INT(CHARGA))*R2R4(INT(CHARGB))
                  RCOVIJ = RCOV(INT(CHARGA)) + RCOV(INT(CHARGB))

                  IF (R2.LT.R2_THR) THEN
                  !IJ = GET_INDX(NCENTB,NCENTA)
                  CALL GET_INDX(NCENTB,NCENTA,IJ)
                  SKIPDFD(IJ) = .FALSE.

! Calculate dC6(i,j)/dCN(i) and dC6(i,j)/dCN(j)
! results in dC6ij
                  CALL GET_dC6_dCNij(MX_CORD,N_ELEM,C6AB,&
     &                               MXC(INT(CHARGA)),&
     &                               MXC(INT(CHARGB)), CN(NCENTA),&
     &                               CN(NCENTB), INT(CHARGA),&
     &                               INT(CHARGB), NCENTA, NCENTB,&
     &                               C6, DC6IJ(NCENTA,NCENTB),&
     &                               DC6IJ(NCENTB,NCENTA) )
! Calculate the damping functions and their derivatives
! Use the Becke--Johnson radius
                  R0BJ = A1*SQRT(3.0E0_realk*R42)+A2   

                  T6 = (R6+R0BJ**6)
                  T8 = (R8+R0BJ**8)

                  DRIJ(IJ) = DRIJ(IJ)&
     &                     - S6*C6*6.0E0_realk*R4*R/(T6*T6)&
     &                     - S8*C6*24.0E0_realk*R42*R6*R/(T8*T8)

! Store part of the dC6 term for later
                  DC6_SAV(IJ) = S6/T6+3.0E0_realk*S8*R42/T8

! Calculate dCN(i)/dr_ij = dCN(j)/dr_ij
                  IF (R2.LT.CN_THRESH) THEN
                    EXPT = exp(-k1*(RCOVIJ/R-1.0E0_realk))
                    DCN(IJ) = -k1*RCOVIJ*EXPT/&
     &                         (R*R*(EXPT+1.0E0_realk)*(EXPT+1.0E0_realk))

! Calculate dC6/dCN * dCN/dr_ij = dC6/dr_ij
                    DC6 = (DC6IJ(NCENTA,NCENTB)+DC6IJ(NCENTB,NCENTA))&
     &                  * DCN(IJ)
                    ! d(C6(IJ))/d(r_ij) 
                    DRIJ(IJ) = DRIJ(IJ) + DC6_SAV(IJ)*DC6
                  ELSE
                    DC6 = 0.0E0_realk
                    DCN(IJ) = 0.0E0_realk
                    SKIPDFD(IJ) = .TRUE.
                  ENDIF

! Loop to calculate dC6(i,k)/dr_ij, dC6(j,k)/dr_ij, dC6(ij)/dr_ik, dC6(ij)/dr_jk
                  DO NCENTC = 1, NCENTB-1

                    !IK = GET_INDX(NCENTA,NCENTC)
                    CALL GET_INDX(NCENTA,NCENTC,IK)
                    !JK = GET_INDX(NCENTB,NCENTC)
                    CALL GET_INDX(NCENTB,NCENTC,JK)


                    DC6=DC6IJ(NCENTA,NCENTC)*DCN(IJ)
                    DRIJ(IJ) = DRIJ(IJ) + DC6_SAV(IK)*DC6

                    DC6 = DC6IJ(NCENTB,NCENTC)*DCN(IJ)
                    DRIJ(IJ) = DRIJ(IJ) + DC6_SAV(JK)*DC6

                    IF (.NOT. SKIPDFD(JK)) THEN
                      DRIJ(JK) =DRIJ(JK) +&
     &                          DC6_SAV(IK)*DC6IJ(NCENTC,NCENTA)*DCN(JK)
                      DRIJ(JK) =DRIJ(JK) +&
     &                          DC6_SAV(IJ)*DC6IJ(NCENTB,NCENTA)*DCN(JK)
                    ENDIF

                    IF (.NOT. SKIPDFD(IK)) THEN
                      DRIJ(IK) =DRIJ(IK) +&
     &                          DC6_SAV(JK)*DC6IJ(NCENTC,NCENTB)*DCN(IK)
                      DRIJ(IK) =DRIJ(IK) +&
     &                          DC6_SAV(IJ)*DC6IJ(NCENTA,NCENTB)*DCN(IK)
                    ENDIF


                  ENDDO  !NCENTC

                  ENDIF

               ENDIF
            ENDDO !NCENTB
         ENDIF
      ENDDO !NCENTA

! All dE/dr_ij calculated. Now calculate dE/dr_ij * dr_ij/dxyz_i

      DO NCENTA = 2, NATOMS 
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE
         IF (CHARGA .LE. N_ELEM .AND. CHARGA .GT. 0) THEN
            ISCOOA = (NCENTA-1)*3
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CENTER(3)

            DO NCENTB = 1, NCENTA-1
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE
               IF (CHARGB .LE. N_ELEM .AND. CHARGB .GT. 0) THEN
                  ISCOOB = (NCENTB-1)*3
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CENTER(3)

!                 distance R between atoms A and B 
!                 and R^6  
                  RX    = CORDAX-CORDBX
                  RY    = CORDAY-CORDBY
                  RZ    = CORDAZ-CORDBZ
                  R2    = RX*RX + RY*RY + RZ*RZ
                  R     = DSQRT(R2)

                  !IJ = GET_INDX(NCENTA,NCENTB)
                  CALL GET_INDX(NCENTA,NCENTB,IJ)
                  GRDX1 = -DRIJ(IJ)*RX/R
                  GRDY1 = -DRIJ(IJ)*RY/R
                  GRDZ1 = -DRIJ(IJ)*RZ/R
                  GRDX2 =  DRIJ(IJ)*RX/R
                  GRDY2 =  DRIJ(IJ)*RY/R
                  GRDZ2 =  DRIJ(IJ)*RZ/R
                  GRADFTD(ISCOOA+1) = GRADFTD(ISCOOA+1) + GRDX1
                  GRADFTD(ISCOOA+2) = GRADFTD(ISCOOA+2) + GRDY1
                  GRADFTD(ISCOOA+3) = GRADFTD(ISCOOA+3) + GRDZ1
                  GRADFTD(ISCOOB+1) = GRADFTD(ISCOOB+1) + GRDX2
                  GRADFTD(ISCOOB+2) = GRADFTD(ISCOOB+2) + GRDY2
                  GRADFTD(ISCOOB+3) = GRADFTD(ISCOOB+3) + GRDZ2
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DFTD3_GRAD_BJ

!---------------------- DFT-D3 SET COORDINATION NUMBERS --------------------------------

      SUBROUTINE SET_COORD_NOS(SETTING,RCOV,CN,N_ELEM,NATOMS)
      use ls_util
      IMPLICIT NONE
      REAL(REALK) :: RCOV(N_ELEM)
      REAL(REALK) :: CN(NATOMS)
      REAL(REALK) :: CN_THRESH, CN_I, RX, RY, RZ, R2, R, RRCOV
      REAL(REALK) :: RR, CNDAMP
      INTEGER     :: IZ, JZ
      INTEGER     :: I, J, N_ELEM,NATOMS
! Parameters k1, k2 and k3
      REAL(REALK) :: k1, k2, k3
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING

      PARAMETER (k1=16.0E0_realk)
      PARAMETER (k2=4.0E0_realk/3.0E0_realk)
      PARAMETER (k3=-4.0E0_realk)
      CN_THRESH = 1600.0E0_realk

      DO I = 1,NATOMS
        CN_I=0.0E0_realk
        DO J = 1, NATOMS
          IF (I.NE.J) THEN
            RX = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(1) - SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)
            RY = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(2) - SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)
            RZ = SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(3) - SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)
            R2    = RX*RX + RY*RY + RZ*RZ
            IF (R2.LT.CN_THRESH) THEN
! Covalent separation in au
              R  = SQRT(R2)
              IZ=INT(SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE)
              JZ=INT(SETTING%MOLECULE(1)%p%ATOM(J)%CHARGE)
              RRCOV = RCOV(IZ) + RCOV(JZ)
              RR = RRCOV/R
! Coord-No Damping Function
              CNDAMP = 1.0E0_realk / (1.0E0_realk+exp(-k1*(RR-1.0E0_realk)))
              CN_I = CN_I + CNDAMP
            ENDIF
          ENDIF
        ENDDO
        CN(I)=CN_I
      ENDDO
      RETURN
      END SUBROUTINE SET_COORD_NOS 


!-------------------- DFT-D3 Set C6 PARAMETERs -----------------------------------------
      SUBROUTINE SET_C6_COEFF(SETTING,MXC,C6AB,CN,NCENTA,NCENTB,C6,&
     &                        N_ELEM,MX_CORD,NATOMS)
      use ls_util
      IMPLICIT NONE
      REAL(REALK) :: C6AB(N_ELEM,N_ELEM,MX_CORD,MX_CORD,3)
      REAL(REALK) :: C6, C6MEM, RSUM, CSUM, R_SAVE, CN_I, CN_J
      REAL(REALK) :: CN1, CN2, R, TMP1
      REAL(REALK) :: CN(NATOMS)
      INTEGER     :: NCENTA, NCENTB, N_ELEM, MX_CORD
      INTEGER     :: I, J, NATOMS
      INTEGER MXC(N_ELEM)
      !  external types
      TYPE(LSSETTING),   INTENT(INOUT) :: SETTING

! Parameters k1, k2 and k3
      REAL(REALK) :: k1, k2, k3
      PARAMETER (k1=16.0E0_realk)
      PARAMETER (k2=4.0E0_realk/3.0E0_realk)
      PARAMETER (k3=-4.0E0_realk)
      REAL(REALK) :: BIGN
      PARAMETER (BIGN = 1.0D99)
      INTEGER     :: CHRGA, CHRGB

      C6MEM  = -BIGN
      RSUM   = 0.0E0_realk
      CSUM   = 0.0E0_realk
      C6     = 0.0E0_realk
      R_SAVE = BIGN

      CHRGA=INT(SETTING%MOLECULE(1)%p%ATOM(NCENTA)%CHARGE)
      CHRGB=INT(SETTING%MOLECULE(1)%p%ATOM(NCENTB)%CHARGE)

! Coord No.s in common block 
      CN_I=CN(NCENTA)
      CN_J=CN(NCENTB)

      DO I = 1, MXC(CHRGA)            !Loop over coord no. refs of atom 1
        DO J = 1, MXC(CHRGB)          !Loop over coord no. refs of atom 2
          C6 = C6AB(CHRGA,CHRGB,I,J,1)
          IF (C6.GT.0) THEN
            CN1 = C6AB(CHRGA,CHRGB,I,J,2) 
            CN2 = C6AB(CHRGA,CHRGB,I,J,3)
            !distance
            R = (CN1-CN_I)**2 + (CN2-CN_J)**2
            IF (R.LT.R_SAVE) THEN
               R_SAVE = R
               C6MEM  = C6
            ENDIF
            TMP1 = exp(k3*R)
            RSUM = RSUM + TMP1
            CSUM = CSUM + TMP1*C6
          ENDIF
        ENDDO
      ENDDO

      IF (RSUM.GT.1.0D-99) THEN
         C6 = CSUM/RSUM
      ELSE
         C6 = C6MEM
      END IF

      RETURN
      END SUBROUTINE SET_C6_COEFF

      SUBROUTINE LOAD_C6AB(MX_CORD,N_ELEM,C6AB,MXCi)
      use ls_util
      IMPLICIT NONE
      INTEGER   :: MX_CORD,N_ELEM,MXCi(N_ELEM),NLINES
      REAL*8    :: C6AB(N_ELEM,N_ELEM,MX_CORD,MX_CORD,3)
      INTEGER   :: IAT,JAT,L,J,K,IL,IADR,JADR,NN,KK
#include "lsdalton_dftd_pars.h"
      !nlines set in the above file
      C6AB=-1
      MXCi=0
      KK=1
      DO NN=1,NLINES
       IAT=INT(PARS(KK+1))
       JAT=INT(PARS(KK+2))
       CALL LIMIT(IAT,JAT,IADR,JADR)
       MXCi(IAT)=MAX(MXCi(IAT),IADR)
       MXCi(JAT)=MAX(MXCi(JAT),JADR)

       C6AB(IAT,JAT,IADR,JADR,1)=PARS(KK)
       C6AB(IAT,JAT,IADR,JADR,2)=PARS(KK+3)
       C6AB(IAT,JAT,IADR,JADR,3)=PARS(KK+4)

       C6AB(JAT,IAT,JADR,IADR,1)=PARS(KK)
       C6AB(JAT,IAT,JADR,IADR,2)=PARS(KK+4)
       C6AB(JAT,IAT,JADR,IADR,3)=PARS(KK+3)
       KK=(NN*5)+1
      ENDDO
      RETURN
      END SUBROUTINE LOAD_C6AB  


      SUBROUTINE LIMIT(IAT,JAT,IADR,JADR)
      IMPLICIT NONE
      INTEGER  :: IAT,JAT,IADR,JADR,I
      IADR=1
      JADR=1
      I=100
 10   IF(IAT.GT.100) THEN
         IAT=IAT-100
         IADR=IADR+1
         GOTO 10
      ENDIf

      I=100
 20   IF(JAT.GT.100) THEN
         JAT=JAT-100
         JADR=JADR+1
         GOTO 20
      ENDIf
      RETURN
      END SUBROUTINE LIMIT


      SUBROUTINE GET_dC6_dCNij(MX_CORD,N_ELEM,C6AB,MXCI,MXCJ,CNI,CNJ,&
     &           IZI,IZJ,IAT,JAT,C6CHECK,DC6I,DC6J)
      IMPLICIT NONE
! Parameters k1, k2 and k3
      REAL(REALK) :: K1, K2, K3
      PARAMETER (K1=16.0E0_realk)
      PARAMETER (K2=4.0E0_realk/3.0E0_realk)
      PARAMETER (K3=-4.0E0_realk)
      INTEGER     :: MX_CORD,N_ELEM
      REAL(REALK) :: C6AB(N_ELEM,N_ELEM,MX_CORD,MX_CORD,3)
      INTEGER     :: MXCI,MXCJ   !mxc(iz(iat))
      REAL(REALK) :: CNI,CNJ
      INTEGER     :: IAT,JAT,IZI,IZJ
      REAL(REALK) :: DC6I,DC6J,C6CHECK

      INTEGER I,J,A,B
      REAL(REALK) :: NUMER,DENOM,DNUMER_I,DDENOM_I,DNUMER_J,DDENOM_J
      REAL(REALK) :: EXPTERM,CN_REFI,CN_REFJ,C6REF,R
      REAL(REALK) :: C6MEM,R_SAVE

      C6MEM=-1.d99
      R_SAVE=9999.0
      NUMER=0.0E0_realk
      DENOM=0.0E0_realk

      DNUMER_I=0.E0_realk
      DDENOM_I=0.E0_realk
      DNUMER_J=0.E0_realk
      DDENOM_J=0.E0_realk


      DO A=1,MXCI
        DO B=1,MXCJ
          C6REF=C6AB(IZI,IZJ,A,B,1)
          IF (C6REF.GT.0) THEN
            CN_REFI=C6AB(IZI,IZJ,A,B,2)
            CN_REFJ=C6AB(IZI,IZJ,A,B,3)
            R=(CN_REFI-CNI)*(CN_REFI-CNI)+(CN_REFJ-CNJ)*(CN_REFJ-CNJ)
            IF (R.LT.R_SAVE) THEN
               R_SAVE=R
               C6MEM=C6REF
            ENDIF
            EXPTERM=EXP(K3*R)
            NUMER=NUMER+C6REF*EXPTERM
            DENOM=DENOM+EXPTERM
            DNUMER_I=DNUMER_I+C6REF*EXPTERM*&
     &             2.E0_realk*K3*(CNI-CN_REFI)
            DDENOM_I=DDENOM_I+EXPTERM*&
     &             2.E0_realk*K3*(CNI-CN_REFI)

            DNUMER_J=DNUMER_J+C6REF*EXPTERM*&
     &             2.E0_realk*K3*(CNJ-CN_REFJ)
            DDENOM_J=DDENOM_J+EXPTERM*&
     &             2.E0_realk*K3*(CNJ-CN_REFJ)
          ENDIF
        ENDDO !b
      ENDDO !a

      IF (DENOM.GT.1.0D-99) THEN
        C6CHECK=NUMER/DENOM
        DC6I=((DNUMER_I*DENOM)-(DDENOM_I*NUMER))&
     &    /(DENOM*DENOM)
        DC6J=((DNUMER_J*DENOM)-(DDENOM_J*NUMER))&
     &    /(DENOM*DENOM)
      ELSE
        C6CHECK=C6MEM
        DC6I=0.0E0_realk
        DC6J=0.0E0_realk
      ENDIF
      END SUBROUTINE GET_dC6_dCNij


      SUBROUTINE SET_R0AB(N_ELEM,ANG_TO_AU,R)
      use ls_util
      IMPLICIT NONE
      INTEGER      :: N_ELEM,I,J,K
      REAL(REALK)  :: R(N_ELEM,N_ELEM),ANG_TO_AU
      REAL(REALK)  :: R0AB(4465)
      R0AB(   1:  70)=(/&
     &   2.1823E0_realk,  1.8547E0_realk,  1.7347E0_realk,  2.9086E0_realk,  2.5732E0_realk,  3.4956E0_realk,  2.3550E0_realk&
     &,  2.5095E0_realk,  2.9802E0_realk,  3.0982E0_realk,  2.5141E0_realk,  2.3917E0_realk,  2.9977E0_realk,  2.9484E0_realk&
     &,  3.2160E0_realk,  2.4492E0_realk,  2.2527E0_realk,  3.1933E0_realk,  3.0214E0_realk,  2.9531E0_realk,  2.9103E0_realk&
     &,  2.3667E0_realk,  2.1328E0_realk,  2.8784E0_realk,  2.7660E0_realk,  2.7776E0_realk,  2.7063E0_realk,  2.6225E0_realk&
     &,  2.1768E0_realk,  2.0625E0_realk,  2.6395E0_realk,  2.6648E0_realk,  2.6482E0_realk,  2.5697E0_realk,  2.4846E0_realk&
     &,  2.4817E0_realk,  2.0646E0_realk,  1.9891E0_realk,  2.5086E0_realk,  2.6908E0_realk,  2.6233E0_realk,  2.4770E0_realk&
     &,  2.3885E0_realk,  2.3511E0_realk,  2.2996E0_realk,  1.9892E0_realk,  1.9251E0_realk,  2.4190E0_realk,  2.5473E0_realk&
     &,  2.4994E0_realk,  2.4091E0_realk,  2.3176E0_realk,  2.2571E0_realk,  2.1946E0_realk,  2.1374E0_realk,  2.9898E0_realk&
     &,  2.6397E0_realk,  3.6031E0_realk,  3.1219E0_realk,  3.7620E0_realk,  3.2485E0_realk,  2.9357E0_realk,  2.7093E0_realk&
     &,  2.5781E0_realk,  2.4839E0_realk,  3.7082E0_realk,  2.5129E0_realk,  2.7321E0_realk,  3.1052E0_realk,  3.2962E0_realk&
     &/)
      R0AB(  71: 140)=(/&
     &   3.1331E0_realk,  3.2000E0_realk,  2.9586E0_realk,  3.0822E0_realk,  2.8582E0_realk,  2.7120E0_realk,  3.2570E0_realk&
     &,  3.4839E0_realk,  2.8766E0_realk,  2.7427E0_realk,  3.2776E0_realk,  3.2363E0_realk,  3.5929E0_realk,  3.2826E0_realk&
     &,  3.0911E0_realk,  2.9369E0_realk,  2.9030E0_realk,  2.7789E0_realk,  3.3921E0_realk,  3.3970E0_realk,  4.0106E0_realk&
     &,  2.8884E0_realk,  2.6605E0_realk,  3.7513E0_realk,  3.1613E0_realk,  3.3605E0_realk,  3.3325E0_realk,  3.0991E0_realk&
     &,  2.9297E0_realk,  2.8674E0_realk,  2.7571E0_realk,  3.8129E0_realk,  3.3266E0_realk,  3.7105E0_realk,  3.7917E0_realk&
     &,  2.8304E0_realk,  2.5538E0_realk,  3.3932E0_realk,  3.1193E0_realk,  3.1866E0_realk,  3.1245E0_realk,  3.0465E0_realk&
     &,  2.8727E0_realk,  2.7664E0_realk,  2.6926E0_realk,  3.4608E0_realk,  3.2984E0_realk,  3.5142E0_realk,  3.5418E0_realk&
     &,  3.5017E0_realk,  2.6190E0_realk,  2.4797E0_realk,  3.1331E0_realk,  3.0540E0_realk,  3.0651E0_realk,  2.9879E0_realk&
     &,  2.9054E0_realk,  2.8805E0_realk,  2.7330E0_realk,  2.6331E0_realk,  3.2096E0_realk,  3.5668E0_realk,  3.3684E0_realk&
     &,  3.3686E0_realk,  3.3180E0_realk,  3.3107E0_realk,  2.4757E0_realk,  2.4019E0_realk,  2.9789E0_realk,  3.1468E0_realk&
     &/)
      R0AB( 141: 210)=(/&
     &   2.9768E0_realk,  2.8848E0_realk,  2.7952E0_realk,  2.7457E0_realk,  2.6881E0_realk,  2.5728E0_realk,  3.0574E0_realk&
     &,  3.3264E0_realk,  3.3562E0_realk,  3.2529E0_realk,  3.1916E0_realk,  3.1523E0_realk,  3.1046E0_realk,  2.3725E0_realk&
     &,  2.3289E0_realk,  2.8760E0_realk,  2.9804E0_realk,  2.9093E0_realk,  2.8040E0_realk,  2.7071E0_realk,  2.6386E0_realk&
     &,  2.5720E0_realk,  2.5139E0_realk,  2.9517E0_realk,  3.1606E0_realk,  3.2085E0_realk,  3.1692E0_realk,  3.0982E0_realk&
     &,  3.0352E0_realk,  2.9730E0_realk,  2.9148E0_realk,  3.2147E0_realk,  2.8315E0_realk,  3.8724E0_realk,  3.4621E0_realk&
     &,  3.8823E0_realk,  3.3760E0_realk,  3.0746E0_realk,  2.8817E0_realk,  2.7552E0_realk,  2.6605E0_realk,  3.9740E0_realk&
     &,  3.6192E0_realk,  3.6569E0_realk,  3.9586E0_realk,  3.6188E0_realk,  3.3917E0_realk,  3.2479E0_realk,  3.1434E0_realk&
     &,  4.2411E0_realk,  2.7597E0_realk,  3.0588E0_realk,  3.3474E0_realk,  3.6214E0_realk,  3.4353E0_realk,  3.4729E0_realk&
     &,  3.2487E0_realk,  3.3200E0_realk,  3.0914E0_realk,  2.9403E0_realk,  3.4972E0_realk,  3.7993E0_realk,  3.6773E0_realk&
     &,  3.8678E0_realk,  3.5808E0_realk,  3.8243E0_realk,  3.5826E0_realk,  3.4156E0_realk,  3.8765E0_realk,  4.1035E0_realk&
     &/)
      R0AB( 211: 280)=(/&
     &   2.7361E0_realk,  2.9765E0_realk,  3.2475E0_realk,  3.5004E0_realk,  3.4185E0_realk,  3.4378E0_realk,  3.2084E0_realk&
     &,  3.2787E0_realk,  3.0604E0_realk,  2.9187E0_realk,  3.4037E0_realk,  3.6759E0_realk,  3.6586E0_realk,  3.8327E0_realk&
     &,  3.5372E0_realk,  3.7665E0_realk,  3.5310E0_realk,  3.3700E0_realk,  3.7788E0_realk,  3.9804E0_realk,  3.8903E0_realk&
     &,  2.6832E0_realk,  2.9060E0_realk,  3.2613E0_realk,  3.4359E0_realk,  3.3538E0_realk,  3.3860E0_realk,  3.1550E0_realk&
     &,  3.2300E0_realk,  3.0133E0_realk,  2.8736E0_realk,  3.4024E0_realk,  3.6142E0_realk,  3.5979E0_realk,  3.5295E0_realk&
     &,  3.4834E0_realk,  3.7140E0_realk,  3.4782E0_realk,  3.3170E0_realk,  3.7434E0_realk,  3.9623E0_realk,  3.8181E0_realk&
     &,  3.7642E0_realk,  2.6379E0_realk,  2.8494E0_realk,  3.1840E0_realk,  3.4225E0_realk,  3.2771E0_realk,  3.3401E0_realk&
     &,  3.1072E0_realk,  3.1885E0_realk,  2.9714E0_realk,  2.8319E0_realk,  3.3315E0_realk,  3.5979E0_realk,  3.5256E0_realk&
     &,  3.4980E0_realk,  3.4376E0_realk,  3.6714E0_realk,  3.4346E0_realk,  3.2723E0_realk,  3.6859E0_realk,  3.8985E0_realk&
     &,  3.7918E0_realk,  3.7372E0_realk,  3.7211E0_realk,  2.9230E0_realk,  2.6223E0_realk,  3.4161E0_realk,  2.8999E0_realk&
     &/)
      R0AB( 281: 350)=(/&
     &   3.0557E0_realk,  3.3308E0_realk,  3.0555E0_realk,  2.8508E0_realk,  2.7385E0_realk,  2.6640E0_realk,  3.5263E0_realk&
     &,  3.0277E0_realk,  3.2990E0_realk,  3.7721E0_realk,  3.5017E0_realk,  3.2751E0_realk,  3.1368E0_realk,  3.0435E0_realk&
     &,  3.7873E0_realk,  3.2858E0_realk,  3.2140E0_realk,  3.1727E0_realk,  3.2178E0_realk,  3.4414E0_realk,  2.5490E0_realk&
     &,  2.7623E0_realk,  3.0991E0_realk,  3.3252E0_realk,  3.1836E0_realk,  3.2428E0_realk,  3.0259E0_realk,  3.1225E0_realk&
     &,  2.9032E0_realk,  2.7621E0_realk,  3.2490E0_realk,  3.5110E0_realk,  3.4429E0_realk,  3.3845E0_realk,  3.3574E0_realk&
     &,  3.6045E0_realk,  3.3658E0_realk,  3.2013E0_realk,  3.6110E0_realk,  3.8241E0_realk,  3.7090E0_realk,  3.6496E0_realk&
     &,  3.6333E0_realk,  3.0896E0_realk,  3.5462E0_realk,  2.4926E0_realk,  2.7136E0_realk,  3.0693E0_realk,  3.2699E0_realk&
     &,  3.1272E0_realk,  3.1893E0_realk,  2.9658E0_realk,  3.0972E0_realk,  2.8778E0_realk,  2.7358E0_realk,  3.2206E0_realk&
     &,  3.4566E0_realk,  3.3896E0_realk,  3.3257E0_realk,  3.2946E0_realk,  3.5693E0_realk,  3.3312E0_realk,  3.1670E0_realk&
     &,  3.5805E0_realk,  3.7711E0_realk,  3.6536E0_realk,  3.5927E0_realk,  3.5775E0_realk,  3.0411E0_realk,  3.4885E0_realk&
     &/)
      R0AB( 351: 420)=(/&
     &   3.4421E0_realk,  2.4667E0_realk,  2.6709E0_realk,  3.0575E0_realk,  3.2357E0_realk,  3.0908E0_realk,  3.1537E0_realk&
     &,  2.9235E0_realk,  3.0669E0_realk,  2.8476E0_realk,  2.7054E0_realk,  3.2064E0_realk,  3.4519E0_realk,  3.3593E0_realk&
     &,  3.2921E0_realk,  3.2577E0_realk,  3.2161E0_realk,  3.2982E0_realk,  3.1339E0_realk,  3.5606E0_realk,  3.7582E0_realk&
     &,  3.6432E0_realk,  3.5833E0_realk,  3.5691E0_realk,  3.0161E0_realk,  3.4812E0_realk,  3.4339E0_realk,  3.4327E0_realk&
     &,  2.4515E0_realk,  2.6338E0_realk,  3.0511E0_realk,  3.2229E0_realk,  3.0630E0_realk,  3.1265E0_realk,  2.8909E0_realk&
     &,  3.0253E0_realk,  2.8184E0_realk,  2.6764E0_realk,  3.1968E0_realk,  3.4114E0_realk,  3.3492E0_realk,  3.2691E0_realk&
     &,  3.2320E0_realk,  3.1786E0_realk,  3.2680E0_realk,  3.1036E0_realk,  3.5453E0_realk,  3.7259E0_realk,  3.6090E0_realk&
     &,  3.5473E0_realk,  3.5327E0_realk,  3.0018E0_realk,  3.4413E0_realk,  3.3907E0_realk,  3.3593E0_realk,  3.3462E0_realk&
     &,  2.4413E0_realk,  2.6006E0_realk,  3.0540E0_realk,  3.1987E0_realk,  3.0490E0_realk,  3.1058E0_realk,  2.8643E0_realk&
     &,  2.9948E0_realk,  2.7908E0_realk,  2.6491E0_realk,  3.1950E0_realk,  3.3922E0_realk,  3.3316E0_realk,  3.2585E0_realk&
     &/)
      R0AB( 421: 490)=(/&
     &   3.2136E0_realk,  3.1516E0_realk,  3.2364E0_realk,  3.0752E0_realk,  3.5368E0_realk,  3.7117E0_realk,  3.5941E0_realk&
     &,  3.5313E0_realk,  3.5164E0_realk,  2.9962E0_realk,  3.4225E0_realk,  3.3699E0_realk,  3.3370E0_realk,  3.3234E0_realk&
     &,  3.3008E0_realk,  2.4318E0_realk,  2.5729E0_realk,  3.0416E0_realk,  3.1639E0_realk,  3.0196E0_realk,  3.0843E0_realk&
     &,  2.8413E0_realk,  2.7436E0_realk,  2.7608E0_realk,  2.6271E0_realk,  3.1811E0_realk,  3.3591E0_realk,  3.3045E0_realk&
     &,  3.2349E0_realk,  3.1942E0_realk,  3.1291E0_realk,  3.2111E0_realk,  3.0534E0_realk,  3.5189E0_realk,  3.6809E0_realk&
     &,  3.5635E0_realk,  3.5001E0_realk,  3.4854E0_realk,  2.9857E0_realk,  3.3897E0_realk,  3.3363E0_realk,  3.3027E0_realk&
     &,  3.2890E0_realk,  3.2655E0_realk,  3.2309E0_realk,  2.8502E0_realk,  2.6934E0_realk,  3.2467E0_realk,  3.1921E0_realk&
     &,  3.5663E0_realk,  3.2541E0_realk,  3.0571E0_realk,  2.9048E0_realk,  2.8657E0_realk,  2.7438E0_realk,  3.3547E0_realk&
     &,  3.3510E0_realk,  3.9837E0_realk,  3.6871E0_realk,  3.4862E0_realk,  3.3389E0_realk,  3.2413E0_realk,  3.1708E0_realk&
     &,  3.6096E0_realk,  3.6280E0_realk,  3.6860E0_realk,  3.5568E0_realk,  3.4836E0_realk,  3.2868E0_realk,  3.3994E0_realk&
     &/)
      R0AB( 491: 560)=(/&
     &   3.3476E0_realk,  3.3170E0_realk,  3.2950E0_realk,  3.2874E0_realk,  3.2606E0_realk,  3.9579E0_realk,  2.9226E0_realk&
     &,  2.6838E0_realk,  3.7867E0_realk,  3.1732E0_realk,  3.3872E0_realk,  3.3643E0_realk,  3.1267E0_realk,  2.9541E0_realk&
     &,  2.8505E0_realk,  2.7781E0_realk,  3.8475E0_realk,  3.3336E0_realk,  3.7359E0_realk,  3.8266E0_realk,  3.5733E0_realk&
     &,  3.3959E0_realk,  3.2775E0_realk,  3.1915E0_realk,  3.9878E0_realk,  3.8816E0_realk,  3.5810E0_realk,  3.5364E0_realk&
     &,  3.5060E0_realk,  3.8097E0_realk,  3.3925E0_realk,  3.3348E0_realk,  3.3019E0_realk,  3.2796E0_realk,  3.2662E0_realk&
     &,  3.2464E0_realk,  3.7136E0_realk,  3.8619E0_realk,  2.9140E0_realk,  2.6271E0_realk,  3.4771E0_realk,  3.1774E0_realk&
     &,  3.2560E0_realk,  3.1970E0_realk,  3.1207E0_realk,  2.9406E0_realk,  2.8322E0_realk,  2.7571E0_realk,  3.5455E0_realk&
     &,  3.3514E0_realk,  3.5837E0_realk,  3.6177E0_realk,  3.5816E0_realk,  3.3902E0_realk,  3.2604E0_realk,  3.1652E0_realk&
     &,  3.7037E0_realk,  3.6283E0_realk,  3.5858E0_realk,  3.5330E0_realk,  3.4884E0_realk,  3.5789E0_realk,  3.4094E0_realk&
     &,  3.3473E0_realk,  3.3118E0_realk,  3.2876E0_realk,  3.2707E0_realk,  3.2521E0_realk,  3.5570E0_realk,  3.6496E0_realk&
     &/)
      R0AB( 561: 630)=(/&
     &   3.6625E0_realk,  2.7300E0_realk,  2.5870E0_realk,  3.2471E0_realk,  3.1487E0_realk,  3.1667E0_realk,  3.0914E0_realk&
     &,  3.0107E0_realk,  2.9812E0_realk,  2.8300E0_realk,  2.7284E0_realk,  3.3259E0_realk,  3.3182E0_realk,  3.4707E0_realk&
     &,  3.4748E0_realk,  3.4279E0_realk,  3.4182E0_realk,  3.2547E0_realk,  3.1353E0_realk,  3.5116E0_realk,  3.9432E0_realk&
     &,  3.8828E0_realk,  3.8303E0_realk,  3.7880E0_realk,  3.3760E0_realk,  3.7218E0_realk,  3.3408E0_realk,  3.3059E0_realk&
     &,  3.2698E0_realk,  3.2446E0_realk,  3.2229E0_realk,  3.4422E0_realk,  3.5023E0_realk,  3.5009E0_realk,  3.5268E0_realk&
     &,  2.6026E0_realk,  2.5355E0_realk,  3.1129E0_realk,  3.2863E0_realk,  3.1029E0_realk,  3.0108E0_realk,  2.9227E0_realk&
     &,  2.8694E0_realk,  2.8109E0_realk,  2.6929E0_realk,  3.1958E0_realk,  3.4670E0_realk,  3.4018E0_realk,  3.3805E0_realk&
     &,  3.3218E0_realk,  3.2815E0_realk,  3.2346E0_realk,  3.0994E0_realk,  3.3937E0_realk,  3.7266E0_realk,  3.6697E0_realk&
     &,  3.6164E0_realk,  3.5730E0_realk,  3.2522E0_realk,  3.5051E0_realk,  3.4686E0_realk,  3.4355E0_realk,  3.4084E0_realk&
     &,  3.3748E0_realk,  3.3496E0_realk,  3.3692E0_realk,  3.4052E0_realk,  3.3910E0_realk,  3.3849E0_realk,  3.3662E0_realk&
     &/)
      R0AB( 631: 700)=(/&
     &   2.5087E0_realk,  2.4814E0_realk,  3.0239E0_realk,  3.1312E0_realk,  3.0535E0_realk,  2.9457E0_realk,  2.8496E0_realk&
     &,  2.7780E0_realk,  2.7828E0_realk,  2.6532E0_realk,  3.1063E0_realk,  3.3143E0_realk,  3.3549E0_realk,  3.3120E0_realk&
     &,  3.2421E0_realk,  3.1787E0_realk,  3.1176E0_realk,  3.0613E0_realk,  3.3082E0_realk,  3.5755E0_realk,  3.5222E0_realk&
     &,  3.4678E0_realk,  3.4231E0_realk,  3.1684E0_realk,  3.3528E0_realk,  3.3162E0_realk,  3.2827E0_realk,  3.2527E0_realk&
     &,  3.2308E0_realk,  3.2029E0_realk,  3.3173E0_realk,  3.3343E0_realk,  3.3092E0_realk,  3.2795E0_realk,  3.2452E0_realk&
     &,  3.2096E0_realk,  3.2893E0_realk,  2.8991E0_realk,  4.0388E0_realk,  3.6100E0_realk,  3.9388E0_realk,  3.4475E0_realk&
     &,  3.1590E0_realk,  2.9812E0_realk,  2.8586E0_realk,  2.7683E0_realk,  4.1428E0_realk,  3.7911E0_realk,  3.8225E0_realk&
     &,  4.0372E0_realk,  3.7059E0_realk,  3.4935E0_realk,  3.3529E0_realk,  3.2492E0_realk,  4.4352E0_realk,  4.0826E0_realk&
     &,  3.9733E0_realk,  3.9254E0_realk,  3.8646E0_realk,  3.9315E0_realk,  3.7837E0_realk,  3.7465E0_realk,  3.7211E0_realk&
     &,  3.7012E0_realk,  3.6893E0_realk,  3.6676E0_realk,  3.7736E0_realk,  4.0660E0_realk,  3.7926E0_realk,  3.6158E0_realk&
     &/)
      R0AB( 701: 770)=(/&
     &   3.5017E0_realk,  3.4166E0_realk,  4.6176E0_realk,  2.8786E0_realk,  3.1658E0_realk,  3.5823E0_realk,  3.7689E0_realk&
     &,  3.5762E0_realk,  3.5789E0_realk,  3.3552E0_realk,  3.4004E0_realk,  3.1722E0_realk,  3.0212E0_realk,  3.7241E0_realk&
     &,  3.9604E0_realk,  3.8500E0_realk,  3.9844E0_realk,  3.7035E0_realk,  3.9161E0_realk,  3.6751E0_realk,  3.5075E0_realk&
     &,  4.1151E0_realk,  4.2877E0_realk,  4.1579E0_realk,  4.1247E0_realk,  4.0617E0_realk,  3.4874E0_realk,  3.9848E0_realk&
     &,  3.9280E0_realk,  3.9079E0_realk,  3.8751E0_realk,  3.8604E0_realk,  3.8277E0_realk,  3.8002E0_realk,  3.9981E0_realk&
     &,  3.7544E0_realk,  4.0371E0_realk,  3.8225E0_realk,  3.6718E0_realk,  4.3092E0_realk,  4.4764E0_realk,  2.8997E0_realk&
     &,  3.0953E0_realk,  3.4524E0_realk,  3.6107E0_realk,  3.6062E0_realk,  3.5783E0_realk,  3.3463E0_realk,  3.3855E0_realk&
     &,  3.1746E0_realk,  3.0381E0_realk,  3.6019E0_realk,  3.7938E0_realk,  3.8697E0_realk,  3.9781E0_realk,  3.6877E0_realk&
     &,  3.8736E0_realk,  3.6451E0_realk,  3.4890E0_realk,  3.9858E0_realk,  4.1179E0_realk,  4.0430E0_realk,  3.9563E0_realk&
     &,  3.9182E0_realk,  3.4002E0_realk,  3.8310E0_realk,  3.7716E0_realk,  3.7543E0_realk,  3.7203E0_realk,  3.7053E0_realk&
     &/)
      R0AB( 771: 840)=(/&
     &   3.6742E0_realk,  3.8318E0_realk,  3.7631E0_realk,  3.7392E0_realk,  3.9892E0_realk,  3.7832E0_realk,  3.6406E0_realk&
     &,  4.1701E0_realk,  4.3016E0_realk,  4.2196E0_realk,  2.8535E0_realk,  3.0167E0_realk,  3.3978E0_realk,  3.5363E0_realk&
     &,  3.5393E0_realk,  3.5301E0_realk,  3.2960E0_realk,  3.3352E0_realk,  3.1287E0_realk,  2.9967E0_realk,  3.6659E0_realk&
     &,  3.7239E0_realk,  3.8070E0_realk,  3.7165E0_realk,  3.6368E0_realk,  3.8162E0_realk,  3.5885E0_realk,  3.4336E0_realk&
     &,  3.9829E0_realk,  4.0529E0_realk,  3.9584E0_realk,  3.9025E0_realk,  3.8607E0_realk,  3.3673E0_realk,  3.7658E0_realk&
     &,  3.7035E0_realk,  3.6866E0_realk,  3.6504E0_realk,  3.6339E0_realk,  3.6024E0_realk,  3.7708E0_realk,  3.7283E0_realk&
     &,  3.6896E0_realk,  3.9315E0_realk,  3.7250E0_realk,  3.5819E0_realk,  4.1457E0_realk,  4.2280E0_realk,  4.1130E0_realk&
     &,  4.0597E0_realk,  3.0905E0_realk,  2.7998E0_realk,  3.6448E0_realk,  3.0739E0_realk,  3.2996E0_realk,  3.5262E0_realk&
     &,  3.2559E0_realk,  3.0518E0_realk,  2.9394E0_realk,  2.8658E0_realk,  3.7514E0_realk,  3.2295E0_realk,  3.5643E0_realk&
     &,  3.7808E0_realk,  3.6931E0_realk,  3.4723E0_realk,  3.3357E0_realk,  3.2429E0_realk,  4.0280E0_realk,  3.5589E0_realk&
     &/)
      R0AB( 841: 910)=(/&
     &   3.4636E0_realk,  3.4994E0_realk,  3.4309E0_realk,  3.6177E0_realk,  3.2946E0_realk,  3.2376E0_realk,  3.2050E0_realk&
     &,  3.1847E0_realk,  3.1715E0_realk,  3.1599E0_realk,  3.5555E0_realk,  3.8111E0_realk,  3.7693E0_realk,  3.5718E0_realk&
     &,  3.4498E0_realk,  3.3662E0_realk,  4.1608E0_realk,  3.7417E0_realk,  3.6536E0_realk,  3.6154E0_realk,  3.8596E0_realk&
     &,  3.0301E0_realk,  2.7312E0_realk,  3.5821E0_realk,  3.0473E0_realk,  3.2137E0_realk,  3.4679E0_realk,  3.1975E0_realk&
     &,  2.9969E0_realk,  2.8847E0_realk,  2.8110E0_realk,  3.6931E0_realk,  3.2076E0_realk,  3.4943E0_realk,  3.5956E0_realk&
     &,  3.6379E0_realk,  3.4190E0_realk,  3.2808E0_realk,  3.1860E0_realk,  3.9850E0_realk,  3.5105E0_realk,  3.4330E0_realk&
     &,  3.3797E0_realk,  3.4155E0_realk,  3.6033E0_realk,  3.2737E0_realk,  3.2145E0_realk,  3.1807E0_realk,  3.1596E0_realk&
     &,  3.1461E0_realk,  3.1337E0_realk,  3.4812E0_realk,  3.6251E0_realk,  3.7152E0_realk,  3.5201E0_realk,  3.3966E0_realk&
     &,  3.3107E0_realk,  4.1128E0_realk,  3.6899E0_realk,  3.6082E0_realk,  3.5604E0_realk,  3.7834E0_realk,  3.7543E0_realk&
     &,  2.9189E0_realk,  2.6777E0_realk,  3.4925E0_realk,  2.9648E0_realk,  3.1216E0_realk,  3.2940E0_realk,  3.0975E0_realk&
     &/)
      R0AB( 911: 980)=(/&
     &   2.9757E0_realk,  2.8493E0_realk,  2.7638E0_realk,  3.6085E0_realk,  3.1214E0_realk,  3.4006E0_realk,  3.4793E0_realk&
     &,  3.5147E0_realk,  3.3806E0_realk,  3.2356E0_realk,  3.1335E0_realk,  3.9144E0_realk,  3.4183E0_realk,  3.3369E0_realk&
     &,  3.2803E0_realk,  3.2679E0_realk,  3.4871E0_realk,  3.1714E0_realk,  3.1521E0_realk,  3.1101E0_realk,  3.0843E0_realk&
     &,  3.0670E0_realk,  3.0539E0_realk,  3.3890E0_realk,  3.5086E0_realk,  3.5895E0_realk,  3.4783E0_realk,  3.3484E0_realk&
     &,  3.2559E0_realk,  4.0422E0_realk,  3.5967E0_realk,  3.5113E0_realk,  3.4576E0_realk,  3.6594E0_realk,  3.6313E0_realk&
     &,  3.5690E0_realk,  2.8578E0_realk,  2.6334E0_realk,  3.4673E0_realk,  2.9245E0_realk,  3.0732E0_realk,  3.2435E0_realk&
     &,  3.0338E0_realk,  2.9462E0_realk,  2.8143E0_realk,  2.7240E0_realk,  3.5832E0_realk,  3.0789E0_realk,  3.3617E0_realk&
     &,  3.4246E0_realk,  3.4505E0_realk,  3.3443E0_realk,  3.1964E0_realk,  3.0913E0_realk,  3.8921E0_realk,  3.3713E0_realk&
     &,  3.2873E0_realk,  3.2281E0_realk,  3.2165E0_realk,  3.4386E0_realk,  3.1164E0_realk,  3.1220E0_realk,  3.0761E0_realk&
     &,  3.0480E0_realk,  3.0295E0_realk,  3.0155E0_realk,  3.3495E0_realk,  3.4543E0_realk,  3.5260E0_realk,  3.4413E0_realk&
     &/)
      R0AB( 981:1050)=(/&
     &   3.3085E0_realk,  3.2134E0_realk,  4.0170E0_realk,  3.5464E0_realk,  3.4587E0_realk,  3.4006E0_realk,  3.6027E0_realk&
     &,  3.5730E0_realk,  3.4945E0_realk,  3.4623E0_realk,  2.8240E0_realk,  2.5960E0_realk,  3.4635E0_realk,  2.9032E0_realk&
     &,  3.0431E0_realk,  3.2115E0_realk,  2.9892E0_realk,  2.9148E0_realk,  2.7801E0_realk,  2.6873E0_realk,  3.5776E0_realk&
     &,  3.0568E0_realk,  3.3433E0_realk,  3.3949E0_realk,  3.4132E0_realk,  3.3116E0_realk,  3.1616E0_realk,  3.0548E0_realk&
     &,  3.8859E0_realk,  3.3719E0_realk,  3.2917E0_realk,  3.2345E0_realk,  3.2274E0_realk,  3.4171E0_realk,  3.1293E0_realk&
     &,  3.0567E0_realk,  3.0565E0_realk,  3.0274E0_realk,  3.0087E0_realk,  2.9939E0_realk,  3.3293E0_realk,  3.4249E0_realk&
     &,  3.4902E0_realk,  3.4091E0_realk,  3.2744E0_realk,  3.1776E0_realk,  4.0078E0_realk,  3.5374E0_realk,  3.4537E0_realk&
     &,  3.3956E0_realk,  3.5747E0_realk,  3.5430E0_realk,  3.4522E0_realk,  3.4160E0_realk,  3.3975E0_realk,  2.8004E0_realk&
     &,  2.5621E0_realk,  3.4617E0_realk,  2.9154E0_realk,  3.0203E0_realk,  3.1875E0_realk,  2.9548E0_realk,  2.8038E0_realk&
     &,  2.7472E0_realk,  2.6530E0_realk,  3.5736E0_realk,  3.0584E0_realk,  3.3304E0_realk,  3.3748E0_realk,  3.3871E0_realk&
     &/)
      R0AB(1051:1120)=(/&
     &   3.2028E0_realk,  3.1296E0_realk,  3.0214E0_realk,  3.8796E0_realk,  3.3337E0_realk,  3.2492E0_realk,  3.1883E0_realk&
     &,  3.1802E0_realk,  3.4050E0_realk,  3.0756E0_realk,  3.0478E0_realk,  3.0322E0_realk,  3.0323E0_realk,  3.0163E0_realk&
     &,  3.0019E0_realk,  3.3145E0_realk,  3.4050E0_realk,  3.4656E0_realk,  3.3021E0_realk,  3.2433E0_realk,  3.1453E0_realk&
     &,  3.9991E0_realk,  3.5017E0_realk,  3.4141E0_realk,  3.3520E0_realk,  3.5583E0_realk,  3.5251E0_realk,  3.4243E0_realk&
     &,  3.3851E0_realk,  3.3662E0_realk,  3.3525E0_realk,  2.7846E0_realk,  2.5324E0_realk,  3.4652E0_realk,  2.8759E0_realk&
     &,  3.0051E0_realk,  3.1692E0_realk,  2.9273E0_realk,  2.7615E0_realk,  2.7164E0_realk,  2.6212E0_realk,  3.5744E0_realk&
     &,  3.0275E0_realk,  3.3249E0_realk,  3.3627E0_realk,  3.3686E0_realk,  3.1669E0_realk,  3.0584E0_realk,  2.9915E0_realk&
     &,  3.8773E0_realk,  3.3099E0_realk,  3.2231E0_realk,  3.1600E0_realk,  3.1520E0_realk,  3.4023E0_realk,  3.0426E0_realk&
     &,  3.0099E0_realk,  2.9920E0_realk,  2.9809E0_realk,  2.9800E0_realk,  2.9646E0_realk,  3.3068E0_realk,  3.3930E0_realk&
     &,  3.4486E0_realk,  3.2682E0_realk,  3.1729E0_realk,  3.1168E0_realk,  3.9952E0_realk,  3.4796E0_realk,  3.3901E0_realk&
     &/)
      R0AB(1121:1190)=(/&
     &   3.3255E0_realk,  3.5530E0_realk,  3.5183E0_realk,  3.4097E0_realk,  3.3683E0_realk,  3.3492E0_realk,  3.3360E0_realk&
     &,  3.3308E0_realk,  2.5424E0_realk,  2.6601E0_realk,  3.2555E0_realk,  3.2807E0_realk,  3.1384E0_realk,  3.1737E0_realk&
     &,  2.9397E0_realk,  2.8429E0_realk,  2.8492E0_realk,  2.7225E0_realk,  3.3875E0_realk,  3.4910E0_realk,  3.4520E0_realk&
     &,  3.3608E0_realk,  3.3036E0_realk,  3.2345E0_realk,  3.2999E0_realk,  3.1487E0_realk,  3.7409E0_realk,  3.8392E0_realk&
     &,  3.7148E0_realk,  3.6439E0_realk,  3.6182E0_realk,  3.1753E0_realk,  3.5210E0_realk,  3.4639E0_realk,  3.4265E0_realk&
     &,  3.4075E0_realk,  3.3828E0_realk,  3.3474E0_realk,  3.4071E0_realk,  3.3754E0_realk,  3.3646E0_realk,  3.3308E0_realk&
     &,  3.4393E0_realk,  3.2993E0_realk,  3.8768E0_realk,  3.9891E0_realk,  3.8310E0_realk,  3.7483E0_realk,  3.3417E0_realk&
     &,  3.3019E0_realk,  3.2250E0_realk,  3.1832E0_realk,  3.1578E0_realk,  3.1564E0_realk,  3.1224E0_realk,  3.4620E0_realk&
     &,  2.9743E0_realk,  2.8058E0_realk,  3.4830E0_realk,  3.3474E0_realk,  3.6863E0_realk,  3.3617E0_realk,  3.1608E0_realk&
     &,  3.0069E0_realk,  2.9640E0_realk,  2.8427E0_realk,  3.5885E0_realk,  3.5219E0_realk,  4.1314E0_realk,  3.8120E0_realk&
     &/)
      R0AB(1191:1260)=(/&
     &   3.6015E0_realk,  3.4502E0_realk,  3.3498E0_realk,  3.2777E0_realk,  3.8635E0_realk,  3.8232E0_realk,  3.8486E0_realk&
     &,  3.7215E0_realk,  3.6487E0_realk,  3.4724E0_realk,  3.5627E0_realk,  3.5087E0_realk,  3.4757E0_realk,  3.4517E0_realk&
     &,  3.4423E0_realk,  3.4139E0_realk,  4.1028E0_realk,  3.8388E0_realk,  3.6745E0_realk,  3.5562E0_realk,  3.4806E0_realk&
     &,  3.4272E0_realk,  4.0182E0_realk,  3.9991E0_realk,  4.0007E0_realk,  3.9282E0_realk,  3.7238E0_realk,  3.6498E0_realk&
     &,  3.5605E0_realk,  3.5211E0_realk,  3.5009E0_realk,  3.4859E0_realk,  3.4785E0_realk,  3.5621E0_realk,  4.2623E0_realk&
     &,  3.0775E0_realk,  2.8275E0_realk,  4.0181E0_realk,  3.3385E0_realk,  3.5379E0_realk,  3.5036E0_realk,  3.2589E0_realk&
     &,  3.0804E0_realk,  3.0094E0_realk,  2.9003E0_realk,  4.0869E0_realk,  3.5088E0_realk,  3.9105E0_realk,  3.9833E0_realk&
     &,  3.7176E0_realk,  3.5323E0_realk,  3.4102E0_realk,  3.3227E0_realk,  4.2702E0_realk,  4.0888E0_realk,  3.7560E0_realk&
     &,  3.7687E0_realk,  3.6681E0_realk,  3.6405E0_realk,  3.5569E0_realk,  3.4990E0_realk,  3.4659E0_realk,  3.4433E0_realk&
     &,  3.4330E0_realk,  3.4092E0_realk,  3.8867E0_realk,  4.0190E0_realk,  3.7961E0_realk,  3.6412E0_realk,  3.5405E0_realk&
     &/)
      R0AB(1261:1330)=(/&
     &   3.4681E0_realk,  4.3538E0_realk,  4.2136E0_realk,  3.9381E0_realk,  3.8912E0_realk,  3.9681E0_realk,  3.7909E0_realk&
     &,  3.6774E0_realk,  3.6262E0_realk,  3.5999E0_realk,  3.5823E0_realk,  3.5727E0_realk,  3.5419E0_realk,  4.0245E0_realk&
     &,  4.1874E0_realk,  3.0893E0_realk,  2.7917E0_realk,  3.7262E0_realk,  3.3518E0_realk,  3.4241E0_realk,  3.5433E0_realk&
     &,  3.2773E0_realk,  3.0890E0_realk,  2.9775E0_realk,  2.9010E0_realk,  3.8048E0_realk,  3.5362E0_realk,  3.7746E0_realk&
     &,  3.7911E0_realk,  3.7511E0_realk,  3.5495E0_realk,  3.4149E0_realk,  3.3177E0_realk,  4.0129E0_realk,  3.8370E0_realk&
     &,  3.7739E0_realk,  3.7125E0_realk,  3.7152E0_realk,  3.7701E0_realk,  3.5813E0_realk,  3.5187E0_realk,  3.4835E0_realk&
     &,  3.4595E0_realk,  3.4439E0_realk,  3.4242E0_realk,  3.7476E0_realk,  3.8239E0_realk,  3.8346E0_realk,  3.6627E0_realk&
     &,  3.5479E0_realk,  3.4639E0_realk,  4.1026E0_realk,  3.9733E0_realk,  3.9292E0_realk,  3.8667E0_realk,  3.9513E0_realk&
     &,  3.8959E0_realk,  3.7698E0_realk,  3.7089E0_realk,  3.6765E0_realk,  3.6548E0_realk,  3.6409E0_realk,  3.5398E0_realk&
     &,  3.8759E0_realk,  3.9804E0_realk,  4.0150E0_realk,  2.9091E0_realk,  2.7638E0_realk,  3.5066E0_realk,  3.3377E0_realk&
     &/)
      R0AB(1331:1400)=(/&
     &   3.3481E0_realk,  3.2633E0_realk,  3.1810E0_realk,  3.1428E0_realk,  2.9872E0_realk,  2.8837E0_realk,  3.5929E0_realk&
     &,  3.5183E0_realk,  3.6729E0_realk,  3.6596E0_realk,  3.6082E0_realk,  3.5927E0_realk,  3.4224E0_realk,  3.2997E0_realk&
     &,  3.8190E0_realk,  4.1865E0_realk,  4.1114E0_realk,  4.0540E0_realk,  3.6325E0_realk,  3.5697E0_realk,  3.5561E0_realk&
     &,  3.5259E0_realk,  3.4901E0_realk,  3.4552E0_realk,  3.4315E0_realk,  3.4091E0_realk,  3.6438E0_realk,  3.6879E0_realk&
     &,  3.6832E0_realk,  3.7043E0_realk,  3.5557E0_realk,  3.4466E0_realk,  3.9203E0_realk,  4.2919E0_realk,  4.2196E0_realk&
     &,  4.1542E0_realk,  3.7573E0_realk,  3.7039E0_realk,  3.6546E0_realk,  3.6151E0_realk,  3.5293E0_realk,  3.4849E0_realk&
     &,  3.4552E0_realk,  3.5192E0_realk,  3.7673E0_realk,  3.8359E0_realk,  3.8525E0_realk,  3.8901E0_realk,  2.7806E0_realk&
     &,  2.7209E0_realk,  3.3812E0_realk,  3.4958E0_realk,  3.2913E0_realk,  3.1888E0_realk,  3.0990E0_realk,  3.0394E0_realk&
     &,  2.9789E0_realk,  2.8582E0_realk,  3.4716E0_realk,  3.6883E0_realk,  3.6105E0_realk,  3.5704E0_realk,  3.5059E0_realk&
     &,  3.4619E0_realk,  3.4138E0_realk,  3.2742E0_realk,  3.7080E0_realk,  3.9773E0_realk,  3.9010E0_realk,  3.8409E0_realk&
     &/)
      R0AB(1401:1470)=(/&
     &   3.7944E0_realk,  3.4465E0_realk,  3.7235E0_realk,  3.6808E0_realk,  3.6453E0_realk,  3.6168E0_realk,  3.5844E0_realk&
     &,  3.5576E0_realk,  3.5772E0_realk,  3.5959E0_realk,  3.5768E0_realk,  3.5678E0_realk,  3.5486E0_realk,  3.4228E0_realk&
     &,  3.8107E0_realk,  4.0866E0_realk,  4.0169E0_realk,  3.9476E0_realk,  3.6358E0_realk,  3.5800E0_realk,  3.5260E0_realk&
     &,  3.4838E0_realk,  3.4501E0_realk,  3.4204E0_realk,  3.3553E0_realk,  3.6487E0_realk,  3.6973E0_realk,  3.7398E0_realk&
     &,  3.7405E0_realk,  3.7459E0_realk,  3.7380E0_realk,  2.6848E0_realk,  2.6740E0_realk,  3.2925E0_realk,  3.3386E0_realk&
     &,  3.2473E0_realk,  3.1284E0_realk,  3.0301E0_realk,  2.9531E0_realk,  2.9602E0_realk,  2.8272E0_realk,  3.3830E0_realk&
     &,  3.5358E0_realk,  3.5672E0_realk,  3.5049E0_realk,  3.4284E0_realk,  3.3621E0_realk,  3.3001E0_realk,  3.2451E0_realk&
     &,  3.6209E0_realk,  3.8299E0_realk,  3.7543E0_realk,  3.6920E0_realk,  3.6436E0_realk,  3.3598E0_realk,  3.5701E0_realk&
     &,  3.5266E0_realk,  3.4904E0_realk,  3.4590E0_realk,  3.4364E0_realk,  3.4077E0_realk,  3.5287E0_realk,  3.5280E0_realk&
     &,  3.4969E0_realk,  3.4650E0_realk,  3.4304E0_realk,  3.3963E0_realk,  3.7229E0_realk,  3.9402E0_realk,  3.8753E0_realk&
     &/)
      R0AB(1471:1540)=(/&
     &   3.8035E0_realk,  3.5499E0_realk,  3.4913E0_realk,  3.4319E0_realk,  3.3873E0_realk,  3.3520E0_realk,  3.3209E0_realk&
     &,  3.2948E0_realk,  3.5052E0_realk,  3.6465E0_realk,  3.6696E0_realk,  3.6577E0_realk,  3.6388E0_realk,  3.6142E0_realk&
     &,  3.5889E0_realk,  3.3968E0_realk,  3.0122E0_realk,  4.2241E0_realk,  3.7887E0_realk,  4.0049E0_realk,  3.5384E0_realk&
     &,  3.2698E0_realk,  3.1083E0_realk,  2.9917E0_realk,  2.9057E0_realk,  4.3340E0_realk,  3.9900E0_realk,  4.6588E0_realk&
     &,  4.1278E0_realk,  3.8125E0_realk,  3.6189E0_realk,  3.4851E0_realk,  3.3859E0_realk,  4.6531E0_realk,  4.3134E0_realk&
     &,  4.2258E0_realk,  4.1309E0_realk,  4.0692E0_realk,  4.0944E0_realk,  3.9850E0_realk,  3.9416E0_realk,  3.9112E0_realk&
     &,  3.8873E0_realk,  3.8736E0_realk,  3.8473E0_realk,  4.6027E0_realk,  4.1538E0_realk,  3.8994E0_realk,  3.7419E0_realk&
     &,  3.6356E0_realk,  3.5548E0_realk,  4.8353E0_realk,  4.5413E0_realk,  4.3891E0_realk,  4.3416E0_realk,  4.3243E0_realk&
     &,  4.2753E0_realk,  4.2053E0_realk,  4.1790E0_realk,  4.1685E0_realk,  4.1585E0_realk,  4.1536E0_realk,  4.0579E0_realk&
     &,  4.1980E0_realk,  4.4564E0_realk,  4.2192E0_realk,  4.0528E0_realk,  3.9489E0_realk,  3.8642E0_realk,  5.0567E0_realk&
     &/)
      R0AB(1541:1610)=(/&
     &   3.0630E0_realk,  3.3271E0_realk,  4.0432E0_realk,  4.0046E0_realk,  4.1555E0_realk,  3.7426E0_realk,  3.5130E0_realk&
     &,  3.5174E0_realk,  3.2884E0_realk,  3.1378E0_realk,  4.1894E0_realk,  4.2321E0_realk,  4.1725E0_realk,  4.1833E0_realk&
     &,  3.8929E0_realk,  4.0544E0_realk,  3.8118E0_realk,  3.6414E0_realk,  4.6373E0_realk,  4.6268E0_realk,  4.4750E0_realk&
     &,  4.4134E0_realk,  4.3458E0_realk,  3.8582E0_realk,  4.2583E0_realk,  4.1898E0_realk,  4.1562E0_realk,  4.1191E0_realk&
     &,  4.1069E0_realk,  4.0639E0_realk,  4.1257E0_realk,  4.1974E0_realk,  3.9532E0_realk,  4.1794E0_realk,  3.9660E0_realk&
     &,  3.8130E0_realk,  4.8160E0_realk,  4.8272E0_realk,  4.6294E0_realk,  4.5840E0_realk,  4.0770E0_realk,  4.0088E0_realk&
     &,  3.9103E0_realk,  3.8536E0_realk,  3.8324E0_realk,  3.7995E0_realk,  3.7826E0_realk,  4.2294E0_realk,  4.3380E0_realk&
     &,  4.4352E0_realk,  4.1933E0_realk,  4.4580E0_realk,  4.2554E0_realk,  4.1072E0_realk,  5.0454E0_realk,  5.1814E0_realk&
     &,  3.0632E0_realk,  3.2662E0_realk,  3.6432E0_realk,  3.8088E0_realk,  3.7910E0_realk,  3.7381E0_realk,  3.5093E0_realk&
     &,  3.5155E0_realk,  3.3047E0_realk,  3.1681E0_realk,  3.7871E0_realk,  3.9924E0_realk,  4.0637E0_realk,  4.1382E0_realk&
     &/)
      R0AB(1611:1680)=(/&
     &   3.8591E0_realk,  4.0164E0_realk,  3.7878E0_realk,  3.6316E0_realk,  4.1741E0_realk,  4.3166E0_realk,  4.2395E0_realk&
     &,  4.1831E0_realk,  4.1107E0_realk,  3.5857E0_realk,  4.0270E0_realk,  3.9676E0_realk,  3.9463E0_realk,  3.9150E0_realk&
     &,  3.9021E0_realk,  3.8708E0_realk,  4.0240E0_realk,  4.1551E0_realk,  3.9108E0_realk,  4.1337E0_realk,  3.9289E0_realk&
     &,  3.7873E0_realk,  4.3666E0_realk,  4.5080E0_realk,  4.4232E0_realk,  4.3155E0_realk,  3.8461E0_realk,  3.8007E0_realk&
     &,  3.6991E0_realk,  3.6447E0_realk,  3.6308E0_realk,  3.5959E0_realk,  3.5749E0_realk,  4.0359E0_realk,  4.3124E0_realk&
     &,  4.3539E0_realk,  4.1122E0_realk,  4.3772E0_realk,  4.1785E0_realk,  4.0386E0_realk,  4.7004E0_realk,  4.8604E0_realk&
     &,  4.6261E0_realk,  2.9455E0_realk,  3.2470E0_realk,  3.6108E0_realk,  3.8522E0_realk,  3.6625E0_realk,  3.6598E0_realk&
     &,  3.4411E0_realk,  3.4660E0_realk,  3.2415E0_realk,  3.0944E0_realk,  3.7514E0_realk,  4.0397E0_realk,  3.9231E0_realk&
     &,  4.0561E0_realk,  3.7860E0_realk,  3.9845E0_realk,  3.7454E0_realk,  3.5802E0_realk,  4.1366E0_realk,  4.3581E0_realk&
     &,  4.2351E0_realk,  4.2011E0_realk,  4.1402E0_realk,  3.5381E0_realk,  4.0653E0_realk,  4.0093E0_realk,  3.9883E0_realk&
     &/)
      R0AB(1681:1750)=(/&
     &   3.9570E0_realk,  3.9429E0_realk,  3.9112E0_realk,  3.8728E0_realk,  4.0682E0_realk,  3.8351E0_realk,  4.1054E0_realk&
     &,  3.8928E0_realk,  3.7445E0_realk,  4.3415E0_realk,  4.5497E0_realk,  4.3833E0_realk,  4.3122E0_realk,  3.8051E0_realk&
     &,  3.7583E0_realk,  3.6622E0_realk,  3.6108E0_realk,  3.5971E0_realk,  3.5628E0_realk,  3.5408E0_realk,  4.0780E0_realk&
     &,  4.0727E0_realk,  4.2836E0_realk,  4.0553E0_realk,  4.3647E0_realk,  4.1622E0_realk,  4.0178E0_realk,  4.5802E0_realk&
     &,  4.9125E0_realk,  4.5861E0_realk,  4.6201E0_realk,  2.9244E0_realk,  3.2241E0_realk,  3.5848E0_realk,  3.8293E0_realk&
     &,  3.6395E0_realk,  3.6400E0_realk,  3.4204E0_realk,  3.4499E0_realk,  3.2253E0_realk,  3.0779E0_realk,  3.7257E0_realk&
     &,  4.0170E0_realk,  3.9003E0_realk,  4.0372E0_realk,  3.7653E0_realk,  3.9672E0_realk,  3.7283E0_realk,  3.5630E0_realk&
     &,  4.1092E0_realk,  4.3347E0_realk,  4.2117E0_realk,  4.1793E0_realk,  4.1179E0_realk,  3.5139E0_realk,  4.0426E0_realk&
     &,  3.9867E0_realk,  3.9661E0_realk,  3.9345E0_realk,  3.9200E0_realk,  3.8883E0_realk,  3.8498E0_realk,  4.0496E0_realk&
     &,  3.8145E0_realk,  4.0881E0_realk,  3.8756E0_realk,  3.7271E0_realk,  4.3128E0_realk,  4.5242E0_realk,  4.3578E0_realk&
     &/)
      R0AB(1751:1820)=(/&
     &   4.2870E0_realk,  3.7796E0_realk,  3.7318E0_realk,  3.6364E0_realk,  3.5854E0_realk,  3.5726E0_realk,  3.5378E0_realk&
     &,  3.5155E0_realk,  4.0527E0_realk,  4.0478E0_realk,  4.2630E0_realk,  4.0322E0_realk,  4.3449E0_realk,  4.1421E0_realk&
     &,  3.9975E0_realk,  4.5499E0_realk,  4.8825E0_realk,  4.5601E0_realk,  4.5950E0_realk,  4.5702E0_realk,  2.9046E0_realk&
     &,  3.2044E0_realk,  3.5621E0_realk,  3.8078E0_realk,  3.6185E0_realk,  3.6220E0_realk,  3.4019E0_realk,  3.4359E0_realk&
     &,  3.2110E0_realk,  3.0635E0_realk,  3.7037E0_realk,  3.9958E0_realk,  3.8792E0_realk,  4.0194E0_realk,  3.7460E0_realk&
     &,  3.9517E0_realk,  3.7128E0_realk,  3.5474E0_realk,  4.0872E0_realk,  4.3138E0_realk,  4.1906E0_realk,  4.1593E0_realk&
     &,  4.0973E0_realk,  3.4919E0_realk,  4.0216E0_realk,  3.9657E0_realk,  3.9454E0_realk,  3.9134E0_realk,  3.8986E0_realk&
     &,  3.8669E0_realk,  3.8289E0_realk,  4.0323E0_realk,  3.7954E0_realk,  4.0725E0_realk,  3.8598E0_realk,  3.7113E0_realk&
     &,  4.2896E0_realk,  4.5021E0_realk,  4.3325E0_realk,  4.2645E0_realk,  3.7571E0_realk,  3.7083E0_realk,  3.6136E0_realk&
     &,  3.5628E0_realk,  3.5507E0_realk,  3.5155E0_realk,  3.4929E0_realk,  4.0297E0_realk,  4.0234E0_realk,  4.2442E0_realk&
     &/)
      R0AB(1821:1890)=(/&
     &   4.0112E0_realk,  4.3274E0_realk,  4.1240E0_realk,  3.9793E0_realk,  4.5257E0_realk,  4.8568E0_realk,  4.5353E0_realk&
     &,  4.5733E0_realk,  4.5485E0_realk,  4.5271E0_realk,  2.8878E0_realk,  3.1890E0_realk,  3.5412E0_realk,  3.7908E0_realk&
     &,  3.5974E0_realk,  3.6078E0_realk,  3.3871E0_realk,  3.4243E0_realk,  3.1992E0_realk,  3.0513E0_realk,  3.6831E0_realk&
     &,  3.9784E0_realk,  3.8579E0_realk,  4.0049E0_realk,  3.7304E0_realk,  3.9392E0_realk,  3.7002E0_realk,  3.5347E0_realk&
     &,  4.0657E0_realk,  4.2955E0_realk,  4.1705E0_realk,  4.1424E0_realk,  4.0800E0_realk,  3.4717E0_realk,  4.0043E0_realk&
     &,  3.9485E0_realk,  3.9286E0_realk,  3.8965E0_realk,  3.8815E0_realk,  3.8500E0_realk,  3.8073E0_realk,  4.0180E0_realk&
     &,  3.7796E0_realk,  4.0598E0_realk,  3.8470E0_realk,  3.6983E0_realk,  4.2678E0_realk,  4.4830E0_realk,  4.3132E0_realk&
     &,  4.2444E0_realk,  3.7370E0_realk,  3.6876E0_realk,  3.5935E0_realk,  3.5428E0_realk,  3.5314E0_realk,  3.4958E0_realk&
     &,  3.4730E0_realk,  4.0117E0_realk,  4.0043E0_realk,  4.2287E0_realk,  3.9939E0_realk,  4.3134E0_realk,  4.1096E0_realk&
     &,  3.9646E0_realk,  4.5032E0_realk,  4.8356E0_realk,  4.5156E0_realk,  4.5544E0_realk,  4.5297E0_realk,  4.5083E0_realk&
     &/)
      R0AB(1891:1960)=(/&
     &   4.4896E0_realk,  2.8709E0_realk,  3.1737E0_realk,  3.5199E0_realk,  3.7734E0_realk,  3.5802E0_realk,  3.5934E0_realk&
     &,  3.3724E0_realk,  3.4128E0_realk,  3.1877E0_realk,  3.0396E0_realk,  3.6624E0_realk,  3.9608E0_realk,  3.8397E0_realk&
     &,  3.9893E0_realk,  3.7145E0_realk,  3.9266E0_realk,  3.6877E0_realk,  3.5222E0_realk,  4.0448E0_realk,  4.2771E0_realk&
     &,  4.1523E0_realk,  4.1247E0_realk,  4.0626E0_realk,  3.4530E0_realk,  3.9866E0_realk,  3.9310E0_realk,  3.9115E0_realk&
     &,  3.8792E0_realk,  3.8641E0_realk,  3.8326E0_realk,  3.7892E0_realk,  4.0025E0_realk,  3.7636E0_realk,  4.0471E0_realk&
     &,  3.8343E0_realk,  3.6854E0_realk,  4.2464E0_realk,  4.4635E0_realk,  4.2939E0_realk,  4.2252E0_realk,  3.7169E0_realk&
     &,  3.6675E0_realk,  3.5739E0_realk,  3.5235E0_realk,  3.5126E0_realk,  3.4768E0_realk,  3.4537E0_realk,  3.9932E0_realk&
     &,  3.9854E0_realk,  4.2123E0_realk,  3.9765E0_realk,  4.2992E0_realk,  4.0951E0_realk,  3.9500E0_realk,  4.4811E0_realk&
     &,  4.8135E0_realk,  4.4959E0_realk,  4.5351E0_realk,  4.5105E0_realk,  4.4891E0_realk,  4.4705E0_realk,  4.4515E0_realk&
     &,  2.8568E0_realk,  3.1608E0_realk,  3.5050E0_realk,  3.7598E0_realk,  3.5665E0_realk,  3.5803E0_realk,  3.3601E0_realk&
     &/)
      R0AB(1961:2030)=(/&
     &   3.4031E0_realk,  3.1779E0_realk,  3.0296E0_realk,  3.6479E0_realk,  3.9471E0_realk,  3.8262E0_realk,  3.9773E0_realk&
     &,  3.7015E0_realk,  3.9162E0_realk,  3.6771E0_realk,  3.5115E0_realk,  4.0306E0_realk,  4.2634E0_realk,  4.1385E0_realk&
     &,  4.1116E0_realk,  4.0489E0_realk,  3.4366E0_realk,  3.9732E0_realk,  3.9176E0_realk,  3.8983E0_realk,  3.8659E0_realk&
     &,  3.8507E0_realk,  3.8191E0_realk,  3.7757E0_realk,  3.9907E0_realk,  3.7506E0_realk,  4.0365E0_realk,  3.8235E0_realk&
     &,  3.6745E0_realk,  4.2314E0_realk,  4.4490E0_realk,  4.2792E0_realk,  4.2105E0_realk,  3.7003E0_realk,  3.6510E0_realk&
     &,  3.5578E0_realk,  3.5075E0_realk,  3.4971E0_realk,  3.4609E0_realk,  3.4377E0_realk,  3.9788E0_realk,  3.9712E0_realk&
     &,  4.1997E0_realk,  3.9624E0_realk,  4.2877E0_realk,  4.0831E0_realk,  3.9378E0_realk,  4.4655E0_realk,  4.7974E0_realk&
     &,  4.4813E0_realk,  4.5209E0_realk,  4.4964E0_realk,  4.4750E0_realk,  4.4565E0_realk,  4.4375E0_realk,  4.4234E0_realk&
     &,  2.6798E0_realk,  3.0151E0_realk,  3.2586E0_realk,  3.5292E0_realk,  3.5391E0_realk,  3.4902E0_realk,  3.2887E0_realk&
     &,  3.3322E0_realk,  3.1228E0_realk,  2.9888E0_realk,  3.4012E0_realk,  3.7145E0_realk,  3.7830E0_realk,  3.6665E0_realk&
     &/)
      R0AB(2031:2100)=(/&
     &   3.5898E0_realk,  3.8077E0_realk,  3.5810E0_realk,  3.4265E0_realk,  3.7726E0_realk,  4.0307E0_realk,  3.9763E0_realk&
     &,  3.8890E0_realk,  3.8489E0_realk,  3.2706E0_realk,  3.7595E0_realk,  3.6984E0_realk,  3.6772E0_realk,  3.6428E0_realk&
     &,  3.6243E0_realk,  3.5951E0_realk,  3.7497E0_realk,  3.6775E0_realk,  3.6364E0_realk,  3.9203E0_realk,  3.7157E0_realk&
     &,  3.5746E0_realk,  3.9494E0_realk,  4.2076E0_realk,  4.1563E0_realk,  4.0508E0_realk,  3.5329E0_realk,  3.4780E0_realk&
     &,  3.3731E0_realk,  3.3126E0_realk,  3.2846E0_realk,  3.2426E0_realk,  3.2135E0_realk,  3.7491E0_realk,  3.9006E0_realk&
     &,  3.8332E0_realk,  3.8029E0_realk,  4.1436E0_realk,  3.9407E0_realk,  3.7998E0_realk,  4.1663E0_realk,  4.5309E0_realk&
     &,  4.3481E0_realk,  4.2911E0_realk,  4.2671E0_realk,  4.2415E0_realk,  4.2230E0_realk,  4.2047E0_realk,  4.1908E0_realk&
     &,  4.1243E0_realk,  2.5189E0_realk,  2.9703E0_realk,  3.3063E0_realk,  3.6235E0_realk,  3.4517E0_realk,  3.3989E0_realk&
     &,  3.2107E0_realk,  3.2434E0_realk,  3.0094E0_realk,  2.8580E0_realk,  3.4253E0_realk,  3.8157E0_realk,  3.7258E0_realk&
     &,  3.6132E0_realk,  3.5297E0_realk,  3.7566E0_realk,  3.5095E0_realk,  3.3368E0_realk,  3.7890E0_realk,  4.1298E0_realk&
     &/)
      R0AB(2101:2170)=(/&
     &   4.0190E0_realk,  3.9573E0_realk,  3.9237E0_realk,  3.2677E0_realk,  3.8480E0_realk,  3.8157E0_realk,  3.7656E0_realk&
     &,  3.7317E0_realk,  3.7126E0_realk,  3.6814E0_realk,  3.6793E0_realk,  3.6218E0_realk,  3.5788E0_realk,  3.8763E0_realk&
     &,  3.6572E0_realk,  3.5022E0_realk,  3.9737E0_realk,  4.3255E0_realk,  4.1828E0_realk,  4.1158E0_realk,  3.5078E0_realk&
     &,  3.4595E0_realk,  3.3600E0_realk,  3.3088E0_realk,  3.2575E0_realk,  3.2164E0_realk,  3.1856E0_realk,  3.8522E0_realk&
     &,  3.8665E0_realk,  3.8075E0_realk,  3.7772E0_realk,  4.1391E0_realk,  3.9296E0_realk,  3.7772E0_realk,  4.2134E0_realk&
     &,  4.7308E0_realk,  4.3787E0_realk,  4.3894E0_realk,  4.3649E0_realk,  4.3441E0_realk,  4.3257E0_realk,  4.3073E0_realk&
     &,  4.2941E0_realk,  4.1252E0_realk,  4.2427E0_realk,  3.0481E0_realk,  2.9584E0_realk,  3.6919E0_realk,  3.5990E0_realk&
     &,  3.8881E0_realk,  3.4209E0_realk,  3.1606E0_realk,  3.1938E0_realk,  2.9975E0_realk,  2.8646E0_realk,  3.8138E0_realk&
     &,  3.7935E0_realk,  3.7081E0_realk,  3.9155E0_realk,  3.5910E0_realk,  3.4808E0_realk,  3.4886E0_realk,  3.3397E0_realk&
     &,  4.1336E0_realk,  4.1122E0_realk,  3.9888E0_realk,  3.9543E0_realk,  3.8917E0_realk,  3.5894E0_realk,  3.8131E0_realk&
     &/)
      R0AB(2171:2240)=(/&
     &   3.7635E0_realk,  3.7419E0_realk,  3.7071E0_realk,  3.6880E0_realk,  3.6574E0_realk,  3.6546E0_realk,  3.9375E0_realk&
     &,  3.6579E0_realk,  3.5870E0_realk,  3.6361E0_realk,  3.5039E0_realk,  4.3149E0_realk,  4.2978E0_realk,  4.1321E0_realk&
     &,  4.1298E0_realk,  3.8164E0_realk,  3.7680E0_realk,  3.7154E0_realk,  3.6858E0_realk,  3.6709E0_realk,  3.6666E0_realk&
     &,  3.6517E0_realk,  3.8174E0_realk,  3.8608E0_realk,  4.1805E0_realk,  3.9102E0_realk,  3.8394E0_realk,  3.8968E0_realk&
     &,  3.7673E0_realk,  4.5274E0_realk,  4.6682E0_realk,  4.3344E0_realk,  4.3639E0_realk,  4.3384E0_realk,  4.3162E0_realk&
     &,  4.2972E0_realk,  4.2779E0_realk,  4.2636E0_realk,  4.0253E0_realk,  4.1168E0_realk,  4.1541E0_realk,  2.8136E0_realk&
     &,  3.0951E0_realk,  3.4635E0_realk,  3.6875E0_realk,  3.4987E0_realk,  3.5183E0_realk,  3.2937E0_realk,  3.3580E0_realk&
     &,  3.1325E0_realk,  2.9832E0_realk,  3.6078E0_realk,  3.8757E0_realk,  3.7616E0_realk,  3.9222E0_realk,  3.6370E0_realk&
     &,  3.8647E0_realk,  3.6256E0_realk,  3.4595E0_realk,  3.9874E0_realk,  4.1938E0_realk,  4.0679E0_realk,  4.0430E0_realk&
     &,  3.9781E0_realk,  3.3886E0_realk,  3.9008E0_realk,  3.8463E0_realk,  3.8288E0_realk,  3.7950E0_realk,  3.7790E0_realk&
     &/)
      R0AB(2241:2310)=(/&
     &   3.7472E0_realk,  3.7117E0_realk,  3.9371E0_realk,  3.6873E0_realk,  3.9846E0_realk,  3.7709E0_realk,  3.6210E0_realk&
     &,  4.1812E0_realk,  4.3750E0_realk,  4.2044E0_realk,  4.1340E0_realk,  3.6459E0_realk,  3.5929E0_realk,  3.5036E0_realk&
     &,  3.4577E0_realk,  3.4528E0_realk,  3.4146E0_realk,  3.3904E0_realk,  3.9014E0_realk,  3.9031E0_realk,  4.1443E0_realk&
     &,  3.8961E0_realk,  4.2295E0_realk,  4.0227E0_realk,  3.8763E0_realk,  4.4086E0_realk,  4.7097E0_realk,  4.4064E0_realk&
     &,  4.4488E0_realk,  4.4243E0_realk,  4.4029E0_realk,  4.3842E0_realk,  4.3655E0_realk,  4.3514E0_realk,  4.1162E0_realk&
     &,  4.2205E0_realk,  4.1953E0_realk,  4.2794E0_realk,  2.8032E0_realk,  3.0805E0_realk,  3.4519E0_realk,  3.6700E0_realk&
     &,  3.4827E0_realk,  3.5050E0_realk,  3.2799E0_realk,  3.3482E0_realk,  3.1233E0_realk,  2.9747E0_realk,  3.5971E0_realk&
     &,  3.8586E0_realk,  3.7461E0_realk,  3.9100E0_realk,  3.6228E0_realk,  3.8535E0_realk,  3.6147E0_realk,  3.4490E0_realk&
     &,  3.9764E0_realk,  4.1773E0_realk,  4.0511E0_realk,  4.0270E0_realk,  3.9614E0_realk,  3.3754E0_realk,  3.8836E0_realk&
     &,  3.8291E0_realk,  3.8121E0_realk,  3.7780E0_realk,  3.7619E0_realk,  3.7300E0_realk,  3.6965E0_realk,  3.9253E0_realk&
     &/)
      R0AB(2311:2380)=(/&
     &   3.6734E0_realk,  3.9733E0_realk,  3.7597E0_realk,  3.6099E0_realk,  4.1683E0_realk,  4.3572E0_realk,  4.1862E0_realk&
     &,  4.1153E0_realk,  3.6312E0_realk,  3.5772E0_realk,  3.4881E0_realk,  3.4429E0_realk,  3.4395E0_realk,  3.4009E0_realk&
     &,  3.3766E0_realk,  3.8827E0_realk,  3.8868E0_realk,  4.1316E0_realk,  3.8807E0_realk,  4.2164E0_realk,  4.0092E0_realk&
     &,  3.8627E0_realk,  4.3936E0_realk,  4.6871E0_realk,  4.3882E0_realk,  4.4316E0_realk,  4.4073E0_realk,  4.3858E0_realk&
     &,  4.3672E0_realk,  4.3485E0_realk,  4.3344E0_realk,  4.0984E0_realk,  4.2036E0_realk,  4.1791E0_realk,  4.2622E0_realk&
     &,  4.2450E0_realk,  2.7967E0_realk,  3.0689E0_realk,  3.4445E0_realk,  3.6581E0_realk,  3.4717E0_realk,  3.4951E0_realk&
     &,  3.2694E0_realk,  3.3397E0_realk,  3.1147E0_realk,  2.9661E0_realk,  3.5898E0_realk,  3.8468E0_realk,  3.7358E0_realk&
     &,  3.9014E0_realk,  3.6129E0_realk,  3.8443E0_realk,  3.6054E0_realk,  3.4396E0_realk,  3.9683E0_realk,  4.1656E0_realk&
     &,  4.0394E0_realk,  4.0158E0_realk,  3.9498E0_realk,  3.3677E0_realk,  3.8718E0_realk,  3.8164E0_realk,  3.8005E0_realk&
     &,  3.7662E0_realk,  3.7500E0_realk,  3.7181E0_realk,  3.6863E0_realk,  3.9170E0_realk,  3.6637E0_realk,  3.9641E0_realk&
     &/)
      R0AB(2381:2450)=(/&
     &   3.7503E0_realk,  3.6004E0_realk,  4.1590E0_realk,  4.3448E0_realk,  4.1739E0_realk,  4.1029E0_realk,  3.6224E0_realk&
     &,  3.5677E0_realk,  3.4785E0_realk,  3.4314E0_realk,  3.4313E0_realk,  3.3923E0_realk,  3.3680E0_realk,  3.8698E0_realk&
     &,  3.8758E0_realk,  4.1229E0_realk,  3.8704E0_realk,  4.2063E0_realk,  3.9987E0_realk,  3.8519E0_realk,  4.3832E0_realk&
     &,  4.6728E0_realk,  4.3759E0_realk,  4.4195E0_realk,  4.3952E0_realk,  4.3737E0_realk,  4.3551E0_realk,  4.3364E0_realk&
     &,  4.3223E0_realk,  4.0861E0_realk,  4.1911E0_realk,  4.1676E0_realk,  4.2501E0_realk,  4.2329E0_realk,  4.2208E0_realk&
     &,  2.7897E0_realk,  3.0636E0_realk,  3.4344E0_realk,  3.6480E0_realk,  3.4626E0_realk,  3.4892E0_realk,  3.2626E0_realk&
     &,  3.3344E0_realk,  3.1088E0_realk,  2.9597E0_realk,  3.5804E0_realk,  3.8359E0_realk,  3.7251E0_realk,  3.8940E0_realk&
     &,  3.6047E0_realk,  3.8375E0_realk,  3.5990E0_realk,  3.4329E0_realk,  3.9597E0_realk,  4.1542E0_realk,  4.0278E0_realk&
     &,  4.0048E0_realk,  3.9390E0_realk,  3.3571E0_realk,  3.8608E0_realk,  3.8056E0_realk,  3.7899E0_realk,  3.7560E0_realk&
     &,  3.7400E0_realk,  3.7081E0_realk,  3.6758E0_realk,  3.9095E0_realk,  3.6552E0_realk,  3.9572E0_realk,  3.7436E0_realk&
     &/)
      R0AB(2451:2520)=(/&
     &   3.5933E0_realk,  4.1508E0_realk,  4.3337E0_realk,  4.1624E0_realk,  4.0916E0_realk,  3.6126E0_realk,  3.5582E0_realk&
     &,  3.4684E0_realk,  3.4212E0_realk,  3.4207E0_realk,  3.3829E0_realk,  3.3586E0_realk,  3.8604E0_realk,  3.8658E0_realk&
     &,  4.1156E0_realk,  3.8620E0_realk,  4.1994E0_realk,  3.9917E0_realk,  3.8446E0_realk,  4.3750E0_realk,  4.6617E0_realk&
     &,  4.3644E0_realk,  4.4083E0_realk,  4.3840E0_realk,  4.3625E0_realk,  4.3439E0_realk,  4.3253E0_realk,  4.3112E0_realk&
     &,  4.0745E0_realk,  4.1807E0_realk,  4.1578E0_realk,  4.2390E0_realk,  4.2218E0_realk,  4.2097E0_realk,  4.1986E0_realk&
     &,  2.8395E0_realk,  3.0081E0_realk,  3.3171E0_realk,  3.4878E0_realk,  3.5360E0_realk,  3.5145E0_realk,  3.2809E0_realk&
     &,  3.3307E0_realk,  3.1260E0_realk,  2.9940E0_realk,  3.4741E0_realk,  3.6675E0_realk,  3.7832E0_realk,  3.6787E0_realk&
     &,  3.6156E0_realk,  3.8041E0_realk,  3.5813E0_realk,  3.4301E0_realk,  3.8480E0_realk,  3.9849E0_realk,  3.9314E0_realk&
     &,  3.8405E0_realk,  3.8029E0_realk,  3.2962E0_realk,  3.7104E0_realk,  3.6515E0_realk,  3.6378E0_realk,  3.6020E0_realk&
     &,  3.5849E0_realk,  3.5550E0_realk,  3.7494E0_realk,  3.6893E0_realk,  3.6666E0_realk,  3.9170E0_realk,  3.7150E0_realk&
     &/)
      R0AB(2521:2590)=(/&
     &   3.5760E0_realk,  4.0268E0_realk,  4.1596E0_realk,  4.1107E0_realk,  3.9995E0_realk,  3.5574E0_realk,  3.5103E0_realk&
     &,  3.4163E0_realk,  3.3655E0_realk,  3.3677E0_realk,  3.3243E0_realk,  3.2975E0_realk,  3.7071E0_realk,  3.9047E0_realk&
     &,  3.8514E0_realk,  3.8422E0_realk,  3.8022E0_realk,  3.9323E0_realk,  3.7932E0_realk,  4.2343E0_realk,  4.4583E0_realk&
     &,  4.3115E0_realk,  4.2457E0_realk,  4.2213E0_realk,  4.1945E0_realk,  4.1756E0_realk,  4.1569E0_realk,  4.1424E0_realk&
     &,  4.0620E0_realk,  4.0494E0_realk,  3.9953E0_realk,  4.0694E0_realk,  4.0516E0_realk,  4.0396E0_realk,  4.0280E0_realk&
     &,  4.0130E0_realk,  2.9007E0_realk,  2.9674E0_realk,  3.8174E0_realk,  3.5856E0_realk,  3.6486E0_realk,  3.5339E0_realk&
     &,  3.2832E0_realk,  3.3154E0_realk,  3.1144E0_realk,  2.9866E0_realk,  3.9618E0_realk,  3.8430E0_realk,  3.9980E0_realk&
     &,  3.8134E0_realk,  3.6652E0_realk,  3.7985E0_realk,  3.5756E0_realk,  3.4207E0_realk,  4.4061E0_realk,  4.2817E0_realk&
     &,  4.1477E0_realk,  4.0616E0_realk,  3.9979E0_realk,  3.6492E0_realk,  3.8833E0_realk,  3.8027E0_realk,  3.7660E0_realk&
     &,  3.7183E0_realk,  3.6954E0_realk,  3.6525E0_realk,  3.9669E0_realk,  3.8371E0_realk,  3.7325E0_realk,  3.9160E0_realk&
     &/)
      R0AB(2591:2660)=(/&
     &   3.7156E0_realk,  3.5714E0_realk,  4.6036E0_realk,  4.4620E0_realk,  4.3092E0_realk,  4.2122E0_realk,  3.8478E0_realk&
     &,  3.7572E0_realk,  3.6597E0_realk,  3.5969E0_realk,  3.5575E0_realk,  3.5386E0_realk,  3.5153E0_realk,  3.7818E0_realk&
     &,  4.1335E0_realk,  4.0153E0_realk,  3.9177E0_realk,  3.8603E0_realk,  3.9365E0_realk,  3.7906E0_realk,  4.7936E0_realk&
     &,  4.7410E0_realk,  4.5461E0_realk,  4.5662E0_realk,  4.5340E0_realk,  4.5059E0_realk,  4.4832E0_realk,  4.4604E0_realk&
     &,  4.4429E0_realk,  4.2346E0_realk,  4.4204E0_realk,  4.3119E0_realk,  4.3450E0_realk,  4.3193E0_realk,  4.3035E0_realk&
     &,  4.2933E0_realk,  4.1582E0_realk,  4.2450E0_realk,  2.8559E0_realk,  2.9050E0_realk,  3.8325E0_realk,  3.5442E0_realk&
     &,  3.5077E0_realk,  3.4905E0_realk,  3.2396E0_realk,  3.2720E0_realk,  3.0726E0_realk,  2.9467E0_realk,  3.9644E0_realk&
     &,  3.8050E0_realk,  3.8981E0_realk,  3.7762E0_realk,  3.6216E0_realk,  3.7531E0_realk,  3.5297E0_realk,  3.3742E0_realk&
     &,  4.3814E0_realk,  4.2818E0_realk,  4.1026E0_realk,  4.0294E0_realk,  3.9640E0_realk,  3.6208E0_realk,  3.8464E0_realk&
     &,  3.7648E0_realk,  3.7281E0_realk,  3.6790E0_realk,  3.6542E0_realk,  3.6117E0_realk,  3.8650E0_realk,  3.8010E0_realk&
     &/)
      R0AB(2661:2730)=(/&
     &   3.6894E0_realk,  3.8713E0_realk,  3.6699E0_realk,  3.5244E0_realk,  4.5151E0_realk,  4.4517E0_realk,  4.2538E0_realk&
     &,  4.1483E0_realk,  3.8641E0_realk,  3.7244E0_realk,  3.6243E0_realk,  3.5589E0_realk,  3.5172E0_realk,  3.4973E0_realk&
     &,  3.4715E0_realk,  3.7340E0_realk,  4.0316E0_realk,  3.9958E0_realk,  3.8687E0_realk,  3.8115E0_realk,  3.8862E0_realk&
     &,  3.7379E0_realk,  4.7091E0_realk,  4.7156E0_realk,  4.5199E0_realk,  4.5542E0_realk,  4.5230E0_realk,  4.4959E0_realk&
     &,  4.4750E0_realk,  4.4529E0_realk,  4.4361E0_realk,  4.1774E0_realk,  4.3774E0_realk,  4.2963E0_realk,  4.3406E0_realk&
     &,  4.3159E0_realk,  4.3006E0_realk,  4.2910E0_realk,  4.1008E0_realk,  4.1568E0_realk,  4.0980E0_realk,  2.8110E0_realk&
     &,  2.8520E0_realk,  3.7480E0_realk,  3.5105E0_realk,  3.4346E0_realk,  3.3461E0_realk,  3.1971E0_realk,  3.2326E0_realk&
     &,  3.0329E0_realk,  2.9070E0_realk,  3.8823E0_realk,  3.7928E0_realk,  3.8264E0_realk,  3.7006E0_realk,  3.5797E0_realk&
     &,  3.7141E0_realk,  3.4894E0_realk,  3.3326E0_realk,  4.3048E0_realk,  4.2217E0_realk,  4.0786E0_realk,  3.9900E0_realk&
     &,  3.9357E0_realk,  3.6331E0_realk,  3.8333E0_realk,  3.7317E0_realk,  3.6957E0_realk,  3.6460E0_realk,  3.6197E0_realk&
     &/)
      R0AB(2731:2800)=(/&
     &   3.5779E0_realk,  3.7909E0_realk,  3.7257E0_realk,  3.6476E0_realk,  3.5729E0_realk,  3.6304E0_realk,  3.4834E0_realk&
     &,  4.4368E0_realk,  4.3921E0_realk,  4.2207E0_realk,  4.1133E0_realk,  3.8067E0_realk,  3.7421E0_realk,  3.6140E0_realk&
     &,  3.5491E0_realk,  3.5077E0_realk,  3.4887E0_realk,  3.4623E0_realk,  3.6956E0_realk,  3.9568E0_realk,  3.8976E0_realk&
     &,  3.8240E0_realk,  3.7684E0_realk,  3.8451E0_realk,  3.6949E0_realk,  4.6318E0_realk,  4.6559E0_realk,  4.4533E0_realk&
     &,  4.4956E0_realk,  4.4641E0_realk,  4.4366E0_realk,  4.4155E0_realk,  4.3936E0_realk,  4.3764E0_realk,  4.1302E0_realk&
     &,  4.3398E0_realk,  4.2283E0_realk,  4.2796E0_realk,  4.2547E0_realk,  4.2391E0_realk,  4.2296E0_realk,  4.0699E0_realk&
     &,  4.1083E0_realk,  4.0319E0_realk,  3.9855E0_realk,  2.7676E0_realk,  2.8078E0_realk,  3.6725E0_realk,  3.4804E0_realk&
     &,  3.3775E0_realk,  3.2411E0_realk,  3.1581E0_realk,  3.1983E0_realk,  2.9973E0_realk,  2.8705E0_realk,  3.8070E0_realk&
     &,  3.7392E0_realk,  3.7668E0_realk,  3.6263E0_realk,  3.5402E0_realk,  3.6807E0_realk,  3.4545E0_realk,  3.2962E0_realk&
     &,  4.2283E0_realk,  4.1698E0_realk,  4.0240E0_realk,  3.9341E0_realk,  3.8711E0_realk,  3.5489E0_realk,  3.7798E0_realk&
     &/)
      R0AB(2801:2870)=(/&
     &   3.7000E0_realk,  3.6654E0_realk,  3.6154E0_realk,  3.5882E0_realk,  3.5472E0_realk,  3.7289E0_realk,  3.6510E0_realk&
     &,  3.6078E0_realk,  3.5355E0_realk,  3.5963E0_realk,  3.4480E0_realk,  4.3587E0_realk,  4.3390E0_realk,  4.1635E0_realk&
     &,  4.0536E0_realk,  3.7193E0_realk,  3.6529E0_realk,  3.5512E0_realk,  3.4837E0_realk,  3.4400E0_realk,  3.4191E0_realk&
     &,  3.3891E0_realk,  3.6622E0_realk,  3.8934E0_realk,  3.8235E0_realk,  3.7823E0_realk,  3.7292E0_realk,  3.8106E0_realk&
     &,  3.6589E0_realk,  4.5535E0_realk,  4.6013E0_realk,  4.3961E0_realk,  4.4423E0_realk,  4.4109E0_realk,  4.3835E0_realk&
     &,  4.3625E0_realk,  4.3407E0_realk,  4.3237E0_realk,  4.0863E0_realk,  4.2835E0_realk,  4.1675E0_realk,  4.2272E0_realk&
     &,  4.2025E0_realk,  4.1869E0_realk,  4.1774E0_realk,  4.0126E0_realk,  4.0460E0_realk,  3.9815E0_realk,  3.9340E0_realk&
     &,  3.8955E0_realk,  2.6912E0_realk,  2.7604E0_realk,  3.6037E0_realk,  3.4194E0_realk,  3.3094E0_realk,  3.1710E0_realk&
     &,  3.0862E0_realk,  3.1789E0_realk,  2.9738E0_realk,  2.8427E0_realk,  3.7378E0_realk,  3.6742E0_realk,  3.6928E0_realk&
     &,  3.5512E0_realk,  3.4614E0_realk,  3.4087E0_realk,  3.4201E0_realk,  3.2607E0_realk,  4.1527E0_realk,  4.0977E0_realk&
     &/)
      R0AB(2871:2940)=(/&
     &   3.9523E0_realk,  3.8628E0_realk,  3.8002E0_realk,  3.4759E0_realk,  3.7102E0_realk,  3.6466E0_realk,  3.6106E0_realk&
     &,  3.5580E0_realk,  3.5282E0_realk,  3.4878E0_realk,  3.6547E0_realk,  3.5763E0_realk,  3.5289E0_realk,  3.5086E0_realk&
     &,  3.5593E0_realk,  3.4099E0_realk,  4.2788E0_realk,  4.2624E0_realk,  4.0873E0_realk,  3.9770E0_realk,  3.6407E0_realk&
     &,  3.5743E0_realk,  3.5178E0_realk,  3.4753E0_realk,  3.3931E0_realk,  3.3694E0_realk,  3.3339E0_realk,  3.6002E0_realk&
     &,  3.8164E0_realk,  3.7478E0_realk,  3.7028E0_realk,  3.6952E0_realk,  3.7669E0_realk,  3.6137E0_realk,  4.4698E0_realk&
     &,  4.5488E0_realk,  4.3168E0_realk,  4.3646E0_realk,  4.3338E0_realk,  4.3067E0_realk,  4.2860E0_realk,  4.2645E0_realk&
     &,  4.2478E0_realk,  4.0067E0_realk,  4.2349E0_realk,  4.0958E0_realk,  4.1543E0_realk,  4.1302E0_realk,  4.1141E0_realk&
     &,  4.1048E0_realk,  3.9410E0_realk,  3.9595E0_realk,  3.8941E0_realk,  3.8465E0_realk,  3.8089E0_realk,  3.7490E0_realk&
     &,  2.7895E0_realk,  2.5849E0_realk,  3.6484E0_realk,  3.0162E0_realk,  3.1267E0_realk,  3.2125E0_realk,  3.0043E0_realk&
     &,  2.9572E0_realk,  2.8197E0_realk,  2.7261E0_realk,  3.7701E0_realk,  3.2446E0_realk,  3.5239E0_realk,  3.4696E0_realk&
     &/)
      R0AB(2941:3010)=(/&
     &   3.4261E0_realk,  3.3508E0_realk,  3.1968E0_realk,  3.0848E0_realk,  4.1496E0_realk,  3.6598E0_realk,  3.5111E0_realk&
     &,  3.4199E0_realk,  3.3809E0_realk,  3.5382E0_realk,  3.2572E0_realk,  3.2100E0_realk,  3.1917E0_realk,  3.1519E0_realk&
     &,  3.1198E0_realk,  3.1005E0_realk,  3.5071E0_realk,  3.5086E0_realk,  3.5073E0_realk,  3.4509E0_realk,  3.3120E0_realk&
     &,  3.2082E0_realk,  4.2611E0_realk,  3.8117E0_realk,  3.6988E0_realk,  3.5646E0_realk,  3.6925E0_realk,  3.6295E0_realk&
     &,  3.5383E0_realk,  3.4910E0_realk,  3.4625E0_realk,  3.4233E0_realk,  3.4007E0_realk,  3.2329E0_realk,  3.6723E0_realk&
     &,  3.6845E0_realk,  3.6876E0_realk,  3.6197E0_realk,  3.4799E0_realk,  3.3737E0_realk,  4.4341E0_realk,  4.0525E0_realk&
     &,  3.9011E0_realk,  3.8945E0_realk,  3.8635E0_realk,  3.8368E0_realk,  3.8153E0_realk,  3.7936E0_realk,  3.7758E0_realk&
     &,  3.4944E0_realk,  3.4873E0_realk,  3.9040E0_realk,  3.7110E0_realk,  3.6922E0_realk,  3.6799E0_realk,  3.6724E0_realk&
     &,  3.5622E0_realk,  3.6081E0_realk,  3.5426E0_realk,  3.4922E0_realk,  3.4498E0_realk,  3.3984E0_realk,  3.4456E0_realk&
     &,  2.7522E0_realk,  2.5524E0_realk,  3.5742E0_realk,  2.9508E0_realk,  3.0751E0_realk,  3.0158E0_realk,  2.9644E0_realk&
     &/)
      R0AB(3011:3080)=(/&
     &   2.8338E0_realk,  2.7891E0_realk,  2.6933E0_realk,  3.6926E0_realk,  3.1814E0_realk,  3.4528E0_realk,  3.4186E0_realk&
     &,  3.3836E0_realk,  3.2213E0_realk,  3.1626E0_realk,  3.0507E0_realk,  4.0548E0_realk,  3.5312E0_realk,  3.4244E0_realk&
     &,  3.3409E0_realk,  3.2810E0_realk,  3.4782E0_realk,  3.1905E0_realk,  3.1494E0_realk,  3.1221E0_realk,  3.1128E0_realk&
     &,  3.0853E0_realk,  3.0384E0_realk,  3.4366E0_realk,  3.4562E0_realk,  3.4638E0_realk,  3.3211E0_realk,  3.2762E0_realk&
     &,  3.1730E0_realk,  4.1632E0_realk,  3.6825E0_realk,  3.5822E0_realk,  3.4870E0_realk,  3.6325E0_realk,  3.5740E0_realk&
     &,  3.4733E0_realk,  3.4247E0_realk,  3.3969E0_realk,  3.3764E0_realk,  3.3525E0_realk,  3.1984E0_realk,  3.5989E0_realk&
     &,  3.6299E0_realk,  3.6433E0_realk,  3.4937E0_realk,  3.4417E0_realk,  3.3365E0_realk,  4.3304E0_realk,  3.9242E0_realk&
     &,  3.7793E0_realk,  3.7623E0_realk,  3.7327E0_realk,  3.7071E0_realk,  3.6860E0_realk,  3.6650E0_realk,  3.6476E0_realk&
     &,  3.3849E0_realk,  3.3534E0_realk,  3.8216E0_realk,  3.5870E0_realk,  3.5695E0_realk,  3.5584E0_realk,  3.5508E0_realk&
     &,  3.4856E0_realk,  3.5523E0_realk,  3.4934E0_realk,  3.4464E0_realk,  3.4055E0_realk,  3.3551E0_realk,  3.3888E0_realk&
     &/)
      R0AB(3081:3150)=(/&
     &   3.3525E0_realk,  2.7202E0_realk,  2.5183E0_realk,  3.4947E0_realk,  2.8731E0_realk,  3.0198E0_realk,  3.1457E0_realk&
     &,  2.9276E0_realk,  2.7826E0_realk,  2.7574E0_realk,  2.6606E0_realk,  3.6090E0_realk,  3.0581E0_realk,  3.3747E0_realk&
     &,  3.3677E0_realk,  3.3450E0_realk,  3.1651E0_realk,  3.1259E0_realk,  3.0147E0_realk,  3.9498E0_realk,  3.3857E0_realk&
     &,  3.2917E0_realk,  3.2154E0_realk,  3.1604E0_realk,  3.4174E0_realk,  3.0735E0_realk,  3.0342E0_realk,  3.0096E0_realk&
     &,  3.0136E0_realk,  2.9855E0_realk,  2.9680E0_realk,  3.3604E0_realk,  3.4037E0_realk,  3.4243E0_realk,  3.2633E0_realk&
     &,  3.1810E0_realk,  3.1351E0_realk,  4.0557E0_realk,  3.5368E0_realk,  3.4526E0_realk,  3.3699E0_realk,  3.5707E0_realk&
     &,  3.5184E0_realk,  3.4085E0_realk,  3.3595E0_realk,  3.3333E0_realk,  3.3143E0_realk,  3.3041E0_realk,  3.1094E0_realk&
     &,  3.5193E0_realk,  3.5745E0_realk,  3.6025E0_realk,  3.4338E0_realk,  3.3448E0_realk,  3.2952E0_realk,  4.2158E0_realk&
     &,  3.7802E0_realk,  3.6431E0_realk,  3.6129E0_realk,  3.5853E0_realk,  3.5610E0_realk,  3.5406E0_realk,  3.5204E0_realk&
     &,  3.5036E0_realk,  3.2679E0_realk,  3.2162E0_realk,  3.7068E0_realk,  3.4483E0_realk,  3.4323E0_realk,  3.4221E0_realk&
     &/)
      R0AB(3151:3220)=(/&
     &   3.4138E0_realk,  3.3652E0_realk,  3.4576E0_realk,  3.4053E0_realk,  3.3618E0_realk,  3.3224E0_realk,  3.2711E0_realk&
     &,  3.3326E0_realk,  3.2950E0_realk,  3.2564E0_realk,  2.5315E0_realk,  2.6104E0_realk,  3.2734E0_realk,  3.2299E0_realk&
     &,  3.1090E0_realk,  2.9942E0_realk,  2.9159E0_realk,  2.8324E0_realk,  2.8350E0_realk,  2.7216E0_realk,  3.3994E0_realk&
     &,  3.4475E0_realk,  3.4354E0_realk,  3.3438E0_realk,  3.2807E0_realk,  3.2169E0_realk,  3.2677E0_realk,  3.1296E0_realk&
     &,  3.7493E0_realk,  3.8075E0_realk,  3.6846E0_realk,  3.6104E0_realk,  3.5577E0_realk,  3.2052E0_realk,  3.4803E0_realk&
     &,  3.4236E0_realk,  3.3845E0_realk,  3.3640E0_realk,  3.3365E0_realk,  3.3010E0_realk,  3.3938E0_realk,  3.3624E0_realk&
     &,  3.3440E0_realk,  3.3132E0_realk,  3.4035E0_realk,  3.2754E0_realk,  3.8701E0_realk,  3.9523E0_realk,  3.8018E0_realk&
     &,  3.7149E0_realk,  3.3673E0_realk,  3.3199E0_realk,  3.2483E0_realk,  3.2069E0_realk,  3.1793E0_realk,  3.1558E0_realk&
     &,  3.1395E0_realk,  3.4097E0_realk,  3.5410E0_realk,  3.5228E0_realk,  3.5116E0_realk,  3.4921E0_realk,  3.4781E0_realk&
     &,  3.4690E0_realk,  4.0420E0_realk,  4.1759E0_realk,  4.0078E0_realk,  4.0450E0_realk,  4.0189E0_realk,  3.9952E0_realk&
     &/)
      R0AB(3221:3290)=(/&
     &   3.9770E0_realk,  3.9583E0_realk,  3.9434E0_realk,  3.7217E0_realk,  3.8228E0_realk,  3.7826E0_realk,  3.8640E0_realk&
     &,  3.8446E0_realk,  3.8314E0_realk,  3.8225E0_realk,  3.6817E0_realk,  3.7068E0_realk,  3.6555E0_realk,  3.6159E0_realk&
     &,  3.5831E0_realk,  3.5257E0_realk,  3.2133E0_realk,  3.1689E0_realk,  3.1196E0_realk,  3.3599E0_realk,  2.9852E0_realk&
     &,  2.7881E0_realk,  3.5284E0_realk,  3.3493E0_realk,  3.6958E0_realk,  3.3642E0_realk,  3.1568E0_realk,  3.0055E0_realk&
     &,  2.9558E0_realk,  2.8393E0_realk,  3.6287E0_realk,  3.5283E0_realk,  4.1511E0_realk,  3.8259E0_realk,  3.6066E0_realk&
     &,  3.4527E0_realk,  3.3480E0_realk,  3.2713E0_realk,  3.9037E0_realk,  3.8361E0_realk,  3.8579E0_realk,  3.7311E0_realk&
     &,  3.6575E0_realk,  3.5176E0_realk,  3.5693E0_realk,  3.5157E0_realk,  3.4814E0_realk,  3.4559E0_realk,  3.4445E0_realk&
     &,  3.4160E0_realk,  4.1231E0_realk,  3.8543E0_realk,  3.6816E0_realk,  3.5602E0_realk,  3.4798E0_realk,  3.4208E0_realk&
     &,  4.0542E0_realk,  4.0139E0_realk,  4.0165E0_realk,  3.9412E0_realk,  3.7698E0_realk,  3.6915E0_realk,  3.6043E0_realk&
     &,  3.5639E0_realk,  3.5416E0_realk,  3.5247E0_realk,  3.5153E0_realk,  3.5654E0_realk,  4.2862E0_realk,  4.0437E0_realk&
     &/)
      R0AB(3291:3360)=(/&
     &   3.8871E0_realk,  3.7741E0_realk,  3.6985E0_realk,  3.6413E0_realk,  4.2345E0_realk,  4.3663E0_realk,  4.3257E0_realk&
     &,  4.0869E0_realk,  4.0612E0_realk,  4.0364E0_realk,  4.0170E0_realk,  3.9978E0_realk,  3.9834E0_realk,  3.9137E0_realk&
     &,  3.8825E0_realk,  3.8758E0_realk,  3.9143E0_realk,  3.8976E0_realk,  3.8864E0_realk,  3.8768E0_realk,  3.9190E0_realk&
     &,  4.1613E0_realk,  4.0566E0_realk,  3.9784E0_realk,  3.9116E0_realk,  3.8326E0_realk,  3.7122E0_realk,  3.6378E0_realk&
     &,  3.5576E0_realk,  3.5457E0_realk,  4.3127E0_realk,  3.1160E0_realk,  2.8482E0_realk,  4.0739E0_realk,  3.3599E0_realk&
     &,  3.5698E0_realk,  3.5366E0_realk,  3.2854E0_realk,  3.1039E0_realk,  2.9953E0_realk,  2.9192E0_realk,  4.1432E0_realk&
     &,  3.5320E0_realk,  3.9478E0_realk,  4.0231E0_realk,  3.7509E0_realk,  3.5604E0_realk,  3.4340E0_realk,  3.3426E0_realk&
     &,  4.3328E0_realk,  3.8288E0_realk,  3.7822E0_realk,  3.7909E0_realk,  3.6907E0_realk,  3.6864E0_realk,  3.5793E0_realk&
     &,  3.5221E0_realk,  3.4883E0_realk,  3.4649E0_realk,  3.4514E0_realk,  3.4301E0_realk,  3.9256E0_realk,  4.0596E0_realk&
     &,  3.8307E0_realk,  3.6702E0_realk,  3.5651E0_realk,  3.4884E0_realk,  4.4182E0_realk,  4.2516E0_realk,  3.9687E0_realk&
     &/)
      R0AB(3361:3430)=(/&
     &   3.9186E0_realk,  3.9485E0_realk,  3.8370E0_realk,  3.7255E0_realk,  3.6744E0_realk,  3.6476E0_realk,  3.6295E0_realk&
     &,  3.6193E0_realk,  3.5659E0_realk,  4.0663E0_realk,  4.2309E0_realk,  4.0183E0_realk,  3.8680E0_realk,  3.7672E0_realk&
     &,  3.6923E0_realk,  4.5240E0_realk,  4.4834E0_realk,  4.1570E0_realk,  4.3204E0_realk,  4.2993E0_realk,  4.2804E0_realk&
     &,  4.2647E0_realk,  4.2481E0_realk,  4.2354E0_realk,  3.8626E0_realk,  3.8448E0_realk,  4.2267E0_realk,  4.1799E0_realk&
     &,  4.1670E0_realk,  3.8738E0_realk,  3.8643E0_realk,  3.8796E0_realk,  4.0575E0_realk,  4.0354E0_realk,  3.9365E0_realk&
     &,  3.8611E0_realk,  3.7847E0_realk,  3.7388E0_realk,  3.6826E0_realk,  3.6251E0_realk,  3.5492E0_realk,  4.0889E0_realk&
     &,  4.2764E0_realk,  3.1416E0_realk,  2.8325E0_realk,  3.7735E0_realk,  3.3787E0_realk,  3.4632E0_realk,  3.5923E0_realk&
     &,  3.3214E0_realk,  3.1285E0_realk,  3.0147E0_realk,  2.9366E0_realk,  3.8527E0_realk,  3.5602E0_realk,  3.8131E0_realk&
     &,  3.8349E0_realk,  3.7995E0_realk,  3.5919E0_realk,  3.4539E0_realk,  3.3540E0_realk,  4.0654E0_realk,  3.8603E0_realk&
     &,  3.7972E0_realk,  3.7358E0_realk,  3.7392E0_realk,  3.8157E0_realk,  3.6055E0_realk,  3.5438E0_realk,  3.5089E0_realk&
     &/)
      R0AB(3431:3500)=(/&
     &   3.4853E0_realk,  3.4698E0_realk,  3.4508E0_realk,  3.7882E0_realk,  3.8682E0_realk,  3.8837E0_realk,  3.7055E0_realk&
     &,  3.5870E0_realk,  3.5000E0_realk,  4.1573E0_realk,  4.0005E0_realk,  3.9568E0_realk,  3.8936E0_realk,  3.9990E0_realk&
     &,  3.9433E0_realk,  3.8172E0_realk,  3.7566E0_realk,  3.7246E0_realk,  3.7033E0_realk,  3.6900E0_realk,  3.5697E0_realk&
     &,  3.9183E0_realk,  4.0262E0_realk,  4.0659E0_realk,  3.8969E0_realk,  3.7809E0_realk,  3.6949E0_realk,  4.2765E0_realk&
     &,  4.2312E0_realk,  4.1401E0_realk,  4.0815E0_realk,  4.0580E0_realk,  4.0369E0_realk,  4.0194E0_realk,  4.0017E0_realk&
     &,  3.9874E0_realk,  3.8312E0_realk,  3.8120E0_realk,  3.9454E0_realk,  3.9210E0_realk,  3.9055E0_realk,  3.8951E0_realk&
     &,  3.8866E0_realk,  3.8689E0_realk,  3.9603E0_realk,  3.9109E0_realk,  3.9122E0_realk,  3.8233E0_realk,  3.7438E0_realk&
     &,  3.7436E0_realk,  3.6981E0_realk,  3.6555E0_realk,  3.5452E0_realk,  3.9327E0_realk,  4.0658E0_realk,  4.1175E0_realk&
     &,  2.9664E0_realk,  2.8209E0_realk,  3.5547E0_realk,  3.3796E0_realk,  3.3985E0_realk,  3.3164E0_realk,  3.2364E0_realk&
     &,  3.1956E0_realk,  3.0370E0_realk,  2.9313E0_realk,  3.6425E0_realk,  3.5565E0_realk,  3.7209E0_realk,  3.7108E0_realk&
     &/)
      R0AB(3501:3570)=(/&
     &   3.6639E0_realk,  3.6484E0_realk,  3.4745E0_realk,  3.3492E0_realk,  3.8755E0_realk,  4.2457E0_realk,  3.7758E0_realk&
     &,  3.7161E0_realk,  3.6693E0_realk,  3.6155E0_realk,  3.5941E0_realk,  3.5643E0_realk,  3.5292E0_realk,  3.4950E0_realk&
     &,  3.4720E0_realk,  3.4503E0_realk,  3.6936E0_realk,  3.7392E0_realk,  3.7388E0_realk,  3.7602E0_realk,  3.6078E0_realk&
     &,  3.4960E0_realk,  3.9800E0_realk,  4.3518E0_realk,  4.2802E0_realk,  3.8580E0_realk,  3.8056E0_realk,  3.7527E0_realk&
     &,  3.7019E0_realk,  3.6615E0_realk,  3.5768E0_realk,  3.5330E0_realk,  3.5038E0_realk,  3.5639E0_realk,  3.8192E0_realk&
     &,  3.8883E0_realk,  3.9092E0_realk,  3.9478E0_realk,  3.7995E0_realk,  3.6896E0_realk,  4.1165E0_realk,  4.5232E0_realk&
     &,  4.4357E0_realk,  4.4226E0_realk,  4.4031E0_realk,  4.3860E0_realk,  4.3721E0_realk,  4.3580E0_realk,  4.3466E0_realk&
     &,  4.2036E0_realk,  4.2037E0_realk,  3.8867E0_realk,  4.2895E0_realk,  4.2766E0_realk,  4.2662E0_realk,  4.2598E0_realk&
     &,  3.8408E0_realk,  3.9169E0_realk,  3.8681E0_realk,  3.8250E0_realk,  3.7855E0_realk,  3.7501E0_realk,  3.6753E0_realk&
     &,  3.5499E0_realk,  3.4872E0_realk,  3.5401E0_realk,  3.8288E0_realk,  3.9217E0_realk,  3.9538E0_realk,  4.0054E0_realk&
     &/)
      R0AB(3571:3640)=(/&
     &   2.8388E0_realk,  2.7890E0_realk,  3.4329E0_realk,  3.5593E0_realk,  3.3488E0_realk,  3.2486E0_realk,  3.1615E0_realk&
     &,  3.1000E0_realk,  3.0394E0_realk,  2.9165E0_realk,  3.5267E0_realk,  3.7479E0_realk,  3.6650E0_realk,  3.6263E0_realk&
     &,  3.5658E0_realk,  3.5224E0_realk,  3.4762E0_realk,  3.3342E0_realk,  3.7738E0_realk,  4.0333E0_realk,  3.9568E0_realk&
     &,  3.8975E0_realk,  3.8521E0_realk,  3.4929E0_realk,  3.7830E0_realk,  3.7409E0_realk,  3.7062E0_realk,  3.6786E0_realk&
     &,  3.6471E0_realk,  3.6208E0_realk,  3.6337E0_realk,  3.6519E0_realk,  3.6363E0_realk,  3.6278E0_realk,  3.6110E0_realk&
     &,  3.4825E0_realk,  3.8795E0_realk,  4.1448E0_realk,  4.0736E0_realk,  4.0045E0_realk,  3.6843E0_realk,  3.6291E0_realk&
     &,  3.5741E0_realk,  3.5312E0_realk,  3.4974E0_realk,  3.4472E0_realk,  3.4034E0_realk,  3.7131E0_realk,  3.7557E0_realk&
     &,  3.7966E0_realk,  3.8005E0_realk,  3.8068E0_realk,  3.8015E0_realk,  3.6747E0_realk,  4.0222E0_realk,  4.3207E0_realk&
     &,  4.2347E0_realk,  4.2191E0_realk,  4.1990E0_realk,  4.1811E0_realk,  4.1666E0_realk,  4.1521E0_realk,  4.1401E0_realk&
     &,  3.9970E0_realk,  3.9943E0_realk,  3.9592E0_realk,  4.0800E0_realk,  4.0664E0_realk,  4.0559E0_realk,  4.0488E0_realk&
     &/)
      R0AB(3641:3710)=(/&
     &   3.9882E0_realk,  4.0035E0_realk,  3.9539E0_realk,  3.9138E0_realk,  3.8798E0_realk,  3.8355E0_realk,  3.5359E0_realk&
     &,  3.4954E0_realk,  3.3962E0_realk,  3.5339E0_realk,  3.7595E0_realk,  3.8250E0_realk,  3.8408E0_realk,  3.8600E0_realk&
     &,  3.8644E0_realk,  2.7412E0_realk,  2.7489E0_realk,  3.3374E0_realk,  3.3950E0_realk,  3.3076E0_realk,  3.1910E0_realk&
     &,  3.0961E0_realk,  3.0175E0_realk,  3.0280E0_realk,  2.8929E0_realk,  3.4328E0_realk,  3.5883E0_realk,  3.6227E0_realk&
     &,  3.5616E0_realk,  3.4894E0_realk,  3.4241E0_realk,  3.3641E0_realk,  3.3120E0_realk,  3.6815E0_realk,  3.8789E0_realk&
     &,  3.8031E0_realk,  3.7413E0_realk,  3.6939E0_realk,  3.4010E0_realk,  3.6225E0_realk,  3.5797E0_realk,  3.5443E0_realk&
     &,  3.5139E0_realk,  3.4923E0_realk,  3.4642E0_realk,  3.5860E0_realk,  3.5849E0_realk,  3.5570E0_realk,  3.5257E0_realk&
     &,  3.4936E0_realk,  3.4628E0_realk,  3.7874E0_realk,  3.9916E0_realk,  3.9249E0_realk,  3.8530E0_realk,  3.5932E0_realk&
     &,  3.5355E0_realk,  3.4757E0_realk,  3.4306E0_realk,  3.3953E0_realk,  3.3646E0_realk,  3.3390E0_realk,  3.5637E0_realk&
     &,  3.7053E0_realk,  3.7266E0_realk,  3.7177E0_realk,  3.6996E0_realk,  3.6775E0_realk,  3.6558E0_realk,  3.9331E0_realk&
     &/)
      R0AB(3711:3780)=(/&
     &   4.1655E0_realk,  4.0879E0_realk,  4.0681E0_realk,  4.0479E0_realk,  4.0299E0_realk,  4.0152E0_realk,  4.0006E0_realk&
     &,  3.9883E0_realk,  3.8500E0_realk,  3.8359E0_realk,  3.8249E0_realk,  3.9269E0_realk,  3.9133E0_realk,  3.9025E0_realk&
     &,  3.8948E0_realk,  3.8422E0_realk,  3.8509E0_realk,  3.7990E0_realk,  3.7570E0_realk,  3.7219E0_realk,  3.6762E0_realk&
     &,  3.4260E0_realk,  3.3866E0_realk,  3.3425E0_realk,  3.5294E0_realk,  3.7022E0_realk,  3.7497E0_realk,  3.7542E0_realk&
     &,  3.7494E0_realk,  3.7370E0_realk,  3.7216E0_realk,  3.4155E0_realk,  3.0522E0_realk,  4.2541E0_realk,  3.8218E0_realk&
     &,  4.0438E0_realk,  3.5875E0_realk,  3.3286E0_realk,  3.1682E0_realk,  3.0566E0_realk,  2.9746E0_realk,  4.3627E0_realk&
     &,  4.0249E0_realk,  4.6947E0_realk,  4.1718E0_realk,  3.8639E0_realk,  3.6735E0_realk,  3.5435E0_realk,  3.4479E0_realk&
     &,  4.6806E0_realk,  4.3485E0_realk,  4.2668E0_realk,  4.1690E0_realk,  4.1061E0_realk,  4.1245E0_realk,  4.0206E0_realk&
     &,  3.9765E0_realk,  3.9458E0_realk,  3.9217E0_realk,  3.9075E0_realk,  3.8813E0_realk,  3.9947E0_realk,  4.1989E0_realk&
     &,  3.9507E0_realk,  3.7960E0_realk,  3.6925E0_realk,  3.6150E0_realk,  4.8535E0_realk,  4.5642E0_realk,  4.4134E0_realk&
     &/)
      R0AB(3781:3850)=(/&
     &   4.3688E0_realk,  4.3396E0_realk,  4.2879E0_realk,  4.2166E0_realk,  4.1888E0_realk,  4.1768E0_realk,  4.1660E0_realk&
     &,  4.1608E0_realk,  4.0745E0_realk,  4.2289E0_realk,  4.4863E0_realk,  4.2513E0_realk,  4.0897E0_realk,  3.9876E0_realk&
     &,  3.9061E0_realk,  5.0690E0_realk,  5.0446E0_realk,  4.6186E0_realk,  4.6078E0_realk,  4.5780E0_realk,  4.5538E0_realk&
     &,  4.5319E0_realk,  4.5101E0_realk,  4.4945E0_realk,  4.1912E0_realk,  4.2315E0_realk,  4.5534E0_realk,  4.4373E0_realk&
     &,  4.4224E0_realk,  4.4120E0_realk,  4.4040E0_realk,  4.2634E0_realk,  4.7770E0_realk,  4.6890E0_realk,  4.6107E0_realk&
     &,  4.5331E0_realk,  4.4496E0_realk,  4.4082E0_realk,  4.3095E0_realk,  4.2023E0_realk,  4.0501E0_realk,  4.2595E0_realk&
     &,  4.5497E0_realk,  4.3056E0_realk,  4.1506E0_realk,  4.0574E0_realk,  3.9725E0_realk,  5.0796E0_realk,  3.0548E0_realk&
     &,  3.3206E0_realk,  3.8132E0_realk,  3.9720E0_realk,  3.7675E0_realk,  3.7351E0_realk,  3.5167E0_realk,  3.5274E0_realk&
     &,  3.3085E0_realk,  3.1653E0_realk,  3.9500E0_realk,  4.1730E0_realk,  4.0613E0_realk,  4.1493E0_realk,  3.8823E0_realk&
     &,  4.0537E0_realk,  3.8200E0_realk,  3.6582E0_realk,  4.3422E0_realk,  4.5111E0_realk,  4.3795E0_realk,  4.3362E0_realk&
     &/)
      R0AB(3851:3920)=(/&
     &   4.2751E0_realk,  3.7103E0_realk,  4.1973E0_realk,  4.1385E0_realk,  4.1129E0_realk,  4.0800E0_realk,  4.0647E0_realk&
     &,  4.0308E0_realk,  4.0096E0_realk,  4.1619E0_realk,  3.9360E0_realk,  4.1766E0_realk,  3.9705E0_realk,  3.8262E0_realk&
     &,  4.5348E0_realk,  4.7025E0_realk,  4.5268E0_realk,  4.5076E0_realk,  3.9562E0_realk,  3.9065E0_realk,  3.8119E0_realk&
     &,  3.7605E0_realk,  3.7447E0_realk,  3.7119E0_realk,  3.6916E0_realk,  4.1950E0_realk,  4.2110E0_realk,  4.3843E0_realk&
     &,  4.1631E0_realk,  4.4427E0_realk,  4.2463E0_realk,  4.1054E0_realk,  4.7693E0_realk,  5.0649E0_realk,  4.7365E0_realk&
     &,  4.7761E0_realk,  4.7498E0_realk,  4.7272E0_realk,  4.7076E0_realk,  4.6877E0_realk,  4.6730E0_realk,  4.4274E0_realk&
     &,  4.5473E0_realk,  4.5169E0_realk,  4.5975E0_realk,  4.5793E0_realk,  4.5667E0_realk,  4.5559E0_realk,  4.3804E0_realk&
     &,  4.6920E0_realk,  4.6731E0_realk,  4.6142E0_realk,  4.5600E0_realk,  4.4801E0_realk,  4.0149E0_realk,  3.8856E0_realk&
     &,  3.7407E0_realk,  4.1545E0_realk,  4.2253E0_realk,  4.4229E0_realk,  4.1923E0_realk,  4.5022E0_realk,  4.3059E0_realk&
     &,  4.1591E0_realk,  4.7883E0_realk,  4.9294E0_realk,  3.3850E0_realk,  3.4208E0_realk,  3.7004E0_realk,  3.8800E0_realk&
     &/)
      R0AB(3921:3990)=(/&
     &   3.9886E0_realk,  3.9040E0_realk,  3.6719E0_realk,  3.6547E0_realk,  3.4625E0_realk,  3.3370E0_realk,  3.8394E0_realk&
     &,  4.0335E0_realk,  4.2373E0_realk,  4.3023E0_realk,  4.0306E0_realk,  4.1408E0_realk,  3.9297E0_realk,  3.7857E0_realk&
     &,  4.1907E0_realk,  4.3230E0_realk,  4.2664E0_realk,  4.2173E0_realk,  4.1482E0_realk,  3.6823E0_realk,  4.0711E0_realk&
     &,  4.0180E0_realk,  4.0017E0_realk,  3.9747E0_realk,  3.9634E0_realk,  3.9383E0_realk,  4.1993E0_realk,  4.3205E0_realk&
     &,  4.0821E0_realk,  4.2547E0_realk,  4.0659E0_realk,  3.9359E0_realk,  4.3952E0_realk,  4.5176E0_realk,  4.3888E0_realk&
     &,  4.3607E0_realk,  3.9583E0_realk,  3.9280E0_realk,  3.8390E0_realk,  3.7971E0_realk,  3.7955E0_realk,  3.7674E0_realk&
     &,  3.7521E0_realk,  4.1062E0_realk,  4.3633E0_realk,  4.2991E0_realk,  4.2767E0_realk,  4.4857E0_realk,  4.3039E0_realk&
     &,  4.1762E0_realk,  4.6197E0_realk,  4.8654E0_realk,  4.6633E0_realk,  4.5878E0_realk,  4.5640E0_realk,  4.5422E0_realk&
     &,  4.5231E0_realk,  4.5042E0_realk,  4.4901E0_realk,  4.3282E0_realk,  4.3978E0_realk,  4.3483E0_realk,  4.4202E0_realk&
     &,  4.4039E0_realk,  4.3926E0_realk,  4.3807E0_realk,  4.2649E0_realk,  4.6135E0_realk,  4.5605E0_realk,  4.5232E0_realk&
     &/)
      R0AB(3991:4060)=(/&
     &   4.4676E0_realk,  4.3948E0_realk,  4.0989E0_realk,  3.9864E0_realk,  3.8596E0_realk,  4.0942E0_realk,  4.2720E0_realk&
     &,  4.3270E0_realk,  4.3022E0_realk,  4.5410E0_realk,  4.3576E0_realk,  4.2235E0_realk,  4.6545E0_realk,  4.7447E0_realk&
     &,  4.7043E0_realk,  3.0942E0_realk,  3.2075E0_realk,  3.5152E0_realk,  3.6659E0_realk,  3.8289E0_realk,  3.7459E0_realk&
     &,  3.5156E0_realk,  3.5197E0_realk,  3.3290E0_realk,  3.2069E0_realk,  3.6702E0_realk,  3.8448E0_realk,  4.0340E0_realk&
     &,  3.9509E0_realk,  3.8585E0_realk,  3.9894E0_realk,  3.7787E0_realk,  3.6365E0_realk,  4.1425E0_realk,  4.1618E0_realk&
     &,  4.0940E0_realk,  4.0466E0_realk,  3.9941E0_realk,  3.5426E0_realk,  3.8952E0_realk,  3.8327E0_realk,  3.8126E0_realk&
     &,  3.7796E0_realk,  3.7635E0_realk,  3.7356E0_realk,  4.0047E0_realk,  3.9655E0_realk,  3.9116E0_realk,  4.1010E0_realk&
     &,  3.9102E0_realk,  3.7800E0_realk,  4.2964E0_realk,  4.3330E0_realk,  4.2622E0_realk,  4.2254E0_realk,  3.8195E0_realk&
     &,  3.7560E0_realk,  3.6513E0_realk,  3.5941E0_realk,  3.5810E0_realk,  3.5420E0_realk,  3.5178E0_realk,  3.8861E0_realk&
     &,  4.1459E0_realk,  4.1147E0_realk,  4.0772E0_realk,  4.3120E0_realk,  4.1207E0_realk,  3.9900E0_realk,  4.4733E0_realk&
     &/)
      R0AB(4061:4130)=(/&
     &   4.6157E0_realk,  4.4580E0_realk,  4.4194E0_realk,  4.3954E0_realk,  4.3739E0_realk,  4.3531E0_realk,  4.3343E0_realk&
     &,  4.3196E0_realk,  4.2140E0_realk,  4.2339E0_realk,  4.1738E0_realk,  4.2458E0_realk,  4.2278E0_realk,  4.2158E0_realk&
     &,  4.2039E0_realk,  4.1658E0_realk,  4.3595E0_realk,  4.2857E0_realk,  4.2444E0_realk,  4.1855E0_realk,  4.1122E0_realk&
     &,  3.7839E0_realk,  3.6879E0_realk,  3.5816E0_realk,  3.8633E0_realk,  4.1585E0_realk,  4.1402E0_realk,  4.1036E0_realk&
     &,  4.3694E0_realk,  4.1735E0_realk,  4.0368E0_realk,  4.5095E0_realk,  4.5538E0_realk,  4.5240E0_realk,  4.4252E0_realk&
     &,  3.0187E0_realk,  3.1918E0_realk,  3.5127E0_realk,  3.6875E0_realk,  3.7404E0_realk,  3.6943E0_realk,  3.4702E0_realk&
     &,  3.4888E0_realk,  3.2914E0_realk,  3.1643E0_realk,  3.6669E0_realk,  3.8724E0_realk,  3.9940E0_realk,  4.0816E0_realk&
     &,  3.8054E0_realk,  3.9661E0_realk,  3.7492E0_realk,  3.6024E0_realk,  4.0428E0_realk,  4.1951E0_realk,  4.1466E0_realk&
     &,  4.0515E0_realk,  4.0075E0_realk,  3.5020E0_realk,  3.9158E0_realk,  3.8546E0_realk,  3.8342E0_realk,  3.8008E0_realk&
     &,  3.7845E0_realk,  3.7549E0_realk,  3.9602E0_realk,  3.8872E0_realk,  3.8564E0_realk,  4.0793E0_realk,  3.8835E0_realk&
     &/)
      R0AB(4131:4200)=(/&
     &   3.7495E0_realk,  4.2213E0_realk,  4.3704E0_realk,  4.3300E0_realk,  4.2121E0_realk,  3.7643E0_realk,  3.7130E0_realk&
     &,  3.6144E0_realk,  3.5599E0_realk,  3.5474E0_realk,  3.5093E0_realk,  3.4853E0_realk,  3.9075E0_realk,  4.1115E0_realk&
     &,  4.0473E0_realk,  4.0318E0_realk,  4.2999E0_realk,  4.1050E0_realk,  3.9710E0_realk,  4.4320E0_realk,  4.6706E0_realk&
     &,  4.5273E0_realk,  4.4581E0_realk,  4.4332E0_realk,  4.4064E0_realk,  4.3873E0_realk,  4.3684E0_realk,  4.3537E0_realk&
     &,  4.2728E0_realk,  4.2549E0_realk,  4.2032E0_realk,  4.2794E0_realk,  4.2613E0_realk,  4.2491E0_realk,  4.2375E0_realk&
     &,  4.2322E0_realk,  4.3665E0_realk,  4.3061E0_realk,  4.2714E0_realk,  4.2155E0_realk,  4.1416E0_realk,  3.7660E0_realk&
     &,  3.6628E0_realk,  3.5476E0_realk,  3.8790E0_realk,  4.1233E0_realk,  4.0738E0_realk,  4.0575E0_realk,  4.3575E0_realk&
     &,  4.1586E0_realk,  4.0183E0_realk,  4.4593E0_realk,  4.5927E0_realk,  4.4865E0_realk,  4.3813E0_realk,  4.4594E0_realk&
     &,  2.9875E0_realk,  3.1674E0_realk,  3.4971E0_realk,  3.6715E0_realk,  3.7114E0_realk,  3.6692E0_realk,  3.4446E0_realk&
     &,  3.4676E0_realk,  3.2685E0_realk,  3.1405E0_realk,  3.6546E0_realk,  3.8579E0_realk,  3.9637E0_realk,  4.0581E0_realk&
     &/)
      R0AB(4201:4270)=(/&
     &   3.7796E0_realk,  3.9463E0_realk,  3.7275E0_realk,  3.5792E0_realk,  4.0295E0_realk,  4.1824E0_realk,  4.1247E0_realk&
     &,  4.0357E0_realk,  3.9926E0_realk,  3.4827E0_realk,  3.9007E0_realk,  3.8392E0_realk,  3.8191E0_realk,  3.7851E0_realk&
     &,  3.7687E0_realk,  3.7387E0_realk,  3.9290E0_realk,  3.8606E0_realk,  3.8306E0_realk,  4.0601E0_realk,  3.8625E0_realk&
     &,  3.7269E0_realk,  4.2062E0_realk,  4.3566E0_realk,  4.3022E0_realk,  4.1929E0_realk,  3.7401E0_realk,  3.6888E0_realk&
     &,  3.5900E0_realk,  3.5350E0_realk,  3.5226E0_realk,  3.4838E0_realk,  3.4594E0_realk,  3.8888E0_realk,  4.0813E0_realk&
     &,  4.0209E0_realk,  4.0059E0_realk,  4.2810E0_realk,  4.0843E0_realk,  3.9486E0_realk,  4.4162E0_realk,  4.6542E0_realk&
     &,  4.5005E0_realk,  4.4444E0_realk,  4.4196E0_realk,  4.3933E0_realk,  4.3741E0_realk,  4.3552E0_realk,  4.3406E0_realk&
     &,  4.2484E0_realk,  4.2413E0_realk,  4.1907E0_realk,  4.2656E0_realk,  4.2474E0_realk,  4.2352E0_realk,  4.2236E0_realk&
     &,  4.2068E0_realk,  4.3410E0_realk,  4.2817E0_realk,  4.2479E0_realk,  4.1921E0_realk,  4.1182E0_realk,  3.7346E0_realk&
     &,  3.6314E0_realk,  3.5168E0_realk,  3.8582E0_realk,  4.0927E0_realk,  4.0469E0_realk,  4.0313E0_realk,  4.3391E0_realk&
     &/)
      R0AB(4271:4340)=(/&
     &   4.1381E0_realk,  3.9962E0_realk,  4.4429E0_realk,  4.5787E0_realk,  4.4731E0_realk,  4.3588E0_realk,  4.4270E0_realk&
     &,  4.3957E0_realk,  2.9659E0_realk,  3.1442E0_realk,  3.4795E0_realk,  3.6503E0_realk,  3.6814E0_realk,  3.6476E0_realk&
     &,  3.4222E0_realk,  3.4491E0_realk,  3.2494E0_realk,  3.1209E0_realk,  3.6324E0_realk,  3.8375E0_realk,  3.9397E0_realk&
     &,  3.8311E0_realk,  3.7581E0_realk,  3.9274E0_realk,  3.7085E0_realk,  3.5598E0_realk,  4.0080E0_realk,  4.1641E0_realk&
     &,  4.1057E0_realk,  4.0158E0_realk,  3.9726E0_realk,  3.4667E0_realk,  3.8802E0_realk,  3.8188E0_realk,  3.7989E0_realk&
     &,  3.7644E0_realk,  3.7474E0_realk,  3.7173E0_realk,  3.9049E0_realk,  3.8424E0_realk,  3.8095E0_realk,  4.0412E0_realk&
     &,  3.8436E0_realk,  3.7077E0_realk,  4.1837E0_realk,  4.3366E0_realk,  4.2816E0_realk,  4.1686E0_realk,  3.7293E0_realk&
     &,  3.6709E0_realk,  3.5700E0_realk,  3.5153E0_realk,  3.5039E0_realk,  3.4684E0_realk,  3.4437E0_realk,  3.8663E0_realk&
     &,  4.0575E0_realk,  4.0020E0_realk,  3.9842E0_realk,  4.2612E0_realk,  4.0643E0_realk,  3.9285E0_realk,  4.3928E0_realk&
     &,  4.6308E0_realk,  4.4799E0_realk,  4.4244E0_realk,  4.3996E0_realk,  4.3737E0_realk,  4.3547E0_realk,  4.3358E0_realk&
     &/)
      R0AB(4341:4410)=(/&
     &   4.3212E0_realk,  4.2275E0_realk,  4.2216E0_realk,  4.1676E0_realk,  4.2465E0_realk,  4.2283E0_realk,  4.2161E0_realk&
     &,  4.2045E0_realk,  4.1841E0_realk,  4.3135E0_realk,  4.2562E0_realk,  4.2226E0_realk,  4.1667E0_realk,  4.0932E0_realk&
     &,  3.7134E0_realk,  3.6109E0_realk,  3.4962E0_realk,  3.8352E0_realk,  4.0688E0_realk,  4.0281E0_realk,  4.0099E0_realk&
     &,  4.3199E0_realk,  4.1188E0_realk,  3.9768E0_realk,  4.4192E0_realk,  4.5577E0_realk,  4.4516E0_realk,  4.3365E0_realk&
     &,  4.4058E0_realk,  4.3745E0_realk,  4.3539E0_realk,  2.8763E0_realk,  3.1294E0_realk,  3.5598E0_realk,  3.7465E0_realk&
     &,  3.5659E0_realk,  3.5816E0_realk,  3.3599E0_realk,  3.4024E0_realk,  3.1877E0_realk,  3.0484E0_realk,  3.7009E0_realk&
     &,  3.9451E0_realk,  3.8465E0_realk,  3.9873E0_realk,  3.7079E0_realk,  3.9083E0_realk,  3.6756E0_realk,  3.5150E0_realk&
     &,  4.0829E0_realk,  4.2780E0_realk,  4.1511E0_realk,  4.1260E0_realk,  4.0571E0_realk,  3.4865E0_realk,  3.9744E0_realk&
     &,  3.9150E0_realk,  3.8930E0_realk,  3.8578E0_realk,  3.8402E0_realk,  3.8073E0_realk,  3.7977E0_realk,  4.0036E0_realk&
     &,  3.7604E0_realk,  4.0288E0_realk,  3.8210E0_realk,  3.6757E0_realk,  4.2646E0_realk,  4.4558E0_realk,  4.2862E0_realk&
     &/)
      R0AB(4411:4465)=(/&
     &   4.2122E0_realk,  3.7088E0_realk,  3.6729E0_realk,  3.5800E0_realk,  3.5276E0_realk,  3.5165E0_realk,  3.4783E0_realk&
     &,  3.4539E0_realk,  3.9553E0_realk,  3.9818E0_realk,  4.2040E0_realk,  3.9604E0_realk,  4.2718E0_realk,  4.0689E0_realk&
     &,  3.9253E0_realk,  4.4869E0_realk,  4.7792E0_realk,  4.4918E0_realk,  4.5342E0_realk,  4.5090E0_realk,  4.4868E0_realk&
     &,  4.4680E0_realk,  4.4486E0_realk,  4.4341E0_realk,  4.2023E0_realk,  4.3122E0_realk,  4.2710E0_realk,  4.3587E0_realk&
     &,  4.3407E0_realk,  4.3281E0_realk,  4.3174E0_realk,  4.1499E0_realk,  4.3940E0_realk,  4.3895E0_realk,  4.3260E0_realk&
     &,  4.2725E0_realk,  4.1961E0_realk,  3.7361E0_realk,  3.6193E0_realk,  3.4916E0_realk,  3.9115E0_realk,  3.9914E0_realk&
     &,  3.9809E0_realk,  3.9866E0_realk,  4.3329E0_realk,  4.1276E0_realk,  3.9782E0_realk,  4.5097E0_realk,  4.6769E0_realk&
     &,  4.5158E0_realk,  4.3291E0_realk,  4.3609E0_realk,  4.3462E0_realk,  4.3265E0_realk,  4.4341E0_realk&
     &/)

      K=0
      DO I=1,N_ELEM
         DO J=1,I
            K=K+1
            R(I,J)=R0AB(K)*ANG_TO_AU
            R(J,I)=R0AB(K)*ANG_TO_AU
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SET_R0AB



      SUBROUTINE SET_R0(N_ELEM_D2,R00)
      IMPLICIT NONE
      !Pass in as R00 space for R0
      !Load from data statement 
      REAL(REALK)  :: R0(86)
      INTEGER      :: N_ELEM_D2,I
      REAL(REALK)  :: R00(N_ELEM_D2)
! the published radii in S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799 (tab 1)
! refer to the following values multiplied by 1.1 (rs6 in this code)
!                        1            2
      DATA      R0/ 0.91E0_realk,0.92E0_realk,&
!             3               4            5            6           7              8           9           10
     &      0.75E0_realk,1.28E0_realk,1.35E0_realk,1.32E0_realk,1.27E0_realk,1.22E0_realk,1.17E0_realk,1.13E0_realk,&
!             11             12           13           14           15           16           17          18
     &      1.04E0_realk,1.24E0_realk,1.49E0_realk,1.56E0_realk,1.55E0_realk,1.53E0_realk,1.49E0_realk,1.45E0_realk,&
!             19             20
     &      1.35E0_realk,1.34E0_realk,&
!             21              22          23            24          25
     &      1.42E0_realk,1.42E0_realk,1.42E0_realk,1.42E0_realk,1.42E0_realk,&
!             26             27          28            29           30
     &      1.42E0_realk,1.42E0_realk,1.42E0_realk,1.42E0_realk,1.42E0_realk,&
!             31             32           33           34           35          36
     &      1.50E0_realk,1.57E0_realk,1.60E0_realk,1.61E0_realk,1.59E0_realk,1.57E0_realk,&
!             37             38
     &      1.48E0_realk,1.46E0_realk,&
!             39              40          41           42           43
     &      1.49E0_realk,1.49E0_realk,1.49E0_realk,1.49E0_realk,1.49E0_realk,&
!             44             45            46          47          48
     &      1.49E0_realk,1.49E0_realk,1.49E0_realk,1.49E0_realk,1.49E0_realk,&
!             49               50            51            52            53     54
     &      1.52E0_realk,1.64E0_realk,1.71E0_realk,1.72E0_realk,1.72E0_realk,1.71E0_realk,&
!             55               56            57            58            59            60           61   
     &      1.638E0_realk,1.602E0_realk,1.564E0_realk,1.594E0_realk,1.594E0_realk,1.594E0_realk,1.594E0_realk,&
!             62               63            64            65            66            67           68   
     &      1.594E0_realk,1.594E0_realk,1.594E0_realk,1.594E0_realk,1.594E0_realk,1.594E0_realk,1.594E0_realk,&
!             69               70            71
     &      1.594E0_realk,1.594E0_realk,1.594E0_realk,&
!             72               73            74            75           76            77            78   
     &      1.625E0_realk,1.611E0_realk,1.611E0_realk,1.611E0_realk,1.611E0_realk,1.611E0_realk,1.611E0_realk,&
!             79
     &      1.611E0_realk,&
!             80              81             82            83            84            85           86   
     &      1.598E0_realk,1.805E0_realk,1.767E0_realk,1.725E0_realk,1.823E0_realk,1.810E0_realk,1.749E0_realk/
      
      DO I = 1,N_ELEM_D2
        R00(I) = R0(I)
      ENDDO

      RETURN
      END SUBROUTINE SET_R0


      SUBROUTINE SET_C6(N_ELEM_D2,C600)
      IMPLICIT NONE 
      !Pass in as C600 space for C6
      !Load from data statement 
      REAL(REALK)   :: C6(86)
      INTEGER       :: N_ELEM_D2,I
      REAL(REALK)   :: C600(N_ELEM_D2)
!                      1            2
      DATA      C6/0.14E0_realk,0.08E0_realk,&
!          3               4           5            6            7            8            9           10
     &   1.61E0_realk,1.61E0_realk,3.13E0_realk,1.75E0_realk,1.23E0_realk,0.70E0_realk,0.75E0_realk,0.63E0_realk,&
!          11             12            13           14            15           16          17           18
     &   5.71E0_realk,5.71E0_realk,10.79E0_realk,9.23E0_realk,7.84E0_realk,5.57E0_realk,5.07E0_realk,4.61E0_realk,&
!          19             20            21          22           23
     &   10.8E0_realk,10.8E0_realk,10.8E0_realk,10.8E0_realk,10.8E0_realk,&
!          24               25        26           27             28          29          30             31
     &   10.8E0_realk,10.8E0_realk,10.8E0_realk,10.8E0_realk,10.8E0_realk,10.8E0_realk,10.8E0_realk,16.99E0_realk,&
!          32               33            34              35             36             37            38    
     &   17.10E0_realk,16.37E0_realk,12.64E0_realk,12.47E0_realk,12.01E0_realk,24.67E0_realk,24.67E0_realk,&
!          39               40            41              42             43             44            45  
     &   24.67E0_realk,24.67E0_realk,24.67E0_realk,24.67E0_realk,24.67E0_realk,24.67E0_realk,24.67E0_realk,&
!          46               47            48              49             50             51            52
     &   24.67E0_realk,24.67E0_realk,24.67E0_realk,37.32E0_realk,38.71E0_realk,38.44E0_realk,31.74E0_realk,&
!          53               54            55              56             57
     &   31.50E0_realk,29.99E0_realk,315.275E0_realk,226.994E0_realk,176.252E0_realk,&
!          58               59            60              61             62             63            64
     &  140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,&
!          65               66            67              68             69             70            71
     &  140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,140.68E0_realk,&
!          72
     &  105.112E0_realk,&
!          73               74            75              76             77             78            79
     &  81.24E0_realk,81.24E0_realk,81.24E0_realk,81.24E0_realk,81.24E0_realk,81.24E0_realk,81.24E0_realk,&
!          80               81            82              83             84             85            86
     &  57.364E0_realk,57.254E0_realk,63.162E0_realk,63.540E0_realk,55.283E0_realk,57.171E0_realk,56.64E0_realk /

      DO I = 1,N_ELEM_D2
        C600(I) = C6(I)
      ENDDO

      RETURN
      END SUBROUTINE SET_C6


      SUBROUTINE SET_RCOV(N_ELEM,RCOV00)
      IMPLICIT NONE 
      !Pass in as RCOV00 space for RCOV
      !Load from data statement 
      REAL(REALK)   ::  RCOV(94)
      INTEGER       ::  N_ELEM,I
      REAL(REALK)   ::  RCOV00(N_ELEM)
!
!     Covalent radii, original values Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197
!     Here scaled by k2 (=4/3) and converted to a_0 (factor = 0.52917726E0_realk)
!     Elements H --> Pu indexed by Z.
!
      DATA RCOV/&
!          1                       2                  3                  4                  5
     & 0.80628308E0_realk, 1.15903197E0_realk, 3.02356173E0_realk, 2.36845659E0_realk, 1.94011865E0_realk,&
!          6                       7                  8                  9                  10
     & 1.88972601E0_realk, 1.78894056E0_realk, 1.58736983E0_realk, 1.61256616E0_realk, 1.68815527E0_realk,&
!          11                      12                 13                 14                 15
     & 3.52748848E0_realk, 3.14954334E0_realk, 2.84718717E0_realk, 2.62041997E0_realk, 2.77159820E0_realk,&
!          16                      17                 18                 19                 20
     & 2.57002732E0_realk, 2.49443835E0_realk, 2.41884923E0_realk, 4.43455700E0_realk, 3.88023730E0_realk,&
!          21                      22                 23                 24                 25
     & 3.35111422E0_realk, 3.07395437E0_realk, 3.04875805E0_realk, 2.77159820E0_realk, 2.69600923E0_realk,&
!          26                      27                 28                 29                 30
     & 2.62041997E0_realk, 2.51963467E0_realk, 2.49443835E0_realk, 2.54483100E0_realk, 2.74640188E0_realk,&
!          31                      32                 33                 34                 35
     & 2.82199085E0_realk, 2.74640188E0_realk, 2.89757982E0_realk, 2.77159820E0_realk, 2.87238349E0_realk,&
!          36                      37                 38                 39                 40
     & 2.94797246E0_realk, 4.76210950E0_realk, 4.20778980E0_realk, 3.70386304E0_realk, 3.50229216E0_realk,&
!          41                      42                 43                 44                 45
     & 3.32591790E0_realk, 3.12434702E0_realk, 2.89757982E0_realk, 2.84718717E0_realk, 2.84718717E0_realk,&
!          46                      47                 48                 49                 50
     & 2.72120556E0_realk, 2.89757982E0_realk, 3.09915070E0_realk, 3.22513231E0_realk, 3.17473967E0_realk,&
!          51                      52                 53                 54                 55
     & 3.17473967E0_realk, 3.09915070E0_realk, 3.32591790E0_realk, 3.30072128E0_realk, 5.26603625E0_realk,&
!          56                      57                 58                 59                 60
     & 4.43455700E0_realk, 4.08180818E0_realk, 3.70386304E0_realk, 3.98102289E0_realk, 3.95582657E0_realk,&
!          61                      62                 63                 64                 65
     & 3.93062995E0_realk, 3.90543362E0_realk, 3.80464833E0_realk, 3.82984466E0_realk, 3.80464833E0_realk,&
!          66                      67                 68                 69                 70
     & 3.77945201E0_realk, 3.75425569E0_realk, 3.75425569E0_realk, 3.72905937E0_realk, 3.85504098E0_realk,&
!          71                      72                 73                 74                 75
     & 3.67866672E0_realk, 3.45189952E0_realk, 3.30072128E0_realk, 3.09915070E0_realk, 2.97316878E0_realk,&
!          76                      77                 78                 79                 80
     & 2.92277614E0_realk, 2.79679452E0_realk, 2.82199085E0_realk, 2.84718717E0_realk, 3.32591790E0_realk,&
!          81                      82                 83                 84                 85
     & 3.27552496E0_realk, 3.27552496E0_realk, 3.42670319E0_realk, 3.30072128E0_realk, 3.47709584E0_realk,&
!          86                      87                 88                 89                 90
     & 3.57788113E0_realk, 5.06446567E0_realk, 4.56053862E0_realk, 4.20778980E0_realk, 3.98102289E0_realk,&
!          91                      92                 93                 94           
     & 3.82984466E0_realk, 3.85504098E0_realk, 3.88023730E0_realk, 3.90543362E0_realk /

      DO I = 1, N_ELEM
         RCOV00(I)=RCOV(I)
      ENDDO

      RETURN
      END SUBROUTINE SET_RCOV


      SUBROUTINE SET_R2R4(N_ELEM,R2R400)
      IMPLICIT NONE 
      !Pass in as R2R400 space for R2R4
      !Load from data statement 
      REAL(REALK)  :: R2R4(94)
      INTEGER      :: N_ELEM,I
      REAL(REALK)  :: R2R400(N_ELEM)
      DATA R2R4 /&
!          1                         2                     3                    4                   5
     & 2.00734898E0_realk,  1.56637132E0_realk,  5.01986934E0_realk,  3.85379032E0_realk,  3.64446594E0_realk,&
!          6                         7                     8                    9                   10
     & 3.10492822E0_realk,  2.71175247E0_realk,  2.59361680E0_realk,  2.38825250E0_realk,  2.21522516E0_realk,&
!          11                        12                    13                   14                  15
     & 6.58585536E0_realk,  5.46295967E0_realk,  5.65216669E0_realk,  4.88284902E0_realk,  4.29727576E0_realk,&
!          16                        17                    18                   19                  20
     & 4.04108902E0_realk,  3.72932356E0_realk,  3.44677275E0_realk,  7.97762753E0_realk,  7.07623947E0_realk,&
!          21                        22                    23                   24                  25
     & 6.60844053E0_realk,  6.28791364E0_realk,  6.07728703E0_realk,  5.54643096E0_realk,  5.80491167E0_realk,&
!          26                        27                    28                   29                  30
     & 5.58415602E0_realk,  5.41374528E0_realk,  5.28497229E0_realk,  5.22592821E0_realk,  5.09817141E0_realk,&
!          31                        32                    33                   34                  35
     & 6.12149689E0_realk,  5.54083734E0_realk,  5.06696878E0_realk,  4.87005108E0_realk,  4.59089647E0_realk,&
!          36                        37                    38                   39                  40
     & 4.31176304E0_realk,  9.55461698E0_realk,  8.67396077E0_realk,  7.97210197E0_realk,  7.43439917E0_realk,&
!          41                        42                    43                   44                  45
     & 6.58711862E0_realk,  6.19536215E0_realk,  6.01517290E0_realk,  5.81623410E0_realk,  5.65710424E0_realk,&
!          46                        47                    48                   49                  50
     & 5.52640661E0_realk,  5.44263305E0_realk,  5.58285373E0_realk,  7.02081898E0_realk,  6.46815523E0_realk,&
!          51                        52                    53                   54                  55
     & 5.98089120E0_realk,  5.81686657E0_realk,  5.53321815E0_realk,  5.25477007E0_realk, 11.02204549E0_realk,&
!          56                        57                    58                   59                  60
     &10.15679528E0_realk,  9.35167836E0_realk,  9.06926079E0_realk,  8.97241155E0_realk,  8.90092807E0_realk,&
!          61                        62                    63                   64                  65
     & 8.85984840E0_realk,  8.81736827E0_realk,  8.79317710E0_realk,  7.89969626E0_realk,  8.80588454E0_realk,&
!          66                        67                    68                   69                  70
     & 8.42439218E0_realk,  8.54289262E0_realk,  8.47583370E0_realk,  8.45090888E0_realk,  8.47339339E0_realk,&
!          71                        72                    73                   74                  75
     & 7.83525634E0_realk,  8.20702843E0_realk,  7.70559063E0_realk,  7.32755997E0_realk,  7.03887381E0_realk,&
!          76                        77                    78                   79                  80
     & 6.68978720E0_realk,  6.05450052E0_realk,  5.88752022E0_realk,  5.70661499E0_realk,  5.78450695E0_realk,&
!          81                        82                    83                   84                  85
     & 7.79780729E0_realk,  7.26443867E0_realk,  6.78151984E0_realk,  6.67883169E0_realk,  6.39024318E0_realk,&
!          86                        87                    88                   89                  90
     & 6.09527958E0_realk, 11.79156076E0_realk, 11.10997644E0_realk,  9.51377795E0_realk,  8.67197068E0_realk,&
!          91                        92                    93                   94           
     & 8.77140725E0_realk,  8.65402716E0_realk,  8.53923501E0_realk,  8.85024712E0_realk /

      DO I = 1,N_ELEM
        R2R400(I)=R2R4(I)
      ENDDO

      RETURN
      END SUBROUTINE SET_R2R4

      END MODULE IIDFTD

