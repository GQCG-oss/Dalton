
module fde_evaluators_dirac


! embed modules
   use fde_types
   use fde_cfg
   use fde_data
   use fde_io
!  use fde_nadd_derv
!  use fde_xcfun_interface
!  use fde_max_block_length

! dirac-specific modules
   use numoper_integrals
   use fde_dirac_matrices_integration
   use electrostatic_potential
   use density_eval
   use ao_eval
   use interface_ao
   use interface_mo_specific

! these are the interfaces to the dirac-specific code
   public fde_calculate_elpot
   public fde_dirac_set_nz
   public fde_dirac_get_nz
   public fde_dirac_get_isymop
   public fde_dirac_set_isymop
   public fde_dirac_get_ihrmop
   public fde_dirac_set_ihrmop

! below we have the routines that actually calculate
   public fde_get_density
   interface fde_get_density
      module procedure fde_dirac_get_density_from_dfcoef
   end interface

   public fde_get_density_matrix

   public fde_qccode_data_interface
   interface fde_qccode_data_interface
      module procedure fde_dirac_data_interface
   end interface

! the matrix for the embedding potential, or its contributions to
! different fock matrices can be calculated in different ways; 
! we take this into account by defining an interface name to the
! code-specific routines
 
   public fde_calculate_emb_pot_mat
   interface fde_calculate_emb_pot_mat
      module procedure fde_dirac_embpot_via_integrator
      module procedure fde_dirac_embpot_via_aoproper
   end interface


   public fde_calculate_emb_linrsp_mat
   interface fde_calculate_emb_linrsp_mat
      module procedure fde_dirac_linrsp_via_integrator
   end interface


   private 

   integer, save :: nba, nba_orb, nba_aux
   integer, save :: dfcoef_unit
   integer, save :: nz = 4
   integer, parameter :: fde_max_nr_mat = 10
   integer, save :: isymop(fde_max_nr_mat)
   integer, save :: ihrmop(fde_max_nr_mat)

   logical, save :: is_dirac_initialized = .false.

   contains

!   ----------------------------------------------------------------------------
      subroutine fde_calculate_elpot(level,grid,gf)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(inout) :: gf
         type(fde_grid), intent(in)         :: grid
         character(len=4)                   :: level
         real(kind=8), allocatable          :: dmat(:)

! dmat is allocated and set in the routine below
         call fde_dirac_get_dmat_from_dfcoef(dmat,level)
! we calculate both the nuclear and electrostatic potentials in one call
! if needed, we can split this up later 
         call fde_dirac_calculate_elpot(dmat,grid%r,gf%elpot,gf%nucpot)

         deallocate(dmat)

      end subroutine fde_calculate_elpot

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_calculate_elpot(dmat,r,vc,vnuc)
!   ----------------------------------------------------------------------------
         real(kind=8), intent(in)  :: dmat(:)
         real(kind=8), intent(in)  :: r(:,:)
         real(kind=8), intent(out) :: vc(:)
         real(kind=8), intent(out) :: vnuc(:)
         integer :: i
         integer :: np
         integer :: irrep = 1
         logical :: nuc_part = .false.    
         logical :: ele_part = .true.    
         
         np = size(r,2)
!
! calculating the hartree potential first
! 
         write (*,*) 'Generating the Hartree potential over the grid'
!
! aspg, 11/10/2012
! important here: get_esp expects to get a density matrix
!                 that is not scaled wih the occupation numbers
!                 (so like what one would get from calling DENMAT).
!                 
         nuc_part = .false.    
         ele_part = .true.    
         call get_esp(np,            &
                      vc,            &
                      irrep,         & 
                      size(dmat, 1), & 
                      dmat,          &
                      r,             &
                      nuc_part,      &
                      ele_part)
!
! now the nuclear potential
!
         write (*,*) 'Generating the nuclear potential over the grid'

         nuc_part = .true.    
         ele_part = .false.    
         call get_esp(np,            &
                      vnuc,          &
                      irrep,         & 
                      size(dmat, 1), & 
                      dmat,          &
                      r,             &
                      nuc_part,      &
                      ele_part)
!
! assemble the final electrostatic potential
!
         vc   =  vnuc + vc 

      end subroutine fde_dirac_calculate_elpot
      
#ifdef PRG_DIRAC
!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_dmat_from_dfcoef(dmat,level)
! the outgoing dmat is not scaled with occupation numbers
!   ----------------------------------------------------------------------------
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "mxcent.h"
#include "nuclei.h"
! Passed/recovered variables
      character(len=4) :: level
      real(kind=8), allocatable, intent(out) :: dmat(:)
! internal variables
      INTEGER :: ierr
      real(kind=8), allocatable :: cmo(:)

      allocate(dmat(n2bbasxq))
      dmat = 0.0d0
      allocate(cmo(ncmotq))
      cmo = 0.0d0

      select case(level)
      case('MP2')
         write(lupri,*) 'Reading CC_DENSITY'
         call get_cc_density(level,'AO',DMAT,ierr)
         if (ierr.ne.0) then
            write(lupri,*) 'Error in reading density matrix'
            call fde_quit('error readin density matrix')
         end if
      case default
!     HF density used
         write(lupri,*) 'Reading DFCOEF'
         call opnfil(LUCOEF,'DFCOEF','UNKNOWN','EMBDRV')
         call reacmo_no_work(lucoef,    &
                             'DFCOEF',  &
                             cmo,       &
                             (/0.0d0/), &
                             (/0/),     &
                             dummy,     &
                             2)
!     SCF construct density matrix
      call denmat(dmat,cmo,0)
      end select

      deallocate(cmo)

      end subroutine fde_dirac_get_dmat_from_dfcoef 


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_density_from_dfcoef(level,grid,gf)
!   ----------------------------------------------------------------------------
      type(grid_function), intent(inout) :: gf
      type(fde_grid), intent(in)         :: grid
      character(len=4), intent(in)       :: level

      INTEGER :: i
      INTEGER :: NDIMENSION
      INTEGER :: NDER=0, DOLND=0

      REAL(KIND=8),ALLOCATABLE :: GAO(:), GAO1(:),GAO2(:),BUF(:),NCNT(:)
      REAL(KIND=8),ALLOCATABLE :: GAB1(:)

      real(kind=8), allocatable :: dmat(:)

      REAL(KIND=8)::dummy,TOLS,TOLOG
      REAL(KIND=8)::gexp=0.0,ncentc=1.0,factor=1.0

      NDIMENSION=nba_orb

!     DMAT does into account the occupation numbers for scf and mp2 ...
      call fde_dirac_get_dmat_from_dfcoef(dmat,level)
!     so we correct it for scf and mp2...
      dmat = 2.0d0*dmat

      allocate(gao (   ndimension))
      allocate(gao1( 3*ndimension))
      allocate(gao2( 6*ndimension))
      allocate(gab1( 3*ndimension))
      allocate(ncnt(   ndimension))
      allocate(buf ( 4*ndimension))

!     Initialize the matrices to zero.

      GAO  = 0.0
      GAO1 = 0.0
      GAO2 = 0.0
      BUF  = 0.0
      NCNT = 0.0
      GAB1 = 0.0

      WRITE(*,*) 'Calculating the density and its derivatives.'
      
      do i = 1, size(grid%r,2)

        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,  &
                grid%r(1,i), grid%r(2,i), grid%r(3,i), &
                    BUF,NDIMENSION,2,DOLND,0)
        CALL GETRHO(gf%n(i),0,GAO,DMAT,BUF)
        CALL GETGRHO(gf%gn(1:3,i),0,GAO,GAO1,DMAT,BUF)
        CALL GETLRHO(gf%hn(1,i),gf%hn(1,i),gf%hn(1,i), &
                     gf%hn(4,i),gf%hn(5,i),gf%hn(6,i), &
                     GAO,GAO1,GAO2,DMAT,BUF)

      END DO

      deallocate (buf)
      deallocate (gao)
      deallocate (gao1)
      deallocate (gao2)
      deallocate (gab1)
      deallocate (ncnt) 
      deallocate (dmat)

      end subroutine fde_dirac_get_density_from_dfcoef



!   ----------------------------------------------------------------------------
      subroutine fde_dirac_embpot_via_aoproper(comb,iprint)
!   ----------------------------------------------------------------------------
!     the following sets up the matrix elements for ground-state embedding 
!     potential, as a so-called "numerical operator" (that is, we have
!     values of the given operator on a set of gridpoints P, and contract
!     that with the value of the orbitals (or their derivatives, whatever)
!     and the integration weight to get the final integral: 
!
!     M_ab = < phi^n_a | operator | phi^m_b > = \sum^P_i w_i phi^n_a(i) operator(i) phi^m_b(i) 
!
!     where i is the i-th grid point, w_i the integration weight, and
!     phi^{n,m}_{a,b}(i) is the value of the {n,m}-th order derivative of
!     the (contracted) basis function {a,b} at point i.
!
!     information about the ao-basis blocks active for a given operator is
!     carried in the variable COMB, as indicated up. In the case of the
!     embedding potential for FDE, we have a diagonal operator just like
!     for the molecular field.
#include "implicit.h"
#include "priunit.h"
            LOGICAL DOINT(4)
            CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
            character*80 string
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"

            integer :: file_unit, iprint, i, irrep_vemb
            real(kind=8), allocatable :: a_oneint(:)
            character(len=60) :: file_name
            type(fde_files)   :: ftmp

            call fde_get_files_info(ftmp)
      
            file_unit = ftmp%embpot%unit
            file_name = ftmp%embpot%name
            call fde_open_file(file_name,file_unit)

 
            label = 'FDEVEMB '
            irrep_vemb = 0  ! embedding potential is totally symmetric

            read(comb,'(4L1)',err=1000) (doint(i),i=1,4)

            allocate(a_oneint(nnbbasx)) 
            a_oneint = 0.0d0
            call NumOper_OneElOpMatrix(a_oneint,file_unit,doint,irrep_vemb)

            IF(IPRINT.GE.4) THEN
              write (string,*) 'Integrals over (numerical) operator: '//label
              call HEADER(string,-1)
              call OUTPAK(a_oneint,NTBAS(0),1,6)
            ENDIF
!
!        Generate integral labels
!
            CALL GETDAT(RTNLBL(1),RTNLBL(2))
            RTNLBL(2)(1:2) = 'SY'
            WRITE(RTNLBL(2)(3:4),'(I2)') 1
            RTNLBL(2)(5:8) = COMB
!
!        Write integrals to file
!
            call WRTPRO(a_oneint,NNBBASX,LABEL,RTNLBL,IPRINT)
!
            call fde_close_file(file_unit)
            deallocate(a_oneint)

            return
 1000       continue
!
!     Not able to read DOINT information from COMB
!
            write(LUPRI,'(A,A)') 'FDE_SaveAOPROPER_StaticEmbPot: '// &
                                 'Not able to read COMB =',COMB
            call fde_quit('Vemb2AOPROPER: Not able to read COMB')

         end subroutine fde_dirac_embpot_via_aoproper

!   ----------------------------------------------------------------------------
!     subroutine fde_dirac_embpot_via_integrator(mat_dim,dmat,fmat,in_nz)
      subroutine fde_dirac_embpot_via_integrator(mat_dim,dmat,fmat)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: mat_dim
         real(kind=8), target, intent(in) :: dmat(*)
         real(kind=8), target  :: fmat(*)

         call fde_dirac_emb_matrices_via_integration(      &
                                  fde_mat_dim   = mat_dim, &
                                  fde_nz        = nz,      &
                                  fde_dmat_0    = dmat,    &
                                  fde_nr_dmat   = 0,       &
                                  fde_nr_fmat   = 1,       &
                                  fde_fmat      = fmat,    &
                                  fde_do_potential = .true.)

      end subroutine fde_dirac_embpot_via_integrator

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_linrsp_via_integrator(mat_dim,dmat0,ndmat,dmat,fmat)
!   ----------------------------------------------------------------------------
         integer, intent(in)  :: mat_dim
         integer, intent(in)  :: ndmat
         real(kind=8), target, intent(in) :: dmat0(*)
         real(kind=8), target :: dmat(*)
         real(kind=8), target :: fmat(*)

         call fde_dirac_emb_matrices_via_integration(       &
                           fde_mat_dim           = mat_dim, &
                           fde_nz                = nz,      &
                           fde_dmat_0            = dmat0,   &
                           fde_nr_dmat           = ndmat,   &
                           fde_nr_fmat           = ndmat,   &
                           fde_dmat              = dmat,    &
                           fde_fmat              = fmat,    &
                           fde_fmat_pg_sym       = isymop,  &
                           fde_dmat_pg_sym       = isymop,  &
                           fde_dmat_ih_sym       = ihrmop,  &
                           fde_response_order_mo = 1)

      end subroutine fde_dirac_linrsp_via_integrator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routines to set/recover values for certain dirac-specific variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_set_isymop(value)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: value(fde_max_nr_mat)
         isymop = value
      end subroutine fde_dirac_set_isymop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_set_ihrmop(value)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: value(fde_max_nr_mat)
         ihrmop = value
      end subroutine fde_dirac_set_ihrmop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_set_nz(value)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: value
         nz = value
      end subroutine fde_dirac_set_nz


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_isymop(value)
!   ----------------------------------------------------------------------------
         integer, intent(out) :: value(fde_max_nr_mat)
         value = isymop
      end subroutine fde_dirac_get_isymop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_ihrmop(value)
!   ----------------------------------------------------------------------------
         integer, intent(out) :: value(fde_max_nr_mat)
         value = ihrmop
      end subroutine fde_dirac_get_ihrmop

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_get_nz(value)
!   ----------------------------------------------------------------------------
        integer, intent(out) :: value
         value = nz
      end subroutine fde_dirac_get_nz


!   ----------------------------------------------------------------------------
      subroutine fde_dirac_data_interface
!   ----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#ifdef MOD_DNF
#include "densfit.h"
#endif
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"

         nba_orb = ntbas(0)
         nba = nba_orb
         dfcoef_unit   = LUCOEF
         mat_dim_quat_as_1d = n2bbasxq

      end subroutine fde_dirac_data_interface

end module fde_evaluators_dirac
