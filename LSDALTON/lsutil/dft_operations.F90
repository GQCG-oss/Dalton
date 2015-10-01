!> dft type module
!> \author T.Kjaergaard
!> \date 2010-02-21
MODULE dft_type
use dft_typetype
use precision
use dft_memory_handling, only: mem_dft_alloc, mem_dft_dealloc
use lsmpi_type, only: LS_MPI_BUFFER
use infpar_module
CONTAINS
subroutine DFT_set_default_config(DFT)
implicit none
type(DFTparam) :: DFT
integer :: i
DFT%IGRID = Grid_Default
DFT%RADIALGRID = 3   !(default TURBO) (1=GC2,2=LMG,3=TURBO) (olddefault LMG)
DFT%PARTITIONING = 5 !(default BLOCKSSF)                    (olddefault becke-original)
!(1=SSF, 2=Becke, 3=Becke-original, 4=block, 5=blockssf, 6=cartesian)
DFT%ZdependenMaxAng = .TRUE.  !default on                   (olddefault off)
DFT%GRIDDONE = 0 !FALSE = GRID NOT DONE YET  
DFT%RADINT = 5.01187E-14_realk                                   !(olddefault 1.0E-11_realk)
DFT%ANGMIN = 5                                             
DFT%ANGINT = 35                                            !(olddefault 31)
DFT%HRDNES = 3 ! hardness of the becke partioning function
DFT%DFTHR0 = 1.0E-9_realk  !not used
DFT%DFTHRI = 2.0E-12_realk
DFT%DFTELS = 1.0E-3_realk
DFT%DFTHRL = 2.0E-10_realk !not used - obsolete
DFT%RHOTHR = 2.0E-15_realk
DFT%NOPRUN = .FALSE. 
DFT%DFTASC = .FALSE.
DFT%DFTPOT = .FALSE.
DFT%DODISP = .FALSE.
DFT%DOORBFREE = .FALSE.
!AMT New Dispersion Defaults
DFT%DO_DFTD2        =       .FALSE.
DFT%L_INP_D2PAR     =       .FALSE.
DFT%DODISP2         =       .FALSE.
DFT%DODISP3         =       .FALSE.
DFT%DO_DFTD3        =       .FALSE.
DFT%DO_BJDAMP       =       .FALSE.
DFT%DO_3BODY        =       .FALSE.
DFT%L_INP_D3PAR     =       .FALSE.
!AMT
DFT%DFTIPT = 1.0E-20_realk      
DFT%DFTBR1 = 1.0E-20_realk
DFT%DFTBR2 = 1.0E-20_realk
DFT%DFTADD = .TRUE.  !disable DFT D
DFT%DISPDONE = .FALSE.
DFT%TURBO = 1 !TURBOMOLE grid                               (olddefault 0)
DFT%NBUFLEN=0
do i=1,80
   DFT%dftfunc(i:i)=' '
enddo
DFT%testNelectrons=.TRUE.                                         !(olddefault FALSE)
DFT%CS00=.FALSE.
DFT%LB94=.FALSE.
DFT%CS00shift=0.0E0_realk
DFT%CS00eHOMO=0.0E0_realk
DFT%CS00ZND1=0.2332E0_realk
DFT%CS00ZND2=0.0116E0_realk !(=0.315eV)
DFT%HFexchangeFac=0.0E0_realk
DFT%XCFUN=.TRUE.
call init_gridObject(dft,DFT%gridObject)
call init_DFTfunc(dft)
end subroutine DFT_set_default_config

subroutine dft_setIntegralSchemeFromInput(schemeDFT,dalton_inpDFT)
type(DFTparam) :: schemeDFT,dalton_inpDFT

schemeDFT=dalton_inpDFT

end subroutine dft_setIntegralSchemeFromInput

subroutine WRITE_FORMATTET_DFT_PARAM(LUPRI,DFT)
implicit none
integer :: LUPRI
type(DFTparam) :: DFT

SELECT CASE(DFT%RADIALGRID)
   case(1); WRITE(LUPRI,'(2X,A35)')'Radial grid: GC2'
   case(2); WRITE(LUPRI,'(2X,A35)')'Radial grid: LMG'
   case(3); WRITE(LUPRI,'(2X,A35)')'Radial grid: TURBO'
CASE DEFAULT
   WRITE (LUPRI,*)'RADIAL GRID TYPE: ',DFT%RADIALGRID,' not recognized'
   CALL lsQUIT('Illegal value of DFT%RADIALGRID',lupri)
END SELECT
SELECT CASE(DFT%PARTITIONING)
   case(1); WRITE(LUPRI,'(2X,A35)')'Using SSF partitioning'
   case(2); WRITE(LUPRI,'(2X,A35)')'Using BECKE partitioning'
   case(3); WRITE(LUPRI,'(2X,A35)')'Using BECKE-ORIGINAL partitioning'
   case(4); WRITE(LUPRI,'(2X,A35)')'Using BLOCK partitioning'
   case(5); WRITE(LUPRI,'(2X,A35)')'Using BLOCKSSF partitioning'
   case(6); WRITE(LUPRI,'(2X,A35)')'Using CARTESIAN partitioning'
CASE DEFAULT
   WRITE (LUPRI,*)'PARTITIONING TYPE: ',DFT%PARTITIONING,' not recognized'
   CALL lsQUIT('Illegal value of DFT%PARTITIONING',lupri)
END SELECT
WRITE(LUPRI,'(2X,A35,L1)') 'ZdependenMaxAng', DFT%ZdependenMaxAng
WRITE(LUPRI,'(2X,A35,I8)')'GRDONE',DFT%GRIDDONE 
WRITE(LUPRI,'(2X,A35,F16.8)') 'RADINT',DFT%RADINT
WRITE(LUPRI,'(2X,A35,I8)')'ANGMIN',DFT%ANGMIN
WRITE(LUPRI,'(2X,A35,I8)')'ANGINT',DFT%ANGINT
WRITE(LUPRI,'(2X,A35,I8)')'Hardness',DFT%HRDNES
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTELS',DFT%DFTELS
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTHR0',DFT%DFTHR0
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTHRI',DFT%DFTHRI
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTHRL',DFT%DFTHRL
WRITE(LUPRI,'(2X,A35,F16.8)') 'RHOTHR',DFT%RHOTHR
WRITE(LUPRI,'(2X,A35,L1)') 'NOPRUN', DFT%NOPRUN
WRITE(LUPRI,'(2X,A35,L1)') 'DFTASC',DFT%DFTASC
WRITE(LUPRI,'(2X,A35,L1)') 'DFTPOT',DFT%DFTPOT
WRITE(LUPRI,'(2X,A35,L1)') 'DISPER', DFT%DODISP
!AMT New Dispersion values
WRITE(LUPRI,'(2X,A35,L1)') 'DO_DFTD2', DFT%DO_DFTD2        
WRITE(LUPRI,'(2X,A35,L1)') 'L_INP_D2PAR', DFT%L_INP_D2PAR     
WRITE(LUPRI,'(2X,A35,L1)') 'DODISP2', DFT%DODISP2         
WRITE(LUPRI,'(2X,A35,L1)') 'DODISP3', DFT%DODISP3         
WRITE(LUPRI,'(2X,A35,L1)') 'DO_DFTD3', DFT%DO_DFTD3        
WRITE(LUPRI,'(2X,A35,L1)') 'DO_BJDAMP', DFT%DO_BJDAMP       
WRITE(LUPRI,'(2X,A35,L1)') 'DO_3BODY', DFT%DO_3BODY        
WRITE(LUPRI,'(2X,A35,L1)') 'L_INP_D3PAR', DFT%L_INP_D3PAR     
!AMT
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTIPT',DFT%DFTIPT
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTBR1',DFT%DFTBR1
WRITE(LUPRI,'(2X,A35,F16.8)') 'DFTBR2',DFT%DFTBR2
WRITE(LUPRI,'(2X,A35,L1)') 'DFTADD', DFT%DFTADD
WRITE(LUPRI,'(2X,A35,L1)') 'DISPDONE',DFT%DISPDONE
WRITE(LUPRI,'(2X,A35,I8)')'TURBO',DFT%TURBO
WRITE(LUPRI,'(2X,A35,I8)')'NBUFLEN',DFT%NBUFLEN
WRITE(LUPRI,'(2X,A75,A)')'dftfunc',DFT%dftfunc
WRITE(LUPRI,'(2X,A35,L1)')'testNelectrons',DFT%testNelectrons
WRITE(LUPRI,'(2X,A35,L1)')'LB94',DFT%LB94
WRITE(LUPRI,'(2X,A35,L1)')'CS00',DFT%CS00
WRITE(LUPRI,'(2X,A35,F16.8)') 'CS00shift',DFT%CS00shift
WRITE(LUPRI,'(2X,A35,F16.8)') 'CS00eHOMO',DFT%CS00eHOMO
WRITE(LUPRI,'(2X,A35,F16.8)') 'CS00ZND1',DFT%CS00ZND1
WRITE(LUPRI,'(2X,A35,F16.8)') 'CS00ZND2',DFT%CS00ZND2
WRITE(LUPRI,'(2X,A35,F16.8)') 'HF exchange Factor',DFT%HFexchangeFac
WRITE(LUPRI,'(2X,A35,L1)')'XCFUN',DFT%XCFUN
end subroutine WRITE_FORMATTET_DFT_PARAM

!> \brief mpi copy the DFTdata structure
!> \author T. Kjaergaard
!> \date 2011
!> \param DFTdata
SUBROUTINE free_DFTdata(DFTdata)
implicit none
type(DFTDATATYPE) :: DFTdata
!
IF(associated(DFTdata%energy))THEN
   call mem_dft_dealloc(DFTDATA%energy)
   nullify(DFTdata%energy)
endif

IF (associated(DFTdata%BMAT))THEN
   call mem_dft_dealloc(DFTDATA%BMAT)
   nullify(DFTdata%BMAT)
ENDIF

IF (associated(DFTdata%FKSM))THEN
   call mem_dft_dealloc(DFTDATA%FKSM)
   nullify(DFTdata%FKSM)
endif

IF (associated(DFTdata%FKSMS))THEN
   call mem_dft_dealloc(DFTDATA%FKSMS)
   nullify(DFTdata%FKSMS)
endif

IF(associated(DFTdata%orb2atom))THEN
   call mem_dft_dealloc(DFTDATA%orb2atom)
   nullify(DFTDATA%orb2atom)   
ENDIF

IF (associated(DFTdata%grad))THEN
   call mem_dft_dealloc(DFTDATA%grad)
   nullify(DFTDATA%grad)   
ENDIF
end subroutine free_DFTdata


subroutine mpicopy_DFTparam(DFT,master)
implicit none
integer(kind=ls_mpik) :: master
type(DFTparam) :: DFT
integer :: i
call LS_MPI_BUFFER(DFT%IGRID,Master)
call LS_MPI_BUFFER(DFT%IDFTtype,Master)
call LS_MPI_BUFFER(DFT%RADIALGRID,Master)
call LS_MPI_BUFFER(DFT%PARTITIONING,Master)
call LS_MPI_BUFFER(DFT%ZdependenMaxAng,Master)
call LS_MPI_BUFFER(DFT%GRIDDONE,Master)
call LS_MPI_BUFFER(DFT%RADINT,Master)
call LS_MPI_BUFFER(DFT%ANGMIN,Master)
call LS_MPI_BUFFER(DFT%ANGINT,Master)
call LS_MPI_BUFFER(DFT%HRDNES,Master)
call LS_MPI_BUFFER(DFT%DFTELS,Master)
call LS_MPI_BUFFER(DFT%DFTHR0,Master)
call LS_MPI_BUFFER(DFT%DFTHRI,Master)
call LS_MPI_BUFFER(DFT%DFTHRL,Master)
call LS_MPI_BUFFER(DFT%RHOTHR,Master)
call LS_MPI_BUFFER(DFT%NOPRUN,Master)
call LS_MPI_BUFFER(DFT%DFTASC,Master)
call LS_MPI_BUFFER(DFT%DFTPOT,Master)
call LS_MPI_BUFFER(DFT%DODISP,Master)
!SR&AB Orbital-free DFT param
call LS_MPI_BUFFER(DFT%DOORBFREE,Master)
IF(DFT%DOORBFREE)THEN
   call LS_MPI_BUFFER(DFT%ORBFREE%numberTSfunc,Master)
   DO i=1,DFT%ORBFREE%numberTSfunc
      call LS_MPI_BUFFER(DFT%ORBFREE%TScoeff(i),Master)
      call LS_MPI_BUFFER(DFT%ORBFREE%TSfunc(i),Master)
   ENDDO
   call LS_MPI_BUFFER(DFT%ORBFREE%kineticFac,Master)
ENDIF
!AMT New Dispersion values
call LS_MPI_BUFFER(DFT%DO_DFTD2,Master)
call LS_MPI_BUFFER(DFT%L_INP_D2PAR,Master)
call LS_MPI_BUFFER(DFT%D2_s6_inp,Master)  
call LS_MPI_BUFFER(DFT%D2_alp_inp,Master)
call LS_MPI_BUFFER(DFT%D2_rs6_inp,Master)
call LS_MPI_BUFFER(DFT%DODISP2,Master)
call LS_MPI_BUFFER(DFT%DODISP3,Master)
call LS_MPI_BUFFER(DFT%DO_DFTD3,Master)
call LS_MPI_BUFFER(DFT%DO_BJDAMP,Master)
call LS_MPI_BUFFER(DFT%DO_3BODY,Master)
call LS_MPI_BUFFER(DFT%L_INP_D3PAR,Master)
call LS_MPI_BUFFER(DFT%D3_s6_inp,Master)
call LS_MPI_BUFFER(DFT%D3_alp_inp,Master)
call LS_MPI_BUFFER(DFT%D3_rs6_inp,Master)
call LS_MPI_BUFFER(DFT%D3_rs18_inp,Master)
call LS_MPI_BUFFER(DFT%D3_s18_inp,Master)
!AMT
call LS_MPI_BUFFER(DFT%DFTIPT,Master)
call LS_MPI_BUFFER(DFT%DFTBR1,Master)
call LS_MPI_BUFFER(DFT%DFTBR2,Master)
call LS_MPI_BUFFER(DFT%DFTADD,Master)
call LS_MPI_BUFFER(DFT%DISPDONE,Master)
call LS_MPI_BUFFER(DFT%TURBO,Master)
call LS_MPI_BUFFER(DFT%NBUFLEN,Master)
call LS_MPI_BUFFER(DFT%dftfunc,Master)
call LS_MPI_BUFFER(DFT%testNelectrons,Master)
call LS_MPI_BUFFER(DFT%LB94,Master)
call LS_MPI_BUFFER(DFT%CS00,Master)
call LS_MPI_BUFFER(DFT%CS00shift,Master)
call LS_MPI_BUFFER(DFT%CS00eHOMO,Master)
call LS_MPI_BUFFER(DFT%CS00ZND1,Master)
call LS_MPI_BUFFER(DFT%CS00ZND2,Master)
call LS_MPI_BUFFER(DFT%HFexchangeFac,Master)
call LS_MPI_BUFFER(DFT%XCFUN,Master)
do i=1,size(DFT%dftfuncObject)
   call LS_MPI_BUFFER(DFT%dftfuncObject(i),Master)
enddo
do i=1,size(DFT%GridObject)
   call LS_MPI_BUFFER(DFT%GridObject(i)%RADIALGRID,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%PARTITIONING,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%ZdependenMaxAng,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%GRIDDONE,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%RADINT,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%ANGMIN,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%ANGINT,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%HRDNES,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%NOPRUN,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%TURBO,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%NBUFLEN,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%Id,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%NBAST,Master)
   call LS_MPI_BUFFER(DFT%GridObject(i)%Numnodes,Master)
enddo

end subroutine mpicopy_DFTparam

subroutine initDFTdatatype(DFTdata)
implicit none
type(DFTDATATYPE) :: DFTdata

nullify(DFTdata%BMAT)
nullify(DFTdata%FKSM)
nullify(DFTdata%FKSMS)
nullify(DFTdata%orb2atom)
nullify(DFTdata%grad)
nullify(DFTdata%Energy)
DFTDATA%LB94=.FALSE.
DFTDATA%CS00=.FALSE.
DFTDATA%CS00shift=0.0E0_realk
DFTDATA%CS00eHOMO=0.0E0_realk
DFTDATA%CS00ZND1=0.0E0_realk
DFTDATA%CS00ZND2=0.0E0_realk
DFTDATA%HFexchangeFac=0.0E0_realk
DFTDATA%nWorkNactBastNblen = 0
DFTDATA%nWorkNactBast  = 0
DFTDATA%nWorkNactBastNactBast  = 0
end subroutine initDFTdatatype

subroutine copyDFTdata(newDFTdata,DFTdata)
implicit none
type(DFTDATATYPE) :: newDFTdata,DFTdata
!
integer :: nbast,ndmat,nbmat,natoms,I,nfmat

nbast = DFTdata%nbast
nfmat = DFTdata%nfmat
ndmat = DFTdata%ndmat
nbmat = DFTdata%nbmat
natoms= DFTdata%natoms
newDFTdata%nbast = nbast
newDFTdata%ndmat = ndmat
newDFTdata%nfmat = nfmat
newDFTdata%nbmat = nbmat
newDFTdata%nWorkNactBastNblen = DFTdata%nWorkNactBastNblen
newDFTdata%nWorkNactBast  = DFTdata%nWorkNactBast
newDFTdata%nWorkNactBastNactBast  = DFTdata%nWorkNactBastNactBast
newDFTdata%LB94 = DFTdata%LB94
newDFTdata%CS00 = DFTdata%CS00
newDFTdata%CS00shift = DFTdata%CS00shift
newDFTdata%CS00eHOMO = DFTdata%CS00eHOMO
newDFTdata%CS00ZND1 = DFTdata%CS00ZND1
newDFTdata%CS00ZND2 = DFTdata%CS00ZND2
newDFTdata%HFexchangeFac = DFTdata%HFexchangeFac
if(associated(DFTdata%Energy))then
   call mem_dft_alloc(newDFTDATA%Energy,ndmat)
   newDFTDATA%Energy=0.0E0_realk
else
   nullify(newDFTdata%energy)
endif

!THE BMAT, FKSM and FKSMS are NOT copied because this routine
!is used to distribute data to a private variable and these should be SHARED 

!if(associated(DFTdata%BMAT))then
!   call mem_dft_alloc(newDFTDATA%BMAT,nbast,nbast,nbmat)
!   CALL DCOPY(nbast*nbast*nbmat,DFTDATA%BMAT,1,newDFTDATA%BMAT,1)
!else
   nullify(newDFTdata%BMAT)
!endif
!if(associated(DFTdata%FKSM))then
!   call mem_dft_alloc(newDFTDATA%FKSM,nbast,nbast,nfmat)
!   CALL LS_DZERO(newDFTDATA%FKSM,nbast*nbast*nfmat)
!else
   nullify(newDFTdata%FKSM)
!endif
!if(associated(DFTdata%FKSMS))then
!   call mem_dft_alloc(newDFTDATA%FKSMS,nbast,nbast,nfmat)
!   CALL LS_DZERO(newDFTDATA%FKSMS,nbast*nbast*nfmat)
!else
   nullify(newDFTdata%FKSMS)
!endif
newDFTdata%dosympart = DFTdata%dosympart
newDFTdata%natoms = natoms

if(associated(DFTdata%orb2atom))then
   call mem_dft_alloc(newDFTDATA%orb2atom,nbast)
   do I=1,nbast
      newDFTDATA%orb2atom(I) = DFTDATA%orb2atom(I)
   enddo
else
   nullify(newDFTDATA%orb2atom)   
endif

if(associated(DFTdata%grad))then
   call mem_dft_alloc(newDFTDATA%grad,3,natoms)   
   CALL LS_DZERO(newDFTDATA%grad,3*natoms)
else
   nullify(newDFTDATA%grad)   
endif

end subroutine copyDFTdata

subroutine DFTdataReduction(inputDFTdata,collectDFTdata)
implicit none
type(DFTDATATYPE) :: inputDFTdata,collectDFTdata
!
integer :: nbast,ndmat,natoms,nfmat,idmat

nbast= inputDFTdata%nbast
ndmat= inputDFTdata%ndmat
nfmat= inputDFTdata%nfmat
natoms= inputDFTdata%natoms

if(associated(inputDFTdata%energy))then
   DO idmat=1,ndmat 
      collectDFTdata%Energy(idmat)=collectDFTdata%Energy(idmat)+inputDFTdata%Energy(idmat)
   ENDDO
   call mem_dft_dealloc(inputDFTDATA%energy)
endif
nullify(inputDFTDATA%energy)   

if(associated(inputDFTdata%FKSM))then 
   !collect result
   CALL DAXPY(nbast*nbast*nfmat,1E0_realk,inputDFTDATA%FKSM,1,collectDFTdata%FKSM,1)
   call mem_dft_dealloc(inputDFTDATA%FKSM)
endif
nullify(inputDFTdata%FKSM)

if(associated(inputDFTdata%FKSMS))then
   !collect result
   CALL DAXPY(nbast*nbast*nfmat,1E0_realk,inputDFTDATA%FKSMS,1,collectDFTdata%FKSMS,1)
   call mem_dft_dealloc(inputDFTDATA%FKSMS)
endif
nullify(inputDFTdata%FKSMS)

if(associated(inputDFTdata%grad))then
   !collect result
   CALL DAXPY(3*natoms,1E0_realk,inputDFTDATA%grad,1,collectDFTdata%grad,1)
   call mem_dft_dealloc(inputDFTDATA%grad)
endif
nullify(inputDFTDATA%grad)   

if(associated(inputDFTdata%BMAT))then
   call mem_dft_dealloc(inputDFTDATA%BMAT)
endif
nullify(inputDFTdata%BMAT)

if(associated(inputDFTdata%orb2atom))then
   call mem_dft_dealloc(inputDFTDATA%orb2atom)   
endif
nullify(inputDFTDATA%orb2atom)   

end subroutine DFTdataReduction

#ifdef VAR_MPI
!> \brief mpi copy the DFTdata structure
!> \author T. Kjaergaard
!> \date 2011
!> \param DFTdata
SUBROUTINE mpicopy_DFTdata(DFTdata,mynum)
implicit none
type(DFTDATATYPE) :: DFTdata
integer(kind=ls_mpik):: mynum
!
integer(kind=ls_mpik) :: master
integer :: nbast,ndmat,nbmat,natoms,I,nfmat
Logical :: TESTassociatedenergy,TESTassociatedBMAT,TESTassociatedFMAT,TESTassociatedFMATS
Logical :: TESTassociatedorb2atom,TESTassociatedgrad,SLAVE
Master= infpar%master
SLAVE = mynum.ne.Master
IF(SLAVE)THEN
   call initDFTdatatype(DFTdata)
ENDIF
!INTEGER
CALL LS_MPI_BUFFER(DFTdata%nbast,Master)
CALL LS_MPI_BUFFER(DFTdata%nfmat,Master)
CALL LS_MPI_BUFFER(DFTdata%ndmat,Master)
CALL LS_MPI_BUFFER(DFTdata%nbmat,Master)
CALL LS_MPI_BUFFER(DFTdata%natoms,Master)
CALL LS_MPI_BUFFER(DFTdata%nWorkNactBastNblen,Master)
CALL LS_MPI_BUFFER(DFTdata%nWorkNactBast,Master)
CALL LS_MPI_BUFFER(DFTdata%nWorkNactBastNactBast,Master)

CALL LS_MPI_BUFFER(DFTdata%LB94,Master)
CALL LS_MPI_BUFFER(DFTdata%CS00,Master)
CALL LS_MPI_BUFFER(DFTdata%CS00shift,Master)
CALL LS_MPI_BUFFER(DFTdata%CS00eHOMO,Master)
CALL LS_MPI_BUFFER(DFTdata%CS00ZND1,Master)
CALL LS_MPI_BUFFER(DFTdata%CS00ZND2,Master)
CALL LS_MPI_BUFFER(DFTdata%HFexchangeFac,Master)

nbast = DFTdata%nbast
nfmat = DFTdata%nfmat
ndmat = DFTdata%ndmat
nbmat = DFTdata%nbmat
natoms= DFTdata%natoms
!DFTdata%electrons = 0.0E0_realk
CALL LS_MPI_BUFFER(DFTdata%dosympart,Master)

! BROADCAST DFTdata%energy
IF (.NOT.SLAVE) THEN
   TESTassociatedenergy = associated(DFTdata%energy)
ENDIF
CALL LS_MPI_BUFFER(TESTassociatedenergy,Master)

if(TESTassociatedenergy)then
   IF (SLAVE) call mem_dft_alloc(DFTDATA%energy,ndmat)
   CALL LS_MPI_BUFFER(DFTDATA%energy,ndmat,Master)
else
   IF (SLAVE)nullify(DFTdata%energy)
endif

! BROADCAST DFTdata%BMAT
IF (.NOT.SLAVE) THEN
   TESTassociatedBMAT = associated(DFTdata%BMAT)
ENDIF
CALL LS_MPI_BUFFER(TESTassociatedBmat,Master)

if(TESTassociatedBmat)then
   IF (SLAVE) call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   CALL LS_MPI_BUFFER(DFTDATA%BMAT,nbast,nbast,nbmat,Master)
else
   IF (SLAVE)nullify(DFTdata%BMAT)
endif

! BROADCAST DFTdata%FKSM
IF (.NOT.SLAVE) THEN
   TESTassociatedFMAT = associated(DFTdata%FKSM)
ENDIF
CALL LS_MPI_BUFFER(TESTassociatedFmat,Master)

if(TESTassociatedFmat)then
   IF (SLAVE)THEN
      call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,nfmat)
      CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*nfmat)
   ENDIF
else
   IF (SLAVE)nullify(DFTdata%FKSM)
endif

! BROADCAST DFTdata%FKSMS
IF (.NOT.SLAVE) THEN
   TESTassociatedFMATS = associated(DFTdata%FKSMS)
ENDIF
CALL LS_MPI_BUFFER(TESTassociatedFmatS,Master)

if(TESTassociatedFmatS)then
   IF (SLAVE)THEN
      call mem_dft_alloc(DFTDATA%FKSMS,nbast,nbast,nfmat)
      CALL LS_DZERO(DFTDATA%FKSMS,nbast*nbast*nfmat)
   ENDIF
else
   IF(SLAVE)nullify(DFTdata%FKSMS)
endif

! BROADCAST DFTdata%orb2atom
IF (.NOT.SLAVE) THEN
   TESTassociatedOrb2atom = associated(DFTdata%orb2atom)
ENDIF
CALL LS_MPI_BUFFER(TESTassociatedOrb2atom,Master)

IF(TESTassociatedOrb2atom)THEN
   IF(SLAVE)call mem_dft_alloc(DFTDATA%orb2atom,nbast)
   CALL LS_MPI_BUFFER(DFTDATA%orb2atom,nbast,Master)
else
   IF(SLAVE)nullify(DFTDATA%orb2atom)   
endif

! BROADCAST DFTdata%grad
IF (.NOT.SLAVE) THEN
   TESTassociatedgrad = associated(DFTdata%grad)
ENDIF
CALL LS_MPI_BUFFER(TESTassociatedgrad,Master)

if(TESTassociatedgrad)then
   IF(SLAVE)THEN
      call mem_dft_alloc(DFTDATA%grad,3,natoms)   
      CALL LS_DZERO(DFTDATA%grad,3*natoms)
   ENDIF
else
   IF(SLAVE)nullify(DFTDATA%grad)   
endif

end subroutine mpicopy_DFTdata

#endif

END MODULE dft_type
