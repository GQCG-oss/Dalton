***********************************************************************

      SUBROUTINE GASCI_old(ISM,ISPC,IPRNT,EREF,IPROCLIST,
     &                     IGROUPLIST,fh_array)
*
* CI optimization in GAS space number ISPC for symmetry ISM
*
*
* Jeppe Olsen, Winter of 1995
*
#ifdef VAR_MPI
      use file_io_model
#endif
      use lucita_energy_types
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLCLBT, KLCLEBT, KLCI1BT, KLCIBT, KLC2B,
     &          KLCBLTP, KBSCLFC, KBLCKIF, KBLCKPD,
     &          KRCCTOS, 
     &          KILUCLIST, KILU1LIST, KILU2LIST, KILU3LIST, KILU4LIST, 
     &          KILU5LIST, KILU6LIST, KILU7LIST, KLBASSPC,
     &          KLEBASC, KLCBASC, KLSPSPCL, KLBLKCLS, KLCLSL, KLCLSC, 
     &          KLCLSE, KLCLSCT, KLCLSET, KLCLSA, KLCLSA2, KLBLKA, 
     &          KLCLSD, KLCLSDT, KVEC1, KVEC2
!               for addressing of WORK
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "comjep.inc"
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) ITEST_OFF
      integer                                   :: file_offset_off
      integer,                      allocatable :: file_offset_fac_i4(:)
      integer(KIND=MPI_OFFSET_KIND),allocatable :: file_offset_fac(:)
      integer(KIND=MPI_OFFSET_KIND),allocatable :: file_offset_array(:)
#endif
#include "parluci.h"
! Common block for communicating with sigma
#include "cands.inc"
!     common block needed for CSF run (variable NCSASM)
      COMMON/SPINFO/MULTSP,MS2P,
     &              MINOP,MAXOP,NTYP,NDPCNT(MXPCTP),NCPCNT(MXPCTP),
     &              NCNATS(MXPCTP,MXPCSM),NDTASM(MXPCSM),NCSASM(MXPCSM),
     &              NCNASM(MXPCSM)
*
      DIMENSION IPROCLIST(LUCI_NMPROC)
      DIMENSION IGROUPLIST(LUCI_NMPROC)
      integer fh_array(nr_files)
      real(8),    allocatable :: tmp_block_scaling_fac(:)
!     general information needed in parallel runs
      integer,    allocatable :: par_dist_block_list(:)
      integer,    allocatable :: block_list(:)
      integer                 :: number_of_dets = 0
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GASCI ')
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
*. Initialize tests on comjeppe.inc
C     MXACJ,MXACIJ,MXAADST
      MXACJ    = 0
      MXACIJ   = 0
      MXAADST  = 0
!     transfer to CANDS
      ICSM           = ISM
      ISSM           = ISM
      ICSPC          = ISPC
      ISSPC          = ISPC

      number_of_dets = XISPSM(ICSM,ICSPC)
      LBLK           = - 1
      IF(ICISTR.EQ.1) LBLK = number_of_dets

      IF(NOCSF.EQ.0)THEN
!       number_of_dets = NCSASM(ICSM)
        call quit('CSF version of LUCITA not active.')
      END  IF
      
      WRITE(LUWRT,'(a,i14)') 
     &     '  @@ Number of determinants/combinations: ',
     &           number_of_dets
      if(number_of_dets .eq. 0)then
      WRITE(LUWRT,*) ' the number of determinants/combinations is zero.'
      WRITE(LUWRT,*) ' I am sure that fascinating discussions about '
      WRITE(LUWRT,*) ' the energy of such a wave function exists, '
      WRITE(LUWRT,*) ' but I am just a dumb program, so I will stop'
      WRITE(LUWRT,*)
      WRITE(LUWRT,*) ' GASCI : Vanishing number of parameters '
      call quit(' GASCI : Vanishing number of parameters ')
      end if
 
!     allocate two CI vector blocks and one resolution vector block
      call get_3_blocks(kvec1,kvec2,kvec3)
!
!     information about block structure ==> memory for partitioning of C vector
      iatp   = 1
      ibtp   = 2
      call z_blkfo(icspc,icsm,iatp,ibtp,klclbt,klclebt,
     &             klci1bt,klcibt,klcbltp,nbatch,nblock)

!     parallel stuff starts here
!     --------------------------
      allocate(block_list(nblock))
      call icopy(nblock,work(klci1bt),1,block_list,1)

C     array for block distribution among cpu's (parallel run)
      allocate(par_dist_block_list(nblock))
      par_dist_block_list = -2

!     allocate array for final sigma-cvec connection storage
      CALL MEMMAN(KRCCTOS,NBLOCK,'ADDS  ',1,'RCCTOS')
      CALL IZERO(WORK(KRCCTOS),NBLOCK)
 
#ifdef VAR_MPI
      IF (LUCI_NMPROC .GT. 1) THEN

      allocate(tmp_block_scaling_fac(nblock))

!     find optimal block distribution among cpu's
!     -------------------------------------------
!     block_list          == list of all blocks containing their length
!     par_dist_block_list == list of all blocks containing their assigned CPU
      CALL BLOCK_DISTR_DRV(NBLOCK,block_list,par_dist_block_list,
     &                     WORK(KRCCTOS),tmp_block_scaling_fac,
     &                     number_of_dets,IPROCLIST)

      deallocate(tmp_block_scaling_fac)

!     ----------------------------------------------------------
!     organize MPI file I/O offsets with the following ordering: 
!     idia, iluc, ilu[2-7], ilu1
!     ----------------------------------------------------------
      allocate(file_offset_fac_i4(nr_files))
      allocate(file_offset_fac(nr_files))
      allocate(file_offset_array(nr_files))
      file_offset_array = 0
      file_offset_fac   = 0
      file_offset_off   = 0

      file_offset_fac_i4(1) = 1
      file_offset_fac_i4(2) = 0
      file_offset_fac_i4(3) = mxciv + nroot
      file_offset_fac_i4(4) = mxciv + nroot
      file_offset_fac_i4(5) = mxciv + nroot
      file_offset_fac_i4(6) = mxciv + nroot
      file_offset_fac_i4(7) = 1
      file_offset_fac_i4(8) = mxciv
      file_offset_fac_i4(9) = nroot

      file_offset_fac(1) = file_offset_fac_i4(1)
      file_offset_fac(2) = file_offset_fac_i4(2)
      file_offset_fac(3) = file_offset_fac_i4(3)
      file_offset_fac(4) = file_offset_fac_i4(4)
      file_offset_fac(5) = file_offset_fac_i4(5)
      file_offset_fac(6) = file_offset_fac_i4(6)
      file_offset_fac(7) = file_offset_fac_i4(7)
      file_offset_fac(8) = file_offset_fac_i4(8)
      file_offset_fac(9) = file_offset_fac_i4(9)
      
      call set_file_io_offset(nr_files,
     &                        file_offset_off,
     &                        file_offset_array,
     &                        file_offset_fac,
     &                        my_vec1_ioff,
     &                        my_vec2_ioff,
     &                        my_act_blk1,
     &                        my_act_blk2,
     &                        my_act_blk_all,
     &                        .false.,
     &                        luci_myproc,
     &                        num_blocks2,
     &                        newcomm_proc,
     &                        igrouplist,
     &                        par_dist_block_list,
     &                        block_list)

!     save output in common block variables
      my_dia_off = file_offset_array(1)
      my_luc_off = file_offset_array(2)
      my_lu2_off = file_offset_array(3)
      my_lu3_off = file_offset_array(4)
      my_lu4_off = file_offset_array(5)
      my_lu5_off = file_offset_array(6)
      my_lu6_off = file_offset_array(7)
      my_lu7_off = file_offset_array(8)
      my_lu1_off = file_offset_array(9)

      deallocate(file_offset_array)
      deallocate(file_offset_fac)

!     length for allocation of file arrays
      iall_luc =         1             * num_blocks2
      iall_lu2 = file_offset_fac_i4(3) * my_act_blk2
      iall_lu3 = file_offset_fac_i4(4) * my_act_blk2
      iall_lu4 = file_offset_fac_i4(5) * my_act_blk2
      iall_lu5 = file_offset_fac_i4(6) * my_act_blk2
      iall_lu6 = file_offset_fac_i4(7) * my_act_blk2
      iall_lu7 = file_offset_fac_i4(8) * my_act_blk2
      iall_lu1 = file_offset_fac_i4(9) * my_act_blk2

      deallocate(file_offset_fac_i4)

!     allocate file arrays
      CALL MEMMAN(KILUCLIST,IALL_LUC,'ADDS  ',1,'LUCLST')
      CALL MEMMAN(KILU1LIST,IALL_LU1,'ADDS  ',1,'LU1LST')
      CALL MEMMAN(KILU2LIST,IALL_LU2,'ADDS  ',1,'LU2LST')
      CALL MEMMAN(KILU3LIST,IALL_LU3,'ADDS  ',1,'LU3LST')
      CALL MEMMAN(KILU4LIST,IALL_LU4,'ADDS  ',1,'LU4LST')
      CALL MEMMAN(KILU5LIST,IALL_LU5,'ADDS  ',1,'LU5LST')
      CALL MEMMAN(KILU6LIST,IALL_LU6,'ADDS  ',1,'LU6LST')
      CALL MEMMAN(KILU7LIST,IALL_LU7,'ADDS  ',1,'LU7LST')
!     initialize ...
      CALL IZERO(WORK(KILUCLIST),IALL_LUC)
      CALL IZERO(WORK(KILU1LIST),IALL_LU1)
      CALL IZERO(WORK(KILU2LIST),IALL_LU2)
      CALL IZERO(WORK(KILU3LIST),IALL_LU3)
      CALL IZERO(WORK(KILU4LIST),IALL_LU4)
      CALL IZERO(WORK(KILU5LIST),IALL_LU5)
      CALL IZERO(WORK(KILU6LIST),IALL_LU6)
      CALL IZERO(WORK(KILU7LIST),IALL_LU7)
!
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif
!     parallel stuff ends here
!     ------------------------

!
      IF(IDIAG.EQ.2) THEN
        LUDIA = LUSC1
      END IF
      IF(.NOT.(IDIAG.EQ.2.AND.IRESTR.EQ.1)) THEN
        ECOREP = 0.0D0
        IF(ICISTR.GE.2) CALL REWINO(LUDIA)
        SHIFT = ECORE_ORIG-ECORE
C       construct H diagonal
        CALL GASDIAT(WORK(KVEC1),LUDIA,SHIFT,ICISTR,I12,
     &               WORK(KLCBLTP),NBLOCK,WORK(KLCIBT),
     &               par_dist_block_list)
        IF(NOCSF.EQ.1.AND.ICISTR.EQ.1) THEN
          CALL REWINO(LUDIA)
          CALL TODSC_LUCI(WORK(KVEC1),number_of_dets,-1,LUDIA)
        END IF
        IF(IPRCIX.GE.2) WRITE(LUWRT,*) ' Diagonal constructed  '
      ELSE
         WRITE(LUWRT,*) ' Diagonal not calculated '
      END IF
#ifdef VAR_MPI
      IF (LUCI_NMPROC .GT. 1) THEN

C     debug printing
      NPTEST_VAR = 00
      IF( NPTEST_VAR .ge. 10 ) THEN
        ITEST_OFF = 0
        ITEST_OFF = MY_DIA_OFF
        DO IIBLK = 1, NUM_BLOCKS
            ILEN = 0
            CALL GET_BLOCK_PROC(par_dist_block_list,IIBLK,IPROC)
            IF( LUCI_MYPROC .eq. IPROC ) THEN
              CALL GET_BLOCK_LENGTH(block_list,IIBLK,ILEN)
              CALL MPI_FILE_READ_AT(IDIA,ITEST_OFF,WORK(KVEC1),ILEN,
     &                              MPI_REAL8,my_STATUS,IERR)
              WRITE(LUWRT,*) ' Diagonal elements',ILEN
              CALL WRTMATMN(WORK(KVEC1),1,ILEN,1,ILEN,LUWRT)
              ITEST_OFF = ITEST_OFF + ILEN
          END IF
        END DO
      ENDIF
      NPTEST_VAR = 00
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif
*
*. Transfer control to optimization routine

C     check for restart without CI --> we want to run an analysis only
      IF( MAXIT .lt. 0) GOTO 404
#ifdef VAR_MPI
      IF (LUCI_NMPROC .GT. 1) THEN
C
C     copy restart vector to MPI-file format
C
      IF(IRESTR .eq. 1)THEN
C
         WRITE(LUWRT,'(1X,A)')
     &   ' restart file LUCITA_CVECS.x will be transformed'//
     &   ' to MPIs file-i/o format ...'
C
         IF( LUCI_MYPROC .eq. LUCI_MASTER ) CALL REWINE(LUC,-1)
         CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU1,MY_LU1_OFF,
     &                         WORK(KILU1LIST),
     &                         par_dist_block_list,block_list,
     &                         NBLOCK,NROOT,irc_save)
        WRITE(LUWRT,'(1X,A)') ' ... done! '
      END IF
C
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif

      CALL CIEIG5(IRESTR,EROOT,WORK(KVEC1),WORK(KVEC2),LUDIA,
     &            LUC,LUHC,LUSC1,LUSC2,LUSC3,LUSC34,LUSC35,LUSC41,
     &            number_of_dets,NBLK,NROOT,MXCIV,MAXIT,LUCIVI,
     &            IPRNT,WORK(KSBEVC),
     &            0,WORK(KH0),WORK(KSBIDT),MXP1,MXP2,MXQ,
     &            WORK(KH0SCR),ECORE,ICISTR,LBLK,IDIAG,
     &            THRES_E,NBATCH,WORK(KLCLBT),WORK(KLCLEBT),
     &            WORK(KLCI1BT),WORK(KLCIBT),INIDEG,E_THRE,C_THRE,
     &            E_CONV,C_CONV,ICLSSEL,WORK(KLBLKCLS),NOCCLS,
     &            WORK(KLCLSC),WORK(KLCLSE), WORK(KLCLSCT),
     &            WORK(KLCLSET),WORK(KLCLSA),WORK(KLCLSL),WORK(KLBLKA),
     &            WORK(KLCLSD),WORK(KLCLSDT),ISKIPEI,WORK(KLC2B),
     &            WORK(KLCLSA2),LBLOCK,IROOTHOMING,WORK(KLBASSPC),
     &            WORK(KLEBASC),WORK(KLCBASC),NCISPC,
     &            ISPC,block_list,par_dist_block_list,WORK(KRCCTOS),
     &            IPROCLIST,IGROUPLIST
#ifdef VAR_MPI
     &            ,WORK(KILU1LIST),WORK(KILU2LIST),WORK(KILU3LIST),
     &            WORK(KILU4LIST),WORK(KILU5LIST),WORK(KILU6LIST),
     &            WORK(KILU7LIST),WORK(KILUCLIST)
#endif
     &            )
C
      EREF = EROOT(IRFROOT)
CTF
*  If TERACI has been chosen, only 1 root is determined.
*  So change NROOTS accordingly, otherwise the program crashes.
 404  if (IDIAG.eq.2) then
        NROOT_SAVE = NROOT
        NROOT = 1
      end if
C
C     rewind C-vector file
C
      CALL REWINE(LUC,-1)
C
#ifdef VAR_MPI

      IF (LUCI_NMPROC .GT. 1) THEN
!     close parallel files that are no longer needed
      call close_file_io_model(nr_files-1,0,fh_array)
C
      IF(IAM_NOT_INV .eq. 0 .or. MAXIT .lt. 0) GOTO 505
C
C     copy c-vectors from nodes and master back to the master
C
C     NOTE: only C-vectors will be saved.
C
      CALL REWINE(LUC,-1)
      DO JROOT = 1, NROOT
        CALL REWINE(LUSC41,-1)
        CALL COPVCD_PAR_BDRIV_REL(ILU1,LUSC41,WORK(KVEC1),
     &                            par_dist_block_list,NUM_BLOCKS2,
     &                            block_list,
     &                            MPI_COMM_WORLD,-1,JROOT,
     &                            WORK(KILU1LIST),MY_LU1_OFF,1)
        IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
          CALL REWINE(LUSC41,-1)
          CALL copvcd(LUSC41,LUC,WORK(KVEC1),0,-1)
        END IF
      END DO
      CALL REWINE(LUC,-1)
      NPTESTVAR = 00
      IF( NPTESTVAR .ge. 10 )THEN
        IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
          DO IVEC = 1, NROOT
            WRITE(LUWRT,*) '  final new C vector',IVEC
            CALL WRTVCD(WORK(KVEC1),LUC,0,-1)
          END DO
        END IF
      END IF
      CALL REWINE(LUC,-1)
      NPTESTVAR = 00
C
 505    CONTINUE
      close(unit=LUSC41,status='DELETE')
*
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif
C
      IF(LUCI_MYPROC.EQ.LUCI_MASTER)THEN
C
C       default: always print NO occ. numbers (IDENSI = 1)
C                which are useful for detailed analysis
        IF( IDENSI .gt. 0 )THEN
          IF( IPRDEN .le. 0 ) IPRDEN = 1 
        END IF
        IPRDEN = 100
        WRITE(luwrt,*) ' IPRDEN raised explicitly in GASCI '
C
      DO JROOT = 1, NROOT
        if (IPRDEN.ge.1) then
          WRITE(LUWRT,'(/2X ,A)')
     &    ' **************************************************'
          WRITE(LUWRT,'(2X ,A,I3)')
     &    ' Natural orbital occupation numbers for ROOT = ',JROOT
          WRITE(LUWRT,'(2X ,A)')
     &    ' **************************************************'
        end if
        IF(ICISTR.EQ.1) THEN
          CALL FRMDSC_LUCI(WORK(KVEC1),number_of_dets,-1,LUC,
     &                     IMZERO,IAMPACK)
          CALL DCOPY(number_of_dets,WORK(KVEC1),1,WORK(KVEC2),1)
          IF(IDENSI.EQ.1) THEN
            KRHO2 = 1
          END IF
          CALL DENSI2(IDENSI,WORK(KRHO1),WORK(KRHO2),
     &           WORK(KVEC1),WORK(KVEC2),0,0,EXPS2)
          if (IPRNCIV.ge.1)
     &    CALL GASANA(WORK(KVEC1),NBLOCK,WORK(KLCIBT),WORK(KLCBLTP),
     &                LUC,ICISTR)
        ELSE
          IF(IDENSI.GE.1.OR.NROOT.GT.1) THEN
            CALL REWINO(LUSC1)
            CALL COPVCD(LUC,LUSC1,WORK(KVEC1),0,-1)
            if(IDENSI .ge. 1)then
              CALL COPVCD(LUSC1,LUHC,WORK(KVEC1),1,-1)
            end if
          END IF
          if(idensi .ge. 1)then
            WRITE(luwrt,*) ' IDENSI in GASCI = ', IDENSI
            WRITE(luwrt,*) ' CALL to STANDARD DENSI2 '
            I12 = IDENSI
            CALL DENSI2(IDENSI,WORK(KRHO1),WORK(KRHO2),
     &                  WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC,EXPS2)
            ISIGDEN = 2
            WRITE(luwrt,*) ' CALL to SIGDEN_CI for density construction'
            CALL SIGDEN_CI(WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC,
     &                     CV_DUMMY,HCV_DUMMY,ISIGDEN
#ifdef VAR_MPI
     &               ,LUCLIST,LUHCLIST,IBLOCKL,NPARBLOCK,RCCTOS,
     &                IGROUPLIST,IPROCLIST
#endif
     &            )

          end if
C
          if (IPRNCIV.ge.1) then
            WRITE(LUWRT,'(/2X ,A)')
     &      ' **************************************************'
            WRITE(LUWRT,'(2X ,A,I3)')
     &      ' Analysis of leading coefficients for ROOT = ',JROOT
            WRITE(LUWRT,'(2X ,A)')
     &      ' **************************************************'
            IF (NROOT.GT.1) THEN
              CALL REWINO(LUSC1)
              CALL GASANA(WORK(KVEC1),NBLOCK,WORK(KLCIBT),WORK(KLCBLTP),
     &                    LUSC1,ICISTR)
            ELSE
              CALL REWINO(LUC)
              CALL GASANA(WORK(KVEC1),NBLOCK,WORK(KLCIBT),WORK(KLCBLTP),
     &                    LUC,ICISTR)
            end if
          end if
C
          IF (IPRNCIV.EQ.2) THEN
             WRITE(LUWRT,*)
             WRITE(LUWRT,*) ' ======================'
             WRITE(LUWRT,*) ' Complete CI expansion '
             WRITE(LUWRT,*) ' ======================'
             WRITE(LUWRT,*)
             IF(NROOT.EQ.1) THEN
               CALL WRTVCD(WORK(KVEC1),LUC  ,1,-1)
             ELSE
               CALL WRTVCD(WORK(KVEC1),LUSC1,1,-1)
             END IF
          END IF
*         ^ End if print of CI vector
        END IF
        IF (NPROP.GT.0) CALL ONE_EL_PROP
      END DO
      END IF
*    /\ end over master
*
CTF
C     If TERACI has been chosen, only 1 root is determined.
C     Reset number of roots.
C
      if (IDIAG.eq.2) NROOT = NROOT_SAVE
C
C     flush local memory
      deallocate(par_dist_block_list)
      deallocate(block_list)

      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUM,'GASCI')
C
C     eliminate scratch units
C     -----------------------
      close(unit=LUSC2, status='DELETE')
      close(unit=LUSC3, status='DELETE')
      close(unit=LUSC34,status='DELETE')
      close(unit=LUSC35,status='DELETE')
      close(unit=LUSC36,status='DELETE')
      close(unit=LUSC37,status='DELETE')
      close(unit=LUSC38,status='DELETE')
      close(unit=LUSC39,status='DELETE')
      close(unit=LUSC40,status='DELETE')
      close(unit=LU91,  status='DELETE')
C
      END
***********************************************************************
*
      SUBROUTINE LUCITA_DRIVER_old(WRK_DALTON,LWRK_DALTON)
*
* L U C I A
*
*
* CI for program for :FCI
*                     RASCI
*                     MRSDCI
*                     GAS GAS GAS GAS GAS GAS
*
* Written by Jeppe Olsen , winter of 1991
*                          GAS version in action summer of 95
*
* Modifications wrt to an MPI implementation by Stefan Knecht, Feb. 06 
*
!     dependencies on F90 modules
!     common modules
      use lucita_cfg
      use lucita_orbital_spaces
      use lucita_mcscf_ci_task
      use lucita_setup
      use lucita_energy_types
!     parallel lucita
#ifdef VAR_MPI
      use communication_model
      use file_io_model
      use sync_coworkers
#endif

      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "files.inc"
*. Parameters for dimensioning
#include "mxpdim.inc"
*.File numbers
#include "clunit.inc"
#include "units.inc"
*.Print flags
#include "cprnt.inc"
#include "lucinp.inc"
#include "cstate.inc"
! contains enviro
#include "crun.inc"
#include "cicisp.inc"
#include "oper.inc"
#include "cgas.inc"
*.Memory
#include "wrkspc.inc"
      real(8) :: XNFLOP    = 0.0d0
      real(8) :: XNCALL    = 0.0d0
      real(8) :: XLCROW    = 0.0d0
      real(8) :: XLCCOL    = 0.0d0
      real(8) :: XLCROWCOL = 0.0d0
      real(8) :: TMULT     = 0.0d0
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT

      DIMENSION WRK_DALTON(LWRK_DALTON)

      integer               :: files_to_close
      integer               :: fh_offset
      integer, allocatable  :: proclist(:)
      integer, allocatable  :: grouplist(:)
      integer, allocatable  :: grouplist_shared_mem(:)
      integer, allocatable  :: fh_array(:)
!     not yet a input parameter
      logical               :: set_common_blocks

      CALL QENTER('LUCITA_MASTER')

      write(lupri,'(/a)') 
     & '   ------------------------------------------------'
      write(lupri,'(a,a12)') '   LUCITA module started with task-id = ',
     &                           lucita_citask_id
      write(lupri,'(a/)') 
     & '   ------------------------------------------------'

!     set LUCITA internal control and task-individual common block information, 
!     insert defaults, cross check input and setup file handles
      call setup_lucita_cb_interface(lucita_citask_id)

!     control variable for invoking the setup procedure inside setup_lucita_orbital_string_cb
!     todo: make it flexible via ci/mcscf input
      set_common_blocks = .true.

!     set LUCITA internal orbital and string common block information 
      call setup_lucita_orbital_string_cb(set_common_blocks,
     &                                    iprorb)

      allocate(proclist(luci_nmproc))
      allocate(grouplist(luci_nmproc))

      proclist  = -1
      grouplist = -1

!     currently we work with 9 scratch files (nr_files in parluci.h)
      nr_files  = 9

#ifdef VAR_MPI
      IF (LUCI_NMPROC .GT. 1) THEN

        allocate(grouplist_shared_mem(luci_nmproc))
        grouplist_shared_mem = -1
      
!       setup the communication model.
!       all important variables (input parameters to module) 
!       are stored on common block /LUPARGROUP/

        call setup_communication_model(nflgrps,
     &                                 iiomod,
     &                                 shared_m,
     &                                 it_shl,
     &                                 ic_shl,
     &                                 grouplist,
     &                                 proclist,
     &                                 grouplist_shared_mem,
     &                                 luci_myproc,
     &                                 luci_nmproc,
     &                                 mpi_comm_world,
     &                                 mynew_id,
     &                                 icomm_id,
     &                                 mynew_id_sm,
     &                                 mynew_id_sm_c,
     &                                 newcomm_proc,
     &                                 icomm_size,
     &                                 newcomm_proc_sm,
     &                                 newcomm_proc_sm_c,
     &                                 mynew_comm,
     &                                 icomm,
     &                                 mynew_comm_sm,
     &                                 mynew_comm_sm_c,
     &                                 my_groupn,
     &                                 n_master,
     &                                 n_master_sm,
     &                                 n_master_sm_c,
     &                                 lupri)

        allocate(fh_array(nr_files))
        fh_array = 0
        call  setup_file_io_model(mynew_comm,
     &                            nr_files,
     &                            fh_array,
     &                            my_groupn,
     &                            newcomm_proc,
     &                            'parci',
     &                            lupri)
  
!       transfer file handles to common block /LUCIAPFILE/ (in parluci.h)
!       ILU1 has to be the last in the row because it is supposed to be
!       closed after we have copied the final solution vector back to
!       sequential format. Before doing that we want to free disk space by
!       deleting all other MPI-I/O files.
        IDIA = fh_array(1)
        ILUC = fh_array(2)
        ILU7 = fh_array(3)
        ILU6 = fh_array(4)
        ILU5 = fh_array(5)
        ILU4 = fh_array(6)
        ILU3 = fh_array(7)
        ILU2 = fh_array(8)
        ILU1 = fh_array(9)
  
!       free currently unused shared-memory group list
        deallocate(grouplist_shared_mem)

      END IF ! LUCI_NMPROC .GT. 1
#endif

!     initialize LUCITA work space
      call setup_lucita_initialize_work_space(wrk_dalton,
     &                                        work,
     &                                        kfree,
     &                                        kfrsav,
     &                                        lwrk_dalton,
     &                                        lfree,
     &                                        mxpwrd,
     &                                        luci_nmproc,
     &                                        nr_files,
     &                                        luwrt)
 
!     set LUCITA internal work space pointers, 
!     generate string and integral information on work space, 
!     and read in integrals.
      call setup_lucita_pointer_string_integrals_work_space(iprorb)

C     ==============================================================
C     Loop over GAS spaces and perform calcalations : CI, PERT, ....
C     ==============================================================
C
      ISKIPEI_INI = ISKIPEI
      I_HAVE_DONE_CC = 0

      DO JCMBSPC = 1, NCMBSPC
        ISKIPEI = 0
*
        WRITE(lupri,'(/15X,A)')
     &  '********************************'
        WRITE(lupri,'(15X,A)')
     &  ' ******************************'
        WRITE(lupri,'(15X,A,I3)')
     &  '   Calculations in space ', JCMBSPC
        WRITE(lupri,'(15X,A)')
     &  ' ******************************'
        WRITE(lupri,'(15X,A/)')
     &  '********************************'
C
        I_EXPAND = 1
C
C       Loop over Calculations in given space
        DO JSEQ = 1,  NSEQCI(JCMBSPC)
C
C         possible expansion of CI VECTORS
C         ================================
          IF(JCMBSPC.NE.1.AND.I_EXPAND.EQ.1)THEN
C
C           Restart from previous spaces ( Assuming a progressing sequence :
C           spaces are just added, not subtracted )
C           ( Used for perturbation calculations and CI calculations and CC )
C
            LUIN  = LUC
            LUOUT = LUSC2
            IF(ICISTR.EQ.1) THEN
              LBLK = XISPSM(IREFSM,JCMBSPC)
            ELSE
              LBLK = -1
            END IF
C           WRITE(luwrt,*) ' LBLK = ', LBLK
            CALL EXPCIV(IREFSM,1,LUIN,JCMBSPC,LUOUT,
     &                  LBLK,LUSC2,
     &                  NROOT,1,IDC,IPRDIA)
C           Last space where vectors were stored
            ISKIPEI = ISKIPEI_INI
C           Expanded vector will be used as initial vector in the
C           zero space calculation. Tell next CI to restart from
C           CI vectors
            IRESTR = 1
          END IF ! Expansion section
C
            IF(JSEQ.EQ.1.AND.JCMBSPC.EQ.2.AND.IRST2.EQ.0) THEN
C             No restart from previous vectors - IRST2 has been set to zero
              IRESTR = 0
              WRITE(lupri,*) ' No restart from previous vectors'
            END IF

!           main driver routine for GAS-CI Davidson calculation
            CALL GASCI(IREFSM,JCMBSPC,IPRDIA,EREF,
     &                 proclist,grouplist,fh_array)

*. Modified one-electron operator in first it
            IF(XLAMBDA.NE.1.0D0 .AND.JCMBSPC.EQ.1) THEN
*. Obtain modified operator for lambda calculations
              WRITE(lupri,*) ' Operator will be modified '
              CALL GENH1(XLAMBDA)
              CALL SCLH2(XLAMBDA)
            END IF
*. Transform CI coeffficients
            IF(ITRACI.NE.0) THEN
              WRITE(lupri,*) ' Control will be transferred to TRACI_CTL'
              CALL TRACI_CTL
            END IF
        END DO ! loop over calculations in a given expansion
*
      END DO ! loop over CI Expansions
*
*. Print info on matrix multiplier
#ifdef LUCI_DEBUG
      IF(LUCI_MYPROC.EQ.LUCI_MASTER)THEN
        CALL PR_MATML_STAT
      END IF
#endif

      IF (LUCI_NMPROC .GT. 1) THEN
!       close all remaining open parallel files (for GASCI: 1 --> ILU1)
        files_to_close = 1
        fh_offset      = 8
#ifdef VAR_MPI
        call close_file_io_model(files_to_close,fh_offset,fh_array)
        deallocate(fh_array)

!       reset communication groups
        call close_communication_model(mynew_comm,icomm,mynew_comm_sm,
     &                                 mynew_comm_sm_c)
#endif
      END IF

!     free memory
      call memrel('lucita.done',wrk_dalton,kfrsav,kfrsav,kfree,lfree)
      deallocate(grouplist)
      deallocate(proclist)
     
      write(lupri,'(/a,a12,a)') ' LUCITA run, task-id = ',
     &      lucita_citask_id(1:12),' finished with no errors.'

      call qexit('LUCITA_MASTER')

      END
***********************************************************************
