!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
***********************************************************************

      subroutine davidson_ci_driver(block_list,par_dist_block_list,
     &                              proclist,grouplist,fh_array,rcctos,
     &                              nblock,iprnt,vec1,vec2,c2,
     &                              docisrdft)
*
* CI optimization in defined GAS spaces
*
*
* Jeppe Olsen, Winter of 1995
*
#ifdef VAR_MPI
      use dalton_mpi
      use file_io_model
      use par_mcci_io
#endif
      use file_type_module
      use lucita_energy_types
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "clunit.inc"
#include "cstate.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "cgas.inc"
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      integer(MPI_INTEGER_KIND):: my_STATUS(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) ITEST_OFF
#endif
#include "parluci.h"
      real(8), intent(inout) :: vec1(*)
      real(8), intent(inout) :: vec2(*)
      real(8), intent(inout) :: c2(*)
      integer, intent(in)    :: proclist(luci_nmproc)
      integer, intent(in)    :: grouplist(luci_nmproc)
      integer, intent(inout) :: fh_array(file_info%active_nr_f_lucipar)
!     general information needed in parallel runs
      integer, intent(in)    :: par_dist_block_list(nblock)
      integer, intent(in)    :: block_list(nblock)
      integer, intent(in)    :: rcctos(nblock)
      logical, intent(in)    :: docisrdft
!-----------------------------------------------------------------------
      integer(8)             :: cvec_file_size
      integer                :: keep_nr_files
!-----------------------------------------------------------------------

      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)

      if(maxit <= 0)then
        write(luwrt,'(/a/)') '  @@ no CI iterations asked for'//
     &                       ' - return! :)'
        return
      end if

      write(luwrt,'(a,i14)') 
     &'  @@ Number of determinants/combinations: ',l_combi

      lblk                 = -1
      if(icistr == 1) lblk = l_combi

!     restart cross check (experimental; should work with ifort/gfortran >= 4.5)
      if(irestr .eq. 0 .and.  LUCI_MYPROC .eq. LUCI_MASTER )then
        cvec_file_size = 0
!       inquire(unit=luc, size=cvec_file_size)
        if(cvec_file_size .gt. 0) irestr = 1
      end if

#ifdef VAR_MPI
      IF(LUCI_NMPROC .GT. 1)THEN

!       update - might have changed because of the above (experimental) cross check
        call dalton_mpi_bcast(irestr,LUCI_MASTER,MPI_COMM_WORLD)

        call izero(file_info%iluxlist(1,2),file_info%max_list_length)
!       copy restart vector to MPI-file format
        IF(IRESTR .eq. 1)THEN
          WRITE(LUWRT,'(A)')
     &    '  restart file LUCITA_CVECS.x will be transformed'//
     &    ' to MPIs file-i/o format ...'
          IF( LUCI_MYPROC .eq. LUCI_MASTER ) CALL REWINE(LUC,-1)
          call mcci_cp_vcd_mpi_2_seq_io_interface(VEC1,LUC,ILU1,
     &                                            MY_LU1_OFF,
     &                                            file_info%
     &                                            iluxlist(1,2),
     &                                            par_dist_block_list,
     &                                            block_list,
     &                                            MPI_COMM_WORLD,
     &                                            NUM_BLOCKS2,NROOT,1,1)
          CALL REWINE(LUC,LBLK)
          WRITE(LUWRT,'(A)') '  ... done! '
        END IF

      end if ! (LUCI_NMPROC .GT. 1) THEN
#endif

!     transfer control to optimization routine
!     ----------------------------------------
      CALL CIEIG5(IRESTR,VEC1,VEC2,C2,LUDIA,
     &            LUC,LUHC,LUSC1,LUSC2,LUSC3,LUSC34,LUSC35,LUSC41,
     &            l_combi,NROOT,MXCIV,MAXIT,
     &            IPRNT,WORK(KSBEVC),0,WORK(KH0),WORK(KSBIDT),
     &            MXP1,MXP2,MXQ,WORK(KH0SCR),ECORE,ICISTR,LBLK,IDIAG,
     &            THRES_E,THRES_E_aux,INIDEG,LBLOCK,IROOTHOMING,NCISPC,
     &            ICSPC,block_list,par_dist_block_list,RCCTOS,
     &            PROCLIST,GROUPLIST,docisrdft,ispnden,idensi
#ifdef VAR_MPI
     &           ,file_info%iluxlist(1,2), file_info%iluxlist(1,3),
     &            file_info%iluxlist(1,5), file_info%iluxlist(1,6),
     &            file_info%iluxlist(1,7), file_info%iluxlist(1,8),
     &            file_info%iluxlist(1,9), file_info%ilublist
#endif
     &            )
 
#ifdef VAR_MPI
!     check for restart without CI
      IF(MAXIT .lt. 0) goto 777

      IF(LUCI_NMPROC .GT. 1)THEN
!       close (and re-open) parallel files that are no longer needed
!       which will free disk space for copyig back the vectors...
        keep_nr_files = 2
        call close_file_io_model(file_info%
     &                           active_nr_f_lucipar-keep_nr_files,
     &                           keep_nr_files,
     &                           fh_array)
!       re-open
        call setup_file_io_model(mynew_comm,file_info%
     &                           active_nr_f_lucipar-keep_nr_files,
     &                           fh_array,keep_nr_files,
     &                           my_groupn,newcomm_proc,
     &                           'parci',luwrt)
!       transfer file handles to common block /LUCIAPFILE/ (in parluci.h)
        ILU2 = fh_array( 3)
        IDIA = fh_array( 4)
        ILUC = fh_array( 5)
        ILU3 = fh_array( 6)
        ILU4 = fh_array( 7)
        ILU5 = fh_array( 8)
        ILU6 = fh_array( 9)
        ILU7 = fh_array(10)

!       copy c-vectors from nodes and master back to the master
        CALL REWINE(LUC,-1)

        call mcci_cp_vcd_mpi_2_seq_io_interface(VEC1,LUC,ILU1,
     &                                          MY_LU1_OFF,
     &                                          file_info%iluxlist(1,2),
     &                                          par_dist_block_list,
     &                                          block_list,
     &                                          MPI_COMM_WORLD,
     &                                          NUM_BLOCKS2,NROOT,1,2)
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      CALL REWINE(LUC,LBLK)
      IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
        DO IVEC = 1, NROOT
          WRITE(LUWRT,*) '  final solution vector for root ==> ',IVEC
          CALL WRTVCD(VEC1,LUC,0,LBLK)
        END DO
      END IF
      CALL REWINE(LUC,LBLK)
#undef LUCI_DEBUG
#endif
      !> save final vectors on different scratch file
      if(docisrdft)then
        CALL REWINE(LUC,LBLK)
        luci_cvec = 99
        open(file='srdft-lucita-final.cvec',unit=luci_cvec,
     &       status='replace',
     &       form='unformatted',action='write',position='rewind')
        do ivec = 1, nroot
          call copvcd(luc,luci_cvec,vec1,0,-1)
        end do
        close(luci_cvec, status='keep')
        CALL REWINE(LUC,LBLK)
      end if


!     eliminate scratch units
!     -----------------------
 777  close(unit=LUSC2, status='DELETE')
      close(unit=LUSC3, status='DELETE')
      close(unit=LUSC34,status='DELETE')
      close(unit=LUSC35,status='DELETE')
      close(unit=LUSC36,status='DELETE')
      close(unit=LUSC37,status='DELETE')
      close(unit=LUSC38,status='DELETE')
      close(unit=LUSC40,status='DELETE')
      close(unit=LUSC41,status='DELETE')
      close(unit=LU91,  status='DELETE')
      if(luci_myproc.eq.luci_master)
     &close(unit=LUSC1, status='DELETE')

!     re-open partially used scratch units (sequential run)
!     -----------------------------------------------------
      if(luci_myproc.eq.luci_master)then
        open(unit=LUSC1,file="LUSC1.0",status='UNKNOWN',
     &       form='UNFORMATTED')
        open(unit=LUSC2,file="LUSC2.0",status='UNKNOWN',
     &       form='UNFORMATTED')
        open(unit=LUSC3,file="LUSC3.0",status='UNKNOWN',
     &       form='UNFORMATTED')
      end if
 
      END
***********************************************************************

      SUBROUTINE CIEIG5(IRESTR,VEC1,VEC2,C2,LUDIA,LU1,LU2,LU3,
     &                  LU4,LU5,LU6,LU7,LU41,NDIM,NROOTSL,MAXVEC,
     &                  MXCIIT,IPRT,PEIGVC,NPRDET,H0,IPNTR,
     &                  NP1,NP2,NQ,H0SCR,EIGSHF,ICISTR,LBLK,IDIAG,
     &                  THRES_E,THRES_E_aux,INIDEG,
     &                  MXLNG,IROOTHOMING,NSPC,ISPC,IBLOCKAR,
     &                  IBLKDSTND,RCCTOS,IPROCLIST,IGROUPLIST,docisrdft,
     &                  ispnden,idensi
#ifdef VAR_MPI
     &                  ,LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                  LU6LIST,LU7LIST,LUCLIST
#endif
     &                  )

!
! Master routine for CI diagonalization
!
#ifdef MOD_SRDFT
      use lucita_mcscf_srdftci_cfg
#endif

      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLXLBT, KLXLEBT, KLXI1BT, KLXIBT
      INTEGER*8 current_free_mem, ksrdftdal_interface
!               for addressing of WORK
      INTEGER   lscr_local
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      integer(MPI_INTEGER_KIND):: my_STATUS(MPI_STATUS_SIZE)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
      DIMENSION LU7LIST(*), LU6LIST(*)
#endif
#include "priunit.h"
#include "parluci.h"
#include "files.inc"
#include "mxpdim.inc"
! Definition of c and sigma
#include "cands.inc"
! NACOB used
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "gasstr.inc"
#include "cgas.inc"
! NSMOB used
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
      LOGICAL CA
      LOGICAL docisrdft
      integer, intent(in) :: ispnden, idensi
C
C
      REAL*8   VEC1(*),VEC2(*),C2(*)
      INTEGER, PARAMETER    :: LLWRK = 500000
      integer, allocatable  :: iscr1(:)
      real(8), allocatable  :: scr1(:)
      real(8), allocatable  :: fsr(:)
      integer               :: nnorbt_internal, nnashx_internal
*. Output from Subspace dagonalization
      DIMENSION H0(*),IPNTR(*),H0SCR(*),PEIGVC(*)
*.
      INTEGER   RCCTOS(*)
      real*8 :: s2value

!     use dynamic instead of static allocation (+ common block)
      allocate(iscr1(llwrk))
      allocate(scr1(llwrk))
*
      IF( IPRT.GT. 1 )  WRITE(luwrt,'(/A)')
     &'          *** information from ci diagonalization  ***'
!     NTEST = 9999 ! NTEST=0
      NTEST = 0
C
C               ====================================
C  1 :               INITIAL VARIATIONAL SUBSPACE
C               ====================================
C
      ninvec       = nrootsl
      ninvec_local = 0
      IF( IRESTR .EQ. 0 ) THEN
        IF(NPRDET .EQ. 0 ) THEN
C         ==================================================
C          Initial guess from lowest elements of CI diagonal
C         ==================================================
C
C         in order treat degeneracies, the lowest 6 * NROOTSL
C         elements are obtained
C
          NFINDM = MIN(NDIM,6*NROOTSL+6)

          ninvec = min(NFINDM,4*nrootsl+6)
!!!! DEBUG
!         NFINDM = 12; ninvec = 12 
!         print *, 'ninvec, nfindm, nrootsl, ndim, maxvec ...', 
!    &              ninvec, nfindm, nrootsl, ndim, maxvec
!!!! DEBUG
          if(ninvec > maxvec) ninvec = maxvec

          IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
              CALL FNDMND_PAR(LUDIA,LBLK,VEC1,NFINDM,NFINDA,
     &                    ISCR1(1+2*NFINDM),SCR1(1+2*NFINDM),ISCR1,
     &                    SCR1,IBLOCKAR,IBLKDSTND,NUM_BLOCKS2,IPRT)
#endif
          ELSE
             CALL FNDMND(LUDIA,LBLK,VEC1,NFINDM,NFINDA,
     &            ISCR1(1+2*NFINDM),SCR1(1+2*NFINDM),ISCR1,SCR1,IPRT)
          END IF
          CALL REWINE(LU1,-1)
          IBASE   = 1
          TEST    = 1.0D-10
!test     inideg  = -1
          irootsi =  1

         write(luwrt,1410)
     &   ninvec,(i,iscr1(I),scr1(i)-eigshf,scr1(i),i=1,ninvec)
 1410 FORMAT(/' (CIEIG5)',I3,' lowest diagonal elements:',
     &      //' Element no. Config.no.    Active energy',
     &        '      Total energy',
     &      //,(I10,' : ',I10,2F18.10))


          DO

            if(irootsi > ninvec) exit
*. Number of degenerate elements
            NDEG = 1
            XVAL = SCR1(IBASE)
   90       CONTINUE
            IF(IBASE-1+NDEG+1 <= NFINDA) THEN
              IF (ABS(SCR1(IBASE-1+NDEG+1)-XVAL) <= TEST) THEN
                NDEG = NDEG + 1
                GOTO 90
              END IF
            END IF

!#define TEST_THIS_VERSION_DEGENERACIES

#ifdef TEST_THIS_VERSION_DEGENERACIES

            IF (INIDEG.EQ.0.AND.NDEG.GT.1) THEN
!             WRITE(luwrt,*) ' NOTE! degenerate initial vectors for CI'
!             WRITE(luwrt,*) ' degree of degeneracy: ',NDEG 
              NDEG = 1
            END IF

!  Initial guess in compressed form in SCR1
            SCALE = 1.0D0/SQRT(DFLOAT(NDEG))
            DO 250 II = 1,NDEG
!  Anti symmetric combination
              IF(INIDEG.EQ.-1) THEN
                SCR1(II) = (-1.0D0)**II * SCALE
!  Symmetric combination
              ELSE IF (INIDEG.EQ.1.OR.INIDEG.EQ.0) THEN
!               SCR1(II) =  SCALE
                SCR1(II) = (-1.0D0)**(II-1) * SCALE
!               write(luwrt,*) 'scr1(ii), ii ==> ',scr1(ii), ii
              END IF
  250       CONTINUE
            IF(IDIAG.EQ.2) THEN
              JPACK = 1
            ELSE
              JPACK = 0
            END IF
            do ii = 1, ndeg

              if(ii > 1) SCR1(II) = (-1.0D0)**(II-1) * SCR1(II)

              IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
                 CALL WRSVCD_PAR(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,
     &                           IBLOCKAR,IBLKDSTND,NUM_BLOCKS2,JPACK,
     &                           IROOTSI,LU1LIST)
#endif
              ELSE
                 CALL WRSVCD(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
     &                       LUDIA,JPACK)
              END IF
            end do

#else
            IF (INIDEG.EQ.0.AND.NDEG.GT.1) THEN
!             WRITE(luwrt,'(a,i3,a,i4,a)') 
!    &        ' degenerate initial vector for CI:'//
!    &        ' vector ',irootsi,' out of ',ninvec,' start vectors'
!             WRITE(luwrt,*) ' degree of degeneracy: ',NDEG 
!             print *, 'plus_combi... ',plus_combi
              if(.not.plus_combi) NDEG = 1
            END IF

*. Initial guess in compressed form in SCR1
            SCALE = 1.0D0/SQRT(DFLOAT(NDEG))
            DO 250 II = 1,NDEG
*. Anti symmetric combination
              IF(INIDEG.EQ.-1) THEN
                SCR1(II) = (-1.0D0)**II * SCALE
*. Symmetric combination
              ELSE IF (INIDEG.EQ.1.OR.INIDEG.EQ.0) THEN
                SCR1(II) =  SCALE
              END IF
  250       CONTINUE
            IF(IDIAG.EQ.2) THEN
              JPACK = 1
            ELSE
              JPACK = 0
            END IF
            IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
               CALL WRSVCD_PAR(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,
     &                      IBLOCKAR,IBLKDSTND,NUM_BLOCKS2,JPACK,
     &                      IROOTSI,LU1LIST)
#endif
            ELSE
               CALL WRSVCD(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
     &                  LUDIA,JPACK)
            END IF
#endif
            do ii = 1,ndeg
              scr1(ii) = 0.0d0
            end do

            ninvec_local = ninvec_local + 1
            ibase        = ibase   + ndeg
            irootsi      = irootsi + ndeg 
          end do
        ELSE
* =====================================
*. Initial approximations are in PEIGVC
* =====================================
          CALL REWINE(LU1,-1)
          IF(IDIAG.EQ.2) THEN
            JPACK = 1
          ELSE
            JPACK = 0
          END IF
          ninvec_local = ninvec
          DO 1984 IROOTSI = 1, NROOTSL
            CALL WRSVCD(LU1,LBLK,VEC1,IPNTR,
     &           PEIGVC((IROOTSI-1)*NPRDET+1),NPRDET,NDIM,LUDIA,JPACK)
 1984     CONTINUE
        END IF
        ninvec = ninvec_local
        if(ninvec < nrootsl) call quit('less startvec than CI roots')
      END IF

!     release scratch memory
      deallocate(iscr1)

!
!     *****************************************************************
!     CI-DFT contributions
!     *****************************************************************
!
      IF (DOCISRDFT) THEN
#ifndef MOD_SRDFT
         call quit('srdft not implemented in this version')
#else

#ifdef VAR_MPI
         call quit('help me stefan... parallel srdft-ci')
#endif

         idum = 0
!        set local marker + allocate space for lucita-dalton integral interface
         call memman(idum,idum,'MARK  ',idum,'srdft1')
         nnorbt_internal = 0
         norbt_internal  = 0
         do i = 1, nirrep
           nnorbt_internal = nnorbt_internal + 
     &                       (nmos_env(i)*(nmos_env(i)+1))/2
           norbt_internal  = norbt_internal + nmos_env(i)
         end do 
         nnashx = (ntoob*(ntoob+1))/2

         allocate(fsr(max(ntoob**2,nnorbt_internal)))
         fsr              = 0.0d0

         call dzero(work(KRHO1_ens),ntoob**2)

         i12_save = i12
         i12      = 1 ! only 1-body terms...
         if(ispnden > 0)then
           i12      = ispnden
         end if

         luci_cvec = 99
         open(file='srdft-lucita.cvec',unit=luci_cvec,status='replace',
     &        form='unformatted',action='write',position='rewind')

!        debug print
!        call wrtvcd(vec1,lu1,-1,-1)

         call rewine(lu1,-1)

! Manu feb. 2014: weights(1) set to 1.0 when one single root is computed
!                 like in GS CI-srDFT calculations for example

         if(nrootsl ==1)then 
         weights(1) = 1.0d0
         end if
         !write(luwrt,*) 'weights(1:nroots)',weights(1:nrootsl)


!        calculate ensemble density matrix
!        ---------------------------------
         do i = 1, nrootsl

           call rewine(lu5,-1)
           call rewine(lu6,-1)
           call copvcd(lu1,lu5,vec1, 0,-1)
           call copvcd(lu5,lu6,vec1,-1,-1)

!          save initial CI vectors on scratch file
           call rewine(lu5,-1)
           call copvcd(lu5,luci_cvec,vec1,0,-1)

           call rewine(lu5,-1)
           call rewine(lu6,-1)

           xxs2dum = 0.0d0
!                                                            1: sigma vector
!                                                            2: density matrix
           call sigden_ci(vec1,vec2,c2,lu5,lu6,cdummy,sdummy,2, 0,
     &                    s2value)                          !   1: rhotype == 1 (regular density matrix)
!                                                           !   3: rhotype == 3 (transition density matrix)
!                                                           !   0: rhotype == 0 (regular density matrix run in CI)

        if(i12 > 1)then
          write(luwrt,'(/a      )')
     &'   ------------------------------------------------'
          write(luwrt,'(a,1f10.3)')
     &'   expectation value of operator <S**2> =', s2value
          write(luwrt,'(a/      )')
     &'   ------------------------------------------------'
        end if

!********************
!!! Odile+Manu debug
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
           CALL HEADER(
     &     'lucita: individual one-el. density matrix, active part',-1)
           CALL OUTPUT(work(krho1),1,ntoob,1,ntoob,ntoob,ntoob,
     &                 -11,luwrt)
           write(luwrt,*) 'weights(i)',weights(i)
#undef LUCI_DEBUG
#endif
!!!!
!********************

           call daxpy(ntoob**2,weights(i),work(krho1),1,
     &                work(KRHO1_ens),1)

         end do ! loop over # roots

!        debug print
!        CALL HEADER(
!    &   'lucita: ensemble one-el. density matrix, active part',-1)
!        CALL OUTPUT(work(KRHO1_ens),1,ntoob,1,ntoob,ntoob,ntoob,
!    &               -11,luwrt)

         !> i12 set to 1 since we only need the 1p-density matrix to save on file
         i12 = 1

         call lucita_putdens_generic(work(krho1_ens),work(krho2),
     &                              fsr,xxdummy,xxdummy,i12,2,1,nrootsl)
!                                                            2: density matrix
!                                                            1: rhotype
#ifdef LUCI_DEBUG
         CALL HEADER(
     &   'lucita: ensemble one-el. density matrix, ap reord   ',-1)
         CALL OUTPAK(work(krho1_ens),ntoob,-1,LUPRI)
#endif

!        reset i12 to allow for 1 and most likely 2-electron contributions
         i12 = i12_save

!        deallocate(fsr)
         close(luci_cvec,status='keep')

!        dump ensemble density on file
         call lucita_putdens_ensemble(work(KRHO1_ens))

!        calculate ensemble srDFT contributions
!        --------------------------------------

!        allocate(fsr(nnorbt_internal))
         fsr = 0.0d0

         IDUM = 0
         CALL MEMMAN(current_free_mem,0,'SFREEM',2,'SHOWMM')
         LSCR_local = current_free_mem - 1000
         call memman(ksrdftdal_interface,lscr_local,'ADDL  ',2,'srdft2')

         
         CALL SRFMAT(FSR,srdft_cmo_lucita,work(KRHO1_ens),
     &               EJCSR_mc2lu,EJVSR_mc2lu,EDSR_mc2lu,EDFT_mc2lu,
     &               emydftaux_mc2lu,UEJCVSR,
     &               work(ksrdftdal_interface),lscr_local,0)
         
         if(luci_myproc == luci_master)then
!          ... Correct EMY/eigshf and CI diagonal
           EMYDFT_mc2lu = 
     &            EJCSR_mc2lu + EJVSR_mc2lu + EDSR_mc2lu + EDFT_mc2lu(1)

!          write(luwrt,*) 'eigshf before',eigshf
!          in lucita: emy is eigshf at this point in the code
           eigshf = eigshf + EMYDFT_mc2lu
!          write(luwrt,*) 'eigshf after',eigshf
!          ... extract and save SR over active indices
           CALL GETAC(fsr,srdft_srac_lucita)

!          debug print
!          call wrtmatmn(srdft_srac_lucita,1,nnashx,1,nnashx,luwrt)
!          call wrtmatmn(work(kint1),1,3,1,3,luwrt)

           call dzero(work(krho1),ntoob**2)

!          add the contributions over the active indices (modify h1)
!          ------------------------------------------------------------
           call lucita_srdft_h1_adapt(srdft_srac_lucita,work(krho1),
     &                                work(kint1))

!          call wrtmatmn(work(kint1),1,3,1,3,luwrt)
!          ... fold SRAC for 
!          -----------------
!           1) Adding it to CI-diagonal to get better pre-conditioning
!           2) Later printing of CIDFT specific energy contributions
           CALL DSPTGE(ntoob,srdft_srac_lucita,fsr)
           CALL DGEFSP(ntoob,fsr,srdft_srac_lucita)

           ESRDV = DDOT(nnashx,work(KRHO1_ens),1,srdft_srac_lucita,1)
           write(luwrt,*) 'ESRDV_ref = ',esrdv
           write(luwrt,*) 'EJVSR_ref = ',EJVSR_mc2lu
         end if

#ifdef VAR_MPI
         if(luci_nmproc > 1)then 
           call mpi_bcast(eigshf,1,mpi_real8,0,mpi_comm_world,irr)
           call mpi_bcast(emydft_mc2lu,1,mpi_real8,0,mpi_comm_world,irr)
           call mpi_bcast(esrdv,1,mpi_real8,0,mpi_comm_world,ierr)
         end if
#endif

!        adding esrdv to H diagonal 
!        ----------------------------
         call manipulate_vcd(ludia,lu5,EMYDFT_mc2lu+ESRDV,
     &                       vec1,-1,-1)

C        eliminate local memory
         IDUM = 0
         CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'srdft1')
         deallocate(fsr)
#endif

!#define TEST_START_VEC_S2
#ifdef TEST_START_VEC_S2
      else ! TEST stefan

         print *, 'bla TEST ',i12,ispnden
         luci_cvec = 99
         open(file='srdft-lucita.cvec',unit=luci_cvec,status='replace',
     &        form='unformatted',action='write',position='rewind')

!        debug print
!        call wrtvcd(vec1,lu1,-1,-1)
         call rewine(lu1,-1)

         do i = 1, nrootsl

           call rewine(lu5,-1)
           call rewine(lu6,-1)
           call copvcd(lu1,lu5,vec1, 0,-1)
           call copvcd(lu5,lu6,vec1,-1,-1)

!          save initial CI vectors on scratch file
           call rewine(lu5,-1)
           call copvcd(lu5,luci_cvec,vec1,0,-1)

           call rewine(lu5,-1)
           call rewine(lu6,-1)

           xxs2dum = 0.0d0
!                                                            1: sigma vector
!                                                            2: density matrix
           call sigden_ci(vec1,vec2,c2,lu5,lu6,cdummy,sdummy,2, 0,
     &                    s2value)                          !   1: rhotype == 1 (regular density matrix)
!                                                           !   3: rhotype == 3 (transition density matrix)
!                                                           !   0: rhotype == 0 (regular density matrix run in CI)

        if(i12 > 1)then
          write(luwrt,'(/a      )')
     &'   ------------------------------------------------'
          write(luwrt,'(a,1f10.3)')
     &'   expectation value of operator <S**2> =', s2value
          write(luwrt,'(a/      )')
     &'   ------------------------------------------------'
        end if

!********************
!!! Odile+Manu debug
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
           CALL HEADER(
     &     'lucita: individual one-el. density matrix, active part',-1)
           CALL OUTPUT(work(krho1),1,ntoob,1,ntoob,ntoob,ntoob,
     &                 -11,luwrt)
#undef LUCI_DEBUG
#endif
!!!!
!********************

         end do ! loop over # roots
         close(luci_cvec,status='keep')

#endif /* TEST_START_VEC_S2 */

      ENDIF ! DOCISRDFT

!
!     initialize arrays
      do i = 1, nrootsl
        eroot(i)          = 0.0d0
        root_residual(i)  = 0.0d0
        root_converged(i) = i
      end do
      nfinal_vec = ninvec
!
!     check for initial CI vector generation only
      IF(MXCIIT .lt. 0) goto 999

!     =================
!      Diagonalization
!     =================
!

      CALL QENTER('CIEIG')
! Inverse iteration modified Davidson with 2 vectors in core
!     write(luwrt,*) ' THRES_E is ',THRES_E
!     THRES_E = 1.0d-12
!     write(luwrt,*) ' bla bla THRES_E reset...'
      IF(IPRT .GE. 5 ) THEN
         WRITE(luwrt,'(A,I3)')
     &   '  Number of roots to be converged..  ',NROOTSL
         WRITE(luwrt,'(A,I3)')
     &   '  Largest allowed number of vectors..',MAXVEC
         WRITE(luwrt,'(A,I3)')
     &   '  Allowed number of CI iterations  ..',MXCIIT
      END IF

*     allocations for SCR1
      KRNRM  = 1
      KEIG   = KRNRM + MXCIIT*ninvec
      KFIN   = KEIG  + MXCIIT*ninvec
      KAPROJ = KFIN + ninvec
      KAVEC  = KAPROJ + MAXVEC*(MAXVEC+1)/2
      KWORK  = KAVEC + MAXVEC ** 2
      KLFREE = KWORK + MAXVEC*(MAXVEC+1)
      IF( IPRT .GE. 5 ) THEN
         WRITE(luwrt,*) ' KRNRM KEIG KFIN KAPROJ KAVEC KWORK KLFREE '
         WRITE(luwrt,'(6I8)')KRNRM,KEIG,KFIN,KAPROJ,KAVEC,KWORK,KLFREE
      END IF
      IF( KLFREE-1 .GT. LLWRK) THEN
         WRITE(luwrt,'(A,2I14)' )
     &   ' Not enough memory in CIEIG5 : neeeded and available ',
     &     KLFREE-1, LLWRK
         WRITE(luwrt,'(A,I14)' )
     &   ' Increase parameter LLWRK in CIEIG5 to   ', KLFREE-1
         call quit( ' insufficient memory in cieig5 ' )
       END IF
*
      IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
C
C     partition CI vector with max. batch length using LBLOCK
C
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CIEIG5')
C
      IATP = 1
      IBTP = 2
C     Arrays for partitioning of rhs/lhs vector (symmetric because the
C     Hamiltonian is totally symmetric)
      NTTS = MXNTTS
      CALL MEMMAN(KLXLBT,   NTTS,'ADDL  ' ,1,'CLBT  ')
      CALL MEMMAN(KLXLEBT,  NTTS,'ADDL  ' ,1,'CLEBT ')
      CALL MEMMAN(KLXI1BT,  NTTS,'ADDL  ' ,1,'CI1BT ')
      CALL MEMMAN(KLXIBT, 8*NTTS,'ADDL  ' ,1,'CIBT  ')

      call Z_BLKFO_partitioning_parallel(icspc,icsm,iatp,ibtp,
     &                                   WORK(KLXLBT),WORK(KLXLEBT),
     &                                   WORK(KLXI1BT),WORK(KLXIBT),
     &                                   nbatch_par,nblock_par,
     &                                   IBLKDSTND)
 
       CALL MICDV6_PAR(VEC1,VEC2,C2,SCR1(KRNRM),SCR1(KEIG),nfinal_vec,
     &                 EROOT,root_converged,root_residual,
     &                 MXCIIT,NDIM,NROOTSL,
     &                 MAXVEC,ninvec,SCR1(KAPROJ),SCR1(KAVEC),
     &                 SCR1(KWORK),IPRT,NPRDET,H0,IPNTR,NP1,NP2,NQ,
     &                 H0SCR,LBLK,EIGSHF,THRES_E,IROOTHOMING,
     &                 luwrt,c_state_of_interest,
     &                 DOCISRDFT.and.NROOTSL>1,weights,
     &                 IBLOCKAR,IBLKDSTND,
     &                 RCCTOS,
     &                 LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                 LU6LIST,LU7LIST,LUCLIST,NBATCH_par,
     &                 WORK(KLXLBT),WORK(KLXLEBT),WORK(KLXI1BT),
     &                 WORK(KLXIBT),IPROCLIST,IGROUPLIST)

#endif
      ELSE ! (LUCI_NMPROC .GT. 1) THEN
!#define blubb
#ifndef blubb
       CALL MICDV6(VEC1,VEC2,C2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &             nfinal_vec,EROOT,root_converged,root_residual,
     &             MXCIIT,NDIM,LU3,LU4,LU5,LU6,LU7,
     &             LUDIA,NROOTSL,
     &             MAXVEC,ninvec,SCR1(KAPROJ),SCR1(KAVEC),
     &             SCR1(KWORK) ,IPRT,NPRDET,H0,IPNTR,NP1,NP2,NQ,
     &             H0SCR,-1,EIGSHF,THRES_E,IROOTHOMING,luwrt,
     &             c_state_of_interest,
     &             DOCISRDFT.and.NROOTSL>1,weights)
!
!!! Manu july 2014: for ensemble CI-srDFT, wavefunctions will be
!converged for a given trial ensemble density. The latter should then be
!updated in the srDFT potential (call SRFMAT ...), and MICDV6 called again until
!convergence is reached.   
!
#else
#ifndef blubb2
       CALL blubb_MICDV4(VEC1,VEC2,C2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &             nfinal_vec,EROOT,root_converged,root_residual,
     &             MXCIIT,NDIM,LU3,LU4,LU5,LU6,LU7,
     &             LUDIA,NROOTSL,
     &             MAXVEC,ninvec,SCR1(KAPROJ),SCR1(KAVEC),
     &             SCR1(KWORK) ,IPRT,NPRDET,H0,IPNTR,NP1,NP2,NQ,
     &             H0SCR,-1,EIGSHF,THRES_E,luwrt)
#else
       CALL MICDV4_enlmd(VEC1,VEC2,C2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &             nfinal_vec,EROOT,root_converged,root_residual,
     &             MXCIIT,NDIM,LU3,LU4,LU5,LU6,LU7,
     &             LUDIA,NROOTSL,
     &             MAXVEC,ninvec,SCR1(KAPROJ),SCR1(KAVEC),
     &             SCR1(KWORK) ,IPRT,NPRDET,H0,IPNTR,NP1,NP2,NQ,
     &             H0SCR,-1,EIGSHF,THRES_E,luwrt)
#endif
#endif
      END IF ! (LUCI_NMPROC .GT. 1) THEN

#ifdef VAR_MPI
      IF (LUCI_NMPROC .GT. 1) THEN
C       eliminate local memory
        IDUM = 0
        CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'CIEIG5')
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif

      CALL QEXIT('CIEIG')

!     release scratch memory
 999  deallocate(scr1)

      END
***********************************************************************
