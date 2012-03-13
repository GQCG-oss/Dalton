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
!
!...  This file gets the AO density matrix from SIRIFC.
!
!...  2012-03-09, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief gets the atomic orbital (AO) density matrix
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param size_ao is the size of AO density matrix
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \return ao_dens contains the AO density matrix
  subroutine get_ao_dens(size_ao, ao_dens, len_work, wrk_space)
    implicit none
    integer, intent(in) :: size_ao
    real(REALK), intent(out) :: ao_dens(size_ao)
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    ! uses NCMOT, NNASHX, N2BASX
#include <inforb.h>
    ! uses LUSIFC
#include <inftap.h>
    integer start_dv_mo     !start of active part of one-electron density matrix (MO) in the workspace
    integer start_dv_ao     !start of active part of one-electron density matrix (AO) in the workspace
    integer start_left_wrk  !start of left workspace
    integer len_left_wrk    !length of left workspace
    logical GETDC           !if calculating DCAO
    logical GETDV           !if calculating DVAO
    integer idummy          !dummy stuff
    logical found           !if found required data from SIRIFC
    ! pushes current subroutine into the stack
    call QENTER("get_ao_dens")
    ! start of active part of one-electron density matrix (MO)
    start_dv_mo = 1+NCMOT
    ! start of active part of one-electron density matrix (AO)
    start_dv_ao = start_dv_mo+1
    ! start of left workspace
    start_left_wrk = start_dv_ao+1
    ! only calculates DCAO
    GETDC = .true.
    GETDV = .false.
    ! lenght of the left workspace
    len_left_wrk = len_work-start_left_wrk+1
    ! checks if the workspace is enough
    if (len_left_wrk<0) call STOPIT("GEN1INT", "get_ao_dens", start_left_wrk-1, len_work)
    ! opens SIRIFC
    if (LUSIFC<=0) &
      call GPOPEN(LUSIFC, "SIRIFC", "OLD", " ", "UNFORMATTED", idummy, .false.)
    rewind(LUSIFC)
    ! reads the molecular orbital coefficients
#if REALK == 4
    call SZERO(wrk_space(1), NCMOT)
#elif REALK == 8
    call DZERO(wrk_space(1), NCMOT)
#else
    call QUIT("Unknown kind of real numbers!")
#endif
    call RD_SIRIFC("CMO", found, wrk_space(1), wrk_space(start_left_wrk), &
                   len_left_wrk)
    if (.not.found) call QUIT("CMO IS NOT FOUND ON SIRIFC!")
    ! reads active part of one-electron density matrix (MO)
    if (GETDV) then
#if REALK == 4
      call SZERO(wrk_space(start_dv_mo), NNASHX)
#elif REALK == 8
      call DZERO(wrk_space(start_dv_mo), NNASHX)
#else
      call QUIT("Unknown kind of real numbers!")
#endif
      call RD_SIRIFC("DV", found, wrk_space(start_dv_mo), &
                     wrk_space(start_left_wrk), len_left_wrk)
      if (.not.found) call QUIT("DV IS NOT FOUND ON SIRIFC!")
#if REALK == 4
      call SZERO(wrk_space(start_dv_ao), N2BASX)
#elif REALK == 8
      call DZERO(wrk_space(start_dv_ao), N2BASX)
#else
      call QUIT("Unknown kind of real numbers!")
#endif
    end if
    ! gets the AO density matrix, using
    !
    ! FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LFRSAV)
    ! Input:
    !   GETDC   if true calculate DCAO
    !   GETDV   if true calculate DVAO
    !   CMO(*)  molecular orbital coefficients
    !   DV(*)   active part of one-electron density matrix (over MO's)
    ! Scratch:
    !   WRK(LFRSAV)
    call FCKDEN(GETDC, GETDV, ao_dens, wrk_space(start_dv_ao), wrk_space(1), &
                wrk_space(start_dv_mo), wrk_space(start_left_wrk), len_left_wrk)
    ! sums DCAO and DVAO
    if (GETDV) then
      start_dv_ao = start_dv_ao-1
      do idummy = 1, N2BASX
        ao_dens(idummy) = ao_dens(idummy)+wrk_space(start_dv_ao+idummy)
      end do
    end if
    ! pops the stack
    call QEXIT("get_ao_dens")
    return
  end subroutine get_ao_dens
