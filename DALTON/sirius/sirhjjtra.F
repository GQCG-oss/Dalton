!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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

!  FILE: DALTON/sirius/sirhjjtra.F
!  (c) Copyright Hans Joergen Aa. Jensen, hjj@sdu.dk (2014)

      subroutine SIR_INTOPEN_HJJTRA()
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
#include "priunit.h"

      write(LUPRI,*) 'SIR_INTOPEN_HJJTRA called'

      end

      subroutine HJJTRA_CTL(KTRLVL, CMO, WORK, KWORK)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer KTRLVL, KWORK
      real*8  CMO(*), WORK(KWORK)
#include "priunit.h"

! Used from include files:
!  inforb.h : ???
#include "inforb.h"

      write(LUPRI,*) 'HJJTRA_CTL called'
      write(LUPRI,*) 'KTRLVL, KWORK:',KTRLVL, KWORK

      CALL QUIT('HJJTRA_CTL not finished yet')

      end

      subroutine HJJTRA_NXTH2M(IC,ID,H2CD,NEEDTP,IDIST)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!  This version is interface routine for .HJJTRA integral transformation.
!
! Purpose:
!    Read next Mulliken two-electron integral distribution (**|cd)
!    where (cd) distribution is needed according to NEEDTP(ITYPCD)
!
! Usage:
! 1) Set IDIST = 0 before first call.
!    Do NOT change IDIST in calling routine
!    until last distribution has been read (signalled by IDIST .eq. -1)
! 2) if IDIST<0, then IC, ID must be set to correspond to the distribution
!    wanted in H2CD
!

      implicit none
#include "priunit.h"

      integer  IC, ID, NEEDTP(*), IDIST
      real*8   H2CD(NORBT,NORBT)

! Used from include files:
!  inforb.h : NORBT, ???
#include "inforb.h"

      write(LUPRI,*) 'HJJTRA_NXTH2M called, IC, ID, IDIST =',
     &   IC, ID, IDIST

      if (IDIST .le. 0) then
         if (IC .gt. 0 .and. IC .le. NORBT .and.
     &       ID .gt. 0 .and. ID .le. NORBT) then
            write (LUPRI,*) 'HJJTRA_NXTH2M: invalid input IC, ID',IC,ID
            call quit('HJJTRA_NXTH2M: invalid input IC, ID')
         end if
      else if (IDIST .eq. 0) then
      else
      end if

      end

      subroutine HJJTRA_NXTH2D(IC,ID,H2CD,NEEDTP,IDIST)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!  This version is interface routine for .HJJTRA integral transformation.
!  
! Purpose:
!    Read next Dirac two-electron integral distribution <**|cd> = (*c|*d)
!    where (cd) distribution is needed according to NEEDTP(ITYPCD)
!
!

      implicit none
#include "priunit.h"

      integer  IC, ID, NEEDTP(*), IDIST
      real*8   H2CD(NORBT,NORBT)

! Used from include files:
!  inforb.h : NORBT, ???
#include "inforb.h"

      call quit('ERROR: HJJTRA_NXTH2D is not implemented')

      end

! -- end of DALTON/sirius/sirhjjtra.F --
