!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
module mlcc_block_import
!
!
!  mlcc types
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: import print units form Dalton
!
   integer, parameter   :: mxcorb = 2400
!  priunit, contains variables used for I/O
#if defined(__CVERSION__)
struct common_priunit {
  int lucmd, lupri, luerr, lustat, luw4, lupot, ninfo, nwarn, iprstatr;
};
extern struct common_priunit priunit_;
#else
!     FILE: priunit.h
      COMMON /PRIUNIT/ LUCMD,                                           &
     &        LUPRI, LUERR, LUSTAT, LUW4, LUPOT, NINFO, NWARN, IPRSTAT, &
     &        LQM3PCM, LPCMQM3
! ---  end of priunit.h ---
#endif
!
!     cc_tap contains information about the integral files
!
       COMMON /CC_TAP/ LUIAJB, LUBFDN, LURES, LUCSOL, LUPRPC,           &
     &                 LUMMMO, LUCCEF, LUCCPO
!
end module mlcc_block_import

