
      module pcm_utils
              implicit none

              public ibtand
              public ibtor
              public ibtshl
              public ibtshr
              public ibtxor
              public getacord

              contains

      function ibtand(i, j)
                
        integer :: i, j
        integer ibtand
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtand = and(i, j)
#else
        ibtand = iand(i, j)
#endif

      end function
      
      function ibtor(i, j)
                
        integer :: i, j
        integer ibtor
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtor = or(i, j)
#else
        ibtor = ior(i, j)
#endif

      end function
      
      function ibtshl(i, j)
                
        integer :: i, j
        integer ibtshl
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtshl = shiftl(i, j)
#else
        ibtshl = ishft(i, j)
#endif

      end function

      function ibtshr(i, j)
                
        integer :: i, j
        integer ibtshr
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtshr = shiftr(i, j)
#else
        ibtshr = ishft(i, -j)
#endif

      end function
      
      function ibtxor(i, j)
                
        integer :: i, j
        integer ibtxor
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtxor = xor(i, j)
#else
        ibtxor = ieor(i, j)
#endif

      end function
      
      subroutine getacord(coora)      
!*****************************************************************************
!
!    getacord : Make list atomic coordinates
!
!               Written oct.2001 by Jesper Kielberg Pedersen
!               Copied here from DIRAC by Roberto Di Remigio, February 2012
!
!*****************************************************************************
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
      real(8), intent(out) :: coora(3,*)

      integer              :: jatom, icent, mulcnt, isymop
!
!     Make the full matrix of cartesian coordinates from CORD(NUCIND)
!
      jatom = 0
      do icent = 1, nucind
         mulcnt = istbnu(icent)
         if (mult(mulcnt) .eq. 1) then
            jatom = jatom + 1
            coora(1,jatom) = cord(1,icent)
            coora(2,jatom) = cord(2,icent)
            coora(3,jatom) = cord(3,icent)
        else
            do isymop = 0, maxopr
            if (iand(isymop,mulcnt) .eq. 0) then
                  jatom = jatom + 1
                  coora(1,jatom) = pt(iand(isymax(1,1),isymop))*cord(1,icent)
                  coora(2,jatom) = pt(iand(isymax(2,1),isymop))*cord(2,icent)
                  coora(3,jatom) = pt(iand(isymax(3,1),isymop))*cord(3,icent)
              end if
            enddo
        end if
      enddo 

      end subroutine
 
      end module pcm_utils
