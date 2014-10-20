!     FILE: include/cbihrs.h
!     control common block for control of FORMSUP in her2sup.F, set under *SUPINT in input /hjaaj
      real*8  THRSUP
      integer IPRSUP
      logical RUNSUP, NOSSUP, ONLY_J
      COMMON /CBIHRS/ THRSUP, IPRSUP, RUNSUP, NOSSUP, ONLY_J
!  -- end of cbihrs.h --
