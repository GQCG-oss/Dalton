! FILE: energy.h for abacus, individual contributions to energy and gradient
!
      REAL*8 ENERKE, ENERNA, ENEREE, ENERNN
      REAL*8 GRADKE, GRADNA, GRADEE, GRADNN
      REAL*8 GRADFS, GRADFT, GRADFTD, PEGRAD

      COMMON /CB_ABACUS_ENERGY/
     &   ENERKE, GRADKE(MXCOOR),
     &   ENERNA, GRADNA(MXCOOR),
     &   ENEREE, GRADEE(MXCOOR),
     &   ENERNN, GRADNN(MXCOOR),
     &           GRADFS(MXCOOR),
     &           GRADFT(MXCOOR),
     &           GRADFTD(MXCOOR),
     &           PEGRAD(MXCOOR)
!  -- end of energy.h --
