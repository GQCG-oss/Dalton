* ==================================================================== *
* ENERGY_DEC.H
* -------------------------------------------------------------------- *
* explicite variable declaration for ENERGY.H
* -------------------------------------------------------------------- *
#if defined (SYS_CRAY)
      REAL ENERKE, ENERNA, ENEREE, ENERNN
      REAL GRADKE, GRADNA, GRADEE, GRADNN, GRADFS
#else
      DOUBLE PRECISION ENERKE, ENERNA, ENEREE, ENERNN
      DOUBLE PRECISION GRADKE, GRADNA, GRADEE, GRADNN, GRADFS
#endif
* -------------------------------------------------------------------- *
