*----------------------------------------------------------------------*
* Define data type conversion constants                                *
* ItoB : integer --> byte                                              *
* RtoB : real*8  --> byte                                              *
* RtoI : real*8  --> integer                                           *
* (units=Bytes)                                                        *
*----------------------------------------------------------------------*
      Integer    ItoB,      RtoB,      RtoI

#ifdef VAR_INT64
      Parameter( ItoB = 8 , RtoB = 8 , RtoI = 1  )
#else
      Parameter( ItoB = 4 , RtoB = 8 , RtoI = 2  )
#endif
