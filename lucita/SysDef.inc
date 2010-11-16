*----------------------------------------------------------------------*
* Define data type conversion constants                                *
* ItoB : integer --> byte                                              *
* RtoB : real*8  --> byte                                              *
* RtoI : real*8  --> integer                                           *
* (units=Bytes)                                                        *
*----------------------------------------------------------------------*
      Integer    ItoB,      RtoB,      RtoI
      Character Arch_Type*13

#if defined(INT_STAR8)
      Parameter( ItoB = 8 , RtoB = 8 , RtoI = 1  )
#else
      Parameter( ItoB = 4 , RtoB = 8 , RtoI = 2  )
#endif

      Data Arch_Type/'IRIX64'/
