

      integer MXPWRD
      COMMON /CMXPWRD/ MXPWRD

!     contains all PARAMETERS defining LUCIA
      integer, parameter :: MXPIRR =    8
!     integer, parameter :: MXPOBS =   35 ! old - with only D2h + subgroups the value below should be fine
      integer, parameter :: MXPOBS =    8
!    
!     max. number of orbs per point group irrep and GA space
!smaa integer, parameter :: MXTSOB =   32
      integer, parameter :: MXTSOB =   60

      integer, parameter :: MXPR4T =   10
      integer, parameter :: MXPR0T =    1
      integer, parameter :: MXPORB =  510
      integer, parameter :: MXPXOT =    9
      integer, parameter :: MXPXST =  100
      integer, parameter :: MXPSHL =  100
      integer, parameter :: MXPL   =   20
      integer, parameter :: MXPXT  =   25
      integer, parameter :: MXPICI =   10
      integer, parameter :: MXPSTT = 1000
      integer, parameter :: MXPCSM =   20
      integer, parameter :: MXPCTP =   30
      integer, parameter :: MXCNSM =    8
      integer, parameter :: MTYP   =   30
!     MXPNGAS = MXPR4T+9
!     max number of GA spaces
      integer, parameter :: MXPNGAS  =  16
      integer, parameter :: MXPNSMST =   8

!     Largest allowed division of space for perturbation operator
      integer, parameter :: MXPPTSPC =  20

!     Largest cycle of spaces in multispace method
      integer, parameter :: MXPCYC =   100
!
      integer, parameter :: MXTIT  =    10
      integer, parameter :: MXSEQ  =     5
      integer, parameter :: MXITEM =    32
