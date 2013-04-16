*comdeck mxpdim.inc $Revision: 1.1.2.1 $ 
*
      INTEGER IMAXLBLKSZ
      COMMON /CMXPWRD/ MXPWRD,IMAXLBLKSZ
*
* contains all PARAMETERS defining LUCIA
      PARAMETER (MXPIRR = 8)
      PARAMETER (MXPOBS = 35)
CTF  Globally included here:
      parameter (MXTSOB = 30)
CSK      parameter (MXTSOB = 50)
*
      PARAMETER (MXPR4T = 10)
      PARAMETER (MXPR0T = 1)
      PARAMETER (MXPORB = 510)
      PARAMETER (MXPXOT = 9)
      PARAMETER (MXPXST = 100)
      PARAMETER (MXPSHL = 100)
      PARAMETER (MXPL = 20)
      PARAMETER (MXPXT = 25)
      PARAMETER (MXPICI = 10)
      PARAMETER (MXPSTT = 1000)
      PARAMETER (MXPCSM = 20)
      PARAMETER (MXPCTP = 30)
      PARAMETER (MXCNSM = 8)
      PARAMETER (MTYP = 30)
*. Note : MXPNGAS = MXPR4T+6 !!
*. Required in order to handle GAS and RAS within /LUCINP/   
      PARAMETER (MXPNGAS = 16)
      PARAMETER(MXPNSMST = 8)
*. Largest allowed division of space for perturbation operator
      PARAMETER (MXPPTSPC=20)
*. Largest cycle of spaces in multispace method
      PARAMETER(MXPCYC = 100)
*
*  For MOLUC
      parameter (MXTIT   =          10)
      parameter (MXSEQ   =           5)
      parameter (MXITEM  =          32)
*     parameter (MXROOT  =          25) ! defined in rasdim.inc
