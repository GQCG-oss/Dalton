C file: include/spinfo.h
C
C MULTS : spin multiplicity
C MS2   : 2 * M_S value of determinants
C
      PARAMETER(MTYP = 30 )
      COMMON/SPINFO/MULTS,MS2,
     &              MINOP,MAXOP,NTYP,NDTFTP(MTYP),NCSFTP(MTYP),
     &              NCNFTP(MTYP,8)
C end of include/spinfo.h
