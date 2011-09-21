!
! File: molinp.h
!
! Contains a copy of the .mol input file in MLINE(1:NMLINE).
! This is done by RDLINE in "abacus/herrdn.F".
! NCLINE(i) points to the line defining center no. i
! NMLAU     points to the line in which it is defined if coordinates
!           are in AU (bohr) or in Angstroem
! NONTYP_QM number of QM types, found by Master in abacus/herrdn.F, used by slaves
!           to only keep the QM centers in their setup (see abacus/herpar.F)
!
      PARAMETER (KMLINE = 2500, max_len_MLINE = 100)
      CHARACTER*(max_len_MLINE) MLINE
      COMMON /MOLINC/ MLINE(KMLINE)
      COMMON /MOLINP/ NMLINE,NCLINE(MXCENT),NMLAU,NONTYP_QM
! -- End of molinp.h --
