!
! File: molinp.h
!
! Contains a copy of the .mol input file in MLINE(1:NMLINE).
! This is done by RDLINE in "abacus/herrdn.F".
! NCLINE(i) points to the line defining atom center no. i
! NMLINE_1  points to "line1" in the .mol input
! NMLINE_basis points to the line where basis set is defined when BASIS in "line1"
! NMLINE_4  points to the line in which it is defined if coordinates
!           are in AU (bohr) or in Angstroem ("line4")
! NONTYP_QM number of QM types for QM3 model, value found by Master in abacus/herrdn.F,
!           used by slaves to only keep the QM centers in their setup (see abacus/herpar.F)
!
      INTEGER    IMLINE, len_MLINE
      PARAMETER (KMLINE = 2500, len_MLINE = 80)
!     note: line length of 80 is hardcoded many places so it is a lot of work to change it!
      CHARACTER*(len_MLINE) MLINE
      COMMON /MOLINC/ MLINE(KMLINE)
      INTEGER         NMLINE,NCLINE        ,NMLINE_1,NMLINE_basis,
     &                NMLINE_4,NONTYP_QM
      COMMON /MOLINP/ NMLINE,NCLINE(MXCENT),NMLINE_1,NMLINE_basis,
     &                NMLINE_4,NONTYP_QM
! -- End of molinp.h --
