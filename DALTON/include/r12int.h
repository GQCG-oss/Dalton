!     File: r12int.h
!
!     The next parameters indicates the number of integer and
!     logical variables to be transferred in a parallel calculation.
!     Must be updated if COMR12 is changed.
!
      INTEGER NR12I, NR12L, IANR12, IAPR12, MAXGAM, XMIADR
      PARAMETER (NR12I = 30, NR12L = 40)
      PARAMETER (MAXGAM = 50)

      REAL*8  GAMAC,  GAMAD,  GAMAB,  GAMAX
!
      LOGICAL R12INT, R12TRA, U12INT, U21INT, R12SQR,                   &
     &        R12CAL, R12NOA, R12NOP, R12NOB, R12OLD, NORXR,            &
     &        R12HYB, COMBSS, R12EIN, R12EOR, R12ECO, R12XXL,           &
     &        V12INT, NOTR12, NOTU12, NOTV12, XXXXXX, CC2R12INT,        &
     &        ANTICO, DIRSCF, U12DIR, R12DIR, LOOPDP, ONEAUX,           &
     &        V12DIR, STEPIV, DCCR12, CUSP12, DUMR12, CCPAIR,           &
     &        NOTONE, NOTTWO, CCSDR12INT,     R12CBS, NOBP,             &
     &        SLATER, MIAU12, AENONB, NOTTRE, CCR12SM           
      LOGICAL R12DIA, R12SVD, R12RST, R12PRP, NATVIR, USEVABKL
!
      INTEGER NTOGAM, NRXR12(8), IANCC2, IAPCC2, LOCFRO(8)
      INTEGER INTGAC, INTGAD, IADV12, IADR12, IADU12, IADU21, NOPP12,   &
     &        COMR12LAST
!
      COMMON /COMR12/ GAMAC,  GAMAD, GAMAB(MAXGAM), GAMAX(MAXGAM),      &
     &        R12INT, R12TRA, U12INT,                                   &
     &        U21INT, R12SQR, R12CAL, R12NOA, R12NOP, R12NOB,           &
     &        R12OLD, NORXR,  R12HYB, COMBSS, R12EIN,                   &
     &        R12EOR, R12ECO, R12XXL, V12INT, NOTR12, NOTU12,           &
     &        NOTV12, ANTICO, DIRSCF, U12DIR, R12DIR, V12DIR,           &
     &        STEPIV, DCCR12, CUSP12, DUMR12,                           &
     &        INTGAC, INTGAD, IADV12, IADR12, IADU12, IADU21, NOPP12,   &
     &        CCPAIR, XXXXXX, CC2R12INT, LOOPDP, ONEAUX,                &
     &        LU21INT, IOFFU21,                                         &
     &        IANR12, IAPR12, NOTONE, NOTTWO,CCSDR12INT, R12CBS, NOBP,  &
     &        XMIADR, MIAU12, AENONB, NRXR12, NOTTRE, IANCC2, IAPCC2,   &
     &        NTOGAM, SLATER, LOCFRO, CCR12SM,                          &
     &        R12DIA, R12SVD, R12RST, R12PRP, NATVIR, USEVABKL  
!
      COMMON /COMR12/ COMR12LAST
      !  Very important !!!
      !  Always keep COMR12LAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
!
!
!
!
!
!
!
      CHARACTER*10    FNVAJKL
      COMMON /R12FNS/ FNVAJKL
!
!
!
!
!
      REAL*8          VCLTHR, SVDTHR, R12LEV
      COMMON /COMVCL/ VCLTHR, SVDTHR, R12LEV
!
!
!
!
!
      INTEGER MBAS1(8), MBAS2(8), MBAS1T, MBAS2T, LU21INT, IOFFU21,     &
     &        NORB1(8), NORB2(8), MBSMAX, NPLSH_R12, CMMMULLAST
      LOGICAL LAUXBS, NOAUXB
      COMMON /CMMMUL/ MBAS1, MBAS2, NORB1, NORB2, MBAS1T, MBAS2T,       &
     &        MBSMAX, NPLSH_R12,                                        &
     &        LAUXBS, NOAUXB
      COMMON /CMMMUL/ CMMMULLAST 
      !  Very important !!!
      !  Always keep CMMMULLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
!
!
!
!
!
!
      CHARACTER*8 LABEL
      COMMON /CMMOLL/ LABEL
      REAL*8  BRASCL,KETSCL
      COMMON /R12SCALE/ BRASCL, KETSCL
! -- end of r12int.h --
