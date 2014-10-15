C --- File: inftra.h ---
C     Global information for the sirius/sirtra.F 2-el integral transforamtion modules
C
      REAL*8  THRP,   THRQ
      INTEGER LSRTAO, IPRTRA, LBUF, NIBUF, NBITS
      LOGICAL USEDRC, HJJTRA, NEWTRA, USE_INTSORT
      COMMON /INFTRA/ THRP,   THRQ,
     &                LSRTAO, IPRTRA, LBUF, NIBUF, NBITS,
     &                USEDRC, HJJTRA, NEWTRA, USE_INTSORT
C
C     Name of file with 2-el integrals in AO basis.
C     Defaults to "AOTWOINT", but can be set to e.g. "AOSR2INT" to transform
C     short range integrals to MO basis.
C
      CHARACTER*8 AO2INTFILE_LABEL
      COMMON /INFTRA_CHAR/ AO2INTFILE_LABEL
C --- end of inftra.h ---
