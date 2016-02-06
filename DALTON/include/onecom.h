!     FILE : onecom.h
!
!     internal information for the one-electron integral code in the her1*.F files
!
      LOGICAL         LDIAG, ONECEN, SPHRA, SPHRB, SPHRAB
      REAL*8 CORAX, CORAY, CORAZ, CORBX, CORBY, CORBZ, SIGNBX, SIGNBY,  &
     &       SIGNBZ, HKAB
      INTEGER NHKTA, NHKTB, KHKTA, KHKTB, KHKTAB, KCKTA, KCKTB, KCKTAB, &
     &        IDENA, IDENB, NCENTA, NCENTB, ICENTA, ICENTB, NUCA, NUCB, &
     &        JSTA, JSTB, MULA, MULB, MAB, JMAX, ISTEPU, ISTEPV, NAHGTF,&
     &        NUMCFA, NUMCFB
      COMMON /ONECOM/ CORAX, CORAY, CORAZ, CORBX, CORBY, CORBZ, SIGNBX, &
     &                SIGNBY, SIGNBZ, HKAB, NHKTA, NHKTB, KHKTA, KHKTB, &
     &                KHKTAB, KCKTA, KCKTB, KCKTAB, IDENA, IDENB,       &
     &                NCENTA, NCENTB, ICENTA, ICENTB, ONECEN, NUCA,     &
     &                NUCB, JSTA, JSTB, MULA, MULB, MAB, LDIAG, JMAX,   &
     &                ISTEPU, ISTEPV, NAHGTF, NUMCFA, NUMCFB, SPHRAB,   &
     &                SPHRA, SPHRB
! --- end of onecom.h ---
