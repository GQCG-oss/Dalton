! FILE: cbihr2.h

! Control variables for her2* : AO 2-electron integrals from Hermit

      REAL*8          THRFAC,    THRTWO
      COMMON /CBRHR2/ THRFAC(2), THRTWO
      INTEGER         IPRTWO, IPRNTA, IPRNTB, IPRNTC, IPRNTD,
     &                ICEDIF, IFTHRS
      COMMON /CBIHR2/ IPRTWO, IPRNTA, IPRNTB, IPRNTC, IPRNTD,
     &                ICEDIF, IFTHRS
      LOGICAL         RUNTWO, RTNTWO, TKTIME, SOFOCK, USRSCR
      COMMON /CBIHR2/ RUNTWO, RTNTWO, TKTIME, SOFOCK, USRSCR
