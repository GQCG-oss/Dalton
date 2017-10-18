! file: cbihr2.h
! Control variables for her2* : AO 2-electron integrals from Hermit
      REAL*8          THRTWO, THRFAC
      INTEGER         IPRTWO, IPRNTA, IPRNTB, IPRNTC, IPRNTD,
     &                ICEDIF, IFTHRS
      LOGICAL         RUNTWO, RTNTWO, TKTIME, SOFOCK, USRSCR
      COMMON /CBIHR2/ THRTWO, THRFAC(2),
     &                IPRTWO, IPRNTA, IPRNTB, IPRNTC, IPRNTD,
     &                ICEDIF, IFTHRS,
     &                RUNTWO, RTNTWO, TKTIME, SOFOCK, USRSCR
! -- end of cbihr2.h --
