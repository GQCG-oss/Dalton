      integer :: NINOB,NACOB,NDEOB,NOCOB,NTOOB
      integer :: NORB0,NORB1,NORB2,NORB3,NORB4
      integer :: NOSPIR,IOSPIR
      integer :: NINOBS,NR0OBS
      integer :: NRSOBS
      integer :: NR4OBS,NACOBS
      integer :: NOCOBS
      integer :: NTOOBS,NDEOBS,NRS4TO
      integer :: IREOTS,IREOST,ISMFTO
      integer :: ITPFSO,IBSO
      integer :: NTSOB,IBTSOB, ITSOB
      integer :: NOBPTS
      integer :: IOBPTS
      integer :: ITOOBS,ITPFTO,ISMFSO
      integer :: NOBPT,NAOS_ENV
      integer :: NMOS_ENV,I_IAD
      COMMON/ORBINP_LUCI/NINOB,NACOB,NDEOB,NOCOB,NTOOB,                 &
     &                   NORB0,NORB1,NORB2,NORB3,NORB4,                 &
     &                   NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),          &
     &                   NINOBS(MXPOBS),NR0OBS(1,MXPOBS),               &
     &                   NRSOBS(MXPOBS,3),                              &
     &                   NR4OBS(MXPOBS,MXPR4T),NACOBS(MXPOBS),          &
     &                   NOCOBS(MXPOBS),                                &
     &                   NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),  &
     &                   IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),  &
     &                   ITPFSO(MXPORB),IBSO(MXPOBS),                   &
     &                   NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),&
     &                   NOBPTS(6+MXPR4T,MXPOBS),                       &
     &                   IOBPTS(6+MXPR4T,MXPOBS),                       &
     &                   ITOOBS(MXPOBS),ITPFTO(MXPORB),ISMFSO(MXPORB),  &
     &                   NOBPT(6+MXPR4T),NAOS_ENV(MXPOBS),              &
     &                   NMOS_ENV(MXPOBS),I_IAD(MXPNGAS)
