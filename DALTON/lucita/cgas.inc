
      integer NGAS,NGSSH, NGSOB, NGSOBT, IGSOCC, IGSINA, IGSDEL
      integer IGSOCCX, NCISPC, NCMBSPC, LCMBSPC, ICMBSPC, NMXOCCLS
      integer IPHGAS, IHPVGAS

      COMMON/CGAS/NGAS,NGSSH(MXPIRR,MXPNGAS),                           &
     &            NGSOB(MXPOBS,MXPNGAS),                                &
     &            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,      &
     &            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC,                     &
     &            NCMBSPC, LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),      &
     &            NMXOCCLS,IPHGAS(MXPNGAS),IHPVGAS(MXPNGAS)
