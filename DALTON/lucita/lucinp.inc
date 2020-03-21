      INTEGER PNTGRP, NIRREP, EXTSPC, MAXML, MAXL, INTSPC, NRS0SH
      INTEGER NRSSH, MNRS1R,MXRS1R,MNRS3R,MXRS3R,LUCI_NACTEL, NSMOB
      INTEGER NRS4SH, MXR4TP, MXHR0,MXER4, NINASH, INTXCI,NDELSH
      INTEGER MNRS10,MXRS30, MNRS1RE,MXRS1RE,MNRS3RE,MXRS3RE, MNRS1ZE
      INTEGER MXRS1ZE,MNRS3ZE,MXRS3ZE, MPORENP,MPORENR,MPORENZ

      COMMON/LUCINP/PNTGRP,NIRREP,       MAXML,MAXL,                    &
     &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),                      &
     &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,LUCI_NACTEL,            &
     &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPIRR,MXPR4T),       &
     &              MXR4TP, MXHR0,MXER4,                                &
     &              NINASH(MXPIRR),                                     &
     &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30,                &
     &              MNRS1RE,MXRS1RE,MNRS3RE,MXRS3RE,                    &
     &              MNRS1ZE,MXRS1ZE,MNRS3ZE,MXRS3ZE,                    &
     &              MPORENP,MPORENR,MPORENZ
