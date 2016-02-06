      REAL*8  NPCORD,  mm_cord,  NPCHRG,  MMCHRG,                       &
     &        NPFPOL,  NPFCAP,  NPFOMG1, NPFGAM1,                       &
     &        NPFOMG2, NPFGAM2, NPFFAC,                                 &
     &        MMFM0,   MMFPOL,                                          &
     &        ENSOLQNP, EESOLQNP, ENSOLMNP, EESOLMNP,                   &
     &        ENSOLQMM, EESOLQMM, ENSOLMMM, EESOLMMM
!
      INTEGER TNPBLK,  TMMBLK,  IPRTLVL,                                &
     &        TNPATM,  TMMATM,  NPATOM, NPFTYP,                         &
     &        MMATOM,  MMFTYP,  TNPFF,  TMMFF,                          &
     &        MMMOL,   TPOLATM, MMSKIP
!
      LOGICAL DONPSUB, DOMMSUB, NPMQGAU,                                &
     &        MQITER,           DONPCAP, DOMMCAP,                       &
     &        DONPPOL, DOMMPOL, NOVDAMP
!
      INTEGER MAXBLK,  MXNPATM, MXMMATM, MXNPFF, MXMMFF
      PARAMETER (MAXBLK = 5)
      PARAMETER (MXNPATM = 10000)
      PARAMETER (MXMMATM = 90000)
      PARAMETER (MXNPFF = 5)
      PARAMETER (MXMMFF = 20)
!
      COMMON /QMNPIN/ NPCORD(3,MXNPATM), mm_cord(3,MXMMATM),            &
     &                NPCHRG(MAXBLK),    MMCHRG(MAXBLK),                &
     &                NPFPOL(MXNPFF),    NPFCAP(MXNPFF),                &
     &                NPFOMG1(MXNPFF),   NPFGAM1(MXNPFF),               &
     &                NPFOMG2(MXNPFF),   NPFGAM2(MXNPFF),               &
     &                NPFFAC(MXNPFF),                                   &
     &                MMFM0(MXMMFF),     MMFPOL(MXMMFF),                &
     &                ENSOLQNP, EESOLQNP, ENSOLMNP, EESOLMNP,           &
     &                ENSOLQMM, EESOLQMM, ENSOLMMM, EESOLMMM,           &
     &                TNPBLK,  TMMBLK,   IPRTLVL,                       &
     &                TNPATM,  TMMATM,   TNPFF,   TMMFF,                &
     &                NPFTYP(MXNPATM),   TPOLATM, MMFTYP(MXMMATM),      &
     &                NPATOM(MAXBLK),    MMATOM(MAXBLK),                &
     &                MMMOL(MXMMATM),    MMSKIP(MXMMATM),               &
     &                DONPSUB, DOMMSUB,  NPMQGAU,                       &
     &                MQITER,            DONPCAP, DOMMCAP,              &
     &                DONPPOL, DOMMPOL,  NOVDAMP
