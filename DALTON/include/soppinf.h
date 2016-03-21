      LOGICAL AOSOP, AORPA, AOHRP, SOPCHK, AOTEST, DCRPA, AOSOC
      LOGICAL TRIPLET
      CHARACTER*7 FNTR1E, FNTR1D, FNTR2E, FNTR2D,
     &            FNRS1E, FNRS1D, FNRS2E, FNRS2D,
     &            FNRO1E, FNRO1D,
     &            FNDIAG, FNDENS, FNFOCK,
     &            FNSAI1, FNSAI2, FNSDA1, FNSDA2,
     &            FNBT1E, FNBT1D, FNBJ1E, FNBJ1D,
     &            FNRI1E, FNRI1D, FNRI2E, FNRI2D,
     &            FNSC1E, FNSC1D, FNSC2E, FNSC2D,
     &            FNSV1E, FNSV1D, FNSV2E, FNSV2D,
     &            FNRV1E, FNRV1D, FNRV2E, FNRV2D,
     &            FNGPVE, FNSOLA, FNGPV1, FNGPV2
      PARAMETER (FNTR1E = 'SO_TR1E', FNTR1D = 'SO_TR1D',
     &           FNTR2E = 'SO_TR2E', FNTR2D = 'SO_TR2D',
     &           FNRS1E = 'SO_RS1E', FNRS1D = 'SO_RS1D',
     &           FNRS2E = 'SO_RS2E', FNRS2D = 'SO_RS2D',
     &           FNRO1E = 'SO_RO1E', FNRO1D = 'SO_RO1D',
     &           FNDIAG = 'SO_DIAG', FNDENS = 'SO_DENS',
     &           FNFOCK = 'SO_FOCK',
     &           FNSAI1 = 'SO_SAI1', FNSAI2 = 'SO_SAI2',
     &           FNSDA1 = 'SO_SDA1', FNSDA2 = 'SO_SDA2',
     &           FNBT1E = 'SO_BT1E', FNBT1D = 'SO_BT1D',
     &           FNBJ1E = 'SO_BJ1E', FNBJ1D = 'SO_BJ1D',
     &           FNRI1E = 'SO_RI1E', FNRI1D = 'SO_RI1D',
     &           FNRI2E = 'SO_RI2E', FNRI2D = 'SO_RI2D',
     &           FNSC1E = 'SO_SC1E', FNSC1D = 'SO_SC1D',
     &           FNSC2E = 'SO_SC2E', FNSC2D = 'SO_SC2D',
     &           FNSV1E = 'SO_SV1E', FNSV1D = 'SO_SV1D',
     &           FNSV2E = 'SO_SV2E', FNSV2D = 'SO_SV2D',
     &           FNRV1E = 'SO_RV1E', FNRV1D = 'SO_RV1D',
     &           FNRV2E = 'SO_RV2E', FNRV2D = 'SO_RV2D',
     &           FNGPVE = 'SO_GPVE', FNSOLA = 'SO_SOLA',
     &           FNGPV1 = 'SO_GPV1', FNGPV2 = 'SO_GPV2')
      PARAMETER (LSOTIM = 43)
      PARAMETER (LSOWTI = 4)
      PARAMETER (LORWCI = 4)
      PARAMETER (LSOSUB = 100)
C
      CHARACTER*16 SOSUBNM
C
      INTEGER    LUTR1E, LUTR1D, LUTR2E, LUTR2D,
     &           LURS1E, LURS1D, LURS2E, LURS2D,
     &           LURO1E, LURO1D,
     &           LUDIAG, LUDENS, LUFOCK,
     &           LUSAI1, LUSAI2, LUSDA1, LUSDA2,
     &           LUBT1E, LUBT1D, LUBJ1E, LUBJ1D,
     &           LURI1E, LURI1D, LURI2E, LURI2D,
     &           LUSC1E, LUSC1D, LUSC2E, LUSC2D,
     &           LUSV1E, LUSV1D, LUSV2E, LUSV2D,
     &           LURV1E, LURV1D, LURV2E, LURV2D,
     &           LUGPVE, LUSOLA, LUGPV1, LUGPV2,
     &           SOMEMO
C
C
      INTEGER SOPPINFlast, SOPPEXClast, RWINFlast
C
CSPAS:10/08-09 including triplet 2p2h vectors
      COMMON /SOPPINF/ SOTIME(LSOTIM), SOWTIM(LSOWTI),
     &                 SOSUBNM(LSOSUB), SOMEMO(LSOSUB), ISOSUB, LWTOTAL,
     &                 SOPCHK,IIJDEN(8,8),
     &                 NIJDEN(8),IABDEN(8,8),NABDEN(8),IAIDEN(8,8),
     &                 NAIDEN(8),N2P2HOP(8),
     &                 NT2AMT1(8), NT2AMT2(8), NT2AMT3(8), NT2AMTT(8),
     &                 IT2AMT1(8,8), IT2AMT2(8,8), IT2AMT3(8,8),
     &                 ISOO(8,8), ITOO(8,8), ISVV(8,8), ITVV(8,8),
     &                 IPRSOP,AOTEST,NSAVMX, NSAVMXORIG,
     &                 AORPA,  AOHRP,  DCRPA,  AOSOP,  AOSOC,
     &                 TRIPLET,
     &                 LUTR1E, LUTR1D, LUTR2E, LUTR2D,
     &                 LURS1E, LURS1D, LURS2E, LURS2D,
     &                 LURO1E, LURO1D,
     &                 LUDIAG, LUDENS, LUFOCK,
     &                 LUSAI1, LUSAI2, LUSDA1, LUSDA2,
     &                 LUBT1E, LUBT1D, LUBJ1E, LUBJ1D,
     &                 LURI1E, LURI1D, LURI2E, LURI2D,
     &                 LUSC1E, LUSC1D, LUSC2E, LUSC2D,
     &                 LUSV1E, LUSV1D, LUSV2E, LUSV2D,
     &                 LURV1E, LURV1D, LURV2E, LURV2D,
     &                 LUGPVE, LUSOLA, LUGPV1, LUGPV2,
     &   SOPPINFlast !  Very important:
      !  Always keep SOPPINFlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
C
      COMMON /SOPPEXC/ NEXCI2(8),THREX2,
     &   SOPPEXClast !  Very important:
      !  Always keep SOPPEXClast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
C
C
      COMMON /RWINF/ SOORWC(LORWCI),
     &   RWINFlast !  Very important:
      !  Always keep RWINFlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

! -- end of soppinf.h
