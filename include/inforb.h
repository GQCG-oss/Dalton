C
C     Parameters NINFI must be updated after changes (for parallelization)
C
C     NOTE: Reals and logicals should appear at the end.
C
      PARAMETER (NINFI = 333)
      COMMON /INFORB/ MULD2H(8,8), NRHF(8),NVIR(8), NFRO(8),
     *       NISH(8),NASH(8),NSSH(8),NOCC(8),NORB(8),NBAS(8),
     *       NNORB(8),NNBAS(8), N2ORB(8),N2BAS(8),
     *       IISH(8),IASH(8),ISSH(8),IOCC(8),IORB(8),IBAS(8),
     *       IIISH(8),IIASH(8),IIORB(8),IIBAS(8),I2ORB(8),I2BAS(8),
     *       ICMO(8), NSYM,
     *       NISHT,NASHT,NSSHT,NOCCT,NORBT,NBAST,NCMOT,NRHFT,NVIRT,
     *       N2ISHX,NNASHX,N2ASHX,NNASHY,NNOCCX,N2OCCX,
     *       NNORBT,NNORBX,N2ORBT,N2ORBX,NNBAST,N2BAST,NNBASX,N2BASX,
     *       NNRHFT,NNRHFX,N2RHFT,N2RHFX,NNVIRT,NNVIRX,N2VIRT,N2VIRX,
     *       NAS1(8),NAS2(8),NAS3(8),NNOCCT,N2OCCT,
     *       NAS1T,NAS2T,NAS3T
C     MXSSYM = maximum number of "super symmetries"
      PARAMETER ( MXSSYM = 100 )
      COMMON /INFOSS/ NSSYM, NORBSS(MXSSYM), IORBSS(MXSSYM),
     *       NINFSS(MXSSYM,3), MXDGSS
