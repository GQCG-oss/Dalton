      PARAMETER (MAXFRO = 60)
      LOGICAL FRORHF,FROVIR
      COMMON /CCORB/ NSYM,NSYMFR,MULD2H(8,8),NLAMDS,N2BAST,
     *               NRHF(8),NVIR(8),NORB(8),NRHFT,NVIRT,NORBT,
     *               NRHFS(8),NVIRS(8),NORBS(8),NRHFTS,NVIRTS,NORBTS,
     *               NBAS(8),NBAST,NNBASX,N2BASX,IORB(8),IBAS(8),
     *               NRHFFR(8),NVIRFR(8),
     *               KFRRHF(MAXFRO,8),KFRVIR(MAXFRO,8),
     *               FRORHF(MAXFRO,8),FROVIR(MAXFRO,8),
     *               NFC,NFV
