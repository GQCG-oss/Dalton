      INTEGER MAXFRO
      INTEGER NSYM,NSYMFR,MULD2H,NLAMDS,N2BAST,                         &
     &        NRHF,NVIR,NORB,NRHFT,NVIRT,NORBT,                         &
     &        NRHFS,NVIRS,NORBS,NRHFTS,NVIRTS,NORBTS,                   &
     &        NBAS,NBAST,NNBASX,N2BASX,IORB,IBAS,                       &
     &        NRHFFR,NVIRFR,                                            &
     &        KFRRHF,KFRVIR,NFC,NFV,                                    &
     &        NRHFA,NRHFSA,NRHFB,NRHFSB,NRHFTB,                         &
     &        CCORBLAST
!
      PARAMETER (MAXFRO = 300)
!
      LOGICAL LGLO
!
      COMMON /CCORB/ NSYM,NSYMFR,MULD2H(8,8),NLAMDS,N2BAST,             &
     &               NRHF(8),NVIR(8),NORB(8),NRHFT,NVIRT,NORBT,         &
     &               NRHFS(8),NVIRS(8),NORBS(8),NRHFTS,NVIRTS,NORBTS,   &
     &               NBAS(8),NBAST,NNBASX,N2BASX,IORB(8),IBAS(8),       &
     &               NRHFFR(8),NVIRFR(8),                               &
     &               KFRRHF(MAXFRO,8),KFRVIR(MAXFRO,8),                 &
     &               NFC,NFV,LGLO,NRHFA(8),NRHFSA(8),NRHFB(8),          &
     &               NRHFSB(8),NRHFTB
!
!
      COMMON /CCORB/ CCORBLAST
      !  Very important !!!
      !  Always keep CCORBLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
