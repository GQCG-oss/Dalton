! --- FILE: qmmm.h ---
!
      LOGICAL SPLDIP,CONMAT,FIXDIP,MMDAMP,MMPROP,                       &
     &        SPLNMR,MMMAT,MMITER,MMDIIS,LCLOSE,                        &
     &        LMMNUL,MMREST,RELMAT,LGFILE, LHFFIL,                      &
     &        ZERODI,ZEROQU,ISOPOL,NOPOL,NOMB,MMTIME,                   &
     &        NEWEXC
      INTEGER MMCENT,NMULT,IPOLTP,NEXLST,EXLIST,                        &
     &        ZEROAL,NNZAL,IDAMP,NQMNUC,MXMMIT,                         &
     &        MXMMDI,NMMAC,IQMMMC,                                      &
     &        NELEME,ELEME,IPQMMM,                                      &
     &        MXMMCT,MXEXCL
!
      PARAMETER(MXMMCT = 9000)
      PARAMETER(MXEXCL = 80)

!
      REAL*8          MMCORD, MUL0MM, MUL1MM, MUL2MM,                   &
     &                POLMM, POLIMM, THRMM, QMPOL, THMMIT,              &
     &                RQMMMC, RCUTMM,DELFLD,                            &
     &                ECHART,EDIPT,EQUADT,                              &
     &                EDELD,EDNUC,EDMULT,ENUMUL
      PARAMETER ( THRMM = 1.0D-10 )

      COMMON /REQMMM/ MMCORD(3,MXMMCT), MUL0MM(MXMMCT),                 &
     &                MUL1MM(3,MXMMCT), MUL2MM(6,MXMMCT),               &
     &                POLMM(6,MXMMCT), POLIMM(MXMMCT),                  &
     &                QMPOL(MXCENT) , THMMIT, RQMMMC, RCUTMM,           &
     &                DELFLD,ECHART,EDIPT,EQUADT,                       &
     &                EDELD,EDNUC,EDMULT,ENUMUL

      COMMON /LOQMMM/ SPLDIP,CONMAT,FIXDIP,MMDAMP,MMPROP,               &
     &                SPLNMR,MMMAT,MMITER,MMDIIS,LCLOSE,LHFFIL,         &
     &                LMMNUL(MXMMCT),MMREST,RELMAT,                     &
     &                LGFILE,ZERODI,ZEROQU,ISOPOL,NOPOL,NOMB,           &
     &                MMTIME,NEWEXC

      COMMON /INQMMM/ MMCENT,NMULT,IPOLTP,NEXLST,                       &
     &                EXLIST(MXEXCL,MXMMCT),ZEROAL(MXMMCT),             &
     &                NNZAL,IDAMP,NQMNUC,MXMMIT,MXMMDI,                 &
     &                NMMAC,IQMMMC,NELEME,ELEME(MXMMCT),IPQMMM

!C --- end of qmmm.h ---
