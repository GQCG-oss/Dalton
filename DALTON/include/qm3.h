! --- FILE: qm3.h ---
      LOGICAL QM3LO1, QM3LO2, LOCLAS
      LOGICAL CCMM, FIXMOM, OLDTG, ONLYOV
      LOGICAL LONUPO, LOELFD, LOSPC, LOEC3
      LOGICAL LOSHAW,REPTST,RELMOM,SLOTH
      LOGICAL LONEPAR, LTWOPAR, LEPSADD, LSIGADD
      LOGICAL SKIPNC, VDWSKP, MYITE, MYMAT, EXPON
      LOGICAL PRFQM3, INTDIR, FORQM3, REDCNT, LGSPOL, RUNQM3
      LOGICAL QMDAMP, NYQMMM, HFFLD, CCFIXF, FFIRST
      LOGICAL MMPCM, LADDMM, FIRST1
!
! ---------------------------------------------------------
! In the present implementation the MXQM3 parameter follows
! the MXCENT_QM parameter in the mxcent.h include file. This is
! crucial for this implementation to work properly!!
! ---------------------------------------------------------
!
      INTEGER ISUBSY,ISUBSI,MXTYP1,NSYSBG,NSYSED
      INTEGER NSISY, ISYTP, NTOTQM3, IQM3PR, ICHRGS
      INTEGER MXDIIT, NUSITE, MXQM3, MXTYPE, NCOMS
      INTEGER NTOTIN, NUALIS, NQMBAS,NMMBA1, NREPMT
      INTEGER ISIGEPS, NSIGEPS, MXQ, NSTATES, ICQM3
      INTEGER NOSIMOLD, NOSIMFIRST, MXITMP
!
      PARAMETER(NMMBA1 = 5000)
      PARAMETER(MXQM3  = 500) ! should be equal to MXCENT_QM in include/mxcent.h
      PARAMETER(MXTYPE = 20)
      PARAMETER(MXQ = MXQM3)
      PARAMETER(NSTATES = 120)
!
      CHARACTER MDLWRD*7
!
      LOGICAL SHAWFC(0:MXTYPE)
      LOGICAL RDFILE(0:MXTYPE), DISMOD(0:MXTYPE)
!
!     ----------------------------------------------
!     IQM3PR takes the role of the IPREAD print flag
!     used in herrdn.F!
!     ----------------------------------------------
!
      REAL*8  QM3CHG,QM3LJA,QM3LJB,ALPIMM
      REAL*8  ALTXX,ALTXY,ALTXZ,ALTYY,ALTYZ,ALTZZ
      REAL*8  ECLPOL,ECLVDW,ECLQM3
      REAL*8  THDISC,ENUQM3,CHAOLD
      REAL*8  EMMPOL,EMMVDW,EMMELC,EMM_MM,EVDWSH,PEDIP1
      REAL*8  ENSQM3,EPOQM3
      REAL*8  XMMQ,  YMMQ,  ZMMQ,  MMQ
      REAL*8  XMMMY, YMMMY, ZMMMY, MMMYX, MMMYY, MMMYZ
      REAL*8  THRSMP,DMMSAVE
      REAL*8  QMCOM, ADAMP

      COMMON /REAQM3/ THDISC,ECLPOL,ECLVDW,ECLQM3,ENUQM3,               &
     &                EMMPOL,EMMVDW,EMMELC,EMM_MM,EVDWSH,               &
     &                PEDIP1,ENSQM3,EPOQM3,THRSMP,DMMSAVE,              &
     &                QMCOM(3),ADAMP

      COMMON /LOGQM3/ RDFILE,DISMOD,QM3LO1,QM3LO2,CCMM,FIXMOM,          &
     &                OLDTG,ONLYOV,LONUPO,LOELFD,LOSPC,LOEC3,NYQMMM,    &
     &                SHAWFC,LOSHAW,REPTST,RELMOM,SLOTH,HFFLD,CCFIXF,   &
     &                LONEPAR,LTWOPAR,LEPSADD,LSIGADD,LOCLAS,           &
     &                SKIPNC,VDWSKP,MYITE,MYMAT,EXPON,PRFQM3,FFIRST,    &
     &                INTDIR, FORQM3, REDCNT, LGSPOL, MMPCM,            &
     &                LADDMM, FIRST1, RUNQM3, QMDAMP

      COMMON /INTQM3/ IQM3PR,ISYTP,NTOTQM3,NUSITE,NCOMS,NTOTIN,         &
     &                MXDIIT, NQMBAS, NREPMT, NSIGEPS, NOSIMOLD,        &
     &                NOSIMFIRST, MXITMP

      COMMON /QM3WRD/ MDLWRD(0:MXTYPE)

      COMMON /QM3GNR/ ISUBSY(MXQM3),ISUBSI(MXQM3),                      &
     &                NSYSBG(0:MXTYPE),NSYSED(0:MXTYPE),                &
     &                NSISY(0:MXTYPE),                                  &
     &                ICHRGS(0:MXTYPE),NUALIS(0:MXTYPE),                &
     &                ISIGEPS(0:MXTYPE),ICQM3(NSTATES)

      COMMON /QM3SYS/ QM3CHG(0:MXTYPE,MXQM3),                           &
     &                QM3LJA(0:MXTYPE,0:MXTYPE),                        &
     &                QM3LJB(0:MXTYPE,0:MXTYPE),                        &
     &                ALPIMM(0:MXTYPE,MXQM3),CHAOLD(MXQM3),             &
     &                ALTXX(0:MXTYPE),ALTXY(0:MXTYPE),                  &
     &                ALTXZ(0:MXTYPE),ALTYY(0:MXTYPE),                  &
     &                ALTYZ(0:MXTYPE),ALTZZ(0:MXTYPE),                  &
     &                XMMQ(MXQ),YMMQ(MXQ),ZMMQ(MXQ),MMQ(MXQ),           &
     &                XMMMY(MXQ),YMMMY(MXQ),ZMMMY(MXQ),                 &
     &                MMMYX(MXQ),MMMYY(MXQ),MMMYZ(MXQ)
! --- end of qm3.h ---
