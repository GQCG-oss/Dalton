! FILE: cbiwlk.h
!
! Info for geometry walk (abawalk.F)
!
! IWKTYP identifies which walk type (copied from abawalk.F):
!     DATA WLKTYP /'trust-region minimization                 ', ! 0
!    &             'mode-following                            ', ! 1
!    &             'gradient extremal                         ', ! 2
!    &             'dynamic walk                              ', ! 3
!    &             'Newton step                               ', ! 4
!    &             'eigenvector step                          ', ! 5
!    &             'numerical differentiation                 ', ! 6
!    &             'IRC path                                  ', ! 7
!    &             'TS by minimization on image surface       '/ ! 8

      PARAMETER (MAXTMP = 20)
      LOGICAL WFPRED, REJECT, KEEPSY, VIBCNV,
     &        START, DOREPW, IMAGE, STRICT, NOORTH, NATCON, V3CAL,
     &        VIBAVE, NMODIF, ECKART, DOTEMP, DOCENT, REUSED, ACCURT,
     &        HESFLW
      COMMON /CBIWLK/ TOLST, TRUSTR, TRUSTI, TRUSTD,                      ! real*8
     &        REJMIN, REJMAX, RTRMIN, RTRGOD,
     &        XMXNUC, ZERGRD, DISPLC,
     &        STRMOM(3*MXCENT), ECKGEO(3,MXCOOR),
     &        SCALCO(3,MXCENT), THRLDP, ANHFAC, TRUMAX,
     &        TEMP(MAXTMP),
     &        ISTMOM(3*MXCENT), NSTMOM, NTEMP,                            ! integer
     &        IPRWLK, IWKTYP, IWKIND, IMODE, ISCTYP,
     &        IPART(MXCENT), NZEROG, IZEROG(MXCOOR),
     &        WFPRED, REJECT, KEEPSY, VIBCNV,                             ! logical
     &        START, DOREPW(0:7), IMAGE, STRICT, NOORTH, NATCON, V3CAL,
     &        VIBAVE, NMODIF, ECKART, DOTEMP, DOCENT, REUSED, ACCURT,
     &        HESFLW
