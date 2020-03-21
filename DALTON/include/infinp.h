! File : infinp.h
!
!  -*- mode:fortran; fortran-continuation-string: "&" -*-
!     this common block contains general Sirius input read in sirius/sirinp.F
!     (specified under **WAVE FUNCTION).  /hjaaj Oct 2003
!     - SCF specific input is in scbrhf.h
!     - orbital specifications are in inforb.h
!
!  MCTYPE     Calculation Type
!   -1     HSROHF
!    0     single configuration (SCF)
!    1     CASSCF
!    2     RASSCF
!   11     APSG
!
      INTEGER         NFLAG, MAXRTS, MXFELT
      PARAMETER (NFLAG = 80, MAXRTS = 150, MXFELT = 20) ! MAXRTS ==> max. number of roots; has to be equal to MAXRTS_LU in lucita/cstate.inc
!
      INTEGER         NFIELD, ISPIN,NMCAVER,ISTATE,LSYM,NACTEL, MCTYPE, &
     &                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,        &
     &                LROOTS,NROOTS,IROOT,                              &
     &                NOROT        ,IMOORD,                             &
     &                IORTO,ICI0,KDEL,ICHECK,NTIT,                      &
     &                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,  &
     &                ITRLVL,ITRFIN,JCHSYM,JCHORB,                      &
     &                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO, N_in_RN,    &
     &                SPINDENS_lvl
      COMMON /INTINP/ NFIELD, ISPIN,NMCAVER,ISTATE,LSYM,NACTEL, MCTYPE, &
     &                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,        &
     &                LROOTS,NROOTS,IROOT(MAXRTS),                      &
     &                NOROT(MXCORB),IMOORD(MXCORB),                     &
     &                IORTO,ICI0,KDEL,ICHECK,NTIT,                      &
     &                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,  &
     &                ITRLVL,ITRFIN,JCHSYM,JCHORB,                      &
     &                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO, N_in_RN,    &
     &                SPINDENS_lvl
!
      LOGICAL         FLAG,        DOSCF,DOMP2,DOCINO,DOCI,DOMC,DORSP,  &
     &                FCVORB,LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP, &
     &                JOLSEN,ABAIPH,INERSI,INERSF,DODFT,DONEVPT,HSROHF, &
     &                BOYORB,PIPORB,DOFCI,DOCISD,DOLUCITA,DOMEP,        &
     &                DO_CUBE, DOAPSG, DO_VIRTRUNC
!     variables for srDFT /hjaaj
      LOGICAL         DOCISRDFT,DOHFSRDFT,DOMCSRDFT,ADDSRI,SRHYBR
      COMMON /LOGINP/ FLAG(NFLAG), DOSCF,DOMP2,DOCINO,DOCI,DOMC,DORSP,  &
     &                FCVORB,LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP, &
     &                JOLSEN,ABAIPH,INERSI,INERSF,DODFT,DONEVPT,HSROHF, &
     &                BOYORB,PIPORB,DOFCI,DOCISD,DOLUCITA,DOMEP,        &
     &                DOCISRDFT,DOHFSRDFT,DOMCSRDFT,ADDSRI,SRHYBR,      &
     &                DO_CUBE, DOAPSG, DO_VIRTRUNC
      LOGICAL         LOGINPlast
      COMMON /LOGINP/ LOGINPlast ! for MPI and getbytespan
      LOGICAL         SUPSYM
      EQUIVALENCE (SUPSYM,FLAG(17))
!
      REAL*8          SPIN, POTNUC, EPSOL,EPSTAT,EPPN,RSOL,             &
     &                THRGRD, THRPWF, THRCI, THRMC, THRCGR,             &
     &                EFIELD, CMAXMO, THROVL,                           &
     &                THRSSY, DEFLVL, WEIGHT_MCAVER, THR_VIRTRUNC
!     variables for srDFT /hjaaj
      REAL*8          THRCIDFT
      COMMON /RELINP/ SPIN, POTNUC, EPSOL,EPSTAT,EPPN,RSOL(3),          &
     &                THRGRD, THRPWF, THRCI, THRMC, THRCGR,             &
     &                EFIELD(MXFELT), CMAXMO, THROVL,                   &
     &                THRSSY, DEFLVL, WEIGHT_MCAVER(MAXRTS),            &
     &                THR_VIRTRUNC, THRCIDFT
!
      character (len=60) :: TITLE(6)
      character (len=72) :: TITMOL(2)
      character (len= 9) :: ci_program
      character (len= 4) :: CENT, TYPE
      character (len= 8) :: LFIELD
      COMMON /CHRINP/ TITLE, TITMOL, ci_program,                        &
     &                CENT(MXCORB), TYPE(MXCORB), LFIELD(MXFELT)
!-- end of infinp.h --
