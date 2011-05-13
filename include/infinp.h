! File : infinp.h
!
!  -*- mode:fortran; fortran-continuation-string: "&" -*-
!     this common block contains general Sirius input read in sirius/sirinp.F
!     (specified under **WAVE FUNCTION).  /hjaaj Oct 2003
!     - SCF specific input is in scbrhf.h
!     - orbital specifications are in inforb.h
!
      INTEGER         NFLAG, MAXRTS, MXFELT
      PARAMETER (NFLAG = 80, MAXRTS = 100, MXFELT = 20)
!
      INTEGER         NFIELD, ISPIN,NMCAVER,ISTATE,LSYM,NACTEL, MCTYPE, &
     &                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,        &
     &                LROOTS,NROOTS,IROOT,                              &
     &                NOROT        ,IMOORD,                             &
     &                IORTO,ICI0,KDEL,ICHECK,NTIT,                      &
     &                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,  &
     &                ITRLVL,ITRFIN,JCHSYM,JCHORB,                      &
     &                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO
      COMMON /INTINP/ NFIELD, ISPIN,ISTATE,LSYM,NACTEL, MCTYPE,         &
     &                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,        &
     &                LROOTS,NROOTS,IROOT(MAXRTS),                      &
     &                NOROT(MXCORB),IMOORD(MXCORB),                     &
     &                IORTO,ICI0,KDEL,ICHECK,NTIT,                      &
     &                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,  &
     &                ITRLVL,ITRFIN,JCHSYM,JCHORB,                      &
     &                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO, NMCAVER
!
      LOGICAL         FLAG,        DOSCF,DOMP2,DOCINO,DOCI,DOMC,DORSP,  &
     &                FCVORB,LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP, &
     &                JOLSEN,ABAIPH,INERSI,INERSF,DODFT,DONEVPT,HSROHF, &
     &                BOYORB,PIPORB,ADDMP2,DOFCI,DOCISD,DOLUCITA,DOMEP
!     variables for srDFT /hjaaj                                         
      LOGICAL         DOCISRDFT,DOHFSRDFT,DOMCSRDFT,ADDSRI,SRHYBR        
      COMMON /LOGINP/ FLAG(NFLAG), DOSCF,DOMP2,DOCINO,DOCI,DOMC,DORSP,  &
     &                FCVORB,LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP, &
     &                JOLSEN,ABAIPH,INERSI,INERSF,DODFT,DONEVPT,HSROHF, &
     &                BOYORB,PIPORB,ADDMP2,DOFCI,DOCISD,DOLUCITA,DOMEP, &
     &                DOCISRDFT,DOHFSRDFT,DOMCSRDFT,ADDSRI,SRHYBR
      LOGICAL         SUPSYM, DORHF
      EQUIVALENCE (SUPSYM,FLAG(17)), (DOSCF,DORHF)
!
      REAL*8          SPIN, POTNUC, EPSOL,EPSTAT,EPPN,RSOL,             &
     &                THRGRD, THRPWF, THRCI, THRMC, THRCGR,             &
     &                EFIELD, CMAXMO, THROVL,                           &
     &                THRSSY, DEFLVL, WEIGHT_MCAVER
!     variables for srDFT /hjaaj
      REAL*8          THRCIDFT
      COMMON /RELINP/ SPIN, POTNUC, EPSOL,EPSTAT,EPPN,RSOL(3),          &
     &                THRGRD, THRPWF, THRCI, THRMC, THRCGR,             &
     &                EFIELD(MXFELT), CMAXMO, THROVL,                   &
     &                THRSSY, DEFLVL, WEIGHT_MCAVER(MAXRTS),            &
     &                THRCIDFT
!
      CHARACTER*60 TITLE(6)
      CHARACTER*72 TITMOL(2)
      CHARACTER*4  CENT,   TYPE
      CHARACTER*8  LFIELD
      COMMON /CHRINP/ TITLE, TITMOL,                                    &
     &                CENT(MXCORB), TYPE(MXCORB), LFIELD(MXFELT)
!-- end of infinp.h --
