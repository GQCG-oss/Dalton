!
!  FILE : infrsp.h
!
!  control variables for the dalton/rsp/ module
!
      INTEGER MFREQ, MCFREQ
      PARAMETER ( MFREQ = 30 , MCFREQ = 30 )

      REAL*8          THCRSP, FREQ,   CFREQ,  ORBSFT
      INTEGER         IREFSY, IPRRSP, MAXIT,  MAXITO,                    &     
     &                NFREQ,  NCFREQ, NCREF,  ISTOCK, MAXOCK,            &
     &                MAXRM,  LPVMAT, NACTT,  NACT,   IACT, JACT
      LOGICAL         OPTORB, INPTES, ANTTES, AVDIA,  TDHF,   RSPCI
      LOGICAL         NOITRA, ORBSPC, ABOCHK, TRPLET, OLSEN,  PHPRES
      LOGICAL         E3TEST, TRPFLG, A2TEST, X2TEST, DIROIT, RSPSUP
      LOGICAL         SOPPA,  HIRPA,  SOPW4,  CCPPA, TDA, CISRPA         !SPAS : 06/11-2009 AOSOP for AO-SOPPA is in soppinf.h
      LOGICAL         DFT_SO, RSPECD, RSPOCD, SOPRPA, TRDQF
      COMMON /INFRSP/ THCRSP, FREQ(MFREQ),    CFREQ(MCFREQ),  ORBSFT,    & ! real*8 variables
     &                IREFSY, IPRRSP, MAXIT,  MAXITO,                    & ! integer variables
     &                NFREQ,  NCFREQ, NCREF,  ISTOCK, MAXOCK,            &
     &                MAXRM,  LPVMAT, NACTT,  NACT(8),IACT(8),JACT(8),   &
     &                OPTORB, INPTES, ANTTES, AVDIA,  TDHF,   RSPCI,     & ! logical variables
     &                NOITRA, ORBSPC, ABOCHK, TRPLET, OLSEN,  PHPRES,    &
     &                E3TEST, TRPFLG, A2TEST, X2TEST, DIROIT, RSPSUP,    &
     &                SOPPA,  HIRPA,  SOPW4,  CCPPA, TDA, CISRPA,        &
     &                DFT_SO, RSPECD, RSPOCD, SOPRPA, TRDQF
! -- end of infrsp.h --
