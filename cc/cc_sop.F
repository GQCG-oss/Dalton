C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
c /* deck CC_SOP */
*=====================================================================*
       SUBROUTINE CC_SOP(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of (frequency-dependent)
*             second-order properties for the Coupled Cluster models
*
*                        CCS, CC2, CCSD, CC3
*
*             and partially for SCF and CIS
*
*     Written by Christof Haettig summer/fall 1999.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclrinf.h>
#include <ccr1rsp.h>
#include <ccrc1rsp.h>
#include <ccroper.h>
#include <cclists.h>
#include <second.h>

* local parameters:
      CHARACTER*(17) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_SOP> ')

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      LOGICAL LADD
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0, TIM1, TIMF, TIMXE, TIMXI, TIMR
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0, TIM1, TIMF, TIMXE, TIMXI, TIMR
#endif

      INTEGER NBSOP, MXTRAN, MXVEC
      INTEGER MXFTRAN, MXXITRAN, MXRTRAN
      INTEGER MXFDOTS, MXXIDOTS, MXRDOTS
      INTEGER MXXETRAN, MXXEDOTS
      INTEGER NFTRAN, NXITRAN, NRTRAN, NXETRAN
      INTEGER KRESULT, KFTRAN, KFDOTS, KFCONS, KXITRAN, KXIDOTS, KXICONS
      INTEGER KXETRAN, KXDOTS, KXCONS, KEDOTS, KECONS
      INTEGER KEND0, LEND0, IOPT, KRCONS, KRDOTS, KRTRAN, IORDER
      INTEGER KEND1, LWRK1, LRESULT

* external functions:
      INTEGER IR1TAMP
      INTEGER IL1ZETA

*---------------------------------------------------------------------*
* print header for second-order property section:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<  OUTPUT FROM COUPLED CLUST',
     &                               'ER LINEAR RESPONSE   >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<     CALCULATION OF SECOND ',
     &                               'ORDER PROPERTIES     >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CC_SOP called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_SOP Workspace:',LWORK
  
      TIM0  = SECOND()

*---------------------------------------------------------------------*
* allocate & initialize work space for property contributions:
*---------------------------------------------------------------------*
      NBSOP = NLROP * NBLRFR

      MXTRAN  = MAX(NLRTLBL,NLRTHFLBL)
      MXVEC   = MAX(NLRTLBL,NLRTHFLBL)

      MXFTRAN  = MXDIM_FTRAN * MXTRAN
      MXXITRAN = 1 * MXTRAN
      MXRTRAN  = 1 * MXTRAN
      MXXETRAN = MXDIM_XEVEC * MXTRAN

      MXFDOTS  = MXVEC * MXTRAN
      MXXIDOTS = MXVEC * MXTRAN
      MXRDOTS  = MXVEC * MXTRAN
      MXXEDOTS = MXVEC * MXTRAN

      LRESULT = 2*NBSOP

      KRESULT = 1
      KFTRAN  = KRESULT + LRESULT
      KFDOTS  = KFTRAN  + MXFTRAN
      KFCONS  = KFDOTS  + MXFDOTS
      KXITRAN = KFCONS  + MXFDOTS
      KXIDOTS = KXITRAN + MXXITRAN
      KXICONS = KXIDOTS + MXXIDOTS
      KRTRAN  = KXICONS + MXXIDOTS
      KRDOTS  = KRTRAN  + MXRTRAN
      KRCONS  = KRDOTS  + MXRDOTS
      KXETRAN = KRCONS  + MXRDOTS
      KXDOTS  = KXETRAN + MXXETRAN
      KXCONS  = KXDOTS  + MXXEDOTS
      KEDOTS  = KXCONS  + MXXEDOTS
      KECONS  = KEDOTS  + MXXEDOTS
      KEND0   = KECONS  + MXXEDOTS
      LEND0   = LWORK   - KEND0

      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_SOP. (1)')
      END IF

      CALL DZERO(WORK(KRESULT),2*NBSOP)

*---------------------------------------------------------------------*
* set up lists for F transformations, ETA{O} and Xi{O} vectors:
*---------------------------------------------------------------------*
      LADD = .FALSE.

      CALL CCLR_SETUP(MXTRAN, MXVEC,
     &           WORK(KFTRAN), WORK(KFDOTS), WORK(KFCONS), NFTRAN,
     &           WORK(KXITRAN),WORK(KXIDOTS),WORK(KXICONS),NXITRAN,
     &           WORK(KRTRAN), WORK(KRDOTS), WORK(KRCONS), NRTRAN,
     &           WORK(KXETRAN),WORK(KXDOTS), WORK(KEDOTS),
     &                         WORK(KXCONS), WORK(KECONS),NXETRAN,
     &           WORK(KRESULT),NBSOP,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IOPT = 5
      CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'L0 ','R1 ',IOPT,'R1 ',
     &                WORK(KFDOTS),WORK(KFCONS),MXVEC,
     &                WORK(KEND0), LEND0)

      TIMF = SECOND() - TIM1

      IF (NFTRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NFTRAN,' F matrix transformations:',TIMF
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate ETA{O} x R1 and Xksi{O} x L1 vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KXCONS),MXVEC*NXETRAN)
      CALL DZERO(WORK(KECONS),MXVEC*NXETRAN)

      IOPT   = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KXETRAN), NXETRAN, IOPT, IORDER, 'L0 ',
     &               'L1 ',WORK(KXDOTS),WORK(KXCONS),
     &               'R1 ',WORK(KEDOTS),WORK(KECONS),
     &               .FALSE.,MXVEC, WORK(KEND0), LEND0 )

      TIMXE = SECOND() - TIM1
      IF (NXETRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NXETRAN,' O1/X1 vector calculation:',TIMXE
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate relaxation X intermediate contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL CCXIMCON(WORK(KXITRAN),NXITRAN,'R1 ',WORK(KXIDOTS),
     &              WORK(KXICONS),MXVEC,WORK(KEND0),LEND0)

      TIMXI = SECOND() - TIM1
      IF (NXITRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NXITRAN,' relaxation contributions:',TIMXI
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate second-order reorthogonalization/relaxation contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL CCREO2CON(WORK(KRTRAN),NRTRAN,WORK(KRDOTS),
     &               WORK(KRCONS),MXVEC,WORK(KEND0),LEND0)

      TIMR = SECOND() - TIM1
      IF (NRTRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NRTRAN ,' 2. order relax. contrbs.:',TIMR
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* collect contributions and sum them up to the final results:
*---------------------------------------------------------------------*
      LADD = .TRUE.

      CALL CCLR_SETUP(MXTRAN, MXVEC,
     &           WORK(KFTRAN), WORK(KFDOTS), WORK(KFCONS), NFTRAN,
     &           WORK(KXITRAN),WORK(KXIDOTS),WORK(KXICONS),NXITRAN,
     &           WORK(KRTRAN), WORK(KRDOTS), WORK(KRCONS), NRTRAN,
     &           WORK(KXETRAN),WORK(KXDOTS), WORK(KEDOTS),
     &                         WORK(KXCONS), WORK(KECONS),NXETRAN,
     &           WORK(KRESULT),NBSOP,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBSOP,' linear response function:', SECOND() - TIM0

*---------------------------------------------------------------------*
* print frequency-dependent polarizabilities:
*---------------------------------------------------------------------*

      KEND1 = KRESULT + LRESULT
      LWRK1 = LWORK   - KEND1 + 1

      CALL  CCSOPPRT(WORK(KRESULT),WORK(KEND1),LWRK1)

      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* print dipole gradients and APT (Cioslowski) population analysis
*---------------------------------------------------------------------*
      IF (IDIGRA .GT. 0) THEN
         IF (CCS) THEN
            CALL AROUND('CCS Dipole Gradients')
         ELSE IF (CC2) THEN
            CALL AROUND('CC2 Dipole Gradients')
         ELSE IF (CC3) THEN
            CALL AROUND('CC3 Dipole Gradients')
         ELSE IF (CCSD) THEN
            CALL AROUND('CCSD Dipole Gradients')
         ELSE
            CALL QUIT('CC_SOP called for unknown Coupled Cluster.')
         END IF
         LEND0 = LWORK - KFTRAN + 1
         CALL CCDIPGRAPR(WORK(KRESULT),WORK(KFTRAN),LEND0)
      ENDIF

C     IF (NQRDISP.EQ.0) THEN
        RETURN
C     END IF
  
      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_SOP                               *
*=====================================================================*
c /* deck ccsopprt */
*=====================================================================*
       SUBROUTINE CCSOPPRT(SOPROP,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: print (frequency-dependent) second-order properties
*
*     Written by Christof Haettig in summer 1999.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
       LOGICAL DIPGRA
#include <cclrinf.h>
#include <ccroper.h>


      CHARACTER*5  BLANKS
      CHARACTER*8  LABA, LABB
      CHARACTER*10 MODEL, RLXA, RLXB
      CHARACTER*80 STRING
      INTEGER LWORK
      INTEGER ISYMA, ISYMB, ISAMA, ISAMB, ISAPROP, IRELAX, M
      INTEGER IFREQ, IOPER, ICMP, ITURN, MXTURN, ITST, ITST2, IPAR
      INTEGER IFRZER, IBLRFR


#if defined (SYS_CRAY)
      REAL SOPROP(NBLRFR,NLROP,2), WORK(LWORK)
      REAL GLEN(3,MBLRFR), GVEL(3,MBLRFR), DELTA(3,MBLRFR)
      REAL GTST(3,MBLRFR)
      REAL SIGN, ANI, ALPHA, PROP, ERROR, ALPHASI, ANISI
      REAL ALFSQ, BETSQ, APAR, APEN, XKAPPA
      REAL TENS(3,3), PVALUES(3), PAXIS(3,3)
      REAL ZERO, HALF, ONE, DPOLAUTSI
#else
      DOUBLE PRECISION SOPROP(NBLRFR,NLROP,2), WORK(LWORK)
      DOUBLE PRECISION GLEN(3,MBLRFR), GVEL(3,MBLRFR), DELTA(3,MBLRFR)
      DOUBLE PRECISION GTST(3,MBLRFR)
      DOUBLE PRECISION SIGN, ANI, ALPHA, PROP, ERROR, ALPHASI, ANISI
      DOUBLE PRECISION ALFSQ, BETSQ, APAR, APEN, XKAPPA
      DOUBLE PRECISION TENS(3,3), PVALUES(3), PAXIS(3,3)
      DOUBLE PRECISION ZERO, HALF, ONE, DPOLAUTSI
#endif
      PARAMETER (ZERO = 0.0d0, HALF = 0.5d0, ONE = 1.0d0)
      PARAMETER (DPOLAUTSI =  1.648778D-41)

*---------------------------------------------------------------------*
* print header for second-order properties: 
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING = ' RESULTS FOR THE SECOND-ORDER PROPERTIES '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:41)//BLANKS ) 
         MODEL = 'CCS'
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:41)//BLANKS )
         MODEL = 'CC2'
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:41)//BLANKS )
         MODEL = 'CC3'
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:41)//BLANKS )
         MODEL = 'CCSD'
      ELSE
         CALL QUIT('CCSOPPRT called for an unknown Coupled '//
     &             'Cluster model.')
      END IF

      IF (IPRSOP.GT.5) THEN
       WRITE(LUPRI,'(/1X,2(1X,A," operator",17X),4X,A,6X,A,/,100("-"))')
     &      'A','B','property ','(asy./sym. Resp.)'
      ELSE
       WRITE(LUPRI,'(/1X,2(1X,A," operator",17X),4X,A,/,72("-"))')
     &      'A','B','property '
      END IF


      DO IOPER = 1, NLROP
         ISYMA = ISYOPR(IALROP(IOPER))
         ISYMB = ISYOPR(IBLROP(IOPER))

         ISAMA = ISYMAT(IALROP(IOPER))
         ISAMB = ISYMAT(IBLROP(IOPER))

         IF (       LALORX(IOPER) ) RLXA = ' (relax.) '
         IF ( .NOT. LALORX(IOPER) ) RLXA = ' (unrel.) '
         IF (       LBLORX(IOPER) ) RLXB = ' (relax.) '
         IF ( .NOT. LBLORX(IOPER) ) RLXB = ' (unrel.) '

         ISAPROP = ISAMA * ISAMB 
         SIGN    = DBLE(ISAPROP)
         IF (ISAPROP.EQ.0) SIGN = +ONE

         IFREQ = 1
         IF (ISYMA.EQ.ISYMB) THEN

           PROP  = -HALF*(SOPROP(IFREQ,IOPER,1) +
     &                     SIGN * SOPROP(IFREQ,IOPER,2))
           ERROR = -HALF*(SOPROP(IFREQ,IOPER,1) -
     &                     SIGN * SOPROP(IFREQ,IOPER,2))

           IF (IPRSOP.GT.5 .OR. ISAPROP.EQ.0) THEN
            IF (ISAPROP.NE.0) THEN
             WRITE(LUPRI,
     &              '(/1X,2(A8,A10,F7.4,3X),G18.10," (",G18.10,")")')
     &          LBLOPR(IALROP(IOPER)),RLXA,ALRFR(IFREQ),
     &          LBLOPR(IBLROP(IOPER)),RLXB,BLRFR(IFREQ), PROP, ERROR
            ELSE
              WRITE(LUPRI,'(/1X,2A)')
     &          'Cannot determine if real or imaginary property...   ',
     &          'print: symmetric (antisymmetric) parts in +/- w'
              WRITE(LUPRI,
     &             '(1X,2(A8,A10,F7.4,3X),G18.10," (",G18.10,")")')
     &          LBLOPR(IALROP(IOPER)),RLXA,ALRFR(IFREQ),
     &          LBLOPR(IBLROP(IOPER)),RLXB,BLRFR(IFREQ), PROP, ERROR
            END IF
           ELSE
            WRITE(LUPRI,'(/1X,2(A8,A10,F7.4,3X),G16.8)')
     &        LBLOPR(IALROP(IOPER)),RLXA,ALRFR(IFREQ),
     &        LBLOPR(IBLROP(IOPER)),RLXB,BLRFR(IFREQ), PROP
           ENDIF
         ELSE
           IF (IPRSOP.GT.5 .OR. ISAPROP.EQ.0) THEN
            WRITE(LUPRI,
     &             '(/1X,2(A8,A10,A7,3X),7X,A,8X," (",9X,A,10X,")")')
     &        LBLOPR(IALROP(IOPER)),RLXA,'   -.-  ',
     &        LBLOPR(IBLROP(IOPER)),RLXB,'   -.-  ',
     &        '---',
     &        '---'
           ELSE IF (IPRSOP.GT.0) THEN
            WRITE(LUPRI,'(/1X,2(A8,A10,A7,3X),6X,A,7X)')
     &        LBLOPR(IALROP(IOPER)),RLXA,'   -.-  ',
     &        LBLOPR(IBLROP(IOPER)),RLXB,'   -.-  ',
     &        '---'
           END IF
         END IF

         DO IFREQ = 2, NBLRFR
           IF (ISYMA.EQ.ISYMB) THEN

             PROP  = -HALF*(SOPROP(IFREQ,IOPER,1) +
     &                      SIGN * SOPROP(IFREQ,IOPER,2))
             ERROR = -HALF*(SOPROP(IFREQ,IOPER,1) -
     &                      SIGN * SOPROP(IFREQ,IOPER,2))

             IF (IPRSOP.GT.5 .OR. ISAPROP.EQ.0) THEN
               WRITE(LUPRI,'(1X,2(18X,F7.4,3X),G18.10," (",G18.10,")")')
     &           ALRFR(IFREQ), BLRFR(IFREQ), PROP, ERROR
             ELSE
               WRITE(LUPRI,'(1X,2(18X,F7.4,3X),G16.8)')
     &           ALRFR(IFREQ), BLRFR(IFREQ), PROP
             END IF
           END IF
         END DO

      END DO

      IF (IPRSOP.GT.5) THEN
       WRITE(LUPRI,'(100("-"),/,/)')
      ELSE
       WRITE(LUPRI,'(72("-"),/,/)')
      END IF



      IF (ALPHA_ISO .OR. ALPHA_ANI) THEN

       WRITE(LUPRI,'(/,/2X,A,36x,A,3X,A,/,80("-"))')
     &       "average  ", "value (au)", "value (C^2 m^2 J^-1)"

       IF (ALPHA_ISO) THEN

          DO IRELAX = 1, 2
            M = OFFALPHA(IRELAX)
            IF (IRELAX.EQ.1) RLXA = ' (relax.) '
            IF (IRELAX.EQ.2) RLXA = ' (unrel.) '
          IF (M.GE.0) THEN

          DO IFREQ = 1, NBLRFR

            IF (CSYM(1:6).EQ.'ATOMIC') THEN
             ALPHA=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
            ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
             ALPHA=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
            ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
             ALPHA=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
     &             -ONE *(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
             ALPHA=ALPHA / 3.0D0
            ELSE IF (CSYM(1:6).EQ.'XYDGEN') THEN
             ALPHA=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
     &             -ONE *(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
             ALPHA=ALPHA / 3.0D0
            ELSE
             ALPHA=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
     &             -HALF*(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
     &             -HALF*(SOPROP(IFREQ,M+5,1)+SIGN*SOPROP(IFREQ,M+5,2))
             ALPHA=ALPHA / 3.0D0
            END IF

            ALPHASI = ALPHA * DPOLAUTSI
            
            IF (IFREQ.EQ.1) THEN
              WRITE(LUPRI,'(/1X,A9,A10,4X,2(F7.4,3X),G16.8,3X,G16.8)')
     &         'alpha_iso',RLXA,ALRFR(IFREQ),BLRFR(IFREQ),ALPHA,ALPHASI
            ELSE
              WRITE(LUPRI,'(/1X,9X,10X,4X,2(F7.4,3X),G16.8,3X,G16.8)')
     &                          ALRFR(IFREQ),BLRFR(IFREQ),ALPHA,ALPHASI
            END IF

          END DO
          END IF
          END DO
       END IF

       MXTURN = 2
       IF (IPRSOP.LT.5) MXTURN = 1

       DO ITURN = 1, MXTURN
        IF (ALPHA_ANI) THEN

        DO IRELAX = 1, 2
          M = OFFALPHA(IRELAX)
          IF (IRELAX.EQ.1) RLXA = ' (relax.) '
          IF (IRELAX.EQ.2) RLXA = ' (unrel.) '
        IF (M.GE.0) THEN

        DO IFREQ = 1, NBLRFR

         IF     (CSYM(1:6).EQ.'ATOMIC') THEN
          TENS(3,3)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(1,1)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(2,2)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(1,3)=ZERO
          TENS(2,3)=ZERO
          TENS(1,2)=ZERO
         ELSE IF(CSYM(1:6).EQ.'SPHTOP') THEN
          TENS(3,3)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(1,1)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(2,2)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(1,3)=ZERO
          TENS(2,3)=ZERO
          TENS(1,2)=ZERO
         ELSE IF(CSYM(1:6).EQ.'LINEAR') THEN
          TENS(3,3)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(1,1)=-HALF*(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
          TENS(2,2)=-HALF*(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
          TENS(1,3)=ZERO
          TENS(2,3)=ZERO
          TENS(1,2)=ZERO
         ELSE IF (CSYM(1:6).EQ.'XYDGEN') THEN
          TENS(3,3)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(1,1)=-HALF*(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
          TENS(2,2)=-HALF*(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
          TENS(1,3)=-HALF*(SOPROP(IFREQ,M+3,1)+SIGN*SOPROP(IFREQ,M+3,2))
          TENS(2,3)=-HALF*(SOPROP(IFREQ,M+3,1)+SIGN*SOPROP(IFREQ,M+3,2))
          TENS(1,2)=-HALF*(SOPROP(IFREQ,M+4,1)+SIGN*SOPROP(IFREQ,M+4,2))
         ELSE
          TENS(3,3)=-HALF*(SOPROP(IFREQ,M+1,1)+SIGN*SOPROP(IFREQ,M+1,2))
          TENS(1,1)=-HALF*(SOPROP(IFREQ,M+2,1)+SIGN*SOPROP(IFREQ,M+2,2))
          TENS(2,2)=-HALF*(SOPROP(IFREQ,M+5,1)+SIGN*SOPROP(IFREQ,M+5,2))
          TENS(1,3)=-HALF*(SOPROP(IFREQ,M+3,1)+SIGN*SOPROP(IFREQ,M+3,2))
          TENS(2,3)=-HALF*(SOPROP(IFREQ,M+6,1)+SIGN*SOPROP(IFREQ,M+6,2))
          TENS(1,2)=-HALF*(SOPROP(IFREQ,M+4,1)+SIGN*SOPROP(IFREQ,M+4,2))
         END IF

         TENS(2,1) = TENS(1,2)
         TENS(3,1) = TENS(1,3)
         TENS(3,2) = TENS(2,3)

         CALL TNSRAN(TENS,PVALUES,PAXIS,ALFSQ,BETSQ,ITST,ITST2,
     &               APAR,APEN,XKAPPA,IPAR)

         ANI   = DSQRT(BETSQ)
         ANISI = ANI * DPOLAUTSI

         IF (ITURN.EQ.1) THEN
              IF (IFREQ.EQ.1) THEN
                WRITE(LUPRI,'(/1X,A9,A10,4X,2(F7.4,3X),G16.8,3X,G16.8)')
     &            'alpha_ani',RLXA,ALRFR(IFREQ),BLRFR(IFREQ), ANI, ANISI
              ELSE
                WRITE(LUPRI,'(/1X,9X,10X,4X,2(F7.4,3X),G16.8,3X,G16.8)')
     &                             ALRFR(IFREQ),BLRFR(IFREQ), ANI, ANISI
              END IF

         ELSE
              WRITE(LUPRI,'(//1X,A6,1X,A6,A,F10.6)') 
     &          RLXA(3:8), MODEL(1:6),
     &          'polarizability for frequency: ',DABS(ALRFR(IFREQ))
              WRITE(LUPRI,'(1X,53("."),/)')
              IF (IPRSOP.GE.6) CALL OUTPUT(TENS,1,3,1,3,3,3,1,LUPRI)
              WRITE(LUPRI,'(2(/,1X,A,20X,F14.6))')
     &                  'alpha^2 Invariant:        ',ALFSQ,
     &                  'beta^2 Invariant:         ',BETSQ,
     &                  'Isotropic Polarizability: ',DSQRT(ALFSQ),
     &                  'Polarizability anisotropy:',DSQRT(BETSQ)
              IF      (ITST.EQ.0 .AND. ITST2.EQ.3) THEN
                 WRITE(LUPRI,'(/1X,A)')
     &            'Polarizability has spherical symmetry.'
              ELSE IF (ITST.EQ.0 .AND. ITST2.EQ.1) THEN
                 WRITE(LUPRI,'(/1X,A)')
     &            'Polarizability has cylinder symmetry.'
                 WRITE(LUPRI,'(2(1X,A,F10.6,A,3X,E15.6,A,/))')
     &            'Parallel component:      ',APAR,' a.u.',
     &                  APAR*DPOLAUTSI,' C^2 m^2 J^-1',
     &            'Perpendicular component: ',APEN,' a.u.',
     &                  APEN*DPOLAUTSI,' C^2 m^2 J^-1'
                 WRITE(LUPRI,'(/,1X,A42,F12.6)')
     &            'Dimensionless polarizability anisotropy:  ',XKAPPA
              ELSE IF (ITST.EQ.0 .AND. ITST2.EQ.0) THEN
                 WRITE(LUPRI,'(/1X,A)')
     &            'Polarizability is diagonal.'
              ELSE
                 WRITE(LUPRI,'(/,1X,A,/)')
     &             'Principal values of diagonalized Polarizability:'
                 WRITE(LUPRI,'(1X,A)')
     &                '        a.u.          C^2 m^2 J^-1'
                 WRITE(LUPRI,'(3(1X,A,F10.6,3X,E15.6,/))')
     &              '1     ',PVALUES(1),PVALUES(1)*DPOLAUTSI,
     &              '2     ',PVALUES(2),PVALUES(2)*DPOLAUTSI,
     &              '3     ',PVALUES(3),PVALUES(3)*DPOLAUTSI
                 WRITE(LUPRI,'(/,1X,A,/)')
     &             'Principal axis of diagonalized Polarizability:'
                 CALL OUTPUT(PAXIS,1,3,1,3,3,3,1,LUPRI) 
              END IF
              WRITE(LUPRI,'(80("-"))')
          END IF

        END DO
        END IF
        END DO
        END IF

        IF (ITURN.EQ.1)  WRITE(LUPRI,'(80("-"))')

        END DO
         
      END IF

C     Optical rotation:
C     Set up diagonal elements of G' tensor and call print routine.
C     -------------------------------------------------------------

      IF (ORLEN .OR. ORVEL .OR. ORMVE) THEN

         CALL DZERO(GLEN,3*NBLRFR)
         CALL DZERO(GVEL,3*NBLRFR)

         IF (ORMVE) THEN
            CALL DZERO(GTST,3*NBLRFR)
            IFRZER = 0
            IBLRFR = 1
            DO WHILE ((IBLRFR.LE.NBLRFR) .AND. (IFRZER.EQ.0))
               IF (ABS(BLRFR(IBLRFR)) .LT. 1.0D-10) IFRZER = IBLRFR
               IBLRFR = IBLRFR + 1
            END DO
            IF (IFRZER .EQ. 0) THEN
               CALL QUIT('Zero freq. not found for .OR MVE')
            END IF
         END IF

         DO IOPER = 1,NLROP

            LABA  = LBLOPR(IALROP(IOPER))
            ISYMA = ISYOPR(IALROP(IOPER))
            ISAMA = ISYMAT(IALROP(IOPER))

            LABB  = LBLOPR(IBLROP(IOPER))
            ISYMB = ISYOPR(IBLROP(IOPER))
            ISAMB = ISYMAT(IBLROP(IOPER))

            IF (LABA(2:4) .NE. 'DIP') GO TO 999
            IF (LABB(2:7) .NE. 'ANGMOM') GO TO 999
            IF (ISYMA .NE. ISYMB) GO TO 999

            ISAPROP = ISAMA * ISAMB
            SIGN    = DBLE(ISAPROP)
            IF (ISAPROP.EQ.0) SIGN = +ONE

            IF (ORLEN) THEN
               IF ((LABA(1:7).EQ.'XDIPLEN') .AND.
     &             (LABB(1:7).EQ.'XANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                            SIGN * SOPROP(IFREQ,IOPER,2))
                     GLEN(1,IFREQ) = -HALF*PROP
                  ENDDO
               ELSE IF ((LABA(1:7).EQ.'YDIPLEN') .AND.
     &                  (LABB(1:7).EQ.'YANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                            SIGN * SOPROP(IFREQ,IOPER,2))
                     GLEN(2,IFREQ) = -HALF*PROP
                  ENDDO
               ELSE IF ((LABA(1:7).EQ.'ZDIPLEN') .AND.
     &                  (LABB(1:7).EQ.'ZANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                            SIGN * SOPROP(IFREQ,IOPER,2))
                     GLEN(3,IFREQ) = -HALF*PROP
                  ENDDO
               ENDIF
            ENDIF

            IF (ORVEL) THEN
               IF ((LABA(1:7).EQ.'XDIPVEL') .AND.
     &             (LABB(1:7).EQ.'XANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        GVEL(1,IFREQ) = -HALF*PROP/BLRFR(IFREQ)
                     ELSE
                        GVEL(1,IFREQ) = 1.0D15
                     ENDIF
                  ENDDO
               ELSE IF ((LABA(1:7).EQ.'YDIPVEL') .AND.
     &                  (LABB(1:7).EQ.'YANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        GVEL(2,IFREQ) = -HALF*PROP/BLRFR(IFREQ)
                     ELSE
                        GVEL(2,IFREQ) = 1.0D15
                     ENDIF
                  ENDDO
               ELSE IF ((LABA(1:7).EQ.'ZDIPVEL') .AND.
     &                  (LABB(1:7).EQ.'ZANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        GVEL(3,IFREQ) = -HALF*PROP/BLRFR(IFREQ)
                     ELSE
                        GVEL(3,IFREQ) = 1.0D15
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (ORMVE) THEN
               IF ((LABA(1:7).EQ.'XDIPVEL') .AND.
     &             (LABB(1:7).EQ.'XANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) -
     &                               SOPROP(IFRZER,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2) -
     &                               SIGN * SOPROP(IFRZER,IOPER,2))
                        GTST(1,IFREQ) = -HALF*PROP/BLRFR(IFREQ)
                     ELSE
                        GTST(1,IFREQ) = 1.0D15
                     ENDIF
                  ENDDO
               ELSE IF ((LABA(1:7).EQ.'YDIPVEL') .AND.
     &                  (LABB(1:7).EQ.'YANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) -
     &                               SOPROP(IFRZER,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2) -
     &                               SIGN * SOPROP(IFRZER,IOPER,2))
                        GTST(2,IFREQ) = -HALF*PROP/BLRFR(IFREQ)
                     ELSE
                        GTST(2,IFREQ) = 1.0D15
                     ENDIF
                  ENDDO
               ELSE IF ((LABA(1:7).EQ.'ZDIPVEL') .AND.
     &                  (LABB(1:7).EQ.'ZANGMOM')) THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) -
     &                               SOPROP(IFRZER,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2) -
     &                               SIGN * SOPROP(IFRZER,IOPER,2))
                        GTST(3,IFREQ) = -HALF*PROP/BLRFR(IFREQ)
                     ELSE
                        GTST(3,IFREQ) = 1.0D15
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

  999       CONTINUE

         ENDDO
         CALL CC_OROUT(GLEN,GVEL,GTST,WORK,LWORK,IPRSOP,MODEL)

      ENDIF

C     Optical activity: Origin analysis.
C     ----------------------------------

      IF (ORGANL) THEN

         CALL DZERO(DELTA,3*NBLRFR)

         DO IOPER = 1,NLROP

            LABA  = LBLOPR(IALROP(IOPER))
            ISYMA = ISYOPR(IALROP(IOPER))
            ISAMA = ISYMAT(IALROP(IOPER))

            LABB  = LBLOPR(IBLROP(IOPER))
            ISYMB = ISYOPR(IBLROP(IOPER))
            ISAMB = ISYMAT(IBLROP(IOPER))

            IF (ISYMA .NE. ISYMB) GO TO 1000

            ISAPROP = ISAMA * ISAMB
            SIGN    = DBLE(ISAPROP)
            IF (ISAPROP.EQ.0) SIGN = +ONE

            IF (LABA(1:7) .EQ. 'XDIPLEN') THEN
               IF (LABB(1:7) .EQ. 'YDIPVEL') THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        DELTA(3,IFREQ) = DELTA(3,IFREQ)
     &                                 + PROP/(6.0D0*BLRFR(IFREQ))
                     ELSE
                        DELTA(3,IFREQ) = DELTA(3,IFREQ) + 1.0D15
                     ENDIF
                  ENDDO
               ELSE IF (LABB(1:7) .EQ. 'ZDIPVEL') THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        DELTA(2,IFREQ) = DELTA(2,IFREQ)
     &                                 - PROP/(6.0D0*BLRFR(IFREQ))
                     ELSE
                        DELTA(2,IFREQ) = DELTA(2,IFREQ) - 1.0D15
                     ENDIF
                  ENDDO
               ENDIF
            ELSE IF (LABA(1:7) .EQ. 'YDIPLEN') THEN
               IF (LABB(1:7) .EQ. 'XDIPVEL') THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        DELTA(3,IFREQ) = DELTA(3,IFREQ)
     &                                 - PROP/(6.0D0*BLRFR(IFREQ))
                     ELSE
                        DELTA(3,IFREQ) = DELTA(3,IFREQ) - 1.0D15
                     ENDIF
                  ENDDO
               ELSE IF (LABB(1:7) .EQ. 'ZDIPVEL') THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        DELTA(1,IFREQ) = DELTA(1,IFREQ)
     &                                 + PROP/(6.0D0*BLRFR(IFREQ))
                     ELSE
                        DELTA(1,IFREQ) = DELTA(1,IFREQ) + 1.0D15
                     ENDIF
                  ENDDO
               ENDIF
            ELSE IF (LABA(1:7) .EQ. 'ZDIPLEN') THEN
               IF (LABB(1:7) .EQ. 'XDIPVEL') THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        DELTA(2,IFREQ) = DELTA(2,IFREQ)
     &                                 + PROP/(6.0D0*BLRFR(IFREQ))
                     ELSE
                        DELTA(2,IFREQ) = DELTA(2,IFREQ) + 1.0D15
                     ENDIF
                  ENDDO
               ELSE IF (LABB(1:7) .EQ. 'YDIPVEL') THEN
                  DO IFREQ = 1,NBLRFR
                     IF (DABS(BLRFR(IFREQ)).GT.ZERO) THEN
                        PROP = HALF*(SOPROP(IFREQ,IOPER,1) +
     &                               SIGN * SOPROP(IFREQ,IOPER,2))
                        DELTA(1,IFREQ) = DELTA(1,IFREQ)
     &                                 - PROP/(6.0D0*BLRFR(IFREQ))
                     ELSE
                        DELTA(1,IFREQ) = DELTA(1,IFREQ) - 1.0D15
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

C           Escape point.
C           -------------

 1000       CONTINUE

         ENDDO

C        Call print routine.
C        -------------------

         CALL CC_ORANL(DELTA,IPRSOP,MODEL)

      ENDIF

C     Optical activity: Calculate and print specific rotation
C                       in length gauge at specified origins.
C     --------------------------------------------------------

      IF (NORGIN .GT. 0) THEN

         IF ((.NOT.ORLEN) .OR. (.NOT.ORGANL)) THEN
            WRITE(LUPRI,'(//,1X,A,I4,A)')
     &      'NORGIN =',NORGIN,' requires ORLEN and ORGANL!'
            CALL QUIT('Inconsistence in CC_SOP')
         ENDIF

         CALL CC_ORGIN(GLEN,DELTA,IPRINT,MODEL)

      ENDIF

      RETURN
      END
C  /* Deck cc_orgin */
      SUBROUTINE CC_ORGIN(GLEN,DELTA,IPRINT,MODEL)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Origin-dependence of length gauge OR.
C
#include <implicit.h>
      DIMENSION GLEN(3,*), DELTA(3,*)
      CHARACTER*10 MODEL
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <codata.h>
#include <nuclei.h>
#include <symmet.h>
#include <cclrinf.h>
#include <cbiher.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_ORGIN')

      PARAMETER (TINY = 1.0D-7)
      PARAMETER (ZERO = 0.0D0, THREE = 3.0D0)
      PARAMETER (HALF = 0.5D0)
      PARAMETER (HUNDRD = 1.0D2)
      PARAMETER (EHTGHZ = XTHZ*1.0D-9)
C FACTOR = 28800 * pi**2 * bohr**4 * Avogadro, with bohr in cm
C        = 288.0D-30 * pi**2 * bohr**4 * Avogadro, with bohr in Angstrom
      PARAMETER
     & (FACTOR = (288.0D-30)*(PI**2)*(XTANG**4)*XFMOL)

      CHARACTER*15 INFITY, UNDEF
      PARAMETER (INFITY = '        oo     ')
      PARAMETER (UNDEF  = '      undefined')
#include <ibtfun.h>

C     Print title.
C     ------------

      IF (NORGIN .GT. 0) THEN
         CALL AROUND('Length Gauge Specific Rotation in'
     &               //' Requested Origins')
         WRITE(LUPRI,'(25X,A,A,/)') 'Coupled cluster model: ',MODEL
      ELSE
         RETURN
      ENDIF

C     Calculate total mass.
C     ---------------------

      DO I = 1,MXCENT
         ISOTOP(I) = 1
      ENDDO

      TMASS = ZERO
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYM = 0, MAXOPR
            IF (IBTAND(ISYM,ISTBNU(IATOM)) .EQ. 0) THEN
               JATOM = JATOM + 1
               NATTYP = NINT(CHARGE(IATOM))
               IF (NATTYP .NE. 0 .AND. .NOT. NOORBT(IATOM)) THEN
                  AMASS = DISOTP(NATTYP,ISOTOP(JATOM),'MASS')
                  TMASS = TMASS + AMASS
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      WRITE(LUPRI,'(5X,A,/,5X,A,/,5X,A)')
     &' Wavelength               Origin                     [alpha]',
     &'      nm        X/bohr      Y/bohr      Z/bohr deg/[dm g/cm**3]',
     &'---------------------------------------------------------------'
      DO IFREQ = 1,NBLRFR
         IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
            BETA  = -(GLEN(1,IFREQ) + GLEN(2,IFREQ) + GLEN(3,IFREQ))
     &             /(THREE*BLRFR(IFREQ))
            FRCM1 = BLRFR(IFREQ)*XTKAYS
            ALPHA = FACTOR*FRCM1*FRCM1*BETA/TMASS
            WVLEN = XTNM/BLRFR(IFREQ)
            DO IORGIN = 0,NORGIN
               ALPHO = ALPHA
               DO I = 1,3
                  ALPHO = ALPHO + ORGIN(I,IORGIN)*DELTA(I,IFREQ)
               ENDDO
           WRITE(LUPRI,'(5X,F11.3,1X,F11.5,1X,F11.5,1X,F11.5,1X,F15.6)')
     &         WVLEN,ORGIN(1,IORGIN),ORGIN(2,IORGIN),ORGIN(3,IORGIN),
     &         ALPHO
            ENDDO
         ELSE
            WRITE(LUPRI,'(5X,A,1X,A,1X,A,1X,A)')
     &      INFITY,UNDEF,UNDEF,UNDEF
         ENDIF
      ENDDO
      WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'

      RETURN
      END
C  /* Deck cc_oranl */
      SUBROUTINE CC_ORANL(DELTA,IPRINT,MODEL)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Origin-dependence of length gauge OR.
C
#include <implicit.h>
      DIMENSION DELTA(3,*)
      CHARACTER*10 MODEL
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <codata.h>
#include <nuclei.h>
#include <symmet.h>
#include <cclrinf.h>
#include <cbiher.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_ORANL')

      PARAMETER (TINY = 1.0D-7)
      PARAMETER (ZERO = 0.0D0, THREE = 3.0D0)
      PARAMETER (HALF = 0.5D0)
      PARAMETER (HUNDRD = 1.0D2)
      PARAMETER (EHTGHZ = XTHZ*1.0D-9)
C FACTOR = 28800 * pi**2 * bohr**4 * Avogadro, with bohr in cm
C        = 288.0D-30 * pi**2 * bohr**4 * Avogadro, with bohr in Angstrom
      PARAMETER
     & (FACTOR = (288.0D-30)*(PI**2)*(XTANG**4)*XFMOL)

      CHARACTER*15 INFITY, UNDEF
      PARAMETER (INFITY = '         oo    ')
      PARAMETER (UNDEF  = '      undefined')
#include <ibtfun.h>

C     Print title.
C     ------------

      IF (ORGANL) THEN
         CALL AROUND('Origin-Dependence of Length Gauge'
     &               //' Specific Rotation')
         WRITE(LUPRI,'(10X,A,A)')
     &   'Coupled cluster model: ',MODEL
         WRITE(LUPRI,'(10X,A,/)')
     &   'Translation: r -> r - a; [alpha] -> [alpha] + a.Delta'
      ELSE
         RETURN
      ENDIF

C     Calculate total mass.
C     ---------------------

      DO I = 1,MXCENT
         ISOTOP(I) = 1
      ENDDO

      TMASS = ZERO
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYM = 0, MAXOPR
            IF (IBTAND(ISYM,ISTBNU(IATOM)) .EQ. 0) THEN
               JATOM = JATOM + 1
               NATTYP = NINT(CHARGE(IATOM))
               IF (NATTYP .NE. 0 .AND. .NOT. NOORBT(IATOM)) THEN
                  AMASS = DISOTP(NATTYP,ISOTOP(JATOM),'MASS')
                  TMASS = TMASS + AMASS
               ENDIF
            ENDIF
         ENDDO
      ENDDO

C     Get and print delta in units of deg [dm g/cm**3]**(-1) per bohr.
C     ----------------------------------------------------------------

      WRITE(LUPRI,'(5X,A,/,5X,A,/,5X,A)')
     & '     Wavelength                       Delta',
     & '         nm             X               Y               Z',
     & '---------------------------------------------------------------'
      DO IFREQ = 1,NBLRFR
         IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
            FRCM1  = BLRFR(IFREQ)*XTKAYS
            FACTOT = FACTOR*FRCM1*FRCM1/TMASS
            DO I = 1,3
               DELTA(I,IFREQ) = FACTOT*DELTA(I,IFREQ)
            ENDDO
            WVLEN = XTNM/BLRFR(IFREQ)
            WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)')
     &      WVLEN,DELTA(1,IFREQ),DELTA(2,IFREQ),DELTA(3,IFREQ)
         ELSE
            WRITE(LUPRI,'(5X,A,1X,A,1X,A,1X,A)')
     &      INFITY,UNDEF,UNDEF,UNDEF
         ENDIF
      ENDDO
      WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'

      CALL FLSHFO(LUPRI)

      RETURN
      END
C  /* Deck cc_orout */
      SUBROUTINE CC_OROUT(GLEN,GVEL,GTST,WORK,LWORK,IPRINT,MODEL)
C
C     Thomas Bondo Pedersen, April 2003.
C     - based on ABACUS routine lnrout.
C
C     Purpose: Optical activity print, length and/or velocity gauge.
C
#include <implicit.h>
      DIMENSION GLEN(3,*), GVEL(3,*), GTST(3,*), WORK(LWORK)
      CHARACTER*10 MODEL
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <codata.h>
#include <nuclei.h>
#include <symmet.h>
#include <cclrinf.h>
#include <cbiher.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_OROUT')

      PARAMETER (INFOM = 0)
      PARAMETER (TINY = 1.0D-7)
      PARAMETER (ZERO = 0.0D0, THREE = 3.0D0)
      PARAMETER (HALF = 0.5D0)
      PARAMETER (HUNDRD = 1.0D2)
      PARAMETER (EHTGHZ = XTHZ*1.0D-9)
C FACTOR = 28800 * pi**2 * bohr**4 * Avogadro, with bohr in cm
C        = 288.0D-30 * pi**2 * bohr**4 * Avogadro, with bohr in Angstrom
      PARAMETER
     & (FACTOR = (288.0D-30)*(PI**2)*(XTANG**4)*XFMOL)

      CHARACTER*15 INFITY, UNDEF
      PARAMETER (INFITY = '        oo     ')
      PARAMETER (UNDEF  = '      undefined')
#include <ibtfun.h>

C     Print title.
C     ------------

      IF (ORLEN .OR. ORVEL .OR. ORMVE) THEN
         CALL AROUND('Optical Rotation')
         WRITE(LUPRI,'(25X,A,A,/)') 'Coupled cluster model: ',MODEL
         CALL FLSHFO(LUPRI)
      ELSE
         RETURN
      ENDIF

C     Calculate total mass.
C     ---------------------

      DO I = 1,MXCENT
         ISOTOP(I) = 1
      ENDDO

      IF (IPRINT .GE. INFOM) CALL HEADER('Isotopic Masses',-1)

      TMASS = ZERO
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYM = 0, MAXOPR
            IF (IBTAND(ISYM,ISTBNU(IATOM)) .EQ. 0) THEN
               JATOM = JATOM + 1
               NATTYP = NINT(CHARGE(IATOM))
               IF (NATTYP .NE. 0 .AND. .NOT. NOORBT(IATOM)) THEN
                  AMASS = DISOTP(NATTYP,ISOTOP(JATOM),'MASS')
                  TMASS = TMASS + AMASS
                  IF (IPRINT .GE. INFOM) THEN
                     WRITE(LUPRI,'(27X,A6,2X,F12.6)')
     &               NAMDEP(JATOM),AMASS
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF (IPRINT .GE. INFOM) THEN
         WRITE (LUPRI,'(/,23X,A,F12.6,A)') 'Total mass: ',TMASS,' amu'
      ENDIF

C     Print frequency table.
C     ----------------------

      CALL HEADER('Frequencies',-1)
      WRITE(LUPRI,'(5X,A,/,5X,A)')
     & '        Hartree           cm-1             GHz             nm',
     & '---------------------------------------------------------------'
      DO IFREQ = 1,NBLRFR
         FINCM1 = BLRFR(IFREQ)*XTKAYS
         FINGHZ = BLRFR(IFREQ)*EHTGHZ
         IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
            FINNM  = XTNM/BLRFR(IFREQ)
            WRITE(LUPRI,'(5X,F15.8,1X,F15.2,1X,F15.2,1X,F15.2)')
     &      BLRFR(IFREQ),FINCM1,FINGHZ,FINNM
         ELSE
            WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,A)')
     &      BLRFR(IFREQ),FINCM1,FINGHZ,INFITY
         ENDIF
      ENDDO
      WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'

C     Print units.
C     ------------

      CALL HEADER('Units',-1)
      WRITE(LUPRI,'(20X,A)')
     & 'Wavelength                : nm'
      WRITE(LUPRI,'(20X,A)')
     & 'G tensor elements and beta: bohr**4'
      WRITE(LUPRI,'(20X,A)')
     & 'Specific Rotation, [alpha]: deg [dm g/cm**3]**(-1)'
      WRITE(LUPRI,'(20X,A)')
     & 'Molar    Rotation, [M]    : deg [mol dm/cm**3]**(-1)'

C     Print diagonal of G tensor.
C     ---------------------------

      IF (ORLEN) THEN
         CALL HEADER('G Tensor Diagonal, Length Gauge',-1)
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     & '     Wavelength         Gxx             Gyy             Gzz',
     & '---------------------------------------------------------------'
         DO IFREQ = 1,NBLRFR
            IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
               WVLEN = XTNM/BLRFR(IFREQ)
               WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)')
     &         WVLEN,GLEN(1,IFREQ),GLEN(2,IFREQ),GLEN(3,IFREQ)
            ELSE
               WRITE(LUPRI,'(5X,A,1X,F15.6,1X,F15.6,1X,F15.6)')
     &         INFITY,GLEN(1,IFREQ),GLEN(2,IFREQ),GLEN(3,IFREQ)
            ENDIF
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'
      END IF

      IF (ORVEL) THEN
         CALL HEADER('G Tensor Diagonal, Velocity Gauge',-1)
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     & '     Wavelength         Gxx             Gyy             Gzz',
     & '---------------------------------------------------------------'
         DO IFREQ = 1,NBLRFR
            IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
               WVLEN = XTNM/BLRFR(IFREQ)
               WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)')
     &         WVLEN,GVEL(1,IFREQ),GVEL(2,IFREQ),GVEL(3,IFREQ)
            ELSE
               WRITE(LUPRI,'(5X,A,1X,A,1X,A,1X,A)')
     &         INFITY,UNDEF,UNDEF,UNDEF
            ENDIF
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'
      END IF

      IF (ORMVE) THEN
         CALL HEADER('G Tensor Diagonal, Modified Velocity Gauge',-1)
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     & '     Wavelength         Gxx             Gyy             Gzz',
     & '---------------------------------------------------------------'
         DO IFREQ = 1,NBLRFR
            IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
               WVLEN = XTNM/BLRFR(IFREQ)
               WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)')
     &         WVLEN,GTST(1,IFREQ),GTST(2,IFREQ),GTST(3,IFREQ)
            ELSE
               WRITE(LUPRI,'(5X,A,1X,A,1X,A,1X,A)')
     &         INFITY,UNDEF,UNDEF,UNDEF
            ENDIF
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'
      END IF

C     Print specific rotation.
C     ------------------------

      IF (ORLEN) THEN
         CALL HEADER('Specific Rotation, Length Gauge',-1)
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     & '     Wavelength         beta            [alpha]         [M]',
     & '---------------------------------------------------------------'
         DO IFREQ = 1,NBLRFR
            IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
               BETA  = -(GLEN(1,IFREQ) + GLEN(2,IFREQ) + GLEN(3,IFREQ))
     &                /(THREE*BLRFR(IFREQ))
               FRCM1 = BLRFR(IFREQ)*XTKAYS
               ALPHA = FACTOR*FRCM1*FRCM1*BETA/TMASS
               BRMBR = ALPHA*TMASS/HUNDRD
               WVLEN = XTNM/BLRFR(IFREQ)
               WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)')
     &         WVLEN,BETA,ALPHA,BRMBR
            ELSE
               WRITE(LUPRI,'(5X,A,1X,A,1X,A,1X,A)')
     &         INFITY,UNDEF,UNDEF,UNDEF
            ENDIF
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'
      ENDIF

      IF (ORVEL) THEN
         CALL HEADER('Specific Rotation, Velocity Gauge',-1)
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     & '     Wavelength         beta            [alpha]         [M]',
     & '---------------------------------------------------------------'
         DO IFREQ = 1,NBLRFR
            IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
               BETA  = -(GVEL(1,IFREQ) + GVEL(2,IFREQ) + GVEL(3,IFREQ))
     &                /(THREE*BLRFR(IFREQ))
               FRCM1 = BLRFR(IFREQ)*XTKAYS
               ALPHA = FACTOR*FRCM1*FRCM1*BETA/TMASS
               BRMBR = ALPHA*TMASS/HUNDRD
               WVLEN = XTNM/BLRFR(IFREQ)
               WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)')
     &         WVLEN,BETA,ALPHA,BRMBR
            ELSE
               WRITE(LUPRI,'(5X,A,1X,A,1X,A,1X,A)')
     &         INFITY,UNDEF,UNDEF,UNDEF
            ENDIF
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'
      ENDIF

      IF (ORMVE) THEN
         CALL HEADER('Specific Rotation, Modified Velocity Gauge',-1)
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     & '     Wavelength         beta            [alpha]         [M]',
     & '---------------------------------------------------------------'
         DO IFREQ = 1,NBLRFR
            IF (DABS(BLRFR(IFREQ)) .GT. TINY) THEN
               BETA  = -(GTST(1,IFREQ) + GTST(2,IFREQ) + GTST(3,IFREQ))
     &                /(THREE*BLRFR(IFREQ))
               FRCM1 = BLRFR(IFREQ)*XTKAYS
               ALPHA = FACTOR*FRCM1*FRCM1*BETA/TMASS
               BRMBR = ALPHA*TMASS/HUNDRD
               WVLEN = XTNM/BLRFR(IFREQ)
               WRITE(LUPRI,'(5X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)')
     &         WVLEN,BETA,ALPHA,BRMBR
            ELSE
               WRITE(LUPRI,'(5X,A,1X,A,1X,A,1X,A)')
     &         INFITY,UNDEF,UNDEF,UNDEF
            ENDIF
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'
      ENDIF

      WRITE(LUPRI,*)
      CALL FLSHFO(LUPRI)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSOPPRT                             *
*---------------------------------------------------------------------*
C  /* Deck ccdipgrapr */
      SUBROUTINE CCDIPGRAPR(SOPROP,WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <abainf.h>
#include <mxcent.h>
#include <nuclei.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccroper.h>

      DIMENSION SOPROP(NBLRFR,NLROP,2)
      DIMENSION WORK(LWORK)

      LOGICAL DOSYMS(8)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CCDIPGRAPR')

C     Return if nothing to do.
C     ------------------------

      IF (IDIGRA .LE. 0) RETURN

C     Set DOSYM for abacus printing.
C     ------------------------------

      DO I = 1,8
         DOSYMS(I) = DOSYM(I)
         DOSYM(I)  = .TRUE.
      END DO

C     Allocation.
C     -----------

      NCOOR  = 3*NUCDEP
      NDIGRA = NCOOR*3

      KDIGRA = 1
      KCSTRA = KDIGRA + NDIGRA*NBLRFR
      KSCTRA = KCSTRA + NCOOR*NCOOR
      KEND   = KSCTRA + NCOOR*NCOOR
      LWRK   = LWORK  - KEND + 1

      KSCR  = KCSTRA
      KEND0 = KSCR  + NDIGRA
      LWRK0 = LWORK - KCSTRA + 1

      LTST = MIN(LWRK,LWRK0)
      IF (LTST .LT. 0) CALL QUIT('Insuf. mem. in '//SECNAM)

C     Extract dipole gradients (APTs).
C     --------------------------------

      CALL CC_GETDGR(SOPROP,WORK(KDIGRA),IDIGRA)

C     Transpose the APTs.
C     -------------------

      CALL CCTRSP12(WORK(KDIGRA),3*NUCDEP,3,NBLRFR,WORK(KSCR),NDIGRA)

C     Print.
C     ------

      DO IFREQ = 1,NBLRFR
         IF (ABS(BLRFR(IFREQ)) .GT. 1.0D-8) THEN
            WRITE(LUPRI,'(/,8X,A)')
     &      'NB: APTs for non-zero frequency follows!!'
            WRITE(LUPRI,'(8X,A,1P,D14.6,A,/)')
     &      'Frequency: ',BLRFR(IFREQ),' Hartree'
         END IF
         KOFF = KDIGRA + NDIGRA*(IFREQ - 1)
         CALL FCPRI(WORK(KOFF),'APT',WORK(KCSTRA),WORK(KSCTRA)) 
         CALL APTPOP(WORK(KOFF),WORK(KEND),
     &               WORK(KCSTRA),WORK(KSCTRA),
     &               WORK(KEND),WORK(KEND),LWRK)
      END DO

C     Restore DOSYM.
C     --------------

      DO I = 1,8
         DOSYM(I) = DOSYMS(I)
      END DO

      RETURN
      END
C  /* Deck cctrsp12 */
      SUBROUTINE CCTRSP12(XMAT,N1,N2,N3,SCR,LSCR)
C
#include <implicit.h>
      DIMENSION XMAT(N1*N2*N3)
      DIMENSION SCR(LSCR)

      N12 = N1*N2
      IF (N12*N3 .LE. 0) RETURN
      IF (LSCR .LT. N12) CALL QUIT('Insufficient memory in CCTRSP12')

      DO I3 = 1,N3

         KOFFX = N12*(I3 - 1) + 1
         CALL DCOPY(N12,XMAT(KOFFX),1,SCR,1)

         DO I2 = 1,N2
            DO I1 = 1,N1
               I12 = N1*(I2 - 1) + I1
               I21 = N2*(I1 - 1) + I2
               XMAT(KOFFX+I21-1) = SCR(I12)
            END DO
         END DO

      END DO

      RETURN
      END
C  /* Deck cc_getdgr */
      SUBROUTINE CC_GETDGR(SOPROP,PROP,IPROP)
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccroper.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <chrxyz.h>
#include <nuclei.h>

      DIMENSION SOPROP(NBLRFR,NLROP,2)
      DIMENSION PROP(3*NUCDEP,3,NBLRFR)
      INTEGER   IDIP(3)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      IOFFAX = 0
      DO IREP = 0,MAXREP

         NAXIS = NAXREP(IREP,1)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'IPTXYZ: ',
     &      (IPTXYZ(I,IREP,1), I=1,NAXIS)
         END IF

         DO I = 1,NAXIS 
            IDIP(IPTXYZ(I,IREP,1)) = IOFFAX + I 
         END DO

         IOFFAX = IOFFAX + NAXIS

         IF (LOCDBG) THEN
            WRITE (LUPRI,'(23X,3(2A,13X))')
     &      ('E',CHRXYZ(-IPTXYZ(I,IREP,1)),I=1,NAXIS)
            WRITE(LUPRI,*) 'Temporary IDIP for sym. ',IREP+1,':'
            WRITE(LUPRI,*) (IDIP(IPTXYZ(I,IREP,1)), I=1,NAXIS)
         END IF

      END DO

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'Final IDIP: ',(IDIP(I),I=1,3)
      ENDIF

      IOPER = IPROP - 1
      DO IAXIS = 1,3
         DO ICOOR = 1,3*NUCDEP

            IOPER = IOPER + 1

            ISYMA = ISYOPR(IALROP(IOPER))
            ISYMB = ISYOPR(IBLROP(IOPER))

            ISAMA = ISYMAT(IALROP(IOPER))
            ISAMB = ISYMAT(IBLROP(IOPER))

            ISAPROP = ISAMA * ISAMB
            SIGN    = DBLE(ISAPROP)
            IF (ISAPROP.EQ.0) SIGN = +ONE

            IF (ISYMA.EQ.ISYMB) THEN
               DO IFREQ = 1,NBLRFR
                  PROP(ICOOR,IDIP(IAXIS),IFREQ) = 
     &            -0.5D0*(SOPROP(IFREQ,IOPER,1)
     &            + SIGN* SOPROP(IFREQ,IOPER,2))
               END DO
            ELSE
               DO IFREQ = 1,NBLRFR
                  PROP(ICOOR,IDIP(IAXIS),IFREQ) = 0.0D0
               END DO
            END IF

            IF (LOCDBG) THEN
               WRITE(LUPRI,*) 'LABELA,LABELB: ',
     &                   LBLOPR(IALROP(IOPER)),',',LBLOPR(IBLROP(IOPER))
               DO IFREQ = 1,NBLRFR
                  WRITE(LUPRI,*) '   ICOOR,IAXIS,FREQ,LR: ',
     &            ICOOR,IAXIS,BLRFR(IFREQ),PROP(ICOOR,IDIP(IAXIS),IFREQ)
               END DO
            END IF

         END DO
      END DO

      RETURN
      END
