C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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

      KRESULT = 1
      KFTRAN  = KRESULT + 2 * NBSOP
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
     &               MXVEC, WORK(KEND0), LEND0 )

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
     &  NBSOP,' quadratic response func.:', SECOND() - TIM0

*---------------------------------------------------------------------*
* print frequency-dependent polarizabilities:
*---------------------------------------------------------------------*

      CALL  CCSOPPRT(WORK(KRESULT))

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
       SUBROUTINE CCSOPPRT(SOPROP)
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
      CHARACTER*10 MODEL, RLXA, RLXB
      CHARACTER*80 STRING
      INTEGER ISYMA, ISYMB, ISAMA, ISAMB, ISAPROP, IRELAX, M
      INTEGER IFREQ, IOPER, ICMP, ITURN, MXTURN, ITST, ITST2, IPAR


#if defined (SYS_CRAY)
      REAL SOPROP(NBLRFR,NLROP,2)
      REAL SIGN, ANI, ALPHA, PROP, ERROR, ALPHASI, ANISI
      REAL ALFSQ, BETSQ, APAR, APEN, XKAPPA
      REAL TENS(3,3), PVALUES(3), PAXIS(3,3)
      REAL ZERO, HALF, ONE, DPOLAUTSI
#else
      DOUBLE PRECISION SOPROP(NBLRFR,NLROP,2)
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
