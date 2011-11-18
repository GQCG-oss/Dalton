!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
c /* deck ccexlr */
*=====================================================================*
       SUBROUTINE CC_EXLR(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Excited state linear response section:
*
*        --  excited state linear response properties
*        --  two-photon transition moments between two excited states
*

*     Written by Christof Haettig summer 1997.
*     Some restructuring and updates for CC3, october 2003, C. Haettig
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <dummy.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccexlrinf.h>
#include <ccroper.h>
#include <ccr1rsp.h>
#include <ccer1rsp.h>
#include <ccel1rsp.h>
#include <ccn2rsp.h>
#include <cco2rsp.h>
#include <cclists.h>
#include <second.h>

* local parameters:
      CHARACTER*(16) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCEXLR> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0, TIM1, TIMG, TIMF, TIMB
      REAL TIMFA, TIMAA, TIMEA, TIMO
      REAL ZERO
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0, TIM1, TIMG, TIMF, TIMB
      DOUBLE PRECISION TIMFA, TIMAA, TIMEA, TIMO
      DOUBLE PRECISION ZERO
#endif
      PARAMETER ( ZERO = 0.0d0 )

      INTEGER NBEXLR, MXTRAN, MXVEC 
      INTEGER MXGTRAN, MXFTRAN, MXF1TRAN, MXFATRAN, MXEATRAN, MXOTRAN
      INTEGER MXGDOTS, MXFDOTS, MXF1DOTS, MXFADOTS, MXEADOTS, MXODOTS
      INTEGER NGTRAN,  NFTRAN,  NF1TRAN,  NFATRAN,  NEATRAN,  NOTRAN
      INTEGER KGTRAN,  KFTRAN,  KF1TRAN,  KFATRAN,  KEATRAN,  KOTRAN
      INTEGER KGDOTS,  KFDOTS,  KF1DOTS,  KFADOTS,  KEADOTS,  KODOTS
      INTEGER KGCONS,  KFCONS,  KF1CONS,  KFACONS,  KEACONS,  KOCONS
      INTEGER NAATRAN, KAATRAN, KAADOTS,  KAACONS
      INTEGER KEND0, LEND0, KEXLRPRP, IOPT, IORDER

* external functions

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<    OUTPUT FROM COUPLED CLU',
     &                               'STER EXCITED STATE   >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<             LINEAR RESPONSE',
     &                                ' SECTION            >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CCEXLR called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CCEXLR Workspace:',LWORK
  
      TIM0 = SECOND()
*---------------------------------------------------------------------*
* allocate & initialize work space for polarizabilities
*---------------------------------------------------------------------*
      NBEXLR = 2 * NEXLRST * NEXLROPER * NEXLRFREQ

      MXTRAN  = NLRTLBL * MAX(NLRTLBL,NER1LBL,NEL1LBL,NQRN2)
      MXVEC   = MAX(NLRTLBL,NER1LBL,NEL1LBL,NO2LBL,NQRN2)

      MXGTRAN  = MXDIM_GTRAN  * MXTRAN
      MXFTRAN  = MXDIM_FTRAN  * MXTRAN
      MXF1TRAN = MXDIM_FTRAN  * MXTRAN
      MXFATRAN = MXDIM_FATRAN * MXTRAN
      MXEATRAN = MXDIM_XEVEC  * MXTRAN
      MXOTRAN  = 1 * MXTRAN

      MXGDOTS  = MXVEC * MXTRAN
      MXFDOTS  = MXVEC * MXTRAN
      MXF1DOTS = MXVEC * MXTRAN
      MXFADOTS = MXVEC * MXTRAN
      MXEADOTS = MXVEC * MXTRAN
      MXODOTS  = MXVEC * MXTRAN

      KEXLRPRP= 1
      KGTRAN  = KEXLRPRP+ 2 * NBEXLR
      KGDOTS  = KGTRAN  + MXGTRAN
      KGCONS  = KGDOTS  + MXGDOTS
      KFTRAN  = KGCONS  + MXGDOTS
      KFDOTS  = KFTRAN  + MXFTRAN
      KFCONS  = KFDOTS  + MXFDOTS
      KF1TRAN = KFCONS  + MXFDOTS
      KF1DOTS = KF1TRAN + MXF1TRAN
      KF1CONS = KF1DOTS + MXF1DOTS
      KFATRAN = KF1CONS + MXF1DOTS
      KFADOTS = KFATRAN + MXFATRAN
      KFACONS = KFADOTS + MXFADOTS
      KAATRAN = KFACONS + MXFADOTS
      KAADOTS = KAATRAN + MXTRAN * MXDIM_XEVEC
      KAACONS = KAADOTS + MXVEC  * MXTRAN
      KEATRAN = KAACONS + MXVEC  * MXTRAN
      KEADOTS = KEATRAN + MXEATRAN
      KEACONS = KEADOTS + MXEADOTS
      KOTRAN  = KEACONS + MXEADOTS
      KODOTS  = KOTRAN  + MXOTRAN
      KOCONS  = KODOTS  + MXODOTS
      KEND0   = KOCONS  + MXODOTS
      LEND0   = LWORK - KEND0

      IF (LEND0.LT.0) THEN
        WRITE (LUPRI,*) 'KEND0,LEND0:',KEND0,LEND0
        CALL QUIT('Insufficient memory in CCEXLR.')
      END IF

      CALL DZERO(WORK,KEND0-1)

*---------------------------------------------------------------------*
* set up lists for G, F, F{A} transformations etc.:
*---------------------------------------------------------------------*
      CALL CCEXLR_SETUP(MXTRAN, MXVEC,
     &            WORK(KGTRAN), WORK(KGDOTS), WORK(KGCONS), NGTRAN,
     &            WORK(KFTRAN), WORK(KFDOTS), WORK(KFCONS), NFTRAN,
     &            WORK(KF1TRAN),WORK(KF1DOTS),WORK(KF1CONS),NF1TRAN,
     &            WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),NFATRAN,
     &            WORK(KAATRAN),WORK(KAADOTS),WORK(KAACONS),NAATRAN,
     &            WORK(KEATRAN),WORK(KEADOTS),WORK(KEACONS),NEATRAN,
     &            WORK(KOTRAN), WORK(KODOTS), WORK(KOCONS), NOTRAN,
     &            WORK(KEXLRPRP),NBEXLR, .FALSE.  )

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IOPT = 5
      CALL CC_GMATRIX('LE ','R1 ','RE ','R1 ',NGTRAN, MXVEC,
     &              WORK(KGTRAN),WORK(KGDOTS),WORK(KGCONS),
     &              WORK(KEND0), LEND0, IOPT )

      TIMG = SECOND() - TIM1
 
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',NGTRAN,' G matrix transformations:',TIMG
      CALL FLSHFO(LUPRI) 
      
*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF (.NOT. USE_EL1) THEN
        IOPT = 5
        CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'LE ','ER1',IOPT,'R1 ',
     &                  WORK(KFDOTS),WORK(KFCONS),MXVEC,
     &                  WORK(KEND0), LEND0)
      ELSE
        IOPT = 5
        CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'EL1','RE ',IOPT,'R1 ',
     &                  WORK(KFDOTS),WORK(KFCONS),MXVEC,
     &                  WORK(KEND0), LEND0)
      END IF

      TIMF = SECOND() - TIM1
 
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',NFTRAN,' F matrix transformations:',TIMF
      CALL FLSHFO(LUPRI) 
      
*---------------------------------------------------------------------*
* calculate more F matrix contributions:
*---------------------------------------------------------------------*
      TIMB = ZERO

      IF (.NOT. USE_O2) THEN
        TIM1 = SECOND()

        IOPT = 5
        CALL CC_FMATRIX(WORK(KF1TRAN),NF1TRAN,'N2 ','R1 ',IOPT,'R1 ',
     &                  WORK(KF1DOTS),WORK(KF1CONS),MXVEC,
     &                  WORK(KEND0), LEND0)

        TIMB = SECOND() - TIM1
 
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &    '>>> Time used for',NF1TRAN,' F matrix transformations:',TIMB
        CALL FLSHFO(LUPRI) 
      END IF

*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL CCQR_FADRV('LE ','o1 ','RE ','R1 ',NFATRAN, MXVEC,
     &                 WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),
     &                 WORK(KEND0), LEND0, 'DOTP' )

      TIMFA = SECOND() - TIM1
 
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NFATRAN,' F{O} matrix transformat.:',TIMFA
      CALL FLSHFO(LUPRI) 

*---------------------------------------------------------------------*
* calculate A{O} matrix contributions:
*---------------------------------------------------------------------*
      TIMAA = ZERO

      IF (.NOT. USE_O2) THEN
        TIM1 = SECOND()

        IOPT   = 5
        IORDER = 1
        CALL CC_XIETA(WORK(KAATRAN), NAATRAN, IOPT, IORDER, 'N2 ',
     &                '---',IDUMMY,        DUMMY,
     &                'R1 ',WORK(KAADOTS),WORK(KAACONS),
     &                .FALSE.,MXVEC, WORK(KEND0), LEND0 )

        TIMAA = SECOND() - TIM1
 
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',NAATRAN,' A{O} matrix transformat.:',TIMAA
        CALL FLSHFO(LUPRI) 
      END IF

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF (.NOT. USE_EL1) THEN
        IOPT   = 5
        IORDER = 1
        CALL CC_XIETA( WORK(KEATRAN), NEATRAN, IOPT, IORDER, 'LE ',
     &                 '---',IDUMMY,       DUMMY,
     &                 'ER1',WORK(KEADOTS),WORK(KEACONS),
     &                 .FALSE.,MXVEC, WORK(KEND0), LEND0 ) 
      ELSE
        IOPT   = 5
        IORDER = 1
        CALL CC_XIETA( WORK(KEATRAN), NEATRAN, IOPT, IORDER, 'EL1',
     &                 '---',IDUMMY,       DUMMY,
     &                 'RE ',WORK(KEADOTS),WORK(KEACONS),
     &                 .FALSE.,MXVEC, WORK(KEND0), LEND0 ) 
      END IF

      TIMEA = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NEATRAN,' ETA{O} vector calculat.: ',TIMEA
      CALL FLSHFO(LUPRI)  
*---------------------------------------------------------------------*
* calculate N2 x O2 dot products:
*---------------------------------------------------------------------*
      TIMO = ZERO

      IF (USE_O2) THEN 
        TIM1 = SECOND()

        CALL CC_DOTDRV('N2 ','O2 ',NOTRAN,MXVEC,
     &                 WORK(KOTRAN), WORK(KODOTS), WORK(KOCONS),
     &                 WORK(KEND0), LEND0 )

        TIMO = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &     '>>> Time used for',NOTRAN,' N2 x O2 dot products: ', TIMO
        CALL FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
* collect contributions and add them excited state LR properties
*---------------------------------------------------------------------*
      CALL CCEXLR_SETUP(MXTRAN, MXVEC,
     &            WORK(KGTRAN), WORK(KGDOTS), WORK(KGCONS), NGTRAN,
     &            WORK(KFTRAN), WORK(KFDOTS), WORK(KFCONS), NFTRAN,
     &            WORK(KF1TRAN),WORK(KF1DOTS),WORK(KF1CONS),NF1TRAN,
     &            WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),NFATRAN,
     &            WORK(KAATRAN),WORK(KAADOTS),WORK(KAACONS),NAATRAN,
     &            WORK(KEATRAN),WORK(KEADOTS),WORK(KEACONS),NEATRAN,
     &            WORK(KOTRAN), WORK(KODOTS), WORK(KOCONS), NOTRAN,
     &            WORK(KEXLRPRP),NBEXLR, .TRUE.  )


*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     & NBEXLR,' excited state linear response func.:', SECOND() - TIM0 

*---------------------------------------------------------------------*
* print output & return:
*---------------------------------------------------------------------*

      CALL  CCEXLRPRT(WORK(KEXLRPRP))

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_HYPPOL                            *
*=====================================================================*

c /* deck ccexlrprt */
*=====================================================================*
       SUBROUTINE CCEXLRPRT(EXLRPRP)
*---------------------------------------------------------------------*
*
*    Purpose: print output for excited state linear response section
*
*
*     Written by Christof Haettig in Juli 1997.
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
#include <ccexlrinf.h>
#include <ccexci.h>
#include <ccroper.h>


      CHARACTER*5  BLANKS
      CHARACTER*80 STRING
      LOGICAL LTWOPHOT
      INTEGER ISYMA, ISYMB, ISYMSI, ISYMSF, ISTATI, ISTATF
      INTEGER IFREQ, IOPER, IDXS, IEXCII, IEXCIF


#if defined (SYS_CRAY)
      REAL EXLRPRP(NEXLRFREQ,NEXLROPER,NEXLRST,2)
      REAL HALF, FREQA, FREQB, EIGVI, EIGVF
#else
      DOUBLE PRECISION EXLRPRP(NEXLRFREQ,NEXLROPER,NEXLRST,2)
      DOUBLE PRECISION HALF, FREQA, FREQB, EIGVI, EIGVF
#endif
      PARAMETER (HALF = 0.5d0)

*---------------------------------------------------------------------*
* initialize flag for two photon transition moments:
*---------------------------------------------------------------------*
      LTWOPHOT = .FALSE.

*---------------------------------------------------------------------*
* print header for excited state polarizabilities:
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING =' RESULTS FOR EXCITED STATES LINEAR RESPONSE PROPERTIES '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:55)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:55)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:55)//BLANKS )
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:55)//BLANKS )
      ELSE
         CALL QUIT('CCEXLRPRT called for an unknown '//
     &        'Coupled Cluster model.')
      END IF

      DO IDXS = 1, NEXLRST
         ISYMSI = IELRSYM(IDXS,1)
         ISYMSF = IELRSYM(IDXS,2)
         ISTATI = IELRSTA(IDXS,1)
         ISTATF = IELRSTA(IDXS,2)
         IEXCII = ISYOFE(ISYMSI) + ISTATI
         IEXCIF = ISYOFE(ISYMSF) + ISTATF
         EIGVI  = EIGVAL(IEXCII)
         EIGVF  = EIGVAL(IEXCIF)

      IF (IEXCII.NE.IEXCIF) THEN

          LTWOPHOT = .TRUE.

      ELSE
C     IF (IEXCII.EQ.IEXCIF) THEN
         WRITE(STRING,'(A,I2,A,I2,3X,A,F12.8,A)') 
     &        ' State number',ISTATI,
     &        ' in symmetry class',ISYMSI,
     &        ' (excitation energy: ',EIGVI,')'
         CALL AROUND(STRING(1:72))

         IF (IPREXLR.GT.5) THEN
           WRITE(LUPRI,'(/1X,2(1X,A,7X),5X,A,10X,A,/,95("-"))')
     &      'A operator','B operator','  alpha','(asy. Resp.)'
         ELSE
           WRITE(LUPRI,'(/1X,2(1X,A,7X),4X,A,/,60("-"))')
     &      'A operator','B operator','  alpha'
         END IF

      DO IOPER = 1, NEXLROPER
         ISYMA = ISYOPR(IAEXLROP(IOPER))
         ISYMB = ISYOPR(IBEXLROP(IOPER))

         IFREQ = 1
         IF (ISYMA.EQ.ISYMB) THEN
          IF (IPREXLR.GT.5) THEN
            WRITE(LUPRI,'(/2X,2(A8,F7.4,3X),G18.10," (",G18.10,")")')
     &        LBLOPR(IAEXLROP(IOPER)),-BEXLRFR(IFREQ),
     &        LBLOPR(IBEXLROP(IOPER)),+BEXLRFR(IFREQ),
     &        -HALF*(EXLRPRP(IFREQ,IOPER,IDXS,1)
     &               +EXLRPRP(IFREQ,IOPER,IDXS,2)),
     &        -HALF*(EXLRPRP(IFREQ,IOPER,IDXS,1)
     &               -EXLRPRP(IFREQ,IOPER,IDXS,2))
          ELSE
            WRITE(LUPRI,'(/2X,2(A8,F7.4,3X),G16.8)')
     &        LBLOPR(IAEXLROP(IOPER)),-BEXLRFR(IFREQ),
     &        LBLOPR(IBEXLROP(IOPER)),+BEXLRFR(IFREQ),
     &        -HALF*(EXLRPRP(IFREQ,IOPER,IDXS,1)
     &               +EXLRPRP(IFREQ,IOPER,IDXS,2))
          ENDIF
         ELSE
          IF (IPREXLR.GT.5) THEN
           WRITE(LUPRI,'(/2X,2(A8,F7.4,3X),7X,A,8X," (",9X,A,10X,")")')
     &       LBLOPR(IAEXLROP(IOPER)),-BEXLRFR(IFREQ),
     &       LBLOPR(IBEXLROP(IOPER)),+BEXLRFR(IFREQ),
     &       '---',
     &       '---'
          ELSE 
           WRITE(LUPRI,'(/2X,2(A8,F7.4,3X),6X,A,7X)')
     &       LBLOPR(IAEXLROP(IOPER)),-BEXLRFR(IFREQ),
     &       LBLOPR(IBEXLROP(IOPER)),+BEXLRFR(IFREQ),
     &       '---'
          END IF
         END IF

         DO IFREQ = 2, NEXLRFREQ
          IF (ISYMA.EQ.ISYMB) THEN
           IF (IPREXLR.GT.5) THEN
            WRITE(LUPRI,'(2X,2(8X,F7.4,3X),G18.10," (",G18.10,")")')
     &        -BEXLRFR(IFREQ), BEXLRFR(IFREQ),
     &        -HALF*(EXLRPRP(IFREQ,IOPER,IDXS,1)
     &               +EXLRPRP(IFREQ,IOPER,IDXS,2)),
     &        -HALF*(EXLRPRP(IFREQ,IOPER,IDXS,1)
     &               -EXLRPRP(IFREQ,IOPER,IDXS,2))
           ELSE
            WRITE(LUPRI,'(2X,2(8X,F7.4,3X),G16.8)')
     &        -BEXLRFR(IFREQ), BEXLRFR(IFREQ),
     &        -HALF*(EXLRPRP(IFREQ,IOPER,IDXS,1)
     &               +EXLRPRP(IFREQ,IOPER,IDXS,2))
           END IF
          END IF
         END DO

      END DO
      END IF
      END DO

      IF (.NOT.LTWOPHOT) RETURN

*---------------------------------------------------------------------*
* print header for two-photon matrix elements between excited states:
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING =' RESULTS FOR TWO-PHOTON TRANSITION MATRIX ELEMENTS '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:51)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:51)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:51)//BLANKS )
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:55)//BLANKS )
      ELSE
         CALL QUIT('CCEXLRPRT called for an unknown '//
     &        'Coupled Cluster model.')
      END IF

      IF ( HALFFR .AND. NEXLRFREQ.NE.1 ) THEN
        WRITE(LUPRI,*) 'error in CCEXLRPRT: HALFFR option is ',
     &             'incompatible with a frequency list.' 
        CALL QUIT('error in CCEXLRPRT.')
      END IF

      WRITE(LUPRI,'(/,/,"+",112("-"),"+")')
      WRITE(LUPRI,'(3A,/,"|",112(" "),"|",/,3A)') 
     &   '|     STATE I          STATE F    ',
     &   '   OPERATOR A       OPERATOR B     ',
     &   '        MOMENTS                             |',
     &   '| SYM IDX  EIGVAL  SYM IDX  EIGVAL',
     &   '  LABEL    FREQ    LABEL    FREQ  ',
     &   '  M{i<-f}(w_B)  M{f<-i}(-w_B)  M{if}*M{fi}   |'
      WRITE(LUPRI,'("+",112("-"),"+")') 

      DO IDXS = 1, NEXLRST
         ISYMSI = IELRSYM(IDXS,1)
         ISYMSF = IELRSYM(IDXS,2)
         ISTATI = IELRSTA(IDXS,1)
         ISTATF = IELRSTA(IDXS,2)
         IEXCII = ISYOFE(ISYMSI) + ISTATI
         IEXCIF = ISYOFE(ISYMSF) + ISTATF
         EIGVI  = EIGVAL(IEXCII)
         EIGVF  = EIGVAL(IEXCIF)


      IF (IEXCII.NE.IEXCIF) THEN
         WRITE(LUPRI,'("|",112(" "),"|")') 

      DO IOPER = 1, NEXLROPER
         ISYMA = ISYOPR(IAEXLROP(IOPER))
         ISYMB = ISYOPR(IBEXLROP(IOPER))

         IFREQ = 1
         FREQB = BEXLRFR(IFREQ)
         IF ( HALFFR ) FREQB = HALF * (EIGVI-EIGVF)
         FREQA = EIGVI - EIGVF - FREQB
         IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMSI,ISYMSF)) THEN
          IF (IOPER.EQ.1) THEN
            WRITE(LUPRI, '("|",2(I3,1X,I3,F9.4,1X),
     &               (A7,F8.4,2X),(A7,F8.4,1X),3(1X,G14.8),"|")')
     &       ISYMSI,ISTATI,EIGVI,ISYMSF,ISTATF,EIGVF,
     &       LBLOPR(IAEXLROP(IOPER)),FREQA,
     &       LBLOPR(IBEXLROP(IOPER)),FREQB,
     &       EXLRPRP(IFREQ,IOPER,IDXS,1),EXLRPRP(IFREQ,IOPER,IDXS,2),
     &       EXLRPRP(IFREQ,IOPER,IDXS,1)*EXLRPRP(IFREQ,IOPER,IDXS,2)
          ELSE
            WRITE(LUPRI, '("|",2(3X,1X,3X,9X,1X),
     &               (A7,F8.4,2X),(A7,F8.4,1X),3(1X,G14.8),"|")')
     &       LBLOPR(IAEXLROP(IOPER)),FREQA,
     &       LBLOPR(IBEXLROP(IOPER)),FREQB,
     &       EXLRPRP(IFREQ,IOPER,IDXS,1),EXLRPRP(IFREQ,IOPER,IDXS,2),
     &       EXLRPRP(IFREQ,IOPER,IDXS,1)*EXLRPRP(IFREQ,IOPER,IDXS,2)
          END IF
         ELSE
          IF (IOPER.EQ.1) THEN
            WRITE(LUPRI, '("|",2(I3,1X,I3,F9.4,1X),
     &               (A7,F8.4,2X),(A7,F8.4,1X),3A14,"   |")')
     &        ISYMSI,ISTATI,EIGVI,ISYMSF,ISTATF,EIGVF,
     &        LBLOPR(IAEXLROP(IOPER)),FREQA,
     &        LBLOPR(IBEXLROP(IOPER)),FREQB,
     &        '      ---     ', '      ---     ','      ---     '
          ELSE
            WRITE(LUPRI, '("|",2(3X,1X,3X,9X,1X),
     &               (A7,F8.4,2X),(A7,F8.4,1X),3A14,"   |")')
     &        LBLOPR(IAEXLROP(IOPER)),FREQA,
     &        LBLOPR(IBEXLROP(IOPER)),FREQB,
     &        '      ---     ', '      ---     ','      ---     '
          END IF
         END IF

         DO IFREQ = 2, NEXLRFREQ
          FREQB = BEXLRFR(IFREQ)
          FREQA = EIGVI - EIGVF - FREQB
          IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMSI,ISYMSF)) THEN
            WRITE(LUPRI, '("| ",32X,
     &            (7X,F8.4,2X),(7X,F8.4,1X),3(1X,G14.8),"|")')
     &          FREQA, FREQB,
     &        EXLRPRP(IFREQ,IOPER,IDXS,1),EXLRPRP(IFREQ,IOPER,IDXS,2),
     &        EXLRPRP(IFREQ,IOPER,IDXS,1)*EXLRPRP(IFREQ,IOPER,IDXS,2)
          END IF
         END DO

      END DO
      END IF
      END DO

      WRITE(LUPRI,'("|",112(" "),"|")') 
      WRITE(LUPRI,'("+",112("-"),"+")') 

      RETURN
      END
*---------------------------------------------------------------------*
*               END OF SUBROUTINE CCEXLRPRT                           *
*---------------------------------------------------------------------*
c /* deck ccexlr_setup */
*=====================================================================*
      SUBROUTINE CCEXLR_SETUP(MXTRAN,  MXVEC,
     &                        IGTRAN,  IGDOTS,  WG,  NGTRAN,
     &                        IFTRAN,  IFDOTS,  WF,  NFTRAN,
     &                        IF1TRAN, IF1DOTS, F1CONS, NF1TRAN,
     &                        IFATRAN, IFADOTS, WFA, NFATRAN,
     &                        IAATRAN, IAADOTS, WAA, NAATRAN,
     &                        IEATRAN, IEADOTS, WEA, NEATRAN,
     &                        IOTRAN,  IODOTS,  WO,  NOTRAN,
     &                        EXLRPRP, MXPROP,  LADD          )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CCEXLR section
*                - list of G matrix transformations 
*                - list of F matrix transformations 
*                - list of F{O} matrix transformations 
*                - list of ETA{O} vector calculations 
*                - list of dot products of N2 and O2 vectors
*
*     LADD = .FALSE.  --> build lists of contributions
*     LADD = .TRUE.   --> add contributions up to properties
*
*     Written by Christof Haettig, july 1997.
*     Some restructuring and updates for CC3, october 2003, C. Haettig
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccexlrinf.h>
#include <ccexci.h>
#include <ccroper.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(22) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCEXLR_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD

      INTEGER MXVEC, MXTRAN, MXPROP

      INTEGER IGTRAN(MXDIM_GTRAN,MXTRAN)
      INTEGER IGDOTS(MXVEC,MXTRAN)

      INTEGER IFTRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)
  
      INTEGER IF1TRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IF1DOTS(MXVEC,MXTRAN)
  
      INTEGER IFATRAN(MXDIM_FATRAN,MXTRAN)
      INTEGER IFADOTS(MXVEC,MXTRAN)

      INTEGER IAATRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IAADOTS(MXVEC,MXTRAN)

      INTEGER IEATRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IEADOTS(MXVEC,MXTRAN)

      INTEGER IOTRAN(MXTRAN)
      INTEGER IODOTS(MXVEC,MXTRAN)

      INTEGER NGTRAN, NFTRAN, NFATRAN, NEATRAN, NOTRAN, NEXLRPROP,
     &        NF1TRAN, NAATRAN

      CHARACTER*(8) LABELA, LABELB

      LOGICAL LPRJ

      INTEGER ISYMA,ISYMB,ISYMSI,ISYMSF,ISTATI,ISTATF,IEXCII,IEXCIF
      INTEGER IFREQ, IOPER, ISIGN, IDXS

      INTEGER IOPA,IOPB,ITA,ITB,IERA,IERB,IER,IEL,IN2,IO2,IELA,IELB
      INTEGER IVEC, ITRAN, I, IDX

      INTEGER MXG, MXF, MXFA, MXEA, MXO, MXF1VEC, MXAA

#if defined (SYS_CRAY)
      REAL EIGVI, EIGVF, FREQA, FREQB
      REAL EXLRPRP(2*MXPROP)
      REAL WG(MXVEC,MXTRAN)
      REAL WF(MXVEC,MXTRAN)
      REAL F1CONS(MXVEC,MXTRAN)
      REAL WFA(MXVEC,MXTRAN)
      REAL WAA(MXVEC,MXTRAN)
      REAL WEA(MXVEC,MXTRAN)
      REAL WO(MXVEC,MXTRAN)
      REAL GCON, FCON1, FCON2, FACON1, FACON2
      REAL EACON1, EACON2, OCON, F1CON, AACON1, AACON2
      REAL HALF, ZERO
#else
      DOUBLE PRECISION EIGVI, EIGVF, FREQA, FREQB
      DOUBLE PRECISION EXLRPRP(2*MXPROP)
      DOUBLE PRECISION WG(MXVEC,MXTRAN)
      DOUBLE PRECISION WF(MXVEC,MXTRAN)
      DOUBLE PRECISION F1CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION WFA(MXVEC,MXTRAN)
      DOUBLE PRECISION WAA(MXVEC,MXTRAN)
      DOUBLE PRECISION WEA(MXVEC,MXTRAN)
      DOUBLE PRECISION WO(MXVEC,MXTRAN)
      DOUBLE PRECISION GCON, FCON1, FCON2, FACON1, FACON2
      DOUBLE PRECISION EACON1, EACON2, OCON, F1CON, AACON1, AACON2
      DOUBLE PRECISION HALF, ZERO
#endif
      PARAMETER ( HALF = 0.5d0, ZERO = 0.0d0 )


* external functions:
      INTEGER IR1TAMP
      INTEGER IER1AMP
      INTEGER IEL1AMP
      INTEGER IN2AMP
      INTEGER IRHSR2


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN
        NGTRAN  = 0
        NFTRAN  = 0
        NF1TRAN  = 0
        NFATRAN = 0
        NAATRAN = 0
        NEATRAN = 0
        NOTRAN  = 0

        DO ITRAN = 1, MXTRAN
          DO I = 1, MXDIM_XEVEC
            IEATRAN(I,ITRAN) = 0
            IAATRAN(I,ITRAN) = 0
          END DO
          IEATRAN(3,ITRAN) = -1
          IEATRAN(4,ITRAN) = -1
          IAATRAN(3,ITRAN) = -1
          IAATRAN(4,ITRAN) = -1
        END DO
      END IF

      MXG  = 0
      MXF  = 0
      MXF1VEC  = 0
      MXFA = 0
      MXAA = 0
      MXEA = 0
      MXO  = 0

      NEXLRPROP = 0
 
      IF ( HALFFR .AND. NEXLRFREQ.NE.1 ) THEN
        WRITE (LUPRI,*) 'error in CCEXLR_SETUP: HALFFR option is',
     &             ' incompatible with a frequency list.' 
        CALL QUIT('error in CCEXLR_SETUP.')
      END IF

*---------------------------------------------------------------------*
* start loop over all excited state linear response properties
*---------------------------------------------------------------------*
 
      DO IOPER = 1, NEXLROPER
        IOPA   = IAEXLROP(IOPER)
        IOPB   = IBEXLROP(IOPER)

        LABELA = LBLOPR(IOPA)
        LABELB = LBLOPR(IOPB)

        ISYMA  = ISYOPR(IOPA)
        ISYMB  = ISYOPR(IOPB)

      DO IDXS = 1, NEXLRST
        ISYMSI = IELRSYM(IDXS,1)
        ISYMSF = IELRSYM(IDXS,2)
        ISTATI = IELRSTA(IDXS,1)
        ISTATF = IELRSTA(IDXS,2)
        IEXCII = ISYOFE(ISYMSI) + ISTATI
        IEXCIF = ISYOFE(ISYMSF) + ISTATF
        EIGVI  = EIGVAL(IEXCII)
        EIGVF  = EIGVAL(IEXCIF)

      IF ( MULD2H(ISYMA,ISYMB) .EQ. MULD2H(ISYMSI,ISYMSF) ) THEN

      DO IFREQ = 1, NEXLRFREQ
        FREQB = BEXLRFR(IFREQ)
        IF (IEXCII.EQ.IEXCIF) THEN
          FREQA  = -FREQB
          LPRJ   = .NOT. NOPROJ
        ELSE
          IF ( HALFFR ) FREQB = HALF * (EIGVI-EIGVF)
          FREQA  = EIGVI - EIGVF -FREQB
          LPRJ   = .FALSE.
        END IF

        NEXLRPROP = NEXLRPROP + 1

      DO ISIGN = +1, -1, -2

       IF (ISIGN.EQ.1) THEN

         ITA  = IR1TAMP(LABELA,.FALSE.,+FREQA,ISYMA)
         ITB  = IR1TAMP(LABELB,.FALSE.,+FREQB,ISYMB)
         IER  = IEXCIF
         IEL  = IEXCII
         IN2  = IN2AMP(IEXCII,-EIGVI,ISYMSI,IEXCIF,EIGVF,ISYMSF)
         IF (USE_O2) THEN
           IO2  = IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                   LABELB,.FALSE.,+FREQB,ISYMB)
         END IF
         IF (.NOT. USE_EL1) THEN
           IERA = IER1AMP(IEXCIF,EIGVF,ISYMSF,LABELA,+FREQA,ISYMA,LPRJ)
           IERB = IER1AMP(IEXCIF,EIGVF,ISYMSF,LABELB,+FREQB,ISYMB,LPRJ)
         ELSE
           IELA = IEL1AMP(IEXCII,EIGVI,ISYMSI,
     &                    LABELA,+FREQA,ISYMA,.FALSE.,LPRJ)
           IELB = IEL1AMP(IEXCII,EIGVI,ISYMSI,
     &                    LABELB,+FREQB,ISYMB,.FALSE.,LPRJ)
         END IF

       ELSE ! switch states indices and signs of the frequencies

         ITA  = IR1TAMP(LABELA,.FALSE.,-FREQA,ISYMA)
         ITB  = IR1TAMP(LABELB,.FALSE.,-FREQB,ISYMB)
         IER  = IEXCII
         IEL  = IEXCIF
         IN2  = IN2AMP(IEXCIF,-EIGVF,ISYMSF,IEXCII,EIGVI,ISYMSI)
         IF (USE_O2) THEN
           IO2  = IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                   LABELB,.FALSE.,-FREQB,ISYMB)
         END IF
         IF (.NOT. USE_EL1) THEN
           IERA = IER1AMP(IEXCII,EIGVI,ISYMSI,LABELA,-FREQA,ISYMA,LPRJ)
           IERB = IER1AMP(IEXCII,EIGVI,ISYMSI,LABELB,-FREQB,ISYMB,LPRJ)
         ELSE
           IELA = IEL1AMP(IEXCIF,EIGVF,ISYMSF,
     &                    LABELA,-FREQA,ISYMA,.FALSE.,LPRJ)
           IELB = IEL1AMP(IEXCIF,EIGVF,ISYMSF,
     &                    LABELB,-FREQB,ISYMB,.FALSE.,LPRJ)
         END IF

       END IF


*---------------------------------------------------------------------*
* set up list of G matrix transformations, 1 permutation
*---------------------------------------------------------------------*
        CALL CC_SETG212(IGTRAN,IGDOTS,MXTRAN,MXVEC,
     &                  IEL,ITA,IER,ITB,ITRAN,IVEC)
        NGTRAN = MAX(NGTRAN,ITRAN)
        MXG    = MAX(MXG,IVEC)
        GCON   = WG(IVEC,ITRAN)

*---------------------------------------------------------------------*
* set up list of F matrix transformations, 2 permutations
*---------------------------------------------------------------------*
        IF (.NOT. USE_EL1) THEN
          CALL CC_SETF12(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                   IEL,IERA,ITB,ITRAN,IVEC)
        ELSE
          CALL CC_SETF12(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                   IELA,IER,ITB,ITRAN,IVEC)
        END IF
        NFTRAN = MAX(NFTRAN,ITRAN)
        MXF    = MAX(MXF,IVEC)
        FCON1  = WF(IVEC,ITRAN)

        IF (.NOT. USE_EL1) THEN
          CALL CC_SETF12(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                   IEL,IERB,ITA,ITRAN,IVEC)
        ELSE
          CALL CC_SETF12(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                   IELB,IER,ITA,ITRAN,IVEC)
        END IF
        NFTRAN = MAX(NFTRAN,ITRAN)
        MXF    = MAX(MXF,IVEC)
        FCON2  = WF(IVEC,ITRAN)

*---------------------------------------------------------------------*
* set up list of F{O} matrix transformations, 2 permutations
*---------------------------------------------------------------------*
        CALL CC_SETFA12(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                  IEL,IOPA,IER,ITB,ITRAN,IVEC)
        NFATRAN = MAX(NFATRAN,ITRAN)
        MXFA   = MAX(MXFA,IVEC)
        FACON1 = WFA(IVEC,ITRAN)

        CALL CC_SETFA12(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                  IEL,IOPB,IER,ITA,ITRAN,IVEC)
        NFATRAN = MAX(NFATRAN,ITRAN)
        MXFA   = MAX(MXFA,IVEC)
        FACON2 = WFA(IVEC,ITRAN)

*---------------------------------------------------------------------*
* set up list of generalized ETA{O} vector calculations, 2 permutations
*---------------------------------------------------------------------*
        IF (.NOT. USE_EL1) THEN
          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IEL,IOPA,0,0,0,0,IERB,ITRAN,IVEC)
        ELSE
          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IELB,IOPA,0,0,0,0,IER,ITRAN,IVEC)
        END IF
        NEATRAN = MAX(NEATRAN,ITRAN)
        MXEA   = MAX(MXEA,IVEC)
        EACON1 = WEA(IVEC,ITRAN)

        IF (.NOT. USE_EL1) THEN
          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IEL,IOPB,0,0,0,0,IERA,ITRAN,IVEC)
        ELSE
          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IELA,IOPB,0,0,0,0,IER,ITRAN,IVEC)
        END IF
        NEATRAN = MAX(NEATRAN,ITRAN)
        MXEA   = MAX(MXEA,IVEC)
        EACON2 = WEA(IVEC,ITRAN)

*---------------------------------------------------------------------*
* set up list of N2 x O2 dot products, 1 permutation
*---------------------------------------------------------------------*
        OCON = ZERO

        IF (USE_O2) THEN
          CALL CC_SETDOT(IOTRAN,IODOTS,MXTRAN,MXVEC,
     &                   IN2,IO2,ITRAN,IVEC)
          NOTRAN = MAX(NOTRAN,ITRAN)
          MXO    = MAX(MXO,IVEC)
          OCON   = WO(IVEC,ITRAN)
        END IF

*---------------------------------------------------------------------*
* set up list of F matrix transformations, 1 permutation
*---------------------------------------------------------------------*
        F1CON = ZERO

        IF (.NOT. USE_O2) THEN
          CALL CC_SETF12(IF1TRAN,IF1DOTS,MXTRAN,MXVEC,
     &                   IN2,ITB,ITA,ITRAN,IVEC)
          NF1TRAN = MAX(NF1TRAN,ITRAN)
          MXF1VEC    = MAX(MXF1VEC,IVEC)
          F1CON   = F1CONS(IVEC,ITRAN)
        END IF

*---------------------------------------------------------------------*
* set up list of generalized Eta{O} vector calculations, 2 permutation
*---------------------------------------------------------------------*
        AACON1 = ZERO
        AACON2 = ZERO
 
        IF (.NOT. USE_O2) THEN
          CALL CC_SETXE('Eta',IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                   IN2,IOPA,0,0,0,0,ITB,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXAA    = MAX(MXAA,IVEC) 
          AACON1  = WAA(IVEC,ITRAN)

          CALL CC_SETXE('Eta',IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                   IN2,IOPB,0,0,0,0,ITA,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXAA    = MAX(MXAA,IVEC) 
          AACON2  = WAA(IVEC,ITRAN)
        END IF

*---------------------------------------------------------------------*
* add contributions up to excited state linear response property:
*---------------------------------------------------------------------*
        IF (LADD) THEN
          IDX =(IDXS-1)*NEXLROPER*NEXLRFREQ+(IOPER-1)*NEXLRFREQ+IFREQ
          IF (ISIGN.EQ.-1) IDX = IDX + NEXLRST*NEXLROPER*NEXLRFREQ

          EXLRPRP(IDX) = GCON + FCON1 + FCON2 + FACON1 + FACON2 +
     &                 EACON1 + EACON2 + OCON + F1CON + AACON1 + AACON2
 
          IF (LOCDBG) THEN
            WRITE (LUPRI,*)
            WRITE (LUPRI,*) MSGDBG, 'IOPER:',IOPER
            WRITE (LUPRI,*) MSGDBG, 'LABELA, LABELB:',LABELA, LABELB
            WRITE (LUPRI,*) MSGDBG, 'FREQA, FREQB:',FREQA,FREQB
            WRITE (LUPRI,*) MSGDBG, 'ISYMSI,ISTATI,EIGVI:',ISYMSI,
     &           ISTATI,EIGVI
            WRITE (LUPRI,*) MSGDBG, 'ISYMSF,ISTATF,EIGVF:',ISYMSF,
     &           ISTATF,EIGVF
            WRITE (LUPRI,*) MSGDBG, 'IDX:  ',IDX
            WRITE (LUPRI,*) MSGDBG, 'EXLRPRP:  ',EXLRPRP(IDX)
            WRITE (LUPRI,*) MSGDBG, 'GCON:',GCON
            WRITE (LUPRI,*) MSGDBG, 'FCON:',FCON1,FCON2
            WRITE (LUPRI,*) MSGDBG, 'FACON:',FACON1,FACON2
            WRITE (LUPRI,*) MSGDBG, 'EACON:',EACON1,EACON2
            WRITE (LUPRI,*) MSGDBG, 'OCON:',OCON
            WRITE (LUPRI,*) MSGDBG, 'F1CON:',F1CON
            WRITE (LUPRI,*) MSGDBG, 'AACON:',AACON1,AACON2
            WRITE (LUPRI,*) MSGDBG, 'SUM:',
     &           GCON+FCON1+FCON2+FACON1+FACON2+EACON1+EACON2+OCON+
     &           F1CON+AACON1+AACON2
            WRITE (LUPRI,*)
          END IF
        END IF

*---------------------------------------------------------------------*
* end loop over all requested excited state linear response properties
*---------------------------------------------------------------------*
      END DO
      END DO
      END IF
      END DO
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      IF (.NOT. LADD) THEN
        WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NEXLRPROP,
     &        ' excited state linear response properties '
        WRITE(LUPRI,'((8X,A,I3,A))') 
     &     ' - ',NGTRAN,  ' generalized G matrix transformations ',
     &     ' - ',NFTRAN,  ' generalized F matrix transformations ',
     &     ' - ',NF1TRAN, ' generalized F matrix transformations ',
     &     ' - ',NFATRAN, ' generalized F{O} matrix transformations ',
     &     ' - ',NAATRAN, ' generalized ETA{O} vecotr calculations ',
     &     ' - ',NEATRAN, ' generalized ETA{O} vector calculations ',
     &     ' - ',NOTRAN,  ' N2 x O2 dot products calculations '
        WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'
      END IF


* G matrix transformations:
      IF (LOCDBG .AND. .NOT.LADD) THEN
        WRITE (LUPRI,*) MSGDBG, 'List of G matrix transformations:'
        DO ITRAN = 1, NGTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IGTRAN(I,ITRAN),I=1,3),(IGDOTS(I,ITRAN),I=1,MXG)
        END DO
        WRITE (LUPRI,*)
      END IF

* F matrix transformations:
      IF (LOCDBG .AND. .NOT.LADD) THEN
        WRITE (LUPRI,*) MSGDBG, 'List of F matrix transformations:'
        DO ITRAN = 1, NFTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IFTRAN(I,ITRAN),I=1,2),(IFDOTS(I,ITRAN),I=1,MXF)
        END DO
        WRITE (LUPRI,*)
      END IF

* more F matrix transformations:
      IF (LOCDBG .AND. .NOT.LADD) THEN
        WRITE (LUPRI,*) MSGDBG, 'List of F matrix transformations:'
        DO ITRAN = 1, NF1TRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IF1TRAN(I,ITRAN),I=1,2),(IF1DOTS(I,ITRAN),I=1,MXF1VEC)
        END DO
        WRITE (LUPRI,*)
      END IF

* F{O} matrix transformations:
      IF (LOCDBG .AND. .NOT.LADD) THEN
        WRITE (LUPRI,*) MSGDBG, 'List of F{O} matrix transformations:'
        DO ITRAN = 1, NFATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (IFATRAN(I,ITRAN),I=1,5),(IFADOTS(I,ITRAN),I=1,MXFA)
        END DO
        WRITE (LUPRI,*)
      END IF

* more ETA{O} vectors calculations:
      IF (LOCDBG .AND. .NOT.LADD) THEN
        WRITE (LUPRI,*) MSGDBG, 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NAATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (IAATRAN(I,ITRAN),I=1,5),(IAADOTS(I,ITRAN),I=1,MXAA)
        END DO
        WRITE (LUPRI,*)
      END IF

* ETA{O} vector calculations:
      IF (LOCDBG .AND. .NOT.LADD) THEN
        WRITE (LUPRI,*) MSGDBG, 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NEATRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IEATRAN(I,ITRAN),I=1,2),(IEADOTS(I,ITRAN),I=1,MXEA)
        END DO
        WRITE (LUPRI,*)
      END IF

* N2 x O2 vector dot products:
      IF (LOCDBG .AND. .NOT.LADD) THEN
        WRITE (LUPRI,*) MSGDBG, 'List of N2 x O2 dot products:'
        DO ITRAN = 1, NOTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &     IOTRAN(ITRAN),(IODOTS(I,ITRAN),I=1,MXO)
        END DO
        WRITE (LUPRI,*)
      END IF


      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCEXLR_SETUP                         *
*---------------------------------------------------------------------*
