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
c /* deck CC_OPA */
*=====================================================================*
       SUBROUTINE CC_OPA(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of first-order transition properties
*             (transition moments and oscillator strengths)
*             for the Coupled Cluster models
*
*                        CCS, CC2, CCSD, CC3
*
*             and partially for SCF and CIS
*
*     Written by Christof Haettig winter 2002/2003.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>
#include <ccopainf.h>
#include <ccsdinp.h>
#include <dummy.h>

* local parameters:
      CHARACTER*(16) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_OPA> ')

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0, TIM1, TIMF, TIMXE1, TIMXE2
      REAL SECOND
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0, TIM1, TIMF, TIMXE1, TIMXE2
      DOUBLE PRECISION SECOND
#endif

      LOGICAL LADD
      INTEGER NBOPA, MXTRAN, MXVEC, NFTRAN, NXE1TRAN, NXE2TRAN,
     &        KRESULT, KFTRAN, KFDOTS, KFCONS, KEND0, LEND0,
     &        KXE1TRAN, KX1DOTS, KX1CONS, KE1DOTS, KE1CONS,
     &        KXE2TRAN, KX2DOTS, KX2CONS,
     &        IOPT, IORDER


* external functions:
      INTEGER IR1TAMP
      INTEGER IL1ZETA

*---------------------------------------------------------------------*
* print header for second-order property section:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '******************************',
     & '*                                   ',
     &                               '                             *',
     & '*<<<<<<    OUTPUT FROM COUPLED CLUST',
     &                               'ER LINEAR RESPONSE    >>>>>>>*',
     & '*                                   ',
     &                               '                             *',
     & '*<<<<<<  CALCULATION OF ONE-PHOTON A',
     &                               'ABSORPTION STRENGHTS  >>>>>>>*',
     & '*                                   ',
     &                               '                             *',
     & '************************************',
     &                               '******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CC_OPA called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_OPA Workspace:',LWORK
  
      TIM0  = SECOND()

*---------------------------------------------------------------------*
* allocate & initialize work space for property contributions:
*---------------------------------------------------------------------*
      ! maximum number of transition moments to compute
      NBOPA   = 2 * NLRSOP * NXLRSST

      ! maximum number of transformations or vector calculations
      MXTRAN   = MAX(NXLRSST,NLRSOP)

      ! maximum number of vectors to dot on 
      MXVEC    = MAX(NXLRSST,NLRSOP)


      KRESULT  = 1
      KEND0    = KRESULT  + NBOPA
               
      KFTRAN   = KEND0
      KFDOTS   = KFTRAN   + MXTRAN * MXDIM_FTRAN
      KFCONS   = KFDOTS   + MXVEC  * MXTRAN
      KEND0    = KFCONS   + MXVEC  * MXTRAN

      KXE1TRAN = KEND0
      KX1DOTS  = KXE1TRAN + MXTRAN * MXDIM_XEVEC
      KX1CONS  = KX1DOTS  + MXVEC  * MXTRAN
      KE1DOTS  = KX1CONS  + MXVEC  * MXTRAN
      KE1CONS  = KE1DOTS  + MXVEC  * MXTRAN
      KEND0    = KE1CONS  + MXVEC  * MXTRAN

      KXE2TRAN = KEND0
      KX2DOTS  = KXE2TRAN + MXTRAN * MXDIM_XEVEC
      KX2CONS  = KX2DOTS  + MXVEC  * MXTRAN
      KEND0    = KX2CONS  + MXVEC  * MXTRAN

      LEND0 = LWORK - KEND0
      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_OPA. (1)')
      END IF

      CALL DZERO(WORK(KRESULT),NBOPA)

*---------------------------------------------------------------------*
* set up lists for F transformations, ETA{O} and Xi{O} vectors:
*---------------------------------------------------------------------*
      LADD = .FALSE.

      CALL CCOPA_SETUP(MXTRAN, MXVEC,
     &                 WORK(KFTRAN),WORK(KFDOTS),WORK(KFCONS),NFTRAN,
     &                 WORK(KXE1TRAN),WORK(KX1DOTS),WORK(KX1CONS),
     &                        WORK(KE1DOTS),WORK(KE1CONS),NXE1TRAN,
     &                 WORK(KXE2TRAN),WORK(KX2DOTS),WORK(KX2CONS),
     &                        NXE2TRAN,
     &                 WORK(KRESULT),NBOPA,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KFCONS),MXVEC*NFTRAN)

      IOPT = 5
      CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'L0 ','RE ',IOPT,'R1 ',
     &                WORK(KFDOTS),WORK(KFCONS),MXVEC,
     &                WORK(KEND0), LEND0)

      TIMF = SECOND() - TIM1

      IF (NFTRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NFTRAN,' F matrix transformations:',TIMF
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate ETA{O} x RE and LE x Xksi{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KX1CONS),MXVEC*NXE1TRAN)
      CALL DZERO(WORK(KE1CONS),MXVEC*NXE1TRAN)

      IOPT   = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KXE1TRAN), NXE1TRAN, IOPT, IORDER, 'L0 ',
     &               'LE ',WORK(KX1DOTS),WORK(KX1CONS),
     &               'RE ',WORK(KE1DOTS),WORK(KE1CONS),
     &               MXVEC, WORK(KEND0), LEND0 )

      TIMXE1 = SECOND() - TIM1
      IF (NXE1TRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NXE1TRAN,' O1/X1 vector calculation:',TIMXE1
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate M1 x Xksi{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KX2CONS),MXVEC*NXE2TRAN)

      IOPT   = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KXE2TRAN), NXE2TRAN, IOPT, IORDER, 'L0 ',
     &               'M1 ',WORK(KX2DOTS),WORK(KX2CONS),
     &               '---',IDUMMY,DUMMY,
     &               MXVEC, WORK(KEND0), LEND0 )

      TIMXE2 = SECOND() - TIM1
      IF (NXE2TRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NXE2TRAN,' O1/X1 vector calculation:',TIMXE2
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* collect contributions and sum them up to the final results:
*---------------------------------------------------------------------*
      LADD = .TRUE.

      CALL CCOPA_SETUP(MXTRAN, MXVEC,
     &                 WORK(KFTRAN),WORK(KFDOTS),WORK(KFCONS),NFTRAN,
     &                 WORK(KXE1TRAN),WORK(KX1DOTS),WORK(KX1CONS),
     &                        WORK(KE1DOTS),WORK(KE1CONS),NXE1TRAN,
     &                 WORK(KXE2TRAN),WORK(KX2DOTS),WORK(KX2CONS),
     &                        NXE2TRAN,
     &                 WORK(KRESULT),NBOPA,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBOPA,' quadratic response func.:', SECOND() - TIM0

*---------------------------------------------------------------------*
* print one-photon absorption properties and return:
*---------------------------------------------------------------------*
      CALL  CCOPAPRT(WORK(KRESULT))

      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_OPA                               *
*=====================================================================*
c /* deck ccopaprt */
*=====================================================================*
      SUBROUTINE CCOPAPRT(TRANMOM)
*---------------------------------------------------------------------*
*
*    Purpose: print transition momens and different one-photon 
*             absorption properties calculated from these
*
*    Written by Christof Haettig in winter 2002/2003.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <pgroup.h>
#include <ccroper.h>
#include <ccopainf.h>
#include <ccsdinp.h>
#include <ccexci.h>
#include <codata.h>


      CHARACTER*5  BLANKS
      CHARACTER*8  LABEL
      CHARACTER*10 MODEL
      CHARACTER*80 STRING


#if defined (SYS_CRAY)
      REAL TRANMOM(2,NLRSOP,NXLRSST)
      REAL DIPLEN(3,2), DIPVEL(3,2), ANGMOM(3,2)
      REAL ZERO, HALF, ONE, TWO, THREE, EIGV, SIGN
      REAL TMF0, TM0F, STRN
      REAL OSCILEN, OSCIVEL, OSCIMIX, ROTALEN, ROTAVEL
#else
      DOUBLE PRECISION TRANMOM(2,NLRSOP,NXLRSST)
      DOUBLE PRECISION DIPLEN(3,2), DIPVEL(3,2), ANGMOM(3,2)
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, THREE, EIGV, SIGN
      DOUBLE PRECISION TMF0, TM0F, STRN
      DOUBLE PRECISION OSCILEN, OSCIVEL, OSCIMIX, ROTALEN, ROTAVEL
#endif
      PARAMETER (ZERO=0.0d0, HALF=0.5d0, ONE=1.0d0, THREE=3.0d0)
      PARAMETER (TWO=2.0d0)

      LOGICAL LDIPLEN,LDIPL(3), LDIPVEL, LDIPV(3), LANGMOM, LANGM(3)
      INTEGER IRSD, ISTATE, ISYME, IMULE, ISTSY, IOPER, ISYMO, ISAMO,
     &        IDX

#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

*---------------------------------------------------------------------*
* print header for second-order properties: 
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING = ' RESULTS FOR ONE-PHOTON ABSORPTION STRENGTHS '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
         MODEL = 'CCS'
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
         MODEL = 'CC2'
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:45)//BLANKS )
         MODEL = 'CC3'
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
         MODEL = 'CCSD'
      ELSE
         CALL QUIT('CCOPAPRT called for an unknown Coupled '//
     &             'Cluster model.')
      END IF

      DO IRSD = 1, NXLRSST
       ISTATE = ILRSST(IRSD)
       ISYME  = ISYEXC(ISTATE)
       IMULE  = 1
       ISTSY  = ISTATE - ISYOFE(ISYME)
       EIGV   = EIGVAL(ISTATE)

       WRITE(LUPRI,'(//5X,A)') 'Transition:'
       WRITE(LUPRI,'(5X,A,I5,2X,"^",I1,A3,X,I1)')
     &   'number, multiplicity, symmetry :',ISTATE,IMULE,REP(ISYME-1)
       WRITE(LUPRI,'(5X,A,F15.10,A,F10.5,A,F10.1,A)')
     &   'frequency :',EIGV, ' a.u.  ',
     &          EIGV*XTEV,' e.V.  ',EIGV*XTKAYS,' cm^-1'

       WRITE(LUPRI,'(3(/5x,a,a))')
     &    '+-----------+-----------------+-----------------+',
     &                                        '---------------------+',
     &    '| operator  |   left moment   |  right moment   |',
     &                                        ' transition strength |',
     &    '+-----------+-----------------+-----------------+',
     &                                        '---------------------+'

       DO IDX = 1, 3
         LDIPL(IDX) = .FALSE.
         LDIPV(IDX) = .FALSE.
         LANGM(IDX) = .FALSE.
       END DO 

       DO IOPER = 1, NLRSOP
         ISYMO = ISYOPR(ILRSOP(IOPER))
         ISAMO = ISYMAT(ILRSOP(IOPER))
         LABEL = LBLOPR(ILRSOP(IOPER))

         SIGN  = DBLE(ISAMO)
         IF (ISAMO.EQ.0) SIGN = +ONE

         TMF0  = TRANMOM(1,IOPER,IRSD)
         TM0F  = TRANMOM(2,IOPER,IRSD)
         STRN  = TM0F * TMF0

         WRITE(LUPRI,
     &     '(5x,"| ",a8,2x,"|",2(x,f15.8,x,"|"),2x,f15.8,4x,"|")')
     &      LABEL,TMF0,TM0F,STRN

         
         ! collect transition moments for special operators

         ! -> here some dipole vectors:
         IF (LABEL(1:1).EQ.'X') IDX = 1
         IF (LABEL(1:1).EQ.'Y') IDX = 2
         IF (LABEL(1:1).EQ.'Z') IDX = 3

         IF      (LABEL(2:7).EQ.'DIPLEN') THEN
            DIPLEN(IDX,1) = TM0F
            DIPLEN(IDX,2) = TMF0
            LDIPL( IDX)   = .TRUE.
         ELSE IF (LABEL(2:7).EQ.'DIPVEL') THEN
            DIPVEL(IDX,1) = TM0F
            DIPVEL(IDX,2) = TMF0
            LDIPV( IDX)   = .TRUE.
         ELSE IF (LABEL(2:7).EQ.'ANGMOM') THEN
            ANGMOM(IDX,1) = TM0F
            ANGMOM(IDX,2) = TMF0
            LANGM( IDX)   = .TRUE.
         END IF

       END DO

       WRITE(LUPRI,'(5x,a,a)')
     &    '+-----------+-----------------+-----------------+',
     &    '---------------------+'

       LDIPLEN = LDIPL(1) .AND. LDIPL(2) .AND. LDIPL(3)
       LDIPVEL = LDIPV(1) .AND. LDIPV(2) .AND. LDIPV(3)
       LANGMOM = LANGM(1) .AND. LANGM(2) .AND. LANGM(3)

       ! dipole oscillator strengths
       IF (LDIPLEN) THEN
         OSCILEN=(TWO*EIGV/THREE)*DDOT(3,DIPLEN(1,1),1,DIPLEN(1,2),1)
         WRITE(LUPRI, '(6x,a,f15.8)')
     &        ' oscillator strength (length gauge)   : ',OSCILEN
       END IF
       IF (LDIPVEL) THEN
         OSCIVEL=(-TWO/(THREE*EIGV))*DDOT(3,DIPVEL(1,1),1,DIPVEL(1,2),1)
         WRITE(LUPRI, '(6x,a,f15.8)')
     &        ' oscillator strength (velocity gauge) : ',OSCIVEL
       END IF
       IF (LDIPLEN .AND. LDIPVEL) THEN
          OSCIMIX = ( TWO / THREE )
     &                 * HALF *( DDOT(3,DIPVEL(1,1),1,DIPLEN(1,2),1)-
     &                           DDOT(3,DIPLEN(1,1),1,DIPVEL(1,2),1)  )
         WRITE(LUPRI, '(6x,a,f15.8)')
     &         ' oscillator strength (mixed gauge)    : ',OSCIMIX
       END IF

       ! dipole rotator strengths
       IF (LDIPLEN .AND. LANGMOM) THEN
          ROTALEN = ( ONE / TWO )
     &       * HALF*( DDOT(3,DIPLEN(1,1),1,ANGMOM(1,2),1) -
     &                DDOT(3,ANGMOM(1,1),1,DIPLEN(1,2),1) )
         WRITE(LUPRI, '(6x,a,f15.8," a.u.")')
     &         ' rotator strength (length gauge)      : ',ROTALEN
       END IF
       IF (LDIPVEL .AND. LANGMOM) THEN
          ROTAVEL = ( -ONE / (TWO*EIGV) )
     &       * HALF*( DDOT(3,DIPVEL(1,1),1,ANGMOM(1,2),1) +
     &                DDOT(3,ANGMOM(1,1),1,DIPVEL(1,2),1) )
         WRITE(LUPRI, '(6x,a,f15.8," a.u.")')
     &         ' rotator strength (velocity gauge)    : ',ROTAVEL
       END IF  

      END DO ! IRSD

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCOPAPRT                             *
*---------------------------------------------------------------------*
c /* deck ccopa_setup */
*=====================================================================*
      SUBROUTINE CCOPA_SETUP(MXTRAN,  MXVEC,
     &                       IFTRAN,  IFDOTS,  FCONS,  NFTRAN,
     &                       IXE1TRAN,IX1DOTS, X1CONS,  
     &                                IE1DOTS, E1CONS, NXE1TRAN,
     &                       IXE2TRAN,IX2DOTS, X2CONS, NXE2TRAN,
     &                       RESULT,  MXOPA,   LADD,   WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC first-order transition moments
*         - list of F matrix transformations with eigenvectors
*         - list of XKSI and ETA vector contractions with eigenvectors
*         - list of XKSI vector contractions with Mbar multipliers
*
*     Written by Christof Haettig, Dec 2002 based on CCLR_SETUP 
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>
#include <ccopainf.h>
#include <ccroper.h>
#include <ccexci.h>
#include <ccsdinp.h>

* local parameters:
      CHARACTER*(21) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCOPA_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD
      INTEGER MXVEC, MXTRAN, MXOPA

      INTEGER IFTRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)
      INTEGER IXE1TRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IX1DOTS(MXVEC,MXTRAN), IE1DOTS(MXVEC,MXTRAN)
      INTEGER IXE2TRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IX2DOTS(MXVEC,MXTRAN)

      INTEGER NFTRAN, NXE1TRAN, NXE2TRAN, LWORK

#if defined (SYS_CRAY)
      REAL RESULT(MXOPA)
      REAL FCONS(MXVEC,MXTRAN)
      REAL X1CONS(MXVEC,MXTRAN), E1CONS(MXVEC,MXTRAN)
      REAL X2CONS(MXVEC,MXTRAN)
      REAL WORK(LWORK)
      REAL ZERO, SIGN, EIGV
      REAL WETARE, WXILE, WXIM1, WF
#else
      DOUBLE PRECISION RESULT(MXOPA)
      DOUBLE PRECISION FCONS(MXVEC,MXTRAN)
      DOUBLE PRECISION X1CONS(MXVEC,MXTRAN), E1CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION X2CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, SIGN, EIGV
      DOUBLE PRECISION WETARE, WXILE, WXIM1, WF
#endif
      PARAMETER (ZERO = 0.0D0)

      CHARACTER LABEL*(8)
      LOGICAL LORX, LPDBS
      INTEGER ITRAN, I, IRSD, ISTATE, ISYME, ISTSY, IOP, IOPER, ISYMO,
     &        IKAP, MXE1VEC, MXE2VEC, IM1VEC, IR1VEC, MFVEC, ITMF0,
     &        ITM0F, IVEC, NBOPA, IDUM

* external functions:
      INTEGER IR1TAMP
      INTEGER ILRMAMP

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      DO ITRAN = 1, MXTRAN
       IXE1TRAN(1,ITRAN)  = 0
       IXE1TRAN(2,ITRAN)  = 0
       IXE1TRAN(3,ITRAN)  = -1
       IXE1TRAN(4,ITRAN)  = -1
       IXE1TRAN(5,ITRAN)  = 0

       IXE2TRAN(1,ITRAN)  = 0
       IXE2TRAN(2,ITRAN)  = 0
       IXE2TRAN(3,ITRAN)  = -1
       IXE2TRAN(4,ITRAN)  = -1
       IXE2TRAN(5,ITRAN)  = 0

       DO I = 1, 3
        IFTRAN(I,ITRAN)  = 0
       END DO

       DO IVEC  = 1, MXVEC
        IFDOTS(IVEC,ITRAN)   = 0
        IX1DOTS(IVEC,ITRAN)  = 0
        IE1DOTS(IVEC,ITRAN)  = 0
        IX2DOTS(IVEC,ITRAN)  = 0
       END DO
      END DO

      NFTRAN   = 0
      NXE1TRAN = 0
      NXE2TRAN = 0

      NBOPA   = 0
      MFVEC   = 0
      MXE2VEC = 0
      MXE1VEC = 0
 

*---------------------------------------------------------------------*
* start loop over all requested transition moments:
*---------------------------------------------------------------------*
      DO IRSD  = 1, NXLRSST
       ISTATE = ILRSST(IRSD)
       ISYME  = ISYEXC(ISTATE)
       ISTSY  = ISTATE - ISYOFE(ISYME)
       EIGV   = EIGVAL(ISTATE)

       DO IOP = 1, NLRSOP
        IOPER = ILRSOP(IOP)
        LORX  = .FALSE.
        ISYMO = ISYOPR(IOPER)
        LABEL = LBLOPR(IOPER)
        LPDBS = LPDBSOP(IOPER)
        IKAP  = 0

        IF (LPDBS) CALL QUIT('perturbation-dependent basis sets not '//
     &              'implemented in CCOPA_SETUP.')


        IF (ISYMO.EQ.ISYME) THEN 

          NBOPA = NBOPA + 1

          IF (NBOPA.GT.MXOPA) THEN
             CALL QUIT('NBOPA out of range in CCOPA_SETUP.')
          END IF

*---------------------------------------------------------------------*
*         in all cases we need Eta{X} x RE 
*---------------------------------------------------------------------*
          CALL CC_SETXE('Eta',IXE1TRAN,IE1DOTS,MXTRAN,MXVEC,
     &                  0,IOPER,IKAP,0,0,0,ISTATE,ITRAN,IVEC)
          NXE1TRAN = MAX(NXE1TRAN,ITRAN)
          MXE1VEC  = MAX(MXE1VEC, IVEC)
          WETARE   = E1CONS(IVEC,ITRAN)

*---------------------------------------------------------------------*
*         in all cases we need LE x Xksi{X} 
*---------------------------------------------------------------------*
          CALL CC_SETXE('Xi ',IXE1TRAN,IX1DOTS,MXTRAN,MXVEC,
     &                  0,IOPER,IKAP,0,0,0,ISTATE,ITRAN,IVEC)
          NXE1TRAN = MAX(NXE1TRAN,ITRAN)
          MXE1VEC  = MAX(MXE1VEC, IVEC)
          WXILE    = X1CONS(IVEC,ITRAN)

*---------------------------------------------------------------------*
*         add M * Xksi{X} or F * RE * R1, depending on LRS2N1
*---------------------------------------------------------------------*
          WXIM1  = ZERO
          WF     = ZERO

          IF (.NOT.CIS) THEN
            IF (LRS2N1) THEN
              IM1VEC = ILRMAMP(ISTATE,EIGV,ISYME)
              CALL CC_SETXE('Xi ',IXE2TRAN,IX2DOTS,MXTRAN,MXVEC,
     &                      0,IOPER,IKAP,0,0,0,IM1VEC,ITRAN,IVEC)
              NXE2TRAN = MAX(NXE2TRAN,ITRAN)
              MXE2VEC  = MAX(MXE2VEC, IVEC)
              WXIM1    = X2CONS(IVEC,ITRAN)
            ELSE
              IR1VEC = IR1TAMP(LABEL,LORX,-EIGV,IDUM)
              CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                       0,ISTATE,IR1VEC,ITRAN,IVEC)
              NFTRAN = MAX(NFTRAN,ITRAN)
              MFVEC  = MAX(MFVEC, IVEC)
              WF     = FCONS(IVEC,ITRAN)
            END IF
          END IF

*---------------------------------------------------------------------*
*          add contributions together:
*---------------------------------------------------------------------*
           IF (LADD) THEN

              ITMF0 = (NLRSOP*(IRSD-1) + IOP-1)*2 + 1
              ITM0F = (NLRSOP*(IRSD-1) + IOP-1)*2 + 2

              RESULT(ITMF0) = WXILE
              RESULT(ITM0F) = WETARE + WXIM1 + WF

              IF (LOCDBG) THEN
                 WRITE (LUPRI,*) 'ISTATE, EIGV:',ISTATE,EIGV
                 WRITE (LUPRI,*) 'OPERATOR:',LABEL
                 WRITE (LUPRI,*) 'IDX = ',NLRSOP*(IRSD-1) + IOP
                 WRITE (LUPRI,*) 'ITMF0,ITM0F:',ITMF0,ITM0F
                 WRITE (LUPRI,*) 'L^f x Xksi{X}:',WXILE
                 WRITE (LUPRI,*) '<L^f|X|CC>:',RESULT(ITMF0)
                 WRITE (LUPRI,*) 'Eta{X} x R^f :',WETARE
                 WRITE (LUPRI,*) 'M^f x Xksi{X}:',WXIM1
                 WRITE (LUPRI,*) 'F x R^f x M^f:',WF
                 WRITE (LUPRI,*) '<Lambda|[X,R^f]|CC>:',RESULT(ITM0F)
              END IF

           END IF

*---------------------------------------------------------------------*
*       end loop over transition moments
*---------------------------------------------------------------------*

        END IF
       END DO
      END DO

      IF      (MFVEC.GT.MXVEC) THEN
         CALL QUIT('MFVEC has been out of bounds in CCOPA_SETUP.')
      ELSE IF (MXE1VEC.GT.MXVEC) THEN
         CALL QUIT('MXE1VEC has been out of bounds in CCOPA_SETUP.')
      ELSE IF (MXE2VEC.GT.MXVEC) THEN
         CALL QUIT('MXE2VEC has been out of bounds in CCOPA_SETUP.')
      ELSE IF (NFTRAN.GT.MXTRAN) THEN
         CALL QUIT('NFTRAN has been out of bounds in CCOPA_SETUP.')
      ELSE IF (NXE1TRAN.GT.MXTRAN) THEN
         CALL QUIT('NXE1TRAN has been out of bounds in CCOPA_SETUP.')
      ELSE IF (NXE2TRAN.GT.MXTRAN) THEN
         CALL QUIT('NXE2TRAN has been out of bounds in CCOPA_SETUP.')
      END IF

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      IF ((.NOT.LADD) .OR. LOCDBG) THEN
       WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NBOPA,
     &      ' transition moments'
       WRITE(LUPRI,'((8X,A,I3,A))') 
     & ' - ',NFTRAN,  ' F matrix transformations with R1 vectors',
     & ' - ',NXE1TRAN,' ETA and XKSI vector calculations ',
     & ' - ',NXE2TRAN,' extra XKSI vector calculations '
       WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'
      END IF

      IF (LOCDBG) THEN

         ! F matrix transformations:
         WRITE(LUPRI,*)'List of F matrix transformations:'
         DO ITRAN = 1, NFTRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IFTRAN(I,ITRAN),I=1,2),(IFDOTS(I,ITRAN),I=1,MFVEC)
         END DO
         WRITE(LUPRI,*)

         ! Xi{O} and ETA{O} vector calculations:
         WRITE(LUPRI,*) 'List of Xi{O} and ETA{O} vector calculations:'
         DO ITRAN = 1, NXE1TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IXE1TRAN(I,ITRAN),I=1,5),(IX1DOTS(I,ITRAN),I=1,MXE1VEC)
           WRITE(LUPRI,'(A,25X,5X,(25I3,20X))') MSGDBG,
     &                               (IE1DOTS(I,ITRAN),I=1,MXE1VEC)
         END DO
         WRITE(LUPRI,*)

         ! extra Xi{O} vector calculations:
         WRITE(LUPRI,*) 'List of extra Xi{O} vector calculations:'
         DO ITRAN = 1, NXE2TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IXE2TRAN(I,ITRAN),I=1,5),(IX2DOTS(I,ITRAN),I=1,MXE2VEC)
         END DO
         WRITE(LUPRI,*)

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCOPA_SETUP                          *
*---------------------------------------------------------------------*
