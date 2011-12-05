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
c /* deck cc_cauchy */
*=====================================================================*
       SUBROUTINE CC_CAUCHY(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: calculation of linear response Cauchy moments
*
*             implemented models: CCS, CC2, CCSD
*
*     Written by Christof Haettig december 1996.
*     modified for ASYMSD option in october 1997, Christof Haettig
*     odd dispersion coefficients in februar 1999, Christof Haettig
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
#include <ccroper.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cclists.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMREF ! symmetry of reference state, assumed total sym.
      PARAMETER (ISYMREF = 1)
    

#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION ZERO, ONE, TWO
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0D0)

      CHARACTER*80 STRING
      CHARACTER*10 MODEL
      CHARACTER*8  LABELA, LABELB, BLANKS
      INTEGER LWORK
      INTEGER ISYM, ITRAN, IVEC, IDUM, IDXSMOM, IOPER, IOPT, IDISP
      INTEGER ISYMA, IRCVECA, ILCVECA, IOPERA, ISAMA, ISAPROP
      INTEGER ISYMB, IRCVECB, ILCVECB, IOPERB, ISAMB, ISACAU
      INTEGER MXTRAN, MXVEC, MXDOTS, NBCAUCHY, LEND0, KEND0, KRESULT
      INTEGER MXFTRAN, KFCTRAN, KFCDOTS, KFCCONS, NFCTRAN
      INTEGER MXXETRAN,KXETRAN,KECDOTS,KECCONS,KXCDOTS,KXCCONS,NXETRAN
      INTEGER KTRK, ISTRT, IORDER, IDISP0

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TRKBAS, SMOM, SIGN, ECONAB, ECONBA, FCON, XCONAB
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TRKBAS, SMOM, SIGN, ECONAB, ECONBA, FCON, XCONAB
#endif

* external functions:
      INTEGER ILRCAMP
      INTEGER ILC1AMP

*---------------------------------------------------------------------*
* print header for Cauchy moment section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(5(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*     COUPLED CLUSTER LINEAR RESPONS',
     &                               'E DISPERSION COEFFICIENTS     *',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************'

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CC_CAUCHY called for unknown Coupled Cluster.')
      END IF

      IF (ASYMSD .AND. CC3)
     &  CALL QUIT('Assymmetric formula not (yet) implemented  for CC3')

* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_CAUCHY Workspace:',LWORK
  
*---------------------------------------------------------------------*
* add C(0),...,C(-4) Cauchy vectors to Cauchy vector list:
*---------------------------------------------------------------------*

* open lists:
      LRC1OPN = .TRUE.
      LLC1OPN = .TRUE.

* add vectors to list:
      DO IOPER = 1, NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        IF (ISYMA.EQ.ISYMB) THEN
          IVEC = ILRCAMP(LABELA, 0,ISYMA)
          IVEC = ILRCAMP(LABELB, 0,ISYMB)
          IF (.NOT.CC3) THEN
          IVEC = ILRCAMP(LABELA,-1,ISYMA)
          IVEC = ILRCAMP(LABELB,-1,ISYMB)
          IVEC = ILRCAMP(LABELA,-2,ISYMA)
          IVEC = ILRCAMP(LABELB,-2,ISYMB)
          IVEC = ILRCAMP(LABELA,-3,ISYMA)
          IVEC = ILRCAMP(LABELB,-3,ISYMB)
          IVEC = ILRCAMP(LABELA,-4,ISYMA)
          IVEC = ILRCAMP(LABELB,-4,ISYMB)
          END IF
          IF (ASYMSD) THEN
            IVEC = ILC1AMP(LABELA, 0,ISYMA)
            IVEC = ILC1AMP(LABELB, 0,ISYMB)
            IF (.NOT.CC3) THEN
            IVEC = ILC1AMP(LABELA,-1,ISYMA)
            IVEC = ILC1AMP(LABELB,-1,ISYMB)
            IVEC = ILC1AMP(LABELA,-2,ISYMA)
            IVEC = ILC1AMP(LABELB,-2,ISYMB)
            IVEC = ILC1AMP(LABELA,-3,ISYMA)
            IVEC = ILC1AMP(LABELB,-3,ISYMB)
            IVEC = ILC1AMP(LABELA,-4,ISYMA)
            IVEC = ILC1AMP(LABELB,-4,ISYMB)
            END IF
          END IF
        ENDIF
      END DO

* close lists:
      LRC1OPN = .FALSE.
      LLC1OPN = .FALSE.

* calculate cauchy vectors for negative orders:
      IF (.NOT.CC3) CALL CCCM_VEC(WORK,LWORK)

*---------------------------------------------------------------------*
* allocate & initialize work space for Cauchy moments
*---------------------------------------------------------------------*

      NBCAUCHY = NLROP * (NLRDISP + 5)

      MXTRAN  = NLRCLBL 
      MXVEC   = NLRCLBL
      MXDOTS  = MXVEC * MXTRAN

      MXFTRAN  = MXDIM_FTRAN * MXTRAN
      MXXETRAN = MXDIM_XEVEC * MXTRAN

      KTRK    = 1
      KRESULT = KTRK    + NLROP
      KFCTRAN = KRESULT + NBCAUCHY
      KFCDOTS = KFCTRAN + MXFTRAN
      KFCCONS = KFCDOTS + MXDOTS
      KXETRAN = KFCCONS + MXDOTS
      KECDOTS = KXETRAN + MXXETRAN
      KECCONS = KECDOTS + MXDOTS
      KXCDOTS = KECCONS + MXDOTS 
      KXCCONS = KXCDOTS + MXDOTS
      KEND0   = KXCCONS + MXDOTS

      LEND0   = LWORK - KEND0

*---------------------------------------------------------------------*
* set up lists F transformations and ETA vectors:
*---------------------------------------------------------------------*
      IDISP0 = -4
      IF (CC3) IDISP0 = 0

      CALL CCCM_SETUP(MXTRAN, MXVEC, IDISP0,
     &                WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,
     &                WORK(KXETRAN),WORK(KECDOTS),WORK(KXCDOTS),NXETRAN)

*---------------------------------------------------------------------*
* calculate F matrix contributions of the type F x RC x RC:
*---------------------------------------------------------------------*
      IOPT = 5 
      CALL CC_FMATRIX(WORK(KFCTRAN),NFCTRAN,'L0 ','RC ',IOPT,'RC ',
     &                WORK(KFCDOTS),WORK(KFCCONS),MXVEC,
     &                WORK(KEND0),LEND0) 

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions of the type ETA{O} x RC
* and the LC x XKSI{O} vector dot products:
*---------------------------------------------------------------------*

      CALL DZERO(WORK(KXCCONS),MXVEC*NXETRAN)
      CALL DZERO(WORK(KECCONS),MXVEC*NXETRAN)

      IOPT = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KXETRAN), NXETRAN, IOPT, IORDER, 'L0 ',
     &               'LC ',WORK(KXCDOTS),WORK(KXCCONS),
     &               'RC ',WORK(KECDOTS),WORK(KECCONS),
     &               .FALSE.,MXVEC, WORK(KEND0), LEND0 )

*=====================================================================*
* collect all contributions and add them up to Cauchy moments:
*=====================================================================*
      CALL DZERO(WORK(KRESULT),NBCAUCHY)

      DO IOPER = 1, NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        ISAMA  = ISYMAT(IOPERA)
        ISAMB  = ISYMAT(IOPERB)

        ! flag: +1: real, -1: imaginary, 0: unknown
        ISAPROP = ISAMA * ISAMB

      IF ( ISYMA .EQ. ISYMB ) THEN

*---------------------------------------------------------------------*
* calculate the Cauchy moments S(+2):
*---------------------------------------------------------------------*
       IF (.NOT.CC3) THEN
       IF (ISAPROP.EQ.1 .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN

        SMOM = ZERO 

        ! Eta{B} x Jacobian x C^A(-4):
        IRCVECB = ILRCAMP(LABELB,-4,IDUM)
        CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KECDOTS),NXETRAN,MXVEC,
     &                0,IOPERA,0,0,0,0,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)

        ! Eta{A} x Jacobian x C^B(-4) (symmetric formula only):
        IF (.NOT.ASYMSD) THEN
          IRCVECA = ILRCAMP(LABELA,-4,IDUM)
          CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KECDOTS),NXETRAN,MXVEC,
     &                  0,IOPERB,0,0,0,0,IRCVECA,ITRAN,IVEC)
          SMOM = SMOM + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        END IF

        ! F matrix contributions (symmetric formula only):
        IF (.NOT.ASYMSD) THEN
          ! F x C^A(-1) x C^B(-3):
          IRCVECA = ILRCAMP(LABELA,-1,IDUM)
          IRCVECB = ILRCAMP(LABELB,-3,IDUM)
          CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                   0,IRCVECA,IRCVECB,ITRAN,IVEC)
          SMOM = SMOM + WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          ! F x C^A(-2) x C^B(-2):
          IRCVECA = ILRCAMP(LABELA,-2,IDUM)
          IRCVECB = ILRCAMP(LABELB,-2,IDUM)
          CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                   0,IRCVECA,IRCVECB,ITRAN,IVEC)
          SMOM = SMOM - WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          
          ! F x C^A(-3) x C^B(-1):
          IRCVECA = ILRCAMP(LABELA,-3,IDUM)
          IRCVECB = ILRCAMP(LABELB,-1,IDUM)
          CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                   0,IRCVECA,IRCVECB,ITRAN,IVEC)
          SMOM = SMOM + WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        END IF

        ! Xksi{A} x Cbar^B(-4) (assymmetric formula only):
        IF (ASYMSD) THEN
          ILCVECB = ILC1AMP(LABELB,-4,IDUM)
          CALL CC_SETXE('Xi ',WORK(KXETRAN),WORK(KXCDOTS),NXETRAN,MXVEC,
     &                  0,IOPERA,0,0,0,0,ILCVECB,ITRAN,IVEC)
          SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        END IF

        ! store the result:
        IDXSMOM = (IOPER-1)*(NLRDISP+5) + 1
        WORK(KRESULT-1+IDXSMOM) = SMOM

        IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,2A,I2,5X,A,I3,A,E20.10)')
     &     LABELA,LABELB,IDXSMOM,'S(',2,')=', WORK(KRESULT-1+IDXSMOM)
        END IF


       END IF 
       END IF 

*---------------------------------------------------------------------*
* calculate zeroth-order Cauchy moments (TRK sum rules):
*---------------------------------------------------------------------*
       IF (.NOT.CC3) THEN
       IF (ISAPROP.EQ.1 .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN
        SMOM = ZERO 

        ! Eta{A} x C^B(-2):
        IRCVECB = ILRCAMP(LABELB,-2,IDUM)
        CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KECDOTS),NXETRAN,MXVEC,
     &                0,IOPERA,0,0,0,0,IRCVECB,ITRAN,IVEC)
        ECONAB = WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        SMOM = SMOM + ECONAB
        IF (LOCDBG) WRITE(LUPRI,*) 'ECONAB:',ECONAB


        ! Eta{B} x C^A(-2) (symmetric formula only):
        IF (.NOT.ASYMSD) THEN
          IRCVECA = ILRCAMP(LABELA,-2,IDUM)
          CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KECDOTS),NXETRAN,MXVEC,
     &                  0,IOPERB,0,0,0,0,IRCVECA,ITRAN,IVEC)
          ECONBA = WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          SMOM = SMOM + ECONBA
        END IF


        ! F x C^A(-1) x C^B(-1) (assymmetric formula only):
        IF (.NOT.ASYMSD) THEN
          IRCVECA = ILRCAMP(LABELA,-1,IDUM)
          IRCVECB = ILRCAMP(LABELB,-1,IDUM)
          CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                   0,IRCVECA,IRCVECB,ITRAN,IVEC)
          FCON = WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          SMOM = SMOM + FCON

        END IF

        ! Xksi{A} x Cbar^B(-2)
        IF (ASYMSD) THEN
          ILCVECB = ILC1AMP(LABELB,-2,IDUM)
          CALL CC_SETXE('Xi ',WORK(KXETRAN),WORK(KXCDOTS),NXETRAN,MXVEC,
     &                  0,IOPERA,0,0,0,0,ILCVECB,ITRAN,IVEC)
          XCONAB = WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          SMOM = SMOM + XCONAB
          IF (LOCDBG) WRITE(LUPRI,*) 'XCONAB:',XCONAB
        END IF

        ! store coupled cluster number:
        IDXSMOM = (IOPER-1)*(NLRDISP+5) + 3
        WORK(KRESULT-1+IDXSMOM) = SMOM

        ! for CCS calc. also pure basis set number ( = RPA result):
        IF (CCS) THEN
          IF (ASYMSD) THEN
            WORK(KTRK-1+IOPER) = TWO * XCONAB
          ELSE
            WORK(KTRK-1+IOPER) = ECONAB + ECONBA + TWO*FCON
          END IF
        END IF

        IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,2A,I2,5X,A,I3,A,E20.10)')
     &     LABELA,LABELB,IDXSMOM,'S(',0,')=', WORK(KRESULT-1+IDXSMOM)
        END IF

       END IF 
       END IF 

*---------------------------------------------------------------------*
* calculate Cauchy moments S(-2*K-2) for K>=0:
*---------------------------------------------------------------------*
      DO IDISP = 0, NLRDISP

       ! flag: +1: even, -1: odd
       ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1

       IF (ISAPROP.EQ.ISACAU .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN

        SMOM = ZERO

c       --------------
c       Eta{A} x RC^B: 
c       --------------
        IRCVECB = ILRCAMP(LABELB,IDISP,IDUM)

        CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KECDOTS),NXETRAN,MXVEC,
     &                0,IOPERA,0,0,0,0,IRCVECB,ITRAN,IVEC)
        ECONAB = WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        SMOM = SMOM + ECONAB
        IF (LOCDBG) WRITE(LUPRI,*) 'ECONAB:',ECONAB

c       --------------
c       Eta{B} x RC^A: 
c       --------------
        IF (.NOT.ASYMSD) THEN
          IRCVECA = ILRCAMP(LABELA,IDISP,IDUM)
          SIGN    = DBLE(ISACAU)

          CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KECDOTS),NXETRAN,MXVEC,
     &                  0,IOPERB,0,0,0,0,IRCVECA,ITRAN,IVEC)
          ECONBA = SIGN * WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          SMOM = SMOM + ECONBA
          IF (LOCDBG) WRITE(LUPRI,*) 'ECONBA:',ECONBA
        END IF


c       ----------------
c       F x RC^A x RC^B: 
c       ----------------
        IF (.NOT.ASYMSD) THEN
          SIGN = ONE
          DO N = 0, IDISP
            IRCVECA = ILRCAMP(LABELA,      N,IDUM)
            IRCVECB = ILRCAMP(LABELB,IDISP-N,IDUM)

            CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                     0,IRCVECA,IRCVECB,ITRAN,IVEC)
            FCON = SIGN * WORK(KFCCONS-1+(ITRAN-1)*MXVEC+IVEC)
            SMOM = SMOM + FCON
            IF (LOCDBG) WRITE(LUPRI,*) 'N,FCON:',N,FCON
            SIGN = -SIGN
          END DO
        END IF

c       ----------------
c       Xi{A} x LC^B: 
c       ----------------
        IF (ASYMSD) THEN
          ILCVECB = ILC1AMP(LABELB,IDISP,IDUM)
          CALL CC_SETXE('Xi ',WORK(KXETRAN),WORK(KXCDOTS),NXETRAN,MXVEC,
     &                  0,IOPERA,0,0,0,0,ILCVECB,ITRAN,IVEC)
          SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        END IF

        IDXSMOM = (IOPER-1)*(NLRDISP+5) + IDISP+5
        WORK(KRESULT-1+IDXSMOM) = SMOM

        IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,2A,I2,5X,A,I3,A,E20.10)') LABELA,LABELB,
     &     IDXSMOM,'S(',-IDISP-2,')=', WORK(KRESULT-1+IDXSMOM)
        END IF

       END IF

      END DO 

*---------------------------------------------------------------------*
* close loop over operator labels
*---------------------------------------------------------------------*
      END IF
      END DO 

*---------------------------------------------------------------------*
* print output & return:
*---------------------------------------------------------------------*

      BLANKS = '     '
      STRING = ' RESULTS FOR DISPERSION COEFFICIENTS'
      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CC_CAUCHY called for an unknown '//
     &             'Coupled Cluster model.')
      END IF

      WRITE (LUPRI,'(3X,2A,/3X,A,/,/3X,A,/3X,A,/,/3X,2A,/3X,A,/)')
     &    ' The dispersion coefficients D_AB(n) for n>0',
     &                        ' are defined through the expansion:',
     &    '                       alpha(w) = sum_n w^n D_AB(n) ',
     &    ' D_AB(-4) corresponds to the S_AB(+2) Cauchy moment,',
     &    ' D_AB(-2) corresponds to the S_AB( 0) Cauchy moment.',
     &    ' Note that the sign of dispersion coefficients is',
     &          ' here defined opposite to',
     &    ' the sign of the linear response functions!'

      WRITE(LUPRI,'(/1X,2(A," operator",3X),2X,A,8X,A,/,88("-"))')
     &    'A','B','n','D_AB'

      DO IOPER = 1, NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        ISAMA  = ISYMAT(IOPERA)
        ISAMB  = ISYMAT(IOPERB)

        ! flag: +1: real, -1: imaginary, 0: unknown
        ISAPROP = ISAMA * ISAMB

        ISTRT = +1
        IF (ISAPROP.EQ.1 .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) ISTRT=-4
        IF (CC3) ISTRT = MAX(ISTRT,0)

      IF (ISYMA.EQ.ISYMB) THEN

        DO IDISP = ISTRT, NLRDISP
          ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1
          IF (ISAPROP.EQ.ISACAU .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN
            IDXSMOM = (IOPER-1)*(NLRDISP+5) + IDISP+5
            SMOM = WORK(KRESULT-1+IDXSMOM)
            IF (IDISP.EQ.ISTRT) THEN
              WRITE(LUPRI,'(/3X,2(A8,4X),I3,5X,G20.12)')
     &           LABELA,LABELB,IDISP,-SMOM
            ELSE
              IF (IDISP.EQ.-2 .AND. CCS) THEN
               TRKBAS = WORK(KTRK-1+IOPER)
               WRITE(LUPRI,'(3X,2(8X,4X),I3,5X,G20.12,A,G20.12,A)') 
     &                IDISP,-SMOM, '(RPA RESULT:',TRKBAS,')'
              ELSE
               WRITE(LUPRI,'(3X,2(8X,4X),I3,5X,G20.12)') IDISP,-SMOM
              END IF
            END IF
          END IF
        END DO

      ELSE

        WRITE(LUPRI,'(/3X,2(A8,4X),3X,5X,A)')
     &     LABELA,LABELB,'  -.---      '

      END IF

      END DO

      WRITE(LUPRI,'(/,88("-"),//)')

      RETURN
      END

*---------------------------------------------------------------------*
*                END OF SUBROUTINE CC_CAUCHY                          *
*---------------------------------------------------------------------*
c /* deck cccm_setup */
*=====================================================================*
      SUBROUTINE CCCM_SETUP(MXTRAN, MXVEC, IDISP0,
     &                      IFCTRAN, IFCDOTS, NFCTRAN,
     &                      IXETRAN, IECDOTS, IXCDOTS, NXETRAN)
*---------------------------------------------------------------------*
*
*    Purpose: set up for CCLR Cauchy moment section
*         - list of F matrix transformations with Cauchy vectors
*         - list of ETA vector calculations 
*         - list of XKSI vector dot products 
*
*     Written by Christof Haettig, december 1996.
*     modified for ASYMSD option, C.H., october 1997.
*     odd dispersion coefficients, C.H., februar 1999.
*     restructured and updated for CC3, C.H., april 2004.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccroper.h>
#include <ccrc1rsp.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCCM_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MXVEC, MXTRAN

      INTEGER IFCTRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IFCDOTS(MXVEC,MXTRAN)
      INTEGER IXETRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IECDOTS(MXVEC,MXTRAN)
      INTEGER IXCDOTS(MXVEC,MXTRAN)

      INTEGER NFRTRAN, NFCTRAN, NXETRAN, NCAUMOM, MXEVEC, MFCVEC

      INTEGER ISYMA,  ISYMB, ITRAN, IVEC, ISYML, IDUM, I, N, IDISP
      INTEGER IRCVECA, IRCVECB, IOPERA, IOPERB, ISAMA, ISAMB, IDISP0
      INTEGER ILCVECA, ILCVECB, INUM, IOPER, ILROP, ISAPROP, ISACAU
      INTEGER IRHSA, IRHSB

      CHARACTER LABELA*(8), LABELB*(8)

* external functions:
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER IRHSR1

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      DO ITRAN = 1, MXTRAN
       DO I = 1, MXDIM_FTRAN
        IFCTRAN(I,ITRAN) = 0
       END DO
       DO I = 1, MXDIM_XEVEC
        IXETRAN(I,ITRAN) = 0
       END DO

       DO IVEC  = 1, MXVEC
        IFCDOTS(IVEC,ITRAN) = 0
        IECDOTS(IVEC,ITRAN) = 0
        IXCDOTS(IVEC,ITRAN) = 0
       END DO
      END DO

      NFCTRAN = 0
      NXETRAN = 0
      NCAUMOM = 0
      MFCVEC  = 0
      MXEVEC  = 0
 
*---------------------------------------------------------------------*
* start loop over all requested Cauchy moments 
*---------------------------------------------------------------------*
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCCM_SETUP: NLROP = ',NLROP
      ENDIF

      DO IOPER = 1, NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        ISAMA  = ISYMAT(IOPERA)
        ISAMB  = ISYMAT(IOPERB)

        ! flag: +1: real, -1 imaginary, 0 unknown
        ISAPROP = ISAMA * ISAMB

      IF (ISYMA.EQ.ISYMB) THEN

        NCAUMOM = NCAUMOM + 1

*---------------------------------------------------------------------*
* set up list of ETA vector calculations:
*---------------------------------------------------------------------*
        DO IDISP = IDISP0, NLRDISP

          ! flag: +1: even, -1: odd 
          ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1

          IF (ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN
     
            ! we need allways Eta{A} x RC^B
            IRCVECB = ILRCAMP(LABELB,IDISP,IDUM)
            CALL CC_SETXE('Eta',IXETRAN,IECDOTS,MXTRAN,MXVEC,
     &                    0,IOPERA,0,0,0,0,IRCVECB,ITRAN,IVEC)
            NXETRAN = MAX(NXETRAN,ITRAN)
            MXEVEC  = MAX(MXEVEC, IVEC)

            ! for symmetric formula also Eta{B} x RC^A is needed
            IF (.NOT.ASYMSD) THEN
              IRCVECA = ILRCAMP(LABELA,IDISP,IDUM)
              CALL CC_SETXE('Eta',IXETRAN,IECDOTS,MXTRAN,MXVEC,
     &                      0,IOPERB,0,0,0,0,IRCVECA,ITRAN,IVEC)
              NXETRAN = MAX(NXETRAN,ITRAN)
              MXEVEC  = MAX(MXEVEC, IVEC)
            END IF

          END IF

        END DO

*---------------------------------------------------------------------*
* set up list of F matrix transformations (symmetric formula only):
*---------------------------------------------------------------------*
        IF (.NOT.ASYMSD) THEN
         DO IDISP = IDISP0, -2
          ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1
          IF (ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN
            DO N = IDISP+1, -1
              IRCVECA = ILRCAMP(LABELA,      N,IDUM)
              IRCVECB = ILRCAMP(LABELB,IDISP-N,IDUM)

              CALL CCQR_SETF(IFCTRAN,IFCDOTS,MXTRAN,MXVEC,
     &                       0,MIN(IRCVECA,IRCVECB),
     &                         MAX(IRCVECA,IRCVECB),ITRAN,IVEC)
              NFCTRAN = MAX(NFCTRAN,ITRAN)
            END DO
          END IF
         END DO

         DO IDISP = MAX(IDISP0,0), NLRDISP
          ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1
          IF (ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN
            DO N = 0, IDISP
              IRCVECA = ILRCAMP(LABELA,      N,IDUM)
              IRCVECB = ILRCAMP(LABELB,IDISP-N,IDUM)

              CALL CCQR_SETF(IFCTRAN,IFCDOTS,MXTRAN,MXVEC,
     &                       0,MIN(IRCVECA,IRCVECB),
     &                         MAX(IRCVECA,IRCVECB),ITRAN,IVEC)
              NFCTRAN = MAX(NFCTRAN,ITRAN)
              MFCVEC  = MAX(MFCVEC, IVEC)
            END DO
          END IF
         END DO
        END IF

*---------------------------------------------------------------------*
* set up list of XKSI x LC dot products (asymmetric formula only):
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN
          DO IDISP = IDISP0, NLRDISP
            ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1
            IF (ISACAU.EQ.ISAPROP.OR.ISAPROP.EQ.0.OR.ALLLRDSPCF) THEN
              ILCVECB = ILC1AMP(LABELB,IDISP,IDUM)
              CALL CC_SETXE('Xi ',IXETRAN,IXCDOTS,MXTRAN,MXVEC,
     &                      0,IOPERA,0,0,0,0,ILCVECB,ITRAN,IVEC)
              NXETRAN = MAX(NXETRAN,ITRAN)
              MXEVEC  = MAX(MXEVEC, IVEC)
            END IF
          END DO
        END IF

*---------------------------------------------------------------------*
* end loop over all requested Cauchy moments
*---------------------------------------------------------------------*
      END IF
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      WRITE (LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NCAUMOM,
     &      ' Cauchy moments'
      WRITE (LUPRI,'((8X,A,I3,A))') 
     & ' - ',NFCTRAN,' F matrix transformations with Cauchy vectors',
     & ' - ',NXETRAN,' ETA and XKSI vector calculations '
      WRITE (LUPRI,'(3X,A,/,/)') 'will be performed.'

      IF (LOCDBG) THEN

        ! list of cauchy vectors:
        WRITE(LUPRI,'(/A)')' LIST OF REQUIRED RIGHT CAUCHY VECTORS:'
        DO I = 1, NLRCLBL
           WRITE(LUPRI,'(I5,3X,A8,I5,2X,1P,I5)') 
     &            I, LRCLBL(I), ISYLRC(I), ILRCAU(I)
        END DO

        ! F matrix transformations:
        WRITE (LUPRI,*)'List of F matrix transformations:'
        DO ITRAN = 1, NFCTRAN
          WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &     (IFCTRAN(I,ITRAN),I=1,2),(IFCDOTS(I,ITRAN),I=1,MFCVEC)
        END DO
        WRITE (LUPRI,*)

        ! Xi{O} and ETA{O} vector calculations:
        WRITE (LUPRI,*) 'List of Xi{O} and ETA{O} vector calculations:'
        DO ITRAN = 1, NXETRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IXETRAN(I,ITRAN),I=1,5),(IXCDOTS(I,ITRAN),I=1,MXEVEC)
           WRITE(LUPRI,'(A,25X,5X,(25I3,20X))') MSGDBG,
     &                               (IECDOTS(I,ITRAN),I=1,MXEVEC)
        END DO
        WRITE (LUPRI,*)

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCM_SETUP                           *
*---------------------------------------------------------------------*
c /* deck CCCM_VEC */
*=====================================================================*
      SUBROUTINE CCCM_VEC(WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: calculate negative-order Cauchy vectors by succesive
*             multiplication with the Jacoby matrix
*
*
*    Written by Christof Haettig March 1997
*    Left Cauchy vectors introduced, October 1997, Christof Haettig
*    Restructured to compute every vector on the list only once and
*    to introduce CC3, March 2004, Christof Haettig
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
#include <ccsdsym.h>
#include <ccroper.h>
#include <cclrinf.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <dummy.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCCM_VEC> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE
#endif
      PARAMETER (ZERO=0.0d0, ONE=1.0d0)

      CHARACTER*(10) MODEL, MODELW
      CHARACTER*(8)  LABELA, LABELB, LABEL
      CHARACTER*(3)  APROXR12
      INTEGER ILROP, IVECR, IVECR1, IVECL, IVECL1, KAXI1, KAXI2
      INTEGER ISYML, IOPT, IOPERA, IOPERB, ISYMA, ISYMB, KAXI0
      INTEGER IOPTW, IDLSTRC, IDLSTLC, NCAU, ISYM, ICAU,
     &        KSING, KDOUB, KORB, KEND1, LEND1, NAMP, IDXLST
    
* external functions:
      INTEGER ILSTSYM
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER IRHSR1
      INTEGER IETA1

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (CCS) THEN
        MODELW = 'CCS       '
        IOPTW  = 1
      ELSE IF (CC2) THEN
        MODELW = 'CC2       '
        IOPTW  = 3
      ELSE IF (CCSD) THEN
        MODELW = 'CCSD      '
        IOPTW  = 3
      ELSE IF (CC3) THEN
        MODELW = 'CC3       '
        IOPTW  = 3
      ELSE
        WRITE (LUPRI,*)'CCCM_VEC> called for unknown coupled '//
     &        'cluster model...'
        CALL QUIT('CCCM_VEC> called for unknown '//
     &            'coupled cluster model...')
      END IF
  

*---------------------------------------------------------------------*
* start loop over right cauchy vectors:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/A)')
     &    ' LIST OF REQUIRED FIRST-ORDER RIGHT CAUCHY VECTORS:'
        DO IDLSTRC = 1, NLRCLBL
           WRITE(LUPRI,'(I5,3X,A8,I5,2X,1P,I5)') 
     &       IDLSTRC,LRCLBL(IDLSTRC),ISYLRC(IDLSTRC),ILRCAU(IDLSTRC)
        END DO
      END IF

      DO ICAU = -1, -4, -1
       DO IDLSTRC = 1, NLRCLBL
        NCAU  = ILRCAU(IDLSTRC)
        ISYM  = ISYLRC(IDLSTRC)
        LABEL = LRCLBL(IDLSTRC)

        IF (ICAU.EQ.NCAU) THEN

          ! allocate work space for one vector:
          KSING = 1                    ! singles
          KDOUB = KSING + NT1AM(ISYM)  ! doubles
          KEND1 = KDOUB + NT2AM(ISYM)
          IF (CCS) KEND1 = KDOUB

          LEND1 = LWORK - KEND1
          IF (LEND1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCCM_VEC.')
          END IF

          KORB = -999999
 
          NAMP = NT1AM(ISYM)
          IF (.NOT. CCS) NAMP = NAMP + NT2AM(ISYM)
        
          IF (NCAU.EQ.-1) THEN
c           -------------
c           C(-1) = -Xksi
c           -------------
            IOPT = 3
            IDXLST = IRHSR1(LABEL,.FALSE.,ZERO,ISYM)
            CALL CC_RDRSP('O1 ',IDXLST,ISYM,IOPT,MODEL,
     &                    WORK(KSING),WORK(KDOUB))
     
            CALL DSCAL(NAMP,-ONE,WORK(KSING),1)

            CALL CC_WRRSP('RC ',IDLSTRC,ISYM,IOPTW,MODELW,WORK(KORB),
     &                    WORK(KSING),WORK(KDOUB),WORK(KEND1),LEND1)

          ELSE IF (NCAU.LT.-1) THEN
c           -------------------
c           C(-n) = A x C(-n+1)
c           -------------------
            IOPT = 3
            IDXLST = ILRCAMP(LABEL,NCAU+1,ISYM)
            CALL CC_RDRSP('RC ',IDXLST,ISYM,IOPT,MODEL,
     &                    WORK(KSING),WORK(KDOUB))
            
            ! Transform XiB vector with Jacobian. Xib is stored at 
            ! the beginning of WORK and will be replaced by the
            ! transformed vector.
            CALL CC_ATRR(ZERO,ISYM,1,WORK(KSING),LWORK,.FALSE.,DUMMY,
     &                   'XXX')
            
            CALL CC_WRRSP('RC ',IDLSTRC,ISYM,IOPTW,MODELW,WORK(KORB),
     &                    WORK(KSING),WORK(KDOUB),WORK(KEND1),LEND1)
          ELSE 
            CALL QUIT('Illegal Cauchy order in CCCM_VEC.')
          END IF

        END IF ! (ICAU.EQ.NCAU)
       END DO ! IDLSTRC
      END DO ! ICAU

*---------------------------------------------------------------------*
* start loop over left cauchy vectors:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
       WRITE (LUPRI,'(/A)')
     &   ' LIST OF REQUIRED FIRST-ORDER LEFT CAUCHY VECTORS:'
       DO IDLSTLC = 1, NLC1LBL
          WRITE(LUPRI,'(I5,3X,A8,I5,2X,1P,I5)') 
     &       IDLSTLC, LBLLC1(IDLSTLC), ISYLC1(IDLSTLC), ILRCAU(IDLSTLC)
       END DO
      END IF

      DO ICAU = -1, -4, -1
       DO IDLSTLC = 1, NLC1LBL
        NCAU  = ILC1CAU(IDLSTLC)
        ISYM  = ISYLC1(IDLSTLC)
        LABEL = LBLLC1(IDLSTLC)

        IF (ICAU.EQ.NCAU) THEN

          ! allocate work space for one vector:
          KSING = 1                    ! singles
          KDOUB = KSING + NT1AM(ISYM)  ! doubles
          KEND1 = KDOUB + NT2AM(ISYM)
          IF (CCS) KEND1 = KDOUB

          LEND1 = LWORK - KEND1
          IF (LEND1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCCM_VEC.')
          END IF

          KORB = -999999
 
          NAMP = NT1AM(ISYM)
          IF (.NOT. CCS) NAMP = NAMP + NT2AM(ISYM)
        
          IF (NCAU.EQ.-1) THEN
c           ---------------
c           Cbar(-1) = +Eta
c           ---------------
            IOPT   = 3
            IDXLST = IETA1(LABEL,.FALSE.,ZERO,ISYM)
            CALL CC_RDRSP('X1 ',IDXLST,ISYM,IOPT,MODEL,
     &                    WORK(KSING),WORK(KDOUB))
     
            CALL CC_WRRSP('LC ',IDLSTLC,ISYM,IOPTW,MODELW,WORK(KORB),
     &                     WORK(KSING),WORK(KDOUB),WORK(KEND1),LEND1)

          ELSE IF (NCAU.LT.-1) THEN
c          -----------------------------------
c          Cbar(-n) = - FC(-n+1) - Cbar(-n+1)A
c          -----------------------------------

            ! calculate - F C(-1):
            IDXLST = ILRCAMP(LABEL,NCAU+1,ISYM)
            CALL CC_FTRAN('L0 ',0,'RC ',IDXLST,WORK(KSING),LWORK)
            CALL DSCAL(NAMP,-ONE,WORK(KSING),1)
            
            ! store as Cbar(-2):
            CALL CC_WRRSP('LC ',IDLSTLC,ISYM,IOPTW,MODELW,WORK(KORB),
     &                     WORK(KSING),WORK(KDOUB),WORK(KEND1),LEND1)
            
            ! read Cbar(-1) and calculate  - Cbar(-1) A:
            IOPT   = 3
            IDXLST = ILC1AMP(LABEL,NCAU+1,ISYM)
            CALL CC_RDRSP('LC ',IDXLST,ISYM,IOPT,MODEL,
     &                     WORK(KSING),WORK(KDOUB))
            
            CALL CC_ATRR(ZERO,ISYM,-1,WORK(KSING),LWORK,.FALSE.,DUMMY,
     &                   'XXX')
            CALL DSCAL(NAMP,-ONE,WORK(KSING),1)
            
            ! add to Cbar(-2) vector on file:
            CALL CC_WARSP('LC',IDLSTLC,ISYM,IOPTW,MODELW,WORK(KORB),
     &                     WORK(KSING),WORK(KDOUB),WORK(KEND1),LEND1)
          ELSE 
            CALL QUIT('Illegal Cauchy order in CCCM_VEC.')
          END IF

        END IF ! (ICAU.EQ.NCAU)
       END DO ! IDLSTLC
      END DO ! ICAU

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCM_VEC                             *
*---------------------------------------------------------------------*
