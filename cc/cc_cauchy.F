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
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
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

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMREF ! symmetry of reference state, assumed total sym.
      PARAMETER (ISYMREF = 1)
    

#if defined (SYS_CRAY)
      REAL ZERO, ONE
#else
      DOUBLE PRECISION ZERO, ONE
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0)

      CHARACTER*80 STRING
      CHARACTER*10 MODEL
      CHARACTER*8  LABELA, LABELB, BLANKS
      INTEGER LWORK
      INTEGER ISYM, ITRAN, IVEC, IDUM, IDXSMOM, IOPER, IOPT, IDISP
      INTEGER ISYMA, IRCVECA, ILCVECA, IOPERA, ISAMA, ISAPROP
      INTEGER ISYMB, IRCVECB, ILCVECB, IOPERB, ISAMB, ISACAU
      INTEGER MXTRAN, MXVEC, MXDOTS, NBCAUCHY, LEND0, KEND0, KRESULT
      INTEGER MXFTRAN, KFCTRAN, KFCDOTS, KFCCONS, NFCTRAN
      INTEGER MXETRAN, KECTRAN, KECDOTS, KECCONS, NECTRAN
      INTEGER MXXTRAN, KXCTRAN, KXCDOTS, KXCCONS, NXCTRAN
      INTEGER KTRK, ISTRT

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TRKBAS, SMOM, SIGN 
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TRKBAS, SMOM, SIGN
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
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD) ) THEN
         CALL QUIT('CC_CAUCHY called for unknown Coupled Cluster.')
      END IF

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
          IVEC = ILRCAMP(LABELA,-1,ISYMA)
          IVEC = ILRCAMP(LABELA,-2,ISYMA)
          IVEC = ILRCAMP(LABELA,-3,ISYMA)
          IVEC = ILRCAMP(LABELA,-4,ISYMA)
          IVEC = ILRCAMP(LABELB, 0,ISYMB)
          IVEC = ILRCAMP(LABELB,-1,ISYMB)
          IVEC = ILRCAMP(LABELB,-2,ISYMB)
          IVEC = ILRCAMP(LABELB,-3,ISYMB)
          IVEC = ILRCAMP(LABELB,-4,ISYMB)
          IF (ASYMSD) THEN
            IVEC = ILC1AMP(LABELA, 0,ISYMA)
            IVEC = ILC1AMP(LABELA,-1,ISYMA)
            IVEC = ILC1AMP(LABELA,-2,ISYMA)
            IVEC = ILC1AMP(LABELA,-3,ISYMA)
            IVEC = ILC1AMP(LABELA,-4,ISYMA)
            IVEC = ILC1AMP(LABELB, 0,ISYMB)
            IVEC = ILC1AMP(LABELB,-1,ISYMB)
            IVEC = ILC1AMP(LABELB,-2,ISYMB)
            IVEC = ILC1AMP(LABELB,-3,ISYMB)
            IVEC = ILC1AMP(LABELB,-4,ISYMB)
          END IF
        ENDIF
      END DO

* close lists:
      LRC1OPN = .FALSE.
      LLC1OPN = .FALSE.

* calculate cauchy vectors for negative orders:
      CALL CCCM_VEC(WORK,LWORK)

*---------------------------------------------------------------------*
* allocate & initialize work space for Cauchy moments
*---------------------------------------------------------------------*

      NBCAUCHY = NLROP * (NLRDISP + 5)

      MXTRAN  = NLRCLBL 
      MXVEC   = NLRCLBL
      MXDOTS  = MXVEC * MXTRAN

      MXFTRAN = 3 * MXTRAN
      MXETRAN = 3 * MXTRAN
      MXXTRAN = 1 * MXTRAN

      KTRK    = 1
      KRESULT = KTRK    + NLROP
      KFCTRAN = KRESULT + NBCAUCHY
      KFCDOTS = KFCTRAN + MXFTRAN
      KFCCONS = KFCDOTS + MXDOTS
      KECTRAN = KFCCONS + MXDOTS
      KECDOTS = KECTRAN + MXETRAN
      KECCONS = KECDOTS + MXDOTS
      KXCTRAN = KECCONS + MXDOTS 
      KXCDOTS = KXCTRAN + MXXTRAN
      KXCCONS = KXCDOTS + MXDOTS
      KEND0   = KXCCONS + MXDOTS
      LEND0   = LWORK - KEND0

*---------------------------------------------------------------------*
* set up lists F transformations and ETA vectors:
*---------------------------------------------------------------------*
      CALL CCCM_SETUP(MXTRAN, MXVEC,
     &                WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,
     &                WORK(KECTRAN),WORK(KECDOTS),NECTRAN,
     &                WORK(KXCTRAN),WORK(KXCDOTS),NXCTRAN )

*---------------------------------------------------------------------*
* calculate F matrix contributions of the type F x RC x RC:
*---------------------------------------------------------------------*
      IOPT = 5 
      CALL CC_FMATRIX(WORK(KFCTRAN),NFCTRAN,'L0 ','RC ',IOPT,'RC ',
     &                WORK(KFCDOTS),WORK(KFCCONS),MXVEC,
     &                WORK(KEND0),LEND0) 

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions of the type ETA{O} x RC:
*---------------------------------------------------------------------*
      CALL CCQR_EADRV('L0 ','o1 ','RC ',NECTRAN, MXVEC,
     &                WORK(KECTRAN),WORK(KECDOTS),WORK(KECCONS),
     &                WORK(KEND0), LEND0, 'DOTP' )

*---------------------------------------------------------------------*
* calculate LC x XKSI{O} vector dot products:
*---------------------------------------------------------------------*
      IF (ASYMSD) THEN
        CALL CC_DOTDRV('LC ','O1 ',NXCTRAN,MXVEC,
     &                 WORK(KXCTRAN), WORK(KXCDOTS), WORK(KXCCONS),
     &                 WORK(KEND0), LEND0 )
      END IF

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
       IF (ISAPROP.EQ.1 .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN

        SMOM = ZERO 

* Eta{B} x Jacobian x C^A(-4):
        IRCVECB = ILRCAMP(LABELB,-4,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERA,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)

* Eta{A} x Jacobian x C^B(-4):
        IRCVECA = ILRCAMP(LABELA,-4,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERB,IRCVECA,ITRAN,IVEC)
        SMOM = SMOM + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)

* F x C^A(-1) x C^B(-3):
        IRCVECA = ILRCAMP(LABELA,-1,IDUM)
        IRCVECB = ILRCAMP(LABELB,-3,IDUM)
        CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                 0,IRCVECA,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

* F x C^A(-2) x C^B(-2):
        IRCVECA = ILRCAMP(LABELA,-2,IDUM)
        IRCVECB = ILRCAMP(LABELB,-2,IDUM)
        CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                 0,IRCVECA,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM - WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

* F x C^A(-3) x C^B(-1):
        IRCVECA = ILRCAMP(LABELA,-3,IDUM)
        IRCVECB = ILRCAMP(LABELB,-1,IDUM)
        CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                 0,IRCVECA,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

* store the result:
        IDXSMOM = (IOPER-1)*(NLRDISP+5) + 1
        WORK(KRESULT-1+IDXSMOM) = SMOM

        IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,2A,I2,5X,A,I3,A,E20.10)')
     &     LABELA,LABELB,IDXSMOM,'S(',2,')=', WORK(KRESULT-1+IDXSMOM)
        END IF

        IF (ASYMSD) THEN

* Eta{A} x C^B(-4):
        IRCVECB = ILRCAMP(LABELB,-4,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERA,IRCVECB,ITRAN,IVEC)
        SMOM = WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        IF (LOCDBG) WRITE (LUPRI,*) 'Eta{A} x C^B:',
     &       WORK(KECCONS-1+(ITRAN-1)*MXVEC+IVEC)

* Xksi{A} x Cbar^B(-4)
        ILCVECB = ILC1AMP(LABELB,-4,IDUM)
        CALL CC_SETDOT(WORK(KXCTRAN),WORK(KXCDOTS),NXCTRAN,MXVEC,
     &                ILCVECB,IOPERA,ITRAN,IVEC)
        SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Xksi{A} x LC^B:',
     &       WORK(KXCCONS-1+(ITRAN-1)*MXVEC+IVEC)
          WRITE(LUPRI,'(5X,3A,I2,5X,A,I3,A,E20.10)') 'ASYM1: ', LABELA,
     &      LABELB, IDXSMOM,'S(',2,')=', SMOM
        END IF


* Eta{B} x C^A(-4):
        IRCVECA = ILRCAMP(LABELA,-4,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERB,IRCVECA,ITRAN,IVEC)
        SMOM = + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)

* Xksi{B} x Cbar^A(-4)
        ILCVECA = ILC1AMP(LABELA,-4,IDUM)
        CALL CC_SETDOT(WORK(KXCTRAN),WORK(KXCDOTS),NXCTRAN,MXVEC,
     &                ILCVECA,IOPERB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

        IF (LOCDBG) WRITE(LUPRI,'(5X,3A,I2,5X,A,I3,A,E20.10)') 
     &   'ASYM2: ', LABELA,LABELB, IDXSMOM,'S(',2,')=', SMOM


        END IF
       END IF 

*---------------------------------------------------------------------*
* calculate zeroth-order Cauchy moments (TRK sum rules):
*---------------------------------------------------------------------*
       IF (ISAPROP.EQ.1 .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN
        SMOM = ZERO 

* Eta{A} x C^B(-2):
        IRCVECB = ILRCAMP(LABELB,-2,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERA,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)


* Eta{B} x C^A(-2):
        IRCVECA = ILRCAMP(LABELA,-2,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERB,IRCVECA,ITRAN,IVEC)
        SMOM = SMOM + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)


* F x C^A(-1) x C^B(-1):
        IRCVECA = ILRCAMP(LABELA,-1,IDUM)
        IRCVECB = ILRCAMP(LABELB,-1,IDUM)
        CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                 0,IRCVECA,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)


* for CCS calculate also pure basis set number ( = RPA result):
        IF (CCS) THEN
          WORK(KTRK-1+IOPER) = 
     &              SMOM + WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        END IF

* store coupled cluster number:
        IDXSMOM = (IOPER-1)*(NLRDISP+5) + 3
        WORK(KRESULT-1+IDXSMOM) = SMOM

        IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,2A,I2,5X,A,I3,A,E20.10)')
     &     LABELA,LABELB,IDXSMOM,'S(',0,')=', WORK(KRESULT-1+IDXSMOM)
        END IF

        IF (ASYMSD) THEN

* Eta{A} x C^B(-2):
        IRCVECB = ILRCAMP(LABELB,-2,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERA,IRCVECB,ITRAN,IVEC)
        SMOM = WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        IF (LOCDBG) WRITE (LUPRI,*) 'Eta{A} x C^B:',
     &       WORK(KECCONS-1+(ITRAN-1)*MXVEC+IVEC)

* Xksi{A} x Cbar^B(-2)
        ILCVECB = ILC1AMP(LABELB,-2,IDUM)
        CALL CC_SETDOT(WORK(KXCTRAN),WORK(KXCDOTS),NXCTRAN,MXVEC,
     &                ILCVECB,IOPERA,ITRAN,IVEC)
        SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Xksi{A} x LC^B:',
     &       WORK(KXCCONS-1+(ITRAN-1)*MXVEC+IVEC)
          WRITE(LUPRI,'(5X,3A,I2,5X,A,I3,A,E20.10)') 
     &       'ASYM1: ', LABELA,LABELB, IDXSMOM,'S(',0,')=', SMOM
        END IF


* Eta{B} x C^A(-2):
        IRCVECA = ILRCAMP(LABELA,-2,IDUM)
        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERB,IRCVECA,ITRAN,IVEC)
        SMOM = + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)

* Xksi{B} x Cbar^A(-2)
        ILCVECA = ILC1AMP(LABELA,-2,IDUM)
        CALL CC_SETDOT(WORK(KXCTRAN),WORK(KXCDOTS),NXCTRAN,MXVEC,
     &                ILCVECA,IOPERB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

        IF (LOCDBG) WRITE(LUPRI,'(5X,3A,I2,5X,A,I3,A,E20.10)') 
     &     'ASYM2: ', LABELA,LABELB, IDXSMOM,'S(',0,')=', SMOM
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

        IRCVECA = ILRCAMP(LABELA,IDISP,IDUM)
        IRCVECB = ILRCAMP(LABELB,IDISP,IDUM)
        SIGN    = DBLE(ISACAU)

        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERA,IRCVECB,ITRAN,IVEC)
        SMOM = SMOM + WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)

        CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                  0,IOPERB,IRCVECA,ITRAN,IVEC)
        SMOM = SMOM + SIGN * WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)


        SIGN = ONE
        DO N = 0, IDISP
          IRCVECA = ILRCAMP(LABELA,      N,IDUM)
          IRCVECB = ILRCAMP(LABELB,IDISP-N,IDUM)

          CALL CCQR_SETF(WORK(KFCTRAN),WORK(KFCDOTS),NFCTRAN,MXVEC,
     &                   0,IRCVECA,IRCVECB,ITRAN,IVEC)
          SMOM = SMOM + SIGN * WORK(KFCCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          SIGN = -SIGN
        END DO

        IDXSMOM = (IOPER-1)*(NLRDISP+5) + IDISP+5
        WORK(KRESULT-1+IDXSMOM) = SMOM

        IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,2A,I2,5X,A,I3,A,E20.10)') LABELA,LABELB,
     &     IDXSMOM,'S(',-IDISP-2,')=', WORK(KRESULT-1+IDXSMOM)
        END IF

        IF (ASYMSD) THEN
          IRCVECA = ILRCAMP(LABELA,IDISP,IDUM)
          IRCVECB = ILRCAMP(LABELB,IDISP,IDUM)
          ILCVECA = ILC1AMP(LABELA,IDISP,IDUM)
          ILCVECB = ILC1AMP(LABELB,IDISP,IDUM)

          CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                    0,IOPERA,IRCVECB,ITRAN,IVEC)
          SMOM = WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          
          CALL CC_SETDOT(WORK(KXCTRAN),WORK(KXCDOTS),NXCTRAN,MXVEC,
     &                  ILCVECB,IOPERA,ITRAN,IVEC)
          SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,3A,I2,5X,A,I3,A,E20.10)') 
     &       'ASYM1: ', LABELA,LABELB, IDXSMOM,'S(',-IDISP-2,')=', SMOM
          END IF


          CALL CCQR_SETEA(WORK(KECTRAN),WORK(KECDOTS),NECTRAN,MXVEC,
     &                    0,IOPERB,IRCVECA,ITRAN,IVEC)
          SMOM = WORK(KECCONS-1 + (ITRAN-1)*MXVEC + IVEC)
          
          CALL CC_SETDOT(WORK(KXCTRAN),WORK(KXCDOTS),NXCTRAN,MXVEC,
     &                  ILCVECA,IOPERB,ITRAN,IVEC)
          SMOM = SMOM + WORK(KXCCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,3A,I2,5X,A,I3,A,E20.10)') 
     &       'ASYM1: ', LABELA,LABELB, IDXSMOM,'S(',-IDISP-2,')=', SMOM
          END IF

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
      ELSE
         CALL QUIT('CC_CAUCHY called for an unkown '//
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
      SUBROUTINE CCCM_SETUP(MXTRAN, MXVEC,
     &                      IFCTRAN, IFCDOTS, NFCTRAN,
     &                      IECTRAN, IECDOTS, NECTRAN,
     &                      IXCTRAN, IXCDOTS, NXCTRAN )
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
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccroper.h>
#include <ccrc1rsp.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCCM_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MXVEC, MXTRAN

      INTEGER IFCTRAN(3,MXTRAN)
      INTEGER IFCDOTS(MXVEC,MXTRAN)
      INTEGER IECTRAN(3,MXTRAN)
      INTEGER IECDOTS(MXVEC,MXTRAN)
      INTEGER IXCTRAN(1,MXTRAN)
      INTEGER IXCDOTS(MXVEC,MXTRAN)

      INTEGER NFRTRAN, NFCTRAN, NXCTRAN, NECTRAN, NCAUMOM

      INTEGER ISYMA,  ISYMB, ITRAN, IVEC, ISYML, IDUM, I, N, IDISP
      INTEGER IRCVECA, IRCVECB, IOPERA, IOPERB, ISAMA, ISAMB
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
       DO I = 1, 3
        IFCTRAN(I,ITRAN) = 0
        IECTRAN(I,ITRAN) = 0
       END DO
       DO I = 1, 1
        IXCTRAN(I,ITRAN) = 0
       END DO

       DO IVEC  = 1, MXVEC
        IFCDOTS(IVEC,ITRAN) = 0
        IECDOTS(IVEC,ITRAN) = 0
        IXCDOTS(IVEC,ITRAN) = 0
       END DO
      END DO

      NFCTRAN = 0
      NECTRAN = 0
      NXCTRAN = 0
      NCAUMOM = 0
 
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
        DO IDISP = -4 , NLRDISP

          ! flag: +1: even, -1: odd 
          ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1

          IF (ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN

            IRCVECA = ILRCAMP(LABELA,IDISP,IDUM)
            IRCVECB = ILRCAMP(LABELB,IDISP,IDUM)
     
            CALL CCQR_SETEA(IECTRAN,IECDOTS,MXTRAN,MXVEC,
     &                      0,IOPERA,IRCVECB,ITRAN,IVEC)
            NECTRAN = MAX(NECTRAN,ITRAN)

            CALL CCQR_SETEA(IECTRAN,IECDOTS,MXTRAN,MXVEC,
     &                      0,IOPERB,IRCVECA,ITRAN,IVEC)
            NECTRAN = MAX(NECTRAN,ITRAN)

          END IF

        END DO

*---------------------------------------------------------------------*
* set up list of F matrix transformations:
*---------------------------------------------------------------------*
        DO IDISP = -4, -2
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

        DO IDISP = 0, NLRDISP
          ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1
          IF (ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLLRDSPCF) THEN
            DO N = 0, IDISP
              IRCVECA = ILRCAMP(LABELA,      N,IDUM)
              IRCVECB = ILRCAMP(LABELB,IDISP-N,IDUM)

              CALL CCQR_SETF(IFCTRAN,IFCDOTS,MXTRAN,MXVEC,
     &                       0,MIN(IRCVECA,IRCVECB),
     &                         MAX(IRCVECA,IRCVECB),ITRAN,IVEC)
              NFCTRAN = MAX(NFCTRAN,ITRAN)
            END DO
          END IF
        END DO

*---------------------------------------------------------------------*
* set up list of XKSI x LC dot products:
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN
          DO IDISP = -4, NLRDISP
            ISACAU = 2 * ( (IDISP/2)*2 - IDISP ) + 1
            IF (ISACAU.EQ.ISAPROP.OR.ISAPROP.EQ.0.OR.ALLLRDSPCF) THEN
              ILCVECA = ILC1AMP(LABELA,IDISP,IDUM)
              ILCVECB = ILC1AMP(LABELB,IDISP,IDUM)
              IRHSA   = IRHSR1(LABELA,.FALSE.,0.0D0,IDUM)
              IRHSB   = IRHSR1(LABELB,.FALSE.,0.0D0,IDUM)
       
              CALL CC_SETDOT(IXCTRAN,IXCDOTS,MXTRAN,MXVEC,
     &                       ILCVECA,IRHSB,ITRAN,IVEC)
              NXCTRAN = MAX(NXCTRAN,ITRAN)

              CALL CC_SETDOT(IXCTRAN,IXCDOTS,MXTRAN,MXVEC,
     &                       ILCVECB,IRHSA,ITRAN,IVEC)
              NXCTRAN = MAX(NXCTRAN,ITRAN)
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
     & ' - ',NECTRAN,' ETA vector calculations '
      WRITE (LUPRI,'(3X,A,/,/)') 'will be performed.'

      IF (LOCDBG) THEN

* list of cauchy vectors:
      WRITE(LUPRI,'(/A)')' LIST OF REQUIRED RIGHT CAUCHY VECTORS:'
      DO I = 1, NLRCLBL
         WRITE(LUPRI,'(I5,3X,A8,I5,2X,1P,I5)') 
     &          I, LRCLBL(I), ISYLRC(I), ILRCAU(I)
      END DO

* F matrix transformations:
      WRITE (LUPRI,*)'List of F matrix transformations:'
      DO ITRAN = 1, NFCTRAN
        WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &   (IFCTRAN(I,ITRAN),I=1,2),(IFCDOTS(I,ITRAN),I=1,MXVEC)
      END DO
      WRITE (LUPRI,*)

* ETA{O} vector calculations:
      WRITE (LUPRI,*) 'List of ETA{O} vector calculations:'
      DO ITRAN = 1, NECTRAN
        WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &   (IECTRAN(I,ITRAN),I=1,2),(IECDOTS(I,ITRAN),I=1,MXVEC)
      END DO
      WRITE (LUPRI,*)

* XKSI{O} x LC dot products:
      WRITE (LUPRI,*) 'List O1 x LC dot products:'
      DO ITRAN = 1, NXCTRAN
        WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &   IXCTRAN(1,ITRAN),(IXCDOTS(I,ITRAN),I=1,MXVEC)
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
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <dummy_dec.h>
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
      LOGICAL INCORE
      PARAMETER (INCORE = .TRUE.)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION WORK(LWORK)
#endif

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABELA, LABELB
      INTEGER ILROP, IVECR, IVECR1, IVECL, IVECL1, KAXI1, KAXI2, KEND1
      INTEGER LEND1, ISYML, IOPT, IOPERA, IOPERB, ISYMA, ISYMB, KAXI0
 
* external functions:
      INTEGER ILSTSYM
      INTEGER ILRCAMP
      INTEGER ILC1AMP

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (CCS) THEN
        MODEL = 'CCS       '
        IOPT  = 1
      ELSE IF (CC2) THEN
        MODEL = 'CC2       '
        IOPT  = 3
      ELSE IF (CCSD) THEN
        MODEL = 'CCSD      '
        IOPT  = 3
      ELSE
        WRITE (LUPRI,*)'CCCM_VEC> called for unkown coupled '//
     &        'cluster model...'
        CALL QUIT('CCCM_VEC> called for unkown '//
     &            'coupled cluster model...')
      END IF
  
      IF (LOCDBG) THEN
       WRITE (LUPRI,'(/A)')
     &   ' LIST OF REQUIRED FIRST-ORDER RIGHT CAUCHY VECTORS:'
       DO I = 1, NLRCLBL
          WRITE(LUPRI,'(I5,3X,A8,I5,2X,1P,I5)') 
     &          I, LRCLBL(I), ISYLRC(I), ILRCAU(I)
       END DO
       WRITE (LUPRI,'(/A)')
     &   ' LIST OF REQUIRED FIRST-ORDER LEFT CAUCHY VECTORS:'
       DO I = 1, NLRCLBL
          WRITE(LUPRI,'(I5,3X,A8,I5,2X,1P,I5)') 
     &          I, LBLLC1(I), ISYLC1(I), ILRCAU(I)
       END DO
      END IF
*---------------------------------------------------------------------*
* start loop over operator labels for A operators:
*---------------------------------------------------------------------*

      DO ILROP = 1, NLROP
        IOPERA = IALROP(ILROP)
        IOPERB = IBLROP(ILROP)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)

C       WRITE (LUPRI,*) MSGDBG, 'LABELA:',LABELA

* allocate work space for one vector:
        KAXI1 = 1                     ! Xi, singles
        KAXI2 = KAXI1 + NT1AM(ISYMA)  ! Xi, doubles
        KEND1 = KAXI2 + NT2AM(ISYMA)
C       IF (CCS) KEND1 = KAXI2
        LEND1 = LWORK - KEND1
      
        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCCM_VEC.')
        END IF

        KAXI0 = -999999
 
*---------------------------------------------------------------------*
*       C(-1) = -Xksi
*---------------------------------------------------------------------*
        CALL CC_XKSI(WORK(KAXI1),LABELA,ISYMA,0,DUMMY,WORK(KEND1),LEND1)
        CALL DSCAL(NT1AM(ISYMA),-1.0d0,WORK(KAXI1),1)
        IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMA),-1.0d0,WORK(KAXI2),1)
     
        IVECR = ILRCAMP(LABELA,-1,ISYMA)
        CALL CC_WRRSP('RC',IVECR,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-1) = Eta
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN
          CALL CC_ETAC(ISYMA,LABELA,WORK(KAXI1),'L0',0,0,DUMMY,
     &                 WORK(KEND1),LEND1)
     
          IVECL = ILC1AMP(LABELA,-1,ISYMA)
          CALL CC_WRRSP('LC',IVECL,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
*       C(-2) = A x C(-1)
*---------------------------------------------------------------------*
*       Transform XiB vector with Jacobian. Xib is stored at the 
*       beginning of WORK and will be replaced by transformed vector.

        IVECR1 = ILRCAMP(LABELA,-1,ISYMA)
        CALL CC_RDRSP('RC',IVECR1,ISYMA,IOPT,MODEL,
     &                WORK(KAXI1),WORK(KAXI2))

        CALL CC_ATRR(ISYMA,1,WORK(KAXI1),LWORK)
 
        IVECR = ILRCAMP(LABELA,-2,ISYMA)
        CALL CC_WRRSP('RC',IVECR,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-2) = - FC(-1) - Cbar(-1)A
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN 

* calculate - F C(-1):
          IVECR = ILRCAMP(LABELA,-1,ISYMA)
          CALL CC_FTRAN('L0 ',0,'RC ',IVECR,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMA),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMA),-1.0d0,WORK(KAXI2),1)

* store as Cbar(-2):
          IVECL = ILC1AMP(LABELA,-2,ISYMA)
          CALL CC_WRRSP('LC',IVECL,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

* read Cbar(-1) and calculate  - Cbar(-1) A:
          IVECL1 = ILC1AMP(LABELA,-1,ISYMA)
          CALL CC_RDRSP('LC',IVECL1,ISYMA,IOPT,MODEL,
     &                   WORK(KAXI1),WORK(KAXI2))
          
          CALL CC_ATRR(ISYMA,-1,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMA),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMA),-1.0d0,WORK(KAXI2),1)

* add to Cbar(-2) vector on file:
          IVECL = ILC1AMP(LABELA,-2,ISYMA)
          CALL CC_WARSP('LC',IVECL,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
*       C(-3) = A x C(-2)
*---------------------------------------------------------------------*
*       Transform XiB vector with Jacobian. Xib is stored at the 
*       beginning of WORK and will be replaced by transformed vector.

        IVECR1 = ILRCAMP(LABELA,-2,ISYMA)
        CALL CC_RDRSP('RC',IVECR1,ISYMA,IOPT,MODEL,
     &                WORK(KAXI1),WORK(KAXI2))

        CALL CC_ATRR(ISYMA,1,WORK(KAXI1),LWORK)
 
        IVECR = ILRCAMP(LABELA,-3,ISYMA)
        CALL CC_WRRSP('RC',IVECR,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-3) = - FC(-2) - Cbar(-2)A  
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN 

* calculate - F C(-2):
          IVECR = ILRCAMP(LABELA,-2,ISYMA)
          CALL CC_FTRAN('L0 ',0,'RC ',IVECR,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMA),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMA),-1.0d0,WORK(KAXI2),1)

* store as Cbar(-3):
          IVECL = ILC1AMP(LABELA,-3,ISYMA)
          CALL CC_WRRSP('LC',IVECL,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

* read Cbar(-2) and calculate  - Cbar(-2) A:
          IVECL1 = ILC1AMP(LABELA,-2,ISYMA)
          CALL CC_RDRSP('LC',IVECL1,ISYMA,IOPT,MODEL,
     &                   WORK(KAXI1),WORK(KAXI2))
          
          CALL CC_ATRR(ISYMA,-1,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMA),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMA),-1.0d0,WORK(KAXI2),1)

* add to Cbar(-3) vector on file:
          IVECL = ILC1AMP(LABELA,-3,ISYMA)
          CALL CC_WARSP('LC',IVECL,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
*       C(-4) = A x C(-3)
*---------------------------------------------------------------------*
*       Transform XiB vector with Jacobian. Xib is stored at the 
*       beginning of WORK and will be replaced by transformed vector.

        IVECR1 = ILRCAMP(LABELA,-3,ISYMA)
        CALL CC_RDRSP('RC',IVECR1,ISYMA,IOPT,MODEL,
     &                WORK(KAXI1),WORK(KAXI2))

        CALL CC_ATRR(ISYMA,1,WORK(KAXI1),LWORK)
 
        IVECR = ILRCAMP(LABELA,-4,ISYMA)
        CALL CC_WRRSP('RC',IVECR,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-4) = - FC(-3) - Cbar(-3)A
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN 

* calculate - F C(-3):
          IVECR = ILRCAMP(LABELA,-3,ISYMA)
          CALL CC_FTRAN('L0 ',0,'RC ',IVECR,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMA),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMA),-1.0d0,WORK(KAXI2),1)

* store as Cbar(-4):
          IVECL = ILC1AMP(LABELA,-4,ISYMA)
          CALL CC_WRRSP('LC',IVECL,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

* read Cbar(-3) and calculate  - Cbar(-3) A:
          IVECL1 = ILC1AMP(LABELA,-3,ISYMA)
          CALL CC_RDRSP('LC',IVECL1,ISYMA,IOPT,MODEL,
     &                   WORK(KAXI1),WORK(KAXI2))
          
          CALL CC_ATRR(ISYMA,-1,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMA),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMA),-1.0d0,WORK(KAXI2),1)

* add to Cbar(-4) vector on file:
          IVECL = ILC1AMP(LABELA,-4,ISYMA)
          CALL CC_WARSP('LC',IVECL,ISYMA,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
      END DO


*---------------------------------------------------------------------*
* start loop over operator labels for B operators:
*---------------------------------------------------------------------*
      DO ILROP = 1, NLROP
        IOPERA = IALROP(ILROP)
        IOPERB = IBLROP(ILROP)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)

C       WRITE (LUPRI,*) MSGDBG, 'LABELB:',LABELB

* allocate work space for one vector:
        KAXI1 = 1                    ! Xi, singles
        KAXI2 = KAXI1 + NT1AM(ISYMB)  ! Xi, doubles
        KEND1 = KAXI2 + NT2AM(ISYMB)
C       IF (CCS) KEND1 = KAXI2
        LEND1 = LWORK - KEND1
      
        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCCM_VEC.')
        END IF
 
*---------------------------------------------------------------------*
*       C(-1) = -Xksi
*---------------------------------------------------------------------*
        CALL CC_XKSI(WORK(KAXI1),LABELB,ISYMB,0,DUMMY,WORK(KEND1),LEND1)
        CALL DSCAL(NT1AM(ISYMB),-1.0d0,WORK(KAXI1),1)
        IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMB),-1.0d0,WORK(KAXI2),1)
     
        IVECR = ILRCAMP(LABELB,-1,ISYMB)
        CALL CC_WRRSP('RC',IVECR,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-1) = Eta
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN
          CALL CC_ETAC(ISYMB,LABELB,WORK(KAXI1),'L0',0,0,DUMMY,
     &                 WORK(KEND1),LEND1)
     
          IVECL = ILC1AMP(LABELB,-1,ISYMB)
          CALL CC_WRRSP('LC',IVECL,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
*       C(-2) = A x C(-1)
*---------------------------------------------------------------------*
*       Transform XiB vector with Jacobian. Xib is stored at the 
*       beginning of WORK and will be replaced by transformed vector.

        IVECR1 = ILRCAMP(LABELB,-1,ISYMB)
        CALL CC_RDRSP('RC',IVECR1,ISYMB,IOPT,MODEL,
     &                WORK(KAXI1),WORK(KAXI2))

        CALL CC_ATRR(ISYMB,1,WORK(KAXI1),LWORK)
 
        IVECR = ILRCAMP(LABELB,-2,ISYMB)
        CALL CC_WRRSP('RC',IVECR,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-2) = - FC(-1) - Cbar(-1)A
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN 

* calculate - F C(-1):
          IVECR = ILRCAMP(LABELB,-1,ISYMB)
          CALL CC_FTRAN('L0 ',0,'RC ',IVECR,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMB),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMB),-1.0d0,WORK(KAXI2),1)

* store as Cbar(-2):
          IVECL = ILC1AMP(LABELB,-2,ISYMB)
          CALL CC_WRRSP('LC',IVECL,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

* read Cbar(-1) and calculate  - Cbar(-1) A:
          IVECL1 = ILC1AMP(LABELB,-1,ISYMB)
          CALL CC_RDRSP('LC',IVECL1,ISYMB,IOPT,MODEL,
     &                   WORK(KAXI1),WORK(KAXI2))
          
          CALL CC_ATRR(ISYMB,-1,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMB),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMB),-1.0d0,WORK(KAXI2),1)

* add to Cbar(-2) vector on file:
          IVECL = ILC1AMP(LABELB,-2,ISYMB)
          CALL CC_WARSP('LC',IVECL,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
*       C(-3) = A x C(-2)
*---------------------------------------------------------------------*
*       Transform XiB vector with Jacobian. Xib is stored at the 
*       beginning of WORK and will be replaced by transformed vector.

        IVECR1 = ILRCAMP(LABELB,-2,ISYMB)
        CALL CC_RDRSP('RC',IVECR1,ISYMB,IOPT,MODEL,
     &                WORK(KAXI1),WORK(KAXI2))

        CALL CC_ATRR(ISYMB,1,WORK(KAXI1),LWORK)
 
        IVECR = ILRCAMP(LABELB,-3,ISYMB)
        CALL CC_WRRSP('RC',IVECR,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-3) = - FC(-2) - Cbar(-2)A 
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN 

* calculate - F C(-2):
          IVECR = ILRCAMP(LABELB,-2,ISYMB)
          CALL CC_FTRAN('L0 ',0,'RC ',IVECR,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMB),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMB),-1.0d0,WORK(KAXI2),1)

* store as Cbar(-3):
          IVECL = ILC1AMP(LABELB,-3,ISYMB)
          CALL CC_WRRSP('LC',IVECL,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

* read Cbar(-2) and calculate  - Cbar(-2) A:
          IVECL1 = ILC1AMP(LABELB,-2,ISYMB)
          CALL CC_RDRSP('LC',IVECL1,ISYMB,IOPT,MODEL,
     &                   WORK(KAXI1),WORK(KAXI2))
          
          CALL CC_ATRR(ISYMB,-1,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMB),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMB),-1.0d0,WORK(KAXI2),1)

* add to Cbar(-3) vector on file:
          IVECL = ILC1AMP(LABELB,-3,ISYMB)
          CALL CC_WARSP('LC',IVECL,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
*       C(-4) = A x C(-3)
*---------------------------------------------------------------------*
*       Transform XiB vector with Jacobian. Xib is stored at the 
*       beginning of WORK and will be replaced by transformed vector.

        IVECR1 = ILRCAMP(LABELB,-3,ISYMB)
        CALL CC_RDRSP('RC',IVECR1,ISYMB,IOPT,MODEL,
     &                WORK(KAXI1),WORK(KAXI2))

        CALL CC_ATRR(ISYMB,1,WORK(KAXI1),LWORK)
 
        IVECR = ILRCAMP(LABELB,-4,ISYMB)
        CALL CC_WRRSP('RC',IVECR,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
*       Cbar(-4) = - FC(-3) - Cbar(-3)A
*---------------------------------------------------------------------*
        IF (ASYMSD) THEN 

* calculate - F C(-3):
          IVECR = ILRCAMP(LABELB,-3,ISYMB)
          CALL CC_FTRAN('L0 ',0,'RC ',IVECR,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMB),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMB),-1.0d0,WORK(KAXI2),1)

* store as Cbar(-4):
          IVECL = ILC1AMP(LABELB,-4,ISYMB)
          CALL CC_WRRSP('LC',IVECL,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)

* read Cbar(-3) and calculate  - Cbar(-3) A:
          IVECL1 = ILC1AMP(LABELB,-3,ISYMB)
          CALL CC_RDRSP('LC',IVECL1,ISYMB,IOPT,MODEL,
     &                   WORK(KAXI1),WORK(KAXI2))
          
          CALL CC_ATRR(ISYMB,-1,WORK(KAXI1),LWORK)
          CALL DSCAL(NT1AM(ISYMB),-1.0d0,WORK(KAXI1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AM(ISYMB),-1.0d0,WORK(KAXI2),1)

* add to Cbar(-4) vector on file:
          IVECL = ILC1AMP(LABELB,-4,ISYMB)
          CALL CC_WARSP('LC',IVECL,ISYMB,IOPT,MODEL,WORK(KAXI0),
     &                   WORK(KAXI1),WORK(KAXI2),WORK(KEND1),LEND1)
        END IF
*---------------------------------------------------------------------*
      END DO

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCM_VEC                             *
*---------------------------------------------------------------------*
