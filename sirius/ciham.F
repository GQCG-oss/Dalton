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
C  /* Deck ciham */
      SUBROUTINE CIHAM(IWAY,MXPRDT,ICSYM,IDET,CIDIA,NCONF,
     *                 ECORE,FCAC,H2AC,XNDXCI,WRK,LWRK,IPRINT)
C
C 29-Jan-1989 hjaaj
C
C Explicit CI Hamiltonian.
C
C Notes:
C  - introduce symtest using IHSYM=0
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IDET(*),CIDIA(*),FCAC(*),H2AC(*),XNDXCI(*),WRK(LWRK)
C
C Used from common blocks:
C   INFINP : MCTYPE
C   INFORB : MULD2H(8,8),NASHT,N2ASHX
C   CIINFO : ICOMBI, IPSIGN
C   DETBAS : KLTSOB
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <ciinfo.h>
#include <detbas.h>
C
C     Set some codes:
C
      PSIGN  = IPSIGN
      NPRDET = 0
      IHSYM  = 1
C     ... Hamiltonian matrix is always symmetric in this version.
      IDODGN = 1
C     ... always do diagonalization in this version
      NROOT  = MIN(IPRINT,NCONF,MXPRDT)
C     ... temporary analyze specification
C
C     Core allocation:
C
      KHAMIL = 1
      IF (IHSYM .EQ. 1) THEN
         KFCAC2 = KHAMIL + (MXPRDT**2 + MXPRDT)/2
      ELSE
         KFCAC2 = KHAMIL + MXPRDT**2
      END IF
      KRJ    = KFCAC2 + N2ASHX
      KRK    = KRJ    + N2ASHX
      KWHAM  = KRK    + N2ASHX
      IF (IDODGN .EQ. 1) THEN
         KEIGVL = KWHAM
         KEIGVC = KEIGVL + MXPRDT
         KWHAM  = KEIGVC + MXPRDT ** 2
      ELSE
         KEIGVL = 999 999 999
         KEIGVC = 999 999 999
      END IF
C
      CALL DSPTSI(NASHT,FCAC,WRK(KFCAC2))
      CALL GETFIJ(WRK(KRJ),WRK(KRK),H2AC)
      IF (MCTYPE .EQ. 2) THEN
C        ... reorder from Sirius to Lunar order for RAS
         CALL REORMT(WRK(KFCAC2),WRK(KWHAM),NASHT,NASHT,
     &               XNDXCI(KLTSOB),XNDXCI(KLTSOB))
         CALL DCOPY(N2ASHX,WRK(KWHAM),1,WRK(KFCAC2),1)
         CALL REORMT(WRK(KRJ),WRK(KWHAM),NASHT,NASHT,
     &               XNDXCI(KLTSOB),XNDXCI(KLTSOB))
         CALL DCOPY(N2ASHX,WRK(KWHAM),1,WRK(KRJ),1)
         CALL REORMT(WRK(KRK),WRK(KWHAM),NASHT,NASHT,
     &               XNDXCI(KLTSOB),XNDXCI(KLTSOB))
         CALL DCOPY(N2ASHX,WRK(KWHAM),1,WRK(KRK),1)
      END IF
      CALL PHPDET(CIDIA,IWAY,MXPRDT,XNDXCI,WRK,KWHAM,LWRK,
     &            ICSYM,IHSYM,WRK(KHAMIL),
     &            WRK(KFCAC2),WRK(KRJ),WRK(KRK),NASHT,MULD2H,
     &            ECORE,ICOMBI,PSIGN,NPRDET,WRK(KEIGVC),WRK(KEIGVL),
     &            IDODGN,NCONF,IDET,NROOT,H2AC,IPRINT)
C     CALL PHPDET(CIDIA,IWAY,MXPRDT,XNDXCI,WORK,KFREE,LFREE,
C    &            ICSYM,IHSYM,HAMIL,ONEBOD,RJ,RK,NORB,SYMPRO,
C    &            ECORE,ICOMBI,PSIGN,NPRDET,EIGVC,EIGVL,
C    &            IDODGN,NDET,IDET,NROOT,FIJKL,NTEST)
C
      IF (NPRDET .NE. MXPRDT) THEN
         IF (IPRINT .GT. 0) THEN
            WRITE (LUPRI,*) 'Number of explicit determinants'
            WRITE (LUPRI,*) '   -- requested maximum',MXPRDT
            WRITE (LUPRI,*) '   -- actually used    ',NPRDET
         END IF
         MXPRDT = NPRDET
      END IF
      RETURN
      END
C  /* Deck phpdet */
      SUBROUTINE PHPDET(CIDIA,IWAY,MXPRDT,XNDXCI,WORK,KFREE,LFREE,
     &                  ICSYM,IHSYM,HAMIL,ONEBOD,RJ,RK,NORB,SYMPRO,
     &                  ECORE,ICOMBI,PSIGN,NPRDET,EIGVC,EIGVL,
     &                  IDODGN,NDET,IDET,NROOT,FIJKL,NTEST)
C
C..  SELECT A SUBSPACE AND CONSTRUCT HAMLITONIAN MATRIX IN THIS SUBSPACE
C
C  PARAMETERS :
C      CIDIA  : CI diagonal (input, DESTROYED if nroot.ne.0 )
C      IWAY   : flag for selecting subspace ( = 1 normal see below)
C      MXPRDT : largest allowed dimension of subspace
C      XNDXCI : string information
C      WORK   : scratch space, length at least ??
C      KFREE  : first free element in work ( = 1 in normal cases )
C      LFREE  : dimension WORK(LFREE)
C      ICSYM  : symmetry of CI space
C      IHSYM  : = 0 : calculate complete subspace hamiltonian matrix
C               = 1 : only lower part of hamiltonian matrix calculated
C      HAMIL  : space for hamilton matrix. DESTROYED if IDOGN .ne. 0
C      ONEBOD : onebody matrix in square form (FI for active orbitals )
C      RJ     : Coulomb integrals in square form
C      RK     : Exchange integrals in square form
C      NORB   : number of active orbitals
C      SYMPRO : 8*8 D2h table
C      ECORE  : core energy
C      ICOMBI : = 0 : use determinant basis
C               = 1 : use combination basis
C      PSIGN  : permutation sign for combinations ( only ICOMBI = 1 )
C      NPRDET : Dimension of subspace used (OUTPUT ! )
C      EIGVC  : eigenvectors stored in column form ( IDOGN .NE. 0 )
C      EIGVL  : eigenvalues of subspace hamiltonian( IDOGN .NE. 0 )
C      IDOGN  : .ne. 0 : diagonalize subspace matrix
C      IDET   : POINTER : IDET(I) gives adress of supspace element i in
C               FULL space ( output ). If IWAY = 3, also input.
C      NDET   : dimension of full space
C      NROOT  : number of roots analyzed with anaci ( IDOGN .ne. 0)
C      FIJKL  : array relayed to H2TUVX routine in DIHDJ
C      NTEST  : print level
C
C
C.. SUBSPACE :
C     THE SUBSPACE IS CHOOSEN AS
C     IWAY = 1 : CHOOSE THE LOWEST VALUES OF
C                THE CI DIAGONAL.THE NUMBER OF DETERMINANTS IS
C                CHOOSEN SO NO DEGENERATE LEVELS ARE SPLITTED .
C                THE NUMBER OF DETERMINANTS USED , NPRDET, CAN
C                THEREFORE BE LOWER THAN MXPRDT.
C     IWAY = 2 : CHOOSE THE FIRST MXPRDT DETERMINANTS
C                STUPID, BUT CONVENIENT FOR TESTING
C     IWAY = 3 : Use input list in IDET(1:MXPRDT). List is not checked.
C     IWAY = 4 : choose the elements with highest absolut value
C                in CIDIA. (e.g. sigma vector for HF determinant).
C
C IDET CONTAINS ADRESSES OF ELEMENTS CHOSEN
C
C IDODGN .GT. 0 : DIAGONALIZE CONSTRUCTED HAMILTON MATRIX.
C                 EIGVL CONTAINS EIGENVALUES ON RETURN
C                 EIGVC CONTAINS EIGENVECTORS(COLUMNS) ON RETURN
C
C JEPPE OLSEN JANUARY 1989
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XNDXCI(*),WORK(*),HAMIL(*),CIDIA(*)
      DIMENSION ONEBOD(*),RJ(*),RK(*),FIJKL(*)
      DIMENSION IDET(MXPRDT), EIGVL(MXPRDT) , EIGVC(MXPRDT ** 2 )
      INTEGER SYMPRO(8,8)
C
C Used from common blocks:
C   STRNUM : ...
C   DETBAS : KLTSOB
C
#include <mxpdim.h>
#include <strnum.h>
#include <detbas.h>
C
      KFREEO = KFREE
C
      IF (NTEST .GT. 0) WRITE (LUPRI,*) '>>> Output from PHPDET <<<'
      IF (NTEST .GT. 10) THEN
         IF( ICOMBI .EQ. 0 ) THEN
             WRITE(LUPRI,*) ' DETERMINANTS IN USE '
         ELSE
             WRITE(LUPRI,*) ' COMBINATIONS IN USE '
         END IF
         WRITE(LUPRI,*) ' MXPRDT ',MXPRDT
      END IF
C
C.. 0 : SELECT SUBSPACE
C
      IF( IWAY .EQ. 1) THEN
C.      FIND NUMBER OF DETERMINANTS LESS OR EQUAL TOP MXPRDT
C       THAT DOES NOT SEPARATE DEGENERATE PAIRS .
C            FNDMN3(VEC,NDIM,MXELMN,IPLACE,NELMN,IPRT,THRES)
        CALL FNDMN3(CIDIA,NDET,MXPRDT,IDET,NPRDET,NTEST,1.0D-3)
      ELSE IF ( IWAY .EQ. 2 ) THEN
        NPRDET = MIN(MXPRDT,NDET)
        CALL ISTVC2(IDET,0,1,NPRDET)
      ELSE IF ( IWAY .EQ. 3 ) THEN
        NPRDET = MXPRDT
C       ... we assume caller supplies valid IDET(1:MXPRDT) and
C           MXPRDT .LE. NDET
      ELSE IF ( IWAY .EQ. 4 ) THEN
        CALL FNDAMX(CIDIA,NDET,MXPRDT,IDET,NPRDET,NTEST,1.0D-3)
      ELSE
        WRITE (LUPRI,*) 'Fatal error in PHPDET, illegal IWAY =',IWAY
        CALL QUIT('Fatal error in PHPDET, illegal IWAY.')
      END IF
      IF (NTEST .GT. 10) WRITE(LUPRI,*) ' IWAY,NPRDET .... ', IWAY,
     &                                    NPRDET
      IF (NTEST .GE. 10) THEN
        WRITE(LUPRI,*) ' IDET IN PHPDET, IWAY = ',IWAY
        CALL IWRTMA(IDET,1,NPRDET,1,NPRDET)
      END IF
C
C.. 1 : LOCAL MEMORY
C
C SPACE FOR STRINGS OF SUBSPACE DETERMINANTS
      CALL MEMADD(KIAST2,NPRDET*NAEL,KFREE,1)
      CALL MEMADD(KIBST2,NPRDET*NBEL,KFREE,1)
C FOR COMBINATIONS AN EXTRA SET IS NEEDED
      IF( ICOMBI.NE. 0 ) THEN
        CALL MEMADD(KIAST3,NPRDET*NAEL,KFREE,1)
        CALL MEMADD(KIBST3,NPRDET*NBEL,KFREE,1)
      ELSE
        KIAST3 = KIAST2
        KIBST3 = KIBST2
      END IF
C SCRATCH SPACE FOR DIHDJ
      LSCR = 4 * NORB + NPRDET
      CALL MEMADD(KSCR, LSCR,KFREE,1)
      IF (KFREE .GT. LFREE) THEN
        WRITE (LUPRI,*) 'Fatal error, insufficient memory in PHPDET'
        WRITE (LUPRI,*) 'Need:',KFREE,'  Available:',LFREE
        CALL QUIT('Fatal error, insufficient memory in PHPDET')
      END IF
C
C.. 3 : OBTAIN STRINGS CORRESPONDING TO SUBSPACE DETERMINANTS
C
      CALL STRFDT(NPRDET,IDET,WORK(KIAST2),WORK(KIBST2),
     &            ICSYM,SYMPRO,XNDXCI,ICOMBI)
       IF( ICOMBI .NE. 0 ) THEN
         CALL ICOPVE(WORK(KIAST2),WORK(KIAST3),NPRDET*NAEL)
         CALL ICOPVE(WORK(KIBST2),WORK(KIBST3),NPRDET*NBEL)
       END IF
C
C.. 4 : OBTAIN CORRESPONDING HAMILTONIAN
C
      CALL DIHDJ(WORK(KIAST2),WORK(KIBST2),NPRDET,WORK(KIAST3),
     &           WORK(KIBST3),NPRDET,NAEL,NBEL,WORK(KSCR),LSCR,
     &           NORB,ONEBOD,RJ,RK,FIJKL,XNDXCI(KLTSOB),
     &           HAMIL,IHSYM,ECORE,ICOMBI,PSIGN)
C
C.. 5 : DIAGONALIZE HAMIL ( ASSUMES IHSYM .NE. 0 USED )
C
      IF( IDODGN .GT. 0 ) THEN
        MROOT = MIN(NPRDET,NROOT)
        CALL EIGEN(HAMIL,EIGVC,NPRDET,0,1)
C.       EXTRACT EIGENVALUES
        CALL XTRCDI(HAMIL,EIGVL,NPRDET,1)
        IF((NTEST.GE.2.AND.NTEST.LT.4) .OR. NTEST.GE.10 ) THEN
         WRITE(LUPRI,'(/A)') ' EIGENVALUES OF SUBSPACE HAMILTONIAN '
         CALL WRTMAT(EIGVL,1,NPRDET,1,NPRDET,0)
         IF( NTEST .GE. 20 ) THEN
           WRITE(LUPRI,'(/A)') ' EIGENVECTORS OF SUBSPACE HAMILTONIAN '
           CALL OUTPUT(EIGVC,1,NPRDET,1,NPRDET,NPRDET,NPRDET,1,LUPRI)
         ELSE IF(NTEST .GE. 10 .AND. MROOT .GT. 0) THEN
           WRITE(LUPRI,'(/A)') ' EIGENVECTORS OF SUBSPACE HAMILTONIAN '
           CALL OUTPUT(EIGVC,1,NPRDET,1,MROOT,NPRDET,NPRDET,1,LUPRI)
         END IF
        END IF
        IF (NTEST.GT.0) WRITE(LUPRI,*) ' RANGE OF EIGENVALUES IN '//
     &               'SUBSPACE ',EIGVL(1),EIGVL(NPRDET)
C
C .. 6 : ANALYZE MROOT APPROXIMATIVE EIGENVECTORS
C
      IF (NTEST .GE. 4) THEN
      DO 1869 IROOT = 1, MROOT
        IOFF = (IROOT-1)*NPRDET + 1
        CALL SETVEC(CIDIA,0.0D0,NDET)
        CALL SCAVEC(CIDIA,EIGVC(IOFF),IDET,NPRDET)
        WRITE(LUPRI,'(/A,I3)') '  Information about root ... ',IROOT
        WRITE(LUPRI,'(A)')     '  ******************************'
        WRITE(LUPRI,'(/A,1P,E15.8)') '   Energy .... ',EIGVL(IROOT)
        CUTOFF = 0.1D0/SQRT(10.0D0)
        CALL ANACIN(CIDIA,ICSYM,CUTOFF,100,XNDXCI,SYMPRO,6,0)
C       CALL ANACIN(CIVEC,ICSYM,THRES,MAXTRM,XNDXCI,SYMPRO,IOUT,INCSFB)
 1869 CONTINUE
      END IF
      END IF
C
      KFREE = KFREEO
      RETURN
      END
C  /* Deck strfdt */
      SUBROUTINE STRFDT(NDET,IDET,IASTR2,IBSTR2,ICSYM,SYMPRO,
     &   XNDXCI,JCOMBI)
C
C NDET DETERMINANTS ARE GIVEN IN ARRAY IDET
C OBTAIN THE CORRESPONDING ALPHA AND BETASTRINGS AND STORE THEM
C IN IASTR2 AND IBSTR2
C
C NOTE ICOMBI IS INTRODUCED THROUGH CIINFO IN THIS IMPLEMENTATION
C
#include <implicit.h>
      INTEGER SYMPRO(8,8)
      DIMENSION IDET(*),IASTR2(*),IBSTR2(*),XNDXCI(*)
C
#include <mxpdim.h>
#include <ciinfo.h>
#include <detbas.h>
#include <strnum.h>
C
         CALL STFDT2(NDET,IDET,IASTR2,IBSTR2,MAXSYM,NOCTPA,NOCTPB,
     &               XNDXCI(KNSSOA),XNDXCI(KNSSOB),XNDXCI(KISSOA),
     &               XNDXCI(KISSOB),XNDXCI(KIOCOC),NAEL,NBEL,
     &               SYMPRO,XNDXCI(KIASTR),XNDXCI(KIBSTR),ICSYM,ICOMBI)
C     SUBROUTINE STFDT2(NDET,IDET,IASTR2,IBSTR2,MAXSYM,NOCTPA,NOCTPB,
C    &                  NBEL,SYMPRO,IASTR,IBSTR,ISYM,ICOMBI)
C
      RETURN
      END
C  /* Deck dihdj */
      SUBROUTINE DIHDJ(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,
     &   IWORK,LWORK,NORB,ONEBOD,RJ,RK,FIJKL,ILTSOB,
     &   HAMIL,ISYM,ECORE,ICOMBI,PSIGN)
C
C A SET OF DETERMINANTS IA DEFINED BY ALPHA AND BETASTRINGS
C IASTR,IBSTR AND ANOTHER SET OF DETERMINATS DEFINED BY STRINGS
C JASTR AND JBSTR ARE GIVEN . OBTAIN CORRESPONDING HAMILTONIAN MATRIX
C
C IF ICOMBI .NE. 0 COMBINATIONS ARE USED FOR ALPHA AND BETA STRING
C THAT DIFFERS :
C   1/SQRT(2) * ( |I1A I2B| + PSIGN * |I2A I1B| )
C
C IF ISYM .EQ. 0 FULL HAMILTONIAN IS CONSTRUCTED
C IF ISYM .NE. 0 LOWER HALF OF HAMILTONIAN IS CONSTRUCTED
C
C JEPPE OLSEN JANUARY 1989
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION JASTR(NAEL,*),JBSTR(NBEL,*)
C
      DIMENSION IWORK(*), HAMIL(*), ONEBOD(NORB,NORB),
     *          RJ(NORB,NORB), RK(NORB,NORB), FIJKL(*), ILTSOB(*)
C
C SCRATCH SPACE
C
C .. 1 : EXPANSION OF ALPHA AND BETA STRINGS OF TYPE I
C
C?    WRITE(LUPRI,*) ' MISTER DIHDJ SPEAKING  ICOMBI = ', ICOMBI
      KLFREE = 1
      KLIAE  = KLFREE
      KLFREE = KLIAE + NORB
      KLIBE  = KLFREE
      KLFREE = KLIBE + NORB
C
      KLJAE  = KLFREE
      KLFREE = KLJAE + NORB
      KLJBE  = KLFREE
      KLFREE = KLJBE + NORB
      IF( ICOMBI .NE. 0 ) THEN
        KLIAB  = KLFREE
        KLFREE = KLFREE + NIDET
      END IF
C?    WRITE (LUPRI,*) ' PSIGN ', PSIGN
C
      IF( ICOMBI .NE. 0 ) THEN
C SET UP ARRAY COMBARING ALPHA AND BETA STRINGS IN IDET LIST
        DO 13 IDET = 1, NIDET
          IAEQIB = 1
          DO 14 IEL = 1, NAEL
            IF(IASTR(IEL,IDET) .NE. IBSTR(IEL,IDET))IAEQIB = 0
   14     CONTINUE
          IWORK(KLIAB-1+IDET) = IAEQIB
   13   CONTINUE
C?      WRITE(LUPRI,*) ' IAEQIB ARRAY : '
C?      CALL IWRTMA(IWORK(KLIAB),1,NIDET,1,NIDET)
      END IF
C
      IF( ISYM .EQ. 0 ) THEN
        LHAMIL = NIDET*NJDET
      ELSE
        LHAMIL = NIDET*(NIDET+1) / 2
      END IF
      CALL SETVEC(HAMIL,0.0D0,LHAMIL)
C
      NTERMS= 0
      NDIF0 = 0
      NDIF1 = 0
      NDIF2 = 0
C.. LOOP OVER J DETERMINANTS
C
      NTEST =   0
      DO 1000 JDET = 1,NJDET
        IF( NTEST .GE. 20 ) WRITE(LUPRI,*) '  ****** 1000 JDET ', JDET
C
C EXPAND JDET
        CALL ISETVC(IWORK(KLJAE),0,NORB)
        CALL ISETVC(IWORK(KLJBE),0,NORB)
C
        IF( ICOMBI .NE. 0 ) THEN
          JAEQJB = 1
          DO 32 IEL = 1, NAEL
            IF(JASTR(IEL,JDET) .NE. JBSTR(IEL,JDET))JAEQJB = 0
   32     CONTINUE
C?        WRITE(LUPRI,*) ' JAEQJB ', JAEQJB
        END IF
C
        DO 40 IAEL = 1, NAEL
          IWORK(KLJAE-1+JASTR(IAEL,JDET) ) = 1
   40   CONTINUE
C
        DO 50 IBEL = 1, NBEL
          IWORK(KLJBE-1+JBSTR(IBEL,JDET) ) = 1
   50   CONTINUE
C
        IF( NTEST .GE. 10 ) THEN
          WRITE(LUPRI,*) ' LOOP 1000 JDET =  ',JDET
C         WRITE(LUPRI,*) ' JASTR AND JBSTR '
C         CALL IWRTMA(JASTR(1,JDET),1,NAEL,1,NAEL)
C         CALL IWRTMA(JBSTR(1,JDET),1,NBEL,1,NBEL)
C         WRITE(LUPRI,*) ' EXPANDED ALPHA AND BETA STRING '
C         CALL IWRTMA(IWORK(KLJAE),1,NORB,1,NORB)
C         CALL IWRTMA(IWORK(KLJBE),1,NORB,1,NORB)
        END IF
C
        IF( ISYM .EQ. 0 ) THEN
          MINI = 1
        ELSE
          MINI = JDET
        END IF
        DO 900 IDET = MINI, NIDET
C?      WRITE(LUPRI,*) '   LOOP 900 IDET .... ',IDET
        IF( ICOMBI .EQ. 0 ) THEN
            NLOOP = 1
        ELSE
          IAEQIB = IWORK(KLIAB-1+IDET)
          IF(IAEQIB+JAEQJB .EQ. 0 ) THEN
            NLOOP = 2
          ELSE
            NLOOP = 1
          END IF
        END IF
C
        DO 899 ILOOP = 1, NLOOP
         NTERMS = NTERMS + 1
C?       WRITE(LUPRI,*) '   899 : ILOOP ' , ILOOP
C
C.. COMPARE DETERMINANTS
C
C SWAP IA AND IB FOR SECOND PART OF COMBINATIONS
        IF( ILOOP .EQ. 2 )
     &  CALL ISWPVE(IASTR(1,IDET),IBSTR(1,IDET),NAEL)
C
        NACM = 0
        DO 61 IAEL = 1, NAEL
          NACM = NACM + IWORK(KLJAE-1+IASTR(IAEL,IDET))
   61   CONTINUE
C
        NBCM = 0
        DO 62 IBEL = 1, NBEL
          NBCM = NBCM + IWORK(KLJBE-1+IBSTR(IBEL,IDET))
   62   CONTINUE
C
        NADIF = NAEL-NACM
        NBDIF = NBEL-NBCM
        IF( NTEST .GE. 10 ) THEN
          WRITE(LUPRI,*) '  LOOP 900 IDET ',IDET
          WRITE(LUPRI,*) ' COMPARISON , NADIF , NBDIF ', NADIF,NBDIF
        END IF
C
        IF(NADIF+NBDIF .GT. 2 ) GOTO 898
C
C
C  FACTOR FOR COMBINATIONS
        IF( ICOMBI .EQ. 0 ) THEN
          CONST = 1.0D0
        ELSE
          IF((JAEQJB +IAEQIB) .EQ.2 ) THEN
            CONST = 1.0D0
          ELSE IF( (JAEQJB+IAEQIB) .EQ. 1 ) THEN
COLD        CONST = 1.0D0/SQRT(2.0D0)*(1.0D0+PSIGN)
            IF (PSIGN .EQ. -1.0D0) GO TO 898
C           ... no contribution between singlet and triplet.
            CONST = SQRT(2.0D0)
          ELSE IF( (JAEQJB+IAEQIB) .EQ. 0 ) THEN
            IF( ILOOP .EQ. 1)  THEN
              CONST = 1.0D0
            ELSE
              CONST = PSIGN
            END IF
          END IF
        END IF
C
C.. FIND DIFFERING ORBITALS AND SIGN FOR PERMUTATION
C
C EXPAND IDET
        CALL ISETVC(IWORK(KLIAE),0,NORB)
        CALL ISETVC(IWORK(KLIBE),0,NORB)
C
          DO 42 IAEL = 1, NAEL
            IWORK(KLIAE-1+IASTR(IAEL,IDET) ) = 1
   42     CONTINUE
C
          DO 52 IBEL = 1, NBEL
            IWORK(KLIBE-1+IBSTR(IBEL,IDET) ) = 1
   52     CONTINUE
C
        IF(NADIF .EQ. 1 ) THEN
          DO 120 IAEL = 1,NAEL
            IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
              IA = IASTR(IAEL,IDET)
              IEL1 = IAEL
              GOTO 121
             END IF
  120     CONTINUE
  121     CONTINUE
C
          DO 130 JAEL = 1,NAEL
            IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
              JA = JASTR(JAEL,JDET)
              JEL1 = JAEL
              GOTO 131
             END IF
  130     CONTINUE
  131     CONTINUE
          SIGNA = (-1)**(JEL1+IEL1)
C?        WRITE(LUPRI,*) ' IA JA SIGNA... ',IA,JA,SIGNA
C
        END IF
C
        IF(NBDIF .EQ. 1 ) THEN
          DO 220 IBEL = 1,NBEL
            IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
              IB = IBSTR(IBEL,IDET)
              IEL1 = IBEL
              GOTO 221
             END IF
  220     CONTINUE
  221     CONTINUE
C
          DO 230 JBEL = 1,NBEL
            IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
              JB = JBSTR(JBEL,JDET)
              JEL1 = JBEL
              GOTO 231
             END IF
  230     CONTINUE
  231     CONTINUE
          SIGNB = (-1)**(JEL1+IEL1)
C?        WRITE(LUPRI,*) ' IB JB SIGNB... ',IB,JB,SIGNB
C
        END IF
        IF(NADIF .EQ. 2 ) THEN
          IDIFF = 0
          DO 320 IAEL = 1,NAEL
            IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
              IF( IDIFF .EQ. 0 ) THEN
                IDIFF = 1
                I1 = IASTR(IAEL,IDET)
                IPERM = IAEL
              ELSE
                I2 = IASTR(IAEL,IDET)
                IPERM = IAEL + IPERM
                GOTO 321
              END IF
            END IF
  320     CONTINUE
  321     CONTINUE
C
          JDIFF = 0
          DO 330 JAEL = 1,NAEL
            IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
              IF( JDIFF .EQ. 0 ) THEN
                JDIFF = 1
                J1 = JASTR(JAEL,JDET)
                JPERM = JAEL
              ELSE
                J2 = JASTR(JAEL,JDET)
                JPERM = JAEL + JPERM
                GOTO 331
              END IF
            END IF
  330     CONTINUE
  331     CONTINUE
          SIGN = (-1)**(IPERM+JPERM)
C
        END IF
C
        IF(NBDIF .EQ. 2 ) THEN
          IDIFF = 0
          DO 420 IBEL = 1,NBEL
            IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
              IF( IDIFF .EQ. 0 ) THEN
                IDIFF = 1
                I1 = IBSTR(IBEL,IDET)
                IPERM = IBEL
              ELSE
                I2 = IBSTR(IBEL,IDET)
                IPERM = IBEL + IPERM
                GOTO 421
               END IF
            END IF
  420     CONTINUE
  421     CONTINUE
C
          JDIFF = 0
          DO 430 JBEL = 1,NBEL
            IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
              IF( JDIFF .EQ. 0 ) THEN
                JDIFF = 1
                J1 = JBSTR(JBEL,JDET)
                JPERM = JBEL
              ELSE
                J2 = JBSTR(JBEL,JDET)
                JPERM = JBEL + JPERM
                GOTO 431
              END IF
            END IF
  430     CONTINUE
  431     CONTINUE
          SIGN = (-1)**(IPERM+JPERM)
C
        END IF
C
C OBTAIN VALUE OF HAMILTONIAN ELEMENT
C
        IF( NADIF .EQ. 2 .OR. NBDIF .EQ. 2 ) THEN
          NDIF2 = NDIF2 + 1
C SIGN * (I1 J1 | I2 J2 ) - ( I1 J2 | I2 J1 )
          XVAL = SIGN*( H2TUVX(I1,J1,I2,J2,FIJKL,ILTSOB)
     *                 -H2TUVX(I1,J2,I2,J1,FIJKL,ILTSOB) )
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 1 ) THEN
          NDIF2 = NDIF2 + 1
C SIGN * (IA JA | IB JB )
C?        WRITE(LUPRI,*) ' IA IB JA JB ', IA,IB,JA,JB
          XVAL = SIGNA*SIGNB* H2TUVX(IA,JA,IB,JB,FIJKL,ILTSOB)
C?        WRITE(LUPRI,*) ' SIGNA SIGNB XVAL ',SIGNA,SIGNB,XVAL
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 0 .OR.
     &           NADIF .EQ. 0 .AND. NBDIF .EQ. 1 )THEN
          NDIF1 = NDIF1 + 1
C SIGN *
C(  H(I1 J1 ) +
C  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 | JORB JORB )
C -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB | JORB J1 ) )
          IF( NADIF .EQ. 1 ) THEN
            I1 = IA
            J1 = JA
            SIGN = SIGNA
          ELSE
            I1 = IB
            J1 = JB
            SIGN = SIGNB
          END IF
C?        WRITE(LUPRI,*) ' ONE DIFF I1 J1 SIGN : ',I1,J1,SIGN
C
          XVAL = ONEBOD(I1,J1)
          DO 520 JAEL = 1, NAEL
            JORB = JASTR(JAEL,JDET)
            XVAL = XVAL + H2TUVX(I1,J1,JORB,JORB,FIJKL,ILTSOB)
  520     CONTINUE
          DO 521 JBEL = 1, NBEL
            JORB = JBSTR(JBEL,JDET)
            XVAL = XVAL + H2TUVX(I1,J1,JORB,JORB,FIJKL,ILTSOB)
  521     CONTINUE
          IF( NADIF .EQ. 1 ) THEN
            DO 522 JAEL = 1, NAEL
              JORB = JASTR(JAEL,JDET)
              XVAL = XVAL - H2TUVX(I1,JORB,JORB,J1,FIJKL,ILTSOB)
  522       CONTINUE
          ELSE
            DO 523 JBEL = 1, NBEL
              JORB = JBSTR(JBEL,JDET)
              XVAL = XVAL - H2TUVX(I1,JORB,JORB,J1,FIJKL,ILTSOB)
  523       CONTINUE
          END IF
          XVAL = XVAL * SIGN
        ELSE IF( NADIF .EQ. 0 .AND. NBDIF .EQ. 0 ) THEN
          NDIF0 = NDIF0 + 1
C SUM(I,J OF JDET) H(I,J) + (I I | J J ) - (I J | J I )
C
          XVAL = ECORE
          DO 650 IAB = 1, 2
            IF(IAB .EQ. 1 ) THEN
              NIABEL = NAEL
            ELSE
              NIABEL = NBEL
            END IF
            DO 640 JAB = 1, 2
              IF(JAB .EQ. 1 ) THEN
                NJABEL = NAEL
              ELSE
                NJABEL = NBEL
              END IF
              DO 630 IEL = 1, NIABEL
                IF( IAB .EQ. 1 ) THEN
                  IORB = IASTR(IEL,IDET)
                ELSE
                  IORB = IBSTR(IEL,IDET)
                END IF
                IF(IAB .EQ. JAB ) XVAL = XVAL + ONEBOD(IORB,IORB)
                IORB = IORB
                DO 620 JEL = 1, NJABEL
                  IF( JAB .EQ. 1 ) THEN
                    JORB = IASTR(JEL,IDET)
                  ELSE
                    JORB = IBSTR(JEL,IDET)
                  END IF
                  XVAL = XVAL + 0.5D0*RJ(IORB,JORB)
                  IF( IAB . EQ. JAB )
     &            XVAL = XVAL - 0.5D0*RK(IORB,JORB)
  620           CONTINUE
  630         CONTINUE
  640       CONTINUE
  650     CONTINUE
        END IF
C
C?      WRITE(LUPRI,*) ' CONST XVAL  ', CONST ,XVAL
        IF( ISYM .EQ. 0 ) THEN
          HAMIL((JDET-1)*NIDET+IDET) =
     &    HAMIL((JDET-1)*NIDET+IDET) + CONST * XVAL
        ELSE
          HAMIL((IDET-1)*IDET/2 + JDET ) =
     &    HAMIL((IDET-1)*IDET/2 + JDET ) + CONST * XVAL
        END IF
C RESTORE ORDER !!!
  898   IF( ILOOP .EQ. 2 )
     &  CALL ISWPVE(IASTR(1,IDET),IBSTR(1,IDET),NAEL)
  899 CONTINUE
  900 CONTINUE
 1000 CONTINUE
C
      IF (NTEST .GE.2) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
     &'  Number of elements differing by 0 excitation.. ',NDIF0
C
         WRITE(LUPRI,*)
     &'  Number of elements differing by 1 excitation.. ',NDIF1
C
         WRITE(LUPRI,*)
     &'  Number of elements differing by 2 excitation.. ',NDIF2
C
         WRITE(LUPRI,*)
     &'  Number of vanishing elments                    ',
     &   NTERMS - NDIF0 - NDIF1 - NDIF2
      END IF
      IF( NTEST .GT. 2 ) THEN
        WRITE(LUPRI,'(/A)') '  HAMILTONIAN MATRIX '
        IF( ISYM .EQ. 0 ) THEN
          CALL WRTMAT(HAMIL,NIDET,NJDET,NIDET,NJDET,0)
        ELSE
          CALL PRSYM(HAMIL,NIDET)
        END IF
      END IF
C
C!    CALL QUIT(' ENFORCED STOP AT END OF DIHDJ ')
      RETURN
      END
C  /* Deck stfdt2 */
      SUBROUTINE STFDT2(NDET,IDET,IASTR2,IBSTR2,MAXSYM,NOCTPA,NOCTPB,
     &                  NSSOA,NSSOB,ISSOA,ISSOB,IOCOC,NAEL,
     &                  NBEL,SYMPRO,IASTR,IBSTR,ISYM,ICOMBI)
C
C EXTRACT ALPHA AND BETASTRINGS FOR NDET DETERMINANTS GIVEN IN
C IDET
C
C RAS VERSION
C
C JEPPE OLSEN , JANUARY 1989
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IDET(NDET),IASTR2(NAEL,NDET),IBSTR2(NBEL,NDET)
C
      DIMENSION NSSOA(NOCTPA,MAXSYM),NSSOB(NOCTPB,MAXSYM)
      DIMENSION ISSOA(NOCTPA,MAXSYM),ISSOB(NOCTPB,MAXSYM)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      INTEGER SYMPRO(8,8)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
C
      NFOUND = 0
      LDET = 0
      DO 1000 IASYM = 1, MAXSYM
        IBSYM = SYMPRO(ISYM,IASYM)
        DO 999 IATP = 1,NOCTPA
          NIA = NSSOA(IATP,IASYM)
          IF (NIA.EQ.0) GO TO 999
          IAOFF = ISSOA(IATP,IASYM)
          DO 901 IBTP = 1, NOCTPB
            IF(IOCOC(IATP,IBTP).LE.0) GOTO 901
            NIB = NSSOB(IBTP,IBSYM)
            IBOFF = ISSOB(IBTP,IBSYM)
            DO 900 IB = 1, NIB
              DO 800 IA = 1, NIA
                IBEFF = IBOFF-1+IB
                IAEFF = IAOFF-1+IA
C?              WRITE(LUPRI,*) ' IAEFF IBEFF ' , IAEFF,IBEFF
C THE FOLLOWING IS NEW LOGICAL ORDERING !
              IF(ICOMBI.NE.0..AND.IAEFF.LT.IBEFF)
     &          GOTO 800
                LDET = LDET + 1
                NEXT  = 0
                DO 500 I =1, NDET
  500           IF(IDET(I) .EQ. LDET ) NEXT = I
                IF(NEXT.NE.0 ) THEN
                  CALL ICOPVE(IASTR(1,IAEFF),IASTR2(1,NEXT),NAEL)
                  CALL ICOPVE(IBSTR(1,IBEFF),IBSTR2(1,NEXT),NBEL)
                  NFOUND = NFOUND + 1
                  IF(NFOUND  .EQ. NDET ) GOTO 1001
                END IF
  800         CONTINUE
  900       CONTINUE
  901     CONTINUE
  999   CONTINUE
 1000 CONTINUE
 1001 CONTINUE
C
      NTEST = 0
      IF( NTEST .GE. 2 ) THEN
         IF( ICOMBI .EQ. 0 ) THEN
           WRITE(LUPRI,*) ' DETERMINANTS SELECTED '
         ELSE
           WRITE(LUPRI,*) ' DETERMINANT COMBINATIONS SELECTED '
         END IF
C
         CALL IWRTMA(IDET,1,NDET,1,NDET)
         WRITE(LUPRI,*) ' SELECTED ALPHA AND BETA STRINGS FROM STRDT1 '
         CALL IWRTMA(IASTR2,NAEL,NDET,NAEL,NDET)
         WRITE(LUPRI,*)
         CALL IWRTMA(IBSTR2,NBEL,NDET,NBEL,NDET)
      END IF
C
      RETURN
      END
C  /* Deck fndmn3 */
      SUBROUTINE FNDMN3(VEC,NDIM,MXELMN,IPLACE,NELMN,IPRT,THRES)
C
C Jeppe Olsen. Last revision 890205/900405 hjaaj.
C
C FIND MXELMN/NELMN LOWEST ELEMENTS IN VEC .
C NELMN IS THE LARGEST NUMBER LOWER THAN MXELMN THAT DOES NOT
C SPLIT DEGENERATE PAIRS
C ORIGINAL PLACES OF THE LOWEST ELEMENTS ARE STORED IN IPLACE
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION VEC(NDIM),IPLACE(*)
C
      PARAMETER ( D1 = 1.0D0 )
C
C. FIRST OCCURANCE OF LOWEST ELEMENT AND LARGEST ELEMENT
      XMIN = VEC(1)
      XMAX = VEC(1)
      IMIN = 1
      IMAX = 1
      DO 100 I = 2,NDIM
        IF( VEC(I) .GT. XMAX) THEN
           XMAX = VEC(I)
           IMAX = I
        ELSE IF( VEC(I) .LT. XMIN) THEN
           XMIN = VEC(I)
           IMIN = I
        END IF
  100 CONTINUE
C
      IF (IPRT.GT.0) WRITE(LUPRI,*) ' >> output from FNDMN3 <<'
      IF(IPRT .GE. 5 ) THEN
         WRITE(LUPRI,*) ' LOWEST VALUE AND PLACE ',XMIN,IMIN
         WRITE(LUPRI,*) ' HIGHST VALUE AND PLACE ',XMAX,IMAX
      END IF
C
      IPLACE(1) = IMIN
      NDEG = 1
      IF(IPRT .GE. 5 ) WRITE(LUPRI,'(A,I8,1P,D15.5,I8,I4)')
     &    'IELMNT XMIN IMIN NDEG ', 1,XMIN,IMIN,NDEG
C
      ITOP = MIN(NDIM,MXELMN+1)
      DO 200 IELMNT = 2,ITOP
        XMINPR = XMIN
        IMINPR = IMIN
        XMIN = XMAX + D1
        IMIN = -1
        DO 150 I = 1,NDIM
          IF (VEC(I) .LT. XMIN) THEN
            IF (VEC(I) .GT. XMINPR) THEN
               IMIN = I
               XMIN = VEC(I)
            ELSE IF (VEC(I) .EQ. XMINPR .AND. I .GT. IMINPR) THEN
               IMIN = I
               XMIN = VEC(I)
               GO TO 151
            END IF
          END IF
  150   CONTINUE
  151   CONTINUE
        IF(XMIN-XMINPR .LE. THRES )THEN
         NDEG = NDEG + 1
        ELSE
         NDEG = 1
        END IF
C
        IF( IELMNT .LE. MXELMN ) IPLACE(IELMNT) = IMIN
        IF( IPRT .GE. 5 ) WRITE(LUPRI,'(A,I8,1P,D15.5,I8,I4)')
     &    'IELMNT XMIN IMIN NDEG ', IELMNT,XMIN,IMIN,NDEG
  200 CONTINUE
C
C CHECK DEGENERACY ON LAST VALUE
      IF(MXELMN .LT. NDIM ) THEN
        NELMN = MXELMN+1-NDEG
      ELSE
        NELMN = NDIM
      END IF
      IF (IPRT .GT. 0)
     &   WRITE(LUPRI,*) ' NUMBER OF ELEMENTS OBTAINED IN FNDMN3 ',NELMN
C
C
      IF( IPRT  .GE. 3 ) THEN
        WRITE(LUPRI,*) ' FROM FNDMN3 : '
        WRITE(LUPRI,*) '   PLACES OF LOWEST ELEMENTS '
        CALL IWRTMA(IPLACE,1,NELMN ,1,NELMN )
      END IF
C
      IF( IPRT .GE. 1 ) THEN
       WRITE(LUPRI,*)
     & ' MIN AND MAX IN SELECTED SUPSPACE ',
     &   VEC(IPLACE(1)),VEC(IPLACE(NELMN))
      END IF
C
      RETURN
      END
C  /* Deck fndamx */
      SUBROUTINE FNDAMX(VEC,NDIM,MXELMN,IPLACE,NELMN,IPRT,THRES)
C
C 890205 hjaaj. Based on FNDMN3 by Jeppe Olsen.
C Last revision 900405 hjaaj.
C
C FIND MXELMN/NELMN ELEMENTS IN VEC with highest absolute value.
C NELMN IS THE LARGEST NUMBER LOWER THAN MXELMN THAT DOES NOT
C SPLIT DEGENERATE PAIRS.
C ORIGINAL PLACES OF THE LOWEST ELEMENTS ARE STORED IN IPLACE.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION VEC(NDIM),IPLACE(*)
C
      PARAMETER ( D1 = 1.0D0 )
C
C. FIRST OCCURANCE OF LOWEST ELEMENT AND LARGEST ELEMENT
      XMIN = ABS(VEC(1))
      XMAX = ABS(VEC(1))
      IMIN = 1
      IMAX = 1
      DO 100 I = 2,NDIM
        AVECI = ABS(VEC(I))
        IF( AVECI .GT. XMAX) THEN
           XMAX = AVECI
           IMAX = I
        ELSE IF( AVECI .LT. XMIN) THEN
           XMIN = AVECI
           IMIN = I
        END IF
  100 CONTINUE
C
      IF (IPRT.GT.0) WRITE(LUPRI,*) ' >> output from FNDAMX <<'
      IF(IPRT .GE. 5 ) THEN
         WRITE(LUPRI,*) ' LOWEST VALUE AND PLACE ',XMIN,IMIN
         WRITE(LUPRI,*) ' HIGHST VALUE AND PLACE ',XMAX,IMAX
      END IF
C
      IPLACE(1) = IMAX
      NDEG = 1
      IF(IPRT .GE. 5 ) WRITE(LUPRI,'(A,I8,1P,D15.5,I8,I4)')
     &    'IELMNT value IMAX NDEG ', 1,VEC(IMAX),IMAX,NDEG
C
      ITOP = MIN(NDIM,MXELMN+1)
      DO 200 IELMNT = 2,ITOP
        XMAXPR = XMAX
        IMAXPR = IMAX
        XMAX   = -D1
        IMAX   = -1
        DO 150 I = 1,NDIM
          AVECI = ABS(VEC(I))
          IF (AVECI .GT. XMAX) THEN
            IF (AVECI .LT. XMAXPR) THEN
               IMAX = I
               XMAX = AVECI
            ELSE IF (AVECI .EQ. XMAXPR .AND. I .GT. IMAXPR) THEN
               IMAX = I
               XMAX = AVECI
               GO TO 151
            END IF
          END IF
  150   CONTINUE
  151   CONTINUE
        IF(XMAXPR-XMAX .LE. THRES )THEN
          NDEG = NDEG + 1
        ELSE
          NDEG = 1
        END IF
C
        IF( IELMNT .LE. MXELMN ) IPLACE(IELMNT) = IMAX
        IF( IPRT .GE. 5 ) WRITE(LUPRI,'(A,I8,1P,D15.5,I8,I4)')
     &    ' IELMNT value IMAX NDEG ',IELMNT,VEC(IMAX),IMAX,NDEG
  200 CONTINUE
C
C CHECK DEGENERACY ON LAST VALUE
      IF(MXELMN .LT. NDIM ) THEN
        NELMN = MXELMN+1-NDEG
      ELSE
        NELMN = NDIM
      END IF
      IF (IPRT .GT. 0)
     &   WRITE(LUPRI,*) ' NUMBER OF ELEMENTS OBTAINED IN FNDAMX ',NELMN
C
C
      IF( IPRT  .GE. 3 ) THEN
        WRITE(LUPRI,*) ' FROM FNDAMX : '
        WRITE(LUPRI,*) '   PLACES OF ELEMENTS of highest absolut value'
        CALL IWRTMA(IPLACE,1,NELMN ,1,NELMN )
      END IF
C
      IF( IPRT .GE. 1 ) THEN
       WRITE(LUPRI,*)
     & ' MAX AND MIN IN SELECTED SUPSPACE ',
     &   VEC(IPLACE(1)),VEC(IPLACE(NELMN))
      END IF
C
      RETURN
      END
C  /* Deck tripk2 */
      SUBROUTINE TRIPK2(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
C
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C              LOWER HALF OF AUTPAK IS STORED IN APAK
C   IWAY = 2 : PACKED TO FULL FORM
C              APAK STORED IN LOWER HALF
C               SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
C.. NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
#include <implicit.h>
#include <priunit.h>
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 J = 1,NDIM
          DO 50  I = J , NDIM
           APAK(IJ+I) = AUTPAK(I,J)
   50     CONTINUE
          IJ = IJ +NDIM-J
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 J = 1,NDIM
          DO 150  I = J,NDIM
           AUTPAK(I,J) = APAK(IJ+I)
           AUTPAK(J,I) = SIGN*APAK(IJ+I)
  150     CONTINUE
          IJ = IJ + NDIM-J
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' AUTPAK AND APAK FROM TRIPK2 '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM,0)
        CALL PRSM2(APAK,NDIM)
      END IF
C
      RETURN
      END
C  /* Deck tribl3 */
      SUBROUTINE TRIBL3(AI,AO,NRC,NCC,LRC,LCC,IBLKI,IWAY,
     &                  IBLKO,IDOBLK,SIGN)
C
C TRANSFER BETWEEN LOWER PACKED FORM AND COMPLETE FORM OF
C A SYMMETRIC OR ANTISYMMETRIC BLOCKED MATRIX .
C
C IWAY = 1 : UNPACKED TO PACKED FORM
C            LOWER HALF OF COMPLETE MATRIX IS PACKED
C
C IWAY = 2 : PACKED TO UNPACKED FORM
C            LOWER HALF OF COMPLETE MATRIX EQUAL PACKED FORM
C            UPPER HALF IS MULTIPLIED WITH SIGN
C IF IDOBLK DIFFERS FROM ZERO IBLK MATRIX IS CONSTRUCTED FOR
C OUTPUT FORM
C
C
C ... JEPPE OLSEN JANUARY 1989
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION AI(*),AO(*)
      DIMENSION LRC(NRC),LCC(NCC)
      DIMENSION IBLKI(NRC,NCC),IBLKO(NRC,NCC)
C
      NTEST = 0
C?    WRITE(LUPRI,*) ' ENTERING TRIBL3 '
      IF( IWAY .EQ. 1 ) THEN
        IF(IDOBLK.NE. 0 ) CALL ISETVC(IBLKO,0,NRC*NCC)
        IOFFO = 1
        DO 10 IRC = 1, NRC
          DO 5 ICC = 1, IRC
C?          WRITE(LUPRI,*) ' IRC ICC IBKLI ', IRC,ICC,IBLKI(IRC,ICC)
            IF(IBLKI(IRC,ICC).GT. 0 ) THEN
              IF(IDOBLK .NE. 0 ) IBLKO(IRC,ICC) = IOFFO
              IF(IRC.EQ.ICC) THEN
                CALL TRIPK2(AI(IBLKI(IRC,ICC)),AO(IOFFO),1,
     &                      LRC(IRC),LCC(ICC),1.0D0)
                IOFFO = IOFFO + LRC(IRC)*(LRC(IRC)+1) / 2
              ELSE
                CALL COPVEC(AI(IBLKI(IRC,ICC)),AO(IOFFO),
     &                      LRC(IRC)*LCC(ICC) )
                IOFFO = IOFFO + LRC(IRC)*LCC(ICC)
              END IF
            END IF
    5     CONTINUE
   10   CONTINUE
C
        IF( NTEST .NE. 0 ) THEN
          WRITE(LUPRI,*)
     &    '   TRIBL3 : UNPACKED MATRIX(INPUT), PACKED MATRIX(OUTPUT)'
          CALL PRBL3(AI,NRC,LRC,NCC,LCC,IBLKI)
          WRITE(LUPRI,*)
C              PRSBL3(A,NRC,LRC,NCC,LCC,IBLK)
          CALL PRSBL3(AO,NRC,LRC,NCC,LCC,IBLKO)
         END IF
      ELSE IF ( IWAY .EQ. 2 ) THEN
        IOFFO = 1
        DO 20 IRC = 1, NRC
          DO 15 ICC = 1, NCC
            IF(IRC.GE.ICC .AND. IBLKI(IRC,ICC) .GT. 0 .OR.
     &         IRC.LT.ICC .AND. IBLKI(ICC,IRC) .GT. 0 ) THEN
               IF(IDOBLK.NE.0 ) IBLKO(IRC,ICC) = IOFFO
               IF(        IRC.GT.ICC) THEN
                 IOFFI = IBLKI(IRC,ICC)
                 CALL COPVEC(AI(IOFFI),AO(IOFFO),LRC(IRC)*LCC(ICC))
               ELSE IF ( IRC .EQ. ICC ) THEN
                 IOFFI = IBLKI(IRC,ICC)
                 CALL TRIPK2(AO(IOFFO),AI(IOFFI),2,LRC(IRC),LRC(IRC),
     &                       SIGN)
               ELSE IF ( IRC .LT. ICC) THEN
                 IOFFI = IBLKI(ICC,IRC)
                 CALL TRPMAT(AI(IOFFI),LCC(ICC),LRC(IRC),AO(IOFFO) )
                 IF(SIGN .EQ. -1.0D0)
     &           CALL SCALVE(AO(IOFFO),SIGN,LRC(IRC)*LCC(ICC) )
               END IF
               IOFFO = IOFFO + LRC(IRC)*LCC(ICC)
            END IF
   15     CONTINUE
   20   CONTINUE
        IF( NTEST .NE. 0 ) THEN
          WRITE(LUPRI,*)
     &    '   TRIBL3 : UNPACKED MATRIX(OUTPUT), PACKED MATRIX(INPUT)'
          CALL PRBL3(AO,NRC,LRC,NCC,LCC,IBLKO)
          WRITE(LUPRI,*)
C              PRSBL3(A,NRC,LRC,NCC,LCC,IBLK)
          CALL PRSBL3(AI,NRC,LRC,NCC,LCC,IBLKI)
         END IF
      END IF
C
      RETURN
      END
C  /* Deck cdiag */
      SUBROUTINE CDIAG(AIN,AOUT,NR,NC,FACTOR,IFLAG,ISYM)
C
C MODIFY DIAGONAL ELEMENTS OF A MATRIX
C
C IFLAG = 1 : MULTIPLY DIAGONAL ELEMENTS WITH FACTOR
C IFLAG = 2 : ADD      DIAGONAL ELEMENTS TO   FACTOR
C
C ISYM = 0 : DIAGONAL BLOCKS ARE STORED IN COMPLETE FORM
C ISYM  .GT. 0 : DIAGONAL BLOCKS ARE STORED IN  PACKED FORM
C                ( LOWER TRIANGULAR FORM;COLUMNWISE PACKED !! )
#include <implicit.h>
      DIMENSION AIN(*),AOUT(*)
C
      IF( ISYM .EQ. 0 ) THEN
        IMAX = MIN(NR,NC)
        CALL COPVEC(AIN,AOUT,NR*NC)
        IF( IFLAG .EQ. 1 ) THEN
          DO 100 I = 1, IMAX
            AOUT((I-1)*NC+I) = AIN((I-1)*NC+I) * FACTOR
  100     CONTINUE
        ELSE
          DO 200 I = 1, IMAX
            AOUT((I-1)*NC+I) = AIN((I-1)*NC+I) + FACTOR
  200     CONTINUE
        END IF
        ELSE IF ( ISYM .GT. 0 ) THEN
C NR = NC FOR SYMMETRIC PACKING
        NDIM = NR
        CALL COPVEC(AIN,AOUT,NDIM*(NDIM+1)/2)
        IF( IFLAG .EQ. 1 ) THEN
          DO 1100 I = 1, NDIM
            II = (I-1)*NDIM -I*(I-1)/2 + I
            AOUT(II) = AIN(II) * FACTOR
 1100     CONTINUE
        ELSE
          DO 1200 I = 1, NDIM
            II = (I-1)*NDIM -I*(I-1)/2 + I
            AOUT(II) = AIN(II) + FACTOR
 1200     CONTINUE
        END IF
      END IF
C
      RETURN
      END
C  /* Deck cdibl3 */
      SUBROUTINE CDIBL3(AI,NRC,NCC,LRC,LCC,IBLK,
     &                  AO,FACTOR,IFLAG,ISYM)
C
C CHANGE DIAGONAL ELEMENTS OF A BLOCKED MATRIX
C IFLAG = 1 : MULTIPLY DIAGONAL ELEMENTS WITH FACTOR
C IFLAG = 2 : ADD    DIAGONAL ELEMENTS WITH FACTOR
C
C ISYM= 0 : DIAGONAL BLOCKS STORED IN COMPLETE FORM
C ISYM.GT.0 : DIAGONAL BLOCKS ARE PACKED IN LOWER TRIANGULAR FORM
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION AI(*),AO(*)
      DIMENSION LRC(NRC),LCC(NCC)
      DIMENSION IBLK(NRC,NCC)
C
      DO 100 IR = 1, NRC
        DO 100 IC = 1, NCC
          IF(IBLK(IR,IC) .GT. 0 ) THEN
            IF( IR .NE. IC ) THEN
              CALL COPVEC(AI(IBLK(IR,IC)),AO(IBLK(IR,IC)),
     &                    LRC(IR)*LCC(IC) )
            ELSE
C                  CDIAG(AIN,AOUT,NR,NC,FACTOR,IFLAG)
              CALL CDIAG(AI(IBLK(IR,IC)),AO(IBLK(IR,IC)),
     &                   LRC(IR),LCC(IC),FACTOR,IFLAG,ISYM)
            END IF
          END IF
  100 CONTINUE
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        IF( IFLAG.EQ.1 ) WRITE(LUPRI,*)
     &  '  DIAGONAL ELEMENTS MULTIPLIES WITH ', FACTOR
        IF( IFLAG.EQ.2 ) WRITE(LUPRI,*)
     &  '  DIAGONAL ELEMENTS ADDED TO ', FACTOR
        WRITE(LUPRI,*) ' CDIABL3 : INPUT AND OUTPUT MATRICES '
C                 PRBL3(A,NRC,LRC,NCC,LCC,IBLK)
        IF( ISYM .EQ. 0 ) THEN
          CALL PRBL3(AI,NRC,LRC,NCC,LCC,IBLK)
          WRITE(LUPRI,*)
          CALL PRBL3(AO,NRC,LRC,NCC,LCC,IBLK)
        ELSEIF( ISYM .GT. 0 ) THEN
          CALL PRBL3(AI,NRC,LRC,NCC,LCC,IBLK)
          WRITE(LUPRI,*)
          CALL PRBL3(AO,NRC,LRC,NCC,LCC,IBLK)
        END IF
      END IF
C
      RETURN
      END
#ifdef UNDEF
/* Comdeck blk_comments */
C
C VERY OFTEN QUANTUM CHEMISTRY CODES INVOLVED A LOT
C OF MANIPULATION OF BLOCKED MATRICES.BLOCKING OF MATRICES
C CAN FOR EXAMPLE RESULT FROM SYMMETRY DIVISION OR DIVISION
C INTO DIFFERENT OCCUPATION TYPES.IT WOULD BE CONVENIENT TO
C HAVE A RATHER GENERAL SET OF ROUTINES FOR BLOCKED MATRICES.
C THE ROUTINES IN THIS FILE IS A BEGINNING TO A START  OF A
C FIRST ATTEMPT ON DEVELOPING THIS SUITE OF CODES .
C
C A MATRIX IS BLOCKED BY DIVIDING THE ROWS AND COLUMNS INTO
C SEVERAL CLASSES.ALL ELEMENTS OF A CLASS ARE CONSECUTIVE ELEMENTS.
C A BLOCKING IS DEFINED BY
C ========================
C
C ROWS :    NRC    :  NUMBER OF ROW CLASSES
C           IBRC(I): ABSOLUTE PLACE OF FIRST ELEMNT IN ROW CLASS I
C           LRC(I) : NUMBER OF ELEMENTS IN ROW CLASS I
C
C COLUMNS : NCC    : NUMBER OF COLUMN CLASSES
C           IBCC(I): ABSOLUTE PLACE OF FIRST ELEMNT IN COLUMN CLASS I
C           LCC(I) : NUMBER OF ELEMENTS IN COLUMN  CLASS I
C
C A BLOCKED MATRIX IS DEFINED BY A MATRIX IBLK(NRC,NCC),
C============================
C            IBLK(IR,IC) GT 0 : BLOCK IR,IC IS INCLUDED.
C                               FIRST ELEMNTS OF BLOCK IS IBLK(IR,IC)
C            ELSE             : BLOCK IR,IC IS NOT INCLUDED
C
C A BLOCKED MATRIX IS STORED AS
C ==============================
C
C     LOOP OVER ROW CLASSES
C       LOOP OVER COLUMN CLASSES ALLOWED FOR THE ROW CLASS
C         LOOP OVER COLUMNS OF BLOCK
C             LOOP OVER ROWS OF BLOCK
C             END OVER LOOP OVER ROWS OF BLOCK
C         END OF LOOP OVER COLUMNS OF BLOCK
C       END OVER LOOPS OF ALLOWED COLUMN BLOCKS FOR GIVEN ROW BLOCK
C     END OF LOOP OVER ROW BLOCKS
C
C
C DEFINITION OF ATTRIBUTES FOR CLASSES AND ELEMENTS WILL SOON FOLLOW
C
C SUBROUTINES WRITTEN
C PRBL3 : PRINT BLOCKED MATRIX A
C PRBL3(A,NRC,LRC,NCC,LCC,IBLK)
C
#endif
C  /* Deck prbl3 */
      SUBROUTINE PRBL3(A,NRC,LRC,NCC,LCC,IBLK)
C
C PRINT BLOCKED MATRIX
C
C JANUARY 1989 , JEPPE OLSEN
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION A(*)
      DIMENSION LRC(NRC),LCC(NCC),IBLK(NRC,NCC)
C
C     WRITE(LUPRI,*) ' ENTERING PRBLM3 : '
      DO 200 IRC = 1, NRC
        DO 100 ICC = 1, NCC
C           WRITE(LUPRI,*) ' IRC ICC IBLK ', IRC,ICC,IBLK(IRC,IRC)
          IF( IBLK(IRC,ICC).GT.0 .AND. LRC(IRC)*LCC(ICC).NE.0) THEN
            WRITE(LUPRI,'(A,2I3)')
     &      '   BLOCK...',IRC,ICC
            WRITE(LUPRI,'(A,2I3)')
     &      '  =================='
            WRITE(LUPRI,*)
            IPNTR = IBLK(IRC,ICC)
            CALL WRTMAT(A(IPNTR),LRC(IRC),LCC(ICC),LRC(IRC),LCC(ICC),0)
            WRITE(LUPRI,*)
          END IF
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck tpbl3 */
      SUBROUTINE TPBL3(AI,NRCI,NCCI,LRCI,LCCI,IBRCI,IBCCI,IBLKI,
     &                 AO,NRCO,NCCO,LRCO,LCCO,IBRCO,IBCCO,IBLKO )
C
C TRANSPOSE BLOCKED MATRIX AI AND STORE IN AO
C ARRAYS FOR TRANSPOSED BLOCKING IS ALSE OBTAINED.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION AI(*),AO(*)
      DIMENSION LRCI(NRCI),LCCI(NCCI),IBRCI(NRCI),IBCCI(NCCI)
      DIMENSION LRCO(NCCI),LCCO(NRCI),IBRCO(NCCI),IBCCO(NRCI)
      DIMENSION IBLKI(NRCI,NCCI),IBLKO(NCCI,NRCI)
C
C SET UP ARRAYS FOR TRANSPOSED BLOCKING
C
      NRCO = NCCI
      NCCO = NRCI
      WRITE(LUPRI,*) ' NRCO AND NCCO ', NRCO,NCCO
C
      CALL ICOPVE(LCCI,LRCO,NRCO)
      CALL ICOPVE(LRCI,LCCO,NCCO)
      CALL ICOPVE(IBCCI,IBRCO,NRCO)
      CALL ICOPVE(IBRCI,IBCCO,NCCO)
      WRITE(LUPRI,*) ' LRCO LCCO IBRCO IBCCO '
      CALL IWRTMA(LRCO,1,NRCO,1,NRCO)
      CALL IWRTMA(LCCO,1,NCCO,1,NCCO)
      CALL IWRTMA(IBRCO,1,NRCO,1,NRCO)
      CALL IWRTMA(IBCCO,1,NCCO,1,NCCO)
C
C TRANSPOSE FILE AND GENERATE IBLKO
C
      CALL ISETVC(IBLKO,0,NRCO*NCCO)
      IPNTR = 1
      DO 200 IRO = 1, NRCO
        DO 100 ICO = 1, NCCO
          IF(IBLKI(ICO,IRO) .GT. 0 ) THEN
            CALL TRPMAT(AI(IBLKI(ICO,IRO)),LRCO(ICO),LCCO(IRO),
     &      AO(IPNTR))
            IBLKO(IRO,ICO) = IPNTR
            IPNTR = IPNTR + LRCO(IRO)*LCCO(ICO)
          END IF
  100   CONTINUE
  200 CONTINUE
C
      NTEST = 1
      IF( NTEST .GE. 1 ) THEN
        WRITE(LUPRI,*) ' INFO FROM TPBL3 : '
        WRITE(LUPRI,*) ' IBKLO ARRAY '
        CALL IWRTMA(IBLKO,NRCO,NCCO,NRCO,NCCO)
        WRITE(LUPRI,*) ' INPUT MATRIX '
        CALL PRBL3(AI,NRCI,LRCI,NCCI,LCCI,IBLKI)
        WRITE(LUPRI,*) ' OUTPUT MATRIX '
        CALL PRBL3(AO,NRCO,LRCO,NCCO,LCCO,IBLKO)
C PRBL3(A,NRC,LRC,NCC,LRC,IBLK)
      END IF
C
      RETURN
      END
      FUNCTION H2TUVX(I,J,K,L,H2AC,ILTSOB)
C
C 29-Jan-1989 hjaaj
C l.r. 18-apr-1990 hjaaj
C
C get (tu|vx) integral where t,u,v,x are in Lunar order
C
#include <implicit.h>
      DIMENSION H2AC(NNASHX,NNASHX), ILTSOB(NASHT)
#include <iratdef.h>
C
C Used from common blocks:
C   INFORB : NASHT,NNASHX
C   INFTAP : LUH2AC
C   CBGETDIS : IADH2
C
#include <inforb.h>
#include <inftap.h>
#include <cbgetdis.h>
C
C     Obtain index in Sirius order
C
      I1 = ILTSOB(I)
      J1 = ILTSOB(J)
      K1 = ILTSOB(K)
      L1 = ILTSOB(L)
      IF (I1 .LT. J1) THEN
         ISWP = I1
         I1   = J1
         J1   = ISWP
      END IF
      IF (K1 .LT. L1) THEN
         ISWP = K1
         K1   = L1
         L1   = ISWP
      END IF
      IJ1 = (I1*I1 - I1)/2 + J1
      KL1 = (K1*K1 - K1)/2 + L1
C
C     Now we can pick up the desired integral
C
      IF (IADH2 .GE. 0) THEN
         CALL READDX(LUH2AC,IADH2+KL1,IRAT*IJ1,H2AC)
         KL1 = 1
      END IF
      H2TUVX = H2AC(IJ1,KL1)
      RETURN
      END
C  /* Deck eigen */
      SUBROUTINE EIGEN(A,R,N,MV,MFKR)
#include <implicit.h>
#include <priunit.h>
      DIMENSION A(*),R(*)
      DATA TESTIT/1.D-20/
      DATA TESTX /1.D-26/
      DATA TESTY /1.D-18/
C
C        PURPOSE
C           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C           MATRIX
C
C        USAGE
C           CALL EIGEN(A,R,N,MV,MFKR)
C
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN ASSCENDING ORDER.
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C       SEQUENCE)
C           MFKR=0 NO SORT
C               =1 SORT
C
C        REMARKS
C           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)
C           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED
C           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN "MATHEMATICAL
C           METHODS FOR DIGITAL COMPUTERS", EDITED BY A. RALSTON AND
C           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7
C
C     ..................................................................
C
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
C     DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
C    1 COSX2,SINCS,RANGE
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT
C        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C        900111-hjaaj: use generic SQRT and ABS
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
      RANGE=1.0D-12
      IF(MV.eq.1) GO TO 25
      IQ=-N
      DO 20 J=1,N
        IQ=IQ+N
      DO 20 I=1,N
        IJ=IQ+I
        R(IJ)=0.0D+00
        IF(I.eq.J) R(IJ)=1.0D+00
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0D+00
      DO 35 J=2,N
      DO 35 I=1,J-1
        IA=I+(J*J-J)/2
        ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(ANORM .LE. 0.0D0) GO TO 165
      ANORM=SQRT(2.0D0*ANORM)
      ANRMX=ANORM*RANGE/DFLOAT(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/DFLOAT(N)
      IF(THR.LT.TESTY)THR=0.D0
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
      IF(ABS(A(LM)).LT.TESTY)A(LM)=0.D0
      IF(ABS(A(LM)).EQ.0.D0.AND.THR.EQ.0.D0)GO TO 130
      IF(ABS(A(LM)) .lt. THR) go to 130
      IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5D+00*(A(LL)-A(MM))
      AJUK=(A(LM)*A(LM)+X*X)
      AJUK=SQRT(AJUK)
      IF(ABS(AJUK).LT.TESTIT) THEN
        WRITE(LUPRI,3000)TESTIT,AJUK,A(LM)
 3000 FORMAT(' ***in EIGEN: DENOMINATOR LT ',D14.6,'. VALUE=',D14.6,
     &'. NUMERATOR=',D14.6)
        Y=0.D0
      ELSE
        Y=-A(LM)/AJUK
      END IF
C     Y=-A(LM)/ SQRT(A(LM)*A(LM)+X*X)
      IF(X .lt. 0.0d0) Y=-Y
      AJUK=(1.D0-Y*Y)
      IF(AJUK.LT.0.D0) THEN
        WRITE(LUPRI,3001) AJUK
 3001 FORMAT(' ***in EIGEN: SQRT OF ',D20.8)
        AJUK=0.D0
      ELSE
        AJUK=SQRT(AJUK)
      END IF
      AJUK=2.D0*(1.D0+AJUK)
      AJUK=SQRT(AJUK)
      SINX=Y/AJUK
C     SINX=Y/ SQRT(2.0D+00*(1.0D+00+( SQRT(1.0D+00-Y*Y))))
      SINX2=SINX*SINX
C     COSX= SQRT(1.0D+00-SINX2)
      AJUK=1.D0-SINX2
      IF(AJUK.LT.TESTX)AJUK=0.D0
      COSX=SQRT(AJUK)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF(I.eq.L .or. I.eq.M) go to 115
      IF(I.lt.M) THEN
        IM=I+MQ
      ELSE
        IM=M+IQ
      END IF
      IF(I.lt.L) THEN
        IL=I+LQ
      ELSE
        IL=L+IQ
      END IF
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF(MV.eq.1) go to 125
      ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0D+00*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 IF (M.ne.N) THEN
        M=M+1
        GO TO 60
      END IF
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
      IF (L.ne.(N-1)) THEN
        L=L+1
        GO TO 55
      END IF
      IF (IND.eq.1) THEN
        IND=0
        GO TO 50
      END IF
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
      IF(THR .gt. ANRMX) go to 45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      IF(MFKR.EQ.0)GO TO 186
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(MM).ge.A(LL)) go to 185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV.eq.1) go to 185
      DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
186   CONTINUE
      RETURN
      END
