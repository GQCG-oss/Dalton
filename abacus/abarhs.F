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
!===========================================================================
!  FILE : abarhs.F
!
!  rhs = "Right Hand Side" of linear response equations for molecular
!        Hessian and more.
!===========================================================================
C  /* Deck rhsinp */
      SUBROUTINE RHSINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      LOGICAL NEWDEF
      PARAMETER (NTABLE = 42)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION NATOM(10)
#include <abainf.h>
#include <cbirhs.h>
#include <inftap.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', 'XXXXXXX',
     &            '.INTSKI', '.INTPRI', '.RETURN',
     &            '.FCKSKI', '.FCKPRI', 'XXXXXXX',
     &            '.NOH1  ', '.NOH2  ', '.FCKTES',
     &            '.SDRSKI', '.TRASKI', '.SDRPRI',
     &            '.TRAPRI', '.SDRTES', 'XXXXXXX',
     &            '.TRATES', '.GDYSKI', '.GDYPRI',
     &            '.NODC  ', '.NODV  ', '.NOPV  ',
     &            '.FSTTES', 'XXXXXXX', '.NOFD  ',
     &            '.NOFS  ', '.NOSSF ', '.SIRPR4',
     &            '.SIRPR6', '.NOORTH', '.GDHAM ',
     &            '.ALLCOM', '.TIME  ', 'XXXXXXX',
     &            '.STOP  ', '.NODDY ', '.NODPTR',
     &            '.PTRPRI', '.SORPRI', '.PTRSKI'/
C
      NEWDEF = WORD .EQ. '*GETSGY'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     ICHANG = ICHANG + 1
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     *                      15,16,17,18,19,20,21,22,23,24,25,
     *                      26,27,28,29,30,31,32,33,34,35,36,
     *                      37,38,39,40,41,42), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in RHSINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in RHSINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, *) IPRALL
                  IF (IPRALL .EQ. IPRDEF) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     IF (IPRINT .EQ. IPRDEF) IPRINT = IPRALL
                     IF (IPRFD  .EQ. IPRDEF) IPRFD  = IPRALL
                     IF (IPRSDR .EQ. IPRDEF) IPRSDR = IPRALL
                     IF (IPRTRA .EQ. IPRDEF) IPRTRA = IPRALL
                     IF (IPRGDY .EQ. IPRDEF) IPRGDY = IPRALL
                  END IF
               GO TO 100
    3             CONTINUE
               GO TO 100
    4             RUNINT = .FALSE.
               GO TO 100
    5             READ (LUCMD, *) IPRINT, IPRNTA, IPRNTB,
     *                                      IPRNTC, IPRNTD
                  IPRSUM = IPRNTA + IPRNTB + IPRNTC + IPRNTD
                  IF (IPRINT .EQ. IPRDEF .AND. IPRSUM .EQ. 0) THEN
                     ICHANG = ICHANG - 1
                  END IF
               GO TO 100
    6             CONTINUE
                  RETUR = .TRUE.
               GO TO 100
    7             RUNFCK = .FALSE.
               GO TO 100
    8             READ (LUCMD,*) IPRFD
                  IF (IPRFD .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    9             CONTINUE
               GO TO 100
   10             NOH1 = .TRUE.
               GO TO 100
   11             NOH2 = .TRUE.
               GO TO 100
   12             FCKTST = .TRUE.
               GO TO 100
   13             RUNSDR = .FALSE.
               GO TO 100
   14             RUNTRA = .FALSE.
               GO TO 100
   15             READ (LUCMD,*) IPRSDR
                  IF (IPRSDR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   16             READ (LUCMD,*) IPRTRA
                  IF (IPRTRA .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   17             SDRTST = .TRUE.
               GO TO 100
   18             CONTINUE
               GO TO 100
   19             TRATST = .TRUE.
               GO TO 100
   20             RUNGDY = .FALSE.
               GO TO 100
   21             READ (LUCMD,*) IPRGDY
                  IF (IPRGDY .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   22             NODC = .TRUE.
               GO TO 100
   23             NODV = .TRUE.
               GO TO 100
   24             NOPV = .TRUE.
               GO TO 100
   25             TESTFS = .TRUE.
               GO TO 100
   26             CONTINUE
               GO TO 100
   27             NOFD = .TRUE.
               GO TO 100
   28             NOFS = .TRUE.
               GO TO 100
   29             NOSSF = .TRUE.
               GO TO 100
   30             READ (LUCMD,*) IPRIN4
                  IF (IPRIN4 .EQ. 0) ICHANG = ICHANG - 1
               GO TO 100
   31             READ (LUCMD,*) IPRIN6
                  IF (IPRIN6 .EQ. 0) ICHANG = ICHANG - 1
               GO TO 100
   32             NOORTH = .TRUE.
               GO TO 100
   33             GDHAM  = .TRUE.
               GO TO 100
   34             RALLCO = .TRUE.
               GO TO 100
   35             TKTIME  = .TRUE.
               GO TO 100
   36             CONTINUE
               GO TO 100
   37             CUT     = .TRUE.
               GO TO 100
   38             NODDY   = .TRUE.
               GO TO 100
   39             NODPTR  = .TRUE.
               GO TO 100
   40             READ (LUCMD,*) IPRPTR
               GO TO 100
   41             READ (LUCMD,*) IPRSOR
               GO TO 100
   42             RUNPTR = .FALSE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in RHSINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in RHSINP.')
            END IF
      END IF
  300 CONTINUE
C
      IPRMAX = MAX(IPRALL,IPRINT,IPRFD,IPRSDR,IPRTRA,IPRGDY)
C
      IF (NEWDEF .AND. (ICHANG .GT. 0) .AND. SKIP) THEN
         CALL HEADER('Changes of defaults for RHSIDE:',0)
         WRITE (LUPRI,'(A)') ' RHSIDE skipped in this run.'
      ELSE IF (NEWDEF .AND. ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for RHSIDE:',0)
         IF (IPRALL .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     *         ' Print level in RHSIDE:',IPRALL
         IF (TKTIME) WRITE (LUPRI,'(A)')
     *         ' Detailed timing for integral calculation '//
     *          'will be provided.'
         IF (NODC) WRITE (LUPRI,'(A)')
     *         ' DC will be neglected in this run.'
         IF (NODV) WRITE (LUPRI,'(A)')
     *         ' DV will be neglected in this run.'
         IF (NOPV) WRITE (LUPRI,'(A)')
     *         ' PV will be neglected in this run.'
         IF (UNDIF) WRITE (LUPRI,'(A)')
     *         ' Undifferentiated integrals will be calculated.'
         IF (NOH1) WRITE (LUPRI,'(A)')
     *         ' One-electron Hamiltonian neglected.'
         IF (NOH2) WRITE (LUPRI,'(A)')
     *         ' Two-electron Hamiltonian neglected.'
         IF (NOORTH) WRITE (LUPRI,'(A)')
     *         ' No orthonormalization in differentiated gradients.'
         IF (GDHAM)  WRITE (LUPRI,'(A)')
     *         ' File LUINTD for GD hamiltonian written.'
         IF (RALLCO) WRITE (LUPRI,'(A)')
     *         ' All components of right-hand sides calculated '//
     *         'simultaneously.'
         IF (MAXSIM .NE. 3) WRITE (LUPRI,'(A,I3)')
     *         ' Maximum number of simultaneous directions:',MAXSIM
         IF (NODDY) WRITE(LUPRI,'(A)')
     *         ' Test calculation of orbital part of Right Hand Side.'
C
C        TWOINT
C
         IF (RUNINT) THEN
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in TWOINT:',IPRINT
            END IF
            IF (IPRNTA + IPRNTB + IPRNTC + IPRNTD .GT. 0) THEN
               WRITE (LUPRI,'(2A,4I3)') ' Extra output for the',
     *            ' following shells:', IPRNTA, IPRNTB, IPRNTC, IPRNTD
               IF (RETUR) WRITE (LUPRI,'(A)')
     *            'Program will exit TWOINT after these shells.'
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' TWOINT skipped in this run.'
         END IF
C
C        GETFD
C
         IF (RUNFCK) THEN
            IF (IPRFD .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in GETFD :',IPRFD
            END IF
            IF (FCKTST) WRITE (LUPRI, '(A)')
     *         ' Testing of Fock matrices.'
            IF (AOMAT) WRITE (LUPRI, '(A)')
     *         ' AO printing of SO Fock matrices.'
            IF (TWOH1) THEN
               WRITE (LUPRI, '(A)')
     *           ' One-electron Hamiltonian multiplied by two in GETFD.'
               WRITE (LUPRI, '(2A)') ' WARNING: TWOH1 not implemented',
     *            ' in this version.'
               CALL QUIT('TWOH1 not implemented.')
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' GETFD skipped in this run.'
         END IF
C
C        GETSD
C
         IF (RUNSDR) THEN
            IF (IPRSDR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in GETSD :',IPRSDR
            END IF
            IF (DIAGTD) WRITE (LUPRI, '(A)')
     *          ' Diagonal overlap matrices (TD(I,I) = 0.25)'
            IF (SDRTST) WRITE (LUPRI, '(A)')
     *         ' AO printing of SO overlap matrices.'
         ELSE
            WRITE (LUPRI,'(A)') ' GETSD skipped in this run.'
         END IF
C
C        DERTRA
C
         IF (RUNTRA) THEN
            IF (IPRTRA .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in DERTRA:',IPRTRA
            END IF
            IF (TRATST) WRITE (LUPRI,'(A)')' Testing on transformation.'
         ELSE
            WRITE (LUPRI,'(A)') ' DERTRA skipped in this run.'
         END IF
C
C        GETGDY
C
         IF (RUNGDY) THEN
            IF (IPRGDY .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in GETGDY:',IPRGDY
            END IF
            IF (TESTFS) WRITE (LUPRI,'(A)')
     *         ' Testing on one-index transformation of Fock matrices'
            IF (NOFD) WRITE (LUPRI,'(A)')
     *         ' FD not included in calculation of Y matrices'
            IF (NOFS) WRITE (LUPRI,'(A)')
     *         ' FS not included in calculation of Y matrices'
            IF (NOSSF) WRITE (LUPRI,'(A)')
     *         ' SSF not included in calculation of Y matrices'
         ELSE
            WRITE (LUPRI,'(A)') ' GETGDY skipped in this run.'
         END IF
         IF (CUT) THEN
            WRITE (LUPRI,'(/,A)') ' Program is stopped after RHSIDE.'
         END IF
C
C        PTRAN
C
         IF (RUNPTR) THEN
            IF (IPRPTR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in PTRAN:',IPRPTR
            END IF
            IF (IPRSOR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in PSORT:',IPRSOR
            END IF
            IF (NODPTR) WRITE(LUPRI,'(A)')
     &         ' Test calculation transformed two-electron densities.'
         ELSE
            WRITE (LUPRI,'(A)') ' PTRAN skipped in this run.'
         END IF
      END IF
      RETURN
      END
C  /* Deck rhsini */
      SUBROUTINE RHSINI
C
C     Initialize /CBIRHS/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbirhs.h>
C
      SKIP   =.NOT.(MOLHES .OR. DIPDER .OR. SHIELD .OR. VCD .OR. MAGSUS
     &                     .OR. SPNSPN .OR. ECD .OR. VROA .OR. MOLGFA
     &                     .OR. SPINRO .OR. MCD .OR. VERDET .OR. OPTROT
     &                     .OR. QPGRAD)
      CUT    = .FALSE.
      RUNINT = .TRUE.
      RUNFCK = .TRUE.
      IPRALL = IPRDEF
      IPRFD  = IPRDEF
      IPRINT = IPRDEF
      IPRPTR = IPRDEF
      IPRSOR = IPRDEF
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      MAXSIM = 3
      DMAGGD(1) = SHIELD .OR. VCD .OR. MAGSUS .OR. ECD .OR. VROA
     &                   .OR. MOLGFA .OR. MCD .OR. VERDET .OR. SPINRO
     &                   .OR. OPTROT
      DMAGGD(2) = SHIELD .OR. VCD .OR. MAGSUS .OR. ECD .OR. VROA
     &                   .OR. MOLGFA .OR. MCD .OR. VERDET .OR. SPINRO
     &                   .OR. OPTROT
      DMAGGD(3) = SHIELD .OR. VCD .OR. MAGSUS .OR. ECD .OR. VROA
     &                   .OR. MOLGFA .OR. MCD .OR. VERDET .OR. SPINRO
     &                   .OR. OPTROT
      TWOH1  = .FALSE.
      NOH1   = .FALSE.
      NOH2   = .FALSE.
      FCKTST = .FALSE.
      RUNSDR = .TRUE.
      RUNTRA = .TRUE.
      RUNPTR = .TRUE.
      IPRSDR = IPRDEF
      IPRTRA = IPRDEF
      IPRIN4 = 0
      IPRIN6 = 0
      SDRTST = .FALSE.
      UNDIF  = .FALSE.
      TRATST = .FALSE.
      RUNGDY = MOLHES .OR. DIPDER .OR. SHIELD .OR. VCD .OR. MAGSUS
     &                .OR. ECD .OR. VROA .OR. MOLGFA .OR. MCD 
     &                .OR. VERDET .OR. SPINRO .OR. OPTROT .OR. QPGRAD
      IPRGDY = IPRDEF
      NODC   = .FALSE.
      NODV   = .FALSE.
      NOPV   = .FALSE.
      TESTFS = .FALSE.
      DIAGTD = .FALSE.
      NOFD   = .FALSE.
      NOFS   = .FALSE.
      NOSSF  = .FALSE.
      NOORTH = .FALSE.
      GDHAM  = .FALSE.
      RETUR  = .FALSE.
      TKTIME = .FALSE.
      AOMAT  = .FALSE.
      NODDY  = .FALSE.
      NODPTR = .FALSE.
      RALLCO = .FALSE.
      RETURN
      END
C  /* Deck rhside */
      SUBROUTINE RHSIDE(DOSFD,DOREAL,DOIMAG,DOBRHS,DOPRHS,HESFS1,
     &                  WORK,LWORK)
C
C     131088 tuh
C
C     This subroutine calculates the differentiated gradients,
C     Y vectors and overlap matrices,
C     for the right hand side of the MCSCF response equations.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <maxorb.h>
#include <aovec.h>
      PARAMETER (D0 = 0.D0, D2 = 2.D0, DP5=0.50D0)
      PARAMETER (MXDMAT = 50)
C
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION EMYD(3),  WORK(LWORK), ISYMDM(MXDMAT), IFCTYP(MXDMAT),
     &          HESFS1(*)
      LOGICAL   DERIV(3), NOORBI, LONDON, DOSFD, DOREAL, DOIMAG,
     &          DOBRHS, DOPRHS, ANTI, DDFOCK, PERTUR, DOTD, FEXIST,
     &          NOBLK, OLDDX, FOUND, PANTI, DIA2SO, ZFS2EL
C
C Used from common blocks:
C   CBIRHS : *
C   GNRINF : PARCAL
C   INFPAR : NODTOT
C
#include <cbirhs.h>
#include <gnrinf.h>
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <dorps.h>
#include <gdvec.h>
#include <energy.h>
#include <suscpt.h>
#include <inforb.h>
#include <inftap.h>
#include <inflin.h>
#include <infpar.h>
#include <cbieri.h>
#include <dftcom.h>
#include "infinp.h"
C
C     For GDHAM true:
C
      DIMENSION JCRD(MXCOOR)
      CHARACTER*8 GDHAM8(4), KEY*4
      DATA GDHAM8 /'********', 'ABACUS  ', '********', 'GDHAMILT'/
C
      IF (SKIP) RETURN
      CALL QENTER('RHSIDE')
C
      KFRSAV = 1
      I2TYP = 0
      DOTD = .TRUE.
C
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRMAX .GT. 0) WRITE (LUPRI,'(A/)')
     *   '  <<<<<<<<<< Output from RHSIDE >>>>>>>>>>'
C
C     Set for RHF
C
      IF (NASHT .EQ. 0) THEN
         NODV   = .TRUE.
         NOPV   = .TRUE.
         RUNTRA = .FALSE.
      END IF
C
C     If GDHAM write LUINTD with derived Hamiltonian and Fock matrices
C     for post-ABACUS programs (e.g. response calculation of non-adia-
C     batic coupling elements). Leading records are written here,
C     matrices for each coordinate are written by WRINTD called from
C     GETGDY. (860401-hjaaj)
C
      IF (GDHAM) THEN
         LUGDH = -1
         CALL GPOPEN(LUGDH,'ABACUS.GDHAM','NEW',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUGDH
         WRITE (LUGDH) GDHAM8
         ICRD = 0
         NCRD = 0
         DO 24 IATOM = 1,NUCIND
            DO 22 ICOOR = 1,3
               ICRD = ICRD + 1
               IF (DCORGD(IATOM,ICOOR,1)) THEN
                  NCRD = NCRD + 1
                  JCRD(NCRD) = ICRD
               END IF
  22        CONTINUE
  24     CONTINUE
         MCRD = MAX(NCRD,4)
         WRITE (LUGDH) NISHT,NASHT,NORBT,NCRD,(JCRD(I),I=1,MCRD)
      END IF
C
C     **********************************************
C     ***** Loop over magnetic field and atoms *****
C     **********************************************
C
      IF ((MCD .OR. VERDET) .AND. .NOT. NOLOND) 
     &     CALL GPOPEN(LUFDC,'FDCMAT','NEW','DIRECT',' ',IRAT*N2ORBX,
     &                 OLDDX)
      KWMAT = 1
      DO 100 IATOM = 0, NUCIND
         IF (IATOM .EQ. 0) THEN
            DERIV(1) = DOBRHS
            DERIV(2) = DOBRHS
            DERIV(3) = DOBRHS
            LONDON = .TRUE.
            KLAST = KFRSAV
         ELSE
            KLAST = KFRSAV
            IF (IATOM .EQ. 1 .AND. DOREAL .AND. (.NOT. NODIFC)) THEN
               KINTRP = KLAST
               KINTAD = KINTRP + (9*MXCENT + 1)/IRAT
               KLAST  = KINTAD + (9*MXCENT + 1)/IRAT
               LWRK   = LWORK - KLAST + 1
               IF (KLAST.GT.LWORK) CALL STOPIT('RHSIDE','GET1IN',
     &                                         KLAST,LWORK)
               NPATOM = 0
               NCOMP  = 0
               CALL GET1IN(DUMMY,'SQHDOR ',NCOMP,WORK(KLAST),LWRK,
     &              LABINT,WORK(KINTRP),WORK(KINTAD),IDUMMY,
     &              .TRUE.,NPATOM,.FALSE.,DUMMY,.FALSE.,DUMMY,IPRGDY)
               IF (NFIELD .GT. 0) THEN
                  NCOMP = 0
                  CALL GET1IN(DUMMY,'DPLGRA ',NCOMP,WORK(KLAST),LWRK,
     &                    LABINT,WORK(KINTRP),WORK(KINTAD),IDUMMY,.TRUE.
     &                    ,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,IPRGDY)
               END IF
            END IF
            DERIV(1) = DOREAL .AND. DCORGD(IATOM,1,1)
            DERIV(2) = DOREAL .AND. DCORGD(IATOM,2,1)
            DERIV(3) = DOREAL .AND. DCORGD(IATOM,3,1)
            LONDON = .FALSE.
            CALL GPOPEN(LU2DER,'AO2DRINT','NEW',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
         END IF
C
C     Restart checks
CKR   These are no longer properly used. For the time being IDORCI used
CKR   used for new spin-spin coupling feature.
C
ckr         ICOUNT = 0
ckr         ITOT   = 0
ckr         DO 101 ISYM = 1, 8
ckr            DO 101 ICOOR = 1, 3
ckr               IF (IATOM .GT. 0) THEN
ckr                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM - 1,1)
ckr                  IF (ISCOOR .GT. 0) THEN
ckr                     ITOT = ITOT + 1
ckr                     IF (IDORCI(ISCOOR,1) .GT. 0) ICOUNT = ICOUNT + 1
ckr                  END IF
ckr               ELSE
ckr                  ISCOOR = 3*NUCDEP + IPTAX(ICOOR,2)
ckr                  ITOT = ITOT + 1
ckr                  IF (IDORCI(ISCOOR,2) .NE.0) ICOUNT = ICOUNT + 1
ckr               END IF
ckr 101     CONTINUE 
ckr         IF (ICOUNT .EQ. ITOT) GO TO 100
C
C     If all directions for this atom calculated in last run, skip
C     the construction of the right-hand sides for this atom.
C     End of restart section
C
         IF (DERIV(1) .OR. DERIV(2) .OR. DERIV(3)) THEN
            IF (IPRMAX .GT. 0) THEN
               WRITE(LUPRI,'(//)')
               IF (IATOM .EQ. 0) THEN
                  CALL HEADER('Differentiation with respect to field',1)
               ELSE
                  WRITE(LUPRI,'(A,I3)')' Differentiation for nucleus',
     &                                   IATOM
                  WRITE(LUPRI,'(A/)')  ' ------------------------------'
               END IF
            END IF
C
C           Allocation for MAKEGD
C
C           Read in from SIRIFC:
C
            KCMO   = KLAST
            KDV    = KCMO   + NCMOT
C
C           Calculated in GETFD:
C           K.Ruud comment: FDC and FDV should be of dimension N2ORBX.
C           However, in order to save memory, they are used as temporary
C           arrays for AO-matrices in GETFD
C
            KDFTFD = KDV    + NNASHX
            IF (DFTADD) THEN
               KFDC   = KDFTFD + 3*N2BASX
            ELSE
               KFDC   = KDFTFD
            END IF
            KFDV   = KFDC   + 3*N2BASX
C
C           Calculated in GETSD:
C
            IF (NODV) THEN
               KTD    = KFDV
            ELSE
               KTD    = KFDV   + 3*N2ORBX
            END IF
C
C           Read in from SIRIFC:
C           KTD should have dimension N2ORBX, but in order to be usable
C           as temporary array, it has been redimensioned to N2BASX
C
            KFOCK = KTD + 3*N2BASX
            KPV   = KFOCK  + N2ORBT
C
            KTEMP = KPV + N2ASHX*N2ASHX
C            IF (SHIELD .OR. SPINRO) THEN
            IF (DOPRHS) THEN
               KWMAT = KTEMP
               KTEMP = KWMAT + 3*N2ORBX
            ELSE
               KWMAT = KTEMP
               KTEMP = KWMAT
            END IF
            IF (MAGSUS .OR. MOLGFA) THEN
               KSYMAT = KTEMP
               KWRKGD  = KSYMAT + 3*N2ORBX
            ELSE
               KSYMAT = KTEMP
               KWRKGD = KSYMAT + 3*N2ORBX
            END IF
            LWRKGD = LWORK  - KWRKGD + 1
            IF (KWRKGD.GT.LWORK)
     &         CALL STOPIT('RHSIDE','MAKEGD',KWRKGD,LWORK)
C
C           Allocation for GETFD
C
            KGETFD = KTD
            LGETFD = LWORK - KGETFD + 1
C
C           Allocation for GETSD
C
            KGETSD = KFOCK
            LGETSD = LWORK -   KGETSD + 1
            IF (IATOM .GT. 0) THEN
               NOORBI = NOORBT(IATOM)
               IF (HELFEY) NOORBI = .TRUE.
            ELSE 
               NOORBI = .FALSE.
            END IF
C
C           *****************************************************
C           ***** Differentiated two-electron AO integrals  *****
C           ***** and two-electron part of susceptibilities *****
C           *****************************************************
C
            IF (RUNINT .AND. .NOT.NOH2 .AND. .NOT.NOORBI
     &         .AND..NOT.(IATOM .EQ. 0 .AND. NOLOND)
     &         .AND..NOT.(FCKDDR.AND.EXPFCK.AND.IATOM.GT.0)) THEN
C
               IF (IATOM .EQ. 0) THEN
                  ITYPE = 5
                  IF (FCKDDR) ITYPE = -5
               ELSE
                  ITYPE = 1
                  IF (FCKDDR.AND..NOT.EXPFCK) ITYPE = 8
               END IF
CTROND
               DDFOCK = (ITYPE.EQ.-5).OR.(ITYPE.EQ.8)
CTROND
               IF ((MAGSUS .OR. MOLGFA) .AND. NASHT .GT. 0 .AND.
     &              RUNPTR) THEN
                  NOBLK = .TRUE.
               ELSE
                  NOBLK = .FALSE.
               END IF
               KJSTRS = KLAST
               KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 + 1)/IRAT
               KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 + 1)/IRAT
               KIORBS = KNCONT + (MXSHEL*MXAOVC*2 + 1)/IRAT
               KJORBS = KIORBS + (MXSHEL*MXAOVC + 1)/IRAT
               KKORBS = KJORBS + (MXSHEL*MXAOVC + 1)/IRAT
               KLAST  = KKORBS + (MXSHEL*MXAOVC + 1)/IRAT
               IF (KLAST .GT. LWORK) CALL STOPIT('TWOEXP','PAOVEC',
     &                                           KLAST,LWORK)
               LWRK   = LWORK - KLAST + 1
               CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                     WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),
     &                     IATOM,NOBLK,IPRALL)
               KLAST = KJORBS
               LWRK = LWORK - KLAST + 1
C
               IF (ITYPE .EQ. 1) THEN
                  NDMAT  = 1
                  NFMAT  = NDMAT
                  KDMAT  = KLAST
                  KFMAT  = KLAST
C                 ... so KDMAT, KFMAT is defined in subroutine calls below,
C                     also when WORK(KFMAT) is not used.
                  MAXDIF = 1
               ELSE IF (ABS(ITYPE) .EQ. 5) THEN
                  NDMAT = 2
                  IF (NODV) NDMAT = 1
                  IF (MAGSUS .OR. MOLGFA) THEN
                     MAXDIF = 2
                  ELSE
                     MAXDIF = 1
                     IF (.NOT.FCKDDR) NDMAT = 1
                  END IF
                  NFMAT  = NDMAT
C
C                 Transform densities
C
                  IF ((MAGSUS .OR. MOLGFA) .AND. NASHT .GT. 0) THEN
                     IF (RUNPTR) THEN
C
C                       Symmetric part
C
                        ANTI =   .FALSE.
                        PANTI =  .FALSE.
                        DIA2SO = .FALSE.
                        ZFS2EL = .FALSE.
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('START ',TIMSTR,TIMEND)
                        CALL PTRAN(NODPTR,WORK(KLAST),LWRK,IPRPTR,ANTI,
     &                             PANTI,DIA2SO,ZFS2EL)
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('PTRAN ',TIMSTR,TIMEND)
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('START ',TIMSTR,TIMEND)
#if defined (SYS_CRAY) || defined (SYS_NEC) || defined (SYS_T90)
                        LWUSE = LWRK
#else 
                        CALL MEMSRT(MEMMIN,WORK(KNCONT))
C     
C     1000000 is just an estimate of what the rest of the sorting 
C     program needs
C
                        LWUSE = MEMMIN + 1000000
                        IF (LWUSE .GT. LWRK)
     &                       CALL STOPIT('TWOEXP','PSORG',LWUSE,LWRK)
#endif
                        CALL GPOPEN(LUPAO,' ','UNKNOWN',' ',' ',IDUMMY,
     &                              .FALSE.)
                        CALL PSORG(WORK,WORK(KLAST),LWUSE,WORK(KNCONT),
     &                             IPRSOR,ANTI,PANTI)
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('PSORT ',TIMSTR,TIMEND)
                        ANTI = .TRUE.
C
C                       Antisymmetric part
C
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('START ',TIMSTR,TIMEND)
                        CALL PTRAN(NODPTR,WORK(KLAST),LWRK,IPRPTR,ANTI,
     &                             PANTI,DIA2SO,ZFS2EL)
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('PTRAN ',TIMSTR,TIMEND)
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('START ',TIMSTR,TIMEND)
                        CALL GPOPEN(LUPAS,' ','UNKNOWN',' ',' ',IDUMMY,
     &                              .FALSE.)
                        CALL PSORG(WORK,WORK(KLAST),LWUSE,WORK(KNCONT),
     &                             IPRSOR,ANTI,PANTI)
                        IF (IPRPTR .GE. 3)
     &                  CALL TIMER('PSORT ',TIMSTR,TIMEND)
                     END IF
                  END IF
C
                  KFMAT = KLAST
                  IF (FCKDDR) THEN
                     NFMAT = 3*NDMAT
                     KLAST = KFMAT + NFMAT*N2BASX
                     CALL DZERO(WORK(KFMAT),NFMAT*N2BASX)
                  END IF
                  KDMAT = KLAST
                  KLAST = KDMAT + N2BASX*NDMAT
                  LWRK  = LWORK - KLAST + 1
                  IF (KLAST .GT. LWORK) THEN
                     CALL STOPIT('RHSIDE','GETDMT',KLAST,LWORK)
                  END IF
                  IF ((MAGSUS .OR. MOLGFA) .OR. FCKDDR) THEN
                     CALL GETDMT(WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,
     &                           NODC,NODV,.TRUE.,IPRINT)
                  END IF
               ELSE IF (ITYPE .EQ. 8) THEN
                  MAXDIF = 1
                  NDMAT = 2
                  IF (NODV) NDMAT = 1
                  NFMAT  = 3*NDMAT*NUCDEG(IATOM)
                  LFMAT = NFMAT*N2BASX
                  KFMAT  = KLAST
                  KDMAT  = KFMAT + LFMAT
                  KLAST  = KDMAT + N2BASX*NDMAT
                  LWRK   = LWORK - KLAST + 1
                  IF (KLAST .GT. LWORK) THEN
                     CALL STOPIT('RHSIDE','ITYPE8',KLAST,LWORK)
                  END IF
                  CALL GETDMT(WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,
     &                        NODC,NODV,.TRUE.,IPRINT)
                  CALL DZERO(WORK(KFMAT),LFMAT)
               END IF
C
               IF (NFMAT .GT. MXDMAT) THEN
                  WRITE (LUPRI,'(/A/A,2I10)')
     &               ' Error in RHSIDE: NFMAT .gt. MXDMAT',
     &               ' NFMAT, MXDMAT =',NFMAT,MXDMAT
                  CALL QUIT('MXDMAT parameter too small in RHSIDE.')
               END IF
               DO I = 1,NFMAT
                  ISYMDM(I) = 0
                  IFCTYP(I) = 13
               END DO
               K2INT  = KLAST
               L2INT  = LWORK - K2INT + 1
C
#if defined (VAR_MPI) || defined (VAR_PVM)
               IF (PARCAL .AND. DDFOCK) THEN
C
                  KNSTAT = KLAST
                  KLAST  = KNSTAT + (NODTOT + 1)/IRAT
                  IF (KLAST .GT. LWRK)
     &               CALL STOPIT('RHSIDE','PARDRV',KLAST,LWRK)
                  LWRK = LWRK - KLAST
                  IPRTYP = 1
                  HFXM0 = HFXMU
                  HFXMU = D0
                  CALL PARDRV(WORK(KFMAT),WORK(KDMAT),NDMAT,ISYMDM,
     &                        IFCTYP,WORK(KLAST),WORK(KNSTAT),DUMMY,
     &                        LWRK,ITYPE,MAXDIF,IATOM,NODV,NOPV,.FALSE.,
     &                        TKTIME,RETUR,IPRINT,IPRTYP,
     &                        IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,
     &                        DUMMY)
                  HFXMU = HFXM0
                  IF (HFXMU.NE.D0) THEN
                     KFCTMP = KLAST
                     KLAST  = KFCTMP + NFMAT*N2BASX
                     IF (KLAST .GT. LWRK)
     &                    CALL STOPIT('TWOEXP','PARDRV',KLAST,LWRK)
                     LWRK = LWRK - NFMAT*N2BASX
                     
                     HFXFC0 = HFXFAC
                     HFXFAC = HFXATT
                     CALL DCOPY(NFMAT*N2BASX,WORK(KFMAT),1,
     &                    WORK(KFCTMP),1)
                     CALL PARDRV(WORK(KFMAT),WORK(KDMAT),NDMAT,ISYMDM,
     &                    IFCTYP,WORK(KLAST),WORK(KNSTAT),DUMMY,
     &                    LWRK,ITYPE,MAXDIF,IATOM,NODV,NOPV,.FALSE.,
     &                    TKTIME,RETUR,IPRINT,IPRTYP,
     &                    IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,
     &                    DUMMY)
                     CALL DAXPY(NFMAT*N2BASX,1D0,WORK(KFCTMP),1,
     &                    WORK(KFMAT),1)
                     HFXFAC = HFXFC0
                  END IF               
                  IF (IPRALL .GT. 2 .AND. (MAGSUS .OR. MOLGFA)) THEN
                     CALL HEADER('Two-electron expectation value '
     &                    //'contribution to susceptibilities',-1)
                     CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
                  END IF
               ELSE
#endif
C                IF (.NOT.DOERIP) THEN
                    HFXM0 = HFXMU
                    HFXMU = D0
                    CALL TWOINT(WORK(K2INT),L2INT,DUMMY,WORK(KFMAT),
     &                          WORK(KDMAT),NDMAT,ISYMDM,IFCTYP,
     &                          DUMMY,IDUMMY,NUMDIS,1,ITYPE,MAXDIF,
     &                          IATOM,NODV,NOPV,.FALSE.,TKTIME,
     &                          IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                          RETUR,IDUMMY,I2TYP,WORK(KJSTRS),
     &                          WORK(KNPRIM),WORK(KNCONT),
     &                          WORK(KIORBS),IDUMMY,IDUMMY,DUMMY,
     &                          DUMMY,DUMMY,DUMMY,.FALSE.)
                    HFXMU = HFXM0
                    IF (HFXMU.NE.D0) THEN
                       HFXFC0 = HFXFAC
                       HFXFAC = HFXATT
                       CALL TWOINT(WORK(K2INT),L2INT,DUMMY,WORK(KFMAT),
     &                             WORK(KDMAT),NDMAT,ISYMDM,IFCTYP,
     &                             DUMMY,IDUMMY,NUMDIS,1,ITYPE,MAXDIF,
     &                             IATOM,NODV,NOPV,.FALSE.,TKTIME,
     &                             IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                             RETUR,IDUMMY,I2TYP,WORK(KJSTRS),
     &                             WORK(KNPRIM),WORK(KNCONT),
     &                             WORK(KIORBS),IDUMMY,IDUMMY,DUMMY,
     &                             DUMMY,DUMMY,DUMMY,.FALSE.)
                       HFXFAC = HFXFC0
                    END IF
C                ELSE
C                   KCCFBT = KLAST
C                   KINDXB = KCCFBT + MXPRIM*MXCONT
C                   KLAST  = KINDXB + 8*MXSHEL*MXCONT
C                   IF (KLAST.GT.LWRK) THEN
C                      CALL STOPIT('TWOEXP','ERIPRO',KLAST,LWRK)
C                   END IF
C                   LWRK = LWRK - KLAST
C                   CALL ERIPRO(DUMMY,DUMMY,NDMAT,ISYMDM,IFCTYP,
C    &                          .FALSE.,.TRUE.,
C    &                          IPRINT,WORK(KCCFBT),WORK(KINDXB),
C    &                          WORK(KLAST),LWRK)
C                END IF
#if defined (VAR_MPI) || defined (VAR_PVM)
               END IF
#endif
               IF(DDFOCK) THEN
                  PERTUR = ITYPE.EQ.8
                  CALL SKLFCK(WORK(KFMAT),DUMMY,WORK(KLAST),LWRK,
     &            IPRINT,.FALSE.,
     &            DDFOCK,.FALSE.,PERTUR,NODV,MAXDIF,LONDON,NDMAT,
     &            ISYMDM,IFCTYP,IATOM,.FALSE.)
               ENDIF

C
               IF ((MAGSUS .OR. MOLGFA) .AND. LONDON .AND. NASHT.GT.0
     &                    .AND. RUNPTR) CALL GPCLOSE(LUPAS,'DELETE')
C
C              Fock matrices for geometrical derivatives
C
               IF (ITYPE .EQ. 8) THEN
                  DO 200 I = 1, 3*NDMAT
                     CALL WRITDX(LUDFCK,I,
     &                           IRAT*NNBASX,WORK(KFMAT + (I-1)*N2BASX))
  200             CONTINUE
               END IF
C
C              Fock matrices for magnetic field derivatives
C
               IF (ITYPE .EQ. -5) THEN
                  IF (MOLHES .OR. DIPDER .OR. QPGRAD .OR. VCD) THEN
                    IREC = 3
                    IF (EXPFCK) IREC = 3*NUCIND
                  ELSE
                    IREC = 0
                  END IF
                  IF (.NOT.NODV) IREC = 2*IREC
                  DO 300 I = 1, 3*NDMAT
                     CALL WRITDX(LUDFCK,IREC + I,
     &                           IRAT*NNBASX,WORK(KFMAT + (I-1)*N2BASX))
  300             CONTINUE
               END IF
            END IF
C
C           *************************************************
C           ***** Noddy routine for construction of RHS *****
C           *************************************************
C
            IF (IATOM .EQ.0 .AND. NODDY) THEN
               CALL OFACHK(RUNSDR,NOH1,NOH2,NOORTH,WORK(K2INT),L2INT,
     &                     IPRALL)
            END IF
C
C           **********************
C           ***** CMO and DV *****
C           **********************
C
            CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KGETFD),LGETFD)
            IF (.NOT.FOUND)
     &         CALL QUIT('RHSIDE error: CMO not found on SIRIFC')
            IF (NASHT .GT. 0) THEN
               IF (NODV) THEN
                  CALL DZERO(WORK(KDV),  NNASHX)
               ELSE
                  CALL RD_SIRIFC('DV',FOUND,WORK(KDV),
     &               WORK(KGETFD),LGETFD)
                  IF (.NOT.FOUND)
     &            CALL QUIT('RHSIDE error: DV not found on SIRIFC')
               END IF
            END IF
C
C           ************************************************************
C           ** AO differentiated inactive and active MO Fock matrices **
C           ************************************************************
C
            IF (RUNFCK) THEN
               IF (IPRFD .GT. 2)
     &            CALL TIMER('START ',TIMSTR,TIMEND)
               CALL GETFD(IATOM,NOORBI,TWOH1,NOH1,NOH2,NODC,NODV,IPRFD,
     &                    DERIV,EMYD,WORK(KCMO),WORK(KDV),WORK(KFDC),
     &                    WORK(KFDV),WORK(KDFTFD),WORK(KGETFD),LGETFD,
     &                    AOMAT,DOREAL,IPRINT)
               IF (IPRFD .GT. 2)
     &            CALL TIMER('GETFD ',TIMSTR,TIMEND)
            END IF
C
C           *********************************************************
C           ***** AO differentiated overlap matrices (MO basis) *****
C           *********************************************************
C
            IF (RUNSDR .AND. .NOT.(IATOM .EQ. 0 .AND. NOLOND) .AND.
     &           .NOT. (IATOM .GT. 0 .AND. HELFEY)) THEN
               IF (IPRSDR .GT. 2)
     &            CALL TIMER('START ',TIMSTR,TIMEND)
               IF (IATOM .EQ. 0) THEN
                  IF (NODIFC) THEN
                     KEY = 'SMAT'
                  ELSE
                     KEY = 'BMAT'
                  END IF
               ELSE
                  IF (NODIFC) THEN
                     KEY = 'OMAT'
                  ELSE
                     KEY = 'DMAT'
                  END IF
               END IF
               CALL GETSD(DERIV,WORK(KCMO),WORK(KTD),
     &                    WORK(KGETSD),LGETSD,IATOM,DIAGTD,
     &                    KEY,DOTD,IPRSDR,IPRINT)
               IF (IPRSDR .GT. 2)
     &            CALL TIMER('GETSD ',TIMSTR,TIMEND)
            ELSE
               CALL DZERO(WORK(KTD),3*N2BASX)
            END IF
C
C           ******************************
C           ***** Read Fock matrix F *****
C           ******************************
C
            CALL RD_SIRIFC('FOCK',FOUND,WORK(KFOCK),WORK(KWRKGD),LWRKGD)
            IF (.NOT.FOUND)
     &      CALL QUIT('RHSIDE error: FOCK not found on SIRIFC')
C
C           ***** Test Fock matrix FOCK and overlap matrices TD *****
C
            IF (RUNSDR .AND. SDRTST .AND. IATOM.GT.0) THEN
               CALL TESTSD(WORK(KTD),WORK(KFOCK),IATOM,DERIV,IPRSDR)
            END IF
C
C           **************************************
C           ***** Construct GD and Y vectors *****
C           **************************************
C
            CALL MAKEGD(DOSFD,EMYD,WORK(KCMO),WORK(KDV),WORK(KFDC),
     &                  WORK(KFDV),WORK(KTD),WORK(KFOCK),
     &                  WORK(KPV),WORK(KWMAT),WORK(KSYMAT),WORK(KDFTFD),
     &                  WORK(KWRKGD),LWRKGD,IATOM,DERIV,MAXSIM,NOORBI,
     &                  NOH1,NOH2,NODC,NODV,NOPV,
     &                  NOORTH,LONDON,RUNTRA,TRATST,FCKTST,
     &                  RUNGDY,LUGDH,IPRMAX,IPRALL,IPRTRA,IPRGDY)
         END IF
C
C        Close and delete AO2MGINT if it has been used
C
         IF (LU2DER .GT. 0) CALL GPCLOSE(LU2DER,'DELETE')
         IF (DOPRHS .AND. IATOM .EQ. 0) THEN
            IF (.NOT. (SHIELD .OR. SPINRO)) THEN
               KWRK = KWMAT
            ELSE
               KWRK = KWMAT + 3*N2ORBX
            END IF
            LWRK = LWORK - KWRK + 1
            CALL PSORHS(WORK(KWMAT),LABINT,WORK(KWRK),LWRK)
         END IF
C
  100 CONTINUE
C
      IF ((MCD .OR. VERDET) .AND. .NOT. NOLOND) 
     &     CALL GPCLOSE(LUFDC,'KEEP')
C
      IF (MOLHES .AND. .NOT. NODIFC .AND. .NOT. HELFEY) THEN
         CALL REORT2(WORK(KTD),WORK(KTD + N2ORBX),
     &               WORK(KTD + 2*N2ORBX),WORK(KFOCK),HESFS1,
     &               IPRGDY,LUPRI)
      END IF
C
C     Close unit for two-electron AO integrals
C
      IF (GDHAM) CALL GPCLOSE(LUGDH,'KEEP')
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after RHSIDE as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in RHSIDE) *****')
      END IF
C
      IF (IPRALL .GE. 2) CALL TIMER('RHSIDE',TIMEIN,TIMOUT)
C
      CALL QEXIT('RHSIDE')
      RETURN
      END
C  /* Deck makegd */
       SUBROUTINE MAKEGD(DOSFD,EMYD,CMO,DV,FDC,FDV,TD,
     &                   FOCK,PV,WMAT,SYMAT,DFTFD,WORK,LWORK,IATOM,
     &                   DERIV,MAXSIM,NOORBI,NOH1,NOH2,NODC,NODV,NOPV,
     &                   NOORTH,LONDON,RUNTRA,TRATST,FCKTST,RUNGDY,
     &                   LUGDH,IPRMAX,IPRALL,IPRTRA,IPRGDY)
C
C     131088 tuh
C
C
C     Purpose:
C
C        Calculate total derivatives of gradient (AO-derivative terms
C        + one-index transformed terms).
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <iratdef.h>
#include <dummy.h>
      PARAMETER (D1 = 1.0D0, DP5 = 0.5D0, DM1 = -1.D0, D0 = 0.D0,
     &           D2 = 2.0D0)
C
      DIMENSION EMYD(3), CMO(NCMOT), WMAT(N2ORBX,3),
     &          DV(NNASHX), PV(*), SYMAT(N2ORBX,*),
     &          FDC(N2BASX,3), FDV(N2ORBX,3), TD(N2BASX,3),
     &          FOCK(N2ORBT), DFTFD(N2BASX,3),WORK(LWORK)
      LOGICAL DERIV(3), NODC, NODV, NOPV, RUNTRA, TRATST,FCKTST, RUNGDY,
     &        NOORBI, NOORTH, NOH1, NOH2 ,LONDON, DOSFD, SUSCEP, ANTI,
     &        FOUND
C
C
#include <abainf.h>
#include <cbisol.h>
#include <cbdtra.h>
#include <symmet.h>
#include <nuclei.h>
#include <dorps.h>
#include <gdvec.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infdim.h>
#include <nodint.h>
#include <inflin.h>
#include <infinp.h>
#include <energy.h>
#include <suscpt.h>
#include <abares.h>
#include <dftcom.h>
C
      CALL QENTER('MAKEGD')
      SUSCEP = IATOM .EQ. 0 .AND. (MAGSUS .OR. MOLGFA)
C
      IF (IPRMAX .GT. 0) WRITE (LUPRI,'(/A/)')
     &   '  <<<<<<<<<< Output from MAKEGD >>>>>>>>>>'
      IF (IPRMAX .GT. 5) THEN
         WRITE (LUPRI,'(1X,A,I5)')  ' IATOM ', IATOM
         WRITE (LUPRI,'(1X,A,3L5)') ' DERIV ', (DERIV(I),I=1,3)
         WRITE (LUPRI,'(1X,A,I5)')  ' MAXSIM ', MAXSIM
         WRITE (LUPRI,'(1X,A,3L5)') ' NODC, NODV, NOPV ',NODC,NODV,NOPV
         WRITE (LUPRI,'(1X,A,L5)')  ' RUNTRA ', RUNTRA
         WRITE (LUPRI,'(1X,A,3F12.6)') ' EMYD ', (EMYD(I),I=1,3)
         IF (IPRMAX .GT. 10) THEN
            CALL HEADER('DV matrix',-1)
            CALL OUTPAK(DV,NASHT,1,LUPRI)
            IF (DERIV(1)) THEN
               CALL HEADER('FDC matrix - x direction',-1)
               CALL OUTPUT(FDC(1,1),1,NORBT,1,NORBT,NBAST,NBAST,1,
     &                     LUPRI)
               IF (.NOT. NODV) THEN
                  CALL HEADER('FDV matrix - x direction',-1)
                  CALL OUTPUT(FDV(1,1),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                        LUPRI)
               END IF
               CALL HEADER('TD matrix - x direction',-1)
               CALL OUTPUT(TD(1,1),1,NORBT,1,NORBT,NBAST,NBAST,1,
     &                     LUPRI)
            END IF
            IF (DERIV(2)) THEN
               CALL HEADER('FDC matrix - y direction',-1)
               CALL OUTPUT(FDC(1,2),1,NORBT,1,NORBT,NBAST,NBAST,1,
     &                     LUPRI)
               IF (.NOT. NODV) THEN
                  CALL HEADER('FDV matrix - y direction',-1)
                  CALL OUTPUT(FDV(1,2),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                        LUPRI)
               END IF
               CALL HEADER('TD matrix - y direction',-1)
               CALL OUTPUT(TD(1,2),1,NORBT,1,NORBT,NBAST,NBAST,1,
     &                     LUPRI)
            END IF
            IF (DERIV(3)) THEN
               CALL HEADER('FDC matrix - z direction',-1)
               CALL OUTPUT(FDC(1,3),1,NORBT,1,NORBT,NBAST,NBAST,1,
     &                     LUPRI)
               IF (.NOT. NODV) THEN
                  CALL HEADER('FDV matrix - z direction',-1)
                  CALL OUTPUT(FDV(1,3),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                        LUPRI)
               END IF
               CALL HEADER('TD matrix - z direction',-1)
               CALL OUTPUT(TD(1,3),1,NORBT,1,NORBT,NBAST,NBAST,1,
     &                     LUPRI)
            END IF
            CALL HEADER('Total Fock matrix',-1)
            DO 50 ISYM = 1, NSYM
               WRITE (LUPRI,'(A,I5)') ' Irrep', ISYM
               NORBI = NORB(ISYM)
               ISTRI = I2ORB(ISYM) + 1
               CALL OUTPUT(FOCK(ISTRI),1,NORBI,1,NORBI,NORBI,NORBI,
     &                     1,LUPRI)
   50       CONTINUE
         END IF
      END IF
C
      IF ((SHIELD .OR. SPINRO) .AND. IATOM .EQ. 0) 
     &   CALL DZERO(WMAT,3*N2ORBX)
C
C     **************************
C     ***** CI information *****
C     **************************
C
C     920217-hjaaj: ICSYM=IHCSYM=LSYMRF means that .TOTSYM
C        always will work for CSFs.  (Both ICSYM and IHCSYM have
C        no effect for determinants: they are overwritten by
C        CALL ABAVAR(ISYM) below.  However, only one CSF <-> det
C        transformation is defined in the current version.)
C
      KCINDX = 1
      IF (NASHT .GT. 1) THEN
         KFREE = KCINDX + LCINDX
         LFREE = LWORK  - KFREE  + 1
         IF (LFREE.GT.LWORK) CALL STOPIT('MAKEGD','GETCIX',KFREE,LWORK)
         CALL GETCIX(WORK(KCINDX),LSYMRF,LSYMRF,WORK(KFREE),LFREE,0)
CSPAS : 3/12-98 : Changes for SOPPA
      ELSE IF (ABASOP) THEN
         KFREE = KCINDX + LCINDX
CKeinSPASmehr
      ELSE
         KFREE = KCINDX
      END IF
      LFREE = LWORK  - KFREE  + 1
C
C     **************************************
C     ***** Get two-electron densities *****
C     **************************************
C
      IF (NASHT .GT. 0) THEN
         IF (LONDON) THEN
            NDIMPV = N2ASHX
         ELSE
            NDIMPV = NNASHX
         END IF
         IF (NOPV) THEN
            CALL DZERO(PV,NDIMPV*NDIMPV)
         ELSE
C
C           Case 1: Densities constructed
C
            IF (LONDON) THEN
               KCREF = KFREE
               KUDV  = KCREF + NCONRF
               KLAST = KUDV  + N2ASHX
               LWRK  = LWORK - KLAST + 1
               IF (KLAST.GT.LWORK) THEN
                  CALL STOPIT('MAKEGD','RSPDM',KLAST,LWORK)
               END IF
               CALL GETREF(WORK(KCREF),NCONRF)
C
               ISPIN1 = 0
               ISPIN2 = 0
               KSTART = 1
               CALL RSPDM(LSYMRF,LSYMRF,NCONRF,NCONRF,WORK(KCREF),
     &                    WORK(KCREF),WORK(KUDV),PV,
     &                    ISPIN1,ISPIN2,.FALSE.,.FALSE.,
     &                    WORK(KCINDX),WORK(KLAST),KSTART,LWRK)
C              CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *                    ISPIN1,ISPIN2,TDM,NORHO2,XNDXCI,WORK,
C    *                    KFREE,LFREE)
C
C           Case 2: Densities read in
C
            ELSE
               CALL RD_SIRIFC('PV',FOUND,PV,WORK(KFREE),LFREE)
               IF (.NOT.FOUND)
     &            CALL QUIT('RHSIDE error: PV not found on SIRIFC')
            END IF
         END IF
         IF (IPRMAX .GT. 10) THEN
            CALL HEADER('PV matrix',-1)
            CALL OUTPUT(PV,1,NDIMPV,1,NDIMPV,NDIMPV,NDIMPV,1,LUPRI)
         END IF
      END IF
C
C     ********************************
C     ***** Loop over symmetries *****
C     ********************************
C
      IPRLIN = IPRGDY
      IROT = 0
      ISCOOR = 0
      DO 100 ISYM = 1, NSYM
      IF (.NOT. DOSYM(ISYM)) GO TO 100
C
C        *********************************
C        ***** Initialize for ABARHS *****
C        *********************************
C
         CALL ABAVAR(ISYM,.FALSE.,IPRMAX,WORK(KFREE),LFREE)
C
C        **********************************
C        ***** Work space allocations *****
C        **********************************
C
C        SQFDTD
C
         KFDCSQ = KFREE
         KFDVSQ = KFDCSQ + N2ORBX
         KTDSQ  = KFDVSQ + N2ORBX
         KDFTSQ = KTDSQ  + N2ORBX 
         IF (DFTADD) THEN
            KUMAT = KDFTSQ + N2ORBX
         ELSE
            KUMAT = KDFTSQ
         END IF
         IF (SOLVNT) THEN
            KWSQ = KUMAT  + N2ORBX
         ELSE
            KWSQ = KUMAT
         END IF
C
C        DERTRA and LONTRA
C
         LH2DAC = N2ASHX*N2ASHX
         KQD   = KWSQ
         KH2D  = KQD + NASHT*NORBT
         KWTRA = KH2D  + LH2DAC
         LWTRA = LWORK - KWTRA + 1
C
C        ABARHS
C
         KQDI   = KWTRA
         KFC    = KQDI   + NASHT*NORBT
         KFV    = KFC    + NNORBT
         KCREF  = KFV    + NNORBT
         KGD    = KCREF  + NCONRF
         KWKRHS = KGD    + NVARPT
CSPAS 3/12-98 : changes for SOPPA
         IF (ABASOP) KWKRHS = KGD    + 2*NVARPT
CKeinSPASmehr
         LWKRHS = LWORK  - KWKRHS + 1
         LAST1  = KWKRHS
C
C        CHKRHS
C
         KH1AO  = KQDI   + NASHT*NORBT
         KH1MO  = KH1AO  + NNBAST
         KW2    = KH1MO  + NNORBT
         LW2    = LWORK  - KW2
         LAST2  = KW2
C
C        GETY
C
         LAST3 = KQDI  + NASHT*NORBT
C
         LAST = MAX(LAST1,LAST2,LAST3)
         IF (LAST.GT.LWORK) CALL STOPIT('MAKEGD',' ',LAST,LWORK)
C
C        ***** Loop over Cartesian directions *****
C
         DO 300 ICOOR = 1, 3
         IF (DERIV(ICOOR)) THEN
            IF (IATOM .NE. 0) THEN
               ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM - 1,1)
               IF (.NOT.DOPERT(ISCOOR,1)) GO TO 300
ckr               IF (IDORCI(ISCOOR,1) .GT. 0) GO TO 300
            ELSE
               IF (ISYMAX(ICOOR,2) .NE. ISYM - 1) GO TO 300
C
C              what happens when no shieldings are calculated?
C              compare open statement for the file in mcorl.
C
               ISCOOR = 3*NUCDEP + IPTAX(ICOOR,2)
ckr               IF (IDORCI(ISCOOR,2) .GT. 0) GO TO 300
            END IF
            IF (IPRMAX .GT. 0) THEN
               WRITE (LUPRI,'(//A,I3,I3,/A/)')
     &          ' Differentiation Direction and symmetry',ICOOR,ISYM,
     &          ' ============================================ '
            END IF
C
C           **********************
C           ***** H2D and QD *****
C           **********************
C
            IF (NASHT .GT. 0) THEN
               CALL DZERO(WORK(KH2D),LH2DAC)
               CALL DZERO(WORK(KQD),NASHT*NORBT)
               IF (NDSINT(ICOOR,ISYM - 1) .GT. 0) THEN
                  IF (RUNTRA .AND. .NOT.(NOH2 .OR. NOORBI)) THEN
                     IF (IATOM .NE. 0) THEN
                        DO 310 I = 1, 3
                           SRTCOO(I) = .FALSE.
 310                    CONTINUE
                        DO 320 I = 1, 8
                           SRTSYM(I) = .FALSE.
 320                    CONTINUE
                        SRTCOO(ICOOR) = .TRUE.
                        SRTSYM(ISYM) = .TRUE.
                        IF (IPRTRA .GT. 2)
     &                       CALL TIMER('START ',TIMSTR,TIMEND)
                        IPRTRD = IPRTRA
                        CALL DERTRA(IATOM,CMO,PV,WORK(KH2D),WORK(KQD),
     &                              WORK(KWTRA),LWTRA,MWORK)
                        IF (IPRTRA .GT. 2)
     &                       CALL TIMER('DERTRA',TIMSTR,TIMEND)
                        IF (TRATST) THEN
                           CALL TESTTR(IATOM,ICOOR,ISYM,PV,WORK(KH2D),
     &                                 WORK(KQD),IPRTRA)
                        END IF
                     ELSE IF (.NOT.NOLOND) THEN
                        IF (IPRTRA .GT. 2)
     &                       CALL TIMER('START ',TIMSTR,TIMEND)
                        CALL LONTRA(ICOOR,ISYM,WORK(KH2D),WORK(KQD),
     &                              CMO,PV,WORK(KWTRA),LWTRA,IPRTRA)
                        IF (IPRTRA .GT. 2)
     &                       CALL TIMER('LONTRA',TIMSTR,TIMEND)
                     END IF
                  END IF
               END IF
               CALL DCOPY(NASHT*NORBT,WORK(KQD),1,WORK(KQDI),1)
            END IF
C
C           ****************************
C           ***** FDC, FDV, and TD *****
C           ****************************
C
            IF (NOORTH) THEN
               CALL DZERO(WORK(KTDSQ),N2ORBX)
            ELSE
               IF (IPRMAX .GT. 10)
     &            WRITE(LUPRI,'(/A)') 'Calling FXFDTD with TD'
               CALL EXFDTD(TD(1,ICOOR),WORK(KTDSQ),ISYM,IPRMAX)
            END IF

            IF (IPRMAX .GT. 10)
     &         WRITE(LUPRI,'(/A)') 'Calling FXFDTD with FDC'
            CALL EXFDTD(FDC(1,ICOOR),WORK(KFDCSQ),ISYM,IPRMAX)

            IF (DFTADD) THEN
               IF (IPRMAX .GT. 10)
     &            WRITE(LUPRI,'(/A)') 'Calling FXFDTD with DFTFD'
               CALL EXFDTD(DFTFD(1,ICOOR),WORK(KDFTSQ),ISYM,IPRMAX)
            END IF
            IF (NASHT .GT. 0) THEN
               IF (IPRMAX .GT. 10)
     &            WRITE(LUPRI,'(/A)') 'Calling FXFDTD with FDV'
               CALL EXFDTD(FDV(1,ICOOR),WORK(KFDVSQ),ISYM,IPRMAX)
            ELSE
               CALL DZERO(WORK(KFDVSQ),N2ORBX)
            END IF
C
C           *************************************
C           ***** Test FDC, FDV, QD and H2D *****
C           *************************************
C
            IF (FCKTST .AND. IATOM.GT.0) THEN
               CALL TESTFD(IATOM,ICOOR,IPRMAX,NOH1,NOH2,NODC,NODV,NOPV,
     &                     WORK(KWTRA),WORK(KFDCSQ),WORK(KFDVSQ),
     &                     WORK(KQD),DV)
            END IF
C
C           **********************
C           ***** GD vectors *****
C           **********************
C
            IF (RUNGDY) THEN
               REWIND LUSIFC
               CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
               READ (LUSIFC)
               READ (LUSIFC)
               READ (LUSIFC)
               IF (NCONRF .GT. 0) THEN
                  CALL READI(LUSIFC,IRAT*NCONRF,WORK(KCREF))
               ELSE
                  READ (LUSIFC)
               END IF
               READ (LUSIFC)
               READ (LUSIFC)
               READ (LUSIFC)
               CALL READI(LUSIFC,IRAT*NNORBT,WORK(KFC))
               IF (NASHT .GT. 0) THEN
                  CALL READI(LUSIFC,IRAT*NNORBT,WORK(KFV))
               ELSE
                  READ (LUSIFC)
                  CALL DZERO(WORK(KFV),NNORBT)
               END IF
               IF (IPRGDY .GT. 10) THEN
                  CALL HEADER('Inactive undiff. Fock matrix',-1)
                  CALL OUTPKB(WORK(KFC),NORB,NSYM,1,LUPRI)
                  IF (NASHT .GT. 0) THEN
                     CALL HEADER('Active undiff. Fock matrix',-1)
                     CALL OUTPKB(WORK(KFV),NORB,NSYM,1,LUPRI)
                  END IF
               END IF
C
C              For solvent, construct differentiated Tlm and Upq
C              =================================================
C
               IF (SOLVNT .AND. IATOM .GT. 0) THEN
                  CALL SOLRHS(CMO,WORK(KTDSQ),DV,WORK(KUMAT),
     &                        WORK(KWKRHS),LWKRHS,ICOOR,ISCOOR,
     &                        ISYM,IPRMAX)
               END IF
C
C
C              *************************************
C              ***** Calculate right-hand side *****
C              *************************************
C
               IF (NVARPT .GT. 0) THEN
                  IF (SOLVNT .AND. IATOM .GT. 0) THEN
                     CALL DAXPY(N2ORBX,D1,WORK(KUMAT),1,WORK(KFDCSQ),1)
                  END IF
                  IF (IPRGDY .GT. 2)
     &               CALL TIMER('START ',TIMSTR,TIMEND)
                  NGD = 1
CSPAS : 3/12-98 : changes for SOPPA
                  IF (ABASOP) THEN
                     CALL ANGRHS(WORK(KFDCSQ),WORK(KGD),WORK(KCINDX),
     &                           WORK(KWKRHS),1,LWKRHS)
                     GRDICT = D0
                     GRDACT = D0
                  ELSE
                     CALL ABARHS(NGD,DV,PV,WORK(KFC),WORK(KFV),
     &                           WORK(KCREF),CMO,WORK(KTDSQ),
     &                           WORK(KFDCSQ),WORK(KFDVSQ),WORK(KQD),
     &                           WORK(KH2D),LH2DAC,GRDICT,GRDACT,
     &                           WORK(KGD),WORK(KCINDX),LONDON,
     &                           WORK(KDFTSQ),WORK(KWKRHS),1,LWKRHS)
                  ENDIF
CKeinSPASmehr
                  IF (IPRGDY .GT. 2)
     &               CALL TIMER('ABARHS',TIMSTR,TIMEND)
                  IF (SOLVNT .AND. IATOM .GT. 0) THEN
                     CALL DAXPY(N2ORBX,-D1,WORK(KUMAT),1,WORK(KFDCSQ),1)
                  END IF
C
C                 ***** Print right-hand side *****
C
                  IF (IPRMAX .GE. 05) THEN
                     WRITE (LUPRI,'(//A,I3,A,I2,A,I2,A)')
     &                 ' Differentiated Gradient for Atom',IATOM,
     &                 ', Coordinate',ICOOR,', and symmetry',ISYM,':'
                     IF (IPRMAX .GT. 20) THEN
                        WRITE(LUPRI,'(/,A,/)')' CSF part of gradient'
                        WRITE(LUPRI,'(6F12.6)')
     &                            (WORK(KGD+J-1), J = 1, NCONST)
                     END IF
                     WRITE(LUPRI,'(//,A)')' Orbital part of gradient'
                     IF (IPRMAX .GE. 20) THEN
                        PRFAC = 0.0D0
                     ELSE
                        PRFAC = 0.1D0
                     END IF
                     CALL PRKAP(NWOPPT,WORK(KGD+NCONST),PRFAC,LUPRI)
                  END IF
C
C                 ***** Average orbital rotations if SUPSYM and sym 1 *****
C
                  IF (ISYM.EQ.1) CALL AVERAG(WORK(KGD+NCONST),NWOPPT,1)
C
C                 Magnetic field
C                 ==============
C
                  IF (IATOM .EQ. 0) THEN
                     IF (ICOOR .EQ. 1) THEN
                        IMGLAB(ISCOOR) = 'MAGFLD.1'
                     ELSE IF (ICOOR .EQ. 2) THEN
                        IMGLAB(ISCOOR) = 'MAGFLD.2'
                     ELSE
                        IMGLAB(ISCOOR) = 'MAGFLD.3'
                     END IF
                     CALL WRITDX(LUGDI,ISCOOR,IRAT*NVARPT,
     &                           WORK(KGD))
                     IF (SHIELD .OR. SPINRO) THEN
                        ISXYZ = IPTAX(ICOOR,2)
                        ANTI  = .TRUE.
                        CALL GETW(WMAT(1,ISXYZ),ANTI,CMO,WORK(KTDSQ),DV,
     &                            ICOOR,WORK(KWKRHS),LWKRHS,IPRMAX)
                     END IF
C
C                 Geometrical distortion
C                 ======================
C
                  ELSE
                     CALL WRITDX (LUGDR,ISCOOR,IRAT*NVARPT,
     &                            WORK(KGD))
#if defined (VAR_OLDGDHAMCODE)
C
C     This has not been corrected. Wrong dimensions on FDC and FDV
C     K.Ruud, june-94
C
                     IF (GDHAM) THEN
                        CALL WRINTD(ICRD,FDC(1,ISIM),FDV(1,ISIM),
     &                              QD(1,ISIM),DV,WRK(KFD),
     &                              H2DAC(IH2DAC(1,1,ISIM)),LUGDH)
                     END IF
#endif
C                    Check right-hand side
C                    950427: works only for symmetric connection
C
                     IF (NODIFC .AND. ISYM .EQ. 1 .AND. NCONST .GT. 0)
     &                  THEN
                        CALL CHKRHS(WORK(KTDSQ),EMYD(ICOOR),GRDACT,
     &                              FDC(1,ICOOR),WORK(KFDCSQ),
     &                              WORK(KH1AO),WORK(KH1MO),CMO,
     &                              WORK(KW2),LW2,IATOM,ICOOR,ISCOOR,
     &                              NOH1,NOH2,NOORTH,IPRGDY)
                     END IF
                  END IF
               END IF
               IF ((MCD .OR. VERDET) .AND. .NOT. NOLOND) THEN
                  CALL WRITDX(LUFDC,ISCOOR-3*NUCDEP,IRAT*N2ORBT,
     &                        WORK(KFDCSQ))
               END IF
C
C              ******************************
C              ***** Construct Y matrix *****
C              ******************************
C
               IF (DOSFD .AND. (SUSCEP .OR. IATOM.GT.0)) THEN
                  IF (NVARPT .EQ. 0 .AND.
     &                 NDSINT(ICOOR,ISYM - 1) .GT. 0) THEN
                     WRITE (LUPRI, '(//A)')
     &                    ' ERROR in MAKEGD for GETY: NVARPT .eq.'//
     &                    ' 0 and NDSINT .gt. 0 not implemented.'
                     CALL QUIT('ERROR in MAKEGD')
                  END IF
                  IX = 1
                  IF (SUSCEP) IX = ICOOR
                  CALL GETY(SYMAT(1,IX),WORK(KTDSQ),
     &                      FOCK,DV,FDC(1,ICOOR),FDV(1,ICOOR),
     &                      WORK(KFDCSQ),WORK(KFDVSQ),WORK(KQDI),
     &                      WORK(KQD),ISYM,IPRGDY,SUSCEP)
C
C                 Write differentiated overlap matrix and Y matrix
C
                  IF (.NOT.SUSCEP) THEN
                     CALL WRITDX(LUSFDA,2*ISCOOR - 1,
     &                           IRAT*N2ORBX,WORK(KTDSQ))
                     CALL WRITDX(LUSFDA,2*ISCOOR,
     &                           IRAT*N2ORBX,SYMAT(1,IX))
                  END IF
C
C     We need the inactive Fock matrix for later, so we write it to disk
C     K.Ruud, June 10 1996
C
               END IF
            END IF
         END IF
         IF (IATOM .EQ. 0) THEN
ckr            IDORCI(ISCOOR,2) = 1
         ELSE
ckr            IDORCI(ISCOOR,1) = 1
         END IF
         CALL ABAWRITRS
  300    CONTINUE
  100 CONTINUE
C
      IF (SUSCEP) THEN
         CALL DZERO(SUSFSY,9)
         DO 400 I = 1, 3
            DO 410 J = 1, 3
               SYELMT = DDOT(N2ORBX,TD(1,I),1,SYMAT(1,J),1)
C
C     Extra correction term, K.Ruud, May 1994
C
               IF (.NOT. NODIFC) THEN
                  CALL DZERO(WORK(KTDSQ),N2ORBX)
                  CALL DGEMM('T','N',NORBT,NORBT,NORBT,1.D0,
     &                       TD(1,I),NORBT,
     &                       TD(1,J),NORBT,0.D0,
     &                       WORK(KTDSQ),NORBT)
                  ISYMPT = MULD2H(ISYMAX(I,2) + 1,ISYMAX(J,2) + 1)
                  DO 116 ISYMP = 1, NSYM
                     ISYMQ = MULD2H(ISYMP,ISYMPT)
                     NORBQ = NORB(ISYMQ)
                     NOCCQ = NOCC(ISYMQ)
                     IF (NOCCQ .GT. 0) THEN
                        IOFFQ1 = IORB(ISYMQ)
                        IOFFRQ = I2ORB(ISYMQ)
                        DO 118 IQ = 1, NOCCQ
                           IQ1 = IOFFQ1 + IQ
                           TERM = D0
                           DO 119 IR = 1, NOCCQ
                              IR1 = (IQ1 - 1)*NORBT + IOFFQ1 + IR - 1
                              IRQ = IOFFRQ + (IR - 1)*NORBQ + IQ
                              TERM = TERM + WORK(KTDSQ+IR1)*FOCK(IRQ)
 119                          CONTINUE
                              SYELMT = SYELMT + TERM
 118                       CONTINUE
                     END IF
 116              CONTINUE
               END IF
               IF (IPRGDY .GT. 10) THEN
                  WRITE(LUPRI,'(1X,A,2I5,F10.5)')
     &                 'Lowest order reorthonormalization', I, J,SYELMT
               END IF
               IS = IPTAX(I,2)
               JS = IPTAX(J,2)
               SUSFSY(IS,JS) = SUSFSY(IS,JS) + SYELMT
               SUSFSY(JS,IS) = SUSFSY(JS,IS) + SYELMT
  410       CONTINUE
  400    CONTINUE
         IF (IPRALL .GT. 2) THEN
            CALL HEADER('Lowest-order reorthonormalization '
     &          //'contribution to susceptibilities',-1)
            CALL OUTPUT(SUSFSY,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
      CALL QEXIT('MAKEGD')
      RETURN
      END
C  /* Deck exfdtd */
      SUBROUTINE EXFDTD(TRMAT,SQMAT,ISYM,IPRINT)
C
C     Extract elements of correct symmetry from square matrix and
C     put them in square matrix.
C
C     131088 tuh, modified for quadratic matrices, june-94
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION TRMAT(NORBT,NORBT), SQMAT(NORBT,NORBT)
#include <inforb.h>
      CALL DZERO(SQMAT,N2ORBX)
      DO 100 ISYM1 = 1, NSYM
         IOFF1 = IORB(ISYM1)
         NORB1 = NORB(ISYM1)
         DO 200 ISYM2 = 1, NSYM
         IF (MULD2H(ISYM1,ISYM2) .EQ. ISYM) THEN
            IOFF2 = IORB(ISYM2)
            NORB2 = NORB(ISYM2)
            DO 300 I = IOFF1 + 1, IOFF1 + NORB1
            DO 300 J = IOFF2 + 1, IOFF2 + NORB2
               SQMAT(I,J) = TRMAT(I,J)
  300       CONTINUE
         END IF
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from EXFDTD',-1)
         WRITE (LUPRI, '(A,I5)') '  Symmetry of perturbation:',ISYM
         CALL OUTPUT(SQMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck wrintd */
      SUBROUTINE WRINTD(ICRD,FDC,FDV,QD,DV,FD,H2DAC,LUGDH)
C
C     1-Apr-1986 Hans Joergen Aa. Jensen
C
C     Write the derivative Hamiltonian and Fock matrices
C     needed for GD vector no. ICRD to LUGDH.
C
C     LUGDH must be positioned by calling program. The structure of
C     LUGDH from ABACUS is:
C
C       rec 1)   label GDHAMILT
C       rec 2)   NISHT,NASHT,NORBT,NCRD,(ICRD(I),I=1,NCRD)
C
C          This routine is supposed to be called NCRD times,
C          it writes:
C
C       rec i.1) FDC   (NNORBT,icrd)
C       rec i.2) FD    (NOCCT,NORBT,icrd)
C       rec i.3) H2DAC (NNASHX,NNASHX,icrd)
C
#include <implicit.h>
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0, D2 = 2.0D0 )
      DIMENSION FDC(NNORBT), FDV(NNORBT), QD(NASHT,NORBT),
     *          DV(*), FD(NOCCT,NORBT), H2DAC(*)
C
C Used from common blocks:
C  INFORB : NISHT,NASHT,NOCCT,NORBT,NNASHX,NNORBT
C  INFIND : IROW(*)
C
#include <inforb.h>
#include <infind.h>
C
C.....rec i.1)
C
      WRITE (LUGDH) FDC
C
C.....rec i.2)
C
      CALL DZERO(FD,(NOCCT*NORBT))
      DO 220 I = 1,NISHT
         DO 210 J = 1,I
            IJ = IROW(I) + J
            VAL = D2*(FDC(IJ) + FDV(IJ))
            FD(I,J) = VAL
            FD(J,I) = VAL
  210    CONTINUE
         DO 215 J = NISHT+1,NORBT
            IJ = IROW(J) + I
            FD(I,J) = D2*(FDC(IJ) + FDV(IJ))
  215    CONTINUE
  220 CONTINUE
      DO 240 I = 1,NASHT
         DO 230 J = 1,NORBT
            VAL = D0
            DO 225 K = 1,NASHT
               IF (K .GT. I) THEN
                  IK = IROW(K) + I
               ELSE
                  IK = IROW(I) + K
               END IF
               IF (K .GT. J) THEN
                  KJ = IROW(K) + J
               ELSE
                  KJ = IROW(J) + K
               END IF
               VAL = VAL + DV(IK)*FDC(KJ)
  225       CONTINUE
            FD(NISHT+I,J) = VAL + QD(I,J)
  230    CONTINUE
  240 CONTINUE
      WRITE (LUGDH) FD
C
C.....rec i.3)
C
      IF (NASHT .GT. 0) THEN
         NH2DAC = NNASHX*NNASHX
         CALL DSCAL(NH2DAC,D2,H2DAC,1)
         JAB = -NNASHX
         DO 350 IA = 1,NASHT
            DO 345 IB = 1,IA
               JAB = JAB + NNASHX
               IF (IA .EQ. IB) THEN
                  CALL DSCAL(NNASHX,DP5,H2DAC(JAB+1),1)
               ELSE
                  ICC = JAB
                  DO 330 IC = 1,NASHT
                     ICC = ICC + IC
                     H2DAC(ICC) = DP5*H2DAC(ICC)
  330             CONTINUE
               END IF
  345       CONTINUE
  350    CONTINUE
         WRITE (LUGDH) (H2DAC(I),I=1,NH2DAC)
      ELSE
         WRITE (LUGDH) D0,D0,D0,D0
      END IF
C
C     End of WRINTD.
C
      RETURN
      END
C  /* Deck tr1tst */
#if defined (VAR_TESTIO)
      SUBROUTINE TR1TST(ICOOR,IPRINT,CMO,DV,TD,F,FSC,FSV,QS,F1,
     *                  H1,H1AO,FS,WRK,LWRK)
C
C     This subroutine tests the one-index transformed Fock matrices
C     and Q matrix by calculating the reorthonormalization term to the
C     molecular gradient (double weight on two-electron terms).
C
#include <implicit.h>
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D2 = 2.0D0)
C
      INTEGER P, Q, T, U, TU, QU
      DIMENSION F1(NOCCT,*), FSC(*), FSV(*), FS(NOCCT,*),
     *          QS(NASHDI,*), H1AO(*), TD(NORBT,NORBT,*),
     *          CMO(*), DV(*), F(NORBT,*), H1(*), WRK(LWRK)
C
C Used from common blocks:
C  INFDIM : NASHDI,?
C
#include <inftap.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
C
      ITRI(I) = I*(I - 1)/2
C
C     *****************************************************
C     ***** Undifferentiated one-electron Hamiltonian *****
C     *****************************************************
C
      CALL RDONEL('ONEHAMIL',.TRUE.,H1AO,NNBAST)
      CALL UTHU(H1AO(1),H1(1),CMO(1),WRK(1),NBAST,NORBT)
C
C     *****************************************
C     ***** Total Fock matrices F1 and FS *****
C     *****************************************
C
C     Inactive block
C
      DO 200 I = 1, NISHT
         DO 210 Q = 1, NORBT
            IF (I .GE. Q) THEN
               IQ = ITRI(I) + Q
            ELSE
               IQ = ITRI(Q) + I
            END IF
            F1(I,Q) = D2*H1(IQ)
            FS(I,Q) = D2*(FSC(IQ) + FSV(IQ))
  210    CONTINUE
  200 CONTINUE
C
C     Active block
C
      DO 300 T = 1, NASHT
         DO 310 Q = 1, NORBT
            F1TQ = D0
            FSTQ = D0
            DO 320 U = 1, NASHT
               IU = NISHT + U
               IF (Q .GE. IU) THEN
                  QU = ITRI(Q) + IU
               ELSE
                  QU = ITRI(IU) + Q
               END IF
               IF (T .GE. U) THEN
                  TU = ITRI(T) + U
               ELSE
                  TU = ITRI(U) + T
               END IF
               F1TQ = F1TQ + DV(TU)*H1(QU)
               FSTQ = FSTQ + DV(TU)*FSC(QU)
  320       CONTINUE
            F1(NISHT + T,Q) = F1TQ
            FS(NISHT + T,Q) = FSTQ + QS(T,Q)
  310    CONTINUE
  300 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('FS MATRIX',-1)
         CALL OUTPUT(FS(1,1),1,NOCCT,1,NORBT,NOCCT,NORBT,1,LUPRI)
      END IF
C
C     **************************************
C     ***** Reorthonormalization terms *****
C     **************************************
C
C     A: Using undifferentiated Fock matrix
C
      TRACE1 = D0
      DO 400 I = 1, NOCCT
         DO 410 J = 1, I
            F1IJ = (F1(I,J) + F1(J,I))*DP5
            TERM = (F(I,J) + F(J,I) - F1IJ)*TD(I,J,ICOOR)
            IF (I .NE. J) TERM = TERM + TERM
            TRACE1 = TRACE1 + TERM
  410    CONTINUE
         DO 420 J = NOCCT + 1, NORBT
            TRACE1 = TRACE1 - F1(I,J)*TD(I,J,ICOOR)
  420    CONTINUE
  400 CONTINUE
      TRACE1 = TRACE1 + TRACE1
C
C     B: Using one-index transformed Fock matrix
C
      TRACE2 = D0
      DO 500 P = 1, NOCCT
         TRACE2 = TRACE2 + FS(P,P)
  500 CONTINUE
C
C     ***************************
C     ***** Compare results *****
C     ***************************
C
      WRITE (LUPRI,'(//,3(A/))')
     *      ' Test on One-Index Transformed Fock Matrices and Q Matrix',
     *      ' --------------------------------------------------------',
     *      '          Trace(SD*F)    Trace(FS)      Difference'
      WRITE (LUPRI,'(5X,3F15.10)') TRACE1, TRACE2, TRACE1 - TRACE2
      RETURN
      END
#endif
C  /* Deck gety */
      SUBROUTINE GETY(YMAT,TD,FOCK,DV,FDC,FDV,FSC,FSV,QD,QS,
     *                ISYMPT,IPRINT,SUSCEP)
C
C     Purpose: Construction of YMAT matrix
C
C     tuh 281088
C     Making Y-mat and TD square, jan.-94, K.Ruud
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxash.h>
#include <maxorb.h>
#include <abainf.h>
      PARAMETER(D0=0.0D0, DP5=0.5D0, D1 = 1.0D0, D2 = 2.0D0)
C
      LOGICAL NOORBI,SUSCEP
      INTEGER P, P1, Q, Q1, Q2, R, R1, U, U1, U2, V, V1, V2, PQ, UQ, VQ,
     &        UV
      DIMENSION YMAT(NORBT,NORBT), TD(NORBT,NORBT), FOCK(N2ORBT),
     *          DV(NNASHX),FDC(NORBT,NORBT),FDV(NORBT,NORBT),
     *          FSC(NORBT,NORBT),FSV(NORBT,NORBT),
     *          QD(NORBT,NASHT),QS(NORBT,NASHT)
C
C     Used from common blocks :
C      INFORB : NORBT,NNORBT,NASHT,...
C      INFIND : IROW(*)
C
#include <inforb.h>
#include <infind.h>
C
      IF (IPRINT .GE. 5) THEN
         CALL TITLER('Output from GETY','*',103)
         WRITE (LUPRI,'(/A,I5)') ' ISYMPT ', ISYMPT
         IF (IPRINT .GE. 10) THEN
            CALL AROUND('Fock matrix in GETY')
            DO 50 ISYM = 1, NSYM
               WRITE (LUPRI,'(A,I5)') ' Irrep', ISYM
               NORBI = NORB(ISYM)
               ISTRI = I2ORB(ISYM) + 1
               CALL OUTPUT(FOCK(ISTRI),1,NORBI,1,NORBI,
     *                     NORBI,NORBI,1,LUPRI)
   50       CONTINUE
            CALL AROUND('DV matrix in GETY')
            CALL OUTPAK(DV,NASHT,1,LUPRI)
            CALL AROUND('FDC matrix in GETY')
            CALL OUTPUT(FDC,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL AROUND('FDV matrix in GETY')
            CALL OUTPUT(FDV,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL AROUND('FSC matrix in GETY')
            CALL OUTPUT(FSC,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL AROUND('FSV matrix in GETY')
            CALL OUTPUT(FSV,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL AROUND('QD matrix in GETY')
            CALL OUTPUT(QD,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            CALL AROUND('QS matrix in GETY')
            CALL OUTPUT(QS,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            CALL AROUND('TD matrix in GETY')
            CALL OUTPUT(TD,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
C
      IF (SUSCEP) THEN
         SUSSGN = D1
      ELSE
         SUSSGN = - D1
      END IF
C
C     ********************
C     ***** Y matrix *****
C     ********************
C
      CALL DZERO(YMAT,N2ORBX)
C
      DO 100 ISYMI = 1, NSYM
         ISYMQ = MULD2H(ISYMI,ISYMPT)
         NISHI = NISH(ISYMI)
         NORBQ = NORB(ISYMQ)
         IF (NISHI.GT.0 .AND. NORBQ.GT.0) THEN
            IOFFI = IORB(ISYMI)
            IOFFQ = IORB(ISYMQ)
            DO 110 I = IOFFI + 1, IOFFI + NISHI
               DO 120 Q = IOFFQ + 1, IOFFQ + NORBQ
                  YMAT(Q,I) = YMAT(Q,I) - (FDC(I,Q) + FSC(I,Q))*D2
                  IF (NASHT .GT. 0) THEN
                     YMAT(Q,I) = YMAT(Q,I) - (FDV(I,Q) + FSV(I,Q))*D2
                  END IF
 120           CONTINUE
 110        CONTINUE
         END IF
 100  CONTINUE
C
      DO 200 ISYMU = 1, NSYM
         ISYMQ = MULD2H(ISYMU,ISYMPT)
         NASHU = NASH(ISYMU)
         NORBQ = NORB(ISYMQ)
         IF (NASHU.GT.0 .AND. NORBQ.GT.0) THEN
            IOFFU1 = IORB(ISYMU) + NISH(ISYMU)
            IOFFU2 = IASH(ISYMU)
            IOFFQ  = IORB(ISYMQ)
            DO 210 U = 1, NASHU
               U1 = IOFFU1 + U
               U2 = IOFFU2 + U
               DO 220 Q = IOFFQ + 1, IOFFQ + NORBQ
                  TERM = D0
                  DO 230 V = 1, NASHU
                     V1 = IOFFU1 + V
                     V2 = IOFFU2 + V
                     IF (U2 .GT. V2) THEN
                        UV = IROW(U2) + V2
                     ELSE
                        UV = IROW(V2) + U2
                     END IF
                     TERM = TERM - DV(UV)*(FSC(V1,Q) + FDC(V1,Q))
 230              CONTINUE
                  YMAT(Q,U1) = YMAT(Q,U1)+TERM+
     &               SUSSGN*(QD(Q,U2)+QS(Q,U2))
 220           CONTINUE
 210        CONTINUE
         END IF
 200  CONTINUE
C
      DO 300 ISYMP = 1, NSYM
         ISYMQ = MULD2H(ISYMP,ISYMPT)
         NORBP = NORB(ISYMP)
         NORBQ = NORB(ISYMQ)
         NOCCQ = NOCC(ISYMQ)
         IF (NORBP.GT.0 .AND. NOCCQ.GT.0) THEN
            IOFFP1 = IORB(ISYMP)
            IOFFQ1 = IORB(ISYMQ)
            IOFFRQ = I2ORB(ISYMQ)
            DO 310 P = 1, NORBP
               P1 = IOFFP1 + P
               DO 320 Q = 1, NOCCQ
                  Q1 = IOFFQ1 + Q
                  TERM = D0
                  DO 330 R = 1, NOCCQ
                     R1  = IOFFQ1 + R
                     IRQ = IOFFRQ + (R - 1)*NORBQ + Q
                     TERM = TERM + TD(P1,R1)*FOCK(IRQ)
 330              CONTINUE
                  IF (NODIFC) THEN
                     YMAT(Q1,P1) = YMAT(Q1,P1) - D2*TERM
                  ELSE
                     YMAT(Q1,P1) = YMAT(Q1,P1) - TERM
                  END IF
 320           CONTINUE
 310        CONTINUE
         END IF
 300  CONTINUE
      IF (IPRINT .GE. 7) THEN
         CALL AROUND('Y matrix in GETY')
         CALL OUTPUT(YMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck chkrhs */
      SUBROUTINE CHKRHS(TD,EMYD,GRDACT,FDC,FSC,H1AO,H1MO,CMO,WRK,LWRK,
     *                  IATOM,ICOOR,ISCOOR,NOH1,NOH2,NOORTH,IPRINT)
C
C     This subroutine retrieves the molecular gradient from the
C     configuration part of the right-hand side of the response
C     equations.
C
C     EMYD is calculated in GETFD
C     GRDACT is calculated in ABARHS
C
C     tuh 101188
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxash.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D2 = 2.0D0)
      PARAMETER (THRS = 1.0D-10)
      LOGICAL NOH1, NOH2, NOORTH
      DIMENSION TD(NORBT,NORBT), FSC(NORBT,NORBT), FDC(NORBT,NORBT),
     *          H1AO(NNBAST), H1MO(NNORBT),
     *          CMO(*), WRK(LWRK)
#include <abainf.h>
#include <inforb.h>
#include <infdim.h>
#include <infind.h>
#include <inftap.h>
#include <energy.h>
C
      IF (IPRINT .GE. 5) THEN
         CALL TITLER('Output from CHKRHS','*',103)
      END IF
      GRDNUC = GRADNN(ISCOOR)
C
C     *******************************************
C     ***** Get one-electron MO Hamiltonian *****
C     *******************************************
C
C     Read in one-electron AO Hamiltonian
C
      CALL RDONEL('ONEHAMIL',.TRUE.,H1AO,NNBAST)
C
C     Transform to MO basis
C
      DO 200 ISYM = 1, NSYM
         NORBI = NORB(ISYM)
         IF (NORBI .GT. 0) THEN
            NBASI = NBAS(ISYM)
            IJ    = IIORB(ISYM)
            CALL UTHU(H1AO(IIBAS(ISYM)+1),H1MO(IJ+1),
     *                CMO(ICMO(ISYM)+1),WRK,NBASI,NORBI)
         END IF
  200 CONTINUE
      IF (IPRINT .GE. 20) THEN
        CALL AROUND('One-electron Hamiltonian in MO basis - CHKRHS')
        CALL OUTPKB(H1MO,NORB,NSYM,1,LUPRI)
      END IF
C
C     ***************************************************************
C     ***** Calculate reorthonormalization of inactive gradient *****
C     ***************************************************************
C
      EMYS = D0
      IF (.NOT.NOORTH) THEN
         DO 400 ISYM = 1, NSYM
            NISHI  = NISH(ISYM)
            NORBI  = NORB(ISYM)
            IORBI  = IORB(ISYM)
            IOFFIJ = IIORB(ISYM) - 2*IORBI
            IF (NISHI .GT. 0) THEN
               DO 500 I = IORBI + 1, IORBI + NISHI
                  II = IROW(I) + I
                  EMYS = EMYS + FSC(I,I) - FDC(I,I)
                  DO 600 J = IORBI + 1, IORBI + NORBI
                     IJMAX = MAX(I,J) - IORBI
                     IJ    = IOFFIJ + (IJMAX*(IJMAX - 3))/2 + I + J
                     EMYS  = EMYS + D2*TD(I,J)*H1MO(IJ)
  600             CONTINUE
  500          CONTINUE
            END IF
  400    CONTINUE
      END IF
      IF (IPRINT .GE. 5) THEN
         WRITE (LUPRI,'(A,F20.10)') ' EMYD ', EMYD
         WRITE (LUPRI,'(A,F20.10)') ' EMYS ', EMYS
      END IF
C
C     *********************************
C     ***** Calculate CI gradient *****
C     *********************************
C
      GRDICT = EMYD + EMYS
      GRADCI = GRDNUC + GRDACT + GRDICT
C
C     *********************************
C     ***** Calculate AO gradient *****
C     *********************************
C
      GRADAO = GRDNUC
      IF (.NOT.NOH1)   GRADAO = GRADAO + GRADKE(ISCOOR) + GRADNA(ISCOOR)
      IF (.NOT.NOH2)   GRADAO = GRADAO + GRADEE(ISCOOR)
      IF (.NOT.NOORTH) GRADAO = GRADAO + GRADFS(ISCOOR)
C
C     *******************
C     ***** Compare *****
C     *******************
C
      DIFFER = GRADCI - GRADAO
      ABSDIF = ABS(DIFFER)
      IF (ABSDIF .GE. THRS) THEN
         CALL HEADER('WARNING: Dif. between mol. CI and AO'//
     *               ' gradients greater than threshold!',0)
         NWNABA = NWNABA + 1
      END IF
      IF ((ABSDIF .GE. THRS) .OR. (IPRINT .GE. 4)) THEN
         WRITE (LUPRI,'(/A,I3,A,I2,A,/A/)')
     *          ' Molecular gradient for atom', IATOM,
     *          ' and coordinate',ICOOR,':',
     *          ' ------------------------------------------------'
         WRITE (LUPRI,'(4X,2A,/)') '    Active(CI)      Inactive',
     *                             '       Nuclear      Gradient'
         WRITE (LUPRI,'(4X,4F14.10)') GRDACT,GRDICT,GRDNUC,GRADCI
         WRITE (LUPRI,'(/,2A,/)')
     *                       '         CI gradient',
     *                       '         AO gradient          Difference'
         WRITE (LUPRI,'(3F20.10)') GRADCI, GRADAO, DIFFER
         IF (NOH1)   WRITE (LUPRI,'(/,4X,2A,/)') ' Note: No',
     *      ' one-electron terms included in this test!'
         IF (NOH2)   WRITE (LUPRI,'(/,4X,2A,/)') ' Note: No',
     *          ' two-electron terms included in this test!'
         IF (NOORTH) WRITE (LUPRI,'(/,4X,2A,/)') ' Note: No',
     *          ' reorthonormalization included in this test!'
      END IF
      RETURN
      END
C  /* Deck ofachk */
      SUBROUTINE OFACHK(RUNSDR,NOH1,NOH2,NOORTH,WORK,LWORK,IPRINT)
C
C     NODDY SUBROUTINE TO FORM ORBITAL PART OF RIGHT HAND SIDE
C     FOR MAGNETIC RESPONSE EQUATIONS (Eq. 51b JCP 84, 6266 (1986))
C     USING THE GENERALISED FOCK MATRIX
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <dummy.h>
C
      LOGICAL RUNSDR, NOH1, NOH2, NOORTH
      DIMENSION WORK(LWORK)
C
#include <nuclei.h>
#include <inforb.h>
C
      IF (NSYM .GT. 1) THEN
         WRITE (LUPRI,'(1X,A)') ' OFACHK only works without symmetry.'
         CALL QUIT('OFACHK does not work with symmetry,')
      END IF
C
C     Allocate core memory
C
      KUMOG  = 1
      KUMOGA = KUMOG  + N2ORBX*N2ORBX
      KOMOGA = KUMOGA + N2ORBX*N2ORBX
      KPDENS = KOMOGA + N2ORBX*N2ORBX
      KDDENS = KPDENS + N2ORBX*N2ORBX
      KDV    = KDDENS + N2ORBX
      KDFCKU = KDV    + N2ASHX
      KDFCKA = KDFCKU + N2ORBX*3
      KOFOCK = KDFCKA + N2ORBX*3
      KSXMAT = KOFOCK + N2ORBX*3
      KFMAT  = KSXMAT + N2ORBX*3
      KYMAT  = KFMAT  + N2ORBX*3
      KCMO   = KYMAT  + N2ORBX*3
      KUMOH  = KCMO   + N2ORBX
      KUMOHA = KUMOH  + N2ORBX
      KOMOHA = KUMOHA + N2ORBX
      KHBMAT = KOMOHA + N2ORBX
      KGWORK = KHBMAT + N2ORBX*3
      KHWORK = KGWORK + N2ORBX*N2ORBX
      KCFKUA = KHWORK + N2ORBX
      KAFKUA = KCFKUA + N2ORBX*3
      KQMTUA = KAFKUA + N2ORBX*3
      KHWRK3 = KQMTUA + N2ORBX*3
      KPSOAO = KHWRK3 + N2ORBX*3
      KPSOMO = KPSOAO + 3*NUCDEP*N2BASX
      KTOP   = KPSOMO + N2ORBX
      IF (KTOP .GT. LWORK) CALL STOPIT('OFACHK',' ',KTOP,LWORK)
C
      LWRK = LWORK - KTOP + 1
C
      CALL OFABOD(WORK(KUMOG),WORK(KUMOGA),WORK(KOMOGA),WORK(KPDENS),
     &            WORK(KDDENS),WORK(KDV),WORK(KDFCKU),WORK(KDFCKA),
     &            WORK(KOFOCK),WORK(KSXMAT),
     &            WORK(KFMAT),WORK(KYMAT),WORK(KCMO),WORK(KUMOH),
     &            WORK(KUMOHA),WORK(KOMOHA),WORK(KGWORK),WORK(KHWORK),
     &            WORK(KCFKUA),WORK(KAFKUA),WORK(KQMTUA),WORK(KHWRK3),
     &            WORK(KHBMAT),WORK(KPSOAO),WORK(KPSOMO),WORK(KTOP),
     &            LWRK,RUNSDR,NOH1,NOH2,NOORTH,IPRINT)
      RETURN
      END
C  /* Deck ofabod */
      SUBROUTINE OFABOD(UMOG,UMOGA,OMOGA,PDENS,DDENS,DV,DFOCKU,DFOCK,
     &                  OFOCK,SXMAT,FMAT,YMAT,CMO,UMOH,UMOHA,OMOHA,
     &                  GWORK,HWORK,CFKUA,AFKUA,QMTUA,HWRK3,HBMAT,
     &                  PSOAO,PSOMO,WORK,LWORK,RUNSDR,NOH1,NOH2,NOORTH,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (DP5 = 0.5D0, D0 = 0.0D0)
      INTEGER TMAT
      LOGICAL RUNSDR, NOH1, NOH2, NOORTH
      INTEGER P, Q, R, S, A, B, C, D, T, U, V, W, X
      DIMENSION  UMOG (NORBT,NORBT,NORBT,NORBT),
     &           UMOGA(NORBT,NORBT,NORBT,NORBT),
     &           OMOGA(NORBT,NORBT,NORBT,NORBT),
     &           UMOH (NORBT,NORBT),
     &           UMOHA(NORBT,NORBT),
     &           OMOHA(NORBT,NORBT),
     &           PDENS(NORBT,NORBT,NORBT,NORBT), DDENS(NORBT,NORBT),
     &           DV(NASHT,NASHT),
     &           DFOCKU(NORBT,NORBT,3),
     &           DFOCK(NORBT,NORBT,3),
     &           OFOCK(NORBT,NORBT,3),
     &           SXMAT(NORBT,NORBT,3),
     &           FMAT (NORBT,NORBT),
     &           YMAT (NORBT,NORBT,3),
     &           CMO  (NBAST,NORBT),
     &           CFKUA(NORBT,NORBT,3),
     &           AFKUA(NORBT,NORBT,3),
     &           QMTUA(NORBT,NORBT,3),
     &           HWRK3(NBAST,NBAST,3),
     &           HWORK(NBAST,NBAST), HBMAT(NBAST,NBAST,3),
     &           GWORK(NBAST,NBAST,NBAST,NBAST),
     &           PSOAO(NBAST,NBAST,*), PSOMO(NORBT,NORBT), WORK(LWORK)
#include <abainf.h>
#include <infinp.h>
#include <nuclei.h>
#include <inforb.h>
#include <efield.h>
#include <chrxyz.h>
      NT(T) = NISHT + T
C
C     Read in MO's
C     ============
C
      CALL OFABMO(CMO,DV,FMAT,WORK,LWORK)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('MO coefficients',-1)
         CALL OUTPUT(CMO,1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
         CALL HEADER('DV matrix',-1)
         CALL OUTPUT(DV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         CALL HEADER('Fock matrix',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Read in derivative overlap matrix
C     =================================
C
C     Read derivative overlap matrix into HWRK3
C     Transform to MO basis in SXMAT
C
      IF (RUNSDR .AND. .NOT.NOLOND) THEN
         NCOMP = 3
         CALL GET1PR(HWRK3,'S1MAG  ',NCOMP,'SQUARE',.TRUE.,WORK,
     &               LWORK,IDUMMY,5)
      ELSE
         CALL DZERO(HWRK3,3*N2ORBX)
      END IF
      CALL DZERO(SXMAT,3*N2BASX)
      DO 50 N=1,3
      DO 51 P=1,NORBT
      DO 51 Q=1,NORBT
      DO 51 A=1,NBAST
      DO 51 B=1,NBAST
51       SXMAT(P,Q,N)=SXMAT(P,Q,N)+HWRK3(A,B,N)*CMO(A,P)*CMO(B,Q)
         IF (IPRINT .GT. 5) THEN
           CALL HEADER(CHRXYZ(N)//' der. overlap matrix (AO basis)',-1)
           CALL OUTPUT(HWRK3(1,1,N),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
           CALL HEADER(CHRXYZ(N)//' der. overlap matrix (MO basis)',-1)
           CALL OUTPUT(SXMAT(1,1,N),1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
         END IF
50    CONTINUE
C
C
C     Read in density matrices
C     ========================
C
C     Density matrices read into DDENS and PDENS
C
      CALL GETDNS(DDENS,PDENS)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('One-electron density matrix',-1)
         CALL OUTPUT(DDENS,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Two-electron density matrix',-1)
         CALL OUTPUT(PDENS,1,N2ORBX,1,N2ORBX,N2ORBX,N2ORBX,1,LUPRI)
      END IF
C
C     Reorthonormalization terms
C     ==========================
C
      NCOMP = 3*NUCDEP
      CALL GET1PR(PSOAO,'PSO    ',NCOMP,'SQUARE',.TRUE.,WORK,LWORK,
     &            IDUMMY,5)
      DO 600 N = 1, 3*NUCDEP
         CALL DZERO(PSOMO,N2ORBX)
         DO 610 P = 1, NORBT
         DO 610 Q = 1, NORBT
         DO 610 A = 1, NBAST
         DO 610 B = 1, NBAST
            PSOMO(P,Q)=PSOMO(P,Q)+PSOAO(A,B,N)*CMO(A,P)*CMO(B,Q)
  610    CONTINUE
         DO 620 X = 1, 3
            ORTO = 0.0D0
            DO 630 P = 1, NORBT
            DO 630 Q = 1, NORBT
            DO 630 R = 1, NORBT
               ORTO = ORTO - DDENS(P,Q)*SXMAT(Q,R,X)*PSOMO(R,P)
  630       CONTINUE
            WRITE (LUPRI,'(A,F12.6,2I5)')
     &         ' Shielding reorthonormalization in OFABOD:',ORTO,X,N
  620    CONTINUE
  600 CONTINUE
C
C     *****************************************************************
C     ***** Read in unmodified AO one- and two-electron integrals *****
C     *****             and transform to MO basis                 *****
C     *****************************************************************
C
C     Read into HWORK and UMOG
C
      CALL OFABHG(HWORK,UMOG,WORK,LWORK)
C
C     Add contribution from electric field, KRu, Aug.-93
C
      IF (NFIELD .GT. 0) THEN
         TMAT  = 1
         LWRK  = TMAT + N2ORBX*3
         IF (LWRK .GT. LWORK) CALL STOPIT('OFABOD',' ',LWRK,LWORK)
         LLEFT = LWORK - LWRK + 1
         DO 555 IFIELD = 1, NFIELD
            IF (EFIELD(IFIELD) .NE. D0) THEN
               IF (LFIELD(IFIELD) .EQ. 'XDIPLEN ') THEN
                  IOFF = 1
               ELSE IF (LFIELD(IFIELD) .EQ. 'YDIPLEN ') THEN
                  IOFF = 1 + N2ORBX
               ELSE IF (LFIELD(IFIELD) .EQ. 'ZDIPLEN ') THEN
                  IOFF = 1 + 2*N2ORBX
               ELSE
                  WRITE (LUPRI,'(/,3A,/)') 'Field type ',
     &                  LFIELD(IFIELD),
     &                  ' not implemented for magnetizabileties'
                  CALL QUIT('Illegal field type for magnetizabileties')
               END IF
               NCOMP = 3
               CALL GET1PR(WORK(TMAT),'DIPLEN ',NCOMP,'SQUARE',
     &                     .TRUE.,WORK(LWRK),LLEFT,IDUMMY,5)
               CALL DAXPY(N2ORBX,EFIELD(IFIELD),WORK(IOFF),1,HWORK,1)
            END IF
 555     CONTINUE
      END IF
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Undifferentiated one-electron AO Hamiltonian',-1)
         CALL OUTPUT(HWORK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL DZERO(UMOH,N2ORBX)
      DO 10 P=1,NORBT
      DO 10 Q=1,NORBT
      DO 10 A=1,NBAST
      DO 10 B=1,NBAST
         UMOH(P,Q)=UMOH(P,Q)+HWORK(A,B)*CMO(A,P)*CMO(B,Q)
10    CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Undifferentiated one-electron MO Hamiltonian',-1)
         CALL OUTPUT(UMOH,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL DZERO(GWORK,N2ORBX*N2ORBX)
      DO 20 P=1,NORBT
      DO 20 A=1,NBAST
      DO 20 B=1,NBAST
      DO 20 C=1,NBAST
      DO 20 D=1,NBAST
20       GWORK(P,B,C,D)=GWORK(P,B,C,D)+UMOG(A,B,C,D)*CMO(A,P)
      CALL DZERO(UMOG,N2ORBX*N2ORBX)
      DO 22 P=1,NORBT
      DO 22 Q=1,NORBT
      DO 22 B=1,NBAST
      DO 22 C=1,NBAST
      DO 22 D=1,NBAST
22       UMOG(P,Q,C,D)=UMOG(P,Q,C,D)+GWORK(P,B,C,D)*CMO(B,Q)
      CALL DZERO(GWORK,N2ORBX*N2ORBX)
      DO 24 P=1,NORBT
      DO 24 Q=1,NORBT
      DO 24 R=1,NORBT
      DO 24 C=1,NBAST
      DO 24 D=1,NBAST
24       GWORK(P,Q,R,D)=GWORK(P,Q,R,D)+UMOG(P,Q,C,D)*CMO(C,R)
      CALL DZERO(UMOG,N2ORBX*N2ORBX)
      DO 26 P=1,NORBT
      DO 26 Q=1,NORBT
      DO 26 R=1,NORBT
      DO 26 S=1,NORBT
      DO 26 D=1,NBAST
26       UMOG(P,Q,R,S)=UMOG(P,Q,R,S)+GWORK(P,Q,R,D)*CMO(D,S)
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Undifferentiated two-electron MO Hamiltonian',-1)
         CALL OUTPUT(UMOG,1,N2ORBX,1,N2ORBX,N2ORBX,N2ORBX,1,LUPRI)
      END IF
C
C     *************************************************************
C     ***** Loop over symmetries and perturbation coordinates *****
C     *************************************************************
C
      NCOMP = 3
      IF (NOLOND) THEN
         CALL GET1PR(HBMAT,'ANGMOM ',NCOMP,'SQUARE',.TRUE.,
     &               WORK,LWORK,IDUMMY,5)
         CALL DSCAL(3*N2ORBX,DP5,HBMAT,1)
      ELSE
         CALL GET1PR(HBMAT,'MAGMOM ',NCOMP,'SQUARE',.TRUE.,
     &               WORK,LWORK,IDUMMY,5)
      END IF
C
C     Add contribution from electric field, KRu, Aug.-93
C
      IF (NFIELD .GT. 0) THEN
         TMAT  = 1
         LWRK  = TMAT + N2ORBX*3
         IF (LWRK .GT. LWORK) CALL STOPIT('OFABOD',' ',LWRK,LWORK)
         LLEFT = LWORK - LWRK + 1
         DO 556 IFIELD = 1, NFIELD
            IF (EFIELD(IFIELD) .NE. D0) THEN
               IF  (LFIELD(IFIELD) .EQ. 'XDIPLEN ') THEN
                  FIELD1 = 'X-FIELD'
               ELSE IF (LFIELD(IFIELD) .EQ. 'YDIPLEN ') THEN
                  FIELD1 = 'Y-FIELD'
               ELSE IF (LFIELD(IFIELD) .EQ. 'ZDIPLEN ') THEN
                  FIELD1 = 'Z-FIELD'
               ELSE
                  WRITE (LUPRI,'(/,3A,/)') 'Field type ',
     &                  LFIELD(IFIELD),
     &                  ' not implemented for magnetizabileties'
                  CALL QUIT('Illegal field type for magnetizabileties')
               END IF
               NCOMP = 3
               CALL GET1PR(WORK(TMAT),'CM1    ',NCOMP,'SQUARE',
     &                     .TRUE.,WORK(LWRK),LLEFT,IDUMMY,5)
               CALL DAXPY(N2ORBX*3,EFIELD(IFIELD),WORK(TMAT),1,HBMAT,1)
            END IF
 556     CONTINUE
      END IF
C
      IF (NOH1) CALL DZERO(HBMAT,3*N2ORBX)
      CALL DZERO(DFOCKU,3*N2ORBX)
      CALL DZERO(DFOCK,3*N2ORBX)
      CALL DZERO(CFKUA,3*N2ORBX)
      CALL DZERO(AFKUA,3*N2ORBX)
      CALL DZERO(QMTUA,3*N2ORBX)
      DO 3 X = 1,3
         WRITE(LUPRI,'(//A,I3,/A)') ' Cartesian direction ',X,
     *                             ' ======================= '
      CALL DZERO(UMOHA,N2ORBX)
      DO 30 P=1,NORBT
      DO 30 Q=1,NORBT
      DO 30 A=1,NBAST
      DO 30 B=1,NBAST
         UMOHA(P,Q)=UMOHA(P,Q)+HBMAT(A,B,X)*CMO(A,P)*CMO(B,Q)
30    CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Differentiated one-electron AO Hamiltonian',-1)
         CALL OUTPUT(HBMAT(1,1,X),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Differentiated one-electron MO Hamiltonian',-1)
         CALL OUTPUT(UMOHA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL DZERO(GWORK,N2ORBX*N2ORBX)
C
C     Read UMOGA
C
      CALL RUMOGA(UMOGA,NORBT,NOH2,NOLOND,WORK,LWORK,X)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Differentiated two-electron AO Hamiltonian',-1)
         CALL OUTPUT(UMOGA,1,N2ORBX,1,N2ORBX,N2ORBX,N2ORBX,1,LUPRI)
      END IF
C
      DO 40 D=1,NBAST
      DO 40 C=1,NBAST
      DO 40 B=1,NBAST
      DO 40 P=1,NORBT
      DO 40 A=1,NBAST
         GWORK(P,B,C,D)=GWORK(P,B,C,D)+UMOGA(A,B,C,D)*CMO(A,P)
   40 CONTINUE
C
      CALL DZERO(UMOGA,N2ORBX*N2ORBX)
      DO 42 D=1,NBAST
      DO 42 C=1,NBAST
      DO 42 Q=1,NORBT
      DO 42 B=1,NBAST
      DO 42 P=1,NORBT
         UMOGA(P,Q,C,D)=UMOGA(P,Q,C,D)+GWORK(P,B,C,D)*CMO(B,Q)
   42 CONTINUE
C
      CALL DZERO(GWORK,N2ORBX*N2ORBX)
      DO 44 D=1,NBAST
      DO 44 R=1,NORBT
      DO 44 C=1,NBAST
      DO 44 Q=1,NORBT
      DO 44 P=1,NORBT
         GWORK(P,Q,R,D)=GWORK(P,Q,R,D)+UMOGA(P,Q,C,D)*CMO(C,R)
   44 CONTINUE
C
      CALL DZERO(UMOGA,N2ORBX*N2ORBX)
      DO 46 D=1,NBAST
      DO 46 S=1,NORBT
      DO 46 R=1,NORBT
      DO 46 Q=1,NORBT
      DO 46 P=1,NORBT
         UMOGA(P,Q,R,S)=UMOGA(P,Q,R,S)+GWORK(P,Q,R,D)*CMO(D,S)
   46 CONTINUE
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Differentiated two-electron MO Hamiltonian',-1)
         CALL OUTPUT(UMOGA,1,N2ORBX,1,N2ORBX,N2ORBX,N2ORBX,1,LUPRI)
      END IF
C
C     ********************************
C     ***** Form {S,h} and {S,g} *****
C     ********************************
C
C     Form {S,h} in HWORK and {S,g} into GWORK
C
      CALL DZERO(HWORK,N2ORBX)
      DO 60 Q=1,NORBT
      DO 60 A=1,NORBT
      DO 60 P=1,NORBT
         HWORK(P,Q)=HWORK(P,Q)-0.5D0*(SXMAT(P,A,X)*UMOH(A,Q)
     +                               -SXMAT(Q,A,X)*UMOH(P,A))
60    CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('One-electron reorthonormalisation term',-1)
         CALL OUTPUT(HWORK,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL DZERO(GWORK,N2ORBX*N2ORBX)
      DO 65 P=1,NORBT
      DO 65 Q=1,NORBT
      DO 65 R=1,NORBT
      DO 65 S=1,NORBT
      DO 65 A=1,NORBT
         GWORK(P,Q,R,S)=GWORK(P,Q,R,S)-0.5D0*(SXMAT(P,A,X)*UMOG(A,Q,R,S)
     +                                     - SXMAT(Q,A,X)*UMOG(P,A,R,S)
     +                                     + SXMAT(R,A,X)*UMOG(P,Q,A,S)
     +                                     - SXMAT(S,A,X)*UMOG(P,Q,R,A))
65    CONTINUE
C
      IF (IPRINT .GT. 15) THEN
         CALL HEADER('Two-electron reorthonormalisation term',-1)
         CALL OUTPUT(GWORK,1,N2ORBX,1,N2ORBX,N2ORBX,N2ORBX,1,LUPRI)
      END IF
C
C     *********************************************
C     ***** Add {S,h} and {S,g} to unmodified *****
C     *****        derivative integrals       *****
C     *********************************************
C
C     OMOHA <--- UMOHA + HWORK     (cf. equation 16a)
C     OMOGA <--- UMOGA + GWORK     (cf. equation 16b)
C
      DO 70 P=1,NORBT
      DO 70 Q=1,NORBT
         OMOHA(P,Q)=UMOHA(P,Q) + HWORK(P,Q)
      DO 70 R=1,NORBT
      DO 70 S=1,NORBT
         OMOGA(P,Q,R,S)=UMOGA(P,Q,R,S) + GWORK(P,Q,R,S)
70    CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Orthonormal one-electron derivative integrals',-1)
         CALL OUTPUT(OMOHA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      IF (IPRINT .GT. 15) THEN
         CALL HEADER('Orthonormal two-electron derivative integrals',-1)
         CALL OUTPUT(OMOGA,1,N2ORBX,1,N2ORBX,N2ORBX,N2ORBX,1,LUPRI)
      END IF
C
C     *****************************************************
C     ***** Contract density matrices with unmodified *****
C     *****    derivative integrals to form total     *****
C     *****       derivative of the Fock matrix       *****
C     *****************************************************
C
C     HWORK <--- DDENS * OMOHA + 2 PDENS * OMOGA   (cf. equation 30)
C
      DO 80 P=1,NORBT
      DO 80 Q=1,NORBT
      DO 80 B=1,NORBT
         DFOCKU(P,Q,X)=DFOCKU(P,Q,X)+DDENS(P,B)*UMOHA(Q,B)
         DFOCK(P,Q,X)=DFOCK(P,Q,X)+DDENS(P,B)*OMOHA(Q,B)
      DO 80 C=1,NORBT
      DO 80 D=1,NORBT
         DFOCKU(P,Q,X)=DFOCKU(P,Q,X)+2D0*PDENS(P,B,C,D)
     +                                  *UMOGA(Q,B,C,D)
         DFOCK(P,Q,X)=DFOCK(P,Q,X)+2D0*PDENS(P,B,C,D)
     +                                  *OMOGA(Q,B,C,D)
80    CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Derivative of the generalised Fock matrix',-1)
         CALL OUTPUT(DFOCK(1,1,X),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     *****************************************
C     ***** Construct orbital part of the *****
C     *****    differentiated gradient    *****
C     *****************************************
C
      DO 90 P=1,NORBT
      DO 90 Q=1,NORBT
90       OFOCK(P,Q,X)=2D0*(DFOCK(P,Q,X)-DFOCK(Q,P,X))
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Orbital part of the differentiated gradient',-1)
         CALL OUTPUT(OFOCK(1,1,X),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     *********************************************************
C     ***** Form unmodified part of the AO differentiated *****
C     *****      inactive and active Fock matrices        *****
C     *****            and auxiliary matrix               *****
C     *********************************************************
C
      DO 220 P=1,NORBT
      DO 220 Q=1,NORBT
         CFKUA(P,Q,X)=UMOHA(P,Q)
      DO 220 I=1,NISHT
220       CFKUA(P,Q,X)=CFKUA(P,Q,X)+2D0*UMOGA(P,Q,I,I)-UMOGA(P,I,I,Q)
C
      DO 230 P=1,NORBT
      DO 230 Q=1,NORBT
      DO 230 T=1,NASHT
      DO 230 U=1,NASHT
230       AFKUA(P,Q,X)=AFKUA(P,Q,X)+DDENS(NT(T),NT(U))
     +            *(UMOGA(P,Q,NT(T),NT(U))-0.5D0*UMOGA(P,NT(T),NT(U),Q))
C
      DO 240 T=1,NASHT
      DO 240 Q=1,NORBT
      DO 240 U=1,NASHT
      DO 240 V=1,NASHT
      DO 240 W=1,NASHT
240       QMTUA(NT(T),Q,X)=QMTUA(NT(T),Q,X)
     +    +2D0*PDENS(NT(T),NT(U),NT(V),NT(W))*UMOGA(Q,NT(U),NT(V),NT(W))
3     CONTINUE
C
C     ******************************
C     ***** Construct Y matrix *****
C     ******************************
C
      CALL DZERO(YMAT,N2ORBX*3)
      DO 100 X=1,3
      DO 110 P=1,NORBT
      DO 110 Q=1,NORBT
         YMAT(P,Q,X)=0.5D0*(DFOCK(P,Q,X) + DFOCKU(P,Q,X))
      DO 110 A=1,NORBT
         YMAT(P,Q,X)=YMAT(P,Q,X)+0.5D0*SXMAT(P,A,X)*FMAT(A,Q)
110   CONTINUE
C
      CALL HEADER('Y matrix for direction '//CHRXYZ(X),-1)
      CALL OUTPUT(YMAT(1,1,X),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
100   CONTINUE
C
C     *******************************************************
C     ***** Construct Y contribution to the static part *****
C     *****       of the molecular derivatives          *****
C     *******************************************************
C
C     Construct second part of Eab in equation (29b)
C
      CALL HEADER('Reorthonormalization terms',-1)
      DO 200 A=1,3
      DO 200 B=1,3
         EAB = 0.D0
      DO 210 P=1,NORBT
      DO 210 Q=1,NORBT
         EAB=EAB-SXMAT(P,Q,A)*YMAT(Q,P,B)-SXMAT(P,Q,B)*YMAT(Q,P,A)
210   CONTINUE
         WRITE (LUPRI,'(1X,2I5,F12.6)') A, B, EAB
200   CONTINUE
C
C     ###########################################
C     ##### Put in call to routine to check #####
C     #####    partitioned Fock matrices    #####
C     ###########################################
C
      CALL HEADER('Partitioned Fock matrix terms',-1)
      CALL PFKCHK(CMO,DDENS,PDENS,SXMAT,UMOH,UMOG,CFKUA,AFKUA,QMTUA,
     +                  WORK,LWORK,IPRINT)
      RETURN
      END
C  /* Deck ofabmo */
      SUBROUTINE OFABMO(CMO,DV,FMAT,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <inforb.h>
      DIMENSION CMO(NCMOT), DV(NASHT,NASHT), FMAT(N2ORBT), WORK(LWORK)
      LOGICAL FOUND
C
      CALL RD_SIRIFC('CMO',FOUND,CMO,WORK,LWORK)
      IF (.NOT.FOUND) CALL QUIT('OFABMO error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         IF (NNASHX .GT. LWORK) CALL STOPIT('OFABMO',' ',LWORK,NNASHX)
         CALL RD_SIRIFC('DV',FOUND,WORK,WORK,LWORK)
         IF (.NOT.FOUND)
     &      CALL QUIT('OFABMO error: DV not found on SIRIFC')
         CALL DSPTSI(NASHT,WORK,DV)
      END IF
      CALL RD_SIRIFC('FOCK',FOUND,FMAT,WORK,LWORK)
      IF (.NOT.FOUND)
     &   CALL QUIT('OFABMO error: FOCK not found on SIRIFC')
      RETURN
      END
C  /* Deck ofabhg */
      SUBROUTINE OFABHG(HWORK,UMOG,WORK,LWORK)
#include <implicit.h>
#include <iratdef.h>
      PARAMETER (LBUF = 600)
#include <inforb.h>
      DIMENSION WORK(LWORK)
      KIBUF = 1
      KWONE = KIBUF + 4*LBUF
      KTOP  = KWONE + NNBAST
      KLEFT = LWORK - KTOP + 1
      IF (KTOP .GT. LWORK) CALL STOPIT('OFABHG',' ',LWORK,KTOP)
      CALL OFABH1(HWORK,UMOG,WORK(KIBUF),WORK(KWONE),
     &            WORK(KTOP),KLEFT)
      RETURN
      END
C  /* Deck ofabh1 */
      SUBROUTINE OFABH1(HWORK,UMOG,IBUF,HWRK,WRK,LWORK)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
#include <inforb.h>
      PARAMETER (LBUF=600)
      INTEGER A, B, C, D
      DIMENSION HWORK(NBAST,NBAST), UMOG(NBAST,NBAST,NBAST,NBAST),
     &          HWRK(NNBAST), IBUF(4,LBUF), MBAS(8), WRK(LWORK)
#include <ibtfun.h>
C
      KFREE = 1
      LFREE = LWORK
C
C     One-electron part
C
      CALL DZERO(HWRK,NNBAST)
      CALL RDONEL('ONEHAMIL',.TRUE.,HWRK,NNBAST)
      CALL DSPTSI(NBAST,HWRK,HWORK)
C
C     Two-electron part
C
      LUINTA = -1
      CALL GPOPEN(LUINTA,'AOTWOINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) MSYM, (MBAS(I),I=1,8), KBUF, NIBUF, NBITS
      IF (KBUF .NE. LBUF) CALL QUIT('OFABH1: KBUF .ne. LBUF')
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('INTE',KINT ,LENINT,WRK,KFREE,LFREE)
      KIINT = KINT + LBUF
C
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
      CALL DZERO(UMOG,N2BASX*N2BASX)
      JBUF = 608
  200 CONTINUE
#if defined (VAR_SPLITFILES)         
         CALL READSI(LUINTA,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LUINTA,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,IBUF,NINT)
         IF (NINT .EQ. 0) GO TO 200
         IF (NINT .LT. 0) GO TO 400
         IF (NINT .GT. 0) THEN
            DO 300 I = 1, NINT
               GINT  = WRK(KINT + I - 1)
               A = IBUF(1,I)
               B = IBUF(2,I)
               C = IBUF(3,I)
               D = IBUF(4,I)
               UMOG(A,B,C,D) = GINT
               UMOG(A,B,D,C) = GINT
               UMOG(B,A,C,D) = GINT
               UMOG(B,A,D,C) = GINT
               UMOG(C,D,A,B) = GINT
               UMOG(C,D,B,A) = GINT
               UMOG(D,C,A,B) = GINT
               UMOG(D,C,B,A) = GINT
  300       CONTINUE
         END IF
      GO TO 200
  400 CONTINUE
      CALL GPCLOSE(LUINTA,'KEEP')
      RETURN
      END
C  /* Deck rumoga */
      SUBROUTINE RUMOGA(UMOGA,NORBT,NOH2,NOLOND,WORK,LWORK,X)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      INTEGER X
      LOGICAL NOH2, NOLOND
      DIMENSION UMOGA(NORBT,NORBT,NORBT,NORBT), WORK(LWORK)
C
      KIBUF  = 1
      KLAST = KIBUF + (600 + 1)/IRAT - 1
      IF (KLAST .GT. LWORK) CALL STOPIT('RUMOGA',' ',KLAST,LWORK)
      CALL RUMOG1(UMOGA,NORBT,NOH2,NOLOND,WORK(KIBUF),X,WORK(KLAST))
      RETURN
      END
C  /* Deck rumog1 */
      SUBROUTINE RUMOG1(UMOGA,NORBT,NOH2,NOLOND,IBUF,X,WORK)
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
      PARAMETER (LBUF = 600)
      LOGICAL NOH2, NOLOND
      INTEGER P,Q,R,S,X
      DIMENSION UMOGA(NORBT,NORBT,NORBT,NORBT), IBUF(600), WORK(*)
#include <inftap.h>      
#include <ibtfun.h>
C
      JBUF = 256
C
      LENINT = IRAT*LBUF + LBUF + 1
      KINT = 1
      KIINT = KINT + LBUF
C
      CALL DZERO(UMOGA,NORBT**4)
      IF (NOH2 .OR. NOLOND) RETURN
#if defined (VAR_INT64)
      call quit('RUMOGA1: NBITS=16 not implemented here yet')
!     Feb 2011 TODO : see her2out.F for how LU2DER is written
#endif
C
      CALL REWSPL(LU2DER)
      CALL MOLLAB('AO2MGINT',LU2DER,LUPRI)
 150  CONTINUE
#if defined (VAR_SPLITFILES)         
         CALL READSI(LU2DER,LENINT,WORK(KINT),JBUF)
#else
         CALL READI(LU2DER,LENINT,WORK(KINT))
#endif
         CALL AOLAB4(WORK(KIINT),LBUF,NIBUF,NBITS,IBUF,LENGTH)
         IF (LENGTH .GT. 0) THEN
            DO 100 I = 1, LENGTH
               LABEL = IBUF(I)
               VALUE = WORK(KINT + I - 1)
               S = IBTAND(LABEL,255)
               IF (S .EQ. 0) THEN
                  ICOMP = IBTAND(IBTSHR(LABEL,8),255)
               ELSE IF (ICOMP .EQ. X) THEN
                  P = IBTAND(IBTSHR(LABEL,24),255)
                  Q = IBTAND(IBTSHR(LABEL,16),255)
                  R = IBTAND(IBTSHR(LABEL, 8),255)
                  UMOGA(P,Q,R,S) = VALUE
                  UMOGA(R,S,P,Q) = VALUE
                  UMOGA(Q,P,S,R) = -VALUE
                  UMOGA(S,R,Q,P) = -VALUE
               END IF
 100        CONTINUE
         ELSE IF (LENGTH .LT. 0) THEN
            GOTO 300
         END IF
         GOTO 150
 300  CONTINUE
      RETURN
      END
C  /* Deck getdns */
      SUBROUTINE GETDNS(DDENS,PDENS)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.D0)
#include <inforb.h>
      DIMENSION DDENS(NORBT,NORBT), PDENS(NORBT,NORBT,NORBT,NORBT)
      CALL DZERO(DDENS,N2ORBX)
      CALL DZERO(PDENS,N2ORBX*N2ORBX)
      DO 100 I = 1, NOCCT
         DDENS(I,I) = D2
  100 CONTINUE
      DO 200 I = 1, NOCCT
         DO 300 J = 1, NOCCT
            PDENS(I,I,J,J) = D2
  300    CONTINUE
  200 CONTINUE
      DO 400 I = 1, NOCCT
         DO 500 J = 1, NOCCT
            PDENS(I,J,J,I) = PDENS(I,J,J,I) - D1
  500    CONTINUE
  400 CONTINUE
      RETURN
      END
C  /* Deck pfkchk */
      SUBROUTINE PFKCHK(CMO,DDENS,PDENS,SXMAT,UMOH,UMOG,CFKUA,AFKUA,
     +                  QMTUA,WORK,LWORK,IPRINT)
C
C     NODDY SUBROUTINE TO CHECK PARTITIONED FOCK MATRICES
C     FOR MAGNETIC RESPONSE EQUATIONS (Eq. 51b JCP 84, 6266 (1986))
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
      DIMENSION  WORK(LWORK),
     &           UMOG(NORBT,NORBT,NORBT,NORBT),
     &           UMOH(NORBT,NORBT),
     &           CFKUA(NORBT,NORBT,3),
     &           AFKUA(NORBT,NORBT,3),
     &           QMTUA(NORBT,NORBT,3),
     &           PDENS(NORBT,NORBT,NORBT,NORBT), DDENS(NORBT,NORBT),
     &           SXMAT(NORBT,NORBT,3),
     &           CMO  (NBAST,NORBT)
C
#include <inforb.h>
C
      IF (NSYM .GT. 1) THEN
         WRITE (LUPRI,'(1X,A)') ' PFKCHK only works without symmetry.'
         CALL QUIT('PFKCHK does not work with symmetry,')
      END IF
C
C
C     Allocate core memory
C
      ICFOK   = 1
      IAFOK   = ICFOK + N2ORBX
      IQMAT   = IAFOK + N2ORBX
      IFMAT   = IQMAT + N2ORBX
      ICFKSA = IFMAT + N2ORBX
      IAFKSA = ICFKSA + N2ORBX*3
      IQMTSA = IAFKSA + N2ORBX*3
      ICWMAT  = IQMTSA + N2ORBX*3
      IAWMAT  = ICWMAT + N2ORBX*3
      ICDFOK  = IAWMAT + N2ORBX*3
      IADFOK  = ICDFOK + N2ORBX*3
      IQDMAT  = IADFOK + N2ORBX*3
      IDFOCK  = IQDMAT + N2ORBX*3
      IWRK1   = IDFOCK + N2ORBX
      IWRK2 = IWRK1 + N2BASX
      KTOP  = IWRK2 + N2BASX
      IF (KTOP .GT. LWORK) CALL STOPIT('PFKCHK',' ',KTOP,LWORK)
C
      LWRK = LWORK - KTOP + 1
C
      CALL PFKBOD(CMO,DDENS,PDENS,SXMAT,UMOH,UMOG,CFKUA,AFKUA,QMTUA,
     +            WORK(ICFOK),WORK(IAFOK),WORK(IQMAT),WORK(IFMAT),
     +            WORK(ICFKSA),WORK(IAFKSA),WORK(IQMTSA),
     +            WORK(ICWMAT),WORK(IAWMAT),WORK(ICDFOK),WORK(IADFOK),
     +            WORK(IQDMAT),WORK(IDFOCK),WORK(IWRK1),WORK(IWRK2),
     +            WORK(KTOP),LWRK,IPRINT)
      RETURN
      END
C  /* Deck pfkbod */
      SUBROUTINE PFKBOD(CMO,DDENS,PDENS,SXMAT,UMOH,UMOG,CFKUA,AFKUA,
     +                  QMTUA,CFOK,AFOK,QMAT,FMAT,CFKSA,AFKSA,QMTSA,
     +                  CWMAT,AWMAT,CDFOK,ADFOK,QDMAT,DFOCK,WRK1,WRK2,
     +                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
      INTEGER P, Q, R, S, A, B, C, D, T, U, V, X
      DIMENSION  WORK(LWORK),
     &           UMOG(NORBT,NORBT,NORBT,NORBT),
     &           UMOH(NORBT,NORBT),
     &           PDENS(NORBT,NORBT,NORBT,NORBT), DDENS(NORBT,NORBT),
     &           SXMAT(NORBT,NORBT,3),
     &           CMO  (NBAST,NORBT),
     &           CFOK(NORBT,NORBT),
     &           AFOK(NORBT,NORBT),
     &           QMAT(NORBT,NORBT),
     &           FMAT(NORBT,NORBT),
     &           CFKUA(NORBT,NORBT,3),
     &           AFKUA(NORBT,NORBT,3),
     &           QMTUA(NORBT,NORBT,3),
     &           CFKSA(NORBT,NORBT,3),
     &           AFKSA(NORBT,NORBT,3),
     &           QMTSA(NORBT,NORBT,3),
     &           CWMAT(NORBT,NORBT,3),
     &           AWMAT(NORBT,NORBT,3),
     &           CDFOK(NORBT,NORBT,3),
     &           ADFOK(NORBT,NORBT,3),
     &           QDMAT(NORBT,NORBT,3),
     &           DFOCK(NORBT,NORBT,3),
     &           WRK1(NBAST,NBAST),
     &           WRK2(NBAST,NBAST)
#include <inforb.h>
#include <chrxyz.h>
      NT(T) = NISHT + T
C     ***************************************************************
C     ***** Noddy subroutine to calculate the AO differentiated *****
C     *****        inactive and active Fock matrices            *****
C     ***************************************************************
C
C     At end hope to be able to compare AO differentiated
C     inactive and active Fock matrices
C
C     **************************************************************
C     ***** Form inactive and active Fock matrices in MO basis *****
C     *****         Back transform if need check on AO         *****
C     **************************************************************
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('One-electron density matrix',-1)
         CALL OUTPUT(DDENS,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Two-electron density matrix',-1)
         CALL OUTPUT(PDENS,1,N2ORBX,1,N2ORBX,N2ORBX,N2ORBX,1,LUPRI)
      END IF
C
C     Form inactive Fock matrix
C
      CALL DZERO(CFOK,N2ORBX)
      DO 10 P=1,NORBT
      DO 10 Q=1,NORBT
         CFOK(P,Q)=UMOH(P,Q)
      DO 10 I=1,NOCCT
10       CFOK(P,Q)=CFOK(P,Q)+2D0*UMOG(P,Q,I,I)-UMOG(P,I,I,Q)
C
C     Form active Fock matrix
C
      CALL DZERO(AFOK,N2ORBX)
      DO 20 P=1,NORBT
      DO 20 Q=1,NORBT
      DO 20 T=1,NASHT
      DO 20 U=1,NASHT
20       AFOK(P,Q)=AFOK(P,Q)+DDENS(NT(T),NT(U))*(UMOG(P,Q,NT(T),NT(U))
     +                            -0.5D0*UMOG(P,NT(T),NT(U),Q))
C
C     Form auxiliary matrix
C
      CALL DZERO(QMAT,N2ORBX)
      DO 30 T=1,NASHT
      DO 30 Q=1,NORBT
      DO 30 U=1,NASHT
      DO 30 V=1,NASHT
      DO 30 X=1,NASHT
30       QMAT(NT(T),Q)=QMAT(NT(T),Q)+2D0*PDENS(NT(T),NT(U),NT(V),NT(X))
     &                                  *UMOG(Q,NT(U),NT(V),NT(X))
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Inactive Fock matrix in MO basis',-1)
         CALL OUTPUT(CFOK,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         CALL HEADER('Active Fock matrix in MO basis',-1)
         CALL OUTPUT(AFOK,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         CALL HEADER('Q matrix in MO basis',-1)
         CALL OUTPUT(QMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL DZERO(WRK1,N2BASX)
      CALL DZERO(WRK2,N2BASX)
      DO 40 M=1,NBAST
      DO 40 N=1,NBAST
      DO 40 P=1,NORBT
      DO 40 Q=1,NORBT
         WRK1(M,N)=WRK1(M,N)+CFOK(P,Q)*CMO(M,P)*CMO(N,Q)
40       WRK2(M,N)=WRK2(M,N)+AFOK(P,Q)*CMO(M,P)*CMO(N,Q)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Inactive Fock matrix in AO basis',-1)
         CALL OUTPUT(WRK1,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Active Fock matrix in AO basis',-1)
         CALL OUTPUT(WRK2,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     ****************************************************************
C     ***** Form the partitions of the Fock matrix and write out *****
C     *****      so have a check on the total Fock matrix        *****
C     ****************************************************************
C
      CALL DZERO(FMAT,N2ORBX)
      DO 50 I=1,NISHT
      DO 50 Q=1,NORBT
50       FMAT(I,Q)=2D0*(CFOK(I,Q)+AFOK(I,Q))
      DO 52 T=1,NASHT
      DO 52 Q=1,NORBT
         FMAT(NT(T),Q)=QMAT(NT(T),Q)
      DO 52 U=1,NASHT
52       FMAT(NT(T),Q)=FMAT(NT(T),Q)+DDENS(NT(T),NT(U))*CFOK(Q,NT(U))
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Total Fock matrix in MO basis',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      ENDIF
C
      CALL DZERO(CFKSA,3*N2ORBX)
      CALL DZERO(AFKSA,3*N2ORBX)
      CALL DZERO(QMTSA,3*N2ORBX)
      CALL DZERO(CWMAT,3*N2ORBX)
      CALL DZERO(AWMAT,3*N2ORBX)
      CALL DZERO(DFOCK,3*N2ORBX)
       DO 1 X=1,3
       WRITE(LUPRI,'(//A,I3,/A)') ' Cartesian direction ',X,
     *                            ' ======================= '
C
C     *********************************************************
C     ***** Form unmodified part of the AO differentiated *****
C     *****      inactive and active Fock matrices        *****
C     *****            and auxiliary matrix               *****
C     *********************************************************
C
      DO 75 M=1,NBAST
      DO 75 N=1,NBAST
      DO 75 P=1,NORBT
      DO 75 Q=1,NORBT
         WRK1(M,N)=WRK1(M,N)+CFKUA(P,Q,X)*CMO(M,P)*CMO(N,Q)
75       WRK2(M,N)=WRK2(M,N)+AFKUA(P,Q,X)*CMO(M,P)*CMO(N,Q)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Inactive u-derivative Fock matrix in AO basis',-1)
         CALL OUTPUT(WRK1,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Active u-derivative Fock matrix in AO basis',-1)
         CALL OUTPUT(WRK2,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     ****************************************************
C     ***** Form reorthonormalisation part of the AO *****
C     *****    differentiated inactive and active    *****
C     *****    Fock matrices and auxiliary matrix    *****
C     ****************************************************
C
      DO 90 P=1,NORBT
      DO 90 Q=1,NORBT
      DO 90 A=1,NORBT
         CFKSA(P,Q,X)=CFKSA(P,Q,X)-0.5D0*(SXMAT(P,A,X)*CFOK(A,Q)
     +                                   -SXMAT(Q,A,X)*CFOK(P,A))
         AFKSA(P,Q,X)=AFKSA(P,Q,X)-0.5D0*(SXMAT(P,A,X)*AFOK(A,Q)
     +                                   -SXMAT(Q,A,X)*AFOK(P,A))
90    CONTINUE
C
      DO 95 P=1,NORBT
      DO 95 Q=1,NORBT
      DO 95 T=1,NASHT
      DO 95 U=1,NASHT
      DO 95 V=1,NASHT
      DO 95 A=1,NORBT
         QMTSA(P,Q,X)=QMTSA(P,Q,X)-0.5D0*PDENS(P,NT(T),NT(U),NT(V))
     +                    *(SXMAT(Q,    A,X)*UMOG(A,NT(T),NT(U),NT(V))
     +                     -SXMAT(NT(T),A,X)*UMOG(Q,A    ,NT(U),NT(V))
     +                     +SXMAT(NT(U),A,X)*UMOG(Q,NT(T),A    ,NT(V))
     +                     -SXMAT(NT(V),A,X)*UMOG(Q,NT(T),NT(U),A))
95    CONTINUE
C
C     Add in the W matrices (cf. equation 42) to CFKSA and AFKSA
C
      DO 100 P=1,NORBT
      DO 100 Q=1,NORBT
      DO 100 I=1,NISHT
      DO 100 R=1,NORBT
         CWMAT(P,Q,X)=CWMAT(P,Q,X)-0.5D0*SXMAT(I,R,X)
     &                     *(UMOG(P,R,I,Q)-UMOG(P,I,R,Q))
100    CONTINUE
C
      DO 105 P=1,NORBT
      DO 105 Q=1,NORBT
      DO 105 R=1,NORBT
      DO 105 T=1,NORBT
      DO 105 U=1,NORBT
         AWMAT(P,Q,X)=AWMAT(P,Q,X)-0.5D0*DDENS(NT(T),NT(U))
     +                    *SXMAT(NT(T),R,X)
     +                    *0.5D0*(UMOG(P,R,NT(U),Q) - UMOG(P,NT(U),R,Q))
105    CONTINUE
C
      DO 110 P=1,NORBT
      DO 110 Q=1,NORBT
         CFKSA(P,Q,X)=CFKSA(P,Q,X)+CWMAT(P,Q,X)
         AFKSA(P,Q,X)=AFKSA(P,Q,X)+AWMAT(P,Q,X)
110    CONTINUE
C
C     **********************************************************
C     ***** Construct total inactive and active derivative *****
C     *****  Fock matrices and derivative auxiliary matrix *****
C     **********************************************************
C
      DO 115 P=1,NORBT
      DO 115 Q=1,NORBT
         CDFOK(P,Q,X)=CFKUA(P,Q,X)+CFKSA(P,Q,X)
         ADFOK(P,Q,X)=AFKUA(P,Q,X)+AFKSA(P,Q,X)
         QDMAT(P,Q,X)=QMTUA(P,Q,X)+QMTSA(P,Q,X)
115    CONTINUE
C
C     **********************************************************
C     ***** Finally construct total derivative Fock matrix *****
C     **********************************************************
C
      DO 120 I=1,NISHT
      DO 120 Q=1,NORBT
120       DFOCK(I,Q,X)=2D0*(CDFOK(I,Q,X)+ADFOK(I,Q,X))
      DO 122 Q=1,NORBT
      DO 122 T=1,NASHT
         DFOCK(NT(T),Q,X)=QDMAT(NT(T),Q,X)
      DO 122 U=1,NASHT
122       DFOCK(NT(T),Q,X)=DFOCK(NT(T),Q,X)+DDENS(NT(T),NT(U))
     +                     *CDFOK(Q,NT(U),X)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Total derivative Fock matrix in MO basis',-1)
         CALL OUTPUT(DFOCK(1,1,X),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      ENDIF
1     CONTINUE
      RETURN
      END
C  /* Deck psorhs */
      SUBROUTINE PSORHS(WMAT,LABINT,WORK,LWORK)
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
      CHARACTER*8 LABINT(*)
      DIMENSION WMAT(*), WORK(LWORK)
#include <nuclei.h>
#include <inforb.h>
#include <infdim.h>
#include <cbirhs.h>
#include <abainf.h>
#include <sigma.h>
#include <spnout.h>
C INFRSP: using LPVMAT
#include <infrsp.h>
C
      IF ((.NOT. (SHIELD .OR. SPNSPN .OR. SPINRO)) .OR. 
     &   ((.NOT. (SHIELD .OR. SPINRO))
     &    .AND. (SPNSPN .AND. (.NOT. DOPSO)))) RETURN
      CALL QENTER('PSORHS')
C
      IF (IPRALL .GE. 3) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPRALL .GE. 4) CALL TITLER('Output from PSORHS','*',103)
C
      IF (SHIELD .OR. SPINRO) CALL DZERO(SIGMAS,9*NUCDEP)
      IF (ABASOP) THEN
         LUDV   = NORBT * NORBT
         LPVX   = LPVMAT
      ELSE
         LUDV   = N2ASHX
         LPVX   = 0
      ENDIF
      KINTAD = 1
      KBINTR = KINTAD + (3*NUCDEP + 1)/IRAT
      KCMO   = KBINTR + (3*NUCDEP + 1)/IRAT
      KUDV   = KCMO   + NCMOT
      KPVX   = KUDV   + LUDV
      KINDX  = KPVX   + LPVX
      KCSTRA = KINDX  + LCINDX
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      LWRK   = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('PSORHS',' ',KLAST,LWORK)
      CALL PSORH1(NODV,WMAT,WORK(KINTAD),LABINT,
     &            WORK(KBINTR),WORK(KCMO),WORK(KUDV),WORK(KPVX),
     &            WORK(KINDX),WORK(KLAST),LWRK,NODDY,IPRALL,
     &            (SHIELD .OR. SPINRO),RALLCO,IPRALL)
C
      IF ((SHIELD .OR. SPINRO) .AND. IPRALL .GE. 4) THEN
         CALL HEADER('Reorthonormalization part of shielding',-1)
         CALL FCPRI(SIGMAS,'SIGMA',WORK(KCSTRA),WORK(KSCTRA))
      END IF
      IF (IPRALL .GE. 3) CALL TIMER('PSORHS',TIMSTR,TIMEND)
      CALL QEXIT('PSORHS')
      RETURN
      END
C  /* Deck psorh1 */
      SUBROUTINE PSORH1(NODV,WMAT,INTADR,LABINT,INTREP,CMO,UDV,PVX,
     &                  XINDX,WORK,LWORK,PSOTST,IPRALL,DOSHI,RALLCO,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <maxorb.h>
#include <maxmom.h>
#include <maxaqn.h>
#include <mxcent.h>
      PARAMETER (DM2 = -2.D0)
      INTEGER X
      LOGICAL NODV, PSOTST, DOSHI, RALLCO
      CHARACTER LABINT(*)*8, LABTMP(3*MXCOOR)*8
      DIMENSION WMAT(*), INTREP(*), INTADR(*), WORK(LWORK),
     &          CMO(*), UDV(*), PVX(*), XINDX(*)
#include <abainf.h>
#include <inftap.h>
#include <nuclei.h>
#include <inforb.h>
#include <infvar.h>
#include <infpri.h>
#include <inflin.h>
#include <infdim.h>
#include <infrsp.h>
#include <gdvec.h>
#include <wrkrsp.h>
#include <symmet.h>
#include <abares.h>
#include <infsop.h>
C
      IPRRSP = IPRALL - 2
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NCMOT,CMO)
      READ (LUSIFC)
      IF (NASHT .GT. 0) THEN
         IF (NODV) THEN
            READ (LUSIFC)
            CALL DZERO(UDV,N2ASHX)
         ELSE
            IF (LWORK.LT.NNASHX)CALL STOPIT('PSORH1','UDV',LWORK,NNASHX)
            CALL READI(LUSIFC,IRAT*NNASHX,WORK)
            CALL DSPTSI(NASHT,WORK,UDV)
         END IF
      ELSE
         READ (LUSIFC)
      END IF
C
      KSYMOP = LSYMRF
C
      TRPLET = .FALSE.
C
      IF (.NOT. ABASOP) THEN
         CALL GETCIX(XINDX,LSYMRF,LSYMRF,WORK,LWORK,0)
      ELSE
C
C       SOPPA :
C
C       Initialize XINDX
C
        CALL DZERO(XINDX,LCINDX)
C
C       Find address array's for SOPPA calculation
C
        CALL SET2SOPPA(XINDX(KABSAD),XINDX(KABTAD),
     *                 XINDX(KIJSAD),XINDX(KIJTAD),
     *                 XINDX(KIJ1AD),XINDX(KIJ2AD),
     *                 XINDX(KIJ3AD),XINDX(KIADR1))
C
        IF (CCPPA) THEN
          CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
        ELSE
          CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
        ENDIF
C
C       reads the MP2 or CCSD correlation coefficients into PV
C
        CALL READT (LUSIFC,LPVMAT,PVX)
C
        IF (IPRINT.GT.10) THEN
          IF (CCPPA) THEN
            WRITE(LUPRI,'(/A)')' PSORH1 : CCSD correlation ',
     &                         'coefficients'
          ELSE
            WRITE(LUPRI,'(/A,A)')' PSORH1 :',
     &                           ' MP2 correlation coefficients'
          ENDIF
          CALL OUTPUT(PVX,1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
        END IF
C
C       reads the MP2 or CCSD second order one particle density matrix 
C
        CALL READT (LUSIFC,NORBT*NORBT,UDV)
C
C       UDV contains the MP2 one-density. Remove the diagonal
C       contribution from the zeroth order. (Added in MP2FAC)
C
        IF (IPRINT.GT.10) THEN
          IF (CCPPA) THEN
            WRITE(LUPRI,'(/A)')' PSORH1 : CCSD density'
          ELSE
            WRITE(LUPRI,'(/A)')' PSORH1 : MP2 density'
          END IF
          CALL OUTPUT(UDV,1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
     &                  LUPRI)
        END IF
C
        CALL SOPUDV(UDV)
C
      END IF
C
C     Construct the appropriate property vectors
C
      IF (DODRCT .AND. .NOT. RALLCO) THEN
         NPATOM = 2
         ISTRT  = 1
         DO 77 I = 1, NUCIND
            NCOMP  = 0
            NPATOM = 2
            CALL GET1IN(DUMMY,'PSO    ',NCOMP,WORK,LWORK,LABTMP,
     &                  INTREP,INTADR,I,.TRUE.,NPATOM,.TRUE.,
     &                  DUMMY,.FALSE.,DUMMY,IPRINT)
            ITYP  = 0
            ITYP2 = 0
            DO 78 J = 1, NSYM
               DO 79 K = 1, NUCIND
                  DO 81 L = 1, 3
                     ISCOOR = IPTCNT(3*(K - 1) + L,J - 1, 2)
                     IF (ISCOOR .GT. 0) THEN
                        ITYP = ITYP + 1
                        IF (I .EQ. K) THEN
                           ITYP2 = ITYP2 + 1
                           LABINT(ITYP) = LABTMP(ITYP2)
                        END IF
                     END IF
 81               CONTINUE
 79            CONTINUE
 78         CONTINUE
 77      CONTINUE
      ELSE
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'PSO    ',NCOMP,WORK,LWORK,LABINT,
     &               INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,
     &               DUMMY,.FALSE.,DUMMY,IPRINT)
      END IF
      NREC = 3*NUCDEP
C
C     Construct right-hand sides and solve response equations and
C     calculate if nescessary relaxation contributions. One symmetry
C     at a time.
C
      DO 100 ISYM = 1, NSYM
         CALL ABAVAR(ISYM,.FALSE.,IPRINT,WORK,LWORK)
         IREFSY = LSYMRF
         NCREF  = NCONRF
         KSYMST = LSYMST
         KSYMOP = LSYMPT
         KZWOPT = NWOPT
         KZCONF = NCONST
         KZVAR  = KZWOPT + KZCONF
         KZYWOP = 2*KZWOPT
         KZYCON = 2*KZCONF
         KZYVAR = 2*KZVAR
         NREC = NGDVEC(ISYM,2)
         DO 200 LREC = 1, NREC
            IREC   = IGDREC(LREC,ISYM,2)
            ISCOOR = IGDCOR(LREC,ISYM,2)
            IF (ISCOOR .GT. 0) THEN
C
               IF (DOSHI) THEN
                  CALL PSOORT(LABINT(ISCOOR),ISCOOR,KSYMOP,WMAT,CMO,
     &                        WORK,LWORK,IPRALL)
               END IF
               IMGLAB(IREC) = LABINT(ISCOOR)
               CALL GETGPV(LABINT(ISCOOR),DUMMY,DUMMY,CMO,UDV,PVX,XINDX,
     &                     ANTSYM,WORK,LWORK)
C
C              Scaled by -2 in accordance with ABARHS
C
               CALL DSCAL(KZVAR,DM2,WORK,1)
               IF (IPRALL .GE. 05) THEN
                  WRITE (LUPRI,'(//A,I3)')
     &              ' Dif. Gradient for symmetry coordinate',ISCOOR
                  IF (IPRALL .GT. 20) THEN
                     WRITE(LUPRI,'(/,A,/)')' CSF part of gradient'
                     WRITE(LUPRI,'(6F12.6)') (WORK(J), J = 1, NCONST)
                  END IF
                  WRITE(LUPRI,'(//,A)')' Orbital part of gradient'
                  IF (IPRALL .GE. 20) THEN
                     PRFAC = 0.0D0
                  ELSE
                     PRFAC = 0.1D0
                  END IF
                  CALL PRKAP(NWOPPT,WORK(1+NCONST),PRFAC,LUPRI)
               END IF
               CALL WRITDX(LUGDI,ISCOOR,IRAT*NVARPT,WORK)
               IF (PSOTST) THEN
                  CALL PSONOD(ISCOOR,NODV,WORK,LWORK,IPRALL)
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck psoort */
      SUBROUTINE PSOORT(LABEL,ISCOOR,KSYMPSO,WMAT,CMO,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
      CHARACTER LABEL*8
      DIMENSION WORK(LWORK), CMO(*), WMAT(*)
#include <inforb.h>
C
      KPSO  = 1
      KLAST = KPSO + N2ORBX
      IF (KLAST .GT. LWORK) CALL STOPIT('PSOORT',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
      CALL PSOOR1(LABEL,ISCOOR,KSYMPSO,WMAT,WORK(KPSO),CMO,
     &            WORK(KLAST),LWRK,IPRINT)
      RETURN
C
      END
C  /* Deck psoor1 */
      SUBROUTINE PSOOR1(LABEL,ISCOOR,KSYMPSO,WMAT,PSOINT,CMO,
     &                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      INTEGER X
      CHARACTER LABEL*8
      DIMENSION PSOINT(NORBT,NORBT), WMAT(N2ORBX,3), CMO(*),
     &          WORK(LWORK)
#include <inforb.h>
#include <sigma.h>
#include <chrxyz.h>
C
      KSYMP = -1
      CALL PRPGET(LABEL,CMO,PSOINT,KSYMP,ANTSYM,WORK,LWORK,IPRINT)
      IF (KSYMP.NE.KSYMPSO)
     &   CALL QUIT('ERROR: unexpected property symmetry in PSOOR1')
      DO 100 X = 1, 3
         ORTO = DDOT(N2ORBX,PSOINT,1,WMAT(1,X),1)
         SIGMAS(X,ISCOOR) = ORTO
         IF (IPRINT .GE. 10) THEN
            CALL HEADER('W matrix in PSOOR1 for dir. '//CHRXYZ(X),-1)
            CALL OUTPUT(WMAT(1,X),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL HEADER('PSO integrals in PSOOR1',-1)
            CALL OUTPUT(PSOINT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE (LUPRI,'(A,F12.6,2I5)')
     &        ' Shielding reorthonormalization in PSOOR1:',ORTO,ISCOOR,X
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck psonod */
      SUBROUTINE PSONOD(ISCOOR,NODV,WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <dummy.h>
C
      LOGICAL NODV
      DIMENSION WORK(LWORK)
C
#include <nuclei.h>
#include <inforb.h>
C
      IF (NSYM .GT. 1) THEN
         WRITE (LUPRI,'(1X,A)') ' PSONOD only works without symmetry.'
         CALL QUIT('PSONOD does not work with symmetry,')
      END IF
C
C     Allocate core memory
C
      KPSOAO = 1
      KPSOMO = KPSOAO + 3*NUCDEP*N2BASX
      KPDENS = KPSOMO + N2ORBX
      KDDENS = KPDENS + N2ORBX*N2ORBX
      KDV    = KDDENS + N2ORBX
      KFOCK  = KDV    + N2ASHX
      KFMAT  = KFOCK  + N2ORBX
      KCMO   = KFMAT  + N2ORBX
      KGRAD  = KCMO   + N2ORBX
      KWRK   = KGRAD  + N2ORBX
      LWRK   = LWORK - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('PSONOD',' ',KWRK,LWORK)
C
      CALL PSONO1(ISCOOR,NODV,WORK(KPSOAO),WORK(KPSOMO),WORK(KPDENS),
     &            WORK(KDDENS),WORK(KDV),WORK(KFOCK),WORK(KFMAT),
     &            WORK(KCMO),WORK(KGRAD),WORK(KWRK),LWRK,IPRINT)
      RETURN
      END
C  /* Deck psono1 */
      SUBROUTINE PSONO1(ISCOOR,NODV,PSOAO,PSOMO,PDENS,DDENS,DV,FOCK,
     &                  FMAT,CMO,GRAD,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      LOGICAL NODV
      INTEGER P, Q, R, S, A, B, C, D
      DIMENSION  PSOAO(NBAST,NBAST,3*NUCDEP),
     &           PSOMO(NORBT,NORBT),
     &           DDENS(NORBT,NORBT),
     &           PDENS(NORBT,NORBT,NORBT,NORBT),
     &           DV(NASHT,NASHT),
     &           FMAT(NORBT,NORBT),
     &           FOCK(NORBT,NORBT),
     &           CMO  (NBAST,NORBT),
     &           GRAD(NORBT,NORBT), WORK(LWORK)
#include <abainf.h>
#include <inforb.h>
C
C     Read in MO's
C     ============
C
      CALL OFABMO(CMO,DV,FMAT,WORK,LWORK)
      IF (NODV) CALL DZERO(DV,N2ASHX)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('MO coefficients',-1)
         CALL OUTPUT(CMO,1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
         CALL HEADER('DV matrix',-1)
         CALL OUTPUT(DV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
C
C
C     Density matrices read into DDENS and PDENS
C
      CALL GETDNS(DDENS,PDENS)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('One-electron density matrix',-1)
         CALL OUTPUT(DDENS,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     PSO integrals in AO basis
C     =========================
C
      NCOMP = 3*NUCDEP
      CALL GET1PR(PSOAO,'PSO    ',NCOMP,'SQUARE',.TRUE.,WORK,LWORK,
     &            IDUMMY,5)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('PSO integrals in AO basis (PSONOD)',-1)
         CALL OUTPUT(PSOAO(1,1,ISCOOR),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &               LUPRI)
      END IF
C
C     PSO integrals in MO basis
C     =========================
C
      CALL DZERO(PSOMO,N2ORBX)
      DO 30 P=1,NORBT
      DO 30 Q=1,NORBT
      DO 30 A=1,NBAST
      DO 30 B=1,NBAST
         PSOMO(P,Q)=PSOMO(P,Q)+PSOAO(A,B,ISCOOR)*CMO(A,P)*CMO(B,Q)
30    CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('PSO integrals in MO basis (PSONOD)',-1)
         CALL OUTPUT(PSOMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Fock matrix
C     ===========
C
      CALL DZERO(FOCK,N2ORBX)
      DO 80 P=1,NORBT
      DO 80 Q=1,NORBT
      DO 80 B=1,NORBT
         FOCK(P,Q)=FOCK(P,Q)+DDENS(P,B)*PSOMO(Q,B)
80    CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Generalised PSO Fock matrix (PSONOD)',-1)
         CALL OUTPUT(FOCK,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Gradient
C     ========
C
      DO 90 P=1,NORBT
      DO 90 Q=1,NORBT
90       GRAD(P,Q)=2D0*(FOCK(P,Q)-FOCK(Q,P))
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Orbital part of PSO gradient (PSONOD)',-1)
         CALL OUTPUT(GRAD(1,1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck getw */
      SUBROUTINE GETW(WMAT,ANTI,CMO,TD,DV,ICOOR,WORK,LWORK,IPRINT)
C
C     tuh 210992 - apres Oui!
C
#include <implicit.h>
#include <priunit.h>
      LOGICAL ANTI
      DIMENSION TD(*), DV(*), WMAT(*), CMO(*), WORK(LWORK)
#include <inforb.h>
C
      KDTOT = 1
      KLAST = KDTOT + N2ORBX
      IF (KLAST .GT. LWORK) CALL STOPIT('GETW',' ',KLAST,LWORK)
      CALL GETW1(CMO,TD,DV,WORK(KDTOT),WMAT,ANTI,ICOOR,IPRINT)
      RETURN
      END
C  /* Deck getw1 */
      SUBROUTINE GETW1(CMO,TD,DV,DTOT,WMAT,ANTI,ICOOR,IPRINT)
C
C     tuh 210992 - apres Oui!
C
#include <implicit.h>
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
      PARAMETER(D2 = 2.D0)
C
      LOGICAL ANTI
      INTEGER U, U1, U2, V, V1, V2, UV
      DIMENSION TD(NORBT,NORBT), DV(NNASHX), WMAT(N2ORBX), CMO(*),
     &          DTOT(NORBT,NORBT)
#include <inforb.h>
#include <infind.h>
C
      IF (IPRINT .GE. 5) THEN
         CALL TITLER('Output from GETW1','*',103)
         WRITE (LUPRI,'(/A,I5)') ' ICOOR ', ICOOR
         IF (IPRINT .GE. 10) THEN
            CALL AROUND('TD matrix in GETW1')
            CALL OUTPUT(TD,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL AROUND('DV matrix in GETW1')
            CALL OUTPAK(DV,NASHT,1,LUPRI)
         END IF
      END IF
C
C     Construct total density matrix
C
      CALL DZERO(DTOT,N2ORBX)
C
      DO 100 ISYMU = 1, NSYM
         NASHU = NASH(ISYMU)
         NORBQ = NORB(ISYMU)
         DO 110 I = 1, NISH(ISYMU)
            INACT = IORB(ISYMU) + I
            DTOT(INACT,INACT) = D2
  110    CONTINUE
         IF (NASHU.GT.0 .AND. NORBQ.GT.0) THEN
            IOFFU1 = IORB(ISYMU) + NISH(ISYMU)
            IOFFU2 = IASH(ISYMU)
            DO 120 U = 1, NASHU
               U1 = IOFFU1 + U
               U2 = IOFFU2 + U
               DO 130 V = 1, NASHU
                  V1 = IOFFU1 + V
                  V2 = IOFFU2 + V
                  IF (U2 .GT. V2) THEN
                     UV = IROW(U2) + V2
                  ELSE
                     UV = IROW(V2) + U2
                  END IF
                  DTOT(U1,V1) = DTOT(U1,V1) + DV(UV)
 130           CONTINUE
 120        CONTINUE
         END IF
 100  CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL AROUND('DTOT in GETW1')
         CALL OUTPUT(DTOT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     WMAT matrix = DTOT*TD
C     But if we use differentiated creation operators we use:
C     WMAT = DTOT*TD(T) -/+ TD*DTOT, K.Ruud, Nov.-93,Feb.-95
C
      CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
     &           DTOT,NORBT,
     &           TD,NORBT,0.D0,
     &           WMAT,NORBT)
      IF (ANTI) THEN
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,-1.D0,
     &              TD,NORBT,
     &              DTOT,NORBT,1.D0,
     &              WMAT,NORBT)
      ELSE
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &              TD,NORBT,
     &              DTOT,NORBT,1.D0,
     &              WMAT,NORBT)
      END IF
      IF (IPRINT .GE. 10) THEN
         CALL AROUND('Final W matrix in GETW1')
         CALL OUTPUT(WMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck solrhs */
      SUBROUTINE SOLRHS(CMO,TD,DV,UMAT,WORK,LWORK,ICOOR,ISCOOR,
     &                  ISYM,IPRMAX)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      DIMENSION CMO(*), TD(*), DV(*), UMAT(*), WORK(LWORK)
#include <cbisol.h>
#include <inforb.h>
      KGLM    = 1
      KTLMD   = KGLM   + LMTOT
      KWMAT   = KTLMD  + LMTOT
      KSYRLM  = KWMAT  + N2ORBX
      KRLMMO  = KSYRLM + (LMNTOT + 1)/IRAT
      KRLMAO  = KRLMMO + N2ORBX
      KLAST   = KRLMAO + NNBASX
      IF (KLAST .GT. LWORK) CALL STOPIT('SOLRHS',' ',KLAST,LWORK)
      LWRK    = LWORK  - KLAST + 1
      CALL SOLRH1(CMO,TD,DV,WORK(KGLM),WORK(KTLMD),WORK(KWMAT),
     &            UMAT,WORK(KRLMAO),WORK(KRLMMO),WORK(KLAST),
     &            LWRK,ICOOR,ISCOOR,ISYM,WORK(KSYRLM),IPRMAX)
      RETURN
      END
C  /* Deck solrh1 */
      SUBROUTINE SOLRH1(CMO,TD,DV,GLM,TLMD,WMAT,UMAT,RLMAO,RLMMO,
     &                  WORK,LWORK,ICOOR,ISCOOR,ISYM,ISYRLM,IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <thrzer.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.D0, THRESH = 1.0D-8)
      LOGICAL ANTI
      DIMENSION GLM(LMTOT), TLMD(LMTOT), TD(NORBT,NORBT),
     &          DV(*), CMO(*), WORK(LWORK), WMAT(NORBT,NORBT),
     &          RLMAO(NNBASX), RLMMO(NORBT,NORBT),
     &          ISYRLM(LMNTOT), UMAT(NORBT,NORBT)
#include <cbisol.h>
#include <orgcom.h>
#include <inforb.h>
#include <infinp.h>
#include <inftap.h>
C
      DIMENSION CAVCNT(3)
C
C     Construct GLM
C
      CALL SOLFL(GLM,EPDIEL,RCAV,LCAVMX)
C
C     Construct WMAT
C
      ANTI = .FALSE.
      CALL GETW(WMAT,ANTI,CMO,TD,DV,ICOOR,WORK,LWORK,IPRINT)
C
C     Read in TLMD
C
      CALL READDX(LUTLM,ISCOOR,IRAT*LMTOT,TLMD)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('TLMD before adding reorthonorm. in SOLRH1',-1)
         CALL OUTPUT(TLMD,1,1,1,LMTOT,1,LMTOT,1,LUPRI)
         CALL HEADER('GLM constructed in SOLRH1',-1)
         CALL OUTPUT(GLM,1,1,1,LMTOT,1,LMTOT,1,LUPRI)
      END IF
C
C     Read LUSOL and check for consistency
C
      CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUPRI)
      READ (LUSOL) LMAXSS, LMTOTX, NNNBAS, CAVCNT
      IF (LMAXSS .LT. LSOLMX) THEN
         WRITE (LUPRI,'(//2A,2(/A,I5))') ' >>> SOLRH1 ERROR,',
     *   ' insufficient number of intgrals on LUSOL',
     *   ' l max from SIRIUS input :',LSOLMX,
     *   ' l max from LUSOL  file  :',LMAXSS
         CALL QUIT('SOLRH1: lmax on LUSOL is too small')
      END IF
      IF ((LMAXSS+1)**2 .NE. LMTOTX) THEN
         WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> SOLRH1 ERROR,',
     *   ' LUSOL file info inconsistent',
     *   ' l_max               :',LMAXSS,
     *   ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *   ' LMTOTX              :',LMTOTX
         CALL QUIT('SOLRH1: LUSOL info not internally consistent')
      END IF
      IF (NNNBAS .NE. NBAST) THEN
         WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> SOLRH1 ERROR,',
     *   ' LUSOL file info inconsistent with SIRIUS input',
     *   ' NBAST - LUSOL       :',NNNBAS,
     *   ' NBAST - SIRIUS      :',NBAST
         CALL QUIT('SOLRH1: LUSOL info inconsistent with SIRIUS input.')
      END IF
C
C     Some machines do not follow IEEE standard
C     and get in trouble with exact comparison 
C     /aug.98, Kenneth Ruud.
C
      CAVDIF = ABS(CAVCNT(1) - CAVORG(1)) +
     &         ABS(CAVCNT(2) - CAVORG(2)) +
     &         ABS(CAVCNT(3) - CAVORG(3)) 
      IF (CAVDIF .GT. THRZER) THEN
         WRITE (LUPRI,'(//2A,3(/A,3F20.15))') ' >>> SOLRH1 ERROR,',
     &   ' LUSOL center of cavity not consistent with ABACUS value',
     &   ' CAVORG(1:3) from LUSOL        :',CAVCNT,
     &   ' CAVORG(1:3) from common block :',CAVORG
         CALL QUIT('SOLRH1: LUSOL cavity center .ne. ABACUS value.')
      END IF
C
      CALL DZERO(UMAT,N2ORBX)
      READ (LUSOL)
C
      LM = 0
      DO 100 L = 0, LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
         DO 110 M = -L,L
            LM = LM + 1
            ISYMOP = ISYRLM(L + M + 1)
            IF (ABS(ISYMOP) .EQ. ISYM) THEN
               IF (ISYMOP .GT. 0) THEN
                  CALL READT(LUSOL,NNBASX,RLMAO)
               ELSE
                  CALL QUIT('SOLRH1 error:requested R(lm) not on LUSOL')
               END IF
               CALL TRAONE(CMO,RLMAO,RLMMO,WORK,LWORK,.FALSE.,ISYM,
     &                     IPRINT)
C
C              Add reorthonormalization term,
C              we get a minus from -TLMED in TLMD = TLMND - TLMED
C
               TLMD(LM) = TLMD(LM) - DDOT(N2ORBX,WMAT,1,RLMMO,1)
               FAC = -D2*GLM(LM)*TLMD(LM)
               CALL DAXPY(N2ORBX,FAC,RLMMO,1,UMAT,1)
               IF (IPRINT .GE. 15) THEN
                  CALL HEADER('RLMAO in SOLRH1',-1)
                  CALL OUTPAK(RLMAO,NBAST,1,LUPRI)
                  CALL HEADER('RLMMO in SOLRH1',-1)
                  CALL OUTPUT(RLMMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
            ELSE IF (ISYMOP .GT. 0) THEN
               READ (LUSOL)
            END IF
  110    CONTINUE
  100 CONTINUE
C
C     Write TLMD
C
      CALL WRITDX(LUTLM,ISCOOR,IRAT*LMTOT,TLMD)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('TLMD after adding reorthonorm. in SOLRH1',-1)
         CALL OUTPUT(TLMD,1,1,1,LMNTOT,1,LMNTOT,1,LUPRI)
         CALL HEADER('UMAT in SOLRH1',-1)
         CALL OUTPUT(UMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      CALL GPCLOSE(LUSOL,'KEEP')
      RETURN
      END
C  /* Deck traone */
      SUBROUTINE TRAONE(CMO,TRSO,SQMO,WORK,LWORK,ANTI,ISYMOP,IPRINT)
#include <implicit.h>
#include <priunit.h>
      LOGICAL ANTI
      DIMENSION CMO(*), TRSO(*), SQMO(*), WORK(LWORK)
#include <inforb.h>
      KSQSO = 1
      KLAST = KSQSO + N2BASX
      IF (KLAST .GT. LWORK) CALL STOPIT('TRAONE',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST + 1
      CALL TRAON1(CMO,SQMO,TRSO,WORK(KSQSO),WORK(KLAST),LWRK,ISYMOP,
     &            ANTI,IPRINT)
      RETURN
      END
C  /* Deck traon1 */
      SUBROUTINE TRAON1(CMO,SQMO,TRSO,SQSO,WORK,LWORK,ISYMOP,ANTI,
     &                  IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <symmet.h>
#include <inforb.h>
#include <infdim.h>
      DIMENSION CMO(NCMOT), WORK(LWORK),
     *          SQMO(NORBT,NORBT), TRSO(NNBASX), SQSO(N2BASX)
      LOGICAL   NONSYM, ANTI
      PARAMETER (D1 =1.0D0)
C
C     ***** Print Section *****
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from TRAON1',-1)
         WRITE (LUPRI,'(A,L5)') ' ANTI   ', ANTI
         WRITE (LUPRI,'(A,I5)') ' ISYMOP ',ISYMOP
         WRITE (LUPRI,'(A,I5)') ' NORBT  ', NORBT
         WRITE (LUPRI,'(A,I5)') ' NBAST  ', NBAST
         WRITE (LUPRI,'(A,2I5)') ' NNBASX, N2BASX ', NNBASX, N2BASX
         WRITE (LUPRI,'(A,2I5)') ' NNORBX, N2ORBX ', NNORBX, N2ORBX
         WRITE (LUPRI,'(A,I5)')  ' NCMOT  ', NCMOT
         WRITE (LUPRI,'(A,I10)') ' LWORK  ', LWORK
      END IF
C
C     ***** Make square matrix *****
C
      IF (.NOT.ANTI) THEN
         CALL DSPTSI(NBAST,TRSO,SQSO)
      ELSE
         CALL DAPTGE(NBAST,TRSO,SQSO)
      END IF
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Triangular SO matrix',-1)
         CALL OUTPAK(TRSO,NBAST,1,LUPRI)
         CALL HEADER('Square SO matrix',-1)
         CALL OUTPUT(SQSO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     ***** Transform matrix from SO to MO basis *****
C
      CALL DZERO(SQMO,N2ORBX)
C
C     Loop over irreps for orbitals
C
      NONSYM = .FALSE.
      DO 100 IREPA = 1,NSYM
      IF (NORB(IREPA) .GT. 0) THEN
         DO 200 IREPB = 1, IREPA
         IF (NORB(IREPB) .GT. 0) THEN
C
C           Check if operator belongs to this symmetry
C
            IREPO = MULD2H(IREPA,IREPB)
            IF (IREPO .EQ. ISYMOP) THEN
               NONSYM = NONSYM .OR. IREPO .GT. 1
C
C              Print symmetries and orbitals
C
               IF (IPRINT.GT.15) THEN
                  WRITE (LUPRI,'(A,3I3)') ' IREPA/B/O',IREPA,IREPB,IREPO
                  WRITE (LUPRI,'(/A,I5,A)')
     *               ' MO coefficients for symmetry', IREPA
                  CALL OUTPUT(CMO(ICMO(IREPA)+1),1,NBAS(IREPA),1,
     *               NORB(IREPA),NBAS(IREPA),NORB(IREPA),1,LUPRI)
                  WRITE (LUPRI,'(/A,I5,A)')
     *               ' MO coefficients for symmetry', IREPB
                  CALL OUTPUT(CMO(ICMO(IREPB)+1),1,NBAS(IREPB),1,
     *               NORB(IREPB),NBAS(IREPB),NORB(IREPB),1,LUPRI)
               END IF
C
C              Transform matrix block(s)
C
               CALL UTHV(CMO(ICMO(IREPA)+1),SQSO,CMO(ICMO(IREPB)+1),
     *                   IREPA,IREPB,NBAS(IREPA),NBAS(IREPB),SQMO,WORK)
C
C              Print transformed matrix thus far
C
               IF (IPRINT.GT.50) THEN
                  WRITE(LUPRI,'(/4A)')' Unfinished matrix in MO basis'
                  CALL OUTPUT(SQMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
            END IF
         END IF
 200     CONTINUE
      END IF
 100  CONTINUE
      IF (NONSYM) THEN
                   FACTOR =   D1
         IF (ANTI) FACTOR = - D1
         CALL TRANSX(SQMO,SQMO,NORBT,NORBT,FACTOR,IPRINT)
      ENDIF
      IF (IPRINT.GT.10) THEN
         CALL HEADER('Matrix in MO basis TRAON1',-1)
         CALL OUTPUT(SQMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck reort2 */
      SUBROUTINE REORT2(TD1,TD2,TEMP,FOCK,HESFS1,IPRINT,LUPRI)
C
C     Extra reorthonormalization term when doing new right-hand sides
C     K.Ruud, june-94
C
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION TD1(N2ORBX), TD2(N2ORBX), TEMP(N2ORBX),
     &          FOCK(N2ORBT), HESFS1(MXCOOR,MXCOOR)
#include <symmet.h>
#include <nuclei.h>
#include <inforb.h>
#include <dorps.h>
#include <inftap.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from REORT2',-1)
         IF (IPRINT .GE. 10) THEN
            CALL AROUND('Fock matrix in REORT2')
            DO 5 ISYM = 1, NSYM
               WRITE (LUPRI,'(A,I5)') ' Irrep', ISYM
               NORBI = NORB(ISYM)
               ISTRI = I2ORB(ISYM) + 1
               CALL OUTPUT(FOCK(ISTRI),1,NORBI,1,NORBI,
     *                     NORBI,NORBI,1,LUPRI)
 5          CONTINUE
         END IF
      END IF
C
      INOFF= 0
      SYELMT = D0
      DO 10 ISYM = 1, NSYM
         INDIR = NCRREP(ISYM - 1,1)
         IF (DOREPS(ISYM - 1) .AND. (INDIR .GT. 0)) THEN
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,2I5)')
     &           ' Symmetry and number of coordinates:',
     &           ISYM, INDIR
            DO 20 ICOOR = 1, INDIR
               IRCD = INOFF + ICOOR
               IF (DOPERT(IRCD,1)) THEN
                  IREC = 2*IRCD - 1
                  CALL READDX(LUSFDA,IREC,IRAT*N2ORBX,TD1)
                  JNOFF = 0
                  DO 30 JSYM = 1, NSYM
                     JNDIR = NCRREP(JSYM - 1,1)
                     IF (DOREPS(JSYM - 1) .AND. (JNDIR .GT. 0)) THEN
                        IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,2I5)')
     &                       ' Symmetry and number of coordinates:',
     &                       JSYM, JNDIR
                        DO 40 JCOOR = 1, JNDIR
                           JRCD = JNOFF + JCOOR
                           IF (DOPERT(JRCD,1)) THEN
                              IREC = 2*JRCD - 1
                              CALL READDX(LUSFDA,IREC,IRAT*N2ORBX,TD2)
                              CALL DGEMM('T','N',NORBT,NORBT,NORBT,1.D0,
     &                                   TD1,NORBT,
     &                                   TD2,NORBT,0.D0,
     &                                   TEMP,NORBT)
                              ISYMPT = MULD2H(ISYM,JSYM)
                              DO 50 ISYMP = 1, NSYM
                                 ISYMQ = MULD2H(ISYMP,ISYMPT)
                                 NORBQ = NORB(ISYMQ)
                                 NOCCQ = NOCC(ISYMQ)
                                 IF (NOCCQ .GT. 0) THEN
                                    IOFFQ1 = IORB(ISYMQ)
                                    IOFFRQ = I2ORB(ISYMQ)
                                    DO 60 IQ = 1, NOCCQ
                                       IQ1 = IOFFQ1 + IQ
                                       TERM = D0
                                       DO 70 IR = 1, NOCCQ
                                          IR1 = (IQ1 - 1)*NORBT
     &                                       + IOFFQ1 + IR
                                          IRQ = IOFFRQ+(IR - 1)*NORBQ+IQ
                                       TERM = TERM + TEMP(IR1)*FOCK(IRQ)
 70                                 CONTINUE
                                    SYELMT = SYELMT + TERM
 60                              CONTINUE
                              END IF
 50                        CONTINUE
                           IF (IPRINT .GT. 10) THEN
                              WRITE(LUPRI,'(1X,A,2I5,F10.5)')
     &                             'Additional correction term',
     &                             IRCD, JRCD, SYELMT
                           END IF
                           HESFS1(IRCD,JRCD) = HESFS1(IRCD,JRCD)
     &                                       + SYELMT
                           SYELMT = D0
                        END IF
 40                  CONTINUE
                  END IF
                  JNOFF = JNOFF + JNDIR
 30               CONTINUE
               END IF
 20         CONTINUE
         END IF
         INOFF = INOFF + INDIR
 10   CONTINUE
      IF (IPRINT .GT. 15) THEN
         CALL AROUND('Extra contribution to HESFS1')
         CALL OUTPUT(HESFS1,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,1,
     &        LUPRI)
      END IF
      RETURN
      END
C  /* Deck angrhs */
      SUBROUTINE ANGRHS(FDC,GP,XINDX,WRK,KFRSAV,LFRSAV)
C
C Dec 1998 Stephan P. A. Sauer
C
C Purpose: Calculate Abacus GP vectors in SOPPA of ANGMOM for right hand
C          side (RHS) of response equations.
C
C Input:
C   FDC    : derivative Fock matrix (property integrals) in MO basis
C
C Output:
C   GP     = right hand side for ABACUS response calculation
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
      DIMENSION FDC(N2ORBX)
      DIMENSION GP(2*NVARPT)
      DIMENSION XINDX(1),    WRK(LFRSAV)
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C
C   ABAINF : IPRALL
C   INFORB : NNASHX,N2ASHX,N2ORBX,NASHT,NORBT
C   INFLIN : NVARPT,LSYMRF,LSYMPT,NCONRF,NCONST,NWOPPT
C   INFRSP : LPVMAT,IPRRSP
C   INFSOP : KABSAD,KABTAD,KIJSAD,KIJTAD,KIJ1AD,KIJ2AD,KIJ3AD,KIADR1
C
#include <abainf.h>
#include <cbirhs.h>
#include <inftap.h>
#include <inforb.h>
#include <infdim.h>
#include <inflin.h>
#include <infrsp.h>
#include <infsop.h>
#include <wrkrsp.h>
C
C
      CALL QENTER('ANGRHS')
C
      IF (IPRALL .GE. 3) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPRALL .GE. 4) CALL TITLER('Output from ANGRHS','*',103)
C
C----------------------------------------------
C     Initialize variables for RESPONSE routines
C----------------------------------------------
C
      IPRRSP = IPRALL - 2
      TRPLET = .FALSE.
      IREFSY = LSYMRF
      NCREF  = NCONRF
      KSYMST = LSYMST
      KSYMOP = LSYMPT
      KZWOPT = NWOPPT
      KZCONF = NCONST
      KZVAR  = KZWOPT + KZCONF
      KZYWOP = 2*KZWOPT
      KZYCON = 2*KZCONF
      KZYVAR = 2*KZVAR
C
C------------------------------------------------------------
C     allocate memory for MP2 density in UDV 
C                 and for MP2 correlation coefficients in PVX
C------------------------------------------------------------
C
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      CALL MEMGET('REAL',KUDV,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVX,LPVMAT,WRK,KFREE,LFREE)
C
C-----------------------------------------------
C     Find address array's for SOPPA calculation
C     and save them in XINDX
C-----------------------------------------------
C
      CALL DZERO(XINDX,LCINDX)
C
      CALL SET2SOPPA(XINDX(KABSAD),XINDX(KABTAD),
     &               XINDX(KIJSAD),XINDX(KIJTAD),
     &               XINDX(KIJ1AD),XINDX(KIJ2AD),
     &               XINDX(KIJ3AD),XINDX(KIADR1))
C
C-----------------------------------------------------------
C     Read MP2/CCSD correlation coefficients in PVX
C     and the unrelaxed density matrix in UDV.
C     Remove the diagonal contribution to the density matrix 
C     from zeroth order, which was added in MP2FAC
C-----------------------------------------------------------
C
      REWIND LUSIFC
      IF (CCPPA) THEN
         CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
      ELSE
         CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
      ENDIF
C
      CALL READT (LUSIFC,LPVMAT,WRK(KPVX))
C
      IF (IPRRSP.GT.10) THEN
         IF (CCPPA) THEN
            WRITE(LUPRI,'(/A,A)')' ANGRHS :',
     &                           ' CCSD correlation coefficients'
         ELSE
            WRITE(LUPRI,'(/A,A)')' ANGRHS :',
     &                           ' MP2 correlation coefficients'
         ENDIF
         CALL OUTPUT(WRK(KPVX),1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
      END IF
C
      CALL READT (LUSIFC,NORBT*NORBT,WRK(KUDV))
C
      IF (IPRRSP.GT.10) THEN
         IF (CCPPA) THEN
            WRITE(LUPRI,'(/A)')' ANGRHS : CCSD density'
         ELSE
            WRITE(LUPRI,'(/A)')' ANGRHS : MP2 density'
         END IF
         CALL OUTPUT(WRK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,LUPRI)
      END IF
C
      CALL SOPUDV(WRK(KUDV))
C
C----------------------------------------------------------------------
C     Allocate memory for the second order ph contribution to GP vector
C----------------------------------------------------------------------
C
      CALL MEMGET('REAL',KPRPSE,N2ORBX,WRK,KFREE,LFREE)
C
C-------------------------
C     Initialize GP vector
C-------------------------
C
      CALL DZERO(GP,KZYVAR)
C
C----------------------------------------     
C     Construct orbital part of GP vector
C----------------------------------------
C
      IF (KZWOPT.GT.0) THEN
C
C        ... make second order contribution in WRK(KPRPSE)
C
         CALL PRPOMP(FDC,WRK(KPRPSE),WRK(KUDV))
C
C        ... add zeroth order contribution
C
         CALL PRPORS(FDC,WRK(KPRPSE),GP)
C
      ENDIF
C
C-------------------------------------
C     Construct 2p2h part of GP vector
C-------------------------------------
C
      IF (KZCONF.GT.0) THEN
C
         CALL PRPCMP(FDC,WRK(KPVX),GP,XINDX(KABSAD),
     &               XINDX(KABTAD),XINDX(KIJSAD),XINDX(KIJTAD),
     &               XINDX(KIADR1),WRK,KFREE,LFREE)
C
      ENDIF
C
C------------------------------------------
C     Scaled by 2 in accordance with ABARHS
C------------------------------------------
C
      CALL DSCAL(KZYVAR,D2,GP,1)

C
C--------------------------------
C     Return the allocated memory
C--------------------------------
C
      CALL MEMREL('ANGRHS',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('ANGRHS')
C
      IF (IPRALL .GE. 3) CALL TIMER('ANGRHS',TIMSTR,TIMEND)
C
      RETURN
      END
! -- end of abarhs.F --
