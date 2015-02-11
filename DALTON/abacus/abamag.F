!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
C  /* Deck maginp */
      SUBROUTINE MAGINP(WORD)
C
C     K.Ruud, April 1992
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      LOGICAL NEWDEF
      PARAMETER (NTABLE = 14)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "abainf.h"
#include "magone.h"
      DATA TABLE /'.SHIELD', '.SPIN-S', '.DIASUS', '.SKIP  ',
     &            '.PRINT ', '.NODC  ', '.NODV  ', '.STOP  ',
     &            '.POINTS', '.NEFIEL', '.QUADRU', '.ALL CO',
     &            '.ELFGRA', '.SECMOM'/
      NEWDEF = WORD .EQ. '*EXPECT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
         READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
            GOTO 100
         ELSE IF (PROMPT .EQ. '.') THEN
            ICHANG = ICHANG + 1
            DO 200 I = 1, NTABLE
               IF (TABLE(I) .EQ. WORD) THEN
                  GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14), I
               END IF
 200        CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GOTO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)') 'Keyword "', WORD,
     &           '" not recognized in EXPECT.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in EXPECT')

 1          CONTINUE !         .SHIELD
               NST    = .TRUE.
               GOTO 100
 2          CONTINUE !         .SPIN-S
               SST    = .TRUE.
               GOTO 100
 3          CONTINUE !         .DIASUS
               MDIASU = .TRUE.
               GOTO 100
 4          CONTINUE !         .SKIP
               MSKIP  = .TRUE.
               GOTO 100
 5          CONTINUE !         .PRINT
               READ (LUCMD,*) MPRINT
               IF (MPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GOTO 100
 6          CONTINUE !         .NODC
               MNODC  = .TRUE.
               GOTO 100
 7          CONTINUE
               MNODV  = .TRUE.
               GOTO 100
 8          CONTINUE
               MCUT   = .TRUE.
               GOTO 100
 9          CONTINUE
               READ (LUCMD,*) NPOINT
               IF (NPOINT .EQ. MPQUAD) THEN
                  ICHANG = ICHANG - 1
               ELSE
                  MPQUAD = NPOINT
                  MPOINT = .TRUE.
               END IF
               GOTO 100
 10         CONTINUE
               MELFLD = .TRUE.
               GOTO 100
 11         CONTINUE
               MQUADR = .TRUE.
               GOTO 100
 12         CONTINUE
               ALLCOM = .TRUE.
               IF (.NOT. DODRCT) ICHANG = ICHANG - 1
               GOTO 100
 13         CONTINUE
               MNQCC  = .TRUE.
               ALLCOM = .TRUE.
               GOTO 100
 14         CONTINUE
               MSECND = .TRUE.
               GOTO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GOTO 300
         ELSE
            WRITE (LUPRI,'(/,3A,/)') 'Prompt "', WORD,
     &           '" not recognized in EXPECT.'
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords', LUPRI)
         CALL QUIT('Illegal prompt in EXPECT')
      END IF
 300  CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for EXPECT:',0)
         IF (MSKIP) THEN
            WRITE (LUPRI,'(A)')
     &          ' No magnetic expectation values calculated in this run'
         ELSE
            IF (MPRINT .NE. IPRDEF) WRITE(LUPRI,'(A,I5)')
     &           ' Print level in EXPECT:',MPRINT
            IF (NST .OR. SST .OR. MDIASU .OR. MASSVE .OR. DARWIN .OR.
     &          MELFLD .OR. MQUADR) THEN
               WRITE(LUPRI,'(/A/)')
     &         ' The following expectation values will be calculated '//
     &          'in this run:'
               IF (NST) WRITE (LUPRI,'(10X,A)')
     &              'Nuclear shielding tensor'
               IF (MDIASU) WRITE (LUPRI,'(10X,A)')
     &                 'Diamagnetic susceptiblity tensor'
               IF (SST) THEN
                  WRITE (LUPRI,'(10X,A)')
     &                 'Spin-spin coupling tensor'
                  WRITE (LUPRI,'(10X,A,I5)')
     &                 'Number of integration points for DSO:', MPQUAD
               END IF
               IF (MELFLD) WRITE (LUPRI,'(10X,A)')
     &             'Nuclear electric fields (diamagnetic spin-rotation)'
               IF (MSECND) WRITE (LUPRI,'(10X,A)')
     &              'Second moments'
               IF (MQUADR) WRITE (LUPRI,'(10X,A)')
     &              'Quadrupole moments'
               IF (MNQCC) WRITE (LUPRI,'(10X,A)')
     &              'Electric field gradients'
            ELSE
               WRITE (LUPRI,'(/,A,/,A,/)')
     &              'No properties demanded in input.',
     &          'No one-electron expectation values calculated in '//
     &          'this run.'
C              MSKIP = .TRUE.
            END IF
            IF (MNODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     &           'density matrix neglected in EXPECT.'
            IF (MNODV) WRITE (LUPRI,'(/,2A)') 'Active one-electron',
     &           'density matrix neglected in EXPECT.'
            IF (ALLCOM .AND. DODRCT .AND. NST) WRITE(LUPRI,'(/,A)')
     &           ' All components of shielding expexctation value '//
     &           'calculated at once'
            IF (ALLCOM .AND. DODRCT .AND. SST) WRITE(LUPRI,'(/,A)')
     &           ' All components of diamagnetic spin-orbit operators'//
     &           ' calculated at once'
            IF (MCUT) WRITE (LUPRI,'(/,A)')
     &           ' Program is stopped after EXPECT.'
         END IF
      END IF
      RETURN
      END
C  /* Deck magini */
      SUBROUTINE MAGINI
C
C     Initialize MAGONE
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "abainf.h"
#include "magone.h"
COBL/AMT CAM magnetizabilities fix
#include "suscpt.h"
C

      MPRINT = IPRDEF
      MSKIP  = .NOT.(SHIELD .OR. MAGSUS .OR. SPNSPN .OR. VCD .OR.
     &               MASSVE .OR. DARWIN .OR. SPINRO .OR. MOLGFA)
      MNODC  = .FALSE.
      MNODV  = .FALSE.
      NST    = (SHIELD .OR. (SPINRO .AND. .NOT. NOLOND))
      SST    =  SPNSPN
      MDIASU = (MAGSUS .OR. (MOLGFA .AND. .NOT. NOLOND))
      TOFILE = .FALSE.
      MCUT   = .FALSE.
      MPOINT = .FALSE.
      MELFLD = (SPINRO .AND. NOLOND)
      MQUADR = QUADRU
      MSECND = SECNDM
      MTHIRD = THIRDM
      MNQCC  = NQCC
      ALLCOM = .NOT. DODRCT
      MPQUAD = 40
COBL/AMT CAM magnetizabilities fix
      DONETWOLOP = .FALSE.
      RETURN
      END
C  /* Deck magint */
      SUBROUTINE MAGINT(SPNDSO,WORK,LWORK,PASS)
C
C     K. Ruud, March 1992
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxmom.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.D0)
      LOGICAL PASS
      CHARACTER*7 WORD
      DIMENSION SPNDSO(*), WORK(LWORK)
#include "magone.h"
#include "nuclei.h"
#include "abainf.h"
#include "inforb.h"
#include "symmet.h"
#include "dorps.h"
#include "spnout.h"
#include "relcor.h"
#include "dipole.h"
#include "pcmlog.h"
      IF (MSKIP) RETURN
      IF (.NOT.DODSO .AND. (SST .AND. .NOT.(NST .OR. MDIASU))) RETURN
      CALL QENTER('MAGINT')
      IF (MPRINT .GT. 1) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (MPRINT .GT. 2) CALL TITLER('Output from MAGINT','*',103)
      DARWN  = D0
      RMASSV = D0
      CALL DZERO(DIPME,3)
C
C     **************************************************************
C     ***** Set up total density and Fock matrices in AO basis *****
C     **************************************************************
C
      KDENMA = 1
      KFOCMA = KDENMA + NNBASX
      KDSO   = KFOCMA + NNBASX
      KFSO   = KDSO   + NNBASX
      KLAST  = KFSO   + NNBASX
      LWRK   = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('MAGINT',' ',KLAST,LWORK)
C 
      CALL DZERO(WORK(KDSO),NNBASX)
      CALL DZERO(WORK(KFSO),NNBASX)
      CALL DSOFSO(WORK(KDSO),WORK(KFSO),WORK(KLAST),LWRK,MPRINT,
     &            MNODC,MNODV)
      CALL DSYM1(WORK(KDENMA),WORK(KFOCMA),WORK(KDSO),WORK(KFSO),
     &           .TRUE.,NBAST,MPRINT)
      KLAST = KDSO
      LWRK = LWRK + 2*NNBASX
C
C     We always start by calculating the dipole moment
C
      NCOMP = 3
      CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &            MPRINT,MNODC,MNODV,'DIPLEN ',NCOMP)
C
C     Local-field corrected dipole moment
C
      IF(PCM.AND.LOCFLD) THEN
         NCOMP = 3
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &        MPRINT,MNODC,MNODV,'DIPLOC ',NCOMP)
      END IF
C
      IF (NST) THEN
         IF (NOLOND) THEN
            WORD = 'NSTCGO '
         ELSE
            WORD = 'NUCSHI '
         END IF
         IF (ALLCOM .OR. .NOT. DODRCT) THEN
            NCOMP = 9*NUCDEP
         ELSE
            NDEGMX = 1
            DO 11 I = 1, NUCIND
               NDEG = NUCDEG(I)
               IF (NDEG .GT. NDEGMX) NDEGMX = NDEG
 11         CONTINUE
            NCOMP = 9*NDEGMX
         END IF
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
         IF (SPINRO .AND. .NOT. NOLOND) THEN
            WORD = 'NSTCGOS'
            NDEGMX = 1
            DO I = 1, NUCIND
               NDEG = NUCDEG(I)
               IF (NDEG .GT. NDEGMX) NDEGMX = NDEG
            END DO
            NCOMP = 9*NDEGMX
            CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),
     &                  LWRK,MPRINT,MNODC,MNODV,WORD,NCOMP)
         END IF
      END IF
      IF (SST .AND. DODSO) THEN
         WORD = 'DSO    '
C
C        Too many components at present
C
         IF (ALLCOM .OR. .NOT. DODRCT) THEN
            NCOMP = (3*NUCDEP*(3*NUCDEP + 1)/2)
         ELSE
            NDEGMX = 1
            DO I = 1, NUCIND
               NDEG = NUCDEG(I)
               IF (NDEG .GT. NDEGMX) NDEGMX = NDEG
            END DO
C
C     In current worst case scenario we need 9(1,2)+6(1,1)+6(2,2) components
C     for the spin-spin coupling constants of a pair of nuclei
C
            NCOMP = 21*NDEGMX*NDEGMX
         END IF
c        Numerical dso integration is faster for large molecules
c        but currently, it is implemented only without symmetry.
c        In case of no symmetry present, both methods should give
c        identical results.
         IF(NSYM.GT.1) THEN
            WRITE(LUPRI,'(/A)')
     &         "Symmetry -> DSO by analytical integration."
            CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),
     &           WORK(KLAST),LWRK,MPRINT,MNODC,MNODV,WORD,NCOMP)
         ELSE
            WRITE(LUPRI,'(/A)')
     &         "No symmetry -> DSO by numerical integration."
            CALL DZERO(SPNDSO,MXCOOR*MXCOOR)
            CALL NUMDSO(SPNDSO,NUCIND,WORK(KLAST),LWRK,MPRINT)
         ENDIF
      END IF
      IF (MDIASU) THEN
         IF (NOLOND) THEN
            WORD = 'DSUSCGO'
         ELSE
            WORD = 'DIASUS '
         END IF
         NCOMP = 6
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
         IF (MOLGFA .OR. (MAGSUS .AND. .NOT. NOLOND)) THEN
            WORD = 'DSUSCGM'
            NCOMP = 6
            CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),
     &                  LWRK,MPRINT,MNODC,MNODV,WORD,NCOMP)
         END IF
      END IF
      IF (DARWIN) THEN
         WORD = 'DARWIN '
         NCOMP = 1
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
      END IF
      IF (MASSVE) THEN
         WORD = 'MASSVEL'
         NCOMP = 1
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
      END IF
      IF (MELFLD) THEN
         WORD = 'NEFIELD'
         NCOMP = 3*NUCDEP
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
      END IF
      IF (MQUADR) THEN
         WORD = 'THETA  '
         NCOMP = 6
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
      END IF
      IF (MSECND) THEN
         WORD = 'SECMOM '
         NCOMP = 6
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
      END IF
      IF (MTHIRD) THEN
         WORD = 'THRMOM '
         NCOMP = 10
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
      END IF
      IF (MNQCC) THEN
         WORD = 'ELFGRDC'
         IF (ALLCOM .OR. .NOT. DODRCT) THEN
            NCOMP = 6*NUCDEP
         ELSE
            NDEGMX = 1
            DO 12 I = 1, NUCIND
               NDEG = NUCDEG(I)
               IF (NDEG .GT. NDEGMX) NDEXMG = NDEG
 12         CONTINUE
            NCOMP = 6*NDEGMX
         END IF
         CALL MAGDRV(SPNDSO,WORK(KDENMA),WORK(KFOCMA),WORK(KLAST),LWRK,
     &               MPRINT,MNODC,MNODV,WORD,NCOMP)
      END IF
      IF (MPRINT .GT. 1) CALL TIMER('MAGINT',TIMSTR,TIMEND)
      PASS = .TRUE.
      IF (MCUT) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Program stopped after MAGINT as requested.'
         CALL QUIT(' ***** End of ABACUS (in MAGINT) *****')
      END IF
      CALL QEXIT('MAGINT')
      RETURN
      END
C  /* Deck magdrv */
      SUBROUTINE MAGDRV(SPNDSO,DENMAT,FOCMAT,WORK,LWORK,IPRINT,NODC,
     &                  NODV,WORD,NCOMP)
C
C     K. Ruud, March 1992
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "maxmom.h"
#include "iratdef.h"
      LOGICAL NODC, NODV
      CHARACTER*7 WORD
      DIMENSION SPNDSO(*), WORK(LWORK), DENMAT(*), FOCMAT(*)
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "inforb.h"
C
       
C
      CALL QENTER('MAGDRV')
C
      MAXTYP = 3*MXCOOR
      IF (WORD .EQ. 'DSO    ') MAXTYP = (3*NUCDEP)**2
      KEXPVA = 1
      KINTRP = KEXPVA + NCOMP
      KINTAD = KINTRP + (9*MXCENT**2+1)/IRAT
      KFMOVL = KINTAD + (9*MXCENT**2+1)/IRAT
      KLAST  = KFMOVL + 6
      LWRK   = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('MAGDRV',' ',KLAST,LWORK)
      CALL MAGDR1(SPNDSO,DENMAT,FOCMAT,
     &            WORK(KEXPVA),WORK(KLAST),LWRK,IPRINT,NODC,NODV,
     &            NNBASX,NNBAST,NBAST,WORD,WORK(KINTRP),WORK(KINTAD),
     &            WORK(KFMOVL),NCOMP)
      CALL QEXIT('MAGDRV')
      RETURN
      END
C  /* Deck magdr1 */
      SUBROUTINE MAGDR1(SPNDSO,DENMAT,FOCMAT,EXPVAL,WORK,LWORK,
     &                  IPRINT,NODC,NODV,NNBASX,NNBAST,NBAST,WORD,
     &                  INTREP,INTADR,FMOVL,NCOMP)
C
C     K.Ruud, April 1992
C
C     This routine controls the calculation of the different magnetic
C     properties integrals. It also constructs the density matrix and
C     calculates the corresponding expectation values.
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "maxmom.h"
#include "iratdef.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      INTEGER TMAT
      LOGICAL NODC, NODV, EXP1VL
      CHARACTER WORD*7, LABINT(3*MXCOOR)*8
      DIMENSION SPNDSO(*),
     &          DENMAT(NNBASX), INTREP(*), INTADR(*),
     &          FOCMAT(NNBASX), EXPVAL(NCOMP), FMOVL(6),WORK(LWORK)
#include "gnrinf.h"
#include "abainf.h"
#include "efield.h"
#include "cbisol.h"
#include "magone.h"
#include "infinp.h"
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
#include "qmmm.h"
cDM added to pass vectors to MAGPCMSEC
#include "pcmdef.h"
#include "pcm.h"
#include "pcmlog.h"
cDM
C
      SIGN = 1.0D0
      EXP1VL = .TRUE.
      KCOMP = 0
      IF (IPRINT .GT. 2) CALL TITLER('Output form MAGDR1','*',103)
      IF (WORD .EQ. 'DSO    ') THEN
         NPATOM = 1
         CALL GET1IN(DUMMY,WORD,KCOMP,WORK(1),LWORK,LABINT,INTREP,
     &               INTADR,MPQUAD,.FALSE.,NPATOM,.TRUE.,EXPVAL,
     &               EXP1VL,DENMAT,IPRINT)
         CALL MAGAVE(SPNDSO,EXPVAL,KCOMP,WORK,LWORK,
     &               IPRINT,WORD,LABINT,SIGN)
      ELSE IF (WORD .EQ. 'NSTCGOS') THEN
         NPATOM = 2
         DO I = 1, NUCIND
            KCOMP = 0
            CALL GET1IN(DUMMY,WORD,KCOMP,WORK(1),LWORK,LABINT,INTREP,
     &                  INTADR,I,.FALSE.,NPATOM,.TRUE.,EXPVAL,
     &                  EXP1VL,DENMAT,IPRINT)
            CALL MAGAVE(SPNDSO,EXPVAL,KCOMP,WORK,LWORK,
     &                  IPRINT,WORD,LABINT,SIGN)
         END DO
      ELSE
         NPATOM = 0
         IF (WORD .EQ. 'DSUSCGM') THEN
            CALL GET1IN(DUMMY,'DSUSCGO',KCOMP,WORK(1),LWORK,LABINT,
     &           INTREP,INTADR,MPQUAD,.FALSE.,NPATOM,.TRUE.,EXPVAL,
     &           EXP1VL,DENMAT,IPRINT)
         ELSE
            CALL GET1IN(DUMMY,WORD,KCOMP,WORK(1),LWORK,LABINT,INTREP,
     &                  INTADR,MPQUAD,.FALSE.,NPATOM,.TRUE.,EXPVAL,
     &                  EXP1VL,DENMAT,IPRINT)
         END IF
C
C     Add contribution from electric field, KRu, Aug.-93
C
         IF (NFIELD .GT. 0 .AND. (WORD .EQ. 'DIASUS ')) THEN
            NCOMP  = 6
            KEXPT  = 1
            NPATOM = 0
            LWRK  = KEXPT + NCOMP
            IF (LWRK .GT. LWORK) CALL STOPIT('MAGDR1',' ',LWRK,LWORK)
            LLEFT = LWORK - LWRK + 1
            DO 555 IFIELD = 1, NFIELD
               IF (EFIELD(IFIELD) .NE. D0) THEN
                  IF (LFIELD(IFIELD) .EQ. 'XDIPLEN ') THEN
                     FIELD2 = 'X-FIELD'
                  ELSE IF (LFIELD(IFIELD) .EQ. 'YDIPLEN ') THEN
                     FIELD2 = 'Y-FIELD'
                  ELSE IF (LFIELD(IFIELD) .EQ. 'ZDIPLEN ') THEN
                     FIELD2 = 'Z-FIELD'
                  ELSE
                    WRITE (LUPRI,'(/,3A,/)') 'Field type ',
     &                  LFIELD(IFIELD), ' not defined for '//
     &                                  'magnetizabileties'
                     CALL QUIT('Illegal field type for '//
     &                         'magnetizabileties')
                  END IF
                  KCOMP = 0
                  CALL GET1IN(DUMMY,'CM2    ',KCOMP,WORK(LWRK),LLEFT,
     &                        LABINT,INTREP,INTADR,MPQUAD,.FALSE.,
     &                        NPATOM,.TRUE.,WORK(KEXPT),EXP1VL,DENMAT,
     &                        IPRINT)
                  CALL DAXPY(KCOMP,EFIELD(IFIELD),WORK(KEXPT),1,
     &                       EXPVAL,1)
               END IF
 555        CONTINUE
         END IF
C
         IF (QM3 .AND. (WORD .EQ. 'DIASUS')) THEN
            NCOMP  = 6
            KSLEXP = 1
            KLAST  = KSLEXP + NCOMP
            LWRK   = LWORK - KLAST + 1
            CALL QM3B2(WORK(KSLEXP),DENMAT,WORK(KLAST),LWRK)
            IF (IPRINT .GT. 5) THEN
               CALL AROUND('Second-order QM/MM contribution')
               CALL OUTPAK(WORK(KSLEXP),3,1,LUPRI)
            END IF
            CALL DAXPY(NCOMP,D1,WORK(KSLEXP),1,EXPVAL,1)
         END IF
C
         IF (QMMM .AND. (WORD .EQ. 'DIASUS')) THEN
            IF (IPQMMM .GT. 1) CALL TIMER('START ',TIMSTR,TIMEND)
            NCOMP  = 6
            KSLEXP = 1
            KLAST  = KSLEXP + NCOMP
            LWRK   = LWORK - KLAST + 1
            CALL QMMMB2(WORK(KSLEXP),DENMAT,WORK(KLAST),LWRK,IPRINT)
            IF (IPRINT .GT. 5) THEN
               CALL AROUND('Second-order QM/MM contribution')
               CALL OUTPAK(WORK(KSLEXP),3,1,LUPRI)
            END IF
            CALL DAXPY(NCOMP,D1,WORK(KSLEXP),1,EXPVAL,1)
            IF (IPQMMM .GT. 1) CALL TIMER('QMMMB2',TIMSTR,TIMEND)
         END IF
C
         IF ((SOLVNT.OR.PCM) .AND. (WORD .EQ. 'DIASUS')) THEN
            NCOMP  = 6
            IF (SOLVNT) THEN
cDM Onsager solvent model contributions
               KSLEXP = 1
               KSLSPH = KSLEXP + NCOMP
               KSLCAR = KSLSPH + (2*LCAVMX + 1)*NCOMP
               KGLM   = KSLCAR + (LCAVMX + 1)*(LCAVMX + 2)*NCOMP/2
               KTLM   = KGLM   + LMTOT
               KTRAMT = KTLM   + LMTOT
               KLAST  =KTRAMT+(LCAVMX + 1)*(LCAVMX + 2)*(2*LCAVMX + 1)/2
               IF (KLAST .GT. LWORK) CALL STOPIT('MAGDR1','MAGSOLEXP',
     &                                           KLAST,LWORK)
               LWRK   = LWORK - KLAST + 1
               CALL MAGSOLEXP(WORK(KSLEXP),WORK(KSLSPH),WORK(KSLCAR),
     &                        INTREP,INTADR,LABINT,WORK(KGLM),
     &                        WORK(KTLM),WORK(KTRAMT),IPRINT,
     &                        DENMAT,WORK(KLAST),LWRK)
            ELSE 
cDM IEF-PCM solvent model contributions
               CALL QUIT('PCM magnetizabilities not rewritten for new'//
     &                   'integral structure')
               KDERMA = 1
               KPOTPR = KDERMA +  NNBASX*NCOMP
               KQSEPR = KPOTPR + NTSIRR*NCOMP
               LWRK = KQSEPR + NTSIRR*NCOMP
               LLEFT=LWORK - LWRK + 1

               CALL MAGPCMSEC(WORK(KDERMA),NCOMP,IPRINT,IPRINT,
     &              WORK(LWRK),LLEFT,DENMAT,WORK(KPOTPR),WORK(KQSEPR))
            END IF
C
            IF (IPRINT.GT.5) THEN
               CALL AROUND('Second-order solvent contribution')
               CALL OUTPAK(WORK(KSLEXP),3,1,LUPRI)
            END IF

cDM add the solvent contribution to perturbed Fock matrix 
	    CALL DAXPY(NCOMP,D1,WORK(KSLEXP),1,EXPVAL,1)
C
         END IF
C
         CALL MAGAVE(SPNDSO,EXPVAL,KCOMP,WORK,LWORK,IPRINT,WORD,
     &               LABINT,SIGN)
      END IF
C
C     Multiplies Fock matrix with second derivative of overlap integrals
C
      IF (WORD .EQ. 'DIASUS ' .OR. WORD .EQ. 'DSUSCGO') THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'S2MAG  ',NCOMP,WORK(1),LWORK,LABINT,INTREP,
     &               INTADR,MPQUAD,.FALSE.,NPATOM,.TRUE.,EXPVAL,EXP1VL,
     &               FOCMAT,IPRINT)
         SIGN = -1.0D0
         CALL MAGAVE(SPNDSO,EXPVAL,NCOMP,WORK,LWORK,IPRINT,
     &               'S2MAG  ',LABINT,SIGN)
      END IF
      RETURN
      END
C  /* Deck magave */
      SUBROUTINE MAGAVE(SPNDSO,EXPVAL,NCOMP,WORK,
     &                  LWORK,IPRINT,WORD,LABINT,SIGN)
C
C     Distribute calculated expectation values into proper arrays
C
#include "implicit.h"
#include "priunit.h"
#include "maxmom.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "inforb.h"
      PARAMETER (DP5 = 0.5D0)
      CHARACTER WORD*7, LABINT(9*MXCENT)*8
      DIMENSION DENMAT(NNBASX),
     &          EXPVAL(NCOMP), SPNDSO(*), WORK(LWORK)
#include "symmet.h"
#include "nuclei.h"
#include "suscpt.h"
#include "quadru.h"
#include "secmom.h"
#include "thimom.h"
#include "molgfa.h"
#include "relcor.h"
#include "dipole.h"
#include "chrxyz.h"
C
      CALL QENTER('MAGAVE')
C
      IF (WORD .EQ. 'DARWIN ') DARWN  = -EXPVAL(1)
      IF (WORD .EQ. 'MASSVEL') RMASSV = -EXPVAL(1)
      IF (WORD .EQ. 'DIPLEN ') THEN
         DIPME(1) = EXPVAL(1)
         DIPME(2) = EXPVAL(2)
         DIPME(3) = EXPVAL(3)
Clf
         KCSTRA = 1
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         CALL DIPNUC(WORK(KCSTRA),WORK(KSCTRA),IPRINT,.FALSE.)
      END IF
Clf Calculation of dipole moment corrected with the local field factor
C for the electronic part we use the same structure as for diplen but
C using the DIPLOC integrals instead
C for the nuclear part we add the local field contribution calculated in 
C DIPLFN to the nuclear part of the dipole which is already present in dipmn
      IF (WORD .EQ. 'DIPLOC ') THEN
         DLFE(1) = EXPVAL(1)
         DLFE(2) = EXPVAL(2)
         DLFE(3) = EXPVAL(3)
         CALL DIPLFN(IPRINT)
         DLFN(1) = DLFN(1) + DIPMN(1)
         DLFN(2) = DLFN(2) + DIPMN(2)
         DLFN(3) = DLFN(3) + DIPMN(3)
      END IF
      IF (WORD .EQ. 'DIASUS ' .OR. WORD .EQ. 'DSUSCGO') THEN
         IJ = 0
         DO 100 I = 1, 3
         DO 100 J = I, 3
            IJ = IJ + 1
            SUSDIA(IPTAX(J,2),IPTAX(I,2)) = -EXPVAL(IJ)
            SUSDIA(IPTAX(I,2),IPTAX(J,2)) = -EXPVAL(IJ)
            GFACDI(IPTAX(J,2),IPTAX(I,2)) = -EXPVAL(IJ)
            GFACDI(IPTAX(I,2),IPTAX(J,2)) = -EXPVAL(IJ)
  100   CONTINUE
      END IF
      IF (WORD .EQ. 'DSUSCGM') THEN
         IJ = 0
         DO 101 I = 1, 3
         DO 101 J = I, 3
            IJ = IJ + 1
            SUSCOM(IPTAX(J,2),IPTAX(I,2)) = -EXPVAL(IJ)
            SUSCOM(IPTAX(I,2),IPTAX(J,2)) = -EXPVAL(IJ)
            GFACDI(IPTAX(J,2),IPTAX(I,2)) =
     &           GFACDI(IPTAX(J,2),IPTAX(I,2)) + EXPVAL(IJ)
            IF (I .NE. J) THEN
               GFACDI(IPTAX(I,2),IPTAX(J,2)) =
     &              GFACDI(IPTAX(I,2),IPTAX(J,2)) + EXPVAL(IJ)
            END IF
 101     CONTINUE
      END IF
      IF (WORD .EQ. 'THETA  ') THEN
         IJ = 0
         DO 110 I = 1, 3
         DO 110 J = I, 3
            IJ = IJ + 1
            QDREL(IPTAX(J,1),IPTAX(I,1)) = -EXPVAL(IJ)
            QDREL(IPTAX(I,1),IPTAX(J,1)) = -EXPVAL(IJ)
 110     CONTINUE
      END IF
      IF (WORD .EQ. 'SECMOM ') THEN
        IJ = 0
        DO 150 I = 1, 3
        DO 150 J = I, 3
           IJ = IJ + 1
           SCMEL(IPTAX(J,1),IPTAX(I,1)) = -EXPVAL(IJ)
           SCMEL(IPTAX(I,1),IPTAX(J,1)) = -EXPVAL(IJ)
 150    CONTINUE   
      END IF
      IF (WORD .EQ. 'THRMOM')THEN
         IJK = 0
         DO 160 I = 1, 3
         DO 160 J = I, 3
         DO 160 K = J, 3          
            IJK = IJK + 1
            THDMEL(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1)) = -EXPVAL(IJK)
            THDMEL(IPTAX(I,1),IPTAX(K,1),IPTAX(J,1)) = -EXPVAL(IJK)
            THDMEL(IPTAX(J,1),IPTAX(I,1),IPTAX(K,1)) = -EXPVAL(IJK)
            THDMEL(IPTAX(J,1),IPTAX(K,1),IPTAX(I,1)) = -EXPVAL(IJK)
            THDMEL(IPTAX(K,1),IPTAX(I,1),IPTAX(J,1)) = -EXPVAL(IJK)
            THDMEL(IPTAX(K,1),IPTAX(J,1),IPTAX(I,1)) = -EXPVAL(IJK)
 160     CONTINUE
      END IF
      IF (WORD .EQ. 'S2MAG  ') THEN
         IJ = 0
         DO 200 I = 1, 3
         DO 200 J = I, 3
            IJ = IJ + 1
            SUSFS(IPTAX(J,2),IPTAX(I,2)) = EXPVAL(IJ)
            SUSFS(IPTAX(I,2),IPTAX(J,2)) = EXPVAL(IJ)
  200   CONTINUE
      END IF
      IF (WORD .EQ. 'DSO    ') THEN
         KCSTRA =  1
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         LWRK   = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL STOPIT('MAGAVE',' ',KLAST,LWORK)
         CALL SPNTRA(SPNDSO,EXPVAL,WORK(KCSTRA),WORK(KSCTRA),NCOMP,
     &               LABINT,IPRINT)
      END IF
      IF (WORD .EQ. 'NUCSHI ' .OR. WORD .EQ. 'NSTCGO ' .OR.
     &    WORD .EQ. 'NSTCGOS') THEN
         NROW = 3*NUCDEP
         KIVEC  = 1
         KCSTRA = KIVEC  + 3*NROW
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         LWRK   = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL STOPIT('MAGAVE',' ',KLAST,LWORK)
         CALL SHITRA(WORK(KIVEC),EXPVAL,WORK(KCSTRA),WORK(KSCTRA),NCOMP,
     &               NROW,LABINT,WORD,IPRINT)
      END IF
      IF (WORD .EQ. 'NEFIELD') THEN
         KIVEC  = 1
         KCSTRA = KIVEC + NCOMP
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         LWRK   = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL STOPIT('MAGAVE',' ',KLAST,LWORK)
         CALL EFNTRA(WORK(KIVEC),EXPVAL,WORK(KCSTRA),WORK(KSCTRA),
     &               NCOMP,LABINT,IPRINT)
      END IF
      IF (WORD .EQ. 'ELFGRDC') THEN
         KIVEC = 1
         KLAST = KIVEC + 9*NUCDEP
         LWRK  = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL STOPIT('MAGAVE',' ',KLAST,LWORK)
         CALL NQCEL(WORK(KIVEC),EXPVAL,NCOMP,IPRINT)
      END IF
C
C     Print section
C
      IF (IPRINT .GT. 3) THEN
         WRITE (LUPRI,'()')
         IF (SIGN .GT. 0.D0) THEN
            WRITE (LUPRI,'(2A)') ' Expectation values for ', WORD
         ELSE
            WRITE (LUPRI,'(/2A)')' Trace of Fock matrix multiplied by ',
     &                           WORD
         END IF
         DO 30 ICOMP = 1, NCOMP
            WRITE (LUPRI,'(A,A,10X,F16.6)') ' Component: ',
     &            LABINT(ICOMP), EXPVAL(ICOMP)
 30      CONTINUE
         WRITE(LUPRI,'(//)')
         IF (WORD .EQ. 'NUCSHI ' .OR. WORD .EQ. 'NSTCGO ') THEN
            IOFF = 0
            DO 40 IATOM = 1, NUCDEP
               WRITE (LUPRI,1000) NAMDEP(IATOM),'x',
     &               (WORK(KIVEC + IOFF +     I*NROW), I = 0, 2)
               WRITE (LUPRI,1000) NAMDEP(IATOM),'y',
     &               (WORK(KIVEC + IOFF + 1 + I*NROW), I = 0, 2)
               WRITE (LUPRI,1000) NAMDEP(IATOM),'z',
     &               (WORK(KIVEC + IOFF + 2 + I*NROW), I = 0, 2)
               IOFF = IOFF + 3
               WRITE (LUPRI,'()')
 40         CONTINUE
         END IF
         IF (WORD .EQ. 'DIASUS ' .OR. WORD .EQ. 'S2MAG  ' .OR.
     &       WORD .EQ. 'DSUSCGO') THEN
            WRITE (LUPRI,'(24X,A,19X,A,19X,A/)') (CHRXYZ(I), I = 1,3)
            WRITE (LUPRI,'(10X,A,4X, F20.10)') CHRXYZ(1), EXPVAL(1)
            WRITE (LUPRI,'(10X,A,4X,2F20.10)') CHRXYZ(2), EXPVAL(2),
     &             EXPVAL(4)
            WRITE (LUPRI,'(10X,A,4X,3F20.10)') CHRXYZ(3), EXPVAL(3),
     &             EXPVAL(5), EXPVAL(6)
         END IF
         IF (WORD .EQ. 'DIPLEN ') THEN
            CALL HEADER('Electronic contributions to dipole moment',-1)
            CALL DP0PRI(DIPME)
         END IF
      END IF
      CALL QEXIT('MAGAVE')
      RETURN
 1000 FORMAT (1X,A6,1X,A1,' :',F17.10,2F24.10)
      END
C  /* Deck get1pr */
      SUBROUTINE GET1PR(SINTMA,WORD,NCOMP,MTFORM,TRIMAT,WORK,LWORK,
     &                  IORDER,IPRINT)
C
C     K.Ruud, March 1992. This routine generates the nescessary input
C     for getting first-order magnetic derivatives from HERMIT
C
#include "implicit.h"
#include "iratdef.h"
#include "dummy.h"
#include "priunit.h"
#include "inftap.h"
#include "mxcent.h"
#include "maxmom.h"
#include "cbiher.h"
#include "nuclei.h"
C
      LOGICAL TOFILE, TRIMAT, DOINT(2,2)
      CHARACTER WORD*7, MTFORM*6, LABINT(3*MXCOOR)*8
      DIMENSION WORK(LWORK), SINTMA(*)
C
      MAXTYP = 3*MXCOOR
      IF (WORD .EQ. 'DSO    ') MAXTYP = (3*NUCDEP)**2
      TOFILE    = .FALSE.
      TRIANG    = TRIMAT
      PROPRI    = IPRINT .GT. 5
C
      DOINT(1,1) = .TRUE.
      DOINT(2,2) = .TRUE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
C
      ALLATM = .TRUE.
C
      KINTRP = 1
      KINTAD = KINTRP + (MAXTYP + 1)/IRAT
      KLAST  = KINTAD + (MAXTYP + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('GET1PR',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      KFREE = 1
      LFREE = LWRK
      CALL PR1IN1(WORK(KLAST),KFREE,LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,WORD,IORDER,IDUMMY,TRIANG,
     &            PROPRI,IPRINT,SINTMA,NCOMP,TOFILE,MTFORM,DOINT,
     &            DUMMY,.FALSE.,DUMMY)
      CALL GPCLOSE(LUPROP,'KEEP')
      RETURN
      END
C  /* Deck get1in */
      SUBROUTINE GET1IN(SINTMA,WORD,NCOMP,WORK,LWORK,LABINT,INTREP,
     &                  INTADR,MPQUAD,TOFILE,KPATOM,TRIMAT,EXPVAL,
     &                  EXP1VL,DENMAT,NPRINT)
C
C     K.Ruud, March 1992. This routine generates the nescessary input
C     for getting first-order magnetic derivatives from HERMIT
C
C     Some comments about the input variables:
C     KPATOM = 1 => Spin-spin couplings
C     KPATOM = 2 => The nuclear shieldings (or electric field gradients)
C                   are calculated for one symmetry
C                   independent nuclei at a time (loop is outside subroutine)
C
C     MPQUAD: For spin-spin couplings: Number of quadrature points
C     MPQUAD: For shieldings and KPATOM=2: which nuclei for which shieldings
C             are to be calculated
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxmom.h"
#include "nuclei.h"
#include "cbiher.h"
#include "dorps.h"
#include "spnout.h"
#include "inftap.h"
#include "maxorb.h"
#include "infpar.h"
C
      LOGICAL TOFILE, TRIMAT, DOINT(2,2), EXP1VL
      CHARACTER WORD*7, LABINT(*)*8
      DIMENSION WORK(LWORK), SINTMA(*), INTREP(*), INTADR(*), EXPVAL(*),
     &          DENMAT(*)
C
      CALL QENTER('GET1IN')
      TRIANG = TRIMAT
      PROPRI = NPRINT .GT. 5
C
      DOINT(1,1) = .TRUE.
      DOINT(2,2) = .TRUE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
C
C     Although not obvious, this is a test of whether we are calculating
C     spin-spin couplings
C
      IF (KPATOM .EQ. 1) THEN
         NPATOM = 0
         DO 10 I = 1, NUCIND
            IF (DOPERT(I,2)) THEN
               NPATOM = NPATOM + 1
               IPATOM(NPATOM) = I
            END IF
 10      CONTINUE
         IF (NPATOM .EQ. 0) THEN
            ALLATM = .TRUE.
            NPATOM = NUCDEP
         ELSE
            ALLATM = .FALSE.
         END IF
         KPATOM = NPATOM
      ELSE IF (KPATOM .EQ. 2) THEN
C
C     We have a loop outside over all symmetry independent nuclei (needed
C     in order to save memory in direct and parallel calculations)
C     K.Ruud, nov.-94
C
         NPATOM = 1
         IPATOM(1) = MPQUAD
         ALLATM = .FALSE.
      ELSE
         ALLATM = .TRUE.
      END IF
      IF (WORD(1:6) .EQ. 'CARMOM') IORCAR = MPQUAD
      IF (WORD(1:6) .EQ. 'LONSOL') THEN
         IORDER = MPQUAD
      ELSE
         IORDER = IDUMMY
      END IF
C
      LUPBKP = LUPROP
      IF (LUPROP .LE. 0 .AND. TOFILE .AND. MYNUM .EQ. MASTER) THEN
         CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     $        'UNFORMATTED',IDUMMY,.FALSE.)
      END IF
      KFREE = 1
      LFREE = LWORK
      CALL PR1IN1(WORK,KFREE,LFREE,INTREP,INTADR,LABINT,
     &            WORD,IORDER,MPQUAD,TRIANG,PROPRI,
     &            NPRINT,SINTMA,NCOMP,TOFILE,'TRIANG',DOINT,EXPVAL,
     &            EXP1VL,DENMAT)
      IF (TOFILE .AND. MYNUM .EQ. MASTER) THEN
         IF (LUPBKP .LE. 0) CALL GPCLOSE(LUPROP,'KEEP')
      END IF
      CALL QEXIT('GET1IN')
      RETURN
      END
C  /* Deck spntra */
      SUBROUTINE SPNTRA(SPNDSO,EXPVAL,CSTRA,SCTRA,
     &                  NCOMP,LABINT,IPRINT)
#include "implicit.h"
#include "maxmom.h"
#include "mxcent.h"
      LOGICAL SAME
      DIMENSION EXPVAL(NCOMP), CSTRA(*), SCTRA(*), SPNDSO(MXCOOR,MXCOOR)
      CHARACTER DSONUM*4, LABINT(9*MXCENT)*8, WORD*8
      ICOMP = 0
      DO 10 ICOMP = 1, NCOMP
         WORD = LABINT(ICOMP)
         DSONUM = WORD(5:)
         CALL SCOMTR(DSONUM,IROW,ICOL)
         SPNDSO(IROW,ICOL) = -EXPVAL(ICOMP)
         SPNDSO(ICOL,IROW) = -EXPVAL(ICOMP)
 10   CONTINUE
      IF (IPRINT .GE. 3) THEN
         CALL HEADER('DSO part of spin-spin coupling tensors',-1)
         CALL PRIHES(SPNDSO,'SPNSPN',CSTRA,SCTRA)
      END IF
      RETURN
      END
C  /* Deck efntra */
      SUBROUTINE EFNTRA(STAVEC,EXPVAL,CSTRA,SCTRA,NCOMP,LABINT,IPRINT)
#include "implicit.h"
#include "maxmom.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      DIMENSION STAVEC(NCOMP), EXPVAL(NCOMP), CSTRA(*), SCTRA(*)
      CHARACTER NSTNUM*2, LABINT(9*MXCENT)*8, WORD*8
#include "nuclei.h"
#include "symmet.h"
      DO 10 ICOMP = 1, NCOMP
         WORD = LABINT(ICOMP)
         NSTNUM = WORD(6:7)
         NUMBER = (ICHAR(NSTNUM(1:1)) - ICHAR('0'))*10
         IROW = NUMBER + ICHAR(NSTNUM(2:2)) - ICHAR('0')
         STAVEC(IROW) = EXPVAL(ICOMP)
 10   CONTINUE
      CALL DIASPR(STAVEC,EXPVAL,CSTRA,SCTRA,IPRINT)
      RETURN
      END
C  /* Deck shitra */
      SUBROUTINE SHITRA(STAVEC,EXPVAL,CSTRA,SCTRA,NCOMP,NROW,LABINT,
     &                  INWORD,IPRINT)
#include "implicit.h"
#include "maxmom.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      DIMENSION STAVEC(NROW,3), EXPVAL(NCOMP), CSTRA(*), SCTRA(*)
      CHARACTER NSTNUM*3, NSTDIR*1, LABINT(9*MXCENT)*8, WORD*8, INWORD*7
#include "sigma.h"
#include "spinro.h"
#include "nuclei.h"
#include "symmet.h"
      DO 10 ICOMP = 1, NCOMP
         WORD = LABINT(ICOMP)
         NSTNUM = WORD(:3)
         NSTDIR = WORD(8:8)
         CALL NCOMTR(NSTNUM,NSTDIR,IROW,MCOMP)
         IF (INWORD .EQ. 'NSTCGOS') THEN
            ELSPRD(IPTAX(MCOMP,2),IROW) = ELSPRD(IPTAX(MCOMP,2),IROW)
     &                                  - EXPVAL(ICOMP)
         ELSE
            SIGMAD(IPTAX(MCOMP,2),IROW) = EXPVAL(ICOMP)
            ELSPRD(IPTAX(MCOMP,2),IROW) = EXPVAL(ICOMP)
         END IF
         STAVEC(IROW,MCOMP) = EXPVAL(ICOMP)
 10   CONTINUE
      IF (IPRINT .GT. 2) THEN
         IF (INWORD .NE. 'NSTCGOS') THEN
            CALL HEADER('Diamagnetic part of shielding',-1)
            CALL FCPRI(SIGMAD,'SIGMA',CSTRA,SCTRA)
         ELSE
            CALL HEADER('Diamagnetic part of spin-rotation',-1)
            CALL FCPRI(ELSPRD,'SIGMA',CSTRA,SCTRA)
         END IF
      END IF
      RETURN
      END
C  /* Deck scomtr */
      SUBROUTINE SCOMTR(DSONUM,IROW,ICOL)
#include "implicit.h"
      CHARACTER DSONUM*4, NUMCAR*1
      DIMENSION NUMBER(4)
C
      IZERO = ICHAR('0')
      DO 10 ICOUNT = 1, 4
         NUMCAR = DSONUM(ICOUNT:ICOUNT)
         NUMBER(ICOUNT) = ICHAR(NUMCAR) - IZERO
 10   CONTINUE
      IROW = 10*NUMBER(1) + NUMBER(2)
      ICOL = 10*NUMBER(3) + NUMBER(4)
      RETURN
      END
C  /* Deck ncomtr */
      SUBROUTINE NCOMTR(NSTNUM,NSTDIR,IROW,MCOMP)
#include "implicit.h"
      CHARACTER NSTNUM*3, NSTDIR*1
C
      NUMBER = (ICHAR(NSTNUM(1:1)) - ICHAR('0'))*100
      NUMBER = NUMBER + (ICHAR(NSTNUM(2:2)) - ICHAR('0'))*10
      IROW   = NUMBER + ICHAR(NSTNUM(3:3)) - ICHAR('0')
      IF (NSTDIR .EQ. 'X') MCOMP = 1
      IF (NSTDIR .EQ. 'Y') MCOMP = 2
      IF (NSTDIR .EQ. 'Z') MCOMP = 3
      RETURN
      END
C  /* Deck sigadd */
      SUBROUTINE SIGADD(AMAT)
#include "implicit.h"
#include "mxcent.h"
      DIMENSION AMAT(3,MXCOOR)
#include "nuclei.h"
#include "sigma.h"
      NCOORD = 3*NUCDEP
      DO 100 I = 1, 3
         DO 200 J = 1, NCOORD
            SIGMAT(I,J) = SIGMAT(I,J) + AMAT(I,J)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck addsus */
      SUBROUTINE ADDSUS(AMAT)
#include "implicit.h"
#include "mxcent.h"
      DIMENSION AMAT(MXCOOR,MXCOOR)
#include "suscpt.h"
#include "nuclei.h"
C
      DO 100 I = 1, 3*NUCDEP
         DO 200 J = 1, 3*NUCDEP
            SUSTOT(J,I) = SUSTOT(J,I) + AMAT(J,I)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck spninp */
      SUBROUTINE SPNINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "abainf.h"
      PARAMETER (NTABLE = 16, D0 = 0.0D0)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IPOINT(MXCENT)
#include "dorps.h"
#include "spnout.h"
C
      DATA TABLE /'.ABUNDA','.PRINT ','.SELECT','.NODSO ','.NOSD  ',
     &            '.NOFC  ','.NOPSO ','.SD+FC ','.SDXFC ','.ISOTOP',
     &            '.COUPLI','.SOS   ','.SOSOCC','.SINGST','.TRIPST',
     &            '.SOSOCS'/
C
      NEWDEF = (WORD .EQ. '*SPIN-S')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),I
                  END IF
 200           CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &               '" not recognized in SPIN-S.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Input keywords in SPNINP')
 1             CONTINUE
                  READ (LUCMD,*) ABUND
                  IF (ABUND .EQ. 1.0D0) ICHANG = ICHANG - 1
               GOTO 100
 2             CONTINUE
                  READ (LUCMD,*) ISPPRI
                  IF (ISPPRI .EQ. 0) ICHANG = ICHANG - 1
               GOTO 100
 3             CONTINUE
                  DO 1030 I = 1, MXCENT
                     DOPERT(I,2) = .FALSE.
 1030             CONTINUE
                  READ (LUCMD,*) NPERT
                  READ (LUCMD, *) (IPOINT(I),I=1,NPERT)
                  CALL HEADER('Only the spin-spin couplings between'//
     &               ' the following nuclei will be calculated:',1)
                  WRITE (LUPRI,*) ' ',(IPOINT(I), I = 1 , NPERT)
                  WRITE (LUPRI,'()')
                  DO 1031 I = 1, NPERT
                     DOPERT(IPOINT(I),2) = .TRUE.
 1031             CONTINUE
                  DOSELE = .TRUE.
               GOTO 100
 4             CONTINUE
                  DODSO = .FALSE.
               GOTO 100
 5             CONTINUE
                  DOSD  = .FALSE.
               GOTO 100
 6             CONTINUE
                  DOFC  = .FALSE.
               GOTO 100
 7             CONTINUE
                  DOPSO = .FALSE.
               GOTO 100
 8             CONTINUE
                  DOSDFC= .TRUE.
                  DOSD  = .FALSE.
                  DOFC  = .FALSE.
               GOTO 100
 9             CONTINUE
                  ANISON = .TRUE.
               GOTO 100
 10            CONTINUE !        .ISOTOP
                  SPNISO = .TRUE.
                  READ (LUCMD, *) (ISOTPS(IS), IS = 1, NATOMS)
                  ABUND = D0
               GOTO 100
 11            CONTINUE !        .COUPLI
                  READ (LUCMD, *)  NUCSPI
                  READ (LUCMD, *) (IPOINT(IS), IS = 1, NUCSPI)
                  IF (NUCSPI .GT. 0) THEN
                     DO IS = 1, MXCENT
                        NCSPNI(IS) = .FALSE.
                     END DO
                  END IF
                  DO IS = 1, NUCSPI
                     IF (IPOINT(IS) .LE. 0 .OR.
     &                   IPOINT(IS) .GT. NUCIND) THEN
                        WRITE (LUPRI,*) 'Non-existent '//
     &                     'symmetry-dependent nucleus specified for'//
     &                     ' .COUPLING NUCLEUS'
                        CALL QUIT('Inconsistent input in *SPIN-SPIN '//
     &                            'module')
                     ELSE
                        NCSPNI(IPOINT(IS)) = .TRUE.
                     END IF
                  END DO
               GOTO 100
CPFP: 27/03-2006: SOS for spin-spin coupling included
 12            CONTINUE
                  SOSSPN = .TRUE.
               GOTO 100
 13            CONTINUE
                  SOSSPN = .TRUE.
                  SOSOCC = .TRUE.
               GOTO 100
 14            CONTINUE
                  READ (LUCMD,*) NSTATS
               GOTO 100
 15            CONTINUE
                  READ (LUCMD,*) NSTATT
               GOTO 100
 16            CONTINUE
                  SOSSPN = .TRUE.
                  SOSOCS = .TRUE.
                  READ (LUCMD,*) NSTATI,NSTATF,NITRST
               GOTO 100
Cend-PFP
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &               '" not recognized in SPIN-S'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT ('Illegal keyword in SPIN-S')
            END IF
      END IF
 300  CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for SPIN-S:',0)
         IF (ABUND .NE. 1.0D0) THEN
            WRITE (LUPRI,'(A,F8.3,A)')
     &            ' Spin-spin couplings printed for atoms with '//
     &            'abundance greater than :', ABUND,' %'
         END IF
         IF (ISPPRI .NE. 0) THEN
            WRITE (LUPRI,'(A,I5)')
     &            ' Print level in spin-spin output routine: ', ISPPRI
         END IF
         IF (.NOT. DODSO) WRITE (LUPRI,'(A)')
     &      ' No diamagnetic spin-orbit contribution calculated'
         IF (.NOT. DOSD)  WRITE (LUPRI,'(A)')
     &      ' No spin-dipole contribution calculated'
         IF (.NOT. DOFC) WRITE (LUPRI,'(A)')
     &      ' No Fermi contact contribution calculated'
         IF (.NOT. DOPSO) WRITE (LUPRI,'(A)')
     &      ' No paramagnetic spin-orbit contribution calculated'
         IF (DOSDFC) WRITE (LUPRI,'(A)')
     &      ' Only the sum of SD and FC will be calculated'
         IF (NUCSPI .GT. 0) 
     &      WRITE (LUPRI,'(A)')' The spin-spin couplings of selected '//
     &           'nuclei with all the others will be calculated'
CPFP: 27/03-2006: SOS for spin-spin coupling included
         IF (SOSSPN) WRITE (LUPRI,'(A)')
     &      ' Orbital analysis of spin-spin couplings'
         IF (SOSOCC) WRITE (LUPRI,'(A,A)')
     &      ' Only contributions from pairs of occupied ',
     &      ' orbitals are printed'
         IF (NSTATS .GT. 0) 
     &      WRITE (LUPRI,'(A,I6,A)')' Only the contributions from the ',
     &           NSTATS,' lowest singlet states are included'
         IF (NSTATT .GT. 0) 
     &      WRITE (LUPRI,'(A,I6,A)')' Only the contributions from the ',
     &           NSTATT,' lowest triplet states are included'
Cend-PFP
      END IF
      RETURN
      END
C  /* Deck spnini */
      SUBROUTINE SPNINI
C
C     Initialize SPNOUT
C
#include "implicit.h"
#include "mxcent.h"
#include "abainf.h"
#include "spnout.h"
C
      ABUND  = 1.0D0
      ISPPRI = IPRDEF
      DODSO  = .TRUE.
      DOSD   = .TRUE.
      DOFC   = .TRUE.
      DOPSO  = .TRUE.
      DOSDFC = .FALSE.
      DOSELE = .FALSE.
      ANISON = .FALSE.
      SPNISO = .FALSE.
      NUCSPI = 0
      DO IS = 1, MXCENT
         NCSPNI(IS) = .TRUE.
      END DO
CPFP: 27/03-2006: SOS for spin-spin coupling included
      SOSSPN= .FALSE.
      SOSOCC= .FALSE.
      NSTATS= 0
      NSTATT= 0
      SOSOCS= .FALSE.
      NSTATI= 0
      NSTATF= 0
      NITRST= 0
Cend-PFP
      RETURN
      END
C  /* Deck trpinp */
      SUBROUTINE TRPINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 12)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "cbitrp.h"
C
C Used from common blocks:
C  gnrinf.h : THR_REDFAC
C
#include "gnrinf.h"
#include "abainf.h"
#include "anrinf.h"
#include "dorps.h"
#include "nuclei.h"
#include "spnout.h"
      DATA TABLE /'.SKIP  ', '.PRINT ','.MAX IT','.THRESH',
     *            '.MAXRED', '.MAXPHP','.NORHS ','.NORSP ',
     *            '.OPTORB', '.INTPRI','.XXXXXX','.STOP  '/
C
      NEWDEF = (WORD .EQ. '*TRPRSP')
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
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/4A/)') ' Keyword "',WORD,
     &            '" not recognized in ',WORD1
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in TRPINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRTRP
                  IF (IPRTRP .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) MAXTRP
                  IF (MAXTRP .EQ. 40) ICHANG = ICHANG - 1
               GO TO 100
    4          CONTINUE ! .THRESH
                  READ (LUCMD,*) THRTRP
                  IF (THR_REDFAC .GT. 0.0D0) THEN
                     WRITE (LUPRI,'(4A,1P,D10.2)') '@ INFO ',
     &               WORD1, WORD,
     &               ' threshold multiplied with general factor',
     &               THR_REDFAC
                     THRTRP = THRTRP*THR_REDFAC
                  END IF
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) MXRM
                  IF (MXRM .EQ. 400) ICHANG = ICHANG - 1
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) MXPHP
                  IF (MXPHP .EQ. 0) ICHANG = ICHANG - 1
               GO TO 100
    7             SKIP = .TRUE.
               GO TO 100
    8          CONTINUE
                  TRPRSP = .FALSE.
               GO TO 100
    9          CONTINUE
                  OOTV   = .TRUE.
               GO TO 100
   10          CONTINUE
               READ (LUCMD,'(I5)') INTPRI
               IF (INTPRI .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   11             CONTINUE
               GO TO 100
   12             CUT    = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in TRPINP'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in TRPINP')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for TRPRSP:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' Triplet response skipped in this run.'
         ELSE
            IF (IPRTRP .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)')
     &              ' Print level in TRPRSP        :',IPRTRP
            END IF
            IF ( INTPRI .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)')
     &              ' Print level of integrals     :', INTPRI
            END IF
            WRITE (LUPRI,'(A,1P,E9.2)')
     &              ' Threshold in triplet response:',THRTRP
            IF (MAXTRP .NE. 40) THEN
               WRITE(LUPRI,'(A,I5)')' Maximum iterations in TRPRSP :',
     &                               MAXTRP
            END IF
            IF (MXRM .NE. 400) THEN
               WRITE (LUPRI,'(A,I5)')' Maximum size of reduced space :',
     &                                MXRM
            END IF
            IF (MXPHP .NE. 0) THEN
               WRITE (LUPRI,'(A,I5)')' Maximum size of explicit '//
     &                               'configuration Hessian :', MXPHP
            END IF
            IF (.NOT. TRPRSP) WRITE (LUPRI,'(A)')
     &           ' No triplet response equations solved in this run.'
            IF (OOTV) WRITE(LUPRI,'(A)')
     &           ' Optimal orbital trial vectors used in this run.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after TRPDRV.'
            END IF
         END IF
      ELSE IF (TRPRSP) THEN
            WRITE (LUPRI,'(A,1P,E9.2)')
     &            ' Threshold in triplet response:',THRTRP
      END IF
      RETURN
      END
C  /* Deck trpini */
      SUBROUTINE TRPINI
C
C     Initialize /CBITRP/
C
C Used from common blocks:
C gnrinf.h : THR_REDFAC
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "gnrinf.h"
#include "abainf.h"
#include "cbitrp.h"
C
      IPRTRP = IPRDEF
      INTPRI = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      THRTRP = 1.D-04
      MAXTRP = 60
      MXRM   = 400
      MXPHP  = 0
      TRPRSP = .TRUE.
      IF (THR_REDFAC .GT. 0.0D0) THEN
         WRITE (LUPRI,'(2A,1P,D10.2)') '@ INFO ',
     &   ' *TRPRES default threshold multiplied with general factor',
     &   THR_REDFAC
         THRTRP = THRTRP*THR_REDFAC
      END IF
      RETURN
      END
C  /* Deck trpdrv */
      SUBROUTINE TRPDRV(WORK,LWORK,PASS)
C
C     Driver routine for the calculation of triplet operator
C     contributions (FC and SD) to spin-spin-relaxation terms
C
C     K.Ruud, Sept. 1992
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      LOGICAL PASS
      CHARACTER*8 SDLAB(9*MXCENT), FCLAB(MXCENT)
      DIMENSION WORK(LWORK)
#include "gnrinf.h"
#include "cbitrp.h"
#include "spnout.h"
#include "symmet.h"
#include "nuclei.h"
#include "infdim.h"
#include "inforb.h"
#include "inflin.h"
#include "infrsp.h"
#include "infvar.h"
#include "wrkrsp.h"
#include "abainf.h"
#include "infsop.h"
C
      IF (SKIP) RETURN
      IF (IPRTRP .GT. 0) CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRTRP .GT. 0) WRITE (LUPRI,'(A/)')
     &    ' <<<<<<<<<< Output from TRPDRV >>>>>>>>>> '
C
C     Allocate necessary work space to TRPRH1
C
C
      IF (ABASOP) THEN
         LUDV   = NORBT * NORBT
         LPVX   = 2*LPVMAT
C The factor 2 is because of triplet properties
      ELSE
         LUDV   = N2ASHX
         LPVX   = 0
      ENDIF
C
C
      KSDREP = 1
      KFCREP = KSDREP + (9*NUCDEP + 1)/IRAT
      KINTAD = KFCREP + (NUCDEP   + 1)/IRAT
      KINDX  = KINTAD + (9*NUCDEP + 1)/IRAT
      KCMO   = KINDX  + LCINDX
      KUDV   = KCMO   + NCMOT
      KPVX   = KUDV   + LUDV
      KLAST  = KPVX   + LPVX
      LWRK   = LWORK  - KLAST + 1
C
C     Allocate necessary work space to TRPRS1
C
      KZERNR = 1
      KNABAT = KZERNR + (10*NUCDEP + 1)/IRAT
      KLAST2   = KNABAT + (10*NUCDEP + 1)/IRAT
      LWRK2 = LWORK - KLAST2 + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('TRPDRV','RHS',KLAST,LWORK)
      IF (KLAST2 .GT. LWORK) CALL STOPIT('TRPDRV','RSP',KLAST2,LWORK)
C
C     Call routine for constructing right-hand sides
C
      CALL TRPRH1(SDLAB,WORK(KSDREP),FCLAB,WORK(KFCREP),
     &            WORK(KINTAD),WORK(KINDX),WORK(KCMO),WORK(KUDV),
     &            WORK(KPVX),WORK(KLAST),LWRK)
C
C     Call routine for solving response equations
C
      IF (TRPRSP) THEN
         CALL TRPRS1(WORK(KZERNR),WORK(KNABAT),
     &               WORK(KLAST2),LWRK2)
      END IF
      IF (IPRTRP .GT. 0) CALL TIMER ('TRPRSP',TIMEIN,TIMOUT)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &        ' Program stopped after TRPRSP as required.'
         CALL QUIT(' ***** End of ABACUS (in TRPRHS) *****')
      END IF
      RETURN
      END
C  /* Deck trprh1 */
      SUBROUTINE TRPRH1(SDLAB,ISDREP,FCLAB,IFCREP,INTADR,XINDX,
     &                  CMO,UDV,PV,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      PARAMETER (DM1 = -1.0D0)
      CHARACTER OPLBL(8)*8, SDLAB(9*MXCENT)*8, FCLAB(MXCENT)*8,
     &          LABEL(9*MXCENT)*8
      LOGICAL   OLDDX, FOUND
      DIMENSION ISDREP(*), IFCREP(*), WORK(LWORK), INTADR(*),
     &          XINDX(*), CMO(*), UDV(*), PV(*)
C
#include "abainf.h"
#include "cbitrp.h"
#include "spnout.h"
#include "nuclei.h"
#include "inftap.h"
#include "infdim.h"
#include "inforb.h"
#include "inflin.h"
#include "infrsp.h"
#include "infvar.h"
#include "symmet.h"
#include "wrkrsp.h"
#include "gdvec.h"
#include "abares.h"
#include "infsop.h"
C
      DATA OPLBL /'EXOPSYM1','EXOPSYM2','EXOPSYM3','EXOPSYM4',
     &            'EXOPSYM5','EXOPSYM6','EXOPSYM7','EXOPSYM8'/
C
C     Allocate work space for necessary information to be read by RSPMC
C
      CALL RD_SIRIFC('CMO',FOUND,CMO,WORK,LWORK)
      IF (.NOT.FOUND) CALL QUIT('TRPRH1 error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         IF (LWORK.LT.NNASHX)CALL STOPIT('TRPRH1','UDV',LWORK,NNASHX)
         CALL RD_SIRIFC('DV',FOUND,WORK,WORK,LWORK)
         IF (.NOT.FOUND)
     &      CALL QUIT('TRPRH1 error: DV not found on SIRIFC')
            CALL DSPTSI(NASHT,WORK,UDV)
      END IF
C
      KSYMOP = LSYMRF
      IPRRSP = IPRTRP
      TRPLET = .TRUE.
C
      IF (.NOT. ABASOP) CALL GETCIX(XINDX,LSYMRF,LSYMRF,WORK,LWORK,0)
C
C
C     SOPPA :
C     reads second order one particle density matrix
C     and doubles correlation coefficients.
C
      IF (ABASOP) THEN
C
C        Initialize XINDX
C
         CALL DZERO(XINDX,LCINDX)
C
C        Find address array's for SOPPA calculation
C
         CALL SET2SOPPA(XINDX(KABSAD),XINDX(KABTAD),
     &                  XINDX(KIJSAD),XINDX(KIJTAD),
     &                  XINDX(KIJ1AD),XINDX(KIJ2AD),
     &                  XINDX(KIJ3AD),XINDX(KIADR1))
C
C
         REWIND (LUSIFC)
         IF (CCPPA) THEN
            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
         CALL READT (LUSIFC,LPVMAT,PV)
C
         CALL TRPKAP(PV(1),PV(1+LPVMAT),XINDX(KIADR1),WORK)
C
         IF (IPRTRP.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')
     &            ' TRPRH1 : CCSD singlet and triplet correlation',
     &            ' coefficients'
            ELSE
               WRITE(LUPRI,'(/A)')
     &            ' TRPRH1 : MP2 singlet and triplet correlation',
     &            ' coefficients'
            ENDIF
            CALL OUTPUT(PV,1,LPVMAT,1,2,LPVMAT,2,1,LUPRI)
         END IF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
         CALL READT (LUSIFC,NORBT*NORBT,UDV)
C
C        UDV contains the MP2/CCSD one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
         IF (IPRTRP.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' TRPRH1 : CCSD density'
            ELSE
               WRITE(LUPRI,'(/A)')' TRPRH1 : MP2 density'
            END IF
            CALL OUTPUT(UDV,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
C
         CALL SOPUDV(UDV)
C
      END IF
C
      NREC = 0
      IF (DOSD) THEN
         IF (DODRCT) THEN
            NPATOM = 2
            ISTRT  = 1
            DO I = 1, NUCIND
               NCOMP  = 0
               CALL GET1IN(DUMMY,'SPIN-DI',NCOMP,WORK,LWORK,LABEL,
     &                     ISDREP,INTADR,I,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &                     .FALSE.,DUMMY,INTPRI)
               ITYP  = 0
               ITYP2 = 0
               DO J = 0, MAXREP
                  DO K = 1, NUCIND
                     DO L = 1, 3
                        ISCOOR = IPTCNT(3*(K - 1) + L,J, 2)
                        IF (ISCOOR .GT. 0) THEN
                           DO M = 1, 3
                              ITYP = ITYP + 1
                              IF (I .EQ. K) THEN
                                 ITYP2 = ITYP2 + 1
                                 SDLAB(ITYP) = LABEL(ITYP2)
                              END IF
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO
            NREC = NREC + 6*NUCDEP
         ELSE
            NPATOM = 0
            NCOMP = 0
            CALL GET1IN(DUMMY,'SPIN-DI',NCOMP,WORK,LWORK,SDLAB,
     &                  ISDREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,
     &                  DUMMY,.FALSE.,DUMMY,INTPRI)
            NREC = NREC + 6*NUCDEP
         END IF
      END IF
      IF (DOFC) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'FERMI C',NCOMP,WORK,LWORK,FCLAB,
     &               IFCREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,
     &               DUMMY,.FALSE.,DUMMY,INTPRI)
         NREC = NREC + NUCDEP
      END IF
      IF (DOSDFC) THEN
         NPATOM = 0
         NCOMP  = 0
         CALL GET1IN(DUMMY,'SDFC   ',NCOMP,WORK,LWORK,SDLAB,
     &               ISDREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &               .FALSE.,DUMMY,INTPRI)
         NREC = NREC + 9*NUCDEP
      END IF
      CALL GPOPEN(LUGDT,ABAGDT,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
chj 7-mar-06: I don't think it matters if old or not ... /hjaaj
chj           and the records are anyway written below
chj   IF (OLDDX) THEN
chj      WRITE (LUPRI,'(/A)') ' Old LUGDT file opened in TRPRH1.'
chj   ELSE
chj      CALL DZERO(WORK,NVARMA)
chj      DO I = 1, NREC
chj         CALL WRITDX(LUGDT,I,IRAT*NVARMA,WORK)
chj      END DO
chj   END IF
C
C     Construct right-hand side
C
      DO 100 ISYM = 1, NSYM
         CALL ABAVAR(ISYM,TRPLET,IPRTRP,WORK,LWORK)
         IREFSY = LSYMRF
         NCREF  = NCONRF
         KSYMST = LSYMST
         KSYMOP = LSYMPT
         KZWOPT = NWOPT
         IF ((NASHT .EQ. 0).AND.(.NOT.ABASOP)) THEN
            KZCONF = 0
         ELSE IF (NASHT .EQ. 1) THEN
            KZCONF = 1
         ELSE
            KZCONF = NCONST
         END IF
         KZVAR  = KZWOPT + KZCONF
         KZYWOP = 2*KZWOPT
         KZYCON = 2*KZCONF
         KZYVAR = 2*KZVAR
C
C     Construct the appropriate right-hand sides for this symmetry
C
         NREC = NTRVEC(ISYM)
         DO 200 LREC = 1, NREC
            IREC  = ITRREC(LREC,ISYM)
            ICOOR = ITRCOR(LREC,ISYM)
            IF (ICOOR .LT. 0) THEN
               LABEL(1) = FCLAB(ABS(ICOOR))
            ELSE
               LABEL(1) = SDLAB(ICOOR)
               ICOOR1 = ICOOR/3
               ICOOR2 = MOD(ICOOR,3)
            END IF
C
            CALL GETGPV(LABEL(1),DUMMY,DUMMY,CMO,UDV,PV,XINDX,
     &                  ANTSYM,WORK,LWORK)
C
            TRPLAB(IREC) = LABEL(1)
            CALL DSCAL(KZVAR,DM1,WORK,1)
            IF (IPRTRP .GE. 05) THEN
               IF (ICOOR .GT. 0) THEN
                  WRITE (LUPRI,'(//A,2I2,A,I2,A)')
     &                 ' SD-Gradient for Coordinate',
     &                 ICOOR1,ICOOR2,
     &                 ', and symmetry',ISYM,':'
               ELSE
                  WRITE (LUPRI,'(//A,I3,A,I2,A)')
     &                 ' FC-Gradient for Atom',
     &                 ABS(ICOOR),', and symmetry',ISYM,':'
               END IF
               IF (IPRTRP .GE. 20) THEN
                  WRITE (LUPRI,'(/,A,/)')
     &                 ' CSF part of gradient'
                  WRITE (LUPRI,'(6F12.6)')
     &                 (WORK(J), J = 1, NCONST)
               END IF
               WRITE (LUPRI,'(//,A)')
     &              ' Orbital part of gradient'
               IF (IPRTRP .GE. 20) THEN
                  PRFAC = 0.0D0
               ELSE
                  PRFAC = 0.1D0
               END IF
               CALL PRKAP(NWOPPT,WORK(NCONST + 1),PRFAC,LUPRI)
            END IF
            CALL WRITDX(LUGDT,IREC,IRAT*NVARPT,WORK)
 200     CONTINUE
 100  CONTINUE
C
      CALL GPCLOSE(LUGDT,'KEEP')
C
      RETURN
      END
C  /* Deck trprs1 */
      SUBROUTINE TRPRS1(ZERNRM,NABATY,WORK,LWORK)
#include "implicit.h"
#include "iratdef.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      LOGICAL   INTTRS, CICLC, HFCLC, TRPCLC, ZERNRM, EXCLC, OLDDX
      PARAMETER (D10 = 10.0D0, DPM5 = -0.50D0)
      DIMENSION WORK(LWORK), NABATY(*), ZERNRM(*)
C
#include "cbitrp.h"
#include "spnout.h"
#include "nuclei.h"
#include "inftap.h"
#include "inflin.h"
#include "infdim.h"
#include "inforb.h"
#include "infvar.h"
#include "gdvec.h"
#include "abares.h"
#include "chrnos.h"
C
      HFCLC  = NASHT .LE. 1
      IF (HFCLC) NCONF = 1
      EXCLC  = .FALSE.
      TRPCLC = .TRUE.
      CICLC  = .FALSE.
      EXVAL  = 0.0D0
      NEXVAL = 1
C
C     Open direct access files
C
      LUGDVE = -1
      LUSOVE = -1
      LUREVE = -1
      CALL GPOPEN(LUGDT,ABAGDT,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
C     LURDT: ...changed NREC to 2*NREC to cover residuals
      CALL GPOPEN(LURDT,ABARDT,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
      DO 50 ISYM = 1, NSYM
         NREC = NTRVEC(ISYM)
         IF (NREC .GT. 0) THEN
            CALL ABAVAR(ISYM,TRPCLC,IPRTRP,WORK,LWORK)
         IF (NVARPT .EQ. 0) GOTO 50
            IF (HFCLC) NCONF = 1
C
C     Collect right-hand sides
C
            DO 10 I = 1, NREC
            NABAOP = 0
               REWIND LUGDVE
               IREC = ITRREC(I,ISYM)
C
C     If this is a restarted run, check that we haven't already solved
C     this response vector
C
               IF (IDORCT(IREC) .NE. 0) GO TO 10
C
               ICOOR = ITRCOR(I,ISYM)
C
C     If we only want to calculate the coupled FC-SD anisotropies we
C     only solve the response equation for the FC,
C     Added by K.Ruud after a suggestion by J.Vaara, March-96
C
               IF (ICOOR .GT. 0 .AND. ANISON) GO TO 10
               CALL READDX (LUGDT,IREC,IRAT*NVARPT,WORK)
               DNORM = DNRM2(NVARPT,WORK,1)
               IF (DNORM .GT. THRTRP/D10) THEN
                  NABAOP = NABAOP + 1
                  ZERNRM(I) = .FALSE.
                  CALL WRITT(LUGDVE,NVARPT,WORK)
               ELSE
                  ZERNRM(I) = .TRUE.
               END IF
               IF (IPRTRP .GT. 5) THEN
                  CALL HEADER('RHS vectors in TRPRS1 for symmetry '//
     &                        CHRNOS(ISYM)//':',-1)
                  WRITE (LUPRI,'(A,I5)') ' Coordinate:',ICOOR
                  WRITE (LUPRI,'(A,I5)') ' Record:    ',IREC
                  WRITE (LUPRI,'(A,D20.14)') ' Norm:      ', DNORM
                  IF (ZERNRM(I)) WRITE (LUPRI,'(A)')
     &             ' *Note: vector excluded because norm .le. thrclc/10'
               END IF
               IF (IPRTRP .GT. 10 .AND. .NOT. ZERNRM(I)) THEN
                  CALL OUTPUT(WORK,1,NVARPT,1,1,NVARPT,1,1,LUPRI)
               END IF
C
C     Solve response equations. In order to improve restart options, we only
C     send one vector at a time to the response program (no loss in
C     efficiency anyway). K.Ruud, Feb.-96
C
               NABATY(1) = 1
               IF (.NOT. ZERNRM(I)) THEN
                  NOCONV = .FALSE.
                  FINRES = 0.0D0
                  CALL ABARSP(CICLC,HFCLC,TRPCLC,OOTV,ISYM,EXCLC,EXVAL,
     &                        NEXVAL,NABATY,1,TRPLAB(IREC),LUGDVE,
     &                        LUSOVE,LUREVE,THRTRP,MAXTRP,IPRTRP,MXRM,
     &                        MXPHP,WORK,LWORK)
                  IF (NOCONV) THEN
                     WRITE (LUPRI,'(/3A,2F10.5)')
     &                  '@ WARNING No convergence for ',TRPLAB(IREC),
     &                  '. Last residual and threshold:', FINRES, THRTRP
                  ELSE
                     WRITE (LUPRI,'(3A,2F10.5)')
     &                  '    Convergence for ',TRPLAB(IREC),
     &                  '. Last residual and threshold:', FINRES, THRTRP
                  END IF
               END IF
C
C     Write solutions and residuals
C
               REWIND LUSOVE
               REWIND LUREVE
               IF (ZERNRM(I)) THEN
                  CALL DZERO(WORK,2*NVARPT)
               ELSE
                  CALL READT(LUSOVE,NVARPT,WORK)
                  CALL READT(LUREVE,NVARPT,WORK(NVARPT+1))
               END IF
C
C     Divide solution by 2 in accordance with ABACUS solver
C
               CALL DSCAL(NVARPT,DPM5,WORK,1)
C
               IREC2 = 2*IREC - 1
               CALL WRITDX(LURDT,IREC2,IRAT*NVARPT,WORK)
               IREC2 = 2*IREC
               CALL WRITDX(LURDT,IREC2,IRAT*NVARPT,WORK(NVARPT+1))
ckr               IDORCT(IREC) = 1
               CALL ABAWRIT_RESTART
               IF (IPRTRP .GT. 5) THEN
                  CALL HEADER('Response and residual vectors in '//
     &                'TRPRS1 for symmetry '//CHRNOS(ISYM)//':',-1)
                  WRITE (LUPRI,'(A,I5)') ' Coordinate:',ICOOR
                  WRITE (LUPRI,'(A,I5)') ' Record:    ',IREC
                  DN1 = DNRM2(NVARPT,WORK,1)
                  DN2 = DNRM2(NVARPT,WORK(1+NVARPT),1)
                  WRITE (LUPRI,'(A,2(D20.14,2X))')
     &                   ' Norms:     ',DN1,DN2
               END IF
               IF (IPRTRP .GT. 10) THEN
                  CALL OUTPUT(WORK,1,NVARPT,1,2,NVARPT,1,2,LUPRI)
               END IF
 10         CONTINUE
         END IF
 50   CONTINUE
      CALL GPCLOSE(LUGDT,'KEEP')
      CALL GPCLOSE(LURDT,'KEEP')
      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
      RETURN
      END
C  /* Deck shires */
      SUBROUTINE SHIRES(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
cLig <> added include abainf.h
#include "abainf.h"
      DIMENSION WORK(LWORK)
#include "nuclei.h"
#include "ctocdcc.h"
cLig >> changed the values of the memory for CTOCD 
      KSPMAT = 1
      IF(CTOCD .OR. CTOMAG)  THEN
        KSMMAT = KSPMAT + 9*MXCENT
        KCMAT  = KSMMAT + 9*MXCENT
      ELSE
        KSMMAT = 1
        KCMAT =  1
      ENDIF
cLig >>
      KDMAT  = KCMAT  + 9*MXCENT
      KPMAT  = KDMAT  + 9*MXCENT
      KAVET  = KPMAT  + 9*MXCENT
      KANIS  = KAVET  + NUCDEP
      KASYM  = KANIS  + NUCDEP
      KAVED  = KASYM  + NUCDEP
      KAVEP  = KAVED  + NUCDEP
      KSPAR  = KAVEP  + NUCDEP
      KAPAR  = KSPAR  + NUCDEP
      KANIS2 = KAPAR  + NUCDEP
      KASYM2 = KANIS2 + NUCDEP
      KSKEW  = KASYM2 + NUCDEP
      KSPAN  = KSKEW  + NUCDEP
      KCSTRA = KSPAN  + NUCDEP
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('SHIRES',' ',KLAST,LWORK)
cLig <> added work(kspmat) and work (ksmmat) for CTOCD
      CALL SHIRE1(WORK(KSPMAT),WORK(KSMMAT),
     &            WORK(KCMAT),WORK(KDMAT),WORK(KPMAT),WORK(KAVET),
     &            WORK(KAVED),WORK(KAVEP),WORK(KANIS),WORK(KASYM),
     &            WORK(KSPAR),WORK(KAPAR),WORK(KANIS2),WORK(KASYM2),
     &            WORK(KSKEW),WORK(KSPAN),WORK(KCSTRA),WORK(KSCTRA),
     &            IPRINT)
      RETURN
      END
C  /* Deck shire1 */
cLig <> added sftpmat and sftmmat
      SUBROUTINE SHIRE1(SFTPMAT,SFTMMAT,
     &                  TMAT,DMAT,PMAT,AVETOT,AVEDIA,AVEPAR,ANIS,ASYM,
     &                  SPAR,APAR,ANIS2,ASYM2,SKEW,SPAN,CSTRA,SCTRA,
     &                  IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      CHARACTER DIAPAR*12
cLig <> added the dimensions of sftpmat and sftmmat
      DIMENSION TMAT(3,3,MXCENT), DMAT(3,3,MXCENT), PMAT(3,3,MXCENT),
     &          AVETOT(NUCDEP), AVEDIA(NUCDEP), AVEPAR(NUCDEP),
     &          ANIS(NUCDEP), ASYM(NUCDEP), SPAR(NUCDEP), APAR(NUCDEP),
     &          ANIS2(NUCDEP), ASYM2(NUCDEP), SKEW(NUCDEP),
     &          SPAN(NUCDEP), CSTRA(*), SCTRA(*),
     &          SFTPMAT(3,3,MXCENT), SFTMMAT(3,3,MXCENT)
cLig >> added abainf.h and trkoor.h
#include "abainf.h"
#include "trkoor.h"
cLig <<
#include "symmet.h"
#include "sigma.h"
#include "cbisol.h"
C GAGORG:
#include "orgcom.h"
C CTOCD
#include "ctocdcc.h"

      LOGICAL PRTSHI(NUCIND)
      LOGICAL ALLATM
      DIMENSION IPATOM(MXCENT)
CTOCD


C

C     Get ipatom stuff (if CC some only atoms can be selected)
C
      IF (CTOMAG) THEN
         CALL SELINF(NPATOM,IPATOM,ALLATM)
      ELSE
         ALLATM = .TRUE.
      END IF

      IF (CTOMAG) THEN
C
C        Print proper title
C
         IF (EXTSOP) THEN
            CALL TITLER(
     &           'CTOCD-DZ SHIELDINGS FROM EXTERNAL LINEAR RESPONSES',
     &           ':',102)
         ELSE IF (ALCCS) THEN
C        IF (ALCCS) THEN
            CALL TITLER('CCS CTOCD-DZ CHEMICAL SHIELDINGS','*',122)
         ELSE IF (ALCC2) THEN
            CALL TITLER('CC2 CTOCD-DZ CHEMICAL SHIELDINGS','*',122)
         ELSE IF (ALCC3) THEN
            CALL TITLER('CC3 CTOCD-DZ CHEMICAL SHIELDINGS','*',122)
         ELSE IF (ALCCSD) THEN
            CALL TITLER('CCSD CTOCD-DZ CHEMICAL SHIELDINGS','*',122)
         END IF
C
C        Initialize DOSYM
C
         DO IREP = 1, MAXREP+1
            DOSYM(IREP) = .TRUE.
         END DO
C
      ELSE
         CALL TITLER('ABACUS - CHEMICAL SHIELDINGS','*',124)
      END IF
C
      IF (CTOCD .AND. (.NOT. CTOMAG)) THEN
cLig >> added CTOCD-DZ results title
        CALL HEADER('CTOCD-DZ results',1)
      ENDIF
      CALL HEADER('Shielding tensors in symmetry coordinates (ppm)',-1)
      IF (CTOCD .OR. CTOMAG) THEN
cLig >> added call to FCPRI for shieldings
         CALL FCPRI(SIGMADZ,'CTOCD',CSTRA,SCTRA)
      ELSE
         CALL FCPRI(SIGMAT,'SIGMAS',CSTRA,SCTRA)
      ENDIF
      NUCINS = NUCIND
      IF (SOLVNT) NUCINS = NUCINS - 1
C
C     Transform to non-symmetry basis
C
      IF(CTOCD .OR. CTOMAG)  THEN
cLig >> added the call to TRADIP ... for CTOCD
         CALL TRADIP(SIGMADZ,DMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
         CALL TRADIP(SIGMAR,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
         IF (DONS) THEN
            CALL DZERO(SFTPMAT,9*NUCDEP)
            CALL DZERO(SFTMMAT,9*NUCDEP)
            CALL TRADIP(SIGMASFTP,SFTPMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
            CALL TRADIP(SIGMASFTM,SFTMMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
            CALL DONSSHIFT(SFTPMAT,SFTMMAT)
            CALL DAXPY(9*NUCDEP,D1,SFTPMAT,1,PMAT,1)
            CALL DAXPY(9*NUCDEP,-D1,SFTPMAT,1,DMAT,1)
         ENDIF
         CALL DCOPY(9*NUCDEP,DMAT,1,TMAT,1)
         CALL DAXPY(9*NUCDEP,D1,PMAT,1,TMAT,1)
      ELSE
        CALL TRADIP(SIGMAT,TMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
        CALL TRADIP(SIGMAR,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
        CALL DCOPY(9*NUCDEP,TMAT,1,DMAT,1)
        CALL DAXPY(9*NUCDEP,-D1,PMAT,1,DMAT,1)
      ENDIF
C
CTOCD
      IF (ALLATM) THEN
         DO I = 1, NUCINS
            PRTSHI(I) = .TRUE.
         END DO
      ELSE
         DO I = 1, NUCINS
            PRTSHI(I) = .FALSE.
         END DO
         DO I = 1,NPATOM
            PRTSHI(IPATOM(I)) = .TRUE.
         END  DO
      END IF
CTOCD
C

      IATOM = 0
      DO 100 I = 1, NUCINS
         DO 200 ISYMOP = 0, MAXOPR
         IF (IAND(ISTBNU(I),ISYMOP).EQ.0) THEN
            IATOM = IATOM + 1
            IF (PRTSHI(I))          ! Selected atom if CTOMAG
     &      CALL SHIANA(TMAT(1,1,IATOM),DMAT(1,1,IATOM),PMAT(1,1,IATOM),
     &                  AVETOT(IATOM),AVEDIA(IATOM),AVEPAR(IATOM),
     &                  ANIS(IATOM),ASYM(IATOM),SPAR(IATOM),APAR(IATOM),
     &                  ANIS2(IATOM),ASYM2(IATOM),SKEW(IATOM),
     &                  SPAN(IATOM),NAMDEP(IATOM),NAMDPX(3*(IATOM-1)+1),
     &                  IPRINT)
         END IF
  200    CONTINUE
  100 CONTINUE
C
      CALL AROUND('Summary of chemical shieldings')
      WRITE (LUPRI,'(A/)') '@1  Definitions from J. Mason, Solid '//
     &       'state Nuc. Magn. Res. 2 (1993), 285'
      IF (NOLOND) THEN
         WRITE(LUPRI,'(A/A,3F20.12)')
     &    '@1 London orbitals (GIAOs) not used.',
     &    '@1 Gauge origin    (bohr):', (GAGORG(I),I=1,3)
      ELSE
         WRITE(LUPRI,'(A)') '@1 London orbitals (GIAOs) has been used.'
      END IF 
      WRITE (LUPRI,'(/A/A)')
     & '@1 atom   shielding       dia      para     skew '/
     &  /'     span     (aniso     asym)',
     & '@1 ----------------------------------------------'/
     &  /'------------------------------'
      IATOM = 0
      DO 300 I = 1, NUCINS
         DO 400 ISYMOP = 0, MAXOPR
         IF (IAND(ISTBNU(I),ISYMOP).EQ.0) THEN
            IATOM = IATOM + 1
            IF (AVETOT(IATOM) .GT. D0) THEN
               DIAPAR = 'paramagnetic'
            ELSE
               DIAPAR = 'diamagnetic '
            END IF
            IF (PRTSHI(I))          ! Selected atom if CTOMAG
     &      WRITE (LUPRI,'(2A,7F10.4)') '@1 ',
     &         NAMDEP(IATOM),AVETOT(IATOM),AVEDIA(IATOM),AVEPAR(IATOM),
     &         SKEW(IATOM),SPAN(IATOM),ANIS2(IATOM),ASYM2(IATOM)
         END IF
  400    CONTINUE
  300 CONTINUE
      WRITE (LUPRI,'()')
C
      CALL AROUND('Summary of chemical shieldings')
      WRITE (LUPRI,'(A/)') '@2  Definitions from Smith, Palke, and'//
     &       ' Grieg, Concepts in Mag. Res. 4 (1992), 107'
      IF (NOLOND) THEN
         WRITE(LUPRI,'(A/A,3F20.12)')
     &    '@2 London orbitals (GIAOs) not used.',
     &    '@2 Gauge origin    (bohr):', (GAGORG(I),I=1,3)
      ELSE
         WRITE(LUPRI,'(A)') '@2 London orbitals (GIAOs) has been used.'
      END IF 
      WRITE (LUPRI,'(/A/A)')
     & '@2 atom   shielding       dia      para     aniso'/
     &  /'      asym        S        A',
     & '@2 ----------------------------------------------'/
     &  /'------------------------------'
      IATOM = 0
      DO 310 I = 1, NUCINS
         DO 410 ISYMOP = 0, MAXOPR
         IF (IAND(ISTBNU(I),ISYMOP).EQ.0) THEN
            IATOM = IATOM + 1
            IF (AVETOT(IATOM) .GT. D0) THEN
               DIAPAR = 'paramagnetic'
            ELSE
               DIAPAR = 'diamagnetic '
            END IF
            IF (PRTSHI(I))          ! Selected atom if CTOMAG
     &      WRITE (LUPRI,'(2A,7F10.4)') '@2 ',
     &         NAMDEP(IATOM),AVETOT(IATOM),AVEDIA(IATOM),AVEPAR(IATOM),
     &         ANIS(IATOM),ASYM(IATOM),SPAR(IATOM),APAR(IATOM)
         END IF
  410    CONTINUE
  310 CONTINUE
      WRITE (LUPRI,'()')
      RETURN
      END
C  /* Deck shiana */
      SUBROUTINE SHIANA(TMAT,DMAT,PMAT,AVETOT,AVEDIA,AVEPAR,ANIS,ASYM,
     &                  SPAR,APAR,ANIS2,ASYM2,SKEW,SPAN,NAME,NAMEX,
     &                  IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "thrzer.h"
#include "inftap.h"
#include "codata.h"
cLig >> added include mxcent abainf and sigma
#include "mxcent.h"
#include "abainf.h"
#include "sigma.h"
cLig <<
      PARAMETER (FACTOR = 1.D6*ALPHA2)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      CHARACTER NAME*6, NAMEX(3)*8, SGN*3, LAPRPC*8, LANAME*8
cmbh
      character lbname*8,midasname*8
cmbh end
      INTEGER Z
      DIMENSION TMAT(3,3), DMAT(3,3), PMAT(3,3), SYM(3,3), ANTI(3,3),
     &          AXES(3,3), PVAL(6), WRK(3), IWRK(3), ABSP(3)
#include "chrxyz.h"
C
CTOCD
#include "ctocdcc.h"
      CHARACTER*10 LABMET
CTOCD
C
      CALL DSCAL(9,FACTOR,TMAT,1)
      CALL DSCAL(9,FACTOR,DMAT,1)
      CALL DSCAL(9,FACTOR,PMAT,1)
C
      AVETOT = (TMAT(1,1) + TMAT(2,2) + TMAT(3,3))/D3
      AVEDIA = (DMAT(1,1) + DMAT(2,2) + DMAT(3,3))/D3
      AVEPAR = (PMAT(1,1) + PMAT(2,2) + PMAT(3,3))/D3
C
C     Parameters calculated in accordance with J.Mason, Solid State
C     Nuc.Magn.Resonance, 2 (1993), 285
C
      DO 100 I = 1, 3
         DO 200 J = 1, 3
            ANTI(I,J) = (TMAT(I,J) - TMAT(J,I))/D2
            SYM (I,J) = (TMAT(I,J) + TMAT(J,I))/D2
  200    CONTINUE
  100 CONTINUE
      CALL DUNIT(AXES,3)
      IJ = 1
      DO 300 I = 1, 3
         DO 310 J = 1, I
            PVAL(IJ) = SYM(I,J)
            IJ = IJ + 1
  310    CONTINUE
  300 CONTINUE
      CALL JACO(PVAL,AXES,3,3,3,WRK,IWRK)
      PVAL(1) = PVAL(1)
      PVAL(2) = PVAL(3)
      PVAL(3) = PVAL(6)
C
      IZ = IDMAX(3,PVAL,1)
      IX = IDMIN(3,PVAL,1)
      DO 320 I = 1, 3
         IF (I .NE. IZ .AND. I .NE. IX) IY = I
 320  CONTINUE
C
      SPAN  = PVAL(IZ) - PVAL(IX)
      IF (SPAN .LE. THRZER) THEN
         SKEW  = D0
         ANIS2 = D0
         ASYM2 = D0
      ELSE
         SKEW  = D3*(AVETOT - PVAL(IY))/SPAN
         ANIS2 = PVAL(IZ) - (PVAL(IX) + PVAL(IY))/D2
         ASYM2 = (PVAL(IY) - PVAL(IX))/(PVAL(IZ) - AVETOT)
      END IF
C
C
      DO 110 I = 1, 3
         SYM(I,I) = SYM(I,I) - AVETOT
 110  CONTINUE
      CALL DUNIT(AXES,3)
      IJ = 1
      DO 325 I = 1, 3
         DO 330 J = 1, I
            PVAL(IJ) = SYM(I,J)
            IJ = IJ + 1
 330     CONTINUE
 325  CONTINUE
      CALL JACO(PVAL,AXES,3,3,3,WRK,IWRK)
      PVAL(1) = PVAL(1)
      PVAL(2) = PVAL(3)
      PVAL(3) = PVAL(6)
C
      IZ = IDAMAX(3,PVAL,1)
      CALL DSWAP(1,PVAL(3),1,PVAL(IZ),1)
      CALL DSWAP(3,AXES(1,3),1,AXES(1,IZ),1)
      IY = IDAMIN(3,PVAL,1)
      CALL DSWAP(1,PVAL(2),1,PVAL(IY),1)
      CALL DSWAP(3,AXES(1,2),1,AXES(1,IY),1)
C
      ANIS = PVAL(3) - (PVAL(2) + PVAL(1))/D2
      IF (PVAL(3) .EQ. 0) THEN
         ASYM = 0.0D0
      ELSE
         ASYM = ABS((PVAL(2) - PVAL(1))/PVAL(3))
      END IF
      APAR = SQRT(DDOT(9,ANTI,1,ANTI,1)/D2)
      SPAR = ABS(ANIS)*SQRT(D1 + (ASYM**2)/D3)
cLig >> added title for CTOCD shieldings
      IF (CTOCD .OR. CTOMAG) THEN
         CALL TITLER('CTOCD-DZ: Chemical shielding for '
     &               //NAME//':','=',1)
      ELSE
cLig <<
C
        CALL TITLER('Chemical shielding for '//NAME//':','=',1)
cLig >> added title for the Gauge Origin
      ENDIF
      IF(DONS) THEN
        CALL HEADER('Gauge Origin set on the atom',1)
      ENDIF
cLig <<
      WRITE (LUPRI,'(3(/,1X,A,F12.4,A),/,2(/,1X,A,F12.4,A))')
     &   ' Shielding constant:',AVETOT,' ppm',
     &   ' Anisotropy:        ',ANIS,  ' ppm',
     &   ' Asymmetry:         ',ASYM,  '    ',
     &   ' S parameter:       ',SPAR,  ' ppm',
     &   ' A parameter:       ',APAR,  ' ppm'
C << JK write results to DALTON.PROP file
      LANAME = NAME//'  '
      LAPRPC = 'SH-CONST'
C
CTOCD
      IF (ALCCS) THEN
         LABMET = 'DZ /CCS   '
      ELSE IF (ALCC2) THEN
         LABMET = 'DZ /CC2   '
      ELSE IF (ALCC3) THEN
         LABMET = 'DZ /CC3   '
      ELSE IF (ALCCSD) THEN
         LABMET = 'DZ/CCSD   '
      ELSE
         LABMET = 'SCF/DFT   '
      END IF
CTOCD
C
CTOCD CALL WRIPRO(AVETOT,"SCF/DFT   ",401,
      CALL WRIPRO(AVETOT, LABMET     ,401,
     *            LANAME,LAPRPC,LAPRPC,LAPRPC,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LAPRPC = 'ANISO   '
CTOCD CALL WRIPRO(ANIS  ,"SCF/DFT   ",401,
      CALL WRIPRO(ANIS  , LABMET     ,401,
     *            LANAME,LAPRPC,LAPRPC,LAPRPC,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LAPRPC = 'ASYMMETR'
CTOCD CALL WRIPRO(ASYM  ,"SCF/DFT   ",401,
      CALL WRIPRO(ASYM  , LABMET     ,401,
     *            LANAME,LAPRPC,LAPRPC,LAPRPC,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LAPRPC = 'S-PARAME'
CTOCD CALL WRIPRO(SPAR  ,"SCF/DFT   ",401,
      CALL WRIPRO(SPAR  , LABMET     ,401,
     *            LANAME,LAPRPC,LAPRPC,LAPRPC,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LAPRPC = 'A-PARAME'
CTOCD CALL WRIPRO(APAR  ,"SCF/DFT   ",401,
      CALL WRIPRO(APAR  , LABMET     ,401,
     *            LANAME,LAPRPC,LAPRPC,LAPRPC,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
C << JK
      CALL HEADER('Total shielding tensor (ppm):',1)
csonia 04/10/95
      IF (LUCME.GT.0) THEN
         WRITE (LUCME,'(2A)') 'Chemical shielding for ',NAME
         WRITE (LUCME,'(A)')  'Total shielding tensor (ppm):'
      END IF
csonia 04/10/95
C << JK write shielding tensor to DALTON.PROP
cmbh change output slightly
      LAPRPC = '_COMP  '
      midasname=NAME//'  '
      call stripblanks(midasname)
      LANAME = 'X'//LAPRPC
      LBNAME = 'X'//LAPRPC
      CALL WRIPRO(TMAT(1,1),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'X'//LAPRPC
      LBNAME = 'Y'//LAPRPC
      CALL WRIPRO(TMAT(2,1),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'X'//LAPRPC
      LBNAME = 'Z'//LAPRPC
      CALL WRIPRO(TMAT(3,1),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'Y'//LAPRPC
      LBNAME = 'X'//LAPRPC
      CALL WRIPRO(TMAT(1,2),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'Y'//LAPRPC
      LBNAME = 'Y'//LAPRPC
      CALL WRIPRO(TMAT(2,2),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'Y'//LAPRPC
      LBNAME = 'Z'//LAPRPC
      CALL WRIPRO(TMAT(3,2),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'Z'//LAPRPC
      LBNAME = 'X'//LAPRPC
      CALL WRIPRO(TMAT(1,3),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'Z'//LAPRPC
      LBNAME = 'Y'//LAPRPC
      CALL WRIPRO(TMAT(2,3),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      LANAME = 'Z'//LAPRPC
      LBNAME = 'Z'//LAPRPC
      CALL WRIPRO(TMAT(3,3),"SCF/DFT   ",402,
     *            LANAME,LBNAME,midasname,midasname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
C << JK
      WRITE (LUPRI,'(20X,3(A,13X),/)') 'Bx', 'By', 'Bz'
      DO 400 ICOOR = 1, 3
         WRITE (LUPRI,'(2X,A8,3F15.8)')
     &          NAMEX(ICOOR), (TMAT(K,ICOOR),K=1,3)
csonia 04/10/95
         IF (LUCME.GT.0) WRITE (LUCME,'(2X,A8,1P,3D23.15)')
     &          NAMEX(ICOOR), (TMAT(K,ICOOR),K=1,3)
csonia 04/10/95
  400 CONTINUE
C
      CALL HEADER('Diamagnetic and paramagnetic contributions (ppm):',1)
      WRITE (LUPRI,'(17X,3(A,9X),3X,3(A,9X),/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
      DO 500 ICOOR = 1, 3
         WRITE (LUPRI,'(2X,A8,3F11.4,3X,3F11.4)') NAMEX(ICOOR),
     &         (DMAT(K,ICOOR),K=1,3), (PMAT(K,ICOOR),K=1,3)
  500 CONTINUE
      WRITE (LUPRI,'(/1X,A,F14.6,8X,A,F14.6)')
     &          ' Diamagnetic contribution:',AVEDIA,
     &          ' Paramagnetic: ',AVEPAR
      CALL HEADER
     &   ('Antisymmetric and traceless symmetric parts (ppm):',1)
      WRITE (LUPRI,'(17X,3(A,9X),3X,3(A,9X),/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
      DO 600 ICOOR = 1, 3
         WRITE (LUPRI,'(2X,A8,3F11.4,3X,3F11.4)') NAMEX(ICOOR),
     &         (ANTI(K,ICOOR),K=1,3), (SYM(K,ICOOR),K=1,3)
  600 CONTINUE
      CALL HEADER('Principal values and axes:',1)
      DO 700 I = 1, 3
         IF (PVAL(I) .GE. D0) THEN
            SGN = '  +'
         ELSE
            SGN = '  -'
         END IF
         WRITE (LUPRI,'(2X,A,2X,F12.6,A,F8.2,A,F7.2,A,3X,3F10.6)')
     &      NAME,PVAL(I) + AVETOT,'  =',AVETOT,SGN,
     &      ABS(PVAL(I)),':', (AXES(J,I),J=1,3)
  700 CONTINUE
      RETURN
      END
C  /* Deck spires */
      SUBROUTINE SPIRES(SPNDSO,SPNPSO,SPNSD,SPNFC,SPSDFC,
     &                  WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION SPNDSO(MXCOOR,MXCOOR), SPNPSO(MXCOOR,MXCOOR),
     &          SPNSD (MXCOOR,MXCOOR), SPNFC (MXCOOR,MXCOOR),
     &          SPSDFC(MXCOOR,MXCOOR)
      DIMENSION WORK(LWORK)
C
C Allocate necessary work space
C
      KCHES  = 1
      KISOT  = KCHES  + 9*6*NUCDEP*NUCDEP
      KANIS  = KISOT  + 5*NUCDEP*(NUCDEP + 1)/2
      KASYM  = KANIS  + NUCDEP*(NUCDEP + 1)/2
      KSPAR  = KASYM  + NUCDEP*(NUCDEP + 1)/2
      KAPAR  = KSPAR  + NUCDEP*(NUCDEP + 1)/2
      KCSTRA = KAPAR  + NUCDEP*(NUCDEP + 1)/2
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('SPIRES',' ',KLAST,LWORK)
      CALL SPIRE1(SPNDSO,SPNPSO,SPNSD,SPNFC,SPSDFC,
     &            WORK(KCHES),WORK(KISOT),WORK(KANIS),WORK(KASYM),
     &            WORK(KSPAR),WORK(KAPAR),WORK(KCSTRA),WORK(KSCTRA))
      RETURN
      END
C  /* Deck spire1 */
      SUBROUTINE SPIRE1(SPNDSO,SPNPSO,SPNSD,SPNFC,SPSDFC,CHESS,
     &                  AVEISO,ANISO,ASYM,SPAR,APAR,CSTRA,SCTRA)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "codata.h"
      PARAMETER (D3 = 3.0D0, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/(4*XFAMU*XFAMU*PMASS*PMASS)
     &                   *XTHZ)
      LOGICAL TEST
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "dorps.h"
#include "cbisol.h"
#include "spnout.h"
#include "cbinum.h"
      DIMENSION SPNDSO(MXCOOR,MXCOOR), SPNPSO(MXCOOR,MXCOOR),
     &          SPNSD (MXCOOR,MXCOOR), SPNFC (MXCOOR,MXCOOR),
     &          SPSDFC(MXCOOR,MXCOOR)
      DIMENSION CHESS(3*NUCDEP,3*NUCDEP,6),
     &          AVEISO(NUCDEP*(NUCDEP + 1)/2,5),
     &          ANISO (NUCDEP*(NUCDEP + 1)/2),
     &          ASYM (NUCDEP*(NUCDEP + 1)/2),
     &          SPAR (NUCDEP*(NUCDEP + 1)/2),
     &          APAR (NUCDEP*(NUCDEP + 1)/2),
     &          COUPLM(3,3), CSTRA(*), SCTRA(*)

C
C Ordering of the properties in the isotropic matrix:
C 1: Total spin-spin coupling matrix
C 2: DSO-contribution
C 3: PSO-contribution
C 4: SD-contribution
C 5: FC-contribution
C 6: FC/SD-contribution
C
C Construction of total spin-spin coupling matrix
C
C
      CALL TITLER
     &    ('ABACUS - INDIRECT NUCLEAR SPIN-SPIN COUPLINGS','*',116)
C
C Transform to non-symmetry basis
C
      CALL DZERO(CHESS,6*9*NUCDEP*NUCDEP)
      CALL TRAHES(SPNDSO,MXCOOR,CHESS(1,1,2),CSTRA,SCTRA,3*NUCDEP,
     &            3*NUCDEP,2)
      CALL TRAHES(SPNPSO,MXCOOR,CHESS(1,1,3),CSTRA,SCTRA,3*NUCDEP,
     &            3*NUCDEP,2)
      CALL TRAHES(SPNSD, MXCOOR,CHESS(1,1,4),CSTRA,SCTRA,3*NUCDEP,
     &            3*NUCDEP,2)
      CALL TRAHES(SPNFC, MXCOOR,CHESS(1,1,5),CSTRA,SCTRA,3*NUCDEP,
     &            3*NUCDEP,2)
      CALL TRAHES(SPSDFC,MXCOOR,CHESS(1,1,6),CSTRA,SCTRA,3*NUCDEP,
     &            3*NUCDEP,2)
C
      DO 20 ICOL = 1, 3*NUCDEP
         DO 10 IROW = 1, 3*NUCDEP
            CHESS(IROW,ICOL,1) = CHESS(IROW,ICOL,2) + CHESS(IROW,ICOL,3)
     &                         + CHESS(IROW,ICOL,4) + CHESS(IROW,ICOL,5)
     &                         + CHESS(IROW,ICOL,6)
 10      CONTINUE
 20   CONTINUE
C
C Print individual contributions, if requested
C
      IF (ISPPRI .GE. 10) THEN
         CALL HEADER('DSO-part of coupling constants',1)
         CALL PRIHES (SPNDSO,'SPNSPN',CSTRA,SCTRA)
         CALL HEADER ('PSO-part of coupling constants',1)
         CALL PRIHES (SPNPSO,'SPNSPN',CSTRA,SCTRA)
         CALL HEADER ('SD-part of coupling constants',1)
         CALL PRIHES (SPNSD,'SPNSPN',CSTRA,SCTRA)
         CALL HEADER ('FC-part of coupling constants',1)
         CALL PRIHES (SPNFC,'SPNSPN',CSTRA,SCTRA)
         CALL HEADER ('SD-FC-crosscoupling part of coupling constants'
     &               ,1)
         CALL PRIHES (SPSDFC,'SPNSPN',CSTRA,SCTRA)
         CALL HEADER('Total spin-spin coupling matrix',1)
         NCOOR = 3*NUCDEP
         CALL PR2DER(CHESS(1,1,1),NCOOR,NCOOR,LUPRI)
      END IF
C
C     *** Final indirect spin-spin coupling analysis. ***
C
      IF (SPNISO) THEN
         CALL ISOSPN(CHESS,AVEISO,ANISO,ASYM,SPAR,APAR)
      ELSE
         CALL FNLSPN(CHESS,AVEISO,ANISO,ASYM,SPAR,APAR)
      END IF
C
C     *****************************************************
C     *** If vibrational averaging we need to write our ***
C     *** results to file.                              ***
C     *****************************************************
C
      IF (NPRPDR) THEN
         CALL WRAVFL(CHESS,3*NUCDEP,3*NUCDEP,6,'SPIN-SPIN',ISPPRI)
      END IF
C
      RETURN
      END
C
C
C  /* Deck isospn */
      SUBROUTINE ISOSPN(CHESS,AVEISO,ANISO,ASYM,SPAR,APAR)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "codata.h"
#include "cbisol.h"
#include "maxaqn.h"
      PARAMETER (D3 = 3.0D0, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/(4*XFAMU*XFAMU*PMASS*PMASS)
     &                   *XTHZ)
      CHARACTER LAPRPC*8, LANAM1*8, LANAM2*8
#include "dorps.h"
#include "nuclei.h"
#include "symmet.h"
#include "spnout.h"
#include "abainf.h"
      LOGICAL TEST
      DIMENSION CHESS(3*NUCDEP,3*NUCDEP,6),
     &                AVEISO(NUCDEP*(NUCDEP + 1)/2,5),
     &                ANISO (NUCDEP*(NUCDEP + 1)/2  ),
     &                ASYM  (NUCDEP*(NUCDEP + 1)/2  ),
     &                SPAR  (NUCDEP*(NUCDEP + 1)/2  ),
     &                APAR  (NUCDEP*(NUCDEP + 1)/2  ),
     &                COUPLM(3,3)

C
C     *** Calculate the various spin parameters. ***
C
      WRITE (LUPRI,'(A/)') ' Definitions from Smith, Palke, and'//
     &       ' Grieg, Concepts in Mag. Res. 4 (1992) 107'
      DO 30 ITYP = 1, 5
         IF (ITYP .EQ. 1) THEN
            CALL SPIANA(CHESS(1,1,ITYP),AVEISO(1,ITYP),ANISO,
     &                  ASYM,SPAR,APAR)
         ELSE
            CALL DZERO(AVEISO(1,ITYP),NUCDEP*(NUCDEP + 1)/2)
            IATIJ = 0
            DO 40 IATOM1 = 1, NUCDEP
               DO 50 IATOM2 = 1, IATOM1
                  IATIJ = IATIJ + 1
                  IF (IATOM1 .GT. IATOM2) THEN
                     CALL DZERO(COUPLM,9)
                     DO 60 ICOOR1 = 1, 3
                     DO 60 ICOOR2 = 1, 3
                        IADR1 = 3*(IATOM1 - 1) + ICOOR1
                        IADR2 = 3*(IATOM2 - 1) + ICOOR2
                        COUPLM(ICOOR1,ICOOR2) = CHESS(IADR1,IADR2,ITYP)
 60                  CONTINUE
                     AVEISO(IATIJ,ITYP) = (COUPLM(1,1) + COUPLM(2,2)
     &                                  +  COUPLM(3,3))/D3
                  ELSE
                     AVEISO(IATIJ,ITYP) = D0
                  END IF
 50            CONTINUE
 40         CONTINUE
         END IF
 30   CONTINUE
C
C Print section
C
      IF (ISPPRI .GT. 5) THEN
         DO 70 ITYP = 1, 5
            IF (ITYP .EQ. 1) THEN
               CALL HEADER('Total spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 2) THEN
               CALL HEADER('DSO contribution to spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 3) THEN
               CALL HEADER('PSO contribution to spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 4) THEN
               CALL HEADER('SD contribution to spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 5) THEN
               CALL HEADER('FC contribution to spin-spin couplings',-1)
            END IF
            CALL HEADER('Isotropic couplings',-1)
            CALL PRIDIS(NAMDEP,AVEISO(1,ITYP),NUCDEP)
            IF (ISPPRI .GT. 10 .AND. ITYP .EQ. 1) THEN
               CALL HEADER('Anisotropic couplings',-1)
               CALL PRIDIS(NAMDEP,ANISO(1),NUCDEP)
               IF (ITYP .GT. 2) THEN
                  CALL HEADER('Asymmetric part of couplings',-1)
                  CALL PRIDIS(NAMDEP,ASYM(1),NUCDEP)
                  CALL HEADER('S-parameter of couplings',-1)
                  CALL PRIDIS(NAMDEP,APAR(1),NUCDEP)
                  CALL HEADER('A-parameter of couplings',-1)
                  CALL PRIDIS(NAMDEP,SPAR(1),NUCDEP)
               END IF
            END IF
 70      CONTINUE
      END IF
C
      CALL AROUND('ABACUS - Final spin-spin couplings')
      IF ((.NOT. DOSD) .OR. (.NOT. DOFC) .OR. (.NOT. DOPSO)
     &                 .OR. (.NOT. DODSO)) THEN
         WRITE (LUPRI,'(A)')
     &         ' >>>WARNING: Total spin-spin couplings not correct '//
     &         'because some'
         WRITE (LUPRI,'(A)')
     &         '             contributions have not been calculated !'
      END IF
      IATOM1 = 0
      NUCINS = NUCIND
      IF (SOLVNT) NUCINS = NUCINS - 1
      DO 80 I1 = 1, NUCINS
         DO 82 ISYM1 = 0, MAXOPR
            IF (IAND(ISTBNU(I1),ISYM1) .EQ. 0) THEN
               IATOM1 = IATOM1 + 1
               IF (DOPERT(I1,2)) THEN
      IATOM2 = 0
      DO 85 I2 = 1, I1
         IF (I2 .EQ. I1) THEN
            MAXSYM = ISYM1 - 1
         ELSE
            MAXSYM = MAXOPR
         END IF
         DO 87 ISYM2 = 0, MAXSYM
            IF (IAND(ISTBNU(I2),ISYM2) .EQ. 0) THEN
               IATOM2 = IATOM2 + 1
               IF (.NOT.(NCSPNI(I1) .OR. NCSPNI(I2))) GOTO 87
               IF (DOPERT(I2,2)) THEN
                  IATIJ = IATOM1*(IATOM1 -1 )/2 + IATOM2
                  NZ1 = IZATOM(I1)
                  NZ2 = IZATOM(I2)
                  CALL TITLER('Indirect spin-spin coupling between '//
     &               NAMDEP(IATOM1)//' and '//
     &               NAMDEP(IATOM2)//':','=',-1)
                  GVAL1 = DISOTP(NZ1,ISOTPS(IATOM1),'GVAL')
                  IF (GVAL1 .NE. 0) THEN
                     GVAL2 = DISOTP(NZ2,ISOTPS(IATOM2),'GVAL')
                     IF (GVAL2 .NE. 0) THEN
                        ABUND1 = DISOTP(NZ1,ISOTPS(IATOM1),'ABUNDANCE')
                        ABUND2 = DISOTP(NZ2,ISOTPS(IATOM2),'ABUNDANCE')
                        NA1    = NINT(DISOTP(NZ1,ISOTPS(IATOM1),'A'))
                        NA2    = NINT(DISOTP(NZ2,ISOTPS(IATOM2),'A'))
                        FACTOR = AUTOHZ*GVAL1*GVAL2
                        WRITE(LUPRI,'(/2X,A,I5,A,F8.3,A,F10.6)')
     &                     'Mass number atom 1:',NA1,';   Abundance:',
     &                     ABUND1,' %;   g factor: ',GVAL1
                        WRITE(LUPRI,'(2X,A,I5,A,F8.3,A,F10.6)')
     &                     'Mass number atom 2:',NA2,';   Abundance:',
     &                     ABUND2,' %;   g factor: ',GVAL2
                        AVETOT = AVEISO(IATIJ,1)*FACTOR
                        ANIS   = ANISO(IATIJ)*FACTOR
                        ASYMM  = ASYM(IATIJ)
                        SPARA  = SPAR(IATIJ)*FACTOR
                        APARA  = APAR(IATIJ)*FACTOR
                        DSO    = AVEISO(IATIJ,2)*FACTOR
                        PSO    = AVEISO(IATIJ,3)*FACTOR
                        SD     = AVEISO(IATIJ,4)*FACTOR
                        FC     = AVEISO(IATIJ,5)*FACTOR
C MBH+JK
                        LANAM1 = NAMDEP(IATOM1)//'  '
                        LANAM2 = NAMDEP(IATOM2)//'  '
                        LAPRPC = 'SP-SPISO'
                        CALL WRIPRO(AVETOT/AUTOHZ,"SCF/DFT   ",410,
     *                  LANAM1,LANAM2,LAPRPC,LAPRPC,
     *                  0.0D0,0.0D0,0.0D0,1,0,0,0)
                        LAPRPC = 'SP-SPDSO'
                        CALL WRIPRO(DSO/AUTOHZ,"SCF/DFT   ",410,
     *                  LANAM1,LANAM2,LAPRPC,LAPRPC,
     *                  0.0D0,0.0D0,0.0D0,1,0,0,0)
                        LAPRPC = 'SP-SPPSO'
                        CALL WRIPRO(PSO/AUTOHZ,"SCF/DFT   ",410,
     *                  LANAM1,LANAM2,LAPRPC,LAPRPC,
     *                  0.0D0,0.0D0,0.0D0,1,0,0,0)
                        LAPRPC = 'SP-SPSD '
                        CALL WRIPRO(SD/AUTOHZ,"SCF/DFT   ",410,
     *                  LANAM1,LANAM2,LAPRPC,LAPRPC,
     *                  0.0D0,0.0D0,0.0D0,1,0,0,0)
                        LAPRPC = 'SP-SPFC '
                        CALL WRIPRO(FC/AUTOHZ,"SCF/DFT   ",410,
     *                  LANAM1,LANAM2,LAPRPC,LAPRPC,
     *                  0.0D0,0.0D0,0.0D0,1,0,0,0)
C                        WRITE(LUPRI,*)'MBH: AUTOHZ VALUE: ',AUTOHZ
C MBH+JK
                        WRITE(LUPRI,'(9(/,2X,A,F12.4,A),/,/)')
     &                       'Isotropic coupling        :',AVETOT,' Hz',
     &                       'Anisotropic coupling      :',ANIS  ,' Hz',
     &                       'Asymmetry                 :',ASYMM ,'   ',
     &                       'S parameter               :',SPARA ,' Hz',
     &                       'A parameter               :',APARA ,' Hz',
     &                       'Isotropic DSO contribution:',DSO   ,' Hz',
     &                       'Isotropic PSO contribution:',PSO   ,' Hz',
     &                       'Isotropic SD contribution :',SD    ,' Hz',
     &                       'Isotropic FC contribution :',FC    ,' Hz'
CSPAS: trying to print the whole tensor
                        DO ICOOR1 = 1, 3
                           DO ICOOR2 = 1, 3
                              IADR1 = 3*(IATOM1 - 1) + ICOOR1
                              IADR2 = 3*(IATOM2 - 1) + ICOOR2
                        COUPLM(ICOOR1,ICOOR2) = CHESS(IADR1,IADR2,1)
     &                                                *FACTOR
                           END DO
                        END DO
                    CALL TITLER('Total Coupling Tensor in Hz','-',-1)
                        WRITE(LUPRI,'(16X,3(A14,1X),/)') 'x','y','z'
                        WRITE(LUPRI,'(16X,A,2X,1P,3(G14.7,1X))')
     &                       'x',(COUPLM(1,J), J=1,3)
                        WRITE(LUPRI,'(16X,A,2X,1P,3(G14.7,1X))')
     &                       'y',(COUPLM(2,J), J=1,3)
                        WRITE(LUPRI,'(16X,A,2X,1P,3(G14.7,1X))')
     &                       'z',(COUPLM(3,J), J=1,3)
CKeinSPASmehr: trying to print the whole tensor
                     END IF
                  END IF
               END IF
            END IF
 87      CONTINUE
 85   CONTINUE
               END IF
            END IF
 82      CONTINUE
 80   CONTINUE
C
      RETURN
      END
C
C
C  /* Deck fnlspn */
      SUBROUTINE FNLSPN(CHESS,AVEISO,ANISO,ASYM,SPAR,APAR)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "codata.h"
#include "cbisol.h"
#include "maxaqn.h"
      PARAMETER (D3 = 3.0D0, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/(4*XFAMU*XFAMU*PMASS*PMASS)
     &                   *XTHZ)
#include "dorps.h"
#include "nuclei.h"
#include "symmet.h"
#include "spnout.h"
      LOGICAL TEST
      DIMENSION CHESS(3*NUCDEP,3*NUCDEP,6),
     &                AVEISO(NUCDEP*(NUCDEP + 1)/2,5),
     &                ANISO (NUCDEP*(NUCDEP + 1)/2  ),
     &                ASYM  (NUCDEP*(NUCDEP + 1)/2  ),
     &                SPAR  (NUCDEP*(NUCDEP + 1)/2  ),
     &                APAR  (NUCDEP*(NUCDEP + 1)/2  ),
     &                COUPLM(3,3)

C
C     *** Calculate the various spin parameters. ***
C
      WRITE (LUPRI,'(A/)') ' Definitions from Smith, Palke, and'//
     &       ' Grieg, Concepts in Mag. Res. 4 (1992) 107'
      DO 30 ITYP = 1, 5
         IF (ITYP .EQ. 1) THEN
            CALL SPIANA(CHESS(1,1,ITYP),AVEISO(1,ITYP),ANISO,
     &                  ASYM,SPAR,APAR)
         ELSE
            CALL DZERO(AVEISO(1,ITYP),NUCDEP*(NUCDEP + 1)/2)
            IATIJ = 0
            DO 40 IATOM1 = 1, NUCDEP
               DO 50 IATOM2 = 1, IATOM1
                  IATIJ = IATIJ + 1
                  IF (IATOM1 .GT. IATOM2) THEN
                     CALL DZERO(COUPLM,9)
                     DO 60 ICOOR1 = 1, 3
                     DO 60 ICOOR2 = 1, 3
                        IADR1 = 3*(IATOM1 - 1) + ICOOR1
                        IADR2 = 3*(IATOM2 - 1) + ICOOR2
                        COUPLM(ICOOR1,ICOOR2) = CHESS(IADR1,IADR2,ITYP)
 60                  CONTINUE
                     AVEISO(IATIJ,ITYP) = (COUPLM(1,1) + COUPLM(2,2)
     &                                  +  COUPLM(3,3))/D3
                  ELSE
                     AVEISO(IATIJ,ITYP) = D0
                  END IF
 50            CONTINUE
 40         CONTINUE
         END IF
 30   CONTINUE
C
C Print section
C
      IF (ISPPRI .GT. 5) THEN
         DO 70 ITYP = 1, 5
            IF (ITYP .EQ. 1) THEN
               CALL HEADER('Total spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 2) THEN
               CALL HEADER('DSO contribution to spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 3) THEN
               CALL HEADER('PSO contribution to spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 4) THEN
               CALL HEADER('SD contribution to spin-spin couplings',-1)
            ELSE IF (ITYP .EQ. 5) THEN
               CALL HEADER('FC contribution to spin-spin couplings',-1)
            END IF
            CALL HEADER('Isotropic couplings',-1)
            CALL PRIDIS(NAMDEP,AVEISO(1,ITYP),NUCDEP)
            IF (ISPPRI .GT. 10 .AND. ITYP .EQ. 1) THEN
               CALL HEADER('Anisotropic couplings',-1)
               CALL PRIDIS(NAMDEP,ANISO(1),NUCDEP)
               IF (ITYP .GT. 2) THEN
                  CALL HEADER('Asymmetric part of couplings',-1)
                  CALL PRIDIS(NAMDEP,ASYM(1),NUCDEP)
                  CALL HEADER('S-parameter of couplings',-1)
                  CALL PRIDIS(NAMDEP,APAR(1),NUCDEP)
                  CALL HEADER('A-parameter of couplings',-1)
                  CALL PRIDIS(NAMDEP,SPAR(1),NUCDEP)
               END IF
            END IF
 70      CONTINUE
      END IF
C
      CALL AROUND('ABACUS - Final spin-spin couplings')
      IF ((.NOT. DOSD) .OR. (.NOT. DOFC) .OR. (.NOT. DOPSO)
     &                 .OR. (.NOT. DODSO)) THEN
         WRITE (LUPRI,'(A)')
     &         ' >>>WARNING: Total spin-spin couplings not correct '//
     &         'because some'
         WRITE (LUPRI,'(A)')
     &         '             contributions have not been calculated !'
      END IF
      IATOM1 = 0
      NUCINS = NUCIND
      IF (SOLVNT) NUCINS = NUCINS - 1
      DO 80 I1 = 1, NUCINS
         DO 82 ISYM1 = 0, MAXOPR
            IF (IAND(ISTBNU(I1),ISYM1) .EQ. 0) THEN
               IATOM1 = IATOM1 + 1
               IF (DOPERT(I1,2)) THEN
      IATOM2 = 0
      DO 85 I2 = 1, I1
         IF (I2 .EQ. I1) THEN
            MAXSYM = ISYM1 - 1
         ELSE
            MAXSYM = MAXOPR
         END IF
         DO 87 ISYM2 = 0, MAXSYM
            IF (IAND(ISTBNU(I2),ISYM2) .EQ. 0) THEN
               IATOM2 = IATOM2 + 1
               IF (.NOT.(NCSPNI(I1) .OR. NCSPNI(I2))) GOTO 87
               IF (DOPERT(I2,2)) THEN
                  IATIJ = IATOM1*(IATOM1 -1 )/2 + IATOM2
                  NZ1 = IZATOM(I1)
                  NZ2 = IZATOM(I2)
                  CALL TITLER('Indirect spin-spin coupling between '//
     &               NAMDEP(IATOM1)//' and '//
     &               NAMDEP(IATOM2)//':','=',-1)
                  TEST = .FALSE.
                  DO 90 ISO1 = 1, 5
                     GVAL1 = DISOTP(NZ1,ISO1,'GVAL')
                     IF (GVAL1 .NE. 0) THEN
                        IF (IATOM1 .EQ. IATOM2) THEN
                           ISOMAX = ISO1
                        ELSE
                           ISOMAX = 5
                        END IF
                        DO 100 ISO2 = 1, ISOMAX
                           GVAL2 = DISOTP(NZ2,ISO2,'GVAL')
                           IF (GVAL2 .NE. 0) THEN
                              ABUND1 = DISOTP(NZ1,ISO1,'ABUNDANCE')
                              ABUND2 = DISOTP(NZ2,ISO2,'ABUNDANCE')
                              IF (     ((ABUND1 .GE. ABUND)
     &                            .AND. (ABUND2 .GE. ABUND))
     &                            .OR.  (.NOT. TEST)) THEN
                              TEST = .TRUE.
                              NA1    = NINT(DISOTP(NZ1,ISO1,'A'))
                              NA2    = NINT(DISOTP(NZ2,ISO2,'A'))
                              FACTOR = AUTOHZ*GVAL1*GVAL2
                              WRITE(LUPRI,'(/2X,A,I5,A,F8.3,A,F10.6)')
     &                           'Mass number atom 1:',NA1,
     &                           ';   Abundance:',ABUND1,
     &                           ' %;   g factor: ',GVAL1
                              WRITE(LUPRI,'(2X,A,I5,A,F8.3,A,F10.6)')
     &                           'Mass number atom 2:',NA2,
     &                           ';   Abundance:',ABUND2,
     &                           ' %;   g factor: ',GVAL2
                              AVETOT = AVEISO(IATIJ,1)*FACTOR
                              ANIS   = ANISO(IATIJ)*FACTOR
                              ASYMM  = ASYM(IATIJ)
                              SPARA  = SPAR(IATIJ)*FACTOR
                              APARA  = APAR(IATIJ)*FACTOR
                              DSO    = AVEISO(IATIJ,2)*FACTOR
                              PSO    = AVEISO(IATIJ,3)*FACTOR
                              SD     = AVEISO(IATIJ,4)*FACTOR
                              FC     = AVEISO(IATIJ,5)*FACTOR
                           WRITE(LUPRI,'(9(/,2X,A,F12.4,A),/,/)')
     &                       'Isotropic coupling        :',AVETOT,' Hz',
     &                       'Anisotropic coupling      :',ANIS  ,' Hz',
     &                       'Asymmetry                 :',ASYMM ,'   ',
     &                       'S parameter               :',SPARA ,' Hz',
     &                       'A parameter               :',APARA ,' Hz',
     &                       'Isotropic DSO contribution:',DSO   ,' Hz',
     &                       'Isotropic PSO contribution:',PSO   ,' Hz',
     &                       'Isotropic SD contribution :',SD    ,' Hz',
     &                       'Isotropic FC contribution :',FC    ,' Hz'
CSPAS: trying to print the whole tensor
                        DO ICOOR1 = 1, 3
                           DO ICOOR2 = 1, 3
                              IADR1 = 3*(IATOM1 - 1) + ICOOR1
                              IADR2 = 3*(IATOM2 - 1) + ICOOR2
                        COUPLM(ICOOR1,ICOOR2) = CHESS(IADR1,IADR2,1)
     &                                                *FACTOR
                           END DO
                        END DO
                    CALL TITLER('Total Coupling Tensor in Hz','-',-1)
                        WRITE(LUPRI,'(16X,3(A14,1X),/)') 'x','y','z'
                        WRITE(LUPRI,'(16X,A,2X,1P,3(G14.7,1X))')
     &                       'x',(COUPLM(1,J), J=1,3)
                        WRITE(LUPRI,'(16X,A,2X,1P,3(G14.7,1X))')
     &                       'y',(COUPLM(2,J), J=1,3)
                        WRITE(LUPRI,'(16X,A,2X,1P,3(G14.7,1X))')
     &                       'z',(COUPLM(3,J), J=1,3)
CKeinSPASmehr: trying to print the whole tensor
                           END IF
                           END IF
 100                    CONTINUE
                     END IF
 90               CONTINUE
               END IF
            END IF
 87      CONTINUE
 85   CONTINUE
               END IF
            END IF
 82      CONTINUE
 80   CONTINUE
C
      RETURN
      END
C  /* Deck spiana */
      SUBROUTINE SPIANA(PRIBAS,AVEISO,ANISO,ASYM,SPAR,APAR)
#include "implicit.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D3 = 3.0D0, D1 = 1.0D0)
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION PRIBAS(3*NUCDEP,3*NUCDEP),
     &          AVEISO(NUCDEP*(NUCDEP + 1)/2),
     &          ANISO(NUCDEP*(NUCDEP + 1)/2),
     &          ASYM(NUCDEP*(NUCDEP + 1)/2),
     &          SPAR(NUCDEP*(NUCDEP + 1)/2),
     &          APAR(NUCDEP*(NUCDEP + 1)/2), COUPLM(3,3),
     &          PVAL(6), AXES(3,3), WRK(3), IWRK(3), SYM(3,3),
     &          ANTI(3,3)
C
      CALL DZERO(AVEISO,NUCDEP*(NUCDEP + 1)/2)
      IATIJ = 0
      DO 10 IATOM1 = 1, NUCDEP
         DO 20 IATOM2 = 1, IATOM1
            IATIJ = IATIJ + 1
            IF (IATOM1 .GT. IATOM2) THEN
               DO 30 ICOOR1 = 1, 3
                  DO 30 ICOOR2 = 1, 3
                     IADR1 = 3*(IATOM1 - 1) + ICOOR1
                     IADR2 = 3*(IATOM2 - 1) + ICOOR2
                     COUPLM(ICOOR1,ICOOR2) = PRIBAS(IADR1,IADR2)
 30            CONTINUE
               AVETOT = (COUPLM(1,1) + COUPLM(2,2) + COUPLM(3,3))/D3
               DO 40 I = 1, 3
                  DO 50 J = 1, 3
                     ANTI(I,J) = (COUPLM(I,J) - COUPLM(J,I))/D2
                     SYM (I,J) = (COUPLM(I,J) + COUPLM(J,I))/D2
 50               CONTINUE
                  SYM(I,I) = SYM(I,I) - AVETOT
 40            CONTINUE
C
               CALL DUNIT(AXES,3)
               IJ = 1
               DO 60 I = 1, 3
                  DO 70 J = 1, I
                     PVAL(IJ) = SYM(I,J)
                     IJ = IJ + 1
 70               CONTINUE
 60            CONTINUE
               CALL JACO(PVAL,AXES,3,3,3,WRK,IWRK)
               PVAL(1) = PVAL(1)
               PVAL(2) = PVAL(3)
               PVAL(3) = PVAL(6)
C
               IZ = IDAMAX(3,PVAL,1)
               CALL DSWAP(1,PVAL(3),1,PVAL(IZ),1)
               CALL DSWAP(3,AXES(1,3),1,AXES(1,IZ),1)
               IY = IDAMIN(3,PVAL,1)
               CALL DSWAP(1,PVAL(2),1,PVAL(IY),1)
               CALL DSWAP(3,AXES(1,2),1,AXES(1,IY),1)
C
               ANISO(IATIJ) = PVAL(3) - (PVAL(2) + PVAL(1))/2
               IF (PVAL(3) .EQ. 0) THEN
                  ASYM(IATIJ) = D0
               ELSE
                  ASYM(IATIJ) = ABS(PVAL(2) - PVAL(1))/PVAL(3)
               END IF
               APAR(IATIJ) = SQRT(DDOT(9,ANTI,1,ANTI,1)/D2)
               SPAR(IATIJ) =
     &             ABS(ANISO(IATIJ))*SQRT(D1+(ASYM(IATIJ)**2)/D3)
               AVEISO(IATIJ) = AVETOT
            ELSE
               ANISO(IATIJ)  = D0
               ASYM(IATIJ)   = D0
               APAR(IATIJ)   = D0
               SPAR(IATIJ)   = D0
               AVEISO(IATIJ) = D0
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
C  /* Deck nucnqc */
      SUBROUTINE NUCNQC(GEOM,LUPRI,IPRINT)
C
C     Nuclear contribution to electric field gradients
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "gnrinf.h"
#include "nuclei.h"
      PARAMETER (D3 = 3.0D0)
      DIMENSION GEOM(3*NATOMS)
#include "symmet.h"
#include "qm3.h"
#include "nqcc.h"

C
      CALL DZERO(UCNNQC,9*MXCENT)
C
      KK = 1
      NATOM = 1
C
      DO 10 I = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               LL = 1
               DO 100 I2 = 1, NUCIND
                  IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
                  DO 110 ISYM = 0, MAXOPR
                     IF (IAND(ISYM,ISTBNU(I2)) .EQ. 0) THEN
                        IF (LL .NE. KK) THEN
                           XCOOR = GEOM(LL) - GEOM(KK)
                           YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                           ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                           R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                           R5 = R2*R2*SQRT(R2)
                           UCNNQC(1,1,NATOM) = UCNNQC(1,1,NATOM)
     &                          + CHARGE(I2)*(D3*XCOOR*XCOOR - R2)/R5
                           UCNNQC(2,2,NATOM) = UCNNQC(2,2,NATOM)
     &                          + CHARGE(I2)*(D3*YCOOR*YCOOR - R2)/R5
                           UCNNQC(3,3,NATOM) = UCNNQC(3,3,NATOM)
     &                          + CHARGE(I2)*(D3*ZCOOR*ZCOOR - R2)/R5
                           UCNNQC(1,2,NATOM) = UCNNQC(1,2,NATOM)
     &                          + CHARGE(I2)*D3*XCOOR*YCOOR/R5
                           UCNNQC(1,3,NATOM) = UCNNQC(1,3,NATOM)
     &                          + CHARGE(I2)*D3*XCOOR*ZCOOR/R5
                           UCNNQC(2,3,NATOM) = UCNNQC(2,3,NATOM)
     &                          + CHARGE(I2)*D3*YCOOR*ZCOOR/R5
                        END IF
                        LL = LL + 3
                     END IF
 110              CONTINUE
 100           CONTINUE
               UCNNQC(2,1,NATOM) = UCNNQC(1,2,NATOM)
               UCNNQC(3,1,NATOM) = UCNNQC(1,3,NATOM)
               UCNNQC(3,2,NATOM) = UCNNQC(2,3,NATOM)
               NATOM = NATOM + 1
               KK = KK + 3
            END IF
 20      CONTINUE
 10   CONTINUE

      IF ( (QM3) .AND. (.NOT. (SKIPNC)) ) THEN
        CALL QM3QCC1(LUPRI,IPRINT)
        IF ( .NOT. (LOSPC) ) CALL QM3QCC2(LUPRI,IPRINT)
      END IF

      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Nuclear contribution to electric field gradient',
     &        -1)
         CALL EFGPRI(UCNNQC,'                                 ')
      END IF
      RETURN
      END
C  /* Deck nucqdr */
      SUBROUTINE NUCQDR(GEOM,ORIGIN,LUPRI,IPRINT)
C
C     Nuclear contribution to quadrupole moments
C     K.Ruud, nov.-94
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
      PARAMETER (DP5 = 0.50D0, D3 = 3.0D0)
      DIMENSION ORIGIN(3), GEOM(3*NATOMS)
#include "symmet.h"
#include "quadru.h"

C
      CALL DZERO(QDRNUC,9)
C
      KK = 1
      NATOM = 1
      IX = IPTAX(1,1)
      IY = IPTAX(2,1)
      IZ = IPTAX(3,1)
C
      DO 10 I = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               XCOOR = GEOM(KK) - ORIGIN(1)
               YCOOR = GEOM(KK + 1) - ORIGIN(2)
               ZCOOR = GEOM(KK + 2) - ORIGIN(3)
               R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
               QDRNUC(IX,IX) = QDRNUC(IX,IX) + DP5*CHARGE(I)*
     &                        (D3*XCOOR*XCOOR - R2)
               QDRNUC(IY,IY) = QDRNUC(IY,IY) + DP5*CHARGE(I)*
     &                        (D3*YCOOR*YCOOR - R2)
               QDRNUC(IZ,IZ) = QDRNUC(IZ,IZ) + DP5*CHARGE(I)*
     &                        (D3*ZCOOR*ZCOOR - R2)
               QDRNUC(IX,IY) = QDRNUC(IX,IY) + DP5*CHARGE(I)*
     &                         D3*XCOOR*YCOOR
               QDRNUC(IX,IZ) = QDRNUC(IX,IZ) + DP5*CHARGE(I)*
     &                         D3*XCOOR*ZCOOR
               QDRNUC(IY,IZ) = QDRNUC(IY,IZ) + DP5*CHARGE(I)*
     &                         D3*YCOOR*ZCOOR
               KK = KK + 3
            END IF
 20      CONTINUE
 10   CONTINUE
      QDRNUC(IY,IX) = QDRNUC(IX,IY)
      QDRNUC(IZ,IX) = QDRNUC(IX,IZ)
      QDRNUC(IZ,IY) = QDRNUC(IY,IZ)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Nuclear contribution to quadrupole moment',-1)
         CALL POLPRI(QDRNUC,'   ',1)
      END IF
      RETURN
      END
C  /* Deck nucscm */
      SUBROUTINE NUCSCM(GEOM,ORIGIN,LUPRI,IPRINT)
C
C     Nuclear contribution to second moments
C     H.Solheim, May -03
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
      PARAMETER (DP5 = 0.50D0, D3 = 3.0D0)
      DIMENSION ORIGIN(3), GEOM(3*NATOMS)
#include "symmet.h"
#include "secmom.h"

C
      CALL DZERO(SCMNUC,9)
      KK = 1
      NATOM = 1
      IX = IPTAX(1,1)
      IY = IPTAX(2,1)
      IZ = IPTAX(3,1)
      DO I = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               XCOOR = GEOM(KK) - ORIGIN(1)
               YCOOR = GEOM(KK + 1) - ORIGIN(2)
               ZCOOR = GEOM(KK + 2) - ORIGIN(3)
               SCMNUC(IX,IX) = SCMNUC(IX,IX) + CHARGE(I)*XCOOR*XCOOR
               SCMNUC(IY,IY) = SCMNUC(IY,IY) + CHARGE(I)*YCOOR*YCOOR
               SCMNUC(IZ,IZ) = SCMNUC(IZ,IZ) + CHARGE(I)*ZCOOR*ZCOOR
               SCMNUC(IX,IY) = SCMNUC(IX,IY) + CHARGE(I)*XCOOR*YCOOR
               SCMNUC(IX,IZ) = SCMNUC(IX,IZ) + CHARGE(I)*XCOOR*ZCOOR
               SCMNUC(IY,IZ) = SCMNUC(IY,IZ) + CHARGE(I)*YCOOR*ZCOOR
               KK = KK + 3
            END IF
         END DO
      END DO
      SCMNUC(IY,IX) = SCMNUC(IX,IY)
      SCMNUC(IZ,IX) = SCMNUC(IX,IZ)
      SCMNUC(IZ,IY) = SCMNUC(IY,IZ)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Nuclear contribution to second moment',-1)
         CALL POLPRI(SCMNUC,'   ',1)
      END IF
      RETURN
      END
C  /* Deck nuc3rd */
      SUBROUTINE NUC3RD(GEOM,ORIGIN,LUPRI,IPRINT)
C
C     Nuclear contribution to third moments
C     H.Solheim, May -03
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
      PARAMETER (DP5 = 0.50D0, D3 = 3.0D0)
      DIMENSION ORIGIN(3), GEOM(3*NATOMS)
#include "symmet.h"
#include "thimom.h"

C
      CALL DZERO(THDNUC,27)
      KK = 1
      NATOM = 1
      IX = IPTAX(1,1)
      IY = IPTAX(2,1)
      IZ = IPTAX(3,1)
      DO I = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               XCOOR = GEOM(KK) - ORIGIN(1)
               YCOOR = GEOM(KK + 1) - ORIGIN(2)
               ZCOOR = GEOM(KK + 2) - ORIGIN(3)
               THDNUC(IX,IX,IX) = THDNUC(IX,IX,IX)
     &              + CHARGE(I)*XCOOR*XCOOR*XCOOR
               THDNUC(IX,IX,IY) = THDNUC(IX,IX,IY)
     &              + CHARGE(I)*XCOOR*XCOOR*YCOOR
               THDNUC(IX,IX,IZ) = THDNUC(IX,IX,IZ)
     &              + CHARGE(I)*XCOOR*XCOOR*ZCOOR
               THDNUC(IX,IY,IY) = THDNUC(IX,IY,IY)
     &              + CHARGE(I)*XCOOR*YCOOR*YCOOR
               THDNUC(IX,IY,IZ) = THDNUC(IX,IY,IZ)
     &              + CHARGE(I)*XCOOR*YCOOR*ZCOOR
               THDNUC(IX,IZ,IZ) = THDNUC(IX,IZ,IZ)
     &              + CHARGE(I)*XCOOR*ZCOOR*ZCOOR
               THDNUC(IY,IY,IY) = THDNUC(IY,IY,IY)
     &              + CHARGE(I)*YCOOR*YCOOR*YCOOR
               THDNUC(IY,IY,IZ) = THDNUC(IY,IY,IZ)
     &              + CHARGE(I)*YCOOR*YCOOR*ZCOOR
               THDNUC(IY,IZ,IZ) = THDNUC(IY,IZ,IZ)
     &              + CHARGE(I)*YCOOR*ZCOOR*ZCOOR
               THDNUC(IZ,IZ,IZ) = THDNUC(IZ,IZ,IZ)
     &              + CHARGE(I)*ZCOOR*ZCOOR*ZCOOR
               KK = KK + 3
            END IF
         END DO
      END DO
      DO I = 1, 3
         DO J = I, 3
            DO K = J, 3
               THDNUC(IPTAX(I,1),IPTAX(K,1),IPTAX(J,1)) =
     &              THDNUC(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
               THDNUC(IPTAX(J,1),IPTAX(I,1),IPTAX(K,1)) =
     &              THDNUC(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
               THDNUC(IPTAX(J,1),IPTAX(K,1),IPTAX(I,1)) =
     &              THDNUC(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
               THDNUC(IPTAX(K,1),IPTAX(I,1),IPTAX(J,1)) =
     &              THDNUC(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
               THDNUC(IPTAX(K,1),IPTAX(K,1),IPTAX(I,1)) =
     &              THDNUC(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
            END DO
         END DO
      END DO
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Nuclear contribution to third moment',-1)
         CALL PRIOCT(THDNUC)
      END IF
      RETURN
      END
C  /* Deck nucmgf */
      SUBROUTINE NUCMGF(GEOM,ORIGIN,LUPRI,IPRINT)
C
C     Calculate nuclear contribution to molecular rotational g-factors
C     K.Ruud, June-94
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
      DIMENSION ORIGIN(3), GEOM(3*(NATOMS+NFLOAT))
#include "symmet.h"
#include "molgfa.h"

C
      CALL DZERO(GFANUC,9)
      KK = 1
      NATOM = 1
      IX = IPTAX(1,2)
      IY = IPTAX(2,2)
      IZ = IPTAX(3,2)
      DO 10 I = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               XCOOR = GEOM(KK) - ORIGIN(1)
               YCOOR = GEOM(KK + 1) - ORIGIN(2)
               ZCOOR = GEOM(KK + 2) - ORIGIN(3)
               R2    = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
               GFANUC(IX,IX) = GFANUC(IX,IX) + CHARGE(I)*(YCOOR*YCOOR
     &                                       + ZCOOR*ZCOOR)
               GFANUC(IY,IY) = GFANUC(IY,IY) + CHARGE(I)*(XCOOR*XCOOR
     &                                       + ZCOOR*ZCOOR)
               GFANUC(IZ,IZ) = GFANUC(IZ,IZ) + CHARGE(I)*(XCOOR*XCOOR
     &                                       + YCOOR*YCOOR)
               GFANUC(IX,IY) = GFANUC(IX,IY) - CHARGE(I)*XCOOR*YCOOR
               GFANUC(IX,IZ) = GFANUC(IX,IZ) - CHARGE(I)*XCOOR*ZCOOR
               GFANUC(IY,IZ) = GFANUC(IY,IZ) - CHARGE(I)*YCOOR*ZCOOR
               KK = KK + 3
            END IF
 20      CONTINUE
 10   CONTINUE
      GFANUC(IY,IX) = GFANUC(IX,IY)
      GFANUC(IZ,IX) = GFANUC(IX,IZ)
      GFANUC(IZ,IY) = GFANUC(IY,IZ)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER(
     &   'Nuclear contribution to molecular rotational g-factor',-1)
         CALL OUTPUT(GFANUC,1,3,1,3,3,3,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck  nucspr */
      SUBROUTINE NUCSPR(GEOM,ORIGIN,IPRINT)
C
C     Calculate nuclear contribution to spin-rotation constants
C     K.Ruud, October-94
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
      PARAMETER (D1 = 1.0D0, D0 = 0.D0)
      DIMENSION GEOM(3*NATOMS), ORIGIN(3)
#include "abainf.h"
#include "symmet.h"
#include "spinro.h"

C
      CALL DZERO(SPRNUC,3*MXCOOR)
      KK = 1
      NATOM1 = 1
      DO 10 IATOM1 = 1, NUCIND
         DO 20 ISYM1 = 0, MAXOPR
            IF (IAND(ISYM1,ISTBNU(IATOM1)) .EQ. 0) THEN
               IF (CHARGE(IATOM1) .GT. D0) THEN
                  LL = 1
                  DO 30 IATOM2 = 1, NUCIND
                     DO 40 ISYM2 = 0, MAXOPR
                        IF (IAND(ISYM2,ISTBNU(IATOM2)) .EQ. 0) THEN
                           ZVAL2 = CHARGE(IATOM2)
                           IF (ZVAL2 .GT. D0) THEN
                              IF (KK .NE. LL) THEN
                                 XDIF = GEOM(KK) - GEOM(LL)
                                 YDIF = GEOM(KK + 1) - GEOM(LL + 1)
                                 ZDIF = GEOM(KK + 2) - GEOM(LL + 2)
                                 PFAC  = ZVAL2/(SQRT(XDIF*XDIF
     &                                + YDIF*YDIF + ZDIF*ZDIF))**3
                                 SPRNUC(1,1,NATOM1) = SPRNUC(1,1,NATOM1)
     &                                + (YDIF*YDIF + ZDIF*ZDIF)*PFAC
                                 SPRNUC(2,2,NATOM1) = SPRNUC(2,2,NATOM1)
     &                                + (XDIF*XDIF + ZDIF*ZDIF)*PFAC
                                 SPRNUC(3,3,NATOM1) = SPRNUC(3,3,NATOM1)
     &                                + (XDIF*XDIF + YDIF*YDIF)*PFAC
                                 SPRNUC(2,1,NATOM1) = SPRNUC(2,1,NATOM1)
     &                                - YDIF*XDIF*PFAC
                                 SPRNUC(3,1,NATOM1) = SPRNUC(3,1,NATOM1)
     &                                - ZDIF*XDIF*PFAC
                                 SPRNUC(1,2,NATOM1) = SPRNUC(1,2,NATOM1)
     &                                - XDIF*YDIF*PFAC
                                 SPRNUC(3,2,NATOM1) = SPRNUC(3,2,NATOM1)
     &                                - ZDIF*YDIF*PFAC
                                 SPRNUC(1,3,NATOM1) = SPRNUC(1,3,NATOM1)
     &                                - XDIF*ZDIF*PFAC
                                 SPRNUC(2,3,NATOM1) = SPRNUC(2,3,NATOM1)
     &                                - YDIF*ZDIF*PFAC
                              END IF
                              LL = LL + 3
                           END IF
                        END IF
 40                  CONTINUE
 30               CONTINUE
                  KK = KK + 3
               END IF
               NATOM1 = NATOM1 + 1
            END IF
 20      CONTINUE
 10   CONTINUE
C
      IF (IPRINT .GE. 2) THEN
        CALL HEADER('Nuclear contribution to spin-rotation constants',1)
        CALL OUTPUT(SPRNUC,1,3,1,3*NUCDEP,3,MXCOOR,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck mgfres */
      SUBROUTINE MGFRES(GEOM,AMASS,WORK,LWORK,KATOM,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "codata.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           THRESH = 1.0D-5)
      LOGICAL PLANAR, LINEAR, D12, D23
      DIMENSION ANGMOM(3), TINERT(9), OMEGA(3,3), WORK(LWORK),
     &          EIGVAL(3), EIGVEC(3,3),
     &          GTRAN1(3,3), GTRAN2(3,3), GTRAN3(3,3), GEOM(KATOM,3),
     &          AMASS(KATOM)
#include "abainf.h"
#include "nuclei.h"
#include "symmet.h"
#include "suscpt.h"
#include "molgfa.h"
#include "orgcom.h"
#include "dftcom.h"

C
      D12 = .FALSE.
      D23 = .FALSE.
C
      IF (IPRINT .GE. 0)
     &     CALL TITLER('ABACUS - MOLECULAR ROTATIONAL g-TENSOR','*',124)
      CALL DZERO(TOTMGF,9)
C
C     Add terms to diamagnetic molecular g factor
C
      IF (.NOT. NOLOND) THEN
         DO 10 I = 1, 3
         DO 10 J = 1, 3
            GFACDI(I,J) = GFACDI(I,J) + SUS2EL(I,J) + SUSFS(I,J)
     &                  + SUSFSY(I,J)
            IF (DFTRUN) THEN
               GFACDI(I,J) = GFACDI(I,J) + SUSDFT(I,J)
            END IF
 10      CONTINUE
      END IF
      DO 200 I = 1, 3
         ANGMOM(I) = D1
         DO 200 J = 1, 3
            ELMGF(I,J)  = D4*ELMGF(I,J)
            GFACDI(I,J) = D4*GFACDI(I,J)
            TOTMGF(I,J) = GFANUC(I,J) + ELMGF(I,J) + GFACDI(I,J)
 200  CONTINUE
C
      IF (IPRINT .GT. 2) THEN
         CALL HEADER(
     &   'Nuclear contribution * proton mass * moment of inertia',-1)
         CALL POLPRI(GFANUC,'   ',-2)
         AVENUC = (GFANUC(1,1) + GFANUC(2,2) + GFANUC(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F16.6)') ' Average value:',AVENUC
C
         CALL HEADER(
     &   'Diamagnetic contribution * proton mass * moment of inertia',
     &   -1)
         CALL POLPRI(GFACDI,'   ',-2)
         AVEDIA = (GFACDI(1,1) + GFACDI(2,2) + GFACDI(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F16.6)') ' Average value:',AVEDIA
C
         CALL HEADER(
     &   'Electronic contribution * proton mass * moment of inertia',-1)
         CALL POLPRI(ELMGF,'   ',-2)
         AVEEL = (ELMGF(1,1) + ELMGF(2,2) + ELMGF(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F16.6)') ' Average value:',AVEEL
C
         CALL HEADER('Total molecular rotational g-tensor'//
     &      ' * proton mass * moment of inertia',-1)
         CALL POLPRI(TOTMGF,'   ',2)
         AVETOT = (TOTMGF(1,1) + TOTMGF(2,2) + TOTMGF(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F16.6)') ' Average value:',AVETOT
      END IF
C
C     Determine principal moments of inertia
C
      JATOM = 1
      DO 100 IATOM = 1, NUCIND
         NATTYP = IZATOM(IATOM)
         NUMIS  = ISOTOP(IATOM)
         DO 110 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
               DO 120 ICOOR = 1, 3
                  GEOM(JATOM,ICOOR) = PT(IAND(ISYMAX(ICOOR,1),ISYMOP))
     &                               *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
 120           CONTINUE
               JATOM = JATOM + 1
            END IF
 110     CONTINUE
 100  CONTINUE
      CALL WLKDIN(GEOM,AMASS,KATOM,ANGMOM,TINERT,OMEGA,EIGVAL,EIGVEC,
     &            .TRUE.,PLANAR,LINEAR)
C
C     Transform molecular rotational g-tensor to principal axis system,
C     and multiply by inverse moment of inertia and proton mass
C
      IF (IPRINT .GE. 0) THEN
         IF (NOLOND) THEN
            WRITE (LUPRI,'(/A)') ' Calculation without London orbitals.'
         ELSE
            WRITE (LUPRI,'(/A)')' London orbitals used.'
         END IF
      END IF
C
C     All components of the moment of inertia not defined for linear molecules
C
      IF (LINEAR) THEN
         IF (ABS(TOTMGF(1,1)) .GT. THRESH) THEN
            AVETOT = TOTMGF(1,1)*PMASS*EIGVAL(3)
            AVEDIA = GFACDI(1,1)*PMASS*EIGVAL(3)
            AVEEL  = ELMGF (1,1)*PMASS*EIGVAL(3)
            AVENUC = GFANUC(1,1)*PMASS*EIGVAL(3)
         ELSE
            AVETOT = TOTMGF(2,2)*PMASS*EIGVAL(3)
            AVEDIA = GFACDI(2,2)*PMASS*EIGVAL(3)
            AVEEL  = ELMGF (2,2)*PMASS*EIGVAL(3)
            AVENUC = GFANUC(2,2)*PMASS*EIGVAL(3)
         END IF
         CALL DZERO(GTRAN,9)
         GTRAN(1,1) = TOTMGF(1,1)*PMASS*EIGVAL(3)
         GTRAN(2,2) = TOTMGF(2,2)*PMASS*EIGVAL(3)
         GTRAN(3,3) = TOTMGF(3,3)*PMASS*EIGVAL(3)
         IF (IPRDEF .GE. 0) THEN
            WRITE (LUPRI,'(A,F16.6)') ' Moment of inertia :',
     &           D1/EIGVAL(3)
            CALL HEADER('Molecular rotational g-factor',1)
            WRITE(LUPRI,'(4(A,F16.8/))')
     &           ' Isotropic g-factor          :',AVETOT,
     &           ' Nuclear part of g-factor    :',AVENUC,
     &           ' Diamagnetic part of g-factor:',AVEDIA,
     &           ' Electronic part of g-factor :',AVEEL
            WRITE (LUPRI,'(A)') ' The molecule is linear'
         END IF
      ELSE
         IF (IPRDEF .GE. 0) THEN
            CALL HEADER
     &    ('Principal moments of inertia (a.u.) and principal axes',0)
            WRITE (LUPRI,'(3X,A,F12.6,6X,3F12.6)')
     &           'IA',D1/EIGVAL(1), (EIGVEC(I,1),I = 1,3),
     &           'IB',D1/EIGVAL(2), (EIGVEC(I,2),I = 1,3),
     &           'IC',D1/EIGVAL(3), (EIGVEC(I,3),I = 1,3)
         END IF
         CALL DZERO(GTRAN,9)
         CALL DZERO(GTRAN1,9)
         CALL DZERO(GTRAN2,9)
         CALL DZERO(GTRAN3,9)
         DO 50 I = 1, 3
            DO 50 J = 1, 3
               DO 60 IX = 1, 3
                  DO 60 IY = 1, 3
                     JX = IPTAX(IX,2)
                     JY = IPTAX(IY,2)
                     CONST = EIGVAL(I)*PMASS
                     GTRAN(I,J)  = GTRAN(I,J) + EIGVEC(IX,I)
     &                            *TOTMGF(JX,JY)*EIGVEC(IY,J)*CONST
                     GTRAN1(I,J) = GTRAN1(I,J) + EIGVEC(IX,I)
     &                            *GFANUC(JX,JY)*EIGVEC(IY,J)*CONST
                     GTRAN2(I,J) = GTRAN2(I,J) + EIGVEC(IX,I)
     &                            *ELMGF(JX,JY)*EIGVEC(IY,J)*CONST
                     GTRAN3(I,J) = GTRAN3(I,J) + EIGVEC(IX,I)
     &                            *GFACDI(JX,JY)*EIGVEC(IY,J)*CONST
 60            CONTINUE
 50      CONTINUE
         DIFF1 = ABS(GTRAN(1,1) - GTRAN(2,2))
         DIFF2 = ABS(GTRAN(2,2) - GTRAN(3,3))
         IF (DIFF1 .LT. THRESH) D12 = .TRUE.
         IF (DIFF2 .LT. THRESH) D23 = .TRUE.
         AVETOT = (GTRAN(1,1)  + GTRAN(2,2)  + GTRAN(3,3) )/D3
         AVENUC = (GTRAN1(1,1) + GTRAN1(2,2) + GTRAN1(3,3))/D3
         AVEEL  = (GTRAN2(1,1) + GTRAN2(2,2) + GTRAN2(3,3))/D3
         AVEDIA = (GTRAN3(1,1) + GTRAN3(2,2) + GTRAN3(3,3))/D3
         IF (IPRINT .GE. 0) THEN
         IF (D12 .AND. D23) THEN
            WRITE (LUPRI,'(4(/A,F16.8))')
     &           ' Molecular rot. g-factor :  ',AVETOT,
     &           ' Nuclear contribution    :  ',AVENUC,
     &           ' Diamagnetic contribution:  ',AVEDIA,
     &           ' Electronic contribution :  ',AVEEL
            WRITE (LUPRI,'(/A)')
     &           ' Molecular rotational g-tensor is spherical.'
         ELSE IF (D12 .OR. D23) THEN
            IF (D12) THEN
               PAR = GTRAN(3,3)
               PER = (GTRAN(1,1) + GTRAN(2,2))/D2
            ELSE
               PAR = GTRAN(1,1)
               PER = (GTRAN(2,2) + GTRAN(3,3))/D2
            END IF
            ANI = PAR - PER
            CALL HEADER('Molecular rotational g-factor',1)
            WRITE (LUPRI,'(6(A,F16.8,/))')
     &           ' Isotropic g-tensor       : ',AVETOT,
     &           ' Nuclear contribution     : ',AVENUC,
     &           ' Electronic contribution  : ',AVEEL,
     &           ' Parallel component:        ',PAR,
     &           ' Perpendicular component:   ',PER,
     &           ' Anisotropy:                ',ANI
            WRITE (LUPRI,'(A)')
     &           ' Molecular rotational g-tensor is cylindrical.'
            CALL HEADER('Molecular rotational g-tensor in principal'//
     &                  ' axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN(I,1),I = 1,3),
     &         'B', (GTRAN(I,2),I = 1,3),
     &         'C', (GTRAN(I,3),I = 1,3)
            CALL HEADER('Nuclear contribution '//
     &                  'in principal axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN1(I,1),I = 1,3),
     &         'B', (GTRAN1(I,2),I = 1,3),
     &         'C', (GTRAN1(I,3),I = 1,3)
            CALL HEADER('Electronic contribution '//
     &                  'in principal axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN2(I,1),I = 1,3),
     &         'B', (GTRAN2(I,2),I = 1,3),
     &         'C', (GTRAN2(I,3),I = 1,3)
            CALL HEADER('Diamagnetic contribution '//
     &                  'in principal axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN3(I,1),I = 1,3),
     &         'B', (GTRAN3(I,2),I = 1,3),
     &         'C', (GTRAN3(I,3),I = 1,3)
         ELSE
            ANI1 = GTRAN(1,1) - (GTRAN(2,2) + GTRAN(3,3))/D2
            ANI2 = GTRAN(2,2) - (GTRAN(1,1) + GTRAN(3,3))/D2
            CALL HEADER('Molecular rotational g-tensor',1)
            WRITE (LUPRI,'(6(1X,A,F16.8,/))')
     &           ' Isotropic g-tensor       : ',AVETOT,
     &           ' Nuclear contribution     : ',AVENUC,
     &           ' Diamagnetic contribution : ',AVEDIA,
     &           ' Electronic contribution  : ',AVEEL,
     &           ' 1st anisotropy:            ',ANI1,
     &           ' 2nd anisotropy:            ',ANI2
            CALL HEADER('Molecular rotational g-tensor in principal'//
     &                  ' axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN(I,1),I = 1,3),
     &         'B', (GTRAN(I,2),I = 1,3),
     &         'C', (GTRAN(I,3),I = 1,3)
            CALL HEADER('Nuclear contribution '//
     &                  'in principal axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN1(I,1),I = 1,3),
     &         'B', (GTRAN1(I,2),I = 1,3),
     &         'C', (GTRAN1(I,3),I = 1,3)
            CALL HEADER('Electronic contribution '//
     &                  'in principal axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN2(I,1),I = 1,3),
     &         'B', (GTRAN2(I,2),I = 1,3),
     &         'C', (GTRAN2(I,3),I = 1,3)
            CALL HEADER('Diamagnetic contribution '//
     &                  'in principal axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &         'A', (GTRAN3(I,1),I = 1,3),
     &         'B', (GTRAN3(I,2),I = 1,3),
     &         'C', (GTRAN3(I,3),I = 1,3)
         END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck sprres */
      SUBROUTINE SPRRES(GEOM,AMASS,GVAL,DIAMAT,PARMAT,WORK,LWORK,
     &                  KATOM,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "codata.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           THRESH = 1.0D-5, DP5 = 1/D2)
      PARAMETER (FAC = ALPHA2/(PMASS*XFAMU*XFAMU)*XTHZ*1.0D-3)
      LOGICAL PLANAR, LINEAR, D12, D23
      DIMENSION ANGMOM(3), TINERT(9), OMEGA(3,3), WORK(LWORK),
     &          EIGVAL(3), EIGVEC(3,3),
     &          GTRANN(3,3), GTRAND(3,3), GEOM(KATOM,3), AMASS(KATOM),
     &          GVAL(KATOM), GTRANP(3,3), PARMAT(3,3,KATOM),
     &          DIAMAT(3,3,KATOM)
#include "abainf.h"
#include "nuclei.h"
#include "symmet.h"
#include "sigma.h"
#include "spinro.h"
#include "orgcom.h"

C
      D12 = .FALSE.
      D23 = .FALSE.
C
      IF (IPRINT .GE. 0)
     &     CALL TITLER('ABACUS - SPIN-ROTATION CONSTANTS','*',104)
      CALL DZERO(TOTSPR,3*MXCOOR)
      CALL DZERO(GTRANT,9*MXCENT)
C
C     We need to transform the dia- and paramagnetic part
C     of the spin-rotation constant
C
      CALL DAXPY(9*NUCDEP,D1,SIGMAS,1,ELSPRD,1)
      KCSTRA = 1
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('SPRRES','TRADIP',KLAST,LWORK)
      IF (NOLOND) THEN
         CALL DCOPY(9*NUCDEP,SPRDNL,1,DIAMAT,1)
         CALL DSCAL(9*NUCDEP,DP5,DIAMAT,1)
      ELSE
         CALL TRADIP(ELSPRD,DIAMAT,WORK(KCSTRA),WORK(KSCTRA),
     &               3*NUCDEP,2,2)
      END IF
      CALL TRADIP(ELSPRP,PARMAT,WORK(KCSTRA),WORK(KSCTRA),3*NUCDEP,2,2)
C
C     We add up the total spin-rotation constants. Also add nuclear
C     and diamagnetic term into one total diamagnetic term
C
      DO 99 IATOM = 1, NUCDEP
         DO 200 I = 1, 3
            ANGMOM(I) = D1
            DO 200 J = 1, 3
               TOTSPR(I,J,IATOM) = DIAMAT(I,J,IATOM) + PARMAT(I,J,IATOM)
     &                           + SPRNUC(I,J,IATOM)/D2
 200        CONTINUE
 99   CONTINUE
C
C     Determine principal moments of inertia
C
      JATOM = 1
      DO 100 IATOM = 1, NUCIND
         NATTYP = IZATOM(IATOM)
         NUMIS  = ISOTOP(IATOM)
         DO 110 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
               GVAL(JATOM)  = DISOTP(NATTYP,NUMIS,'GVAL')
               DO 120 ICOOR = 1, 3
                  GEOM(JATOM,ICOOR) = PT(IAND(ISYMAX(ICOOR,1),ISYMOP))
     &                               *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
 120           CONTINUE
               JATOM = JATOM + 1
            END IF
 110     CONTINUE
 100  CONTINUE
      CALL WLKDIN(GEOM,AMASS,KATOM,ANGMOM,TINERT,OMEGA,EIGVAL,EIGVEC,
     &            .TRUE.,PLANAR,LINEAR)
C
C     Transform each spin-rotation constant to principal axis system,
C     and multiply by inverse moment of inertia and various constants
C
C     All components of the moment of inertia not defined for linear molecules
C
      IF (LINEAR) THEN
         IF (IPRINT .GE. 0) THEN
            WRITE (LUPRI,'(A)') ' The molecule is linear'
            WRITE (LUPRI,'(A,F16.6)')' Moment of inertia :',
     &                                  D1/EIGVAL(3)
         END IF
         DO 98 IATOM = 1, NUCDEP
            IF (GVAL(IATOM) .NE. 0.0D0) THEN
               IF (IPRINT .GE. 0) THEN
                  CALL HEADER('Spin-rotation constants (kHz) for   '//
     &                        NAMDEP(IATOM),-1)
                  WRITE (LUPRI,'(2X,A18,F6.2,/)')
     &                 'Nuclear g-value  :',GVAL(IATOM)
               END IF
               FACTOT = EIGVAL(3)*GVAL(IATOM)*FAC
               IF (ABS(TOTSPR(1,1,IATOM)) .GT. THRESH) THEN
                  AVETOT = TOTSPR(1,1,IATOM)*FACTOT
                  AVENUC = SPRNUC(1,1,IATOM)*FACTOT/D2
                  AVEPAR = PARMAT(1,1,IATOM)*FACTOT
                  AVEDIA = DIAMAT(1,1,IATOM)*FACTOT
               ELSE
                  AVETOT = TOTSPR(2,2,IATOM)*FACTOT
                  AVENUC = SPRNUC(2,2,IATOM)*FACTOT/D2
                  AVEPAR = PARMAT(2,2,IATOM)*FACTOT
                  AVEDIA = DIAMAT(2,2,IATOM)*FACTOT
               END IF
               CALL DZERO(GTRANT(1,1,IATOM),9)
               GTRANT(1,1,IATOM) = AVETOT
               GTRANT(2,2,IATOM) = AVETOT
               GTRANT(3,3,IATOM) = AVETOT
               IF (IPRINT .GE. 0) WRITE(LUPRI,'(4(1X,A,F12.4,/))')
     &              ' Total spin-rotation constant       :',AVETOT,
     &              ' Nuclear part of spin-rotation      :',AVENUC,
     &              ' Diamagnetic part of spin-rotation  :',AVEDIA,
     &              ' Paramagnetic part of spin-rotation :',AVEPAR
            END IF
 98      CONTINUE
      ELSE
         IF (IPRINT .GE. 0) THEN
            CALL HEADER
     &    ('Principal moments of inertia (a.u.) and principal axes',0)
            WRITE (LUPRI,'(3X,A,F12.6,6X,3F12.6)')
     &           'IA',D1/EIGVAL(1), (EIGVEC(I,1),I = 1,3),
     &           'IB',D1/EIGVAL(2), (EIGVEC(I,2),I = 1,3),
     &           'IC',D1/EIGVAL(3), (EIGVEC(I,3),I = 1,3)
         END IF
         DO 45 IATOM = 1, NUCDEP
            IF (GVAL(IATOM) .NE. 0.0D0) THEN
               FACTOT = GVAL(IATOM)*FAC
               CALL DZERO(GTRANN,9)
               CALL DZERO(GTRAND,9)
               CALL DZERO(GTRANP,9)
               DO 50 I = 1, 3
               DO 50 J = 1, 3
                  DO 60 IX = 1, 3
                     DO 60 IY = 1, 3
                        JX = IPTAX(IX,2)
                        JY = IPTAX(IY,2)
                        CONST = EIGVAL(I)*FACTOT
                        GTRANT(I,J,IATOM)=GTRANT(I,J,IATOM)+EIGVEC(IX,I)
     &                       *TOTSPR(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
                        GTRANN(I,J) = GTRANN(I,J) + EIGVEC(IX,I)
     &                       *SPRNUC(JX,JY,IATOM)*EIGVEC(IY,J)*CONST/D2
                        GTRAND(I,J) = GTRAND(I,J) + EIGVEC(IX,I)
     &                       *DIAMAT(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
                        GTRANP(I,J) = GTRANP(I,J) + EIGVEC(IX,I)
     &                       *PARMAT(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
 60               CONTINUE
 50            CONTINUE
               IF (IPRINT .GE. 0) THEN
                  CALL HEADER('Spin rotation constants (kHz) for '//
     &                      NAMDEP(IATOM),-1)
                  WRITE (LUPRI,'(2X,A18,F6.2,/)')
     &                 'Nuclear g-value  :',GVAL(IATOM)
                  WRITE (LUPRI,'(2X,A,/)') 'Components in the '//
     &              'principal axis system'
                  IF (IPRINT .GT. 2) THEN
                     CALL HEADER('Nuclear contribution',-1)
                     WRITE (LUPRI,'(3X,A,6X,3F16.6)')
     &                    'A', (GTRANN(I,1),I = 1,3),
     &                    'B', (GTRANN(I,2),I = 1,3),
     &                    'C', (GTRANN(I,3),I = 1,3)
                     CALL HEADER('Diamagnetic contribution',-1)
                     WRITE (LUPRI,'(3X,A,6X,3F16.6)')
     &                    'A', (GTRAND(I,1),I = 1,3),
     &                    'B', (GTRAND(I,2),I = 1,3),
     &                    'C', (GTRAND(I,3),I = 1,3)
                     CALL HEADER('Paramagnetic contribution',-1)
                     WRITE (LUPRI,'(3X,A,6X,3F16.6)')
     &                    'A', (GTRANP(I,1),I = 1,3),
     &                    'B', (GTRANP(I,2),I = 1,3),
     &                    'C', (GTRANP(I,3),I = 1,3)
                     CALL HEADER('Total spin rotation constants',-1)
                  END IF
                  WRITE (LUPRI,'(3X,A,6X,3F16.6)')
     &                 'A', (GTRANT(I,1,IATOM),I = 1,3),
     &                 'B', (GTRANT(I,2,IATOM),I = 1,3),
     &                 'C', (GTRANT(I,3,IATOM),I = 1,3)
               END IF
            END IF
 45      CONTINUE
      END IF
      RETURN
      END
C  /* Deck nqcel */
      SUBROUTINE NQCEL(TEMP,EXPVAL,NCOMP,IPRINT)
C
C     Collect and arrange electronic contribution to electric field gradients
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nqcc.h"
#include "nuclei.h"
      DIMENSION TEMP(3,3,NUCDEP), EXPVAL(NCOMP)
#include "symmet.h"

C
      ITYP = 0
C
      CALL DZERO(TEMP,9*NUCDEP)
C
      DO 100 IATOM = 1, NUCIND
         DO 200 ICOOR1 = 1, 3
            DO 210 ICOOR2 = ICOOR1, 3
               ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
               DO 300 IREPC = 0, MAXREP
                  IF (IAND(ISTBNU(IATOM),
     &                 IEOR(IREPC,ISYMIJ)).EQ.0) THEN
                     ITYP = ITYP + 1
                     TEMP(ICOOR1,ICOOR2,IATOM) =
     &                    TEMP(ICOOR1,ICOOR2,IATOM) + EXPVAL(ITYP)
     &                                               /NUCDEG(IATOM)
                  END IF
 300           CONTINUE
 210        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      NATOM = 1
C
      DO 110 IATOM = 1, NUCIND
         DO 310 IREPC = 0, MAXREP
            IF (IAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
               DO 220 I = 1, 3
                  DO 230 J = I, 3
                     ELNQC(I,J,NATOM) = TEMP(I,J,IATOM)*
     &                    PT(IAND(IEOR(ISYMAX(I,1),
     &                    ISYMAX(J,1)),IREPC))
                     IF (I .NE. J) ELNQC(J,I,NATOM) = ELNQC(I,J,NATOM)
 230              CONTINUE
 220           CONTINUE
               NATOM = NATOM + 1
            END IF
 310     CONTINUE
 110  CONTINUE
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Electronic contribution to electric field '//
     &               'gradient',-1)
         CALL EFGPRI(ELNQC,'                                 ')
      END IF
      RETURN
      END
C  /* Deck diaspr */
      SUBROUTINE DIASPR(EXPVAL,ELFLD,CSTRA,SCTRA,IPRINT)
C
C     Routine for calculating the diamagnetic contribution to the
C     spin-rotation constant, K.Ruud - Oct.-95
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
      DIMENSION EXPVAL(3*NUCDEP), ELFLD(3,NUCDEP),
     &          CSTRA(3*NUCDEP,3*NUCDEP), SCTRA(3*NUCDEP,3*NUCDEP)
#include "abainf.h"
#include "symmet.h"
#include "spinro.h"
#include "orgcom.h"
#include "helfey.h"

C
      CALL DZERO(SPRDNL,3*MXCOOR)
      CALL DZERO(HFGEL,MXCOOR)
C
C     We start by sorting the electric field results
C
      CALL TRACOR(CSTRA,SCTRA,1,3*NUCDEP,0)
      CALL DGEMM('N','T',1,3*NUCDEP,3*NUCDEP,1.D0,
     &           EXPVAL,1,
     &           SCTRA(1,1),3*NUCDEP,0.D0,
     &           ELFLD,1)
C
C     Then we construct the diamagnetic contribution
C
      NATOM = 1
      DO 10 IATOM = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               COORX = PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
     &               - CMXYZ(1)
               COORY = PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
     &               - CMXYZ(2)
               COORZ = PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
     &               - CMXYZ(3)
               SPRDNL(1,1,NATOM) = (ELFLD(2,NATOM)*COORY
     &                           +  ELFLD(3,NATOM)*COORZ)
               SPRDNL(2,2,NATOM) = (ELFLD(1,NATOM)*COORX
     &                           +  ELFLD(3,NATOM)*COORZ)
               SPRDNL(3,3,NATOM) = (ELFLD(1,NATOM)*COORX
     &                           +  ELFLD(2,NATOM)*COORY)
               SPRDNL(1,2,NATOM) =  ELFLD(2,NATOM)*COORX
               SPRDNL(1,3,NATOM) =  ELFLD(3,NATOM)*COORX
               SPRDNL(2,1,NATOM) =  ELFLD(1,NATOM)*COORY
               SPRDNL(2,3,NATOM) =  ELFLD(3,NATOM)*COORY
               SPRDNL(3,1,NATOM) =  ELFLD(1,NATOM)*COORZ
               SPRDNL(3,2,NATOM) =  ELFLD(2,NATOM)*COORZ
               NATOM = NATOM + 1
            END IF
 20      CONTINUE
 10   CONTINUE
C
      IF (IPRINT .GE. 2) THEN
        CALL HEADER('Diamagnetic contribution to '//
     &              'spin-rotation constants',1)
        CALL OUTPUT(SPRDNL,1,3,1,3*NUCDEP,3,MXCOOR,1,LUPRI)
      END IF
      RETURN
 1000 FORMAT (1X,A6,F17.10,2F24.10)
      END
C    
C  /* Deck selinf */    
      SUBROUTINE SELINF(MPATOM,JPATOM,ATMALL)
C        
C 
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "cbiher.h"
C    
      DIMENSION JPATOM(MXCENT)
      LOGICAL ATMALL
C    
C    
      ATMALL = ALLATM
C    
      MPATOM = NPATOM
      DO I = 1,MXCENT
         JPATOM(I) = IPATOM(I)
      END DO
C           
      RETURN
      END  
