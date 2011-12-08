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
C  /* Deck troinp */
      SUBROUTINE TROINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 9)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <cbitro.h>
      LOGICAL NEWDEF
C
      DATA TABLE /'.SKIP  ', '.PRINT ', '.STOP  ', '.THRESH', '.COMPAR',
     *            'XXXXXXX', 'xGD    ', 'xRD    ', 'xNOROT '/
C    .GD, .RD, and .NOROT disabled because not working currently.
C    Nov 2004/hjaaj
C
      NEWDEF = (WORD .EQ. '*TROINV')
      THRSHD = THRESH
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
                     GO TO (1,2,3,4,5,6,7,8,9), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in TROINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in TROINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3             CONTINUE
                  CUT    = .TRUE.
               GO TO 100
    4             READ (LUCMD,*) THRESH
                  IF (THRESH .EQ. THRSHD) ICHANG = ICHANG - 1
               GO TO 100
    5             COMPAR = .TRUE.
               GO TO 100
    6          GO TO 100
    7             GDTRO = .TRUE.
                  GDALL = .TRUE.
               GO TO 100
    8             RDTRO = .TRUE.
C                 -- currently RDTRO only possible if NOROT
                  NOROT1 = .TRUE.
               GO TO 100
    9             NOROT1 = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in TROINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in TROINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for TROINV:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' TROINV skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in TROINV:',IPRINT
            END IF
            IF (THRESH .NE. THRSHD) THEN
               WRITE(LUPRI,'(A,D12.6)')
     *            ' Threshold for linear dependence: ',THRESH
            END IF
            IF (NOROT1) WRITE (LUPRI,'(/,A,/)')
     *         ' Rotational symmetry not used.'
            IF (COMPAR) WRITE (LUPRI,'(/,2A)') ' Translational and',
     *         ' rotational symmetry will be used for comparison.'
            IF (GDTRO) WRITE(LUPRI,'(/,2A)')' Differentiated gradients',
     *         ' will be calculated using symmetry.'
            IF (RDTRO) WRITE (LUPRI,'(/,2A)') ' Solution vectors',
     *         ' will be calculated using symmetry.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after TROINV.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck troini */
      SUBROUTINE TROINI
C
C     Initialize /CBITRO/
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <nuclei.h>
! nuclei.h needed for NUCDEP
#include <abainf.h>
#include <cbitro.h>
#include <infinp.h>
C
      IPRINT = IPRDEF
      THRESH = 0.1D00
      COMPAR = VCD .OR. QPGRAD
      SKIP   = (.NOT. ((MOLGRD .OR. MOLHES .OR. DIPDER .OR. QPGRAD) 
     &         .AND.   DOSYM(1))) .OR. NUCDEP .GE. 10
      CUT    = .FALSE.
      NOROT1  = .FALSE.
      HESTRO = .TRUE.
      GDALL  = .FALSE.
      GDTRO  = .FALSE.
      RDTRO  = .FALSE.
      TROGRD = MOLGRD .AND. .NOT. (NFIELD .GT. 0)
      TROHES = MOLHES
      TRODIP = DIPDER
      RETURN
      END
C  /* Deck troinv */
      SUBROUTINE TROINV(WORK,LWORK)
C
C     This subroutine calculates the full molecular gradient, molecular
C     Hessian and dipole derivatives using translational and rotational
C     symmetries.
C
C     The method is described in Acta Chem. Scand. A42 (1988) 515
C
C     March 1985 tuh
C     Dipole derivatives added June 1985 tuh
C     Fully analytic expression for Hessian December 1986 tuh
C     Rewritten for symmetry Nov/Dec 1988 tuh
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
C
      DIMENSION WORK(LWORK)
C
#include <cbitro.h>
#include <abainf.h>
#include <energy.h>
#include <taymol.h>
#include <nuclei.h>
#include <dorps.h>
#include <infdim.h>
#include <infvar.h>
#include <exeinf.h>
#include <trkoor.h>
C
C     The variable FIRST has been replaced with FTRONV in /EXEINF/,
C     which can be controlled from outside. This is necessary because
C     TROINV must be executed twice in each iterations, the first time
C     being an initialization.
C
      IF (.NOT.FTRONV .AND. SKIP) RETURN
      CALL QENTER('TROINV')
C
C     Several variables are stored in /EXEINF/.
C
      IF (FTRONV) THEN
         GTRONV = GDTRO
         HTRONV = HESTRO
         RTRONV = NOROT1
      ELSE
         GDTRO  = GTRONV
         HESTRO = HTRONV
         NOROT1  = RTRONV
      END IF
C
      IF (IPRINT .GT. 0) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL TITLER('Output from TROINV','*',103)
      END IF
C
C     Reset SKIP if necessary
C
      IF (FTRONV) THEN
         NCOOR = 3*NUCDEP
         DO 100 I = 1, NCOOR
            SKIP = SKIP .OR. .NOT.DOPERT(I,1)
  100    CONTINUE
         IF (SKIP .AND. GDALL) THEN
            WRITE (LUPRI,'(/A/)') ' >>> Skip of TROINV is ignored '/
     *                /'because all GD vectors needed.'
            SKIP   = .FALSE.
            HESTRO = .FALSE.
            GDTRO  = .TRUE.
            RDTRO  = .FALSE.
            NOROT1  = .TRUE.
         END IF
      END IF
      IF (SKIP) THEN
         CALL SETPER(IPRINT)
         KTMAT = 1
         KAMAT = KTMAT + 6*NCOOR
         KTVEC = KAMAT + 6*NCOOR
         KLAST = KTVEC +   NCOOR
         IF (KLAST.GT.LWORK) CALL STOPIT('TROINV','GETTR',KLAST,LWORK)
         CALL GETTR(WORK(KTMAT),WORK(KAMAT),WORK(KTVEC))
      ELSE
C
C        Initialization
C
         IF (FTRONV) THEN
            KTMAT = 1
            KWRK  = KTMAT + 6*NCOOR
            LWRK  = LWORK - KWRK + 1
            CALL TRCOOR(WORK(KTMAT),WORK(KWRK),LWRK)
C
C        Calculation
C
         ELSE
            KHESTR = 1
            KDIPTR = KHESTR + MXCOOR*MXCOOR
            KWRK   = KDIPTR + 3*MXCOOR
            LWRK   = LWORK  - KWRK + 1
            IF (HESTRO) THEN
               CALL TRAROT(WORK(KHESTR),WORK(KDIPTR),WORK(KWRK),LWRK)
            END IF
C
C           GDALL and GDTRO are equivalent, but GDALL may be set outside
C           HESTRO module (851126/HJAAJ)
C
C           GDTRO = GDALL
            GDTRO = .FALSE.
C
            IF (RDTRO .OR. GDTRO) THEN
               IF (NOROT1 .OR. (GDALL .AND. .NOT.RDTRO)) THEN
                  IVCIND = 1
                  IVCDEP = IVCIND + NVAR
                  WRITE (LUPRI,'(//2A//)')
     *               ' Subroutine VECTRA has not been symmetrized yet.',
     *               ' Calculation cannot proceed.'
                  CALL QUIT('VECTRA not working.')
C                 CALL VECTRA(WORK(IVCIND),WORK(IVCDEP))
               ELSE
                  WRITE (LUPRI,'(//2A//)')
     *               ' Use of rotational symmetry for GD and/or RD',
     *               ' not implemented in this version.'
               END IF
            END IF
         END IF
      END IF
      CALL SETDCR('ABACUS')
      CALL SETGD
      FTRONV = .FALSE.
      IF (IPRINT .GT. 2) CALL TIMER('TROINV',TIMSTR,TIMEND)
      IF (CUT) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Program stopped after TROINV as required.'
         CALL QUIT(' ***** End of ABACUS (in TROINV) *****')
      END IF
      CALL QEXIT('TROINV')
      RETURN
      END
C  /* Deck setper */
      SUBROUTINE SETPER(IPRINT)
C
C     Sets DCORD, DCORGD, NGDVEC, IGDREC, IGDCOR, NGDTOT
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0)
C
#include <cbitrp.h>
#include <spnout.h>
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <dorps.h>
#include <gdvec.h>
C
#include <chrnos.h>
#include <ibtfun.h>
C
C     Some initialization for spin-spin couplings
C
      IF (.NOT. DOSELE) THEN
         DO 10 I = 1, MXCENT
            DOPERT(I,2) = .FALSE.
            NUCCHA = IZATOM(I)
            IF (NUCCHA .NE. 0) THEN
               DO 20 ISTP = 1, 5
                  GVAL   = DISOTP(NUCCHA,ISTP,'GVAL')
                  ABUND1 = DISOTP(NUCCHA,ISTP,'ABUNDANCE')
                  IF (.NOT. ((GVAL .EQ. D0) .OR.
     &                (ABUND1 .LT. ABUND))) DOPERT(I,2) = .TRUE.
  20           CONTINUE
            END IF
  10     CONTINUE
      END IF
C
C     DCORD & DCORGD
C
      DO 100 I = 1, 2
C        ... I = 1 geometric pertubations, I = 2 magnetic perturbations
         DOPERT(0,I) = .FALSE.
         DO 110 JATOM = 1, NUCIND
            DO 120 ICOOR = 1, 3
               JCOOR = 3*(JATOM - 1) + ICOOR
               DCORD (JATOM,ICOOR,I) = .FALSE.
               DCORGD(JATOM,ICOOR,I) = .FALSE.
               DO 130 IREP = 0, MAXREP
                  ISYMCR = IPTCNT(JCOOR,IREP,I)
                  IF (ISYMCR.NE.0) THEN
                     IF ((I .EQ. 1) .OR.
     &                   (I .EQ. 2) .AND. (VCD .OR. SHIELD)) THEN
                        DCORD (JATOM,ICOOR,I) = DCORD (JATOM,ICOOR,I)
     &                                          .OR. DOPERT(ISYMCR,I)
                        DCORGD(JATOM,ICOOR,I) = DCORGD(JATOM,ICOOR,I)
     &                                          .OR. DOPERT(ISYMCR,I)
                     ELSE
                        DCORD (JATOM,ICOOR,I) = DCORD (JATOM,ICOOR,I)
     &                                          .OR. DOPERT(JATOM,I)
                        DCORGD(JATOM,ICOOR,I) = DCORGD(JATOM,ICOOR,I)
     &                                          .OR. DOPERT(JATOM,I)
                     END IF
                  END IF
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 30) THEN
         CALL HEADER('DCORD',-1)
         DO 200 JATOM = 1, NUCIND
            WRITE (LUPRI,'(2(1X,3L5))')((DCORD(JATOM,I,J),I=1,3),J=1,2)
  200    CONTINUE
         CALL HEADER('DCORGD',-1)
         DO 210 JATOM = 1, NUCIND
            WRITE (LUPRI,'(2(1X,3L5))')((DCORGD(JATOM,I,J),I=1,3),J=1,2)
  210    CONTINUE
      END IF
C
C     Number of coordinates in each symmetry for which linear equations
C     will be solved
C
C     NGDVEC(8)             - # GD vectors in this symmetry
C     IGDCOR(3*NUCIND,IREP) - points from GD vector to coordinate
C     IGDREC(3*NUCIND,IREP) - points from GD vector to record
C     IDOREC(80*MXCENT)     - Restart information
C
      IF (.NOT. RESTAR) THEN
         CALL IZERO(IDORCT,80*MXCENT)
         CALL IZERO(IDORCI,48*(MXCENT + 1))
         CALL IZERO(IGDCOR,16*(MXCOOR + 3))
         CALL IZERO(IGDREC,16*(MXCOOR + 3))
         CALL IZERO(ITRCOR,80*MXCENT)
         CALL IZERO(ITRREC,80*MXCENT)
         CALL IZERO(NGDVEC,16)
      END IF
      DO 300 I = 1, 2
         NGDTOT(I) = 0
         IF (((I .EQ. 1) .AND. (MOLHES .OR. DIPDER .OR. QPGRAD)) .OR. 
     &       ((I .EQ. 2) .AND. ((SPNSPN .AND. DOPSO) .OR. SHIELD
     &                                       .OR. SPINRO))) THEN
            DO 310 IREP = 0, MAXREP
               NSYMCR = 0
               DO 320 JATOM = 1, NUCIND
                  DO 330 ICOOR = 1, 3
                     JCOOR = 3*(JATOM - 1) + ICOOR
                     ISCOOR = IPTCNT(JCOOR,IREP,I)
                     IF (ISCOOR .NE. 0) THEN
                        IF (I .EQ. 1) THEN
                           IF (DOPERT(ISCOOR,I)) THEN
                              NSYMCR = NSYMCR + 1
                              IGDCOR(NSYMCR,IREP + 1,I) = ISCOOR
                              IGDREC(NSYMCR,IREP + 1,I) = ISCOOR
                           END IF
                        ELSE
                           IF (DOPERT(JATOM,I) .OR. SHIELD 
     &                          .OR. SPINRO) THEN
                              NSYMCR = NSYMCR + 1
                              IGDCOR(NSYMCR,IREP + 1,I) = ISCOOR
                              IGDREC(NSYMCR,IREP + 1,I) = ISCOOR
                              IF (.NOT.NCSPNI(JATOM)) IDORCI(ISCOOR,2)=2
                           END IF
                        END IF
                     END IF
 330              CONTINUE
 320           CONTINUE
               NGDVEC(IREP + 1,I) = NSYMCR
 310        CONTINUE
            NGDTOT(I) = ISUM(MAXREP+1,NGDVEC(1,I),1)
         END IF
 300  CONTINUE
C     
C     Dipole derivatives
C
      DO 400 I = 1, 2
         IF ((I.EQ.1 .AND. (DIPDER .OR. POLAR  .OR. VCD)) .OR.
     &       (I.EQ.2 .AND. (SHIELD .OR. MAGSUS .OR. VCD .OR.
     &                      MOLGFA .OR. SPINRO .OR. VROA .OR.
     &                      OPTROT))) THEN
            IAXSYM = 0
            IAXREC = 3*NUCDEP
            DO 410 IREP = 0, MAXREP
               NSYMCR = NGDVEC(IREP + 1,I)
               DO 420 IAX = 1, NAXREP(IREP,I)
                  IAXSYM = IAXSYM + 1
                  IAXREC = IAXREC + 1
                  NSYMCR = NSYMCR + 1
                  IGDCOR(NSYMCR,IREP + 1,I) = - IAXSYM
                  IGDREC(NSYMCR,IREP + 1,I) =   IAXREC
 420           CONTINUE
               NGDVEC(IREP + 1,I) = NSYMCR
 410        CONTINUE
            NGDTOT(I) = ISUM(MAXREP+1,NGDVEC(1,I),1)
         END IF
 400  CONTINUE
C
C     Quadrupole derivatives
C
      IF (QPGRAD) THEN
         IAXREC = 3*NUCDEP + 3
         DO IX = 1, 3
            DO IY = IX , 3
               IREP = IBTXOR(ISYMAX(IX,1),ISYMAX(IY,1))
               NSYMCR = NGDVEC(IREP + 1,1)
               IAXREC = IAXREC + 1
               NSYMCR = NSYMCR + 1
               IAXSYM = IX*10 + IY
               IGDCOR(NSYMCR,IREP + 1,1) = - IAXSYM
               IGDREC(NSYMCR,IREP + 1,1) =   IAXREC
               NGDVEC(IREP + 1,1) = NSYMCR
            END DO
         END DO
         NGDTOT(1) = ISUM(MAXREP + 1,NGDVEC(1,1),1)
      END IF

C
C     Number of triplet operators in each symmetry for which linear equations
C     will be solved
C
C     We first take the spin-dipole operators, then the Fermi contact
C
      IATOM  = 0
      ISYMCR = 0
      DO 421 IREP1 = 0, MAXREP
         NSYMCR = 0
         IF (DOSD .OR. DOSDFC) THEN
            DO 423 IREP2 = 0, MAXREP
               DO 425 JATOM = 1, NUCIND
                  IF (DOPERT(JATOM,2)) THEN
                     DO 427 ICOOR1 = 1, 3
                        ISCOR1 = IPTCNT(3*(JATOM - 1) + ICOOR1,IREP2,2)
                        IF (ISCOR1 .GT. 0) THEN
                        DO 429 ICOOR2 = 1, ICOOR1
                           ISYM = IBTXOR(IREP2,ISYMAX(ICOOR2,2))
                           IF (ISYM .EQ. IREP1) THEN
                              NSYMCR = NSYMCR + 1
                              ISYMCR = ISYMCR + 1
                              ISCOR2 = 3*(ISCOR1 - 1) + ICOOR2
                              ITRCOR(NSYMCR,IREP1 + 1) = ISCOR2
                              ITRREC(NSYMCR,IREP1 + 1) = ISYMCR
                              IF (.NOT.NCSPNI(JATOM)) IDORCT(ISYMCR) = 1
                           END IF
 429                    CONTINUE
                        END IF
 427                 CONTINUE
                  END IF
 425           CONTINUE
 423        CONTINUE
         END IF
         IF (DOFC) THEN
            DO 431 JATOM = 1, NUCIND
               IF (IBTAND(IREP1,ISTBNU(JATOM)).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOPERT(JATOM,2)) THEN
                     NSYMCR = NSYMCR + 1
                     ISYMCR = ISYMCR + 1
                     ITRCOR(NSYMCR,IREP1 + 1) = - IATOM
                     ITRREC(NSYMCR,IREP1+1) = ISYMCR
                     IF (.NOT.NCSPNI(JATOM)) IDORCT(ISYMCR) = 1
                  END IF
               END IF
 431        CONTINUE
         END IF
         NTRVEC(IREP1 + 1) = NSYMCR
 421  CONTINUE
C
      IF (IPRINT .GT. 30) THEN
         CALL AROUND('Information in COMMON/GDVEC/')
         DO 500 I = 1, 2
            IF (I .EQ. 1) THEN
               WRITE (LUPRI,'(//1X,A)')
     &             ' ... for geometrical and electric perturbations:'
            ELSE
               WRITE (LUPRI,'(//1X,A)')
     &             ' ... for magnetic perturbations:'
            END IF
            WRITE (LUPRI,'(//1X,A, I5)') ' NGDTOT ',NGDTOT(I)
            DO 510 ISYM = 1, MAXREP + 1
               CALL HEADER('Symmetry '//CHRNOS(ISYM),1)
               WRITE (LUPRI,'(1X,A,I5)') ' NGDVEC(ISYM,I) ',
     &                                     NGDVEC(ISYM,I)
               WRITE (LUPRI,'(1X,/A,(10I5))') ' IGDCOR(*,ISYM,I) ',
     &                         (IGDCOR(J,ISYM,I),J=1,NGDVEC(ISYM,I))
               WRITE (LUPRI,'(1X,/A,(10I5))') ' IGDREC(*,ISYM,I) ',
     &                         (IGDREC(J,ISYM,I),J=1,NGDVEC(ISYM,I))
 510        CONTINUE
 500     CONTINUE
      END IF
C
      DO 600 I = 1, 2
         NDCORD(I) = 0
         DO 610 JATOM = 1,NUCIND
            DO 620 ICOOR = 1, 3
               IF (DCORD(JATOM,ICOOR,I)) THEN
                  NDCORD(I) = NDCORD(I) + 1
               END IF
  620       CONTINUE
  610    CONTINUE
  600 CONTINUE
      RETURN
      END
C  /* Deck trcoor */
      SUBROUTINE TRCOOR(TMAT,WORK,LWORK)
C
C     This subroutine determines a set of independent translations and
C     rotations, and a corresponding set of Cartesian coordinates.
C
C     tuh Dec/Nov 1988
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      LOGICAL TRCORD(MXCENT,3)
      DIMENSION TMAT(6,NCOOR), WORK(LWORK)
C
#include <abainf.h>
#include <cbitro.h>
#include <nuclei.h>
#include <dorps.h>
#include <symmet.h>
#include <trkoor.h>
C
C     *************************************************
C     ***** Determine trarot coordinates and TMAT *****
C     *************************************************
C
      KAMAT = 1
      KTVEC = KAMAT + 6*NCOOR
      KLAST = KTVEC +   NCOOR
      IF (KLAST.GT.LWORK)
     *     CALL STOPIT('TRCOOR','before GETTR',KLAST,LWORK)
      CALL GETTR(TMAT,WORK(KAMAT),WORK(KTVEC))
C
C     ********************************************************
C     ***** Determine Cartesian directions to be skipped *****
C     ********************************************************
C
      KAMAT = 1
      KDONE = KAMAT + 6*NCOOR
      KPREF = KDONE +   NCOOR
      KTDEP = KPREF +   NCOOR
      KLAST = KTDEP +   NCOOR
      IF (KLAST .GT. LWORK)
     *     CALL STOPIT('TRCOOR','before SETDEP',KLAST,LWORK)
      CALL SETDEP(TMAT,WORK(KAMAT),WORK(KDONE),WORK(KPREF),WORK(KTDEP))
      IF (COMPAR) THEN
         DO 100 ISCOOR = 1, NCOOR
            DOPERT(ISCOOR,1) = .TRUE.
  100    CONTINUE
      ELSE
         DO 200 ISCOOR = 1, NCOOR
            DOPERT(ISCOOR,1) = .NOT.DEPEND(ISCOOR)
  200    CONTINUE
      END IF
      CALL SETPER(IPRINT)
      IF (IPRINT .GE. 1) THEN
        WRITE (LUPRI,'(//)')
        CALL HEADER('Dependent and independent Cartesian coordinates',1)
         DO 600 IREP = 0, MAXREP
            IF (NCRREP(IREP,1) .GT. 0) THEN
               WRITE (LUPRI,'(/1X,A,I1/)') ' Symmetry ',IREP+1
               DO 610 IATOM = 1, NUCIND
                  DO 620 ICOOR = 1, 3
                     ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                     IF (ISCOOR .GT. 0) THEN
                        IF (DEPEND(ISCOOR)) THEN
                           WRITE (LUPRI,'(8X,A6,A)')
     &                      NAMEX(IPTCOR(ISCOOR,1)),'   -    dependent '
                        ELSE
                           WRITE (LUPRI,'(8X,A6,A)')
     &                      NAMEX(IPTCOR(ISCOOR,1)),'   -  independent '
                        END IF
                     END IF
  620             CONTINUE
  610          CONTINUE
            END IF
  600    CONTINUE
      END IF
      RETURN
      END
C  /* Deck setdep */
      SUBROUTINE SETDEP(TMAT,AMAT,DONE,NPREF,ITBDEP)
C
C     Determines a set of Cartesian coordinates whose gradient and
C     Hessian elements can be determined by translational and rotational
C     symmetry.
C
C     Dec 1988 tuh
C
#include <implicit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL OK, OKALL, ONCE, DONE(NCOOR)
      DIMENSION TMAT(6,NCOOR), AMAT(6,NCOOR), NPREF(NCOOR),
     *          ITBDEP(NCOOR)
C
#include <cbitro.h>
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
C
      IF (IPRINT .GE. 5) THEN
         CALL TITLER('Output from SETDEP','*',103)
      END IF
      NSET = 0
      DO 100 I = 1, NCOOR
         DONE(I) = .FALSE.
  100 CONTINUE
      CALL DZERO(AMAT,6*NCOOR)
C
C     (A) Eliminate if possible one atom completely
C     =============================================
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Coordinates for atom completely eliminated',-1)
      END IF
C
C     Select atom of multiplicity one and highest charge
C
      MATOM  = 0
      CHAMAX = D0
      DO 200 IATOM = 1, NUCIND
         IF (MULT(ISTBNU(IATOM)) .EQ. 1) THEN
            CHARGI = CHARGE(IATOM)
            IF (CHARGI .GT. CHAMAX) THEN
               MATOM  = IATOM
               CHAMAX = CHARGI
            END IF
         END IF
  200 CONTINUE
C
C     Set coordinates of this atom - if any - as dependent
C
      IF (MATOM .GT. 0) THEN
         DO 210 ICOOR = 1, 3
            DO 220 IREP = 0, MAXREP
               ISCOOR = IPTCNT(3*(MATOM - 1) + ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  CALL DCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
                  IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *                 ' Calling LINDEP for coordinate ',ISCOOR
                  CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
                  IF (OK) THEN
                     NSET = NSET + 1
                     ITBDEP(NSET) = ISCOOR
                     DONE(ISCOOR) = .TRUE.
                     IF (NSET .EQ. NCDEP) GO TO 1000
                  END IF
               END IF
  220       CONTINUE
  210    CONTINUE
      END IF
      IF (IPRINT .GE. 5) THEN
         WRITE (LUPRI,'(/A,I5)')
     *      ' Number of dependent coordinates found so far: ',NSET
         IF (NSET .GT. 0) WRITE (LUPRI,'(A,6I5)') ' Coordinates found:',
     *                                              (ITBDEP(I),I=1,NSET)
      END IF
C
C     (B) Eliminate if possible one coordinate for all symmetries
C     ===========================================================
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Coordinates eliminated for all symmetries',-1)
      END IF
C
C     Identify Cartesian directions which may be eliminated for
C     all symmetries
C
      NDIR = 0
      IDIR = 0
      DO 300 IATOM = 1, NUCIND
         DO 310 ICOOR = 1, 3
            IDIR = IDIR + 1
            OKALL = .TRUE.
            ONCE  = .FALSE.
            DO 320 IREP = 0, MAXREP
               ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  IF (.NOT.DONE(ISCOOR)) THEN
                     ONCE = .TRUE.
                     CALL DCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
                     IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *                    ' Calling LINDEP for coordinate ',ISCOOR
                     CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
                     OKALL = OKALL .AND. OK
                  END IF
               END IF
  320       CONTINUE
            IF (ONCE .AND. OKALL) THEN
               NDIR = NDIR + 1
               NPREF(NDIR) = IDIR
            END IF
  310    CONTINUE
  300 CONTINUE
      IF (IPRINT .GE. 5) THEN
         IF (NDIR .GT. 0) THEN
            CALL HEADER('Unsorted list of directions',1)
            WRITE (LUPRI,'(12I5)') (NPREF(I),I=1,NDIR)
         ELSE
            WRITE (LUPRI,'(A)') ' No coordinates found.'
         END IF
      END IF
C
C     Sort these coordinates by charge of their atoms
C
      DO 330 I = 1, NDIR - 1
         IDIR = NPREF(I)
         CHARGI = CHARGE((IDIR + 2)/3)
         DO 340 J = I + 1, NDIR
            JDIR = NPREF(J)
            CHARGJ = CHARGE((JDIR + 2)/3)
            IF (CHARGJ .GT. CHARGI) THEN
               NPREF(I) = JDIR
               NPREF(J) = IDIR
            END IF
  340    CONTINUE
  330 CONTINUE
      IF (IPRINT .GE. 5) THEN
         IF (NDIR .GT. 0) THEN
            CALL HEADER('Sorted list of directions',1)
            WRITE (LUPRI,'(12I5)') (NPREF(I),I=1,NDIR)
         END IF
      END IF
C
C     Set these coordinates - if any - as dependent
C
      DO 350 IDIR = 1, NDIR
         DO 360 IREP = 0, MAXREP
            ISCOOR = IPTCNT(NPREF(IDIR),IREP,1)
            IF (ISCOOR .GT. 0) THEN
               IF (.NOT.DONE(ISCOOR)) THEN
                  CALL DCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
                  IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *                 ' Calling LINDEP for coordinate ',ISCOOR
                  CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
                  IF (OK) THEN
                     NSET = NSET + 1
                     ITBDEP(NSET) = ISCOOR
                     DONE(ISCOOR) = .TRUE.
                     IF (NSET .EQ. NCDEP) GO TO 1000
                  END IF
               END IF
            END IF
  360    CONTINUE
  350 CONTINUE
      IF (IPRINT .GE. 5) THEN
         WRITE (LUPRI,'(/A,I5)')
     *      ' Number of dependent coordinates found so far: ',NSET
         IF (NSET .GT. 0) WRITE (LUPRI,'(A,6I5)') ' Coordinates found:',
     *                                              (ITBDEP(I),I=1,NSET)
      END IF
C
C     (D) Set remaining coordinates until all have been set
C     =====================================================
C
      IF (IPRINT .GE. 5) CALL HEADER('Remaining coordinates',-1)
C
C     Sort coordinates by charge of their atoms
C
      NDIR = 0
      DO 400 I = 1, NCOOR
         IF (.NOT.DONE(I)) THEN
            NDIR = NDIR + 1
            NPREF(NDIR) = I
         END IF
  400 CONTINUE
      DO 410 I = 1, NDIR
         DO 420 J = I + 1, NDIR
            ISCOOR = NPREF(I)
            CHARGI = CHARGE((IPTCOR(ISCOOR,1) + 2)/3)
            JSCOOR = NPREF(J)
            CHARGJ = CHARGE((IPTCOR(JSCOOR,1) + 2)/3)
            IF (CHARGJ .GT. CHARGI) THEN
               NPREF(I) = JSCOOR
               NPREF(J) = ISCOOR
            END IF
  420    CONTINUE
  410 CONTINUE
      IF (IPRINT .GE. 5) THEN
         IF (NDIR .GT. 0) THEN
            CALL HEADER('Sorted list of directions',1)
            WRITE (LUPRI,'(12I5)') (NPREF(I),I=1,NDIR)
         END IF
      END IF
      DO 430 IDIR = 1, NDIR
         ISCOOR = NPREF(IDIR)
         CALL DCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
         IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *       ' Calling LINDEP for coordinate ',ISCOOR
         CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
         IF (OK) THEN
            NSET = NSET + 1
            ITBDEP(NSET) = ISCOOR
            DONE(ISCOOR) = .TRUE.
            IF (NSET .EQ. NCDEP) GO TO 1000
         END IF
  430 CONTINUE
 1000 CONTINUE
      IF (NSET .NE. NCDEP) THEN
         WRITE (LUPRI,'(/,2A,/,A,/)')
     &         ' ERROR in SETDEP - less than 6 (5)',
     &         ' independent Cartesian directions',' found.'
         WRITE (LUPRI,'(A,E12.6)')
     &         ' Current threshold for dependency: ',THRESH
         WRITE (LUPRI,'(2A)')
     &         ' Attempt calculation with decreased',
     &         ' threshold .THRESH under *TROINV' 
         WRITE (LUPRI,'(A)')
     &         ' Alternatively use .SKIP under *TROINV'
         WRITE (LUPRI,'(2A)')
     &         ' Note: These directives must be used',
     &         ' under **START, **EACH STEP and **PROPERTIES.'
         CALL QUIT('ERROR in SETDEP.')
      END IF
C
C     Set up DEPEND(NCOOR)
C
      ICOOR = 0
      DO 500 ISCOOR = 1, NCOOR
         DEPEND(ISCOOR) = .FALSE.
  500 CONTINUE
      DO 510 IDEP = 1, NCDEP
         DEPEND(ITBDEP(IDEP)) = .TRUE.
  510 CONTINUE
      IF (IPRINT .GE. 5) THEN
         WRITE (LUPRI,'(/A,I5)')
     *      ' Total number of dependent coordinates found: ',NSET
         IF (NSET .GT. 0) WRITE (LUPRI,'(A,6I5)') ' Coordinates found:',
     *                                              (ITBDEP(I),I=1,NSET)
         CALL HEADER('DEPEND in SETDEP',1)
         WRITE (LUPRI,'(2X,35L2)') (DEPEND(I),I=1,NCOOR)
      END IF
      RETURN
      END
C  /* Deck trarot */
      SUBROUTINE TRAROT(HESTR,DIPTR,WORK,LWORK)
C
C     tuh Dec 1988
C
C     This subroutine determines the full molecular gradient and
C     Hessian using translational and rotational symmetry.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      DIMENSION WORK(LWORK), HESTR(MXCOOR,MXCOOR), DIPTR(3,MXCOOR)
#include <abainf.h>
#include <cbitro.h>
#include <taymol.h>
#include <nuclei.h>
#include <moldip.h>
#include <trkoor.h>
#include <symmet.h>
#include <chrxyz.h>
      IF (HESTRO) THEN
         CALL DZERO(HESTR,MXCOOR*MXCOOR)
         DO 50 I = 1, NCOOR
            IF (.NOT.DEPEND(I)) THEN
               DO 60 J = 1, I
                  IF (.NOT.DEPEND(J)) THEN
                     HESSIJ = HESMOL(I,J)
                     HESTR(I,J) = HESSIJ
                     HESTR(J,I) = HESSIJ
                  END IF
   60          CONTINUE
            END IF
   50    CONTINUE
      END IF
      IF (TRODIP) CALL DZERO(DIPTR,3*MXCOOR)
C
C     Loop over symmetries
C
      IOFFCR = 0
      IOFFAX = 0
      IRPEND = MAXREP
      IF (.NOT.(TROHES.OR.TRODIP)) IRPEND = 0
      DO 100 IREP = 0, IRPEND
      NCR  = NCRREP(IREP,1)
      NDEP = NTRREP(IREP)
      NAX  = NAXREP(IREP,1)
      IF (DOSYM(IREP + 1)) THEN
         NIND = NCR - NDEP
         IF (IPRINT .GE. 5) THEN
            WRITE (LUPRI,'(//A,I5)') ' Symmetry:                  ',IREP
            WRITE (LUPRI,'(A,I5)') ' Number of dep. coordinates:',NDEP
            WRITE (LUPRI,'(A,I5)') ' Number of ind. coordinates:',NIND
         END IF
         IF (NDEP .GT. 0) THEN
C
C           ********************
C           ***** T matrix *****
C           ********************
C
            KPTCOL = 1
            KTDD   = KPTCOL + NCOOR
            KTDINV = KTDD   + NDEP*NDEP
            KTDI   = KTDINV + NDEP*NDEP
            KPVT   = KTDI   + NDEP*NIND
            KSCR   = KPVT   + NDEP
            KLAST  = KSCR   + NDEP
            IF (KLAST .GT. LWORK)
     *          CALL STOPIT('TRAROT','T matrix',KLAST,LWORK)
C
C           TDD
C
            CALL GETTRO(WORK(KTDD),WORK(KPTCOL),NDEP,NDEP,
     *                  'DEP','TMAT','TC',IREP,IPRINT)
            IF (IPRINT .GE. 10) THEN
               CALL AROUND('TDD - TRAROT')
               CALL OUTPUT(WORK(KTDD),1,NDEP,1,NDEP,NDEP,NDEP,1,LUPRI)
            END IF
C
C           TDDINV
C
            CALL DGEINV(NDEP,WORK(KTDD),WORK(KTDINV),
     *                  WORK(KPVT),WORK(KSCR),INFO)
            IF (INFO .NE. 0) THEN
               WRITE (LUPRI,'(//,A,I5,A,/)') ' INFO =', INFO,
     *                                       ' from DGEINV in TRAROT '
               CALL QUIT('ERROR in DGESOL.')
            END IF
            IF (IPRINT .GE. 10) THEN
               CALL AROUND('TDDINV - TRAROT')
               CALL OUTPUT(WORK(KTDINV),1,NDEP,1,NDEP,NDEP,NDEP,1,LUPRI)
            END IF
C
C           TDI
C
            IF (NIND .GT. 0) THEN
               CALL GETTRO(WORK(KTDI),WORK(KPTCOL),NDEP,NIND,
     *                     'IND','TMAT','TC',IREP,IPRINT)
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('TDI - TRAROT')
                  CALL OUTPUT(WORK(KTDI),1,NDEP,1,NIND,
     *                        NDEP,NIND,1,LUPRI)
               END IF
            END IF
C
C           ******************************
C           ***** Molecular Gradient *****
C           ******************************
C
            IF (TROGRD .AND. IREP.EQ.0) THEN
               IF (IPRINT .GE. 5) CALL HEADER('Gradient Section',-1)
               KFI    = KTDI   + NDEP*NIND
               KTDIFI = KFI    + NIND
               KFD    = KTDIFI + NDEP
               KLAST  = KFD    + NDEP
               IF (KLAST .GT. LWORK)
     *             CALL STOPIT('TRAROT','gradient',KLAST,LWORK)
C
C              FI vector
C
               IIND = 0
               DO 200 ISCOOR = 1, NCR
                  IF (.NOT.DEPEND(ISCOOR)) THEN
                     WORK(KFI + IIND) = GRDMOL(ISCOOR)
                     IIND = IIND + 1
                  END IF
  200          CONTINUE
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('FI - TRAROT')
                  CALL OUTPUT(WORK(KFI),1,NIND,1,1,NIND,1,1,LUPRI)
               END IF
C
C              FD vector
C
               CALL DZERO(WORK(KFD),NDEP)
               IF (NIND .GT. 0) THEN
                  CALL DGEMM('N','N',NDEP,1,NIND,1.D0,
     &                       WORK(KTDI),NDEP,
     &                       WORK(KFI),NIND,0.D0,
     &                       WORK(KTDIFI),NDEP)
                  CALL DGEMM('N','N',NDEP,1,NDEP,-1.D0,
     &                       WORK(KTDINV),NDEP,
     &                       WORK(KTDIFI),NDEP,1.D0,
     &                       WORK(KFD),NDEP)
                  IF (IPRINT .GE. 10) THEN
                     CALL AROUND('FD - TRAROT')
                     CALL OUTPUT(WORK(KFD),1,NDEP,1,1,NDEP,1,1,LUPRI)
                  END IF
C
C              Compare FD with explicitly calculated gradient
C     
                  CALL CMPGRD(WORK(KFD),NDEP,IPRINT)
               END IF
            END IF

C           *****************************
C           ***** Molecular Hessian *****
C           *****************************
C
            IF (TROHES) THEN
               IF (IPRINT .GE. 5) CALL HEADER('Hessian Section',-1)
               KBDD  = KTDI + NDEP*NIND
               KBDI  = KBDD + NDEP*NDEP
               KHII  = KBDI + NDEP*NIND
               KHDI  = KHII + NIND*NIND
               KHDD  = KHDI + NDEP*NIND
               KSCR  = KHDD + NDEP*NDEP
               KLAST = KSCR + NDEP*NDEP
               IF (KLAST .GT. LWORK)
     *             CALL STOPIT('TRAROT','Hessian',KLAST,LWORK)
C
C              BDD
C
               CALL GETTRO(WORK(KBDD),WORK(KPTCOL),NDEP,NDEP,
     *                     'DEP','HESS','TC',IREP,IPRINT)
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('BDD - TRAROT')
                  CALL OUTPUT(WORK(KBDD),1,NDEP,1,NDEP,
     *                        NDEP,NDEP,1,LUPRI)
               END IF
C
C              HDD = TDD-1*BDD
C
               CALL DGEMM('N','N',NDEP,NDEP,NDEP,1.D0,
     &                    WORK(KTDINV),NDEP,
     &                    WORK(KBDD),NDEP,0.D0,
     &                    WORK(KHDD),NDEP)
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('HDD (DD contributions only) - TRAROT')
                  CALL OUTPUT(WORK(KHDD),1,NDEP,1,NDEP,
     *                        NDEP,NDEP,1,LUPRI)
               END IF
               IF (NIND .GT. 0) THEN
C
C                 BDI
C
                  CALL GETTRO(WORK(KBDI),WORK(KPTCOL),NDEP,NIND,
     *                        'IND','HESS','TC',IREP,IPRINT)
                  IF (IPRINT .GE. 10) THEN
                     CALL AROUND('BDI - TRAROT')
                     CALL OUTPUT(WORK(KBDI),1,NDEP,1,NIND,
     *                           NDEP,NIND,1,LUPRI)
                  END IF
C
C                 HII
C
                  IIND = 0
                  DO 300 I = IOFFCR + 1, IOFFCR + NCR
                     IF (.NOT.DEPEND(I)) THEN
                        IIND = IIND + 1
                        JIND = 0
                        DO 310 J = IOFFCR + 1, I
                           IF (.NOT.DEPEND(J)) THEN
                              JIND = JIND + 1
                              IJ   = (IIND - 1)*NIND + JIND - 1
                              JI   = (JIND - 1)*NIND + IIND - 1
                              HESSIJ = HESMOL(I,J)
                              WORK(KHII + IJ) = HESSIJ
                              WORK(KHII + JI) = HESSIJ
                           END IF
  310                   CONTINUE
                     END IF
  300             CONTINUE
                  IF (IPRINT .GE. 10) THEN
                     CALL AROUND('HII - TRAROT')
                     CALL OUTPUT(WORK(KHII),1,NIND,1,NIND,
     *                           NIND,NIND,1,LUPRI)
                  END IF
C
C                 HDI
C
                  CALL DGEMM('N','N',NDEP,NIND,NIND,-1.D0,
     &                       WORK(KTDI),NDEP,
     &                       WORK(KHII),NIND,1.D0,
     &                       WORK(KBDI),NDEP)
                  CALL DGEMM('N','N',NDEP,NIND,NDEP,1.D0,
     &                       WORK(KTDINV),NDEP,
     &                       WORK(KBDI),NDEP,0.D0,
     &                       WORK(KHDI),NDEP)
                  IF (IPRINT .GE. 10) THEN
                     CALL AROUND('HDI - TRAROT')
                     CALL OUTPUT(WORK(KHDI),1,NDEP,1,NIND,
     *                           NDEP,NIND,1,LUPRI)
                  END IF
C
C                 HDD
C
                  CALL DGEMM('N','T',NDEP,NDEP,NIND,1.D0,
     &                       WORK(KTDI),NDEP,
     &                       WORK(KHDI),NDEP,0.D0,
     &                       WORK(KSCR),NDEP)
                  CALL DGEMM('N','N',NDEP,NDEP,NDEP,-1.D0,
     &                       WORK(KTDINV),NDEP,
     &                       WORK(KSCR),NDEP,1.D0,
     &                       WORK(KHDD),NDEP)
               END IF
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('HDD - TRAROT')
                  CALL OUTPUT(WORK(KHDD),1,NDEP,1,NDEP,
     *                        NDEP,NDEP,1,LUPRI)
               END IF
C
C              Add calculated elements to HESTR
C
               IDEP = 0
               DO 400 I = IOFFCR + 1, IOFFCR + NCR
                  IF (DEPEND(I)) THEN
                     IDEP = IDEP + 1
                     JDEP = 0
                     JIND = 0
                     DO 410 J = IOFFCR + 1, IOFFCR + NCR
                        IF (DEPEND(J)) THEN
                           JDEP = JDEP + 1
                           IJ   = (JDEP - 1)*NDEP + IDEP - 1
                           HESTR(I,J) = WORK(KHDD + IJ)
                        ELSE
                           JIND = JIND + 1
                           IJ   = (JIND - 1)*NDEP + IDEP - 1
                           HESTR(I,J) = WORK(KHDI + IJ)
                           HESTR(J,I) = WORK(KHDI + IJ)
                        END IF
  410                CONTINUE
                  END IF
  400          CONTINUE
            END IF
C
C           ***************************
C           ***** Dipole Gradient *****
C           ***************************
C
            IF (TRODIP .AND. NAX.GT.0) THEN
               IF (IPRINT .GE. 5) THEN
                  CALL HEADER('Dipole Section',-1)
                  WRITE (LUPRI,'(1X,A,3(A,1X))')
     *                  'Components of this symmetry: ',
     *                  (CHRXYZ(-IPTXYZ(I,IREP,1)),I=1,NAX)
               END IF
C
C              Work space allocations
C
               KDGI   = KTDI   + NDEP*NIND
               KBMAT  = KDGI   + NIND*NAX
               KDGD   = KBMAT  + NDEP*NAX
               KLAST  = KDGD   + NDEP*NAX
               IF (KLAST .GT. LWORK)
     *            CALL STOPIT('TRAROT','dipole',KLAST,LWORK)
C
C              DGI vector
C
               IIND = 0
               DO 500 IAX = IOFFAX + 1, IOFFAX + NAX
                  DO 510 ISCOOR = IOFFCR + 1, IOFFCR + NCR
                     IF (.NOT.DEPEND(ISCOOR)) THEN
                        WORK(KDGI + IIND) = DIP1(IAX,ISCOOR)
                        IIND = IIND + 1
                     END IF
  510             CONTINUE
  500          CONTINUE
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('DGI - TRAROT')
                  CALL OUTPUT(WORK(KDGI),1,NIND,1,NAX,NIND,NAX,1,LUPRI)
               END IF
C
C              B matrix
C
               CALL GETTRO(WORK(KBMAT),WORK(KPTCOL),NDEP,NAX,
     *                     'DEP','DIPOLE','TC',IREP,IPRINT)
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('B matrix - TRAROT')
                  CALL OUTPUT(WORK(KBMAT),1,NDEP,1,NAX,NDEP,NAX,1,LUPRI)
               END IF
C
C              DGD vector
C
               IF (NIND .GT. 0) CALL DGEMM('N','N',NDEP,NAX,NIND,-1.D0,
     &                       WORK(KTDI),NDEP,
     &                       WORK(KDGI),NIND,1.D0,
     &                       WORK(KBMAT),NDEP)
               CALL DGEMM('N','N',NDEP,NAX,NDEP,1.D0,
     &                    WORK(KTDINV),NDEP,
     &                    WORK(KBMAT),NDEP,0.D0,
     &                    WORK(KDGD),NDEP)
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('DGD - TRAROT')
                  CALL OUTPUT(WORK(KDGD),1,NDEP,1,NAX,NDEP,NAX,1,LUPRI)
               END IF
C
C              Construct DIPTR
C
               IJ = 0
               DO 600 I = IOFFAX + 1, IOFFAX + NAX
                  DO 610 J = IOFFCR + 1, IOFFCR + NCR
                     IF (DEPEND(J)) THEN
                        DIPTR(I,J) = WORK(KDGD + IJ)
                        IJ   = IJ   + 1
                     ELSE
                        DIPTR(I,J) = DIP1(I,J)
                     END IF
  610             CONTINUE
  600          CONTINUE
            END IF
         END IF
      END IF
      IOFFCR = IOFFCR + NCR
      IOFFAX = IOFFAX + NAX
  100 CONTINUE
      KCSTRA = 1
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('TRAROT','TRACOR',KLAST,LWORK)
      IF (TROHES) THEN
         IF (IPRINT .GE. 5) THEN
            CALL HEADER('Total molecular Hessian - TRAROT',-1)
            CALL PRIHES(HESTR,'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (COMPAR) THEN
            CALL CMPHES(HESTR,IPRINT,WORK(KCSTRA),WORK(KSCTRA))
            DO 700 I = 1, NCOOR
            DO 700 J = 1, I - 1
               HESMOL(J,I) = HESMOL(I,J)
  700       CONTINUE
         ELSE
            CALL DCOPY(MXCOOR*MXCOOR,HESTR,1,HESMOL,1)
         END IF
      END IF
      IF (TRODIP) THEN
         IF (IPRINT .GE. 5) THEN
            CALL HEADER('Total dipole gradient - TRAROT',-1)
            CALL FCPRI(DIPTR,'APT',WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (COMPAR) THEN
            CALL CMPDIP(DIPTR,IPRINT,WORK(KCSTRA),WORK(KSCTRA))
         ELSE
            CALL DCOPY(3*NCOOR,DIPTR,1,DIP1,1)
         END IF
      END IF
      RETURN
      END
C  /* Deck cmpgrd */
      SUBROUTINE CMPGRD(FD,NDEP,IPRINT)
C
C     Compares the elements of the molecular gradient calculated by
C     translational and rotational symmetry with those calculated
C     directly.
C
C     tuh 071288
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      DIMENSION FD(NDEP)
      LOGICAL   LPRINT
#include <abainf.h>
#include <taymol.h>
#include <nuclei.h>
#include <trkoor.h>
#include <symmet.h>
      PARAMETER (D0 = 0.0D0)
      IF (IPRINT .GT. 0) THEN
         LPRINT = .TRUE.
      ELSE
         LPRINT = .FALSE.
      END IF
   10 IF (LPRINT) THEN
         CALL AROUND(
     *    'Gradient elements calculated directly and from symmetry')
         CALL HEADER(
     * 'Coordinate    Direct calc.     Tra.-rot. sym.   Difference  ',7)
      END IF
      IIND = 0
      DIFMAX = D0
      DO 100 ISCOOR = 1, NCRREP(0,1)
         IF (DEPEND(ISCOOR)) THEN
            IIND = IIND + 1
            GRDML = GRDMOL(ISCOOR)
            GRDTR = FD(IIND)
            DIFFER = GRDML - GRDTR
            ABSDIF = ABS(DIFFER)
            IF (ABSDIF .GE. DIFMAX) THEN
               DIFMAX = ABSDIF
               IDIFMX = ISCOOR
            END IF
            IF (LPRINT) THEN
               WRITE (LUPRI,'(8X,A6,3X,3F17.10)')
     &                NAMEX(IPTCOR(ISCOOR,1)), GRDML, GRDTR, DIFFER
            END IF
         END IF
  100 CONTINUE
      IF (DIFMAX.GT.1.D-06) THEN
         IF (.NOT.LPRINT) THEN
            LPRINT = .TRUE.
            GO TO 10
         END IF
         WRITE (LUPRI,'(///A)') ' >>> WARNING from TRAROT:'
         NWNABA = NWNABA + 1
      END IF
      IF (LPRINT .OR. DIFMAX.GT.1.D-12)
     &   WRITE(LUPRI,'(/A,I3,A,1P,E10.2)')
     *      ' Largest deviation from tra-rot symmetry is in'//
     *      ' gradient element',IDIFMX,':', DIFMAX
      RETURN
      END
C  /* Deck cmphes */
      SUBROUTINE CMPHES(HESTR,IPRINT,CSTRA,SCTRA)
C
C     Compares the elements of the molecular Hessian calculated by
C     translational and rotational symmetry with those calculated
C     directly.
C
C     tuh 071288
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION HESTR(MXCOOR,MXCOOR), CSTRA(*), SCTRA(*)
#include <abainf.h>
#include <taymol.h>
#include <trkoor.h>
      PARAMETER (D0 = 0.0D0)
      DIFMAX = D0
      DO 100 I = 1, NCOOR
         DO 200 J = 1, I
            HESM = HESMOL(I,J)
            HEST = HESTR (I,J)
            DIFF = HESM - HEST
            HESTR(I,J) = DIFF
            ABSDIF = ABS(DIFF)
            IF (ABSDIF .GE. DIFMAX) THEN
               DIFMAX = ABSDIF
               IDIFMX = I
               JDIFMX = J
            END IF
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 2 .OR. DIFMAX .GE. 1.0D-06) THEN
         IF (DIFMAX .GE. 1.0D-06) THEN
            WRITE (LUPRI,'(///7X,A)') ' >>> WARNING from TRAROT:'
            NWNABA = NWNABA + 1
         END IF
         CALL AROUND('Diff. between Hessian calculated directly'//
     *      ' and using tra.-rot. symmetry')
         CALL PRIHES(HESTR,'CENTERS',CSTRA,SCTRA)
      END IF
      WRITE (LUPRI,'(7X,2A,2(I2,A),1P,E8.2)')
     *  ' Largest difference in',
     *  ' Hessian for element (',IDIFMX,',',JDIFMX,'): ',
     *    DIFMAX
      RETURN
      END
C  /* Deck cmpdip */
      SUBROUTINE CMPDIP(DIPTR,IPRINT,CSTRA,SCTRA)
C
C     Compares the elements of the dipole gradient calculated by
C     translational and rotational symmetry with those calculated
C     directly.
C
C     tuh 171289
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION DIPTR(3,MXCOOR), CSTRA(*), SCTRA(*)
#include <abainf.h>
#include <moldip.h>
#include <trkoor.h>
      DIFMAX = D0
      DO 100 I = 1, 3
         DO 200 J = 1, NCOOR
            DIPF = DIP1(I,J)
            DIPT = DIPTR (I,J)
            DIFF = DIPF - DIPT
            DIPTR(I,J) = DIFF
            ABSDIF = ABS(DIFF)
            IF (ABSDIF .GE. DIFMAX) THEN
               DIFMAX = ABSDIF
               IDIFMX = I
               JDIFMX = J
            END IF
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 2 .OR. DIFMAX .GE. 1.0D-06) THEN
         IF (DIFMAX .GE. 1.0D-06) THEN
            WRITE (LUPRI,'(///7X,A)') ' >>> WARNING from TRAROT:'
            NWNABA = NWNABA + 1
         END IF
         CALL AROUND('Diff. between APTs calculated'//
     *      ' directly and using tra.-rot. symmetry')
         CALL FCPRI(DIPTR,'APT',CSTRA,SCTRA)
      END IF
      WRITE (LUPRI,'(7X,2A,2(I2,A),1P,E8.2)')
     *  ' Largest difference in',
     *  ' dipole gradient for element (',IDIFMX,',',JDIFMX,'): ',
     *    DIFMAX
      RETURN
      END
C  /* Deck gettro */
      SUBROUTINE GETTRO(AMAT,IPTCOL,NROW,NCOL,KEY1,KEY2,KEY3,IREP,
     *                  IPRINT)
C
C     Constructs partial derivatives with respect to tra.-rot.
C     coordinates, in particular
C
C     - partial derivatives of Cartesian coordinates ('TMAT')
C     - partial derivatives of molecular gradient    ('HESS')
C     - partial derivatives of dipole moment         ('DIPOLE')
C
C     tuh 071288
C     tuh 171289 - modified to include dipole moment derivatives
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
      LOGICAL DEPMAT, BOTH, DOTMAT, DOHMAT, DOORTH, TC, DODMAT
      CHARACTER*(*) KEY1, KEY2, KEY3
      DIMENSION AMAT(NROW,NCOL), IPTCOL(NCOOR)
C
#include <taymol.h>
#include <moldip.h>
#include <nuclei.h>
#include <trkoor.h>
#include <symmet.h>
C
      IF (IPRINT.GE.5) THEN
         CALL TITLER('Output from GETTRO','*',103)
         WRITE (LUPRI,'(/A,2I5)') ' NROW, NCOL ', NROW, NCOL
         WRITE (LUPRI,'(2A)') ' KEY1 ', KEY1
         WRITE (LUPRI,'(2A)') ' KEY2 ', KEY2
         WRITE (LUPRI,'(2A)') ' KEY3 ', KEY3
         WRITE (LUPRI,'(A,I2)') ' Symmetry ',IREP + 1
      END IF
C
C     Process key words
C
      IF (KEY1 .EQ. 'DEP') THEN
         DEPMAT = .TRUE.
         BOTH   = .FALSE.
      ELSE IF (KEY1 .EQ. 'IND') THEN
         DEPMAT = .FALSE.
         BOTH   = .FALSE.
      ELSE IF (KEY1 .EQ. 'BOTH') THEN
         DEPMAT = .FALSE.
         BOTH   = .TRUE.
      ELSE
         WRITE (LUPRI,'(/2A)') ' First keyword illegal in GETTRO:', KEY1
         CALL QUIT('Illegal keyword in GETTRO.')
      END IF
C
      IF (KEY2 .EQ. 'TMAT') THEN
         DOTMAT = .TRUE.
         DOHMAT = .FALSE.
         DOORTH = .FALSE.
         DODMAT = .FALSE.
      ELSE IF (KEY2 .EQ. 'TORTHO') THEN
         DOTMAT = .TRUE.
         DOHMAT = .FALSE.
         DOORTH = .TRUE.
         DODMAT = .FALSE.
      ELSE IF (KEY2 .EQ. 'HESS') THEN
         DOTMAT = .FALSE.
         DOHMAT = .TRUE.
         DOORTH = .FALSE.
         DODMAT = .FALSE.
      ELSE IF (KEY2 .EQ. 'DIPOLE') THEN
         DOTMAT = .FALSE.
         DOHMAT = .FALSE.
         DOORTH = .FALSE.
         DODMAT = .TRUE.
      ELSE
         WRITE (LUPRI,'(/2A)') ' Second keyword illegal in GETTRO:',KEY2
         CALL QUIT('Illegal second keyword in GETTRO.')
      END IF
C
      IF (KEY3 .EQ. 'TC') THEN
         TC = .TRUE.
         NTR = NROW
      ELSE IF (KEY3 .EQ. 'CT') THEN
         TC = .FALSE.
         NTR = NCOL
      ELSE
         WRITE (LUPRI,'(/2A)') ' Third keyword illegal in GETTRO:',KEY3
         CALL QUIT('Illegal third keyword in GETTRO.')
      END IF
C
C     Set IPTCOL
C
      IF (DODMAT) THEN
         ICOL = 0
         DO 100 ICOOR = 1, 3
            IF (ISYMAX(ICOOR,1) .EQ. IREP) THEN
               ICOL = ICOL + 1
               IPTCOL(IPTAX(ICOOR,1)) = ICOL
            END IF
  100    CONTINUE
      ELSE
         ICOL = 0
         DO 110 IATOM = 1, NUCIND
            DO 120 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  IF (BOTH .OR. (DEPEND(ISCOOR) .EQV. DEPMAT)) THEN
                     ICOL = ICOL + 1
                     IPTCOL(ISCOOR) = ICOL
                  END IF
               END IF
  120       CONTINUE
  110    CONTINUE
      END IF
      IF (IPRINT .GE. 15) THEN
         WRITE (LUPRI,'(A,L5)') ' DEPMAT ', DEPMAT
         WRITE (LUPRI,'(A,L5)') ' BOTH   ', BOTH
         WRITE (LUPRI,'(A,L5)') ' DOTMAT ', DOTMAT
         WRITE (LUPRI,'(A,L5)') ' DOHMAT ', DOHMAT
         CALL HEADER('IPTCOL',1)
         WRITE (LUPRI,'(12I5)') (IPTCOL(I),I=1,NCOOR)
      END IF
C
C     Construct matrix
C
      CALL DZERO(AMAT,NROW*NCOL)
      IF (DODMAT) THEN
         MAXATM = 1
         CALL MOLCHR(ICHRGE)
      ELSE
         MAXATM = NUCIND
      END IF
      DO 200 ITR = 1, NTR
         ITRO = IPTTRO(ITR,IREP+1)
         DO 300 IATOM = 1, MAXATM
            IF (DOORTH) THEN
               FACTOR = SQRT(FMULT(ISTBNU(IATOM)))
            ELSE
               FACTOR = D1
            END IF
C
C           Translation
C
            IF (ITRO .LT. 4) THEN
               IF (DOTMAT) THEN
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ITRO,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     IF (BOTH .OR. (DEPEND(ISCOOR) .EQV. DEPMAT)) THEN
                        IF (TC) THEN
                           AMAT(ITR,IPTCOL(ISCOOR)) = FACTOR
                        ELSE
                           AMAT(IPTCOL(ISCOOR),ITR) = FACTOR
                        END IF
                     END IF
                  END IF
               ELSE IF (DODMAT) THEN
                  IF (ISYMAX(ITRO,1) .EQ. IREP) THEN
                     IF (TC) THEN
                        AMAT(ITR,IPTCOL(IPTAX(ITRO,1))) = ICHRGE
                     ELSE
                        AMAT(IPTCOL(IPTAX(ITRO,1)),ITR) = ICHRGE
                     END IF
                  END IF
               END IF
C
C           Rotation
C
            ELSE
               IA   = MOD(ITRO - 2,3) + 1
               IB   = MOD(ITRO - 3,3) + 1
               IF (DODMAT) THEN
                  IF (ISYMAX(IA,1) .EQ. IREP) THEN
                     IF (TC) THEN
                       AMAT(ITR,IPTCOL(IPTAX(IA,1))) = DIP0(IB)
                     ELSE
                       AMAT(IPTCOL(IPTAX(IA,1)),ITR) = DIP0(IB)
                     END IF
                  END IF
                  IF (ISYMAX(IB,1) .EQ. IREP) THEN
                     IF (TC) THEN
                        AMAT(ITR,IPTCOL(IPTAX(IB,1))) = - DIP0(IA)
                     ELSE
                        AMAT(IPTCOL(IPTAX(IB,1)),ITR) = - DIP0(IA)
                     END IF
                  END IF
               ELSE
                  IF (DOTMAT) THEN
                     FACA = FACTOR*CORD(IA,IATOM)
                     FACB = FACTOR*CORD(IB,IATOM)
                  ELSE
                     ISCORA = IPTCNT(3*(IATOM - 1) + IA,0,1)
                     ISCORB = IPTCNT(3*(IATOM - 1) + IB,0,1)
                     IF (ISCORA .GT. 0) THEN
                        FACA = FACTOR*GRDMOL(ISCORA)
                     ELSE
                        FACA = 0.0D0
                     END IF
                     IF (ISCORB .GT. 0) THEN
                        FACB = FACTOR*GRDMOL(ISCORB)
                     ELSE
                        FACB = 0.0D0
                     END IF
                  END IF
                  ISCORA = IPTCNT(3*(IATOM - 1) + IA,IREP,1)
                  ISCORB = IPTCNT(3*(IATOM - 1) + IB,IREP,1)
                  IF (ISCORA .GT. 0) THEN
                     IF (BOTH .OR. (DEPEND(ISCORA).EQV.DEPMAT)) THEN
                        IF (TC) THEN
                          AMAT(ITR,IPTCOL(ISCORA)) = FACB
                        ELSE
                          AMAT(IPTCOL(ISCORA),ITR) = FACB
                        END IF
                     END IF
                  END IF
                  IF (ISCORB .GT. 0) THEN
                     IF (BOTH .OR. (DEPEND(ISCORB).EQV.DEPMAT)) THEN
                        IF (TC) THEN
                           AMAT(ITR,IPTCOL(ISCORB)) = - FACA
                        ELSE
                           AMAT(IPTCOL(ISCORB),ITR) = - FACA
                        END IF
                     END IF
                  END IF
               END IF
            END IF
  300    CONTINUE
  200 CONTINUE
      IF (IPRINT .GE. 15) THEN
         CALL AROUND('AMAT - GETTRO')
         CALL OUTPUT(AMAT,1,NROW,1,NCOL,NROW,NCOL,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck vectra */
      SUBROUTINE VECTRA(VECIND,VECDEP)
C
C     This subroutine determines the dependent differentiated gradients
C     or CP SCF solutions using translational symmetry.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.D00)
C
      LOGICAL OLDDX, INDCOR
      DIMENSION VECIND(*), VECDEP(*)
C
C Used from common blocks:
C  NUCLEI: DCORGD, DCORD
C
#include <abainf.h>
#include <cbitro.h>
#include <nuclei.h>
#include <dorps.h>
#include <infvar.h>
#include <inftap.h>
C
#include <chrxyz.h>
C
C     Run over vector types (gradient or solution vector)
C
      DO 100 ITYPE = 1, 2
         IF (ITYPE .EQ. 1) THEN
            IF (.NOT. GDTRO) GO TO 100
            LUX = LUGDR
            CALL GPOPEN(LUX,ABAGDR,'OLD',' ','DIRECT',IRAT*NVAR,OLDDX)
            IF (COMPAR .OR. IPRINT .GE. 5) THEN
               CALL HEADER('Translational symmetry for GD vectors',1)
            END IF
         ELSE
            IF (.NOT. RDTRO) GO TO 100
            LUX = LURDR
            CALL GPOPEN(LUX,ABARDR,'OLD',' ','DIRECT',IRAT*NVAR,OLDDX)
            IF (COMPAR .OR. IPRINT .GE. 5) THEN
               CALL HEADER('Translational symmetry for RD vectors',1)
            END IF
         END IF
C
C        Run over coordinates
C
         DO 200 ICOOR = 1, 3
            IF (COMPAR .OR. IPRINT .GE. 5) THEN
               CALL HEADER('Translational symmetry for coordinate '//
     *                     CHRXYZ(ICOOR),1)
            END IF
C
C           Calculate dependent vector
C
            CALL DZERO(VECDEP,NVAR)
            DO 300 IATOM = 1, NUCIND
               IF (ITYPE .EQ. 1) THEN
                  INDCOR = DCORGD(IATOM,ICOOR,1)
               ELSE
                  INDCOR = DCORD(IATOM,ICOOR,1)
               END IF
               IF (IPRINT .GE. 5) THEN
                  WRITE (LUPRI,'(/A,I5)') ' ATOM NUMBER: ', IATOM
                  IF (INDCOR) THEN
                     WRITE (LUPRI,'(/A)') ' INDEPENDENT ATOM '
                  ELSE
                     WRITE (LUPRI,'(/A)') ' DEPENDENT ATOM '
                  END IF
               END IF
               IF (INDCOR) THEN
                  ICRIND = 3*(IATOM - 1) + ICOOR
                  CALL READDX (LUX,ICRIND,IRAT*NVAR,VECIND)
                  IF (IPRINT.GE.20) THEN
                     CALL HEADER('VECIND',-1)
                     CALL OUTPUT(VECIND,1,1,1,NVAR,1,NVAR,1,LUPRI)
                  END IF
                  DO 400 I = 1, NVAR
                     VECDEP(I) = VECDEP(I) - VECIND(I)
  400             CONTINUE
               ELSE
                  ICRDEP = 3*(IATOM - 1) + ICOOR
               END IF
  300       CONTINUE
C
C           Compare (check to see that sum is zero)
C
            IF (COMPAR) THEN
               DIFMAX = ZERO
               DIF2 = ZERO
               DO 500 I = 1, NVAR
                  DIFFER = ABS(VECDEP(I))
                  DIF2 = DIF2 + DIFFER*DIFFER
                  IF (DIFFER .GT. DIFMAX) THEN
                     DIFMAX = DIFFER
                     MAXELE = I
                  END IF
  500          CONTINUE
               DIF2 = SQRT(DIF2)
               IF (DIF2 .GE. 1.0D-10) THEN
                  WRITE (LUPRI,'(///7X,A)') ' >>> WARNING from VECTRA:'
                  NWNABA = NWNABA + 1
               END IF
               WRITE (LUPRI,'(/,A,D12.6)') ' Square root deviation '//
     *            ' between vectors: ',  DIF2
               WRITE (LUPRI,'(/,A,D12.6,A,I5)') ' Maximum difference ',
     *                           DIFMAX, ' for element ', MAXELE
C
C           Write calculated vector
C
            ELSE
               IF (IPRINT .GE. 15) THEN
                  CALL HEADER('VECDEP',-1)
                  CALL OUTPUT(VECDEP(1),1,1,1,NVAR,1,NVAR,1,LUPRI)
               END IF
               CALL WRITDX (LUX,ICRDEP,IRAT*NVAR,VECDEP)
            END IF
  200    CONTINUE
  100 CONTINUE
      IF (GDTRO) CALL GPCLOSE(LUGDR,'KEEP')
      IF (RDTRO) CALL GPCLOSE(LURDR,'KEEP')
      RETURN
C     END OF VECTRA
      END
C  /* Deck gettr */
      SUBROUTINE GETTR(TMAT,AMAT,TVEC)
C
C     This subroutine constructs the full matrix of derivatives of
C     Cartesian coordinates with respect to tra.-rot. coordinates.
C     Independent translations and rotations are identified.
C
C     tuh 071288
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.D00)
C
      LOGICAL OK, DONE
      DIMENSION TMAT(6,NCOOR), TVEC(NCOOR), AMAT(NCOOR,6)
C
#include <cbitro.h>
#include <nuclei.h>
#include <trkoor.h>
#include <symmet.h>
#include <chrxyz.h>
#include <ibtfun.h>
C
C     (I) Set up TMAT for 6 trarot coordinates
C
      IF (IPRINT .GE. 2)
     *   CALL HEADER('Translational and rotational coordinates',1)
      NCDEP = 0
      ITR = 0
      ICR = 1
      DO 100 IREP = 0, MAXREP
         NTR = 0
         NCR = NCRREP(IREP,1)
         DONE = .FALSE.
         DO 200 IA = 1, 3
C
C           Translation
C
            IF (ISYMAX(IA,1) .EQ. IREP) THEN
               IF (IPRINT .GE. 2) THEN
                  IF (.NOT.DONE) WRITE (LUPRI,'(/1X,A,I2/)')
     *              ' Symmetry ', IREP + 1
                  DONE = .TRUE.
                  WRITE (LUPRI,'(7X,3A)')
     *               ' Translation along ',CHRXYZ(IA),' axis'
               END IF
               CALL DZERO(TVEC,NCOOR)
               DO 300 IATOM = 1, NUCIND
                  ISCOOR = IPTCNT(3*(IATOM - 1) + IA,IREP,1)
                  IF (ISCOOR .GT. 0) TVEC(ISCOOR) = D1
  300          CONTINUE
               CALL DCOPY(NCOOR,TVEC,1,AMAT(1,NTR+1),1)
               CALL LINDEP(AMAT,THRESH,OK,NCOOR,6,NTR,1,IPRINT)
               IF (OK) THEN
                  ITR = ITR + 1
                  NTR = NTR + 1
                  IPTTRO(NTR,IREP + 1) = IA
                  CALL DCOPY(NCOOR,TVEC,1,TMAT(ITR,1),6)
               END If
            END IF
            DO 400 IB = 1, IA - 1
C
C              Rotation
C
               IF (IBTXOR(ISYMAX(IA,1),ISYMAX(IB,1)).EQ.IREP) THEN
                  IC = IBTXOR(IA,IB)
                  IF (IPRINT .GE. 2) THEN
                     IF (.NOT.DONE) WRITE (LUPRI,'(/1X,A,I2/)')
     *              ' Symmetry ', IREP + 1
                     DONE = .TRUE.
                     WRITE (LUPRI,'(7X,3A)')
     *                  ' Rotation around ',CHRXYZ(IC),' axis'
                  END IF
                  CALL DZERO(TVEC,NCOOR)
                  DO 500 IATOM = 1, NUCIND
                     ISCOOR = IPTCNT(3*(IATOM - 1) + IA,IREP,1)
                     IF (ISCOOR.GT.0) TVEC(ISCOOR) = CORD(IB,IATOM)
                     ISCOOR = IPTCNT(3*(IATOM - 1) + IB,IREP,1)
                     IF (ISCOOR.GT.0) TVEC(ISCOOR) = - CORD(IA,IATOM)
  500             CONTINUE
                  CALL DCOPY(NCOOR,TVEC,1,AMAT(1,NTR+1),1)
                  CALL LINDEP(AMAT,THRESH,OK,NCOOR,6,NTR,1,IPRINT)
                  IF (OK) THEN
                     ITR = ITR + 1
                     NTR = NTR + 1
                     IPTTRO(NTR,IREP + 1) = 3 + IC
                     CALL DCOPY(NCOOR,TVEC,1,TMAT(ITR,1),6)
                  ELSE
                     IF (IPRINT .GE. 2) THEN
                        WRITE (LUPRI,'(7X,A)') '      - not independent'
                     END IF
                  END If
               END IF
  400       CONTINUE
  200    CONTINUE
         NTRREP(IREP) = NTR
         NCDEP = NCDEP + NTR
  100 CONTINUE
      IF (IPRINT .GE. 15) THEN
         CALL AROUND('T matrix - GETTR')
         IOFFTR = 0
         IOFFCR = 0
         DO 600 IREP = 0, MAXREP
            WRITE (LUPRI,'(//A,I2)') '  Symmetry ', IREP + 1
            NTR = NTRREP(IREP)
            NCR = NCRREP(IREP,1)
            CALL OUTPUT(TMAT,IOFFTR + 1,IOFFTR + NTR,
     *                       IOFFCR + 1,IOFFCR + NCR,
     *                       6,NCOOR,1,LUPRI)
            IOFFTR = IOFFTR + NTR
            IOFFCR = IOFFCR + NCR
  600    CONTINUE
         CALL AROUND('IPTTRO - GETTR')
         DO 610 ITR = 1, 6
            WRITE(LUPRI,'(4X,8I4)') (IPTTRO(ITR,ISYM),ISYM=1,MAXREP+1)
  610    CONTINUE
      END IF
      RETURN
      END
C  /* Deck lindep */
      SUBROUTINE LINDEP(AMAT,THRESH,OK,NROW,NCOL,NOLD,NNEW,IPRINT)
C
C     Checks input vector for linear dependency against old vectors
C     and adds input vector to old ones if independent.
C
C     tuh Nov 1988
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL OK
      DIMENSION AMAT(NROW,NCOL)
C
      IF (NOLD + NNEW .GT. NCOL) THEN
         OK = .FALSE.
      ELSE
         OK = .TRUE.
         DO 100 INEW = NOLD + 1, NOLD + NNEW
            DO 200 IOLD = 1, INEW - 1
               COEF = - DDOT(NROW,AMAT(1,IOLD),1,AMAT(1,INEW),1)
               CALL DAXPY(NROW,COEF,AMAT(1,IOLD),1,AMAT(1,INEW),1)
  200       CONTINUE
            PROD = DDOT(NROW,AMAT(1,INEW),1,AMAT(1,INEW),1)
            PROD = SQRT(PROD)
            IF (PROD .GT. THRESH) THEN
               FACTOR = D1/PROD
            ELSE
               FACTOR = D0
            END IF
            CALL DSCAL(NROW,FACTOR,AMAT(1,INEW),1)
            OK = OK .AND. PROD .GT. THRESH
  100    CONTINUE
      END IF
      IF (IPRINT .GE. 15) THEN
         CALL TITLER('Output from LINDEP','*',103)
         WRITE (LUPRI,'(A,2I5)') ' NROW, NCOL  ', NROW, NCOL
         WRITE (LUPRI,'(A,2I5)') ' NOLD, NNEW  ', NOLD, NNEW
         WRITE (LUPRI,'(A,L5)') ' OK ', OK
         CALL AROUND('AMAT - LINDEP')
         CALL OUTPUT(AMAT,1,NROW,1,NOLD+NNEW,NROW,NCOL,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck setgd */
      SUBROUTINE SETGD
C
C     Set DCORGD, indicating which GD vectors to calculate
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
C Used from common blocks:
C  ABAINF: GDALL
C  CBITRO: COMPAR
C  NUCLEI: DCORD, DCORGD, NTRACO, ITRACO
C
#include <abainf.h>
#include <cbitro.h>
#include <nuclei.h>
#include <dorps.h>
C
C     DCORGD is initialized to all .true. in HESINP subroutine/
C
      IF (GDALL .AND. .NOT. COMPAR) THEN
         DO 100 J = 1,NTRACO
            ICOOR = ITRACO(J)
            INUC  = (ICOOR-1)/3 + 1
            ICOOR = ICOOR - (INUC-1)*3
            DCORGD(INUC,ICOOR,1) = .FALSE.
  100    CONTINUE
      ELSE
         DO 300 I = 1, MXCENT
            DO 200 J = 1, 3
               DCORGD(I,J,1) = DCORD(I,J,1)
  200       CONTINUE
  300    CONTINUE
      END IF
      RETURN
      END
C  /* Deck setdcr */
      SUBROUTINE SETDCR(KEY)
C
C     tuh Mar 08 90
C
C     This subroutine sets up array logicals DOREPS(0:7) and
C     DOCOOR(3,MXCENT).
C
C     DOCOOR(IXYZ,IATOM) is false if the corresponding symmetry independent
C     atom IATOM and Cartesian direction IXYX does not contribute to any
C     symmetry coordinate.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
C
      CHARACTER*(*) KEY
C
#include <abainf.h>
#include <cbitro.h>
#include <nuclei.h>
#include <symmet.h>
#include <dorps.h>
#include <trkoor.h>
C
      IF (KEY .EQ. 'HERMIT') THEN
         DO 100 IREP = 0, MAXREP
            DOREPS(IREP) = .TRUE.
  100    CONTINUE
      ELSE
         DO 200 IREP = 0, MAXREP
            IF (COMPAR) THEN
               DOREPS(IREP) = DOSYM(IREP+1)
            ELSE
               DOREPS(IREP) = DOSYM(IREP+1) .AND.
     &          (((MOLGRD.OR.MOLHES).AND.
     &           (NCRREP(IREP,1).GT.NTRREP(IREP)))
     &            .OR. ((DIPDER.OR.POLAR) .AND.(NTRREP(IREP).GT.0)))
            END IF
  200    CONTINUE
         ICOOR = 0
         DO 300 IATOM = 1, NUCIND
            DO 400 IXYZ = 1, 3
               ICOOR = ICOOR + 1
               DOCOOR(IXYZ,IATOM) = .FALSE.
               DO 500 IREP = 0, MAXREP
                  IF (DOREPS(IREP).AND.(IPTCNT(ICOOR,IREP,1).GT.0)) THEN
                     DOCOOR(IXYZ,IATOM) = .TRUE.
                  END IF
  500          CONTINUE
  400       CONTINUE
  300    CONTINUE
         IF (IPRINT .GE. 10) THEN
           CALL AROUND('DOREPS in SETDCR')
           WRITE (LUPRI,'(1X,8L5)') (DOREPS(I),I=0,MAXREP)
           CALL AROUND
     &      ('Contributing sym. unique coordinates (DOCOOR in SETDCR)')
           CALL HEADER('      atom        x         y         z     ',1)
           DO 600 IATOM = 1, NUCIND
              WRITE (LUPRI,'(1X,I10,3L10)')IATOM,(DOCOOR(I,IATOM),I=1,3)
  600      CONTINUE
         END IF
      END IF
      RETURN
      END
