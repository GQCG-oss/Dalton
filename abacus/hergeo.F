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
#ifdef REV_LOG
c===========================================================================
c942508-kr: As P.Taylor pointed out, IPRINT never used, and thus removed from
cCBIANA.
c===========================================================================
#endif
C  /* Deck anainp */
      SUBROUTINE ANAINP(WORD)
C
C  5-Jul-1985 Hans Jorgen Aa. Jensen
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 4)
      PARAMETER (MAXANG = 20)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
      LOGICAL SKIP
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
      DATA TABLE /'.SKIP  ', '.XXXXXX', '.ANGLES', '.DIHEDR'/
      DATA MANG/0/, MDIHED/0/
C
      CALL ANAINI
C
      NEWDEF = (WORD .EQ. '*GEOANA')
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
                     GO TO (1,2,3,4), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ANAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ANAINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) NANG
                  MANG = MIN(MAXANG,NANG)
                  DO 310 I = 1,MANG
                     READ(LUCMD,*) (IANG(J,I),J=1,3)
  310             CONTINUE
                  MANG = NANG - MANG
                  DO 320 I = 1,MANG
                     READ(LUCMD,'()')
  320             CONTINUE
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) NDIHED
                  MDIHED = MIN(MAXANG,NDIHED)
                  DO 410 I = 1,MDIHED
                     READ(LUCMD,*) (IDIHED(J,I),J=1,4)
  410             CONTINUE
                  MDIHED = NDIHED - MDIHED
                  DO 420 I = 1,MDIHED
                     READ(LUCMD,'()')
  420             CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ANAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ANAINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for GEOANA:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' GEOANA skipped in this run.'
         ELSE
            IF (NANG .GT. 0) THEN
               WRITE (LUPRI,'(/A/)')
     *            ' Following angles will be calculated:'
               DO 1310 I = 1,NANG
                  WRITE (LUPRI,'(I10,A,4I5)') I,' : ',(IANG(J,I),J=1,3)
 1310          CONTINUE
               IF (MANG .GT. 0) THEN
                  WRITE (LUPRI,'(/A,I3,A)') ' The last',MANG,
     *               ' angles specified go beyond current maximum',
     *               ' and will not be printed.'
               END IF
            END IF
            IF (NDIHED .GT. 0) THEN
               WRITE (LUPRI,'(/A/)')
     *            ' Following dihedral angles will be calculated:'
               DO 1410 I = 1,NDIHED
                  WRITE (LUPRI,'(I10,A,4I5)')I,' : ',(IDIHED(J,I),J=1,4)
 1410          CONTINUE
               IF (MDIHED .GT. 0) THEN
                  WRITE (LUPRI,'(/A,I3,A)') ' The last',MDIHED,
     *               ' dihedral angles specified go beyond current',
     *               ' maximum and will not be printed.'
               END IF
            END IF
         END IF
         WRITE (LUPRI,'(/)')
      END IF
      RETURN
      END
C  /* Deck anaini */
      SUBROUTINE ANAINI
C
C     Initialize /CBIANA/
C
#include <implicit.h>
      PARAMETER (MAXANG = 20)
      LOGICAL SKIP
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
      NANG   = 0
      NDIHED = 0
      SKIP   = .FALSE.
      RETURN
      END
C  /* Deck geoana */
      SUBROUTINE GEOANA(COORD,PRINT,DIF,NBONDS,PUNCH,LUIP,WORK,LWORK)
C
C 30-Jun-1985 Hans Jorgen Aa. Jensen
C Modified for symmetry 25-Sep-1989 tuh
C Modified for differential geomtries 18-Oct-1989 tuh
C
#include <implicit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      LOGICAL PRINT, DIF, PUNCH
      DIMENSION WORK(LWORK)
#include <nuclei.h>
C
      CALL QENTER('GEOANA')
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KVEC,3*NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDIST2,NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KBOND ,NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KCHRG ,NUCDEP,       WORK,KFREE,LFREE)
      CALL MEMGET('INTE',IPAIR,2*NUCDEP*NUCDEP,WORK,KFREE,LFREE)
C
      CALL GEOAN1(COORD,PRINT,DIF,NBONDS,PUNCH,LUIP,WORK(KVEC),
     &            WORK(KDIST2),WORK(KBOND),WORK(KCHRG),WORK(IPAIR))
C
      CALL MEMREL('GEOANA',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('GEOANA')
      RETURN
      END
C  /* Deck geoan1 */
      SUBROUTINE GEOAN1(COORD,PRINT,DIF,NBONDS,PUNCH,LUIP,VEC,DIST2,
     &                  BONDED,ICHARG,IPAIRS)
C
C Modified for more selective printing of bonded atoms,
C     Jan-1995 Hanne Heiberg
C     Added punching atom bonds for Gamess graphic output, K.Ruud-95
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <facang.h>
      PARAMETER (MAXANG = 20)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL SKIP, DIF, PRINT, PUNCH
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
#include <gnrinf.h>
#include <inftap.h>
#include <qm3.h>
#include <nuclei.h>
#include <symmet.h>
C
      DIMENSION COORD(3,*),
     *          DIST(MXCENT*(MXCENT+1)/2),
     *          ANGLE(MAXANG), DIHED(MAXANG), ICHARG(NUCDEP)
      DIMENSION VEC(3,NUCDEP,NUCDEP), DIST2(NUCDEP,NUCDEP), 
     &          IPAIRS(2,NUCDEP*NUCDEP)
      LOGICAL   BONDED(NUCDEP,NUCDEP)
      CHARACTER*6 NUCNAM(4)
      SAVE DIST, ANGLE, DIHED
      ARCCOS(ARG) = FACANG*ACOS(ARG)
#include <ibtfun.h>
C
      NBONDS = 0
      IF (NUCDEP .EQ. 1) RETURN
C
C     set up bond vectors
C
      IATOMA = 0
      DO 100 ICENTA = 1, NUCIND
         DO 110 IA = 0, MAXOPR
            IF (IBTAND(IA,ISTBNU(ICENTA)) .EQ. 0) THEN
               IATOMA = IATOMA + 1
               CXA = PT(IBTAND(ISYMAX(1,1),IA))*COORD(1,ICENTA)
               CYA = PT(IBTAND(ISYMAX(2,1),IA))*COORD(2,ICENTA)
               CZA = PT(IBTAND(ISYMAX(3,1),IA))*COORD(3,ICENTA)
C
               ICHARG(IATOMA) = IZATOM(ICENTA)
C
               IATOMB = 0
               DO 200 ICENTB = 1, NUCIND
                  DO 210 IB = 0, MAXOPR
                     IF (IBTAND(IB,ISTBNU(ICENTB)) .EQ. 0) THEN
                        IATOMB = IATOMB + 1
                        IF (IATOMB .GT. IATOMA) GO TO 110
C                       ... next IATOMA, only IATOMB .le. IATOMA needed
                          CXB=PT(IBTAND(ISYMAX(1,1),IB))*COORD(1,ICENTB)
                          CYB=PT(IBTAND(ISYMAX(2,1),IB))*COORD(2,ICENTB)
                          CZB=PT(IBTAND(ISYMAX(3,1),IB))*COORD(3,ICENTB)
                          VEC(1,IATOMB,IATOMA) = XTANG*(CXA - CXB)
                          VEC(2,IATOMB,IATOMA) = XTANG*(CYA - CYB)
                          VEC(3,IATOMB,IATOMA) = XTANG*(CZA - CZB)
                          VEC(1,IATOMA,IATOMB) = -VEC(1,IATOMB,IATOMA)
                          VEC(2,IATOMA,IATOMB) = -VEC(2,IATOMB,IATOMA)
                          VEC(3,IATOMA,IATOMB) = -VEC(3,IATOMB,IATOMA)
                     END IF
  210             CONTINUE
  200          CONTINUE
            END IF
  110    CONTINUE
  100 CONTINUE
C
C     set up distance matrix
C
      DISTMX = D0
      IDISTMX = 0
      JDISTMX = 0
      IJ = 0
      DO 400 I = 1,NUCDEP
         DO 300 J = 1,I
            IJ = IJ + 1
            DISTAN = VEC(1,J,I)*VEC(1,J,I) + VEC(2,J,I)*VEC(2,J,I)
     *             + VEC(3,J,I)*VEC(3,J,I)
            DISTAN = SQRT(DISTAN)
            DIST2(I,J) = DISTAN
            DIST2(J,I) = DISTAN
            IF (.NOT.DIF) THEN
               DIST(IJ) = DISTAN
            ELSE
               DIST(IJ) = DISTAN - DIST(IJ)
            END IF
            IF (ABS(DIST(IJ)) .GT. DISTMX) THEN
               DISTMX  = DIST(IJ)
               IDISTMX = I
               JDISTMX = J
            END IF
  300    CONTINUE
  400 CONTINUE
C
      IF (.NOT. QM3) THEN
        IF (PRINT) THEN
          IF (NUCDEP .LE. 12*(IPRUSR+1) .AND. DISTMX .NE. D0) THEN
C hjaaj Oct 2003: this output is only useful for small molecules ...
            IF (DIF) THEN
              CALL HEADER
     *        ('Differential interatomic separations (in Angstroms):',2)
            ELSE
              CALL HEADER('Interatomic separations (in Angstroms):',2)
            END IF
            CALL PRIDIS(NAMDEP,DIST,NUCDEP)
          END IF
          IF (DISTMX .NE. D0) THEN
          IF (DIF) THEN
            WRITE(LUPRI,'(/A,F12.6,A/5A)')
     &         '  Max differential change in interatomic separation is',
     &         DISTMX,' Angstroms','  between atoms "',
     &         NAMDEP(IDISTMX),'" and "',NAMDEP(JDISTMX),'".'
          ELSE
            WRITE(LUPRI,'(/A,F10.4,A/5A)')
     &         '  Max interatomic separation is',
     &         DISTMX,' Angstroms','  between atoms "',
     &         NAMDEP(IDISTMX),'" and "',NAMDEP(JDISTMX),'".'
          END IF
          END IF
        END IF
      END IF
C
      IF (PRINT .AND. .NOT. DIF) THEN

        DO 10 J= 1,NUCDEP
          RADJ = RADIUS(ICHARG(J))
          IF (RADJ .LT. 0.0D0 .AND. ICHARG(J) .GT. 0) THEN
C             do not print if cavity center or floating orbital /hjaaj
              WRITE(LUPRI,*)
     &        'RADIUS FOR ATOM WITH ATOMIC NUMBER ',
     &        ICHARG(J),' IS UNAVAILABLE, USING 2.0 AA'
            RADJ = 2.0D0
          END IF
          DO 20 I= 1, J-1
            RADI = RADIUS(ICHARG(I))
            IF (RADI .LT. 0.0D0 .AND. ICHARG(I) .GT. 0) THEN
              RADI = 2.0D0
            END IF
            IF (RADI .LT. 0 .OR. RADJ .LT. 0) THEN
C             not bonded if cavity center or floating orbital /hjaaj
              BONDED(I,J) = .FALSE.
              BONDED(J,I) = .FALSE.
            ELSE IF (DIST2(I,J) .LT. (1.2D0 * (RADI + RADJ))) THEN
              NBONDS = NBONDS + 1
              IPAIRS(1,NBONDS) = I
              IPAIRS(2,NBONDS) = J
              BONDED(I,J) = .TRUE.
              BONDED(J,I) = .TRUE.
            ELSE
              BONDED(I,J) = .FALSE.
              BONDED(J,I) = .FALSE.
            END IF
 20       CONTINUE
          BONDED(J,J) = .FALSE.
 10     CONTINUE
        IF (.NOT. QM3) THEN
          CALL HEADER('Bond distances (angstroms):',1)
          WRITE (LUPRI,'(18X,A/18X,A)')
     $          'atom 1     atom 2       distance',
     $          '------     ------       --------'
          DO 500 I = 1, NUCDEP
            DO 510 J = 1, I-1
              IF (BONDED(I,J)) THEN
                NUCNAM(1) = NAMDEP(I)
                NUCNAM(2) = NAMDEP(J)
                WRITE(LUPRI,'(A,2X,A6,5X,A6,F15.6)')
     $               '  bond distance:',
     &               NUCNAM(1), NUCNAM(2), DIST2(I,J)
              END IF
 510        CONTINUE
 500      CONTINUE
C
          IF (NUCDEP .GT. 2 .AND. NANG .LE. 0) THEN
            CALL HEADER('Bond angles (degrees):',1)
            WRITE (LUPRI,'(18X,A/18X,A)')
     $            'atom 1     atom 2     atom 3         angle',
     $            '------     ------     ------         -----'
C
            IJK = 0
            DO 40, I= 1,NUCDEP
              DO 50, J= 1, NUCDEP - 1
                DO 60, K= J + 1, NUCDEP
                  IF (BONDED(I,J) .AND. BONDED(I,K)) THEN
                    IJK = IJK + 1
                    NUCNAM(1) = NAMDEP(J)
                    NUCNAM(2) = NAMDEP(I)
                    NUCNAM(3) = NAMDEP(K)
                    ANG = VECANG(VEC(1,I,J),VEC(1,I,K))
                    ANG = ANG*FACANG
                    WRITE(LUPRI,'(A,5X,A6,5X,A6,5X,A6,F14.3)')
     *                    '  bond angle:',NUCNAM(1),NUCNAM(2),
     *                    NUCNAM(3),ANG
                  END IF
 60             CONTINUE
 50           CONTINUE
 40         CONTINUE
            IF (IJK .EQ. 0) WRITE(LUPRI,'(5X,A)') 'No angles found'
          END IF
        END IF
      END IF
C
C     Punch bonding information in Gamess output format
C
      IF (PUNCH) THEN
         IF(NBONDS.LE.6) THEN
            WRITE(LUIP,8010) (IPAIRS(1,I),IPAIRS(2,I),I=1,NBONDS)
         ELSE
            WRITE(LUIP,8020) (IPAIRS(1,I),IPAIRS(2,I),I=1,6)
            WRITE(LUIP,8030) (IPAIRS(1,I),IPAIRS(2,I),I=7,NBONDS)
         END IF
      END IF
C         
      IF (NANG .GT. 0) THEN
         IF (PRINT) THEN
            CALL HEADER('Angles according to input list:',2)
            WRITE (LUPRI,'(A/A)')
     *       '    atom 1     atom 2     atom 3         angle (degrees)',
     *       '    ------     ------     ------         ---------------'
         END IF
         DO 1000 I = 1,NANG
            I1 = IANG(1,I)
            I2 = IANG(2,I)
            I3 = IANG(3,I)
            IMX = MAX(I1,I2,I3)
            IF (IMX .GT. NUCDEP) THEN
               IF (PRINT) WRITE (LUPRI,'(/A/)')
     &            ' *GEOANA input error for .ANGLES: non-existent atom'
               GO TO 1000
            END IF
            NUCNAM(1) = NAMDEP(I1)
            NUCNAM(2) = NAMDEP(I2)
            NUCNAM(3) = NAMDEP(I3)
            IF (I1 .NE. I2 .AND. I2 .NE. I3) THEN
               ANG = VECANG(VEC(1,I2,I1),VEC(1,I2,I3))
               ANG = ANG*FACANG
               IF (.NOT.DIF) THEN
                  ANGLE(I) = ANG
               ELSE
                  ANGLE(I) = ANG - ANGLE(I)
               END IF
               IF (PRINT) WRITE (LUPRI,'(4X,A6,5X,A6,5X,A6,F20.3)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),ANGLE(I)
            ELSE
               IF (PRINT) WRITE (LUPRI,'(4X,A6,5X,A6,5X,A6,10X,A)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),'undefined'
            END IF
 1000    CONTINUE
      END IF
C
      IF (NDIHED .GT. 0) THEN
          IF (PRINT) WRITE (LUPRI,'(//A/A)')
     *       '    atom 1     atom 2     atom 3     atom 4'//
     *       '    dihedral angle (degrees)',
     *       '    ------     ------     ------     ------'//
     *       '    ------------------------'
         DO 2000 I = 1,NDIHED
            I1 = IDIHED(1,I)
            I2 = IDIHED(2,I)
            I3 = IDIHED(3,I)
            I4 = IDIHED(4,I)
            IMX = MAX(I1,I2,I3,I4)
            IF (IMX .GT. NUCDEP) THEN
               IF (PRINT) WRITE (LUPRI,'(/A/)')
     &            ' *GEOANA input error for .DIHEDR: non-existent atom'
               GO TO 2000
            END IF
            NUCNAM(1) = NAMDEP(I1)
            NUCNAM(2) = NAMDEP(I2)
            NUCNAM(3) = NAMDEP(I3)
            NUCNAM(4) = NAMDEP(I4)
            X1 = VEC(2,I2,I1)*VEC(3,I2,I3) - VEC(2,I2,I3)*VEC(3,I2,I1)
            X2 = VEC(3,I2,I1)*VEC(1,I2,I3) - VEC(3,I2,I3)*VEC(1,I2,I1)
            X3 = VEC(1,I2,I1)*VEC(2,I2,I3) - VEC(1,I2,I3)*VEC(2,I2,I1)
            Y1 = VEC(2,I3,I2)*VEC(3,I3,I4) - VEC(2,I3,I4)*VEC(3,I3,I2)
            Y2 = VEC(3,I3,I2)*VEC(1,I3,I4) - VEC(3,I3,I4)*VEC(1,I3,I2)
            Y3 = VEC(1,I3,I2)*VEC(2,I3,I4) - VEC(1,I3,I4)*VEC(2,I3,I2)
            Z1 = X2*Y3 - X3*Y2
            Z2 = X3*Y1 - X1*Y3
            Z3 = X1*Y2 - X2*Y1
            SENSE = Z1*VEC(1,I2,I3) + Z2*VEC(2,I2,I3) + Z3*VEC(3,I2,I3)
            SENSE = SIGN(D1,SENSE)
            ANG = X1*Y1 + X2*Y2 + X3*Y3
            DDD = (X1*X1 + X2*X2 + X3*X3) * (Y1*Y1 + Y2*Y2 + Y3*Y3)
            IF (DDD .GT. 1.D-10) THEN
               ANG = ANG / SQRT(DDD)
               IF (ABS(ANG) .GT. D1) ANG = SIGN(D1,ANG)
               ANG = SENSE*ARCCOS(ANG)
               IF (.NOT.DIF) THEN
                  DIHED(I) = ANG
               ELSE
                  DIHED(I) = ANG - DIHED(I)
               END IF
               IF (PRINT) WRITE(LUPRI,'(4X,A6,5X,A6,5X,A6,5X,A6,F20.3)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),NUCNAM(4),DIHED(I)
            ELSE
               IF (PRINT) WRITE(LUPRI,'(4X,A6,5X,A6,5X,A6,5X,A6,10X,A)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),NUCNAM(4),'undefined'
            END IF
 2000    CONTINUE
      END IF
C
      IF (PRINT) WRITE (LUPRI,'(/)')
      RETURN
 8010 FORMAT('BONDATOMS ',6(I4,I4,2X))
 8020 FORMAT('BONDATOMS ',6(I4,I4,2X),' >')
 8030 FORMAT(7(I4,I4,2X),:,' >')
      END
C  /* Deck pridis */
      SUBROUTINE PRIDIS (NAMDEP,DISMAT,NROW)
C
C 30-Jun-1985 Hans Jorgen Aa. Jensen
C (based on OUTPAK by Nelson H.F. Beebe)
C
C Print bond distance matrix (or other matrix over atoms)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (KCOL=6)
      CHARACTER*6 NAMDEP(*)
      DIMENSION DISMAT(*)
      INTEGER BEGIN
C
      LAST = MIN(NROW,KCOL)
      BEGIN = 1
 1050 NCOL = 1
      WRITE (LUPRI,1000) (NAMDEP(I),I = BEGIN,LAST)
      WRITE (LUPRI,1000) ('------' ,I = BEGIN,LAST)
      DO 40 K = BEGIN,NROW
         KTOTAL = (K*(K-1))/2 + BEGIN - 1
         WRITE (LUPRI,2000) NAMDEP(K),
     *      (DISMAT(KTOTAL+J),J = 1,NCOL)
         IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40 CONTINUE
      WRITE (LUPRI,'()')
      LAST = MIN(LAST+KCOL,NROW)
      BEGIN = BEGIN+NCOL
      IF (BEGIN.LE.NROW) GO TO 1050
      RETURN
 1000 FORMAT (8X,6(4X,A6,2X))
 2000 FORMAT (1X,A6,':',6F12.6)
      END
C  /* Deck radius */
      FUNCTION RADIUS(NCHARGE)
#include <implicit.h>
#include <priunit.h>
#include <qm3.h>
C     Based on covalent radii and metallic radii in Angstrom.
C     Returns -1 where data is inavailable
C     Oct 2006 hjaaj: changed Hydrogen from 30 to 40 pm,
C              such that H2 is printed as bonded ;-) .
      DIMENSION RAD(100)
      DATA (RAD(I), I = 1, 100)/
     &        40.,  155.,  160.,  110.,
     & 90.,   80.,   70.,   68.,   65.,
     &154.,  190.,  160.,  140.,  110.,
     &110.,  105.,  105.,  190.,  238.,
     &200.,  165.,  145.,  135.,  130.,
     &125.,  125.,  125.,  125.,  125.,
     &140.,  140.,  130.,  120.,  120.,
     &120.,  200.,  255.,  215.,  180.,
     &160.,  145.,  140.,  135.,  130.,
     &130.,  135.,  140.,  155.,  160.,
     &160.,  140.,  140.,  140.,  220.,
     &270.,  220.,  185.,  180.,  180.,
     &180.,  180.,  180.,  200.,  180.,
     &175.,  175.,  175.,  175.,  170.,
     &170.,  170.,  155.,  145.,  140.,
     &135.,  135.,  135.,  135.,  145.,
     &155.,  170.,  175.,  170.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100./
C
      IF (NCHARGE .GE. -2 .AND. NCHARGE .LE. 0) THEN
Chj      =  0: solvent cavity center or floating orbital
Chj      = -1: multiple basis for r12 methods
Chj      = -2: point charges
C
         RADIUS = -1.0D0
      ELSE IF (NCHARGE .LT. 1 .OR. NCHARGE .GT. 100) THEN
        IF (.NOT. QM3) THEN
          WRITE (LUPRI,*) 
     &          'ERROR, RADIUS called with NCHARGE =',NCHARGE
          CALL QUIT('RADIUS called with illegal NCHARGE')
        END IF
      ELSE
        RADIUS = 0.01D0 * RAD(NCHARGE)
      END IF
      RETURN
      END
C  /* Deck vdwrad */
      FUNCTION VDWRAD(NCHARGE)
#include <implicit.h>
#include <priunit.h>
C     Based on van der Waals radii in Angstrom.
C     Returns -1 where data is inavailable
      DIMENSION RAD(100)
      DATA (RAD(I), I = 1, 100)/
     &       110.,  220.,  122.,   63.,
     &155.,  155.,  140.,  135.,  130.,
     &154.,  190.,  160.,  140.,  110.,
     &202.,  220.,  150.,  150.,  220.,
     &188.,  181.,  175.,  277.,  239.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100./
C
      IF (NCHARGE .GE. -2 .AND. NCHARGE .LE. 0) THEN
Chj      =  0: solvent cavity center or floating orbital
Chj      = -1: multiple basis for r12 methods
Chj      = -2: point charges
C
         VDWRAD = -1.0D0
      ELSE IF (NCHARGE .LT. 1 .OR. NCHARGE .GT. 100) THEN
         WRITE (LUPRI,*) 'ERROR, VDWRAD called with NCHARGE =',NCHARGE
         CALL QUIT('VDWRAD called with illegal NCHARGE')
      ELSE
         VDWRAD = 0.01D0 * RAD(NCHARGE)
         IF (VDWRAD .LT. 0.0D0) THEN
C           if no table value, use covalent radius plus 0.65 AA
C           (as for B-F above) /Dec. 2006, hjaaj
            WDWRAD = RADIUS(NCHARGE) 
            IF (VDWRAD .GT. 0.0D0) VDWRAD = VDWRAD + 0.65D0
         END IF
      END IF
      RETURN
      END
