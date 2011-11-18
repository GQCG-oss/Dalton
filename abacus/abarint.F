#define OLD_DIAINT
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
#include "single.h"
C
C File:  abacus/abarint.F
C
C New module taking care of redundant internal coordinates.
C The procedures are called from abaopt.F /960818-vebjornb
C
C  /* Deck inired */
      SUBROUTINE INIRED(MXRCRD,MX2CRD,WILBMT,BMTRAN,BMTINV,PJINMT,
     &                  TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,
     &                  WRK,LWRK)
C
C     Initializes everything that has to do with redundant internal
C     coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
      DIMENSION WILBMT(MXRCRD*MX2CRD), BMTRAN(MXRCRD*MXRCRD)
      DIMENSION BMTINV(MXRCRD*MX2CRD), PJINMT(MXRCRD*MXRCRD)
      DIMENSION TMPMAT(MX2CRD,MX2CRD), TMPMT2(MX2CRD,MX2CRD)
      DIMENSION TMPMT3(MX2CRD,MX2CRD), TMPMT4(MX2CRD,MX2CRD)
      DIMENSION TMPMT5(MX2CRD,MX2CRD), TMPMT6(MX2CRD,MX2CRD)
      DIMENSION WRK(LWRK)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     First we have to determine the redundant internal coordinates.
C
      CALL FNDRED(TMPMAT,TMPMT2)
      IF (REDRED) THEN
         CALL RREDUN
         IF (IPRINT .GE. IPRMIN) WRITE(LUPRI,'(/A,I6)')
     &        ' Number of internal coordinates after removals :',IINTCR
      END IF
C
C     Then we construct Wilsons B matrix
C
      CALL GETWIL(MXRCRD,MX2CRD,TMPMAT,WILBMT,BMTRAN,TMPMT2)
C
C     ... and the derivative of the B matrix
C
      IF (IPRINT .GE. IPRDBG)
     &   CALL GETDWL(MXRCRD,TMPMAT,TMPMT2,TMPMT3,WILBMT)
C
C     We also need the inverse of the rectangular B matrix.
C
      CALL GTBINV(MXRCRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT,
     &     BMTRAN,BMTINV,PJINMT,TMPMT5,TMPMT6,WRK,LWRK)
      RETURN
      END

C  /* Deck fndred */
      SUBROUTINE FNDRED(ATMARR,IBOND)
C
C     Finds natural redundant internal coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
      DIMENSION ATMARR(MXCENT,8), IBOND(MXCENT,MXCENT)
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      CALL QENTER('FNDRED')
      IPRSAV = IPRINT
      IF (FINDRE) IPRINT = MAX(IPRINT,IPRMIN)
      CALL DZERO(ATMARR,8*MXCENT)
      CALL IZERO(INTCRD,48*MXCENT)
      CALL IZERO(IBOND,MXCENT*MXCENT)
      IINTCR = 0
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from FNDRED','*',103)
      END IF
C
C     First we initialize the ATMARR array.
C
      CALL ATMINI(ATMARR,IATOM,.FALSE.)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Expanded array of atoms (Aangstrom)',-1)
         CALL OUTPUT(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      ICRTCR = 3*IATOM
C
C     We find all bonds.
C
      CALL FNDBND(ATMARR,IBOND,IATOM)
C
C     Then all angles
C
      CALL FNDANG(ATMARR,IBOND,IATOM)
C
C     Finally all dihedral angles
C
      IF (.NOT. NODIHE) CALL FNDIHD(ATMARR,IBOND,IATOM)
C
C     Determine value of the various coordinates
C
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
      ITMP = IREDIC
      IREDIC = IINTCR
      CALL GETINT(IATOM,8*MXCENT,ATMARR,CRDIN1)
      IREDIC = ITMP
C
C     Output
C
      IF (IPRINT .GE. IPRMIN) CALL PRNRED(IINTCR,CRDIN1)
C
C     Check that the number of internal coordinates do
C     not exceed the limit of 8*MXCENT.
C     If it does, try giving intelligent suggestions.
C
      IF ((IINTCR .GT. 8*MXCENT) .OR.
     &     ((RATFUN .OR. GDIIS) .AND. (IINTCR+1 .GT. 8*MXCENT))) THEN
         WRITE(LUPRI,'(//A//)') '***** CRITICAL ERROR IN OPTIMIZE ' //
     &        'MODULE *****'
         WRITE(LUPRI,'(A,I6,A)') 'The specified geometry requires',
     &        IINTCR,' internal coordinates.'
         WRITE(LUPRI,'(A,I6,A)') 'This exceeds the current limit of',
     &        8*MXCENT,' coordinates.'
         IF (IINTCR .GT. 15*ICRTCR) THEN
            WRITE(LUPRI,'(/A)') 'WARNING: The molecular geometry has' //
     &           ' probably been specified in units of Aangstrom'
            WRITE(LUPRI,'(A)') 'without setting the Aangstrom mark ' //
     &           'in the input file (see manual).'
         ELSE IF (CARTCO .AND. INMDHS) THEN
            WRITE(LUPRI,'(/A)') 'Run the calculation with an ' //
     &           'initial diagonal Hessian, for example:'
            WRITE(LUPRI,'(A)') '          .INITEV'
            WRITE(LUPRI,'(A/)') '          0.6D0'
         ELSE
            WRITE(LUPRI,'(/A)') 'Try running the calculation in ' //
     &           'Cartesian coordinates.'
         END IF
         CALL QUIT('Too many internal coordinates!')
      END IF
      IC = 1
 20   CONTINUE
      IF (IC .LE. IINTCR) THEN
         IF (ICNSTR(IC) .EQ. 2) THEN
            DO 30 I = IC+1, IINTCR
               DO 32 J = 1, 6
                  INTCRD(I-1,J) = INTCRD(I,J)
 32            CONTINUE
               ICNSTR(I-1) = ICNSTR(I)
 30         CONTINUE
            IINTCR = IINTCR - 1
         ELSE
            IC = IC + 1
         END IF
         GOTO 20
      END IF
      IPRINT = IPRSAV
      CALL QEXIT('FNDRED')
      RETURN
      END

C  /* Deck prnred */
      SUBROUTINE PRNRED(ICRD,VALINT)
C
C     Prints out information about redundant internal coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "codata.h"
      DIMENSION VALINT(ICRD)
      CHARACTER CRDTXT*26, MARK*3, UNTTXT*5
      LOGICAL REMCRD
C
      CALL HEADER('Redundant internal coordinates',-1)
      IF (IPRINT .LT. 3) WRITE(LUPRI,'(11X,A/)')
     &   '(Auxiliary bonds are not printed at this print level.)'
      REMCRD = .FALSE.
      IREMCRD = 0
      DO 10 I = 1, ICRD
         MARK = '   '
         IF (CONOPT .AND. (ICNSTR(I) .EQ. 1)) MARK = '*  '
         IF (ICNSTR(I) .EQ. 2) THEN
            MARK = '#  '
            REMCRD = .TRUE.
            IREMCRD = IREMCRD + 1
         END IF
         CRDTXT = '                          '
         IF (INTCRD(I,1) .LT. 10) THEN
            UNTTXT = ' Ang.'
            IF (INTCRD(I,1) .EQ. 1) CRDTXT='Regular bond'
            IF (INTCRD(I,1) .EQ. 3) CRDTXT='Interfragment bond'
            IF (INTCRD(I,1) .EQ. 4) CRDTXT='Additional interfrag. bond'
            IF (INTCRD(I,1) .EQ. 5) CRDTXT='Hydrogen bond'
            IF (INTCRD(I,1) .EQ. 7) CRDTXT='Auxiliary bond'
            IF (INTCRD(I,1) .EQ. 9) CRDTXT='Manual bond'
            IF ((INTCRD(I,1) .NE. 7) .OR. (IPRINT .GE. 3))
     &           WRITE(LUPRI,'(I8,2A,2I4,F22.5,A5)') I,
     &           MARK,CRDTXT,INTCRD(I,2),INTCRD(I,3),
     &           XTANG*VALINT(I),UNTTXT
         ELSE IF (INTCRD(I,1) .GT. 20) THEN
            UNTTXT = ' deg.'
            IF (INTCRD(I,1) .EQ. 21) CRDTXT='Dihedral angle'
            IF (INTCRD(I,1) .EQ. 23) CRDTXT='Extra dihedral angle'
            IF (INTCRD(I,1) .EQ. 25)
     &           CRDTXT='Interfrag. dihedral angle'
            IF (INTCRD(I,1) .EQ. 29) CRDTXT='Manual dihedral angle'
            WRITE(LUPRI,'(I8,2A,4I4,F14.3,A5)') I,
     &           MARK,CRDTXT,INTCRD(I,2),INTCRD(I,3),
     &           INTCRD(I,4),INTCRD(I,5),1.8D2*VALINT(I)/PI,UNTTXT
         ELSE
            UNTTXT = ' deg.'
            IF (INTCRD(I,1) .EQ. 11) CRDTXT='Regular angle'
            IF (INTCRD(I,1) .EQ. 12) CRDTXT='Additional angle'
            IF (INTCRD(I,1) .EQ. 19) CRDTXT='Manual angle'
            WRITE(LUPRI,'(I8,2A,3I4,F18.3,A5)') I,
     &           MARK,CRDTXT,INTCRD(I,2),INTCRD(I,3),
     &           INTCRD(I,4),180.D0*VALINT(I)/PI,UNTTXT
         END IF
 10   CONTINUE
      IF (CONOPT .OR. REMCRD) THEN
         WRITE(LUPRI,*)
         IF (CONOPT) WRITE(LUPRI,*)
     &        '*) Constrained (fixed) coordinate.'
         IF (REMCRD) WRITE(LUPRI,*)
     &        '#) Coordinate will be removed.'
      END IF
      WRITE(LUPRI,'(/A,I6)')
     &     ' Total number of redundant internal coordinates:',ICRD
      IF (IREMCRD .GT. 0) WRITE(LUPRI,'(/A,I6)')
     &     ' Total number of removed internal coordinates:  ',IREMCRD
      RETURN
      END

C  /* Deck redvib */
      SUBROUTINE REDVIB(NCORD,NMODE,MXRCRD,MX2CRD,EVEC,IFRQCM,IMAGIN,
     &     ATMARR,CRDRIN,WILBMT,BMTRAN,TMPMAT,TMPMT2)
C
C     Converts normal coordinates to redundant internal coordinates.
C
#include "implicit.h"
C
#include "codata.h"
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "priunit.h"
#include "optinf.h"
#include "symmet.h"
#include "nuclei.h"
      DIMENSION EVEC(NCORD,NCORD), IFRQCM(NCORD)
      DIMENSION IMAGIN(NCORD), ATMARR(MXCENT,8)
      DIMENSION CRDRIN(MXRCRD), WILBMT(MXRCRD,MXCOOR)
      DIMENSION BMTRAN(MXRCRD,MXRCRD), TMPMAT(MX2CRD,MX2CRD)
      DIMENSION TMPMT2(MX2CRD,MX2CRD)
      LOGICAL LOGDEL
      CHARACTER*4 CRDTXT
      CHARACTER*1 CHRIMG(0:1)
      DATA CHRIMG /' ','i'/
#include "ibtfun.h"
C
      SXFAMU = SQRT(XFAMU)
      IF (.NOT. (REDINT .OR. DELINT)) RETURN
      ITMP = IREDIC
      LOGDEL = DELINT
      IC = IINTCR
      IF (DELINT) IC = IREDIC
      DELINT = .FALSE.
      CALL DZERO(TMPMT2,MX2CRD*MX2CRD)
      CALL HEADER('Analysis in redundant internal coordinates',1)
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
C
C     We have to subtract the components of the recently
C     determined step.
C
      CALL WLKCOR(STPSYM,TMPMT2,ICRTCR,MXCOOR,-1)
      DO 10 ICENT = 1, NUCIND
         DO 15 J = 1, 3
            CORD(J,ICENT) = CORD(J,ICENT) - TMPMT2(3*(ICENT-1)+J,1)
 15      CONTINUE
 10   CONTINUE
C
C     We have to calculate the Wilson B matrix and the
C     values of all redundant internal coordinates.
C
      CALL GETINT(IATOM,MXRCRD,ATMARR,CRDRIN)
      IP = IPRINT
      IPRINT = -1
      CALL GETWIL(MXRCRD,MX2CRD,ATMARR,WILBMT,BMTRAN,TMPMAT)
      IPRINT = IP
      CALL HEADER('Expanded array of atoms (Bohr)',-1)
      CALL OUTPUT(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      CALL PRNRED(IC,CRDRIN)
C
C     The molecular coordinates are restored.
C
      DO 20 ICENT = 1, NUCIND
         DO 25 J = 1, 3
            CORD(J,ICENT) = CORD(J,ICENT) + TMPMT2(3*(ICENT-1)+J,1)
 25      CONTINUE
 20   CONTINUE
C
C     Finally the eigenvectors are transformed
C
      CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT2,MX2CRD*MX2CRD)
      DO 40 IMODE = 1, NMODE
         DO 45 II = 1, IC
            DO 47 ICOOR = 1, 3*IATOM
               TMPMAT(II,IMODE) = TMPMAT(II,IMODE) +
     &              SXFAMU*EVEC(ICOOR,IMODE)*WILBMT(II,ICOOR)
 47         CONTINUE
 45      CONTINUE
 40   CONTINUE
      CALL HEADER('Normal coordinates (internal)',-1)

      ISTR = 1
      NBATCH = (NMODE + 4)/5
      DO 100 IBATCH = 1, NBATCH
         IEND = MIN(ISTR + 4,NMODE)
         NUMB = IEND - ISTR + 1
         WRITE (LUPRI,'(/A12,5(I5,A2,I4,A1))') '            ',
     &      (II,'  ',IFRQCM(II),CHRIMG(IMAGIN(II)), II = ISTR,IEND)
         LENH = 10 + NUMB*12
         WRITE (LUPRI,'(2X,70A1)') ('-', II = 1,LENH)
         WRITE (LUPRI,'()')
         DO 110 ICOOR = 1, IINTCR
            CRDTXT = 'bend'
            IF (INTCRD(ICOOR,1) .LT. 10) CRDTXT = 'stre'
            IF (INTCRD(ICOOR,1) .GE. 20) CRDTXT = 'tors'
            IF (INTCRD(ICOOR,1) .NE. 7)
     &           WRITE (LUPRI,'(I5,1X,A4,(T13,5F12.6))') ICOOR,CRDTXT,
     &           (TMPMAT(ICOOR,II),II=ISTR,IEND)
 110     CONTINUE
         ISTR = ISTR + 5
 100  CONTINUE
      IREDIC = ITMP
      DELINT = LOGDEL
      RETURN
      END

C  /* Deck fndbnd */
      SUBROUTINE FNDBND(ATMARR,IBOND,IATOM)
C
C     Finds all bonds.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
#include "pi.h"
      DIMENSION ATMARR(MXCENT,8), VEC1(3), VEC2(3), IBOND(MXCENT,MXCENT)
      LOGICAL POSSIB
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We have to find all regular bonds.
C
      IBEFOR = IINTCR
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Regular bonds',-1)
      END IF
      DO 10 I = 1, IATOM - 1
         DO 15 J = I + 1, IATOM
            RADI = ATMARR(I,5)
            RADJ = ATMARR(J,5)            
            VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
            VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
            VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
            DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
C
C     Regular bonds have dist. less than 1.3 times sum of cov. radii
C
            IF (DIST .LE. 1.3D0*(RADI+RADJ)) THEN
               IBOND(I,J) = 1
               IBOND(J,I) = 1
               IINTCR = IINTCR + 1
               INTCRD(IINTCR,1) = 1
               INTCRD(IINTCR,2) = I
               INTCRD(IINTCR,3) = J
               IF (IPRINT .GE. IPRMAX) THEN
                  WRITE(LUPRI,'(A,2I4)')
     &                 ' Bond between atoms      : ',i,j
               END IF
            END IF
 15      CONTINUE
 10   CONTINUE
      IF ((IPRINT .GE. IPRMAX) .AND. (IBEFOR .EQ. IINTCR)) THEN
         WRITE(LUPRI,'(A)') ' None were found.'
      END IF
C
C     Next we have to check if single atoms or fragments are too far
C     away to be bonded by normal bonds, in which case interfragment
C     bond(s) must be assigned.
C
C     The first step is to assign fragment numbers to all atoms, this
C     process is done iteratively.
C
      IBEFOR = IINTCR
      IFRAG = 1
      IAT = 1
 20   CONTINUE
      ATMARR(IAT,6) = IFRAG*1.0D0
 25   CONTINUE
      ICHANG = 0
      DO 30 II = 1,IATOM
         IF (NINT(ATMARR(II,6)) .EQ. IFRAG) THEN
            DO 35 JJ = 1, IATOM
               IF ((IBOND(II,JJ) .GE. 1) .AND.
     &              (NINT(ATMARR(JJ,6)) .EQ. 0)) THEN
                  ATMARR(JJ,6) = IFRAG*1.0D0
                  ICHANG = 1
               END IF
 35         CONTINUE
         END IF
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Assigning atoms to fragments',-1)
         CALL OUTPUT(ATMARR,1,IATOM,6,6,MXCENT,8,1,LUPRI)
      END IF
      IF (ICHANG .GT. 0) GOTO 25
C
      II = 1
 32   CONTINUE
      IF ((II .LE. IATOM) .AND. (ATMARR(II,6) .GT. 0.5D0)) THEN
         II = II + 1
         GOTO 32
      END IF
      IF (II .LE. IATOM) THEN
         IFRAG = IFRAG + 1
         IAT = II
         ATMARR(IAT,6) = IFRAG*1.0D0
         GOTO 20
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,'(/A,I4)') 
     &        ' Number of separate fragments found: ', IFRAG
      END IF
C
C     Two fragments get connected through the shortest interfragment distance.
C     All distances shorter than the minimum of 2.0 Ang. and this distance
C     are also marked as interfragment bonds.
C
      IF (IFRAG .GT. 1) THEN
         DSTMIN = 1.0D10
         DO 40 II = 1, IATOM
            DO 45 JJ = 1, IATOM
               IF (NINT(ATMARR(II,6)) .NE. NINT(ATMARR(JJ,6))) THEN
                  VEC1(1) = ATMARR(II,2)-ATMARR(JJ,2)
                  VEC1(2) = ATMARR(II,3)-ATMARR(JJ,3)
                  VEC1(3) = ATMARR(II,4)-ATMARR(JJ,4)
                  DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                  IF (DIST .LT. DSTMIN) THEN
                     DSTMIN = DIST
                     IATM1 = II
                     IFRG1 = NINT(ATMARR(II,6))
                     IATM2 = JJ
                     IFRG2 = NINT(ATMARR(JJ,6))
                  END IF
               END IF
 45         CONTINUE
 40      CONTINUE
         IBOND(IATM1,IATM2) = 3
         IBOND(IATM2,IATM1) = 3
         IINTCR = IINTCR + 1
         INTCRD(IINTCR,1) = 3
         INTCRD(IINTCR,2) = IATM1
         INTCRD(IINTCR,3) = IATM2
         BNDLIM = MIN(1.3D0*DSTMIN, 2.0D0)
C
C     Additional interfragment bonds
C
         DO 50 II = 1, IATOM
            IF (NINT(ATMARR(II,6)) .EQ. IFRG1) THEN
               DO 55 JJ = 1, IATOM
                  IF ((IBOND(II,JJ) .EQ. 0) .AND.
     &                 (NINT(ATMARR(JJ,6)) .EQ. IFRG2)) THEN
                     VEC1(1) = ATMARR(II,2)-ATMARR(JJ,2)
                     VEC1(2) = ATMARR(II,3)-ATMARR(JJ,3)
                     VEC1(3) = ATMARR(II,4)-ATMARR(JJ,4)
                     DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     IF (DIST .LE. BNDLIM) THEN
                        BNDTYP = 4
C
C     We make a distinction between interfragment bonds that are
C     approximately as short as the shortest one and additional
C     longer interfragment bonds.
C
                        IF (DIST .LE. DSTMIN*1.05D0) BNDTYP = 3
                        IBOND(II,JJ) = BNDTYP
                        IBOND(JJ,II) = BNDTYP
                        IINTCR = IINTCR + 1
                        INTCRD(IINTCR,1) = BNDTYP
                        IF (II .LT. JJ) THEN
                           INTCRD(IINTCR,2) = II
                           INTCRD(IINTCR,3) = JJ
                        ELSE
                           INTCRD(IINTCR,2) = JJ
                           INTCRD(IINTCR,3) = II
                        END IF
                     END IF
                  END IF
 55            CONTINUE
            END IF
 50      CONTINUE
C
C     The process of connecting fragments runs iteratively, until all
C     atoms are connected, i.e. the number of fragments is one.
C
         IF (IFRAG .GT. 1) THEN
            DO 57 I = 1, IATOM
               ATMARR(I,6) = D0
 57         CONTINUE
            IAT = 1
            IFRAG = 1
            GOTO 20
         END IF
      END IF
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Interfragment bonds',-1)
         DO 60 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,2I4)')
     &           ' Interfragment bond between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3)
 60      CONTINUE
      END IF
C
C     Then we have to find all hydrogen bonds in the system.
C     X --- H ... Y (X,Y = N, O, F, P, S, Cl)
C
      IBEFOR = IINTCR
      DO 70 I = 1, IATOM
         IF (NINT(ATMARR(I,1)) .EQ. 1) THEN
            POSSIB = .FALSE.
            DO 72 J = 1, IATOM
               ID = NINT(ATMARR(J,1))
               IF ((IBOND(I,J) .EQ. 1) .AND. ((ID .EQ. 7) .OR.
     &              (ID .EQ. 8) .OR. (ID .EQ. 9) .OR. (ID .EQ. 15)
     &              .OR. (ID .EQ. 16) .OR. (ID .EQ. 17))) THEN
                  POSSIB = .TRUE.
                  VEC2(1) = ATMARR(J,2)-ATMARR(I,2)
                  VEC2(2) = ATMARR(J,3)-ATMARR(I,3)
                  VEC2(3) = ATMARR(J,4)-ATMARR(I,4)
               END IF
 72         CONTINUE
            IF (POSSIB) THEN
               VDWR1 = VDWRAD(1)
               DO 75 J = 1, IATOM
                  ID = NINT(ATMARR(J,1))
                  IF ((IBOND(I,J) .EQ. 0) .AND. ((ID .EQ. 7) .OR.
     &                 (ID .EQ. 8) .OR. (ID .EQ. 9) .OR. (ID .EQ. 15)
     &                  .OR. (ID .EQ. 16) .OR. (ID .EQ. 17))) THEN
                     VEC1(1) = ATMARR(J,2)-ATMARR(I,2)
                     VEC1(2) = ATMARR(J,3)-ATMARR(I,3)
                     VEC1(3) = ATMARR(J,4)-ATMARR(I,4)
                     DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     ICHG = NINT(ATMARR(J,1))
C
C     We require the distance to be less than 0.9 times the sum of the
C     van der Waals radii and larger than the covalent radii. Finally
C     we say that the angle (X-H-Y) has to be larger than 90 degrees.
C
                     VDWRI = VDWRAD(ICHG)
                     IF ((DIST .LE. 0.9D0*(VDWR1+VDWRI))
     &                .AND. (DIST .GE. (RADIUS(1)+RADIUS(ICHG)))) THEN
                        BNDANG = VECANG(VEC1,VEC2)
                        IF (ABS(BNDANG) .GE. 0.49D0*PI) THEN
                           IBOND(I,J) = 5
                           IBOND(J,I) = 5
                           IINTCR = IINTCR + 1
                           INTCRD(IINTCR,1) = 5
                           IF (I .LT. J) THEN
                              INTCRD(IINTCR,2) = I
                              INTCRD(IINTCR,3) = J
                           ELSE
                              INTCRD(IINTCR,2) = J
                              INTCRD(IINTCR,3) = I
                           END IF
                        END IF
                     END IF
                  END IF
 75            CONTINUE
            END IF
         END IF
 70   CONTINUE
C
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Hydrogen bonds',-1)
         DO 77 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,2I4)')
     &           ' Hydrogen bond between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3)
 77      CONTINUE
      END IF
C
C     Bond coordinates that have been manually specified are added
C     (if not already present)
C
      IBEFOR = IINTCR
      IF (ADDCRD .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Manually added bonds',-1)
      END IF
      DO 700 I = 1, IADDCR(0,1)
         IF ((IADDCR(I,3) .EQ. 0) .AND. (IADDCR(I,4) .EQ. 0)) THEN
            IF ((IADDCR(I,1) .GT. IATOM)
     &           .OR. (IADDCR(I,2) .GT. IATOM)) THEN
               WRITE(LUPRI,'(A,I4)')
     &              ' Critical error in additional coordinate #', I
               WRITE(LUPRI,'(A)') ' Chosen atom number '
     &              // 'is larger than the number of atoms.'
               CALL QUIT(' Error detected in the specification of '
     &              // 'additional bond coordinate(s)!')
            END IF
            IF (IBOND(IADDCR(I,1),IADDCR(I,2)) .EQ. 0) THEN
               IINTCR = IINTCR + 1
               INTCRD(IINTCR,1) = 9
               INTCRD(IINTCR,2) = IADDCR(I,1)
               INTCRD(IINTCR,3) = IADDCR(I,2)
            END IF
         END IF
 700  CONTINUE
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         DO 710 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,2I4)')
     &           ' Manual bond between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3)
 710     CONTINUE
      END IF
C
C     We also determine "extra bonds" between atoms. These will
C     not generate angles and dihedral angles, but they seem to
C     make geometry optimization in redundant coordinates more
C     efficient. Only used for approx. Hessians.
C
      IBEFOR = IINTCR
C      IF (.NOT. NOAUX) THEN
      IF (.NOT. NOAUX .AND. (MODHES .OR. CMBMOD .OR. INMDHS)) THEN
         IF (IPRINT .GE. IPRMAX) CALL HEADER('Auxiliary bonds',-1)
         DO 80 I = 1, IATOM - 1
            DO 85 J = I + 1, IATOM
               RADI = ATMARR(I,5)
               RADJ = ATMARR(J,5)            
               VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
               VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
               VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
               DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
C
C     Extra bonds have dist. less than 2.5 times sum of cov. radii
C
               IF ((DIST .LE. 2.5D0*(RADI+RADJ)) .AND.
     &              (IBOND(I,J) .EQ. 0)) THEN
                  IINTCR = IINTCR + 1
                  INTCRD(IINTCR,1) = 7
                  INTCRD(IINTCR,2) = I
                  INTCRD(IINTCR,3) = J
                  IF (IPRINT .GE. IPRMAX) WRITE(LUPRI,'(A,2I4)')
     &                 ' Auxiliary bond between atoms: ',i,j
               END IF
 85         CONTINUE
 80      CONTINUE
         IF ((IPRINT .GE. IPRMAX) .AND. (IBEFOR .EQ. IINTCR)) 
     &        WRITE(LUPRI,'(A)') ' None were added.'
      END IF
      RETURN
      END

C  /* Deck fndang */
      SUBROUTINE FNDANG(ATMARR,IBOND,IATOM)
C
C     Finds all angles.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
#include "pi.h"
      DIMENSION ATMARR(MXCENT,8), VEC1(3), VEC2(3), IBOND(MXCENT,MXCENT)
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We have to find all angles between bonded atoms
C
      IBEFOR = IINTCR
      DO 10 I = 1, IATOM
         NUMBND = 0
         DO 15 J = 1, IATOM
            IF ((IBOND(I,J) .GT. 0) .AND. (IBOND(I,J) .NE. 9))
     &           NUMBND = NUMBND + 1
 15      CONTINUE
         IF (NUMBND .GT. 1) THEN
            DO 20 J = 1, IATOM - 1
               IF ((IBOND(I,J) .GT. 0) .AND. (IBOND(I,J) .NE. 4) 
     &              .AND. (IBOND(I,J) .NE. 9)) THEN
                  VEC1(1) = ATMARR(J,2)-ATMARR(I,2)
                  VEC1(2) = ATMARR(J,3)-ATMARR(I,3)
                  VEC1(3) = ATMARR(J,4)-ATMARR(I,4)
                  DO 25 K = J + 1, IATOM
                     IF ((IBOND(I,K) .GT. 0) .AND.
     &                    (IBOND(I,K) .NE. 4) .AND.
     &                    (IBOND(I,K) .NE. 9)) THEN
                        IINTCR = IINTCR + 1
                        INTCRD(IINTCR,1) = 11
                        INTCRD(IINTCR,2) = J
                        INTCRD(IINTCR,3) = I
                        INTCRD(IINTCR,4) = K
C
C     We have to check for linear or near linear bonds (these get two
C     internal orthogonal bond coordinates). The "critical" bond angle
C     is 175 degrees.
C
                        VEC2(1) = ATMARR(K,2)-ATMARR(I,2)
                        VEC2(2) = ATMARR(K,3)-ATMARR(I,3)
                        VEC2(3) = ATMARR(K,4)-ATMARR(I,4)
                        BNDANG = VECANG(VEC1,VEC2)
                        IF ((ABS(BNDANG) .GE. 175D0*PI/180D0) .AND.
     &                       (.NOT. NOADDA)) THEN
                           IINTCR = IINTCR + 1
                           INTCRD(IINTCR,1) = 12
                           INTCRD(IINTCR,2) = J
                           INTCRD(IINTCR,3) = I
                           INTCRD(IINTCR,4) = K
                        END IF
                     END IF
 25               CONTINUE
               END IF
 20         CONTINUE
         END IF
 10   CONTINUE
C
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Angles',-1)
         DO 30 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,3I4)')
     &           ' Angles between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3),INTCRD(I,4)
 30      CONTINUE
      END IF
C
C     Angle coordinates that have been manually specified are added
C     (if not already present)
C
      IBEFOR = IINTCR
      IF (ADDCRD .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Manually added angles',-1)
      END IF
      DO 50 I = 1, IADDCR(0,1)
         IF ((IADDCR(I,3) .GT. 0) .AND. (IADDCR(I,4) .EQ. 0)) THEN
            IF ((IADDCR(I,1) .GT. IATOM) .OR. (IADDCR(I,2) .GT. IATOM) 
     &           .OR. (IADDCR(I,3) .GT. IATOM))  THEN
               WRITE(LUPRI,'(A,I4)')
     &              ' Critical error in additional coordinate #', I
               WRITE(LUPRI,'(A)') ' Chosen atom number '
     &              // 'is larger than the number of atoms.'
               CALL QUIT(' Error detected in the specification of '
     &              // 'additional angle coordinate(s)!')
            END IF
            IBND1 = IBOND(IADDCR(I,1),IADDCR(I,2))
            IBND2 = IBOND(IADDCR(I,2),IADDCR(I,3))
C     Check if the three atoms should already have produced an angle coordinate
            IF (.NOT.(((IBND1.GT.0).AND.(IBND1.NE.4).AND.(IBND1.NE.9))
     &      .AND.((IBND2.GT.0).AND.(IBND2.NE.4).AND.(IBND2.NE.9)))) THEN
               IINTCR = IINTCR + 1
               INTCRD(IINTCR,1) = 19
               INTCRD(IINTCR,2) = IADDCR(I,1)
               INTCRD(IINTCR,3) = IADDCR(I,2)
               INTCRD(IINTCR,4) = IADDCR(I,3)
            END IF
         END IF
 50   CONTINUE
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         DO 710 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,3I4)')
     &           ' Manual angle between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3),INTCRD(I,4)
 710     CONTINUE
      END IF
      RETURN
      END

C  /* Deck fndihd */
      SUBROUTINE FNDIHD(ATMARR,IBOND,IATOM)
C
C     Finds all dihedral angles.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
#include "pi.h"
      DIMENSION ATMARR(MXCENT,8), IBOND(MXCENT,MXCENT)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3), IBNDCN(16)
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (DEG175 = 175D0*PI/180D0)
      LOGICAL FRGBND, FOUND
C
C     We have to find all dihedral angles between bonded atoms
C
      IF (IATOM .LE. 3) RETURN
      FRGBND = .FALSE.
      INODIH = IINTCR
      IBEFOR = IINTCR
      DO 10 I = 1, IATOM
         DO 15 J = 1, IATOM
            IF (IBOND(I,J) .EQ. 3) FRGBND = .TRUE.
            IF ((IBOND(I,J) .GT. 0) .AND. (IBOND(I,J) .NE. 4)) THEN
               DO 20 K = 1, IATOM
                  IF ((IBOND(J,K) .GT. 0) .AND. (K .NE. I) .AND.
     &                 (IBOND(J,K) .NE. 4)) THEN 
                     DO 25 L = I, IATOM
                        IF ((IBOND(K,L) .GT. 0) .AND. (L .NE. J) .AND.
     &                       (L .NE. I) .AND. (IBOND(K,L) .NE. 4)) THEN
C
C     Check for linearity.
C
                           VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
                           VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
                           VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
                           VEC2(1) = ATMARR(K,2)-ATMARR(J,2)
                           VEC2(2) = ATMARR(K,3)-ATMARR(J,3)
                           VEC2(3) = ATMARR(K,4)-ATMARR(J,4)
                           BNDAN1  = VECANG(VEC1,VEC2)
                           VEC1(1) = -VEC2(1)
                           VEC1(2) = -VEC2(2)
                           VEC1(3) = -VEC2(3)
                           VEC2(1) = ATMARR(L,2)-ATMARR(K,2)
                           VEC2(2) = ATMARR(L,3)-ATMARR(K,3)
                           VEC2(3) = ATMARR(L,4)-ATMARR(K,4)
                           BNDAN2  = VECANG(VEC1,VEC2)
                           IF ((ABS(BNDAN1) .LT. DEG175)
     &                          .AND. (ABS(BNDAN2) .LT. DEG175)) THEN
                              IINTCR = IINTCR + 1
                              INTCRD(IINTCR,1) = 21
                              INTCRD(IINTCR,2) = I
                              INTCRD(IINTCR,3) = J
                              INTCRD(IINTCR,4) = K
                              INTCRD(IINTCR,5) = L
                           END IF
                        END IF
 25                  CONTINUE
                  END IF
 20            CONTINUE
            END IF
 15      CONTINUE
 10   CONTINUE
C
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Dihedral angles',-1)
         DO 30 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,4I4)')
     &           ' Dihedral angles between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3),INTCRD(I,4),INTCRD(I,5)
 30      CONTINUE
      END IF
C
C     For weakly bound systems, we want to make sure that we add
C     dihedrals across the interfragment bonds (if not already present)
C
      IF (FRGBND) THEN
         DO 40 IFBCRD = 1, IINTCR
            IF (INTCRD(IFBCRD,1) .EQ. 3) THEN
               I = INTCRD(IFBCRD,2)
               J = INTCRD(IFBCRD,3)
               FOUND = .FALSE.
C
C     Check if dihedral(s) already span this interfragment bond
C
               DO 41 ICRD = 1, IINTCR
                  IF (INTCRD(ICRD,1) .EQ. 21) THEN
                     IF ((INTCRD(ICRD,3) .EQ. I) .AND.
     &                    (INTCRD(ICRD,4) .EQ. J))
     &                    FOUND = .TRUE.
                     IF ((INTCRD(ICRD,3) .EQ. J) .AND.
     &                    (INTCRD(ICRD,4) .EQ. I))
     &                    FOUND = .TRUE.
                  END IF
 41            CONTINUE
C
C     If not, we have to do a wider search. What we are looking for are
C     dihedrals A-B-C-D, where the actual bonding situation is A--B...C--X--D,
C     i.e. BC is the interfragment bond while C and D are both bonded
C     to the same atom X
C
               IF (.NOT. FOUND) THEN
                  DO 42 ICRD = 1, IINTCR
                     IF (INTCRD(ICRD,1) .EQ. 1) THEN
                        II1 = 0
                        IF (INTCRD(ICRD,2) .EQ. I) THEN
                           II1 = INTCRD(ICRD,3)
                           II2 = I
                           II3 = J
                        ELSE IF (INTCRD(ICRD,3) .EQ. I) THEN
                           II1 = INTCRD(ICRD,2)
                           II2 = I
                           II3 = J
                        ELSE IF (INTCRD(ICRD,2) .EQ. J) THEN
                           II1 = INTCRD(ICRD,3)
                           II2 = J
                           II3 = I
                        ELSE IF (INTCRD(ICRD,3) .EQ. J) THEN
                           II1 = INTCRD(ICRD,2)
                           II2 = J
                           II3 = I
                        END IF
                        IF (II1 .GT. 0) THEN
                           DO 43 IAT1 = 1, IATOM
                              IF (IBOND(IAT1,II3) .EQ. 1) THEN
                                 DO 44 IAT2 = 1, IATOM
                                    IF ((IBOND(IAT1,IAT2) .EQ. 1) .AND.
     &                                   (IAT2 .NE. II3)) THEN
                                       II4 = IAT2
C     Check for linearity.
                                       VEC1(1) = ATMARR(II1,2)
     &                                      - ATMARR(II2,2)
                                       VEC1(2) = ATMARR(II1,3)
     &                                      - ATMARR(II2,3)
                                       VEC1(3) = ATMARR(II1,4)
     &                                      - ATMARR(II2,4)
                                       VEC2(1) = ATMARR(II3,2)
     &                                      - ATMARR(II2,2)
                                       VEC2(2) = ATMARR(II3,3)
     &                                      - ATMARR(II2,3)
                                       VEC2(3) = ATMARR(II3,4)
     &                                      - ATMARR(II2,4)
                                       BNDAN1  = VECANG(VEC1,VEC2)
                                       VEC1(1) = -VEC2(1)
                                       VEC1(2) = -VEC2(2)
                                       VEC1(3) = -VEC2(3)
                                       VEC2(1) = ATMARR(II4,2)
     &                                      - ATMARR(II3,2)
                                       VEC2(2) = ATMARR(II4,3)
     &                                      - ATMARR(II3,3)
                                       VEC2(3) = ATMARR(II4,4)
     &                                      - ATMARR(II3,4)
                                       BNDAN2  = VECANG(VEC1,VEC2)
                                       IF ((ABS(BNDAN1).LT.DEG175).AND.
     &                                    (ABS(BNDAN2).LT.DEG175)) THEN
                                          IINTCR = IINTCR + 1
                                          INTCRD(IINTCR,1) = 25
                                          INTCRD(IINTCR,2) = II1
                                          INTCRD(IINTCR,3) = II2
                                          INTCRD(IINTCR,4) = II3
                                          INTCRD(IINTCR,5) = II4
                                       END IF
                                    END IF
 44                              CONTINUE
                              END IF
 43                        CONTINUE
                        END IF
                     END IF
 42               CONTINUE
               END IF
            END IF
 40      CONTINUE
      END IF
C
C     Systems with more than three atoms where no regular dihedral
C     angles has been found, need addition of dihedral angles to take
C     care of out of plane bending. First we find the central atom.
C
      ICNT = 0
      IBNDNR = 0
      IF ((IATOM .GE. 4) .AND. (IBEFOR .EQ. IINTCR)) THEN
         DO 50 I = 1, IATOM
            INUMBN = 0
            DO 55 J = 1, IATOM
               IF (IBOND(I,J) .GT. 0) INUMBN = INUMBN + 1
 55         CONTINUE
            IF (INUMBN .GE. 3) THEN
               ICNT   = I
               IBNDNR = INUMBN
            END IF
 50      CONTINUE
      END IF
      IF ((ICNT .GT. 0) .AND. (IBNDNR .LE. 16)) THEN
         CALL IZERO(IBNDCN,16)
         IBEFOR = IINTCR
         II = 1
         DO 60 I = 1, IATOM
            IF (IBOND(I,ICNT) .GT. 0) THEN
               IBNDCN(II) = I
               II = II + 1
            END IF
 60      CONTINUE
         DO 65 II = 1, IBNDNR
            VEC1(1) = ATMARR(IBNDCN(II),2)-ATMARR(ICNT,2)
            VEC1(2) = ATMARR(IBNDCN(II),3)-ATMARR(ICNT,3)
            VEC1(3) = ATMARR(IBNDCN(II),4)-ATMARR(ICNT,4)
            DO 67 JJ = 1, IBNDNR
               IF (JJ .NE. II) THEN
                  VEC2(1) = ATMARR(IBNDCN(JJ),2)-ATMARR(ICNT,2)
                  VEC2(2) = ATMARR(IBNDCN(JJ),3)-ATMARR(ICNT,3)
                  VEC2(3) = ATMARR(IBNDCN(JJ),4)-ATMARR(ICNT,4)
                  BNDAN1  = VECANG(VEC1,VEC2)
                  IF (ABS(BNDAN1) .LT. DEG175) THEN
                     VEC3(1) = -VEC2(1)
                     VEC3(2) = -VEC2(2)
                     VEC3(3) = -VEC2(3)
                     DO 69 KK = 1, IBNDNR
                        IF ((KK .NE. II) .AND. (KK .NE. JJ)
     &                       .AND. (IINTCR .EQ.IBEFOR)) THEN
                           VEC4(1) = ATMARR(IBNDCN(KK),2)
     &                          -ATMARR(IBNDCN(JJ),2)
                           VEC4(2) = ATMARR(IBNDCN(KK),3)
     &                          -ATMARR(IBNDCN(JJ),3)
                           VEC4(3) = ATMARR(IBNDCN(KK),4)
     &                          -ATMARR(IBNDCN(JJ),4)
                           BNDAN2  = VECANG(VEC3,VEC4)
                           IF (ABS(BNDAN2) .LT. DEG175) THEN
C
C     We need well-balanced dihedral angles, and once we have found one
C     suitable angle, we rely on permutational symmetry.
C
C     First we have to do some sorting.
C
                              I = II
                              J = JJ
                              K = KK
                              ITMP = IBNDCN(1)
                              IBNDCN(1) = IBNDCN(I)
                              IBNDCN(I) = ITMP
                              IF (J .EQ. 1) J = I
                              IF (K .EQ. 1) K = I
                              I = 1
                              ITMP = IBNDCN(2)
                              IBNDCN(2) = IBNDCN(J)
                              IBNDCN(J) = ITMP
                              IF (K .EQ. 2) K = J
                              J = 2
                              ITMP = IBNDCN(3)
                              IBNDCN(3) = IBNDCN(K)
                              IBNDCN(K) = ITMP
                              K = 3
C
C     Then we find all the dihedral angles
C
                              DO 100 LL = 0, IBNDNR - 1
                                 IINTCR = IINTCR + 1
                                 INTCRD(IINTCR,1) = 23
                                 INTCRD(IINTCR,2) =
     &                                IBNDCN(MOD(LL,IBNDNR)+1)
                                 INTCRD(IINTCR,3) = ICNT
                                 INTCRD(IINTCR,4) =
     &                                IBNDCN(MOD(LL+1,IBNDNR)+1)
                                 INTCRD(IINTCR,5) =
     &                                IBNDCN(MOD(LL+2,IBNDNR)+1)
 100                          CONTINUE

 123                          CONTINUE

                           END IF
                        END IF
 69                  CONTINUE
                  END IF
               END IF
 67         CONTINUE
 65      CONTINUE
         IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
            CALL HEADER('Extra dihedral angles',-1)
            DO 70 I = IBEFOR + 1, IINTCR
               WRITE(LUPRI,'(A,4I4)')
     &              ' Dihedral angles between atoms: ',
     &              INTCRD(I,2),INTCRD(I,3),INTCRD(I,4),INTCRD(I,5)
 70         CONTINUE
         END IF
      END IF
C
C     Dihedral angle coordinates that have been manually specified are added
C     (if not already present)
C
      IBEFOR = IINTCR
      IF (ADDCRD .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Manually added dihedral angles',-1)
      END IF
      DO 90 I = 1, IADDCR(0,1)
         IF ((IADDCR(I,3) .GT. 0) .AND. (IADDCR(I,4) .GT. 0)) THEN
            IF ((IADDCR(I,1) .GT. IATOM) .OR. (IADDCR(I,2) .GT. IATOM)
     &           .OR. (IADDCR(I,3) .GT. IATOM)
     &           .OR. (IADDCR(I,4) .GT. IATOM)) THEN
               WRITE(LUPRI,'(A,I4)')
     &              ' Critical error in additional coordinate #', I
               WRITE(LUPRI,'(A)') ' Chosen atom number '
     &              // 'is larger than the number of atoms.'
               CALL QUIT(' Error detected in the specification of '
     &              // 'additional dihedral coordinate(s)!')
            END IF
            IBND1 = IBOND(IADDCR(I,1),IADDCR(I,2))
            IBND2 = IBOND(IADDCR(I,2),IADDCR(I,3))
            IBND3 = IBOND(IADDCR(I,3),IADDCR(I,4))
C     Check if the four atoms should already have produced a dihedral coordinate
            IF (.NOT.(((IBND1.GT.0).AND.(IBND1.NE.4).AND.(IBND1.NE.9))
     &      .AND.((IBND2.GT.0).AND.(IBND2.NE.4).AND.(IBND2.NE.9))
     &      .AND.((IBND3.GT.0).AND.(IBND3.NE.4).AND.(IBND3.NE.9)))) THEN
               IINTCR = IINTCR + 1
               INTCRD(IINTCR,1) = 29
               INTCRD(IINTCR,2) = IADDCR(I,1)
               INTCRD(IINTCR,3) = IADDCR(I,2)
               INTCRD(IINTCR,4) = IADDCR(I,3)
               INTCRD(IINTCR,5) = IADDCR(I,4)
            END IF
         END IF
 90   CONTINUE
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         DO 95 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,4I4)')
     &           ' Manual dihedral between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3),INTCRD(I,4),INTCRD(I,5)
 95      CONTINUE
      END IF
      RETURN
      END

C  /* Deck rredun */
      SUBROUTINE RREDUN
C
C     Reduces the redundancy of the redundant internal coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      IBEFOR = IINTCR
      NDIHD = 0
      DO 10 I = 1, IINTCR
         IF (INTCRD(I,1) .GE. 20) NDIHD = NDIHD + 1
 10   CONTINUE
C
C     The number of dihedral angles should be larger than 10 and
C     comprise more than 1/3 of the total number of coordinates.
C
      IF ((NDIHD .LT. 10) .OR. (NDIHD .LT. IINTCR/3)) RETURN
C
      IFIRST = IINTCR - NDIHD + 1
      IC = IFIRST
C
C     We pick out every fifth dihedral angle.
C
 20   CONTINUE
      IF (MOD(IC*1.0D0, 10.0D0) .GT. 1.0D-3) THEN
         IC = IC + 1
         GOTO 20
      END IF
C
C     Then we mark all the internal coordinates that have the same value.
C
 25   CONTINUE
      VAL = ABS(CRDIN1(IC))
      IF (VAL .LT. 1.0D6) THEN
         DO 30 I = IFIRST, IINTCR
            IF (ABS(ABS(CRDIN1(I))-VAL) .LT. 1.0D-8) THEN
               CRDIN1(I) = 1.1D6
            END IF
 30      CONTINUE
      END IF
      IF (IC+5 .LE. IINTCR) THEN
         IC = IC + 10
         GOTO 25
      END IF
C
C     Then all marked coordinates are removed.
C
      IC = IINTCR - NDIHD + 1
 35   CONTINUE
      IF (IC .LE. IINTCR) THEN
         IF (CRDIN1(IC) .GT. 1.0D6) THEN
            DO 40 I = IC, IINTCR-1
               DO 42 J = 1, 6
                  INTCRD(I,J) = INTCRD(I+1,J)
 42            CONTINUE
               CRDIN1(I) = CRDIN1(I+1)
 40         CONTINUE
            IINTCR = IINTCR - 1
            GOTO 35
         ELSE
            IC = IC + 1
            GOTO 35
         END IF
      END IF
      NDIHD = IINTCR - IFIRST + 1
      IF ((NDIHD .GT. 10) .AND. (NDIHD .GT. IINTCR/5)) THEN
         IC = IINTCR - NDIHD + 1
         GOTO 20
      END IF
      IF (IPRINT .GE. IPRMIN) THEN
         WRITE(LUPRI,'(/I10,A/)') IBEFOR-IINTCR,
     &        ' dihedral angles were removed to reduce redundancy!'
      END IF
      RETURN
      END

C  /* Deck getwil */
      SUBROUTINE GETWIL(MXRCRD,MX2CRD,ATMARR,WILBMT,BMTRAN,TMPMAT)
C
C     Constructs the Wilson B matrix used for transformations between
C     Cartesian and redundant natural internal coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
#include "pi.h"
      DIMENSION ATMARR(MXCENT,8), WILBMT(MXRCRD,MXCOOR)
      DIMENSION BMTRAN(MXRCRD,MXRCRD), TMPMAT(MX2CRD,MX2CRD)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3), VEC6(3)
      CHARACTER TXT*32
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (DEG175 = 175D0*PI/180D0)
C
      CALL DZERO(WILBMT,MXRCRD*MXCOOR)
      CALL DZERO(ATMARR,8*MXCENT)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from GETWIL','*',103)
      END IF
C
C     First we initialize the ATMARR array.
C
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Expanded array of atoms (Bohr)',-1)
         CALL OUTPUT(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      ICRTCR = 3*IATOM
C
      IF (REBILD) THEN
         IINTCR = IREDIC
         IREDIC = -1
      END IF
      NRIC = IINTCR
      IF (DELINT .AND. (IREDIC .GT. 0)) NRIC = IREDIC
      DO 20 IC = 1, NRIC
C
C     The connection between Cartesian coordinates and bonds:
C     -------------------------------------------------------
C
         IF (INTCRD(IC,1) .LT. 10) THEN
            VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
            VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
            VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
            CALL NRMLVC(VEC1)
            DO 25 II = 1, 3
               WILBMT(IC,(INTCRD(IC,2)-1)*3+II) =  VEC1(II)
               WILBMT(IC,(INTCRD(IC,3)-1)*3+II) = -VEC1(II)
 25         CONTINUE
C
            IF (IPRINT .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
               VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
               VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
               BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results for coordinate',IC
               WRITE(LUPRI,*)
               DO 27 II = 1, IATOM
                  DO 28 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     BNDL2 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     CHANGE = (BNDL2-BNDL1)/1.0D-6
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0D-6
                     WRITE(LUPRI,'(2I5,F16.6)') II,JJ,CHANGE
 28               CONTINUE
 27            CONTINUE
            END IF
C
C     The connection between Cartesian coordinates and angles:
C     --------------------------------------------------------
C
         ELSE IF (INTCRD(IC,1) .LT. 20) THEN
            VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
            VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
            VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
            VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
            VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
            VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
            BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
            BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
            CALL NRMLVC(VEC1)
            CALL NRMLVC(VEC2)
C
C     All regular angles
C
            IF (INTCRD(IC,1) .EQ. 11) THEN
               CALL VECPRD(VEC1,VEC2,VEC3)
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0D-8) THEN
                  VEC3(1) =  VEC1(2)+VEC1(3)
                  VEC3(2) = -VEC1(1)+VEC1(3)
                  VEC3(3) = -VEC1(1)-VEC1(2)
               END IF
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0D-8) THEN
                  VEC3(1) =  VEC1(2)-VEC1(3)
                  VEC3(2) = -VEC1(1)-VEC1(3)
                  VEC3(3) =  VEC1(1)+VEC1(2)
               END IF
               CALL NRMLVC(VEC3)
C
C     Second coordinate of angles larger than 175 degrees.
C
            ELSE
               CALL VECPRD(VEC1,VEC3,VEC4)
               CALL NRMLVC(VEC4)
               VEC3(1) = VEC4(1)
               VEC3(2) = VEC4(2)
               VEC3(3) = VEC4(3)
            END IF
C
            CALL VECPRD(VEC1,VEC3,VEC4)
            CALL VECPRD(VEC3,VEC2,VEC5)
            CALL NRMLVC(VEC4)
            CALL NRMLVC(VEC5)
            DO 30 II = 1, 3
               WILBMT(IC,(INTCRD(IC,2)-1)*3+II) = VEC4(II)/BNDL1
               WILBMT(IC,(INTCRD(IC,3)-1)*3+II) =
     &                          -VEC4(II)/BNDL1 - VEC5(II)/BNDL2
               WILBMT(IC,(INTCRD(IC,4)-1)*3+II) = VEC5(II)/BNDL2
 30         CONTINUE
C
            IF (IPRINT .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
               VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
               VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
               VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
               VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
               VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
               CALL VECPRD(VEC3,VEC1,VEC4)
               CALL NRMLVC(VEC4)
               ANGLE = VECANG(VEC1,VEC2)
               IF ((INTCRD(IC,1) .EQ. 12) .OR. (INTCRD(IC+1,1) .EQ. 12))
     &              ANGLE = VECANG(VEC1,VEC4) + VECANG(VEC4,VEC2)
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results for coordinate',IC
               WRITE(LUPRI,*)
               DO 47 II = 1, IATOM
                  DO 48 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     ANGLE2 = VECANG(VEC1,VEC2)
                     IF ((INTCRD(IC,1) .EQ. 12) .OR.
     &                  (INTCRD(IC+1,1) .EQ. 12))
     &                    ANGLE2 = VECANG(VEC1,VEC4) + VECANG(VEC4,VEC2)
                     CHANGE = (ANGLE2-ANGLE)/1.0D-6
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0D-6
                     WRITE(LUPRI,'(2I5,F16.6)') II,JJ,CHANGE
 48               CONTINUE
 47            CONTINUE
            END IF
C
C     The connection between Cartesian coordinates and dihedral angles:
C     -----------------------------------------------------------------
C
         ELSE
            VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
            VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
            VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
            VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
            VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
            VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
            VEC3(1) = ATMARR(INTCRD(IC,5),2) - ATMARR(INTCRD(IC,4),2)
            VEC3(2) = ATMARR(INTCRD(IC,5),3) - ATMARR(INTCRD(IC,4),3)
            VEC3(3) = ATMARR(INTCRD(IC,5),4) - ATMARR(INTCRD(IC,4),4)
            BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
            BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
            BNDL3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
            CALL NRMLVC(VEC1)
            CALL NRMLVC(VEC2)
            CALL NRMLVC(VEC3)
            UW =  DDOT(3,VEC1,1,VEC2,1)
            VW =  DDOT(3,VEC3,1,VEC2,1)
            CALL VECPRD(VEC1,VEC2,VEC4)
            CALL VECPRD(VEC3,VEC2,VEC5)
C
C     Check for linearity. The dihedral angle is undefined if three
C     or four of the atoms are linear.
C
            IF ((DDOT(3,VEC4,1,VEC4,1) .GT. 1.0D-16) .AND.
     &           (DDOT(3,VEC5,1,VEC5,1) .GT. 1.0D-16)) THEN
               DO 50 II = 1, 3
C
C     We construct a few "building blocks".
C
                  CMPU1 = VEC4(II)/(BNDL1*(1.0D0-UW*UW))
                  CMPV1 = VEC5(II)/(BNDL3*(1.0D0-VW*VW))
                  CMPU2 = VEC4(II)*UW/(BNDL2*(1.0D0-UW*UW))
                  CMPV2 = VEC5(II)*VW/(BNDL2*(1.0D0-VW*VW))
                  WILBMT(IC,(INTCRD(IC,2)-1)*3+II) =  CMPU1
                  WILBMT(IC,(INTCRD(IC,3)-1)*3+II) = -CMPU1+CMPU2-CMPV2
                  WILBMT(IC,(INTCRD(IC,4)-1)*3+II) =  CMPV1+CMPV2-CMPU2
                  WILBMT(IC,(INTCRD(IC,5)-1)*3+II) = -CMPV1
 50            CONTINUE
            END IF
C
            IF (IPRINT .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
               VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
               VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
               VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
               VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
               VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
               VEC3(1) = ATMARR(INTCRD(IC,5),2) - ATMARR(INTCRD(IC,4),2)
               VEC3(2) = ATMARR(INTCRD(IC,5),3) - ATMARR(INTCRD(IC,4),3)
               VEC3(3) = ATMARR(INTCRD(IC,5),4) - ATMARR(INTCRD(IC,4),4)
               CALL NRMLVC(VEC2)
               CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
               VEC1(1) = VEC1(1) - CMPNT1*VEC2(1)
               VEC1(2) = VEC1(2) - CMPNT1*VEC2(2)
               VEC1(3) = VEC1(3) - CMPNT1*VEC2(3)
               CALL NRMLVC(VEC1)
               CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
               VEC3(1) = VEC3(1) - CMPNT2*VEC2(1)
               VEC3(2) = VEC3(2) - CMPNT2*VEC2(2)
               VEC3(3) = VEC3(3) - CMPNT2*VEC2(3)
               CALL NRMLVC(VEC3)
               DIHED = VECANG(VEC1,VEC3)
               IF (ABS(DIHED) .GT. DEG175) THEN
                  CALL VECPRD(VEC1,VEC2,VEC4)
                  CALL NRMLVC(VEC4)
                  DIHED = VECANG(VEC1,VEC4)+VECANG(VEC4,VEC3)
               END IF
               VEC5(1) = VEC1(1)
               VEC5(2) = VEC1(2)
               VEC5(3) = VEC1(3)
               VEC6(1) = VEC3(1)
               VEC6(2) = VEC3(2)
               VEC6(3) = VEC3(3)
               WRITE(LUPRI,'(/A,I6/)')
     &         'Numerical results for dihedral coordinate:',IC
               FAC = -D1
               DO 97 II = 1, IATOM
                  DO 98 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC3(1) = ATMARR(INTCRD(IC,5),2)
     &                    - ATMARR(INTCRD(IC,4),2)
                     VEC3(2) = ATMARR(INTCRD(IC,5),3)
     &                    - ATMARR(INTCRD(IC,4),3)
                     VEC3(3) = ATMARR(INTCRD(IC,5),4)
     &                    - ATMARR(INTCRD(IC,4),4)
                     CALL NRMLVC(VEC2)
                     CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
                     VEC1(1) = VEC1(1) - CMPNT1*VEC2(1)
                     VEC1(2) = VEC1(2) - CMPNT1*VEC2(2)
                     VEC1(3) = VEC1(3) - CMPNT1*VEC2(3)
                     CALL NRMLVC(VEC1)
                     CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
                     VEC3(1) = VEC3(1) - CMPNT2*VEC2(1)
                     VEC3(2) = VEC3(2) - CMPNT2*VEC2(2)
                     VEC3(3) = VEC3(3) - CMPNT2*VEC2(3)
                     CALL NRMLVC(VEC3)
                     IF (FAC .LT. D0) THEN
                        FAC = D1
                        CALL VECPRD(VEC1,VEC3,VEC4)
                        IF (DDOT(3,VEC4,1,VEC2,1) .LT. D0) FAC = -D1
                     END IF
                     DIHED2 = VECANG(VEC1,VEC3)
                     CH1 = FAC*(VECANG(VEC1,VEC6)
     &                    -VECANG(VEC5,VEC6))/1.0D-6
                     CH2 = FAC*(VECANG(VEC3,VEC5)
     &                    -VECANG(VEC6,VEC5))/1.0D-6
                     IF (ABS(DIHED2) .GT. DEG175) THEN
                        CALL VECPRD(VEC1,VEC2,VEC4)
                        CALL NRMLVC(VEC4)
                        DIHED2 = VECANG(VEC1,VEC4)+VECANG(VEC4,VEC3)
                     END IF
                     CHANGE = FAC*(DIHED2-DIHED)/1.0D-6
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0D-6
                     WRITE(LUPRI,'(2I5,3F16.6)') II,JJ,CHANGE,CH1,CH2
 98               CONTINUE
 97            CONTINUE
            END IF
         END IF
C
         IF (IPRINT .GE. IPRMAX) THEN
            TXT = 'Internal coordinate number:     '
            WRITE(TXT(28:32),'(I5)') IC
            CALL HEADER(TXT,-1)
            DO 99 K = 1, IATOM
               WRITE(LUPRI,'(I10,3(5X,A,F16.6))')
     &              K,'x',WILBMT(IC,(K-1)*3+1),
     &                'y',WILBMT(IC,(K-1)*3+2),
     &                'z',WILBMT(IC,(K-1)*3+3)
 99         CONTINUE
         END IF
 20   CONTINUE
C
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Wilson B matrix [B(ij) = dq(i)/dx(j)]',-1)
         CALL OUTPUT(WILBMT,1,NRIC,1,ICRTCR,MXRCRD,MXCOOR,1,LUPRI)
      END IF
C
C     If we're using delocalized internal coordinates, we make
C     the transformation to the active space.
C
      IF (DELINT .AND. (IREDIC .GT. 0)) THEN
         CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
         DO 100 I = 1, NRIC
            DO 102 J = 1, ICRTCR
               TMPMAT(I,J) = WILBMT(I,J)
 102        CONTINUE
 100     CONTINUE
         CALL DZERO(WILBMT,MXRCRD*MXCOOR)
         DO 120 I = 1, IINTCR
            DO 122 J = 1, ICRTCR
               DO 124 K = 1, IREDIC
                  WILBMT(I,J) = WILBMT(I,J) + BMTRAN(K,I)*TMPMAT(K,J)
 124           CONTINUE
 122        CONTINUE
 120     CONTINUE
         IF (IPRINT .GE. IPRMIN) THEN
            CALL HEADER('Active space B matrix',-1)
            CALL OUTPUT(WILBMT,1,IINTCR,1,ICRTCR,
     &           MXRCRD,MXCOOR,1,LUPRI)
         END IF
      END IF
      RETURN
      END

C  /* Deck getdwl */
      SUBROUTINE GETDWL(MXRCRD,DWILBM,ATMARR,TMPMAT,WILBMT)
C
C     Constructs the derivative of the Wilson B matrix used for transformations
C     between Cartesian and redundant natural internal coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
      DIMENSION DWILBM(MXCOOR,MXCOOR), WILBMT(MXRCRD,MXCOOR)
      DIMENSION ATMARR(MXCENT,8), TMPMAT(MXCOOR,MXCOOR)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3), VEC6(3)
      CHARACTER TXT*32
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      CALL DZERO(ATMARR,8*MXCENT)
      IF (IPRINT .GE. IPRMAX) THEN
         CALL TITLER('Output from GETDWL','*',103)
      END IF
C
C     First we initialize the ATMARR array.
C
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Expanded array of atoms (Bohr)',-1)
         CALL OUTPUT(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      ICRTCR = 3*IATOM
C
      DO 20 IC = 1, IINTCR
         CALL DZERO(TMPMAT,MXCOOR*MXCOOR)
         CALL DZERO(DWILBM,MXCOOR*MXCOOR)
         IF (INTCRD(IC,1) .LT. 10) THEN
C
C     Numerically:
C
            IF (IPRINT .GE. IPRDBG) THEN
               DO 30 II = 2, 3
                  DO 35 JJ = 1, 3
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) + 1.0D-6
                     VEC1(1)= ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2)= ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3)= ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     CALL NRMLVC(VEC1)
                     CALL DZERO(VEC2,3)
                     DO 40 J = 1, 3
                        VEC2(J) = 1.0D0
                        COMPNT = DDOT(3,VEC1,1,VEC2,1)
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,3)-1)*3+J) = -COMPNT
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,2)-1)*3+J) =  COMPNT
                        VEC2(J) = 0.0D0
 40                  CONTINUE
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) - 1.0D-6
 35               CONTINUE
 30            CONTINUE
               DO 50 II = 2, 3
                  DO 51 JJ = 1, 3
                     DO 53 KK = 2, 3
                        DO 54 LL = 1, 3
                           DWILBM((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) =
     &                          (TMPMAT((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) -
     &                          WILBMT(IC,(INTCRD(IC,II)-1)*3+JJ))
     &                          /1.0D-6
 54                     CONTINUE
 53                  CONTINUE
 51               CONTINUE
 50            CONTINUE
               DO 56 II = 1, ICRTCR
                  DO 58 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0D0
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0D-6)
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 58               CONTINUE
 56            CONTINUE
               CALL HEADER('Numerically diff. Wilson B matrix',-1)
               WRITE(LUPRI,*) 'Internal coordinate no.',IC
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
C     Analytically:
C
            CALL GTDWLM(MXRCRD,IC,DWILBM,ATMARR,WILBMT)
            IF (IPRINT .GE. IPRMAX) THEN
               CALL HEADER('Analytically diff. Wilson B matrix',-1)
               WRITE(LUPRI,*) 'Internal coordinate no.',IC
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
         ELSE IF (INTCRD(IC,1) .LT. 20) THEN
C
C     Numerically:
C
            IF (IPRINT .GE. IPRDBG) THEN
               DO 60 II = 2, 4
                  DO 65 JJ = 1, 3
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
                     CALL NRMLVC(VEC1)
                     CALL NRMLVC(VEC2)
C
C     All regular angles
C
                     IF (INTCRD(IC,1) .EQ. 11) THEN
                        CALL VECPRD(VEC1,VEC2,VEC3)
                        VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
                        IF (VNRM .LE. 1.0D-8) THEN
                           VEC3(1) =  VEC1(2)+VEC1(3)
                           VEC3(2) = -VEC1(1)+VEC1(3)
                           VEC3(3) = -VEC1(1)-VEC1(2)
                        END IF
                        CALL NRMLVC(VEC3)
C
C     Second coordinate of angles larger than 175 degrees.
C
                     ELSE
                        CALL VECPRD(VEC1,VEC3,VEC4)
                        CALL NRMLVC(VEC4)
                        VEC3(1) = VEC4(1)
                        VEC3(2) = VEC4(2)
                        VEC3(3) = VEC4(3)
                     END IF
C
                     CALL VECPRD(VEC1,VEC3,VEC4)
                     CALL NRMLVC(VEC4)
                     CALL VECPRD(VEC3,VEC2,VEC5)
                     CALL NRMLVC(VEC5)
C
                     DO 70 J = 1, 3
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,2)-1)*3+J) =  (VEC4(J)/BNDL1)
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,3)-1)*3+J) = -(VEC4(J)/BNDL1)
     &                       - (VEC5(J)/BNDL2)
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,4)-1)*3+J) =  (VEC5(J)/BNDL2)
 70                  CONTINUE
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) - 1.0D-6
 65               CONTINUE
 60            CONTINUE
C
               DO 90 II = 2, 4
                  DO 91 JJ = 1, 3
                     DO 93 KK = 2, 4
                        DO 94 LL = 1, 3
                           DWILBM((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) =
     &                          (TMPMAT((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) -
     &                          WILBMT(IC,(INTCRD(IC,II)-1)*3+JJ))
     &                          /1.0D-6
 94                     CONTINUE
 93                  CONTINUE
 91               CONTINUE
 90            CONTINUE
               DO 96 II = 1, ICRTCR
                  DO 98 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0D0
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0D-6)
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 98               CONTINUE
 96            CONTINUE
               CALL HEADER('Numerically diff. Wilson B matrix',-1)
               WRITE(LUPRI,*) 'Internal coordinate no.',IC
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
C     Analytically:
C
            CALL GTDWLM(MXRCRD,IC,DWILBM,ATMARR,WILBMT)
            IF (IPRINT .GE. IPRMAX) THEN
               CALL HEADER('Analytically diff. Wilson B matrix',-1)
               WRITE(LUPRI,*) 'Internal coordinate no.',IC
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
         ELSE
C
C     Numerically:
C
            IF (IPRINT .GE. IPRDBG) THEN
               DO 100 II = 2, 5
                  DO 105 JJ = 1, 3
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC3(1) = ATMARR(INTCRD(IC,5),2)
     &                    - ATMARR(INTCRD(IC,4),2)
                     VEC3(2) = ATMARR(INTCRD(IC,5),3)
     &                    - ATMARR(INTCRD(IC,4),3)
                     VEC3(3) = ATMARR(INTCRD(IC,5),4)
     &                    - ATMARR(INTCRD(IC,4),4)
                     VNRM1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     VNRM2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
                     VNRM3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
                     CALL NRMLVC(VEC1)
                     CALL NRMLVC(VEC2)
                     CALL NRMLVC(VEC3)
                     UW =  DDOT(3,VEC1,1,VEC2,1)
                     VW =  DDOT(3,VEC3,1,VEC2,1)
                     CALL VECPRD(VEC1,VEC2,VEC4)
                     CALL VECPRD(VEC2,VEC3,VEC5)
                     CALL VECPRD(VEC3,VEC2,VEC6)
C
C     Check for linearity. The dihedral angle is undefined if three
C     or four of the atoms are linear.
C
                     IF ((DDOT(3,VEC4,1,VEC4,1) .GT. 1.0D-16) .AND.
     &                    (DDOT(3,VEC5,1,VEC5,1) .GT. 1.0D-16)) THEN
                        DO 120 J = 1, 3
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,2)-1)*3+J)
     &                          = VEC4(J)/(VNRM1*(1.0D0-UW*UW))
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,3)-1)*3+J) =
     &                          (VEC4(J)*(VNRM1*UW-VNRM2))/
     &                          (VNRM1*VNRM2*(1.0D0-UW*UW))
     &                          - (VW*VEC6(J))/(VNRM2*(1.0D0-VW*VW))
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,4)-1)*3+J) =
     &                          (VEC6(J)*(VNRM3*VW+VNRM2))/
     &                          (VNRM3*VNRM2*(1.0D0-VW*VW))
     &                          - (UW*VEC4(J))/(VNRM2*(1.0D0-UW*UW))
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,5)-1)*3+J) =
     &                          VEC5(J)/(VNRM3*(1.0D0-VW*VW))
 120                    CONTINUE
                     END IF
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) - 1.0D-6
 105              CONTINUE
 100           CONTINUE
C
               DO 200 II = 2, 5
                  DO 210 JJ = 1, 3
                     DO 220 KK = 2, 5
                        DO 230 LL = 1, 3
                           DWILBM((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) =
     &                          (TMPMAT((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) -
     &                          WILBMT(IC,(INTCRD(IC,II)-1)*3+JJ))
     &                          /1.0D-6
 230                    CONTINUE
 220                 CONTINUE
 210              CONTINUE
 200           CONTINUE
               DO 250 II = 1, ICRTCR
                  DO 255 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0D0
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0D-6)
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 255              CONTINUE
 250           CONTINUE
               CALL HEADER('Numerically diff. Wilson B matrix',-1)
               WRITE(LUPRI,*) 'Internal coordinate no.',IC
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
C     Analytically:
C
            CALL GTDWLM(MXRCRD,IC,DWILBM,ATMARR,WILBMT)
            IF (IPRINT .GE. IPRMAX) THEN
               CALL HEADER('Analytically diff. Wilson B matrix',-1)
               WRITE(LUPRI,*) 'Internal coordinate no.',IC
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
         END IF
 20   CONTINUE
      RETURN
      END

C  /* Deck gtdwl0 */
      SUBROUTINE GTDWL0(MXRCRD,IC,DWILBM,ATMARR,WILBMT,BMTRAN,
     &     WORK,LWORK)
C
C     Ths subroutine determines the coordinate system in work, and
C     allocates memory if needed.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
      DIMENSION DWILBM(MXCOOR,MXCOOR), ATMARR(MXCENT,8)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION WORK(LWORK)
      IF (.NOT. DELINT) THEN
         CALL GTDWLM(MXRCRD,IC,DWILBM,ATMARR,WILBMT)
      ELSE
         KMTRED = 1
         KMTDEL = KMTRED + IREDIC*ICRTCR*ICRTCR
         KWRK   = KMTDEL + IINTCR*ICRTCR*ICRTCR
         IF (KWRK .GT. LWORK) CALL STOPIT('GTDWL0',' ',KWRK,LWORK)
         CALL GTDWL1(MXRCRD,IC,DWILBM,ATMARR,WILBMT,BMTRAN,
     &        WORK(KMTRED),WORK(KMTDEL),IREDIC,IINTCR,ICRTCR)
      END IF
      RETURN
      END

C  /* Deck gtdwl1 */
      SUBROUTINE GTDWL1(MXRCRD,IC,DWILBM,ATMARR,WILBMT,BMTRAN,REDMAT,
     &     DELMAT,IRED,IDEL,ICRT)
C
C     Ths subroutine acts as a buffer for the subroutine GTDWLM when
C     delocalized internal coordinates are used.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
      DIMENSION DWILBM(MXCOOR,MXCOOR), ATMARR(MXCENT,8)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION REDMAT(IRED,ICRT,ICRT)
      DIMENSION DELMAT(IDEL,ICRT,ICRT)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C     We only calculate the redundant derivatives once.
C
      IF (IC .EQ. 1) THEN
         CALL DZERO(REDMAT,IRED*ICRT*ICRT)
         CALL DZERO(DELMAT,IDEL*ICRT*ICRT)
         DO 10 IR = 1, IRED
            CALL GTDWLM(MXRCRD,IR,DWILBM,ATMARR,WILBMT)
            DO 20 I = 1, ICRT
               DO 22 J = 1, ICRT
                  REDMAT(IR,I,J) = DWILBM(I,J)
 22            CONTINUE
 20         CONTINUE
 10      CONTINUE
C
C     Then we transform it to delocalized coordinates.
C
         DO 30 I1 = 1, IDEL
            DO 32 J1 = 1, IRED
               DO 34 K1 = 1, ICRT
                  DO 36 K2 = 1, ICRT
                     DELMAT(I1,K1,K2) = DELMAT(I1,K1,K2)
     &                    + BMTRAN(J1,I1)*REDMAT(J1,K1,K2)
 36               CONTINUE
 34            CONTINUE
 32         CONTINUE
 30      CONTINUE
      END IF
      CALL DZERO(DWILBM,MXCOOR*MXCOOR)
      DO 50 I = 1, ICRT
         DO 52 J = 1, ICRT
            DWILBM(I,J) = DELMAT(IC,I,J)
 52      CONTINUE
 50   CONTINUE
      RETURN
      END

C  /* Deck gtdwlm */
      SUBROUTINE GTDWLM(MXRCRD,IC,DWILBM,ATMARR,WILBMT)
C
C     Returns the derivative of the Wilson B matrix for _one_
C     redundant internal coordinate (IC).
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
      DIMENSION DWILBM(MXCOOR,MXCOOR), ATMARR(MXCENT,8)
      DIMENSION WILBMT(MXRCRD,MXCOOR)
      DIMENSION VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC5(3)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
      CALL DZERO(DWILBM,MXCOOR*MXCOOR)
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
C
C     The derivative for bonds:
C     -------------------------
C
      IF (INTCRD(IC,1) .LT. 10) THEN
         VEC1(1)=ATMARR(INTCRD(IC,2),2)-ATMARR(INTCRD(IC,3),2)
         VEC1(2)=ATMARR(INTCRD(IC,2),3)-ATMARR(INTCRD(IC,3),3)
         VEC1(3)=ATMARR(INTCRD(IC,2),4)-ATMARR(INTCRD(IC,3),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         CALL NRMLVC(VEC1)
         DO 10 II = 1, 3
            DO 15 JJ = 1, 3
C
C     The variable FACIJ makes the distinction between "diagonal" and
C     "off-diagonal" elements, where "diagonal" refers to Cartesian
C     displacements (II = JJ) and not the atoms.
C
               FACIJ = D0
               IF (II .EQ. JJ) FACIJ = D1
               CMPUU = (FACIJ - VEC1(II)*VEC1(JJ))/BNDL1
               DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,2)-1)*3+JJ)
     &              =  CMPUU
               DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &              = -CMPUU
               DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &          = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
               DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &              =  CMPUU
 15         CONTINUE
 10      CONTINUE
C
C     The derivative for angles:
C     --------------------------
C
      ELSE IF (INTCRD(IC,1) .LT. 20) THEN
         VEC1(1)=ATMARR(INTCRD(IC,2),2)-ATMARR(INTCRD(IC,3),2)
         VEC1(2)=ATMARR(INTCRD(IC,2),3)-ATMARR(INTCRD(IC,3),3)
         VEC1(3)=ATMARR(INTCRD(IC,2),4)-ATMARR(INTCRD(IC,3),4)
         VEC2(1)=ATMARR(INTCRD(IC,4),2)-ATMARR(INTCRD(IC,3),2)
         VEC2(2)=ATMARR(INTCRD(IC,4),3)-ATMARR(INTCRD(IC,3),3)
         VEC2(3)=ATMARR(INTCRD(IC,4),4)-ATMARR(INTCRD(IC,3),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
         CALL NRMLVC(VEC1)
         CALL NRMLVC(VEC2)
         COSA = DDOT(3,VEC1,1,VEC2,1)
         SINA = SQRT(1-COSA*COSA)
C
C     We can only calculate the derivative, if sin(q) is non-zero, as
C     all expressions contains a division by sin(q).
C
         IF (ABS(SINA) .GT. 1.0D-4) THEN
            DO 25 II = 1, 3
               DADMI = WILBMT(IC,(INTCRD(IC,2)-1)*3+II)
               DADNI = WILBMT(IC,(INTCRD(IC,4)-1)*3+II)
               DADOI = WILBMT(IC,(INTCRD(IC,3)-1)*3+II)
               DO 27 JJ = 1, 3
                  DADMJ = WILBMT(IC,(INTCRD(IC,2)-1)*3+JJ)
                  DADNJ = WILBMT(IC,(INTCRD(IC,4)-1)*3+JJ)
                  DADOJ = WILBMT(IC,(INTCRD(IC,3)-1)*3+JJ)
C
C     The variable FACIJ makes the distinction between "diagonal" and
C     "off-diagonal" elements, where "diagonal" refers to Cartesian
C     displacements (II = JJ) and not the atoms.
C
                  FACIJ = D0
                  IF (II .EQ. JJ) FACIJ = D1
C
C     We only need a few "building blocks" to make all expressions.
C
                  CMPUU  = (VEC1(II)*VEC2(JJ) + VEC1(JJ)*VEC2(II)
     &                 - 3.0D0*COSA*VEC1(II)*VEC1(JJ)
     &                 + FACIJ*COSA)/(SINA*BNDL1*BNDL1)
                  CMPVV  = (VEC2(II)*VEC1(JJ) + VEC2(JJ)*VEC1(II)
     &                 - 3.0D0*COSA*VEC2(II)*VEC2(JJ)
     &                 + FACIJ*COSA)/(SINA*BNDL2*BNDL2)
                  CMPWW1 = (VEC1(II)*VEC1(JJ) + VEC2(II)*VEC2(JJ)
     &                 - COSA*VEC1(II)*VEC2(JJ)
     &                 - FACIJ)/(SINA*BNDL1*BNDL2)
                  CMPWW2 = (VEC2(II)*VEC2(JJ) + VEC1(II)*VEC1(JJ)
     &                 - COSA*VEC2(II)*VEC1(JJ)
     &                 - FACIJ)/(SINA*BNDL1*BNDL2)
C
C     All the different elements are constructed.
C
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,2)-1)*3+JJ)
     &             = CMPUU - DADMI*DADMJ*COSA/SINA
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &             = -CMPUU - CMPWW1 - DADMI*DADOJ*COSA/SINA
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &             = CMPWW1 - DADMI*DADNJ*COSA/SINA
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &             = CMPUU + CMPVV + CMPWW1 + CMPWW2
     &             - DADOI*DADOJ*COSA/SINA
                  DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &             = -CMPVV - CMPWW2 - DADNI*DADOJ*COSA/SINA
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,4)-1)*3+II)
     &             = DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &             = CMPVV - DADNI*DADNJ*COSA/SINA
 27            CONTINUE
 25         CONTINUE
C
C     If the second derivative is undefined (sin(q) ~ 0), we set the
C     matrix equal to the zero matrix.
C
         ELSE
            DO 30 I = 1, ICRTCR
               DO 32 J = 1, ICRTCR
                  DWILBM(I,J) = D0
 32            CONTINUE
 30         CONTINUE
         END IF
C
C     The derivative for dihedral angles:
C     -----------------------------------
C
      ELSE 
         VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
         VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
         VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
         VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
         VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
         VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
         VEC3(1) = ATMARR(INTCRD(IC,5),2) - ATMARR(INTCRD(IC,4),2)
         VEC3(2) = ATMARR(INTCRD(IC,5),3) - ATMARR(INTCRD(IC,4),3)
         VEC3(3) = ATMARR(INTCRD(IC,5),4) - ATMARR(INTCRD(IC,4),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
         BNDL3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
         CALL NRMLVC(VEC1)
         CALL NRMLVC(VEC2)
         CALL NRMLVC(VEC3)
         CALL VECPRD(VEC1,VEC2,VEC4)
         CALL VECPRD(VEC3,VEC2,VEC5)
         UW   = DDOT(3,VEC1,1,VEC2,1)
         VW   = DDOT(3,VEC3,1,VEC2,1)
         UWVW = DDOT(3,VEC4,1,VEC5,1)
         UVW  = DDOT(3,VEC1,1,VEC5,1)
C
C     We can only calculate the derivative, if sin(q) (or (u * (v x w)))
C     is non-zero, as all expressions contains a division by sin(q).
C
         IF (ABS(UVW) .GT. 1.0D-6) THEN
            DO 40 II = 1, 3
               DADMI = WILBMT(IC,(INTCRD(IC,2)-1)*3+II)
               DADNI = WILBMT(IC,(INTCRD(IC,5)-1)*3+II)
               DO 45 JJ = 1, 3
                  DADMJ = WILBMT(IC,(INTCRD(IC,2)-1)*3+JJ)
                  DADNJ = WILBMT(IC,(INTCRD(IC,5)-1)*3+JJ)
C
C     The variable FACIJ makes the distinction between "diagonal" and
C     "off-diagonal" elements, where "diagonal" refers to Cartesian
C     displacements (II = JJ) and not the atoms. The sign of this
C     factor follows this table:
C
C                                 i
C                              1  2  3
C                             ---------
C                           1| 0  +  -
C                         j 2| -  0  +
C                           3| +  -  0
C
                  FACIJ = (JJ-II)*(-0.5D0)**(ABS(JJ-II))
C
C     For some of the elements we need the third Cartesian direction,
C     that is k =/ i,j.
C
                  KK = MAX(1,MIN(6-II-JJ,3))
C
C     With the following "building blocks" all expressions can be constructed:
C
                  CMPUU  = (VEC4(II)*(UW*VEC2(JJ)-VEC1(JJ))
     &                   + VEC4(JJ)*(UW*VEC2(II)-VEC1(II)))/
     &                     (BNDL1*BNDL1*(1.0D0-UW*UW)**2)
                  CMPVV  = (VEC5(II)*(VW*VEC2(JJ)-VEC3(JJ))
     &                   + VEC5(JJ)*(VW*VEC2(II)-VEC3(II)))/
     &                     (BNDL3*BNDL3*(1.0D0-VW*VW)**2)
                  CMPUW  = 0.5D0*(VEC4(II)*(VEC2(JJ)-2.0D0*UW*VEC1(JJ)
     &                   + UW*UW*VEC2(JJ))+VEC4(JJ)*(VEC2(II)
     &                   - 2.0D0*UW*VEC1(II)+UW*UW*VEC2(II)))/
     &                     (BNDL1*BNDL2*(1.0D0-UW*UW)**2)
                  CMPVW  = 0.5D0*(VEC5(II)*(VEC2(JJ)-2.0D0*VW*VEC3(JJ)
     &                   + VW*VW*VEC2(JJ))+VEC5(JJ)*(VEC2(II)
     &                   - 2.0D0*VW*VEC3(II)+VW*VW*VEC2(II)))/
     &                     (BNDL3*BNDL2*(1.0D0-VW*VW)**2)
                  CMPWW1 = 0.5D0*(VEC4(II)*(VEC1(JJ)+UW*UW*VEC1(JJ)
     &                   - 3.0D0*UW*VEC2(JJ)+UW*UW*UW*VEC2(JJ))
     &                   + VEC4(JJ)*(VEC1(II)+UW*UW*VEC1(II)
     &                   - 3.0D0*UW*VEC2(II)+UW*UW*UW*VEC2(II)))/
     &                     (BNDL2*BNDL2*(1.0D0-UW*UW)**2)
                  CMPWW2 = 0.5D0*(VEC5(II)*(VEC3(JJ)+VW*VW*VEC3(JJ)
     &                   - 3.0D0*VW*VEC2(JJ)+VW*VW*VW*VEC2(JJ))
     &                   + VEC5(JJ)*(VEC3(II)+VW*VW*VEC3(II)
     &                   - 3.0D0*VW*VEC2(II)+VW*VW*VW*VEC2(II)))/
     &                     (BNDL2*BNDL2*(1.0D0-VW*VW)**2)
                  CMPKK1 = (UW*VEC2(KK)-VEC1(KK))/
     &                     (BNDL1*BNDL2*(1.0D0-UW*UW))
                  CMPKK2 = (VW*VEC2(KK)-VEC3(KK))/
     &                     (BNDL3*BNDL2*(1.0D0-VW*VW))
C
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,2)-1)*3+JJ)
     &                 = CMPUU
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &                 = -CMPUU + CMPUW + FACIJ*CMPKK1
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = -CMPUW - FACIJ*CMPKK1
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,5)-1)*3+JJ)
     &                 = -CMPVV
                  DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &                 = -CMPVW - FACIJ*CMPKK2
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,5)-1)*3+II)
     &             = DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = CMPVV + CMPVW + FACIJ*CMPKK2
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,5)-1)*3+II)
     &             = DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &                 = CMPUU - 2.0D0*CMPUW - CMPWW1 + CMPWW2
                  DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = CMPUW + CMPVW + CMPWW1 - CMPWW2
     &                 + FACIJ*(CMPKK1-CMPKK2)
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,3)-1)*3+II)
     &             = DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = -CMPVV - 2.0D0*CMPVW + CMPWW2 - CMPWW1
 45            CONTINUE
 40         CONTINUE
C
C     If the second derivative is undefined (sin(q) ~ 0), we set the
C     matrix equal to the identity matrix.
C
         ELSE
            DO 50 I = 1, ICRTCR
               DWILBM(I,I) = D1
 50         CONTINUE
         END IF
      END IF
      RETURN
      END

C  /* Deck gtbinv */
      SUBROUTINE GTBINV(MXRCRD,UMAT,WMAT,VMAT,VAL,WILBMT,BMTRAN,
     &     BMTINV,PJINMT,VMATRED,TMPMAT,WRK,LWRK)
C
C     Decomposes the transpose of Wilson's B matrix, then one can
C     construct the inverse of the rectangular matrix.
C
C     7/5-2009 Optimization of code, mostly by replacing matrix multiplies 
C              with dgemm (and by introducing VMATRED). /Simen Reine
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
      DIMENSION UMAT(MXCOOR,MXRCRD), WMAT(MXRCRD,MXRCRD)
      DIMENSION VMAT(MXRCRD,MXRCRD), VAL(MXRCRD)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION BMTINV(MXRCRD,MXCOOR), PJINMT(MXRCRD,MXRCRD)
      DIMENSION VMATRED(MXRCRD,MXRCRD), TMPMAT(MXRCRD,MXRCRD)
      DIMENSION WRK(LWRK)
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      IF (IPRINT .GE. IPRMED) CALL TITLER('Output from GTBINV','*',103)
      IDIM1 = ICRTCR
      IDIM2 = IINTCR
#ifdef VAR_LINSCA
C     New, more efficient verison of DECOMP
C     VMATRED is used to store VMATTR in DECOMP temporarily
      CALL DECOMP(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,WILBMT,TMPMAT,
     &     VMATRED,VAL,WRK,LWRK)
#else
C     The old version of DECOMP
      CALL DCMPOL(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,VAL,WILBMT,
     &     BMTINV)
#endif
C
C     Some output
C
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('U matrix in GTBINV',-1)
         CALL OUTPUT(UMAT,1,IDIM1,1,IDIM2,MXCOOR,MXRCRD,1,LUPRI)
         CALL HEADER('W matrix in GTBINV',-1)
         CALL OUTPUT(WMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
         CALL HEADER('V matrix in GTBINV',-1)
         CALL OUTPUT(VMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
      END IF
C
C     Based on the singular value decomposition of the B matrix,
C     we can construct a non-redundant B matrix mixing the various
C     internal primitives. The eigenvectors are stored in BMTRAN.
C
C     Not optimized time-wise! /SR
C
      IF (DELINT .AND. (IREDIC .LT. 0)) THEN
C
C     The redundant B matrix is copied to BMTRAN.
C
         CALL DZERO(BMTRAN,MXRCRD*MXRCRD)
         DO 540 I = 1, IDIM2
            DO 545 J = 1, IDIM1
               BMTRAN(I,J) = WILBMT(I,J)
 545        CONTINUE
 540     CONTINUE
C
C     We remove all vectors in V corresponding to singular values.
C
Chjaaj   WMAT(IDIM2+1,IDIM2+1) = 1.1D6
Chjaaj-Oct07: bugfix: this goes out of bounds!!!!!
         IDIM3 = IDIM2
         IC = 1
 550     CONTINUE
         IF (ABS(WMAT(IC,IC)) .LT. 1.0D-8) THEN
            DO 560 I = IC, IDIM3 - 1
               DO 565 J = 1, IDIM2
                  VMAT(J,I) = VMAT(J,I+1)
 565           CONTINUE
               DO 567 J = 1, IDIM1
                  UMAT(J,I) = UMAT(J,I+1)
 567           CONTINUE
               WMAT(I,I) = WMAT(I+1,I+1)
 560        CONTINUE
Chjaaj-Oct07: new code follows next 5 lines
            IF (IDIM3 .EQ. IDIM2) THEN
               WMAT(IDIM3,IDIM3) = 1.1D6
            ELSE
               WMAT(IDIM3,IDIM3) = WMAT(IDIM3+1,IDIM3+1)
            END IF
            IDIM3 = IDIM3 - 1
            GOTO 550
         ELSE IF (ABS(WMAT(IC,IC)) .LT. 1.0D6) THEN
            IC = IC + 1
            IF (IC .LE. IDIM2) GOTO 550
         END IF
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Reduced U matrix',-1)
            CALL OUTPUT(UMAT,1,IDIM1,1,IDIM3,MXCOOR,MXRCRD,1,LUPRI)
            CALL HEADER('Reduced W matrix',-1)
            CALL OUTPUT(WMAT,1,IDIM3,1,IDIM3,MXRCRD,MXRCRD,1,LUPRI)
            CALL HEADER('Reduced V matrix (non-redund. eigenvectors)'
     &           ,-1)
            CALL OUTPUT(VMAT,1,IDIM2,1,IDIM3,MXRCRD,MXRCRD,1,LUPRI)
         END IF
C
C     We calculate the weights of the different primitives
C
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Primitive weights in final space',-1)
            DO 600 IPRIM = 1, IDIM2
               WG = 0.0D0
               DO 610 I = 1, IDIM3
                  WG = WG + VMAT(IPRIM,I)*VMAT(IPRIM,I)
 610           CONTINUE
               WRITE(LUPRI,'(A,I4,F24.6)') '        ',IPRIM,WG
 600        CONTINUE
         END IF
C
C     We construct the B matrix for the delocalized
C     internal coordinates.
C
         CALL DZERO(WILBMT,MXRCRD*MXCOOR)
         DO 620 I = 1, IDIM3
            DO 622 J = 1, IDIM1
               DO 624 K = 1, IDIM2
                  WILBMT(I,J) = WILBMT(I,J) + VMAT(K,I)*BMTRAN(K,J)
 624           CONTINUE
 622        CONTINUE
 620     CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('B matrix of active space',-1)
            CALL OUTPUT(WILBMT,1,IDIM3,1,IDIM1,MXRCRD,MXCOOR,1,LUPRI)
         END IF
C
C     The eigenvectors (V matrix) is stored in BMTRAN
C
         CALL DZERO(BMTRAN,MXRCRD*MXRCRD)
         DO 630 I = 1, IDIM2
            DO 635 J = 1, IDIM3
               BMTRAN(I,J) = VMAT(I,J)
 635        CONTINUE
 630     CONTINUE
         IREDIC = IINTCR
         IINTCR = IDIM3
         IDIM2 = IDIM3
C
C     Then we start the decomposition again with this matrix
C     (which now should have no singular values).
C
#ifdef VAR_LINSCA
         CALL DECOMP(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,WILBMT,TMPMAT,
     &        VMATRED,VAL,WRK,LWRK)
#else
         CALL DCMPOL(MXRCRD,IDIM1,IDIM3,UMAT,WMAT,VMAT,VAL,
     &        WILBMT,BMTINV)
#endif
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Non-redundant U matrix in GTBINV',-1)
            CALL OUTPUT(UMAT,1,IDIM1,1,IDIM2,MXCOOR,MXRCRD,1,LUPRI)
            CALL HEADER('Non-redundant W matrix in GTBINV',-1)
            CALL OUTPUT(WMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
            CALL HEADER('Non-redundant V matrix in GTBINV',-1)
            CALL OUTPUT(VMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
         END IF
      END IF
C
C     The inverse is constructed:
C
C           (B^t)^-1 = V (1/w_ii) U^t
C
C     All non-negative values of W are inverted, values close to zero
C     are set to zero.
C
      NPROJ = 0
      DO 130 I = 1, IDIM2
         IF (ABS(WMAT(I,I)) .LE. 1.0D-10) THEN
            WMAT(I,I) = D0
         ELSE
            WMAT(I,I) = 1.0D0/WMAT(I,I)
            NPROJ = NPROJ + 1
         END IF
 130  CONTINUE
C
C     A projection operator that will be used later,
C     is also constructed: P = V W W^-1 V^t
C
C     Set up reduced V
      NK = 0
      DO K=1,IDIM2
        IF (ABS(WMAT(K,K)) .GE. 1.0D-10) THEN
          NK = NK + 1
          DO I=1,IDIM2
            VMATRED(I,NK) = VMAT(I,K)
          ENDDO
        ENDIF
      ENDDO
C     Calculate P
      CALL DZERO(PJINMT,MXRCRD*MXRCRD)
      CALL DGEMM('N','T',IDIM2,IDIM2,NK,1.0d0,VMATRED,MXRCRD,
     &           VMATRED,MXRCRD,0.0d0,PJINMT,MXRCRD)

      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Projection matrix',-1)
         CALL OUTPUT(PJINMT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
      END IF
C
C     We construct (1/w_ii) U^t
C
      CALL DGEMM('N','T',IDIM2,IDIM1,IDIM2,1.0d0,WMAT,MXRCRD,
     &           UMAT,MXCOOR,0.0d0,BMTINV,MXRCRD)
C
C     Then V (1/w_ii) U^t
C
      CALL DGEMM('N','N',IDIM2,IDIM1,IDIM2,1.0d0,VMAT,MXRCRD,
     &           BMTINV,MXRCRD,0.0d0,WMAT,MXRCRD)

      CALL DZERO(BMTINV,MXRCRD*MXCOOR)
      DO 180 J = 1, IDIM1
         DO 182 I = 1, IDIM2
            BMTINV(I,J) = WMAT(I,J)
 182     CONTINUE
 180  CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Inverse of B^t',-1)
         CALL OUTPUT(BMTINV,1,IDIM2,1,IDIM1,MXRCRD,MXCOOR,1,LUPRI)
      END IF
      IF (IPRINT .GE. IPRMAX) THEN
         WRITE(LUPRI,'(//A,I8/A,I8/)')
     &      ' Total number of internal coordinates:',IINTCR,
     &      ' Number of non-redundant coordinates :',NPROJ
      END IF
C
C     Put number of coordinates to be projected away in NPROJ
C
      NPROJ = IINTCR - NPROJ
      RETURN
      END

#ifdef VAR_LINSCA
C  /* Deck decomp */
      SUBROUTINE DECOMP(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,WILBMT,
     &                  WBMTTR,VMATTR,SVAL,WRK,LWRK)
C
C     Decomposes the transpose of Wilson's B matrix.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
      DIMENSION UMAT(MXCOOR,MXRCRD), WMAT(MXRCRD,MXRCRD)
      DIMENSION VMAT(MXRCRD,MXRCRD)
      DIMENSION WILBMT(MXRCRD,MXCOOR),WRK(LWRK)
      DIMENSION WBMTTR(MXCOOR,MXRCRD)
      DIMENSION VMATTR(MXRCRD,MXRCRD)
      DIMENSION SVAL(MXCOOR)
      PARAMETER (D0 = 0.0D0)
C
C     The transpose of Wilsons B matrix is decomposed to:
C
C           B^t = U W V^t  ;   U,V orthogonal, W diagonal
C
      CALL DZERO(UMAT,MXRCRD*MXCOOR)
      CALL DZERO(WMAT,MXRCRD*MXRCRD)
      CALL DZERO(VMAT,MXRCRD*MXRCRD)
      CALL DZERO(WBMTTR,MXCOOR*MXRCRD)
      CALL DZERO(SVAL,MXCOOR)
C
C     The transpose is placed in WBMTTR
C
      DO I = 1, IDIM2
         DO J = 1, IDIM1
            WBMTTR(J,I) = WILBMT(I,J)
         ENDDO
      ENDDO
C
C     The transpose of Wilsons B matrix is decomposed to:
C
C           B^t = U W V^t  ;   U,V orthogonal, W diagonal
C
      CALL DGESVD('A','A',IDIM1,IDIM2,WBMTTR,MXCOOR,SVAL,UMAT,MXCOOR,
     &            VMATTR,MXRCRD,WRK,LWRK,INFO)
      IF (INFO.NE.0) THEN
        WRITE(*,*) 'Error in DECOMP! Info = ',INFO
        CALL QUIT('Error in DECOMP!')
      ENDIF
      DO I=1,IDIM2
        WMAT(I,I) = SVAL(I)
      ENDDO
      DO I=1,IDIM2
        DO J=1,IDIM2
          VMAT(J,I) = VMATTR(I,J)
        ENDDO
      ENDDO
C
C     The transpose is placed in WBMTTR
C
      DO I = 1, IDIM2
         DO J = 1, IDIM1
            UMAT(J,I) = WILBMT(I,J)
         ENDDO
      ENDDO
      BMTNRM = D0
      SCL    = D0
      SGNFAC = D0
      DO 20 I = 1, IDIM2
         SVAL(I) = SCL*SGNFAC
         SCL    = D0
         SGNFAC = D0
         SUM    = D0
         II = I + 1
         IF (I .LE. IDIM1) THEN
            DO 25 J = I, IDIM1
               SCL = SCL + ABS(UMAT(J,I))
 25         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 30 J = I, IDIM1
                  UMAT(J,I) = UMAT(J,I)/SCL
                  SUM = SUM + UMAT(J,I)*UMAT(J,I)
 30            CONTINUE
               DIAG = UMAT(I,I)
               SGNFAC = -SIGN(SQRT(SUM),DIAG)
               FAC = SGNFAC*DIAG - SUM
               UMAT(I,I) = DIAG - SGNFAC
               DO 35 J = II, IDIM2
                  SUM = D0
                  DO 37 K = I, IDIM1
                     SUM = SUM + UMAT(K,I)*UMAT(K,J)
 37               CONTINUE
                  FAC2 = SUM/FAC
                  DO 39 K = I, IDIM1
                     UMAT(K,J) = UMAT(K,J)+FAC2*UMAT(K,I)
 39               CONTINUE
 35            CONTINUE
               DO 40 J = I, IDIM1
                  UMAT(J,I) = SCL*UMAT(J,I)
 40            CONTINUE
            END IF
         END IF
         WMAT(I,I) = SCL*SGNFAC
         SUM    = D0
         SGNFAC = D0
         SCL    = D0
         IF ((I .LE. IDIM1) .AND. (I .NE. IDIM2)) THEN
            DO 45 J = II, IDIM2
               SCL = SCL + ABS(UMAT(I,J))
 45         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 50 J = II, IDIM2
                  UMAT(I,J) = UMAT(I,J)/SCL
                  SUM = SUM + UMAT(I,J)*UMAT(I,J)
 50            CONTINUE
               FAC = UMAT(I,II)
               SGNFAC = -SIGN(SQRT(SUM),FAC)
               FAC2 = SGNFAC*FAC - SUM
               UMAT(I,II) = FAC - SGNFAC
               DO 52 J = II, IDIM2
                  SVAL(J) = UMAT(I,J)/FAC2
 52            CONTINUE
               DO 54 J = II, IDIM1
                  SUM = D0
                  DO 56 K = I+1, IDIM2
                     SUM = SUM + UMAT(J,K)*UMAT(I,K)
 56               CONTINUE
                  DO 58 K = II, IDIM2
                     UMAT(J,K) = UMAT(J,K) + SUM*SVAL(K)
 58               CONTINUE
 54            CONTINUE
               DO 60 J = II, IDIM2
                  UMAT(I,J) = SCL*UMAT(I,J)
 60            CONTINUE
            END IF
         END IF
         BMTNRM = MAX(BMTNRM,(ABS(WMAT(I,I)) + ABS(SVAL(I))))
 20   CONTINUE
      DO 70 I = IDIM2, 1, -1
         IF (I .LT. IDIM2) THEN
            IF (SGNFAC .NE. D0) THEN
               DO 72 J = II, IDIM2
                  VMAT(J,I) = UMAT(I,J)/(SGNFAC*UMAT(I,II))
 72            CONTINUE
               DO 74 J = II, IDIM2
                  SUM = D0
                  DO 76 K = II, IDIM2
                     SUM = SUM + UMAT(I,K)*VMAT(K,J)
 76               CONTINUE
                  DO 78 K = II, IDIM2
                     VMAT(K,J) = VMAT(K,J) + SUM*VMAT(K,I)
 78               CONTINUE
 74            CONTINUE
            END IF
            DO 79 J = II, IDIM2
               VMAT(I,J) = D0
               VMAT(J,I) = D0
 79         CONTINUE
         END IF
         VMAT(I,I) = 1.0D0
         SGNFAC = SVAL(I)
         II = I
 70   CONTINUE
      DO 80 I = MIN(IDIM1,IDIM2), 1, -1
         II = I + 1
         SGNFAC = WMAT(I,I)
         DO 82 J = II, IDIM2
            UMAT(I,J) = D0
 82      CONTINUE
         IF (SGNFAC .NE. D0) THEN
            SGNFAC = 1.0D0/SGNFAC
            DO 84 J = II, IDIM2
               SUM = D0
               DO 86 K = II, IDIM1
                  SUM = SUM + UMAT(K,I)*UMAT(K,J)
 86            CONTINUE
               FAC = SGNFAC*(SUM/UMAT(I,I))
               DO 88 K = I, IDIM1
                  UMAT(K,J) = UMAT(K,J) + FAC*UMAT(K,I)
 88            CONTINUE
 84         CONTINUE
            DO 90 J = I, IDIM1
               UMAT(J,I) = SGNFAC*UMAT(J,I)
 90         CONTINUE
         ELSE
            DO 92 J = I, IDIM1
               UMAT(J,I) = D0
 92         CONTINUE
         END IF
         UMAT(I,I) = UMAT(I,I) + 1.0D0
 80   CONTINUE
      DO 100 J = IDIM2, 1, -1
         DO 105 ITRS = 1, 50
            DO 110 K = J, 1, -1
               KK = K - 1
               IF ((ABS(SVAL(K))+BMTNRM) .EQ. BMTNRM) GOTO 902
               IF ((ABS(WMAT(KK,KK))+BMTNRM) .EQ. BMTNRM) GOTO 901
 110        CONTINUE
 901        CONTINUE
            FC   = D0
            SUM = 1.0D0
            DO 112  I = K, J
               FAC = SUM*SVAL(I)
               SVAL(I) = FC*SVAL(I)
               IF ((ABS(FAC)+BMTNRM) .EQ. BMTNRM) GOTO 902
               SGNFAC = WMAT(I,I)
               FAC2 = SQRT(FAC*FAC + SGNFAC*SGNFAC)
               WMAT(I,I) = FAC2
               FAC2 = 1.0D0/FAC2
               FC  = SGNFAC*FAC2
               SUM =   -FAC*FAC2
               DO 114 L = 1, IDIM1
                  FY = UMAT(L,KK)
                  FZ = UMAT(L,I)
                  UMAT(L,KK) = FY*FC + FZ*SUM
                  UMAT(L,I)  = FZ*FC - FY*SUM
 114           CONTINUE
 112        CONTINUE
 902        CONTINUE
            FZ = WMAT(J,J)
            IF (K .EQ. J) THEN
               IF (FZ .LT. D0) THEN
                  WMAT(J,J) = -FZ
                  DO 116 L = 1, IDIM2
                     VMAT(L,J) = -VMAT(L,J)
 116              CONTINUE
               END IF
               GOTO 903
            END IF
            IF (ITRS .EQ. 50) CALL QUIT('No convergence in GTBINV.')
            FX = WMAT(K,K)
            JJ = J - 1
            FY = WMAT(JJ,JJ)
            SGNFAC = SVAL(JJ)
            FAC2 = SVAL(J)
            FAC = ((FY-FZ)*(FY+FZ) + (SGNFAC-FAC2)*(SGNFAC+FAC2))/
     &           (2.0D0*FAC2*FY)
            SGNFAC = SQRT(FAC*FAC + 1.0D0)
            FAC = ((FX-FZ)*(FX+FZ)+FAC2*((FY/(FAC+SIGN(SGNFAC,FAC)))
     &           -FAC2))/FX
            FC  = 1.0D0
            SUM = 1.0D0
            DO 120 L = K, JJ
               I = L + 1
               SGNFAC = SVAL(I)
               FY = WMAT(I,I)
               FAC2   = SGNFAC*SUM
               SGNFAC = SGNFAC*FC
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               SVAL(L) = FZ
               FC  = FAC/FZ
               SUM = FAC2/FZ
               FAC = (FX*FC) + (SGNFAC*SUM)
               SGNFAC = -(FX*SUM) + (SGNFAC*FC)
               FAC2 = FY*SUM
               FY = FY*FC
               DO 122 M = 1, IDIM2
                  FX = VMAT(M,L)
                  FZ = VMAT(M,I)
                  VMAT(M,L) = (FX*FC) + (FZ*SUM)
                  VMAT(M,I) =-(FX*SUM) + (FZ*FC)
 122           CONTINUE
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               WMAT(L,L) = FZ
               IF (FZ .NE. D0) THEN
                  FZ = 1.0D0/FZ
                  FC = FAC*FZ
                  SUM = FAC2*FZ
               END IF
               FAC = (FC*SGNFAC) + (SUM*FY)
               FX = -(SUM*SGNFAC) + (FC*FY)
               DO 124 M = 1, IDIM1
                  FY = UMAT(M,L)
                  FZ = UMAT(M,I)
                  UMAT(M,L) = (FY*FC) + (FZ*SUM)
                  UMAT(M,I) = -(FY*SUM) + (FZ*FC)
 124           CONTINUE
 120        CONTINUE
            SVAL(K) = D0
            SVAL(J) = FAC
            WMAT(J,J) = FX
 105     CONTINUE
 903     CONTINUE
 100  CONTINUE
      RETURN
      END
#endif

C  /* Deck dcmpol */
      SUBROUTINE DCMPOL(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,VAL,
     &     WILBMT,BMTINV)
C
C     Old version of subroutine to decompose the transpose of Wilson's B matrix.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
      DIMENSION UMAT(MXCOOR,MXRCRD), WMAT(MXRCRD,MXRCRD)
      DIMENSION VMAT(MXRCRD,MXRCRD), VAL(MXRCRD)
      DIMENSION WILBMT(MXRCRD,MXCOOR)
      DIMENSION BMTINV(MXRCRD,MXCOOR)
      PARAMETER (D0 = 0.0D0)
C
      CALL DZERO(BMTINV,MXRCRD*MXCOOR)
      CALL DZERO(UMAT,MXRCRD*MXCOOR)
      CALL DZERO(WMAT,MXRCRD*MXRCRD)
      CALL DZERO(VMAT,MXRCRD*MXRCRD)
      CALL DZERO(VAL,MXRCRD)
C
C     The transpose is placed in UMAT
C
      DO 10 I = 1, IDIM2
         DO 12 J = 1, IDIM1
            UMAT(J,I) = WILBMT(I,J)
 12      CONTINUE
 10   CONTINUE
C
C     The transpose of Wilsons B matrix is decomposed to:
C
C           B^t = U W V^t  ;   U,V orthogonal, W diagonal
C
      BMTNRM = D0
      SCL    = D0
      SGNFAC = D0
      DO 20 I = 1, IDIM2
         VAL(I) = SCL*SGNFAC
         SCL    = D0
         SGNFAC = D0
         SUM    = D0
         II = I + 1
         IF (I .LE. IDIM1) THEN
            DO 25 J = I, IDIM1
               SCL = SCL + ABS(UMAT(J,I))
 25         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 30 J = I, IDIM1
                  UMAT(J,I) = UMAT(J,I)/SCL
                  SUM = SUM + UMAT(J,I)*UMAT(J,I)
 30            CONTINUE
               DIAG = UMAT(I,I)
               SGNFAC = -SIGN(SQRT(SUM),DIAG)
               FAC = SGNFAC*DIAG - SUM
               UMAT(I,I) = DIAG - SGNFAC
               DO 35 J = II, IDIM2
                  SUM = D0
                  DO 37 K = I, IDIM1
                     SUM = SUM + UMAT(K,I)*UMAT(K,J)
 37               CONTINUE
                  FAC2 = SUM/FAC
                  DO 39 K = I, IDIM1
                     UMAT(K,J) = UMAT(K,J)+FAC2*UMAT(K,I)
 39               CONTINUE
 35            CONTINUE
               DO 40 J = I, IDIM1
                  UMAT(J,I) = SCL*UMAT(J,I)
 40            CONTINUE
            END IF
         END IF
         WMAT(I,I) = SCL*SGNFAC
         SUM    = D0
         SGNFAC = D0
         SCL    = D0
         IF ((I .LE. IDIM1) .AND. (I .NE. IDIM2)) THEN
            DO 45 J = II, IDIM2
               SCL = SCL + ABS(UMAT(I,J))
 45         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 50 J = II, IDIM2
                  UMAT(I,J) = UMAT(I,J)/SCL
                  SUM = SUM + UMAT(I,J)*UMAT(I,J)
 50            CONTINUE
               FAC = UMAT(I,II)
               SGNFAC = -SIGN(SQRT(SUM),FAC)
               FAC2 = SGNFAC*FAC - SUM
               UMAT(I,II) = FAC - SGNFAC
               DO 52 J = II, IDIM2
                  VAL(J) = UMAT(I,J)/FAC2
 52            CONTINUE
               DO 54 J = II, IDIM1
                  SUM = D0
                  DO 56 K = I+1, IDIM2
                     SUM = SUM + UMAT(J,K)*UMAT(I,K)
 56               CONTINUE
                  DO 58 K = II, IDIM2
                     UMAT(J,K) = UMAT(J,K) + SUM*VAL(K)
 58               CONTINUE
 54            CONTINUE
               DO 60 J = II, IDIM2
                  UMAT(I,J) = SCL*UMAT(I,J)
 60            CONTINUE
            END IF
         END IF
         BMTNRM = MAX(BMTNRM,(ABS(WMAT(I,I)) + ABS(VAL(I))))
 20   CONTINUE
      DO 70 I = IDIM2, 1, -1
         IF (I .LT. IDIM2) THEN
            IF (SGNFAC .NE. D0) THEN
               DO 72 J = II, IDIM2
                  VMAT(J,I) = UMAT(I,J)/(SGNFAC*UMAT(I,II))
 72            CONTINUE
               DO 74 J = II, IDIM2
                  SUM = D0
                  DO 76 K = II, IDIM2
                     SUM = SUM + UMAT(I,K)*VMAT(K,J)
 76               CONTINUE
                  DO 78 K = II, IDIM2
                     VMAT(K,J) = VMAT(K,J) + SUM*VMAT(K,I)
 78               CONTINUE
 74            CONTINUE
            END IF
            DO 79 J = II, IDIM2
               VMAT(I,J) = D0
               VMAT(J,I) = D0
 79         CONTINUE
         END IF
         VMAT(I,I) = 1.0D0
         SGNFAC = VAL(I)
         II = I
 70   CONTINUE
      DO 80 I = MIN(IDIM1,IDIM2), 1, -1
         II = I + 1
         SGNFAC = WMAT(I,I)
         DO 82 J = II, IDIM2
            UMAT(I,J) = D0
 82      CONTINUE
         IF (SGNFAC .NE. D0) THEN
            SGNFAC = 1.0D0/SGNFAC
            DO 84 J = II, IDIM2
               SUM = D0
               DO 86 K = II, IDIM1
                  SUM = SUM + UMAT(K,I)*UMAT(K,J)
 86            CONTINUE
               FAC = SGNFAC*(SUM/UMAT(I,I))
               DO 88 K = I, IDIM1
                  UMAT(K,J) = UMAT(K,J) + FAC*UMAT(K,I)
 88            CONTINUE
 84         CONTINUE
            DO 90 J = I, IDIM1
               UMAT(J,I) = SGNFAC*UMAT(J,I)
 90         CONTINUE
         ELSE
            DO 92 J = I, IDIM1
               UMAT(J,I) = D0
 92         CONTINUE
         END IF
         UMAT(I,I) = UMAT(I,I) + 1.0D0
 80   CONTINUE
      DO 100 J = IDIM2, 1, -1
         DO 105 ITRS = 1, 50
            DO 110 K = J, 1, -1
               KK = K - 1
               IF ((ABS(VAL(K))+BMTNRM) .EQ. BMTNRM) GOTO 902
               IF ((ABS(WMAT(KK,KK))+BMTNRM) .EQ. BMTNRM) GOTO 901
 110        CONTINUE
 901        CONTINUE
            FC   = D0
            SUM = 1.0D0
            DO 112  I = K, J
               FAC = SUM*VAL(I)
               VAL(I) = FC*VAL(I)
               IF ((ABS(FAC)+BMTNRM) .EQ. BMTNRM) GOTO 902
               SGNFAC = WMAT(I,I)
               FAC2 = SQRT(FAC*FAC + SGNFAC*SGNFAC)
               WMAT(I,I) = FAC2
               FAC2 = 1.0D0/FAC2
               FC  = SGNFAC*FAC2
               SUM =   -FAC*FAC2
               DO 114 L = 1, IDIM1
                  FY = UMAT(L,KK)
                  FZ = UMAT(L,I)
                  UMAT(L,KK) = FY*FC + FZ*SUM
                  UMAT(L,I)  = FZ*FC - FY*SUM
 114           CONTINUE
 112        CONTINUE
 902        CONTINUE
            FZ = WMAT(J,J)
            IF (K .EQ. J) THEN
               IF (FZ .LT. D0) THEN
                  WMAT(J,J) = -FZ
                  DO 116 L = 1, IDIM2
                     VMAT(L,J) = -VMAT(L,J)
 116              CONTINUE
               END IF
               GOTO 903
            END IF
            IF (ITRS .EQ. 50) CALL QUIT('No convergence in GTBINV.')
            FX = WMAT(K,K)
            JJ = J - 1
            FY = WMAT(JJ,JJ)
            SGNFAC = VAL(JJ)
            FAC2 = VAL(J)
            FAC = ((FY-FZ)*(FY+FZ) + (SGNFAC-FAC2)*(SGNFAC+FAC2))/
     &           (2.0D0*FAC2*FY)
            SGNFAC = SQRT(FAC*FAC + 1.0D0)
            FAC = ((FX-FZ)*(FX+FZ)+FAC2*((FY/(FAC+SIGN(SGNFAC,FAC)))
     &           -FAC2))/FX
            FC  = 1.0D0
            SUM = 1.0D0
            DO 120 L = K, JJ
               I = L + 1
               SGNFAC = VAL(I)
               FY = WMAT(I,I)
               FAC2   = SGNFAC*SUM
               SGNFAC = SGNFAC*FC
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               VAL(L) = FZ
               FC  = FAC/FZ
               SUM = FAC2/FZ
               FAC = (FX*FC) + (SGNFAC*SUM)
               SGNFAC = -(FX*SUM) + (SGNFAC*FC)
               FAC2 = FY*SUM
               FY = FY*FC
               DO 122 M = 1, IDIM2
                  FX = VMAT(M,L)
                  FZ = VMAT(M,I)
                  VMAT(M,L) = (FX*FC) + (FZ*SUM)
                  VMAT(M,I) =-(FX*SUM) + (FZ*FC)
 122           CONTINUE
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               WMAT(L,L) = FZ
               IF (FZ .NE. D0) THEN
                  FZ = 1.0D0/FZ
                  FC = FAC*FZ
                  SUM = FAC2*FZ
               END IF
               FAC = (FC*SGNFAC) + (SUM*FY)
               FX = -(SUM*SGNFAC) + (FC*FY)
               DO 124 M = 1, IDIM1
                  FY = UMAT(M,L)
                  FZ = UMAT(M,I)
                  UMAT(M,L) = (FY*FC) + (FZ*SUM)
                  UMAT(M,I) = -(FY*SUM) + (FZ*FC)
 124           CONTINUE
 120        CONTINUE
            VAL(K) = D0
            VAL(J) = FAC
            WMAT(J,J) = FX
 105     CONTINUE
 903     CONTINUE
 100  CONTINUE
      RETURN
      END

C  /* Deck gq2gx */
      SUBROUTINE GQ2GX(MXRCRD,GRADQ,GRADX,WILBMT)
C
C     Transforms the gradient in redundant internal coordinates to
C     Cartesian coordinates:
C                              g_x = B^t g_q
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
      DIMENSION GRADQ(MXRCRD), GRADX(MXCOOR), WILBMT(MXRCRD,MXCOOR)
      CALL DZERO(GRADX,MXCOOR)
      DO 10 I = 1, ICRTCR
         DO 12 J = 1, IINTCR
            GRADX(I) = GRADX(I) + WILBMT(J,I)*GRADQ(J)
 12      CONTINUE
 10   CONTINUE
      RETURN
      END

C  /* Deck gx2gq */
      SUBROUTINE GX2GQ(MXRCRD,GRADX,GRADQ,BMTINV)
C
C     Transforms the gradient in Cartesian coordinates to redundant
C     internal coordinates:
C                              g_q = (B^t)^-1 g_x
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
      DIMENSION GRADX(MXCOOR), GRADQ(MXRCRD), BMTINV(MXRCRD,MXCOOR)
      CALL DZERO(GRADQ,MXRCRD)
      DO 10 I = 1, IINTCR
         DO 12 J = 1, ICRTCR
            GRADQ(I) = GRADQ(I) + BMTINV(I,J)*GRADX(J)
 12      CONTINUE
 10   CONTINUE
      RETURN
      END

C  /* Deck hq2hx */
      SUBROUTINE HQ2HX(MXRCRD,MX2CRD,ATMARR,TMPMAT,TMPMT2,HESSQ,
     &     GRADQ,HESSX,WILBMT,BMTRAN,WORK,LWORK)
C
C     Transforms the Hessian in redundant internal coordinates to
C     Cartesian coordinates:
C                              H_x = B^t H_q B + B'^t g_q
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
      DIMENSION ATMARR(MXCENT,8), TMPMAT(MX2CRD,MX2CRD)
      DIMENSION TMPMT2(MX2CRD,MX2CRD), HESSQ(MXRCRD,MXRCRD)
      DIMENSION GRADQ(MXRCRD), HESSX(MXCOOR,MXCOOR)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION WORK(LWORK)
      CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
      DO 10 IC = 1, IINTCR
         CALL GTDWL0(MXRCRD,IC,HESSX,ATMARR,WILBMT,BMTRAN,WORK,LWORK)
         DO 12 I = 1, ICRTCR
            DO 14 J = 1, ICRTCR
               TMPMAT(I,J) = TMPMAT(I,J) + HESSX(I,J)*GRADQ(IC)
 14         CONTINUE
 12      CONTINUE
 10   CONTINUE
      CALL DZERO(TMPMT2,MX2CRD*MX2CRD)
      DO 20 I = 1, IINTCR
         DO 22 J = 1, ICRTCR
            DO 24 K = 1, IINTCR
               TMPMT2(I,J) = TMPMT2(I,J) + HESSQ(I,K)*WILBMT(K,J)
 24         CONTINUE
 22      CONTINUE
 20   CONTINUE
      CALL DZERO(HESSX,MXCOOR*MXCOOR)
      DO 30 I = 1, ICRTCR
         DO 32 J = 1, ICRTCR
            DO 34 K = 1, IINTCR
               HESSX(I,J) = HESSX(I,J) + WILBMT(K,I)*TMPMT2(K,J)
 34         CONTINUE
            HESSX(I,J) = HESSX(I,J) + TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
C
C     We make sure the Cartesian Hessian is symmetric.
C
      DO 40 I = 1, ICRTCR
         DO 42 J = 1, I
            HESSX(J,I) = HESSX(I,J)
 42      CONTINUE
 40   CONTINUE
      RETURN
      END

C  /* Deck hx2hq */
      SUBROUTINE HX2HQ(MXRCRD,MX2CRD,ATMARR,TMPMAT,TMPMT2,TMPMT3,
     &     HESSX,GRADQ,HESSQ,WILBMT,BMTINV,BMTRAN,WORK,LWORK)
C
C     Transforms the Hessian in Cartesian coordinates to redundant
C     internal coordinates:
C                         H_q = (B^t)^-1 (H_x - B'^t g_q) B^-1
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
      DIMENSION ATMARR(MXCENT,8), TMPMAT(MXCOOR,MXCOOR)
      DIMENSION TMPMT2(MX2CRD,MX2CRD), TMPMT3(MX2CRD,MX2CRD)
      DIMENSION HESSX(MXCOOR,MXCOOR), GRADQ(MXRCRD)
      DIMENSION HESSQ(MXRCRD,MXRCRD)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTINV(MXRCRD,MXCOOR)
      DIMENSION BMTRAN(MXRCRD,MXRCRD), WORK(LWORK)
      LOGICAL INCLGT
      CALL DZERO(HESSQ,MXRCRD*MXRCRD)
      CALL DZERO(TMPMAT,MXCOOR*MXCOOR)
      CALL DZERO(TMPMT2,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT3,MX2CRD*MX2CRD)
C
C     INCLGT determines wether the gradient term in the transformation
C     is included. There seems to be no reason NOT to include it!!
C     (in some articles this term has been accused of ruining the
C     eigenvalue structure of the Hessian).
C
      INCLGT = .TRUE.
C     
      IF (INCLGT) THEN
         DO 10 IC = 1, IINTCR
            CALL GTDWL0(MXRCRD,IC,TMPMAT,ATMARR,WILBMT,BMTRAN,
     &           WORK,LWORK)
            DO 12 I = 1, ICRTCR
               DO 14 J = 1, ICRTCR
                  TMPMT2(I,J) = TMPMT2(I,J) + TMPMAT(I,J)*GRADQ(IC)
 14            CONTINUE
 12         CONTINUE
 10      CONTINUE
      END IF
      DO 16 J = 1, ICRTCR
         DO 18 I = 1, ICRTCR
            TMPMT2(I,J) = HESSX(I,J) - TMPMT2(I,J)
 18      CONTINUE
 16   CONTINUE
      DO 20 I = 1, ICRTCR
         DO 22 J = 1, IINTCR
            DO 24 K = 1, ICRTCR
               TMPMT3(I,J) = TMPMT3(I,J) + TMPMT2(I,K)*BMTINV(J,K)
 24         CONTINUE
 22      CONTINUE
 20   CONTINUE
      DO 30 I = 1, IINTCR
         DO 32 J = 1, IINTCR
            DO 34 K = 1, ICRTCR
               HESSQ(I,J) = HESSQ(I,J) + BMTINV(I,K)*TMPMT3(K,J)
 34         CONTINUE
 32      CONTINUE
 30   CONTINUE
      RETURN
      END

C  /* Deck cghint */
      SUBROUTINE CGHINT(MXRCRD,MX2CRD,SCLVEC,HESTMP,TMPMAT,TMPMT2,
     &     TMPMT3,TMPMT4,WILBMT,BMTINV,BMTRAN,HESINT,WORK,LWORK)
C
C     Transforms Cartesian gradient and Hessian to redundant internal
C     coordinates.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "optinf.h"
#include "taymol.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION SCLVEC(MXCOOR), HESTMP(MXCOOR,MXCOOR)
      DIMENSION TMPMAT(MX2CRD,MX2CRD), TMPMT2(MX2CRD,MX2CRD)
      DIMENSION TMPMT3(MX2CRD,MX2CRD), TMPMT4(MX2CRD,MX2CRD)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTINV(MXRCRD,MXCOOR)
      DIMENSION BMTRAN(MXRCRD,MXRCRD)
      DIMENSION HESINT(MXRCRD,MXRCRD), WORK(LWORK)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We initialize the scaling vector (the inverse of the
C     normalization vector in WLKCGH).
C
      CALL DZERO(SCLVEC,MXCOOR)
C
C Latest version of g77 has problems with loop-unrolling here,
C loop rewritten to avoid problems. VB Oct 2003
C
      DO 1 IREP = 0, MAXREP
         DO 2 ICRD = 1, 3*NUCIND
#include "scalar.h"
            ISCOOR = IPTCNT(ICRD,IREP,1)
            IF (ISCOOR .GT. 0) THEN
               ICENT = (ICRD-1)/3 + 1
               SCLVEC(ISCOOR) = SQRT(FMULT(ISTBNU(ICENT)))
            END IF
 2       CONTINUE
 1    CONTINUE
C
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from CGHINT','*',103)
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Scaling vector',-1)
            CALL OUTPUT(SCLVEC,1,1,1,NCRTOT,1,MXCOOR,1,LUPRI)
         END IF
      END IF
      IF (MAXREP .GT. 0) THEN
         CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
         CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,NCRREP(0,1),3*NUCDEP)
         CALL GX2GQ(MXRCRD,TMPMAT,GRDINT,BMTINV)
         GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian gradient',-1)
            CALL OUTPUT(TMPMAT,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
         END IF
         CALL DZERO(HESTMP,MXCOOR*MXCOOR)
         CALL TRAHES(HESMOL,MXCOOR,HESTMP,TMPMT2,TMPMT3,MXCOOR,
     &               3*NUCDEP,1)
         CALL HX2HQ(MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,HESTMP,
     &        GRDINT,HESINT,WILBMT,BMTINV,BMTRAN,WORK,LWORK)
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian Hessian',-1)
            CALL OUTPUT(HESTMP,1,ICRTCR,1,ICRTCR,MXCOOR,MXCOOR,1,LUPRI)
         END IF
      ELSE
         CALL GX2GQ(MXRCRD,GRDMOL,GRDINT,BMTINV)
         GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian gradient',-1)
            CALL OUTPUT(GRDMOL,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
         END IF
         CALL HX2HQ(MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,HESMOL,
     &        GRDINT,HESINT,WILBMT,BMTINV,BMTRAN,WORK,LWORK)
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian Hessian',-1)
            CALL OUTPUT(HESMOL,1,ICRTCR,1,ICRTCR,MXCOOR,MXCOOR,1,LUPRI)
         END IF
      END IF
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Internal gradient',-1)
         CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
         CALL HEADER('Internal Hessian',-1)
         CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck prjint */
      SUBROUTINE PRJINT(MXRCRD,NDIM,PJINMT,CONMAT,HESINT,TMPMT1,TMPMT2,
     &     TMPMT3,TMPMT4,WORK,LWORK)
C
C     Because of possible redundancies, we have to project
C     both the gradient and Hessian. If a constrained optimization
C     has been requested, we have to construct a matrix of constraints
C     and modify the projection matrix.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
      DIMENSION PJINMT(MXRCRD,MXRCRD), CONMAT(MXRCRD,MXRCRD)
      DIMENSION HESINT(MXRCRD,MXRCRD), TMPMT1(NDIM,NDIM)
      DIMENSION TMPMT2(NDIM,NDIM), TMPMT3(NDIM,NDIM)
      DIMENSION TMPMT4(NDIM,NDIM), WORK(LWORK)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      NDIM2 = NDIM*NDIM
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from PRJINT','*',103)
      END IF
C
C     The matrix of constraints is constructed if necessary.
C
      IF (CONOPT) THEN
         CALL DZERO(CONMAT,MXRCRD*MXRCRD)
         CALL DZERO(TMPMT1,NDIM2)
         CALL DZERO(TMPMT2,NDIM2)
         CALL DZERO(TMPMT3,NDIM2)
         CALL DZERO(TMPMT4,NDIM2)
         DO 10 I = 1, NDIM
            IF (ICNSTR(I) .EQ. 1) THEN
               CONMAT(I,I) = 1.0D0
               TMPMT1(I,I) = PJINMT(I,I)
               DO 15 J = I+1, NDIM
C
C     At the same time we construct the product CP'C,
C     Where C is the matrix of constraints and P' the original
C     projection matrix.
C
                  TMPMT1(I,I) = PJINMT(I,I)
                  IF (ICNSTR(J) .EQ. 1) THEN
                     TMPMT1(I,J) = PJINMT(I,J)
                     TMPMT1(J,I) = PJINMT(I,J)
                  END IF
 15            CONTINUE
            END IF
 10      CONTINUE
C
C     We find the inverse of CP'C
C
         CALL DGEINV(NDIM,TMPMT1,TMPMT2,TMPMT3,WORK,INFO)
C
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Matrix of constraints',-1)
            CALL OUTPUT(CONMAT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
            CALL HEADER('CP''C',-1)
            CALL OUTPUT(TMPMT1,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
            CALL HEADER('Inverse of CP''C',-1)
            CALL OUTPUT(TMPMT2,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
         END IF
C
C     We find CP' and P'C
C
         CALL DZERO(TMPMT1,NDIM2)
         CALL DZERO(TMPMT3,NDIM2)
         DO 20 I = 1, NDIM
            IF (ICNSTR(I) .EQ. 1) THEN
               DO 22 J = 1, NDIM
                  DO 24 K = 1, NDIM
                     TMPMT1(I,J) = TMPMT1(I,J) + CONMAT(I,K)*PJINMT(K,J)
                     TMPMT3(J,I) = TMPMT3(J,I) + PJINMT(J,K)*CONMAT(K,I)
 24               CONTINUE
 22            CONTINUE
            END IF
 20      CONTINUE
C
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('CP''',-1)
            CALL OUTPUT(TMPMT1,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
            CALL HEADER('P''C',-1)
            CALL OUTPUT(TMPMT3,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
         END IF
C
C     We calculate (CP'C)^-1 CP'
C
         DO 30 I = 1, NDIM
            IF (ICNSTR(I) .EQ. 1) THEN
               DO 32 J = 1, NDIM
                  DO 34 K = 1, NDIM
                     TMPMT4(I,J) = TMPMT4(I,J) + TMPMT2(I,K)*TMPMT1(K,J)
 34               CONTINUE
 32            CONTINUE
            END IF
 30      CONTINUE
C
C     then we find P'C (CP'C)^-1 CP'
C
         CALL DZERO(TMPMT1,NDIM2)
         DO 40 I = 1, NDIM
            IF (ICNSTR(I) .EQ. 1) THEN
               DO 42 J = 1, NDIM
                  DO 44 K = 1, NDIM
                     TMPMT1(I,J) = TMPMT1(I,J) + TMPMT3(I,K)*TMPMT4(K,J)
 44               CONTINUE
 42            CONTINUE
            END IF
 40      CONTINUE
C
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('(CP''C)^-1 CP''',-1)
            CALL OUTPUT(TMPMT4,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
            CALL HEADER('P''C (CP''C)^-1 CP''',-1)
            CALL OUTPUT(TMPMT1,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
         END IF
C
C     The projection matrix is then modified to include the constraints.
C
         DO 50 I = 1, NDIM
            DO 52 J = 1, NDIM
               PJINMT(I,J) = PJINMT(I,J) - TMPMT1(I,J)
 52         CONTINUE
 50      CONTINUE
      END IF
C
      IF (IPRINT .GE. IPRMED) THEN
         IF (CONOPT) THEN
            CALL HEADER('Projection matrix w/constraints',-1)
         ELSE
            CALL HEADER('Projection matrix',-1)
         END IF
         CALL OUTPUT(PJINMT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Unprojected gradient',-1)
         CALL OUTPUT(GRDINT,1,1,1,NDIM,1,MXRCRD,1,LUPRI)
         CALL HEADER('Unprojected Hessian',-1)
         CALL OUTPUT(HESINT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
      END IF
C
C     First we project the gradient
C
      CALL DZERO(TMPMT1,NDIM2)
      DO 70 I = 1, NDIM
         DO 72 J = 1, NDIM
            TMPMT1(I,1) = TMPMT1(I,1) + PJINMT(I,J)*GRDINT(J)
 72      CONTINUE
 70   CONTINUE
      DO 75 I = 1, NDIM
         GRDINT(I) = TMPMT1(I,1)
 75   CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Projected gradient',-1)
         CALL OUTPUT(GRDINT,1,1,1,NDIM,1,MXRCRD,1,LUPRI)
      END IF
C
C     Then the Hessian
C
      CALL DZERO(TMPMT1,NDIM2)
      DO 80 I = 1, NDIM
         DO 82 J = 1, NDIM
            DO 84 K = 1, NDIM
               TMPMT1(I,J) = TMPMT1(I,J) + PJINMT(I,K)*HESINT(K,J)
 84         CONTINUE
 82      CONTINUE
 80   CONTINUE
      CALL DZERO(HESINT,MXRCRD*MXRCRD)
      DO 85 I = 1, NDIM
         DO 87 J = 1, NDIM
            DO 89 K = 1, NDIM
               HESINT(I,J) = HESINT(I,J) + TMPMT1(I,K)*PJINMT(K,J)
 89         CONTINUE
 87      CONTINUE
 85   CONTINUE
C
C     The projected Hessian is "stabilized".
C     1.0D4 is an arbitrary (high) value.
C     This section has been commented out because this
C     stabilization seems unnecessary, it just messes up some
C     of the eigenvalues.
C
C      DO 90 I = 1, NDIM
C         DO 95 J = 1, NDIM
C            HESINT(I,J) = HESINT(I,J) - 1.0D4*PJINMT(I,J)
C 95      CONTINUE
C         HESINT(I,I) = HESINT(I,I) + 1.0D4
C 90   CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Projected Hessian',-1)
         CALL OUTPUT(HESINT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck diaint */
      SUBROUTINE DIAINT(MXRCRD,MX2CRD,NCRDHS,EVEC,EVCTMP,HESPCK,
     &     TMPHES,TMPMAT,THRIND,HESINT,DIDENT,WRK,LWRK)
C
C     The Hessian in redundant internal coordinates is diagonalized.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
      DIMENSION EVEC(MX2CRD,MX2CRD), EVCTMP(NCRDHS,NCRDHS)
      DIMENSION HESPCK(NCRDHS*NCRDHS), TMPHES(NCRDHS,NCRDHS)
      DIMENSION TMPMAT(MX2CRD*MX2CRD), HESINT(MXRCRD,MXRCRD)
      DIMENSION DIDENT(MXRCRD,MXRCRD), WRK(LWRK)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0D0)
#ifdef OLD_DIAINT
      DO 10 J = 1, IINTCR
         DO 12 I = 1, IINTCR
            TMPHES(I,J) = HESINT(I,J)
 12      CONTINUE
 10   CONTINUE
      IF (RATFUN .AND. (.NOT. SADDLE)) THEN
         DO 15 I = 1, IINTCR
            TMPHES(I,NCRDHS) = GRDINT(I)
            TMPHES(NCRDHS,I) = GRDINT(I)
 15      CONTINUE
         TMPHES(NCRDHS,NCRDHS) = D0
      END IF
      CALL DZERO(HESPCK,NCRDHS*NCRDHS)
      CALL DSITSP(NCRDHS,TMPHES,HESPCK)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from DIAINT','*',103)
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Internal gradient',-1)
            CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
               CALL HEADER('Internal Hessian',-1)
               CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,
     &              MXRCRD,MXRCRD,1,LUPRI)
            IF (RATFUN .AND. (.NOT. SADDLE)) THEN
               CALL HEADER('Augmented internal Hessian',-1)
               CALL OUTPUT(TMPHES,1,NCRDHS,1,NCRDHS,
     &              NCRDHS,NCRDHS,1,LUPRI)
            END IF
            CALL HEADER('Packed Hessian',-1)
            CALL OUTPUT(HESPCK,1,1,1,(NCRDHS*(NCRDHS+1))/2,
     &           1,NCRDHS*NCRDHS,1,LUPRI)
         END IF
      END IF
      CALL DUNIT(EVCTMP,NCRDHS)
      CALL JACO(HESPCK,EVCTMP,NCRDHS,NCRDHS,NCRDHS,
     &     TMPMAT(1),TMPMAT(MX2CRD+1))
      INDTOT = 0
      DO 20 J = 1, NCRDHS
         EVAL(J) = HESPCK(J*(J+1)/2)
         GRDDIA(J) = DDOT(NCRDHS,GRDINT,1,EVCTMP(1,J),1)
         DO 22 I = 1, NCRDHS
            EVEC(I,J) = EVCTMP(I,J)
 22      CONTINUE
         IF (EVAL(J) .LT. -THRIND) INDTOT = INDTOT + 1
 20   CONTINUE
#else
C
C     Code for circumventing the time-demanding call to JACO. /SR 4/5-2009
C
      DO J = 1, IINTCR
         DO I = 1, IINTCR
            TMPHES(I,J) = HESINT(I,J)
         ENDDO
      ENDDO
      IF (RATFUN .AND. (.NOT. SADDLE)) THEN
         DO I = 1, IINTCR
            TMPHES(I,NCRDHS) = GRDINT(I)
            TMPHES(NCRDHS,I) = GRDINT(I)
         ENDDO
         TMPHES(NCRDHS,NCRDHS) = D0
      END IF

      CALL DUNIT(DIdent,NCRDHS)

#ifdef VAR_DEBUG
! March 2010 - debugging gfortran problem with segmentation error in
! dsygv /hjaaj
      write (lupri,*) 'ncrdhs,lwrk',ncrdhs,lwrk
      write (lupri,*) 'TMPHES'
      call output(tmphes,1,ncrdhs,1,ncrdhs,ncrdhs,ncrdhs,-1,lupri)
      write (lupri,*) 'DIDENT'
      call output(dident,1,ncrdhs,1,ncrdhs,ncrdhs,ncrdhs,-1,lupri)
      write (lupri,*) 'eval'
      call output(eval,1,1,1,ncrdhs,1,ncrdhs,-1,lupri)
      write (lupri,*) 'WRK(1:40)',WRK(1:40)
      call flshfo(lupri)
#endif
      
      CALL DSYGV(1,'V','U',NCRDHS,TMPHES,NCRDHS,DIdent,NCRDHS,EVAL,
     &          WRK,LWRK,INFO)
      INDTOT = 0
      DO J = 1, NCRDHS
         GRDDIA(J) = DDOT(NCRDHS,GRDINT,1,TMPHES(1,J),1)
         DO I = 1, NCRDHS
            EVEC(I,J) = TMPHES(I,J)
         ENDDO
         IF (EVAL(J) .LT. -THRIND) INDTOT = INDTOT + 1
      ENDDO
#endif
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Eigenvalues in DIAINT',-1)
         CALL OUTPUT(EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
         CALL HEADER('Eigenvectors in DIAINT',-1)
         CALL OUTPUT(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
         CALL HEADER('Gradient (diagonal rep.) in DIAINT',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
      END IF
C
C     The eigenvalues are sorted
C
      DO 25 I = 1, NCRDHS
         IF (ABS(EVAL(I)) .LT. 1.0D-6) EVAL(I) = EVAL(I) +  1.0D5
 25   CONTINUE
      DO 30 I = 1, NCRDHS
         JMIN = I
         EMIN = EVAL(I)
         DO 35 J = (I + 1), NCRDHS
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               JMIN = J
            END IF
 35      CONTINUE
         IF (JMIN .NE. I) THEN
            CALL DSWAP(1,EVAL(I),1,EVAL(JMIN),1)
            CALL DSWAP(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
            CALL DSWAP(1,GRDDIA(I),1,GRDDIA(JMIN),1)
         END IF
 30   CONTINUE
      IF (RATFUN .AND. (.NOT. SADDLE) .AND. (EVAL(1) .LT. -THRIND)
     &     .AND. INDTOT .GT. 0) INDTOT = INDTOT - 1
      DO 40 I = 1, NCRDHS
         IF (ABS((ABS(EVAL(I)) - 1.0D5)) .LT. 1.0D-3)
     &        EVAL(I) = EVAL(I) -  1.0D5
 40   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Sorted eigenvalues in WLKEI1',-1)
         CALL OUTPUT(EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
         CALL HEADER('Sorted eigenvectors in WLKEI1',-1)
         CALL OUTPUT(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
         CALL HEADER('Gradient (sorted) in WLKEI1',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck fnstin */
      SUBROUTINE FNSTIN(MXRCRD,MX2CRD,NCRDHS,SCLVEC,HESINT,EVEC,
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,CSTEP,WILBMT,BMTRAN,
     &     BMTINV,GRDARR,STPARR,ACTIVE,EMOD,VECMOD,STPLIN)
C
C     We determine the step in redundant internal coordinates,
C     then we find the corresponding step vector in
C     Cartesian coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "pi.h"
#include "cbiwlk.h"
      DIMENSION SCLVEC(MXCOOR)
      DIMENSION HESINT(MXRCRD,MXRCRD), EVEC(MX2CRD,MX2CRD)
      DIMENSION TMPMAT(MX2CRD*MX2CRD), TMPMT2(MX2CRD*MX2CRD)
      DIMENSION TMPMT3(MX2CRD*MX2CRD*2), TMPMT4(MX2CRD,MX2CRD)
      DIMENSION TMPMT5(MX2CRD,MX2CRD)
      DIMENSION CSTEP(MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION BMTINV(MXRCRD,MXCOOR)
      DIMENSION GRDARR(25,MXRCRD), STPARR(MXRCRD,MXRCRD)
      DIMENSION STPLIN(MXRCRD), VECMOD(MXCOOR)
      LOGICAL   INSIDE, ACTIVE, DOSCAL
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      NVEC = IINTCR - NPROJ
      IF (LNSRCH .AND. (.NOT. RATFUN) .AND. (ITRNMR .GT. 0))
     &     CALL LINSRC(IINTCR,MXRCRD,GRDINT,GRDARR(1,1),STPLIN,
     &     STPARR(1,1),TMPMAT,TMPMT2,ACTIVE,EMOD)
      IF (ACTIVE) THEN
         DO 5 J = 1, IINTCR
            DO 7 I = 1, KEPTIT
               STPARR(J,I) = STPARR(J,I) - STPLIN(J)
 7          CONTINUE
            IF (.NOT. RATFUN)
     &           GRDDIA(J) = DDOT(IINTCR,GRDINT,1,EVEC(1,J),1)
 5       CONTINUE
         IF (.NOT. RATFUN) THEN
            IF (IPRINT .GT. 5) THEN
               CALL HEADER('Diagonal interpolated gradient',1)
               CALL OUTPUT(GRDDIA,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
            END IF
         END IF
      END IF
C
      IF (IPRINT .GE. IPRMED) CALL TITLER('Output from FNSTIN','*',103)
C
C     First comes the trust region method
C
      IF (TRSTRG .OR. (GDIIS .AND. (KEPTIT .LT. 3))) THEN
C
C     For saddle point optimizations, we construct the image function.
C
         IF (SADDLE) THEN
            IMODE = NSPMOD
C
C     We can follow a specific eigenvector if needed...
C
            IF (NSPMOD .GT. 0) THEN
               CALL FNDMOD(.TRUE.,MXRCRD,EVEC,WILBMT,VECMOD,
     &              TMPMT2,TMPMT3,IMODE)
            ELSE
C
C     ... otherwise we have to find the first mode with a non-zero
C     gradient-element.
C
               IMODE = 1
 3             CONTINUE
               IF ((ABS(GRDDIA(IMODE)) .LE. 1.0D-10)
     &              .AND. (IMODE .LE. NVEC)) THEN
                  IMODE = IMODE + 1
                  GOTO 3
               ELSE IF (IMODE .GT. NVEC) THEN
                  IMODE = 1
               END IF
            END IF
            IF (IPRINT .GE. IPRMED) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Making image function by changing ' //
     &              'the sign of mode ',IMODE
               WRITE(LUPRI,*)
            END IF
            CALL MAKIMG(IINTCR,NVEC,MXRCRD,EVAL,GRDDIA,
     &           STPDIA,IMODE,.FALSE.)
         END IF
C
C     First we find the internal step in diagonal representation.
C
         CALL DZERO(STPDIA,MXRCRD)
         DO 10 I = 1, IINTCR
C
C     The eigenvalue threshold was changed from 1.0D-8 due to
C     problems with optimization in delocalized internals.
C
            IF (ABS(EVAL(I)) .LE. 1.0D-6) THEN
               STPDIA(I) = D0
            ELSE
               STPDIA(I) = -GRDDIA(I)/EVAL(I)
            END IF
 10      CONTINUE
         STPNRM = SQRT(DDOT(IINTCR,STPDIA,1,STPDIA,1))
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Newton step',-1) 
            CALL OUTPUT(STPDIA,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
         END IF
C
C     If Newton step is larger that trust radius, we take a step
C     to the boundary. If the Hessian index is larger than zero,
C     the level-shifted step will also be employed, provided the
C     Newton step is larger than 0.5D-3. For saddle points we
C     employ the level-shift when the index is different from 1.
C
         IF (((STPNRM .GT. TRSTRA) .AND. (.NOT. NOTRST)) .OR.
     &        ((.NOT. SADDLE) .AND. (INDTOT .GT. 0) .AND.
     &        (STPNRM .GE. 0.5D-3)) .OR. (SADDLE .AND.
     &        (INDTOT .NE. 1))) THEN
            IF (IPRINT .GE. IPRMED) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm of Newton step:', STPNRM
               WRITE(LUPRI,'(A,F15.10/)')' Trust radius       :', TRSTRA
            END IF
            INSIDE = .FALSE.
            IF (STPNRM .LT. TRSTRA) INSIDE = .TRUE.
            CALL LSHFT0(IINTCR,IINTCR,EVAL,GRDDIA,STPDIA,
     &           MIN(TRSTRA,STPNRM),RNU,.FALSE.,ZERGRD,INSIDE,IPRINT)
            BNDNRM = TRSTRA
            STPNRM = SQRT(DDOT(IINTCR,STPDIA,1,STPDIA,1))
            IF (IPRINT .GE. IPRMED) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm, boundary step:', STPNRM
            END IF
         END IF
C
C     For saddle point optimizations, we restore the original function.
C
         IF (SADDLE) THEN
            DOSCAL = (.NOT. NEWTON)
            IF (INITHS .AND. (ITRNMR .EQ. 0)) DOSCAL = .FALSE.
            IF (DOSCAL) THEN
               CMPLIM = MAX(TRSTRA*0.67D0, 0.5D0)
               DO 35 I = 1, IINTCR
                  IF (ABS(STPDIA(I)) .GT. CMPLIM)
     &                 STPDIA(I) = SIGN(CMPLIM,STPDIA(I))
 35            CONTINUE
            END IF
            CALL MAKIMG(IINTCR,NVEC,MXRCRD,EVAL,GRDDIA,
     &           STPDIA,IMODE,.TRUE.)
         END IF
C
C     Energy is predicted, will be used later to update trust radius.
C
         ERGPRD = DDOT(IINTCR,GRDDIA,1,STPDIA,1)
     &        + 0.5D0*DV3DOT(IINTCR,STPDIA,EVAL,STPDIA)
C
C     If the predicted energy is positive, it means the Newton step is
C     towards a maximum/saddle point. We then simply reverse the
C     step direction (a bit dirty, but seems to work).
C
         IF ((.NOT. SADDLE) .AND. (ERGPRD .GT. 0.0D0)) THEN
            WRITE(LUPRI,'(A,F25.15)')
     &           ' Original pred. energy change',ERGPRD
            WRITE(LUPRI,*) 'Reversing step!'
            DO 40 I = 1, IINTCR
               STPDIA(I) = -STPDIA(I)
 40         CONTINUE
            ERGPRD = DDOT(IINTCR,GRDDIA,1,STPDIA,1)
     &           + 0.5D0*DV3DOT(IINTCR,STPDIA,EVAL,STPDIA)
            WRITE(LUPRI,'(A,F25.15)')
     &           ' New pred. energy change',ERGPRD
         END IF
         IF (IPRINT .GE. IPRMIN) THEN
            IF (IPRINT .GE. IPRDBG) THEN
               CALL HEADER('Internal diagonal gradient',-1)
               CALL OUTPUT(GRDDIA,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
               CALL HEADER('Eigenvalues',-1)
               CALL OUTPUT(EVAL,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
            END IF
            IF (IPRINT .GE. IPRMED) THEN
               CALL HEADER('Internal diagonal step',-1)
               CALL OUTPUT(STPDIA,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
            END IF
            WRITE (LUPRI,'(/A,F25.15)') ' Predicted energy change',
     &           ERGPRD
         END IF
C
C     The diagonal step is converted to ordinary internal
C     coordinates.
C
         CALL DZERO(STPINT,MXRCRD)
         DO 20 I = 1, IINTCR
            DO 22 J = 1, IINTCR
               STPINT(I) = STPINT(I) + STPDIA(J)*EVEC(I,J)
 22         CONTINUE
C
C     For angles and dihedral angles we have to avoid step components
C     giving multiples of 2*pi.
C
            IF (REDINT .AND. INTCRD(I,1) .GT. 10)
     &           STPINT(I) = MOD(STPINT(I),2.0D0*PI)
            IF (ABS(STPINT(I)) .LE. 1.0D-8) STPINT(I) = D0
 20      CONTINUE
         IF (IPRINT .GE. IPRMIN) THEN
            CALL HEADER('Step in internal coordinates',-1)
            CALL OUTPUT(STPINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
         END IF
C
C     The rational function method
C
      ELSE IF (RATFUN) THEN
         IF (SADDLE) THEN
            CALL PRFSTI(MXRCRD,MX2CRD,NCRDHS,HESINT,EVEC,
     &           TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,VECMOD)
         ELSE
            CALL RFSTP(MX2CRD,NCRDHS,MXRCRD,IINTCR,EVEC,STPINT,GRDINT,
     &           TMPMAT,HESINT)
         END IF
         IF (IPRINT .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)')
     &           ' Predicted energy change', ERGPRD
         END IF
C
C     The Geometrical DIIS method
C
      ELSE IF (GDIIS) THEN
         CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
         CALL DZERO(TMPMT4,MX2CRD*MX2CRD)
C
C     First we have to construct the inverse Hessian.
C
         DO 210 I = 1, IINTCR
            IF (ABS(EVAL(I)) .GE. 1.0D-6) THEN
               DO 212 J = 1, IINTCR
                  TMPMT4(I,J) = EVEC(J,I)/EVAL(I)
 212           CONTINUE
            END IF
 210     CONTINUE
         DO 215 I = 1, IINTCR
            DO 217 J = 1, IINTCR
               DO 219 K = 1, IINTCR
                  TMPMAT(I+(J-1)*IINTCR) = TMPMAT(I+(J-1)*IINTCR) +
     &                 EVEC(I,K)*TMPMT4(K,J)
 219           CONTINUE
 217        CONTINUE
 215     CONTINUE
C
C     Then the DIIS-step is determined
C
         CALL GDISTP(MXRCRD,IINTCR,MXRCRD,MX2CRD,STPDIA,GRDINT,HESINT,
     &        TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,GRDARR,STPARR)
         IF (IPRINT .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)')
     &           ' Predicted energy change', ERGPRD
         END IF
         DO 250 I = 1, IINTCR
            IF (ABS(STPDIA(I)) .GE. 1.0D-6) THEN
               STPINT(I) = STPDIA(I)
            ELSE
               STPINT(I) = D0
            END IF
 250     CONTINUE
      END IF
      IF (ACTIVE) THEN
         DO 300 I = 1, IINTCR
            STPINT(I) = STPINT(I) + STPLIN(I)
 300     CONTINUE
         ERGPRD = ERGPRD + (EMOD-ENERGY)
         WRITE (LUPRI,'(/A,F25.15)')
     &        ' Modified energy prediction due to line search', ERGPRD
      END IF
      STPNRM = SQRT(DDOT(IINTCR,STPINT,1,STPINT,1))
C
C     Find Cartesian step vector.
C
      CALL STPI2C(MXRCRD,SCLVEC,TMPMAT(1),TMPMAT(8*MXCENT+1),
     &     TMPMAT(16*MXCENT+1),TMPMT2(1),
     &     TMPMT2(MX2CRD+1),TMPMT2(2*MX2CRD+1),TMPMT2(3*MX2CRD+1),
     &     TMPMT2(4*MX2CRD+1),CSTEP,BMTRAN,BMTINV,TMPMT3(1),
     &     TMPMT3(1+MXCOOR*MXCOOR),.FALSE.)
C
C     If the optimization is constrained, we take a small step to
C     reimpose the constraints.
C
      IF (CONOPT) THEN
         CALL STPI2C(MXRCRD,SCLVEC,TMPMAT(1),TMPMAT(8*MXCENT+1),
     &        TMPMAT(16*MXCENT+1),TMPMT2(1),
     &        TMPMT2(MX2CRD+1),TMPMT2(2*MX2CRD+1),TMPMT2(3*MX2CRD+1),
     &        TMPMT2(4*MX2CRD+1),CSTEP,BMTRAN,BMTINV,TMPMT3(1),
     &        TMPMT3(1+MXCOOR*MXCOOR),.TRUE.)
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Internal values before step',-1)
            CALL OUTPUT(TMPMT2(1),1,1,1,IINTCR,1,MX2CRD,1,LUPRI)
         END IF
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Internal values after step',-1)
            CALL OUTPUT(TMPMT2(2*MX2CRD+1),1,1,1,IINTCR,1,
     &           MX2CRD,1,LUPRI)
         END IF
      ELSE
         IF (IPRINT .GE. IPRMAX) THEN
            NRIC = IINTCR
            CALL HEADER('Internal values before step',-1)
            CALL OUTPUT(TMPMT2(MX2CRD+1),1,1,1,IINTCR,1,
     &           MX2CRD,1,LUPRI)
         END IF
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Internal values after step',-1)
            CALL OUTPUT(TMPMT2(2*MX2CRD+1),1,1,1,IINTCR,1,
     &           MX2CRD,1,LUPRI)
         END IF
      END IF
      RETURN
      END
      
C  /* Deck getint */
      SUBROUTINE GETINT(IATOM,MXRCRD,ATMCRD,VALINT)
C
C     Determines the value of all redundant internal coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "priunit.h"
#include "pi.h"
      DIMENSION ATMCRD(MXCENT,8), VALINT(MXRCRD)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3)
      PARAMETER (D0 = 0.0D0, DEG179 = 179D0*PI/180D0)
      CALL DZERO(VALINT,MXRCRD)
      NRIC = IINTCR
      IF (DELINT) NRIC = IREDIC
      DO 10 IC = 1, NRIC
         IF (INTCRD(IC,1) .LT. 10) THEN
            VEC1(1) = ATMCRD(INTCRD(IC,2),2) - ATMCRD(INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(INTCRD(IC,2),3) - ATMCRD(INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(INTCRD(IC,2),4) - ATMCRD(INTCRD(IC,3),4)
            VALINT(IC) = SQRT(DDOT(3,VEC1,1,VEC1,1))
         ELSE IF (INTCRD(IC,1) .LT. 20) THEN
            VEC1(1) = ATMCRD(INTCRD(IC,2),2) - ATMCRD(INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(INTCRD(IC,2),3) - ATMCRD(INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(INTCRD(IC,2),4) - ATMCRD(INTCRD(IC,3),4)
            VEC2(1) = ATMCRD(INTCRD(IC,4),2) - ATMCRD(INTCRD(IC,3),2)
            VEC2(2) = ATMCRD(INTCRD(IC,4),3) - ATMCRD(INTCRD(IC,3),3)
            VEC2(3) = ATMCRD(INTCRD(IC,4),4) - ATMCRD(INTCRD(IC,3),4)
C
C     Regular angles
C
            IF (INTCRD(IC,1) .EQ. 11) THEN
               CALL VECPRD(VEC1,VEC2,VEC3)
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0D-8) THEN
                  VEC3(1) =  VEC1(2)+VEC1(3)
                  VEC3(2) = -VEC1(1)+VEC1(3)
                  VEC3(3) = -VEC1(1)-VEC1(2)
               END IF
               CALL NRMLVC(VEC3)
C
C     Second coordinate of angles larger than 175 degrees.
C
            ELSE
               CALL VECPRD(VEC1,VEC3,VEC4)
               CALL NRMLVC(VEC4)
               VEC3(1) = VEC4(1)
               VEC3(2) = VEC4(2)
               VEC3(3) = VEC4(3)
            END IF
            CALL VECPRD(VEC3,VEC1,VEC4)
            CALL NRMLVC(VEC4)
            VALINT(IC) = VECANG(VEC1,VEC2)
            IF ((INTCRD(IC,1) .EQ. 12) .OR. (INTCRD(IC+1,1) .EQ. 12))
     &           VALINT(IC) = VECANG(VEC1,VEC4) + VECANG(VEC4,VEC2)
C
C     Dihedral angles
C
         ELSE
            VEC1(1) = ATMCRD(INTCRD(IC,2),2) - ATMCRD(INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(INTCRD(IC,2),3) - ATMCRD(INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(INTCRD(IC,2),4) - ATMCRD(INTCRD(IC,3),4)
            VEC2(1) = ATMCRD(INTCRD(IC,4),2) - ATMCRD(INTCRD(IC,3),2)
            VEC2(2) = ATMCRD(INTCRD(IC,4),3) - ATMCRD(INTCRD(IC,3),3)
            VEC2(3) = ATMCRD(INTCRD(IC,4),4) - ATMCRD(INTCRD(IC,3),4)
            VEC3(1) = ATMCRD(INTCRD(IC,5),2) - ATMCRD(INTCRD(IC,4),2)
            VEC3(2) = ATMCRD(INTCRD(IC,5),3) - ATMCRD(INTCRD(IC,4),3)
            VEC3(3) = ATMCRD(INTCRD(IC,5),4) - ATMCRD(INTCRD(IC,4),4)
            CALL NRMLVC(VEC2)
            CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
            VEC1(1) = VEC1(1) - CMPNT1*VEC2(1)
            VEC1(2) = VEC1(2) - CMPNT1*VEC2(2)
            VEC1(3) = VEC1(3) - CMPNT1*VEC2(3)
            CALL NRMLVC(VEC1)
            CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
            VEC3(1) = VEC3(1) - CMPNT2*VEC2(1)
            VEC3(2) = VEC3(2) - CMPNT2*VEC2(2)
            VEC3(3) = VEC3(3) - CMPNT2*VEC2(3)
            CALL NRMLVC(VEC3)
            IF ((DDOT(3,VEC1,1,VEC1,1) .GT. 1.0D-16) .AND.
     &           (DDOT(3,VEC3,1,VEC3,1) .GT. 1.0D-16)) THEN
               VALINT(IC) = VECANG(VEC1,VEC3)
               IF (ABS(VALINT(IC)) .GT. DEG179) THEN
                  CALL VECPRD(VEC1,VEC2,VEC4)
                  CALL NRMLVC(VEC4)
                  DCMP2 = VECANG(VEC4,VEC3)
                  IF (DCMP2 .GE. PI/2.D0) DCMP2 = PI - DCMP2
                  VALINT(IC) = VECANG(VEC1,VEC4) + DCMP2
               END IF
               CALL VECPRD(VEC1,VEC3,VEC4)
               IF (DDOT(3,VEC2,1,VEC4,1) .LT. 0.0D0)
     &              VALINT(IC) = -VALINT(IC)
            ELSE
               VALINT(IC) = D0
            END IF
         END IF
 10   CONTINUE
      RETURN
      END

C  /* Deck stpi2c */
      SUBROUTINE STPI2C(MXRCRD,SCLVEC,CRDORG,CRDOLD,CRDNEW,VALORG,
     &     VALOLD,VALNEW,DIFF,TMPVEC,CSTEP,BMTRAN,BMTINV,
     &     CSTRA,SCTRA,CORREC)
C
C     Transforms internal step to Cartesian step iteratively.
C     The logical parameter CORREC indicates if the internal step
C     is a small correctional step to reimpose constraints in a
C     constrained optimization.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "optinf.h"
#include "pi.h"
      DIMENSION SCLVEC(MXCOOR), CRDORG(MXCENT,8), CRDOLD(MXCENT,8)
      DIMENSION CRDNEW(MXCENT,8), VALOLD(MXRCRD), VALORG(MXRCRD)
      DIMENSION VALNEW(MXRCRD), DIFF(MXRCRD)
      DIMENSION TMPVEC(MXCOOR), CSTEP(MXCOOR)
      DIMENSION BMTRAN(MXRCRD,MXRCRD), BMTINV(MXRCRD,MXCOOR)
      DIMENSION CSTRA(MXCOOR*MXCOOR), SCTRA(MXCOOR*MXCOOR)
      LOGICAL CORREC, ADJUST
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (ITRLIM = 25)
      CHARACTER STPTXT*16
C
      ADJUST = .FALSE.
      NRIC = IINTCR
      NNNRIC = 0
      IF (DELINT) THEN
         NRIC = IREDIC
         NNNRIC = IINTCR
      END IF
C
C     To be able to properly control the step we must use
C     the primitive (redundant) internal coordinates. If we use
C     delocalized internals, we transform the step to redundant
C     internals. VALORG is used to store the non-redundant step.
C
      IF (DELINT) THEN
         CALL DZERO(VALORG,MXRCRD)
         DO 100 I = 1, NNNRIC
            VALORG(I) = STPINT(I)
 100     CONTINUE
         CALL DZERO(STPINT,MXRCRD)
         DO 105 I = 1, NRIC
            DO 107 J = 1, NNNRIC
               STPINT(I) = STPINT(I) + BMTRAN(I,J)*VALORG(J)
 107        CONTINUE
 105     CONTINUE
      END IF
C
C     We have to find the value of all internal coordinates for the
C     old geometry. If a correctional step is requested, we copy
C     the values of the internal coordinates after the major step.
C
      IF (.NOT. CORREC) THEN
         CALL ATMINI(CRDOLD,IATOM,.TRUE.)
         CALL GETINT(IATOM,MXRCRD,CRDOLD,VALOLD)
         RMSLIM = 1.0D-6
      ELSE
         IATOM = ICRTCR / 3
C        ... BUGFIX hjaaj Oct 07, IATOM was not defined.
C            Is used for TMPVEC below, which has dimension ICRTCR,
C            based on this I divided by three to get IATOM (see DO 15 loop)
         DO 3 I = 1, NRIC
            VALORG(I) = VALOLD(I)
            VALOLD(I) = VALNEW(I)
 3       CONTINUE
         DO 5 I = 1, MXCENT
            DO 7 J = 1, 4
               CRDORG(I,J) = CRDOLD(I,J)
               CRDOLD(I,J) = CRDNEW(I,J)
 7          CONTINUE
 5       CONTINUE
         RMSLIM = 1.0D-9
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Previous coordinates',-1)
         CALL OUTPUT(CRDOLD,1,IATOM,1,4,MXCENT,8,1,LUPRI)
         IF (DELINT) THEN
            CALL HEADER('Previous primitive internal values',-1)
         ELSE
            CALL HEADER('Previous internal values',-1)
         END IF
         CALL OUTPUT(VALOLD,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
         IF (DELINT) THEN
            CALL HEADER('Step to take (primitives)',-1)
         ELSE
            CALL HEADER('Step to take',-1)
         END IF
         CALL OUTPUT(STPINT,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
      END IF
C
C     First estimate of Cartesian step.
C
      CALL DZERO(TMPVEC,MXCOOR)
      IF (.NOT. CORREC) THEN
         IF (DELINT) THEN
            DO 110 I = 1, ICRTCR
               DO 111 J = 1, NNNRIC
                  TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*VALORG(J)
 111           CONTINUE
 110         CONTINUE
          ELSE
             DO 10 I = 1, ICRTCR
               DO 11 J = 1, NRIC
                  TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*STPINT(J)
 11            CONTINUE
 10         CONTINUE
         END IF
      ELSE
C
C     SCTRA is used for temporary storage
C
         DO 12 IC = 1, NRIC
            IF (ICNSTR(IC) .GT. 0) THEN
               SCTRA(IC) = CRDIN1(IC) - VALOLD(IC)
               IF (INTCRD(IC,1) .GT. 10)
     &              SCTRA(IC) = MOD(SCTRA(IC),2.0D0*PI)
               IF ((INTCRD(IC,1) .GT. 20) .AND.
     &              (ABS(SCTRA(IC)) .GT. PI)) THEN
                  IF (SCTRA(IC) .GT. 0.0D0) THEN
                     SCTRA(IC) = SCTRA(IC) - 2.0D0*PI
                  ELSE
                     SCTRA(IC) = SCTRA(IC) + 2.0D0*PI
                  END IF
               END IF
            ELSE
               SCTRA(IC) = 0.0D0
            END IF
 12      CONTINUE
         DO 13 I = 1, ICRTCR
            DO 14 J = 1, NRIC
               TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*SCTRA(J)
 14         CONTINUE
 13      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Projected step',-1)
            CALL OUTPUT(SCTRA,1,1,1,NRIC,1,MXCOOR,1,LUPRI)
         END IF
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Inverse of B^t',-1)
         CALL OUTPUT(BMTINV,1,IINTCR,1,ICRTCR,MXRCRD,MXCOOR,1,LUPRI)
         CALL HEADER('First estimate of Cartesian step',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
      END IF
      ITRCRD = 1
      RMSOLD = D0
      RMS1ST = -D1
 123  CONTINUE
      DO 15 I = 1, IATOM
         CRDNEW(I,1) = CRDOLD(I,1)
         DO 17 J = 1, 3
            CRDNEW(I,J+1) = CRDOLD(I,J+1) + TMPVEC((I-1)*3+J)
 17      CONTINUE
 15   CONTINUE
      CALL GETINT(IATOM,MXRCRD,CRDNEW,VALNEW)
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('New coordinates',-1)
         CALL OUTPUT(CRDNEW,1,IATOM,1,4,MXCENT,8,1,LUPRI)
         IF (DELINT) THEN
            CALL HEADER('New primitive internal values',-1)
         ELSE
            CALL HEADER('New internal values',-1)
         END IF
         CALL OUTPUT(VALNEW,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
      END IF
      CALL DZERO(DIFF,MXRCRD)
      RMSINT = D0
      IF (.NOT. CORREC) THEN
         DO 18 I = 1, NRIC
            DIFF(I) = STPINT(I) - (VALNEW(I) - VALOLD(I))
            IF (INTCRD(I,1) .GT. 10)
     &           DIFF(I) = MOD(DIFF(I),2.0D0*PI)
            IF ((INTCRD(I,1) .GT. 20) .AND.
     &           (ABS(DIFF(I)) .GT. PI)) THEN
               IF (DIFF(I) .GT. 0.0D0) THEN
                  DIFF(I) = DIFF(I) - 2.0D0*PI
               ELSE
                  DIFF(I) = DIFF(I) + 2.0D0*PI
               END IF
            END IF
            IF (ABS(DIFF(I)) .LT. 1.0D-14) DIFF(I) = D0
 18      CONTINUE
         RMSINT = DDOT(NRIC,DIFF,1,DIFF,1)
      ELSE
         DO 19 I = 1, NRIC
            DIFF(I) = (CRDIN1(I) - VALNEW(I))*(ICNSTR(I)*1.0D0)
            IF (INTCRD(I,1) .GT. 10) DIFF(I) = MOD(DIFF(I),2.0D0*PI)
            IF ((INTCRD(I,1) .GT. 20) .AND. (ABS(DIFF(I)) .GT. PI)) THEN
               IF (DIFF(I) .GT. 0.0D0) THEN
                  DIFF(I) = DIFF(I) - 2.0D0*PI
               ELSE
                  DIFF(I) = DIFF(I) + 2.0D0*PI
               END IF
            END IF
            IF (ABS(DIFF(I)) .LT. 1.0D-14) DIFF(I) = D0
            RMSINT = RMSINT + DIFF(I)*DIFF(I)
 19      CONTINUE
      END IF
C
C     The difference wich will be used for the next iteration is
C     transformed back to delocalized internals.
C
      IF (IPRINT .GE. IPRDBG) THEN
         IF (DELINT) THEN
            CALL HEADER('Difference to wanted step (primitives)',-1)
         ELSE
            CALL HEADER('Difference to wanted step',-1)
         END IF
         CALL OUTPUT(DIFF,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
      END IF
      IF (DELINT) THEN
         CALL DZERO(SCTRA,MXCOOR*MXCOOR)
         RMSINT = 0.0D0
         DO 200 I = 1, NRIC
            SCTRA(I) = DIFF(I)
 200     CONTINUE
         CALL DZERO(DIFF,MXRCRD)
         DO 205 I = 1, NNNRIC
            DO 207 J = 1, NRIC
               DIFF(I) = DIFF(I) + BMTRAN(J,I)*SCTRA(J)
 207        CONTINUE
            RMSINT = RMSINT + DIFF(I)*DIFF(I)
 205     CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Difference to wanted step',-1)
            CALL OUTPUT(DIFF,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
         END IF
      END IF
      RMSINT = SQRT(RMSINT/(1.0D0*NRIC))
      IF (RMS1ST .LT. D0) RMS1ST = RMSINT
      DIFRMS = ABS(RMSOLD-RMSINT)
      RMSOLD = RMSINT
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Root-mean-square of difference',-1)
         WRITE(LUPRI,'(A,G16.6)') 'Value:    ',RMSINT
         CALL HEADER('Change in root-mean-square of difference',-1)
         WRITE(LUPRI,'(A,G16.6)') 'Value:    ',DIFRMS
      END IF
      IF ((RMSINT .GE. RMSLIM) .AND. (RMSINT .LE. 1.0D2)
     &     .AND. (DIFRMS .GE. 1.0D-12) .AND. (ITRCRD .LE. ITRLIM)) THEN
         DO 20 I = 1, ICRTCR
            DO 22 J = 1, IINTCR
               TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*DIFF(J)
 22         CONTINUE
 20      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            STPTXT = 'Updated step #XX'
            WRITE(STPTXT(15:16),'(I2)') ITRCRD
            CALL HEADER(STPTXT,-1)
            CALL OUTPUT(TMPVEC,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
         END IF
         ITRCRD = ITRCRD + 1
         GOTO 123
C      ELSE IF (((ITRCRD .GE. ITRLIM) .AND. (RMS1ST .LE. 1.0D2*RMSINT))
C     &        .OR. ((ITRCRD .LT. ITRLIM) .AND.
C     &        (RMS1ST .LE. 1.0D2*RMSINT) .AND.
C     &        (RMSINT .GE. 1.0D2*RMSLIM))
C     &        .OR. (RMSINT .GE. 1.0D0)) THEN
      ELSE IF ((ITRCRD .GE. ITRLIM) .OR. ((RMS1ST .LE. 1.0D2*RMSINT)
     &        .AND. (RMSINT .GE. 1.1D0*RMSLIM)) .OR.
     &        (RMSINT .GE. 1.0D0)) THEN
C      ELSE IF ((ITRCRD .GE. ITRLIM) .OR. (RMSINT .GE. 1.0D2)) THEN
         IF (IPRINT .GE. IPRMIN) THEN
            WRITE(LUPRI,*)
     &           'Step does not converge, reverting to first estimate.'
            WRITE(LUPRI,*)
         END IF
         CALL DZERO(TMPVEC,MXCOOR)
         IF (.NOT. CORREC) THEN
            IF (DELINT) THEN
               DO 130 I = 1, ICRTCR
                  DO 132 J = 1, NNNRIC
                     TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*VALORG(J)
 132              CONTINUE
 130           CONTINUE
            ELSE
               DO 30 I = 1, ICRTCR
                  DO 32 J = 1, NRIC
                     TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*STPINT(J)
 32               CONTINUE
 30            CONTINUE
            END IF
         ELSE
            DO 35 I = 1, ICRTCR
               DO 37 J = 1, NRIC
                  TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*SCTRA(J)
 37            CONTINUE
 35         CONTINUE
         END IF
         DO 40 I = 1, IATOM
            DO 42 J = 1, 3
               CRDNEW(I,J+1) = CRDOLD(I,J+1) + TMPVEC((I-1)*3+J)
 42         CONTINUE
 40      CONTINUE
         CALL GETINT(IATOM,MXRCRD,CRDNEW,VALNEW)
         DO 45 I = 1, NRIC
            STPINT(I) = VALNEW(I) - VALOLD(I)
            IF (INTCRD(I,1) .GT. 10) STPINT(I) = MOD(STPINT(I),2.0D0*PI)
            IF ((INTCRD(I,1) .GT. 20)
     &           .AND. (ABS(STPINT(I)) .GT. PI)) THEN
               IF (STPINT(I) .GT. 0.0D0) THEN
                  STPINT(I) = STPINT(I) - 2.0D0*PI
               ELSE
                  STPINT(I) = STPINT(I) + 2.0D0*PI
               END IF
            END IF
 45      CONTINUE
C
C     If we use delocalized internal coordinates, we have to
C     transform the values. SCTRA is used for temporary storage.
C
         IF (DELINT) THEN
            ADJUST = .TRUE.
            CALL DZERO(SCTRA,MXRCRD)
            DO 300 I = 1, NRIC
               SCTRA(I) = STPINT(I)
 300        CONTINUE
            CALL DZERO(STPINT,MXRCRD)
            DO 305 I = 1, NNNRIC
               DO 307 J = 1, NRIC
                  STPINT(I) = STPINT(I) + BMTRAN(J,I)*SCTRA(J)
 307           CONTINUE
 305        CONTINUE
         END IF
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Adjusted step',-1)
            CALL OUTPUT(STPINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
         END IF
      END IF
      IF (CORREC) THEN
         CALL GETINT(IATOM,MXRCRD,CRDNEW,VALNEW)
         DO 50 I = 1, NRIC
            STPINT(I) = VALNEW(I) - VALORG(I)
            IF (INTCRD(I,1) .GT. 10) STPINT(I) = MOD(STPINT(I),2.0D0*PI)
            IF ((INTCRD(I,1) .GT. 20)
     &           .AND. (ABS(STPINT(I)) .GT. PI)) THEN
               IF (STPINT(I) .GT. 0.0D0) THEN
                  STPINT(I) = STPINT(I) - 2.0D0*PI
               ELSE
                  STPINT(I) = STPINT(I) + 2.0D0*PI
               END IF
            END IF
 50      CONTINUE
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Adjusted step',-1)
            CALL OUTPUT(STPINT,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
         END IF
         DO 60 I = 1, IATOM
            DO 62 J = 1, 3
               TMPVEC((I-1)*3+J) = CRDNEW(I,J+1) - CRDORG(I,J+1)
 62         CONTINUE
 60      CONTINUE
      END IF
C
C     The step in delocs is moved back to STPINT
C
      IF (DELINT .AND. (.NOT. ADJUST)) THEN
         CALL DZERO(STPINT,MXRCRD)
         DO 800 I = 1, IINTCR
            STPINT(I) = VALORG(I)
 800     CONTINUE
      END IF
C
C     We symmetrize the step, and scale it appropriately.
C
      CALL DZERO(STPSYM,MXRCRD)
C
      CALL TRACOR(CSTRA,SCTRA,1,ICRTCR,0)
      CALL DGEMM('N','N',ICRTCR,1,ICRTCR,1.D0,
     &           CSTRA,ICRTCR,
     &           TMPVEC,MXCOOR,0.D0,
     &           STPSYM,MXCOOR)
      DO 90 I = 1, ICRTCR
         STPSYM(I) = STPSYM(I)/SCLVEC(I)
 90   CONTINUE
      CALL WLKCOR(STPSYM,CSTEP,ICRTCR,MXCOOR,IPRINT)
C
C     Experimental extra stabilization by freezing atoms.
C     All Cartesian step vector components are simply set to zero
C     for the selected frozen atoms. Note that this is only applied
C     to the correctional step
C
      IF (CORREC .AND. (IFREEZ(0) .GT. 0)) THEN
         DO 510 I = 1, IFREEZ(0)
            DO 515 J = 1, 3
               CSTEP(3*(IFREEZ(I)-1)+J) = 0.D0
 515        CONTINUE
 510     CONTINUE
C     Calculate new internal coordinates
         DO 530 I = 1, IATOM
            DO 535 J = 1, 3
               CRDNEW(I,J+1) = CRDORG(I,J+1)+ CSTEP(3*(I-1)+J)
 535        CONTINUE
 530     CONTINUE
         CALL GETINT(IATOM,MXRCRD,CRDNEW,VALNEW)
         DO 550 I = 1, NRIC
            STPINT(I) = VALNEW(I) - VALORG(I)
            IF (INTCRD(I,1) .GT. 10) STPINT(I) = MOD(STPINT(I),2.0D0*PI)
            IF ((INTCRD(I,1) .GT. 20)
     &           .AND. (ABS(STPINT(I)) .GT. PI)) THEN
               IF (STPINT(I) .GT. 0.0D0) THEN
                  STPINT(I) = STPINT(I) - 2.0D0*PI
               ELSE
                  STPINT(I) = STPINT(I) + 2.0D0*PI
               END IF
            END IF
 550     CONTINUE
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Adjusted step after freezing',-1)
            CALL OUTPUT(STPINT,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
         END IF
C     If .FRZITR has been specified, the freezing is turned off after
C     the requested number of iterations
         IF ((ITRFRZ .GT. 0) .AND. ((ITRNMR+1) .GE. ITRFRZ))
     &        IFREEZ(0) = 0
      END IF
      RETURN
      END

C  /* Deck dobrki */
      SUBROUTINE DOBRKI(MXRCRD,MX2CRD,EXHER,EXSIR,EXABA,GEINFO,EVEC,
     &     BMTINV,CSTEP,TMPVEC,TMPMAT,NONT)
C
C     This procedure does the actual breaking of symmetry.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "molinp.h"
#include "gnrinf.h"
#include "optinf.h"
#include "cbiwlk.h"
#include "cbirea.h"
#include "nuclei.h"
#include "symmet.h"
#include "trkoor.h"
#include "huckel.h"
      LOGICAL EXHER, EXSIR, EXABA, BRKALS, DOCART, DOOWN, AUTOSY, NOSYM
      CHARACTER FILENM*10, TMPLN*80, TMPTXT*11, CRT*1, BSNM*80, ID3*1
      CHARACTER*1 KASYM(3,3)
      DIMENSION GEINFO(0:ITRMAX,6), NONT(MXCENT), JCO(MXAQN,MXCENT)
      DIMENSION EVEC(MX2CRD,MX2CRD), BMTINV(MXRCRD,MXCOOR)
      DIMENSION CSTEP(MXCOOR), Q(MXCENT), IQM(MXCENT)
      DIMENSION TMPVEC(IINTCR), TMPMAT(IINTCR*IINTCR)
      DIMENSION RADIUS_PCM(MXCENT), ALPHA_PCM(MXCENT)
      SAVE INDOLD
      DATA INDOLD /999/
#include "ibtfun.h"
      CALL DZERO(STPDIA,MXRCRD)
      CALL DZERO(CSTEP,MXCOOR)
      CALL DZERO(TMPVEC,IINTCR)
      CALL DZERO(TMPMAT,IINTCR*IINTCR)
      ISYMBR = MAXREP
      START = .TRUE.
      BRKALS = .FALSE.
      CALL IZERO(NUCNUM, MXCENT*8)
C
C     We find the highest symmetry number with a non-zero index.
C     Some test runs indicated this as the most effective way to
C     decrease the energy.
C
      IF ((ITRBRK .EQ. ITRNMR-1) .AND. (INDOLD .EQ. INDTOT))
     &     CALL QUIT
     &     ('*** DOBRK *** Breaking of symmetry was unsuccesful.')
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)'***** NOTE! *****'
      WRITE(LUPRI,*) 'Due to non-zero index of total ' //
     &     'Hessian, the geometry has to be altered!'
      WRITE(LUPRI,*)
C
C     The eigenvectors of the symmetry to be broken, are copied
C     to TMPMAT.
C
      NCR = 0
      IJ = 1
      DO 20 J = 1, IINTCR
         IF (EVAL(J) .LT. -1.0D-6) THEN
            DO 22 I = 1, IINTCR
               TMPMAT(IJ) = EVEC(I,J)
               IJ = IJ + 1
 22         CONTINUE
            NCR = NCR + 1
         END IF
 20   CONTINUE
      DO 27 I = 1, NCR
         STPDIA(I) = 1.0D0/SQRT(NCR*1.0D0)
 27   CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Geometry before break of symmetry',1)
         CALL PRIGEO(CORD)
         CALL HEADER('Step in diagonal representation',1)
         CALL OUTPUT(STPDIA,1,1,1,NCR,1,MXRCRD,1,LUPRI)
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Eigenvector basis',1)
         CALL OUTPUT(TMPMAT,1,IINTCR,1,NCR,IINTCR,IINTCR,1,LUPRI)
      END IF
      IJ = 1
      DO 30 II = 1, NCR
         DO 32 I = 1, IINTCR
            TMPVEC(I) = TMPVEC(I) + STPDIA(II)*TMPMAT(IJ)
            IJ = IJ + 1
 32      CONTINUE
 30   CONTINUE
C
C     The symmetry-breaking step is scaled.
C
      DO 35 I = 1, IINTCR
         TMPVEC(I) = TMPVEC(I)*1.0D0
 35   CONTINUE
C
      STPNRM = SQRT(DDOT(IINTCR,TMPVEC,1,TMPVEC,1))
      CALL DZERO(CSTEP,MXCOOR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Internal step vector',1)
         CALL OUTPUT(TMPVEC,1,1,1,IINTCR,1,IINTCR,1,LUPRI)
         WRITE(LUPRI,'(/A,F15.10/)') ' Norm of internal step:', STPNRM
      END IF
C
C     Transform internal step to Cartesian coordinates
C
      CALL DZERO(CSTEP,MXCOOR)
      DO 37 I = 1, ICRTCR
         DO 39 J = 1, IINTCR
            CSTEP(I) = CSTEP(I) + BMTINV(J,I)*TMPVEC(J)
 39      CONTINUE
 37   CONTINUE
C
C     Occasionally the symmetry breaking step is zero.
C
      TMP = 0.0D0
      DO 31 I = 1, ICRTCR
         TMP = TMP + CSTEP(I)*CSTEP(I)
 31   CONTINUE
      IF (TMP .LE. ZERGRD) CALL QUIT
     &     ('*** DOBRK *** Breaking of symmetry was unsuccesful.')
C
C     The new geometry is calculated.
C
 77   CONTINUE
      IJ = 1
      DO 45 J = 1, NUCIND
         DO 47 I = 1, 3
            CORD(I,J) = CORD(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 47      CONTINUE
 45   CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Cartesian step vector in non-symmetry basis',1)
         CALL PRIGEO(CSTEP)
         CALL HEADER('Geometry after break of symmetry',1)
         CALL PRIGEO(CORD)
      END IF
      CALL WLKMOL(CORD)
C
C     We remove the number of symmetry operations, so that new symmetry
C     will be added in the next iteration. If we have broken the symmetry
C     earlier but ends up with the same energy, we're having trouble.
C     This is brutally resolved by removing all symmetry, that is the
C     molecule has to be minimized withim the C1 point group.
C
      TMPLN = MLINE(1)
      BASIS = (TMPLN(1:5) .EQ. 'BASIS')
      I = 4
      IF (BASIS) I = 5
      TMPLN = MLINE(I)
      II = 0
      DO 49 J = 0, ITRNMR-1
         IF ((GEINFO(J,3) .LT. -0.5D0) .AND.
     &        (ABS(ENERGY-GEINFO(J,1)) .LT. THRERG*10)) II=II+1
 49   CONTINUE
      IF (II .GT. 0) BRKALS = .TRUE.
      IF (INDEX(TMPLN,'ATO') .NE. 0) THEN
         CALL LINE4(TMPLN,NONTYP,NSYMOP,CRT,KCHARG,THRS,ADDSYM,KASYM,
     &              ID3,DOCART,DOOWN)
         AUTOSY = .TRUE.
         NOSYM = .FALSE.
         IF (BRKALS) NOSYM = .TRUE.
         ID3 = ' '
         CALL LINE4W(TMPLN,NONTYP,NSYMOP,KCHARG,THRS,AUTOSY,NOSYM,
     &               KASYM,ID3,DOCART,DOOWN)
      ELSE
         READ(TMPLN(1:5),'(I5)') NONTYP 
         WRITE(TMPLN(10:20), '(A11)') '           '
         IF (BRKALS) WRITE(TMPLN(10:20), '(A11)') '0          '
      END IF
      MLINE(I) = TMPLN
C
      IATOM = 1
      I = NCLINE(IATOM)-1
      DO 50 J = 1, NONTYP
         IF (INDEX(MLINE(I),'CHA') .NE. 0) THEN
            CALL LINE5R(MLINE(I),Q(J),NONT(J),MBSI,IQM(J),JCO(1,J),
     &                  MXAQN,BASIS,ATOMBA,LMULBS,BSNM,
     &                  RADIUS_PCM(J), ALPHA_PCM(J))
         ELSE
            READ(MLINE(I),'(BN,6X,F4.1,I5)') Q, NONT(J)
         END IF
         IATOM = IATOM + NONT(J)
         I = NCLINE(IATOM)-1
 50   CONTINUE
C
C     We run over all atom types and all symmetry independent centers and
C     and expand them to all atoms. The molecule input is modified
C     according to this.
C
      IATOM = 1
      DO 60 ITYP = 1, NONTYP
         IXTRA = 0
         DO 70 I = 1, NONT(ITYP)
            MULCNT = ISTBNU(IATOM)
            DO 80 ISYMOP = 1, MAXOPR
               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  COOX=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
                  COOY=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
                  COOZ=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
                  INSLIN = NCLINE(IATOM)
                  DO 90 J = NMLINE, INSLIN, -1
                     MLINE(J+1) = MLINE(J)
 90               CONTINUE
                  DO 95 J = (IATOM+1), NUCIND
                     NCLINE(J) = NCLINE(J) + 1
 95               CONTINUE
                  TMPLN = MLINE(INSLIN + 1)
                  IPOS = INDEX(TMPLN,'Isotope=')
                  IF (IPOS .EQ. 0) THEN
                     TMPTXT = '           '
                  ELSE
                     TMPTXT = TMPLN(IPOS:IPOS+11)
                  END IF
                  IF (ABS(COOX).ge.100.0D0 .OR.
     &                ABS(COOY).ge.100.0D0 .OR.
     &                ABS(COOZ).ge.100.0D0) THEN
                     WRITE(TMPLN(5:76),'(3F20.10,1X,A11)')
     &                    COOX,COOY,COOZ,TMPTXT
                  ELSE
                     WRITE(TMPLN(5:76),'(3F20.15,1X,A11)')
     &                    COOX,COOY,COOZ,TMPTXT
                  END IF
                  MLINE(INSLIN + 1) = TMPLN
                  NMLINE = NMLINE + 1
                  IXTRA = IXTRA + 1
               END IF
 80         CONTINUE
            IATOM = IATOM + 1
 70      CONTINUE
         TMPLN = MLINE(NCLINE(IATOM - NONT(ITYP)) - 1)
         IF (INDEX(TMPLN,'CHA') .NE. 0) THEN
            CALL LINE5W(TMPLN,Q(ITYP),NONT(ITYP) + IXTRA,MBSI,BASIS,
     &           ATOMBA,LMULBS,BSNM,IQM(ITYP),JCO(1,ITYP),MXAQN,
     &           RADIUS_PCM(ITYP), ALPHA_PCM(ITYP))
          ELSE
            WRITE(TMPLN(11:15),'(I5)') NONT(ITYP) + IXTRA
         END IF
         MLINE(NCLINE(IATOM - NONT(ITYP)) - 1) = TMPLN
 60   CONTINUE
C
C     Write updated geometry to files.
C
      CALL PNCMOL(ITRNMR,IPRINT)
C
C     To mark that symmetry was broken, the index of the hessian is
C     given a negative sign. This is interpreted in PRIINF.
C
      GEINFO(ITRNMR,3) = -ABS(GEINFO(ITRNMR,3))
C
C     Several variables and arrays has to be modified/reset to be able
C     to continue calculation with new symmetry.
C
      BRKSYM = .FALSE.
      ITRBRK = ITRNMR
      INDOLD = INDTOT
      GECONV = .FALSE.
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      NEWSYM = .TRUE.
      DOHUCKEL = .TRUE.
      RSTARR = .TRUE.
      KEEPHE = .FALSE.
      ERGOLD =  ENERGY
      CALL IZERO(NUCNUM, MXCENT*8)
      CALL IZERO(NCRREP, 16)
      CALL IZERO(IPTCNT, MXCENT*48)
      CALL IZERO(NAXREP, 16)
      CALL IZERO(INDHES, 8)
      RETURN
      END

C  /* Deck prfsti */
      SUBROUTINE PRFSTI(MXRCRD,MX2CRD,NCRDHS,HESINT,EVEC,
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,VECMOD)
C
C     Controls saddle point optimization in redundant internal
C     coordinates using the partitioned rational function approach.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
      DIMENSION HESINT(MXRCRD,MXRCRD),EVEC(MX2CRD,MX2CRD)
      DIMENSION TMPMAT(MX2CRD*MX2CRD),TMPMT2(MX2CRD*MX2CRD)
      DIMENSION TMPMT3(MX2CRD*MX2CRD),TMPMT4(MX2CRD,MX2CRD)
      DIMENSION TMPMT5(MX2CRD,MX2CRD)
      DIMENSION VECMOD(MXCOOR)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     For saddle point optimizations, we can follow a specific eigenvector.
C     Due to the fact that we are separating one mode for maximization,
C     NCRDHS is temporarily reduced by one.
C
      IMODE = NSPMOD
      NCRDHS = NCRDHS-1
      IF (NSPMOD .GT. 0) THEN
         CALL FNDMOD(.FALSE.,MXRCRD,EVEC,TMPMAT,VECMOD,
     &        TMPMT2,TMPMT3,IMODE)
C
C     If the lowest mode has a gradient element of zero, we have to pick
C     another mode for maximization (or we will end up in a minimum!).
C
      ELSE
         IMODE = 1
 50      CONTINUE
         IF ((ABS(GRDDIA(IMODE)) .LT. 1.0D-10) .AND.
     &        (IMODE .LT. NCRDHS)) THEN
            IMODE = IMODE + 1
            GOTO 50
C     
C     If we find no such mode, we just set IMODE = 1, because we must
C     be at a stationary point.
C     
         ELSE IF (ABS(GRDDIA(IMODE)) .LT. 1.0D-10) THEN
            IMODE = 1
         END IF
      END IF
      IF (IPRINT .GE. IPRMAX) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Mode ',IMODE,' will be partitioned ' //
     &        'out and maximized.'
         WRITE(LUPRI,*)
      END IF
C
C     The selected mode is placed at the very end.
C
      CALL DZERO(TMPMT5,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT4,MX2CRD*MX2CRD)
      DO 400 I = 1, NCRDHS
         DO 402 J = 1, IMODE-1
            TMPMT5(I,J) = EVEC(I,J)
 402     CONTINUE
         DO 403 J = IMODE, NCRDHS-1
            TMPMT5(I,J) = EVEC(I,J+1)
 403     CONTINUE
         TMPMT4(1,I) = EVAL(I)
         TMPMT4(2,I) = GRDDIA(I)
         TMPMT5(I,NCRDHS) = EVEC(I,IMODE)
 400  CONTINUE
      TMPVAL = TMPMT4(1,IMODE)
      DO 406 I = IMODE, NCRDHS-1
         TMPMT4(1,I) = TMPMT4(1,I+1)
         TMPMT4(2,I) = TMPMT4(2,I+1)
 406  CONTINUE
      TMPMT4(1,NCRDHS) = TMPVAL
      TMPMT4(2,NCRDHS) = GRDDIA(IMODE)
C
C     We then make the augmented Hessian that will be minimized.
C
      CALL DZERO(TMPMT2,NCRDHS*NCRDHS)
      DO 410 I = 1, NCRDHS-1
         TMPMT2(I+(I-1)*NCRDHS) = TMPMT4(1,I)
         TMPMT2(I+(NCRDHS-1)*NCRDHS) = TMPMT4(2,I)
         TMPMT2(NCRDHS+(I-1)*NCRDHS) = TMPMT4(2,I)
 410  CONTINUE
      TMPMT2(NCRDHS+(NCRDHS-1)*NCRDHS) = D0
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Augmented Hessian',-1)
         CALL OUTPUT(TMPMT2,1,NCRDHS,1,NCRDHS,NCRDHS,
     &        NCRDHS,1,LUPRI)
      END IF
      CALL DZERO(TMPMT3,MX2CRD*MX2CRD)
      CALL DSITSP(NCRDHS,TMPMT2,TMPMT3)
      CALL DUNIT(TMPMT2,NCRDHS)
      CALL JACO(TMPMT3,TMPMT2,NCRDHS,NCRDHS,NCRDHS,
     &     TMPMAT(1),TMPMAT(MX2CRD+1))
      DO 420 J = 1, NCRDHS
         EVAL(J) = TMPMT3(J*(J+1)/2)
         DO 425 I = 1, NCRDHS
            EVEC(I,J) = TMPMT2(I+(J-1)*NCRDHS)
 425     CONTINUE
 420  CONTINUE
C
C     We add 1.0D5 to all eigenvalues that are essentially zero
C     for the sorting.
C
      DO 427 I = 1, NCRDHS
         IF (ABS(EVAL(I)) .LE. 1.0D-8) EVAL(I) = EVAL(I) + 1.0D5
 427  CONTINUE
      DO 430 I = 1, NCRDHS
         JMIN = I
         EMIN = EVAL(I)
         DO 435 J = (I + 1), NCRDHS
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               JMIN = J
            END IF
 435     CONTINUE
         IF (JMIN .NE. I) THEN
            CALL DSWAP(1,  EVAL  (I),1,EVAL  (JMIN),1)
            CALL DSWAP(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
C     CALL DSWAP(1,GRDDIA(I),1,GRDDIA(JMIN),1)
         END IF
 430  CONTINUE
      DO 440 I = 1, NCRDHS
         IF (ABS(ABS(EVAL(I))-1.0D5) .LT. 1.0D-3)
     &        EVAL(I) = EVAL(I) - 1.0D5
 440  CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('RF-eigenvalues',-1)
         CALL OUTPUT(EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
         CALL HEADER('RF-eigenvectors',-1)
         CALL OUTPUT(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
      END IF
      CALL PRFSTP(MX2CRD,NCRDHS,MXRCRD,EVEC,STPINT,GRDINT,
     &     TMPMAT,HESINT,1)
C
C     In the case of saddle point optimization, we also need
C     to take care of the second partition and combine the two.
C
      CMPLIM = MAX(TRSTRA*0.67D0, 0.30D0)
      DO 500 I = 1, NCRDHS-1
         IF (ABS(STPINT(I)) .GT. CMPLIM)
     &        STPINT(I) = SIGN(CMPLIM,STPINT(I))
         TMPMT4(3,I) = STPINT(I)
 500  CONTINUE
      NCRDHS = NCRDHS + 1
C
C     We then make the augmented Hessian that will be maximized.
C
      CALL DZERO(TMPMT2,NCRDHS*NCRDHS)
      TMPMT2(1) = TMPMT4(1,NCRDHS-1)
      TMPMT2(2) = TMPMT4(2,NCRDHS-1)
      TMPMT2(3) = TMPMT4(2,NCRDHS-1)
      TMPMT2(4) = D0
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Augmented Hessian',-1)
         CALL OUTPUT(TMPMT2,1,2,1,2,2,2,1,LUPRI)
      END IF
      CALL DZERO(TMPMT3,MX2CRD*MX2CRD)
      CALL DSITSP(2,TMPMT2,TMPMT3)
      CALL DUNIT(TMPMT2,2)
      CALL JACO(TMPMT3,TMPMT2,2,2,2,TMPMAT(1),TMPMAT(1+MX2CRD))
      DO 510 J = 1, 2
         EVAL(J) = TMPMT3(J*(J+1)/2)
         DO 515 I = 1, 2
            EVEC(I,J) = TMPMT2(I+(J-1)*2)
 515     CONTINUE
 510  CONTINUE
      DO 517 I = 1, 2
         IF (ABS(EVAL(I)) .LE. 1.0D-8) EVAL(I) = EVAL(I) + 1.0D5
 517  CONTINUE
C
C     The eigenvalues are sorted
C
      IF (EVAL(1) .GT. EVAL(2)) THEN
         CALL DSWAP(1,  EVAL  (1),1,EVAL  (2),1)
         CALL DSWAP(MX2CRD,EVEC(1,1),1,EVEC(1,2),1)
      END IF
      DO 520 I = 1, 2
         IF (ABS(ABS(EVAL(I))-1.0D5) .LT. 1.0D-3)
     &        EVAL(I) = EVAL(I) - 1.0D5
 520  CONTINUE
      CALL PRFSTP(MX2CRD,2,IINTCR,EVEC,STPINT,GRDINT,
     &     TMPMAT,HESINT,2)
      TMPVAL = STPINT(1)
      IF (ABS(TMPVAL) .GT. CMPLIM)
     &     TMPVAL = SIGN(CMPLIM,TMPVAL)
      CALL DZERO(STPSYM,MXCOOR)
      DO 530 I = 1, NCRDHS-2
         STPSYM(I) = TMPMT4(3,I)
 530  CONTINUE
      STPSYM(NCRDHS-1) = TMPVAL
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Diagonal RF-step',-1)
         CALL OUTPUT(STPSYM,1,1,1,NCRDHS-1,1,MXRCRD,1,LUPRI)
      END IF
C
C     The symmetry step is constructed from the original eigenvectors
C     (of the normal Hessian) and the diagonal RF-step.
C
      CALL DZERO(STPINT,MXRCRD)
      CALL DZERO(EVEC,MX2CRD*MX2CRD)
      DO 540 I = 1, NCRDHS-1
         DO 545 J = 1, NCRDHS-1
            EVEC(I,J) = TMPMT5(I,J)
 545     CONTINUE
         EVAL(I) = TMPMT4(1,I)
 540  CONTINUE
      DO 550 I = 1, NCRDHS-1
         DO 555 J = 1, NCRDHS-1
            STPINT(I) = STPINT(I) + EVEC(I,J)*STPSYM(J)
 555     CONTINUE
 550  CONTINUE
C
C     The final call to RFSTP to do scaling with respect to the
C     trust radius.
C
      CALL PRFSTP(MX2CRD,NCRDHS,MXRCRD,EVEC,STPINT,GRDINT,
     &     TMPMAT,HESINT,3)
      RETURN
      END

C  /* Deck rhoij */
      FUNCTION RHOIJ(ATMARR,I,J,ALFA,DIST,ORIG_LINDH)
C
C     Calculates the function
C        rho_ij = exp[alfa_ij (r_ref,ij^2 - r_ij^2)]
C     for the model Hessian by Roland Lindh
C
C     Revised Nov 2002 hjaaj for using covalent/metallic radii
C     instead of simplified values in DIST
C     (the problem is that they do NOT work for heavier
C     nuclei - as Roland also notes in his paper;
C     with covalent/metallic radii we should get reasonable
C     r_ref for all bonds)
C
C     for atom no. I:
C     ATMARR(I,1) = charge
C     ATMARR(I,2:4) = (x,y,z) in Bohr
C     ATMARR(I,5) = covalent/metallic radius in Angstrom
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "codata.h"
      DIMENSION ATMARR(MXCENT,8)
      DIMENSION ALFA(3,3), DIST(3,3), VEC(3)
      LOGICAL   ORIG_LINDH
      IROWI = 2
      IROWJ = 2
      IF (ATMARR(I,1) .LE. 2) IROWI = 1
      IF (ATMARR(J,1) .LE. 2) IROWJ = 1
      IF (ATMARR(I,1) .GE. 11) IROWI = 3
      IF (ATMARR(J,1) .GE. 11) IROWJ = 3
      VEC(1) = ATMARR(I,2) - ATMARR(J,2)
      VEC(2) = ATMARR(I,3) - ATMARR(J,3)
      VEC(3) = ATMARR(I,4) - ATMARR(J,4)
      V2 = VEC(1)*VEC(1) + VEC(2)*VEC(2) + VEC(3)*VEC(3)

      IF (ORIG_LINDH) THEN
         RREFIJ = DIST(IROWI,IROWJ)
      ELSE
         RREFIJ = (ATMARR(I,5) + ATMARR(J,5))/XTANG
      END IF

      RHOIJ = EXP( ALFA(IROWI,IROWJ)*(RREFIJ*RREFIJ - V2) )
c     write (lupri,*) 'i,j,rrefij,v2,rhoij',i,j,rrefij,v2,rhoij
      RETURN
      END

C  /* Deck bldhes */
      SUBROUTINE BLDHES(MXRCRD,ATMARR,HESINT)
C
C     Builds a simple model Hessian in redundant internal coordinates.
C     As described by Roland Lindh.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "optinf.h"
      DIMENSION ATMARR(MXCENT,8), HESINT(MXRCRD,MXRCRD)
C
C     Here we assign all the necessary parameters
C
      PARAMETER (STRET = 0.450D0, ROTAT = 0.150D0, TORSN = 0.005D0)
      DIMENSION ALFA(3,3), DIST(3,3)
      SAVE      ALFA, DIST
      DATA ALFA / 1.0000D0,  0.3949D0,  0.3949D0,
     &            0.3949D0,  0.2800D0,  0.2800D0,
     &            0.3949D0,  0.2800D0,  0.2800D0/
      DATA DIST / 1.35D0, 2.10D0, 2.53D0,
     &            2.10D0, 2.87D0, 3.40D0,
     &            2.53D0, 3.40D0, 3.40D0/
C
      CALL DZERO(HESINT,MXRCRD*MXRCRD)
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
C
C     We loop over all the internal coordinates and build up the
C     diagonal Hessian.
C
      NRIC = IINTCR
      IF (DELINT) NRIC = IREDIC
      DO 10 IC = 1, NRIC
        IF (INTCRD(IC,1) .LT. 10) THEN
          FAC = STRET
          RIJ = RHOIJ(ATMARR,INTCRD(IC,2),INTCRD(IC,3),ALFA,DIST,LINDHD)
          RJK = 1.0D0
          RKL = 1.0D0
        ELSE IF (INTCRD(IC,1) .LT. 20) THEN
          FAC = ROTAT
          RIJ = RHOIJ(ATMARR,INTCRD(IC,2),INTCRD(IC,3),ALFA,DIST,LINDHD)
          RJK = RHOIJ(ATMARR,INTCRD(IC,3),INTCRD(IC,4),ALFA,DIST,LINDHD)
          RKL = 1.0D0
        ELSE
          FAC = TORSN
          RIJ = RHOIJ(ATMARR,INTCRD(IC,2),INTCRD(IC,3),ALFA,DIST,LINDHD)
          RJK = RHOIJ(ATMARR,INTCRD(IC,3),INTCRD(IC,4),ALFA,DIST,LINDHD)
          RKL = RHOIJ(ATMARR,INTCRD(IC,4),INTCRD(IC,5),ALFA,DIST,LINDHD)
        END IF
C
C     Feb 2009 - VB:
C     Attempt to scale interfragment bonds, no definite improvement, thus turned off
C
C        IF ((INTCRD(IC,1) .EQ. 3) .OR.(INTCRD(IC,1) .EQ. 4)) THEN
C           FAC = FAC*0.5D0
C        END IF
        HESINT(IC,IC) = FAC*RIJ*RJK*RKL
 10   CONTINUE
      RETURN
      END

C  /* Deck nrmlvc */
      SUBROUTINE NRMLVC(VEC)
C
C     Normalizes a three-dimensional vector.
C
#include "implicit.h"
      DIMENSION VEC(3)
      VECNRM = SQRT(DDOT(3,VEC,1,VEC,1))
      VEC(1) = VEC(1)/VECNRM
      VEC(2) = VEC(2)/VECNRM
      VEC(3) = VEC(3)/VECNRM
      RETURN
      END

C  /* Deck nrmlvx */
      SUBROUTINE NRMLVX(ICRD,VEC)
C
C     Normalizes any vector.
C
#include "implicit.h"
      DIMENSION VEC(ICRD)
      VECNRM = SQRT(DDOT(ICRD,VEC,1,VEC,1))
      IF (VECNRM .LE. 1.0D-12) VECNRM = 1.0D0
      DO 10 I = 1, ICRD
         VEC(I) = VEC(I)/VECNRM
 10   CONTINUE
      RETURN
      END

C  /* Deck vecprd */
      SUBROUTINE VECPRD(VEC1,VEC2,VPRD)
C
C     Finds the vector product of two three-dimensional vectors.
C
#include "implicit.h"
      DIMENSION VEC1(3),VEC2(3),VPRD(3)
      VPRD(1) =  VEC1(2)*VEC2(3) - VEC2(2)*VEC1(3)
      VPRD(2) = -VEC1(1)*VEC2(3) + VEC2(1)*VEC1(3)
      VPRD(3) =  VEC1(1)*VEC2(2) - VEC2(1)*VEC1(2)
      RETURN
      END
C --- end of abarint.F ---
