!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C FILE    : hersym.F
C
#ifdef OLD_LOG
========================================================================
950522-kenneth
Added solvent in FIND_SYMMETRY,
as well as using NUCIND instead ATMARR in loops
950208-vebjorn
SYMADD: Completely new module to determine molecular symmetry, find the
        symmetry elements, remove symmetry-dependent atoms and finally
        return a string containing the full point group.
========================================================================
#endif
C  /* Deck symadd */
      SUBROUTINE SYMADD(WORK,LWORK,NONTYP,NSYMOP,NONT,KATOM,KASYM,
     &                  CLASS,TOLLRN_INPUT,IPRSYM_INPUT)
C*****************************************************************************
C
C     This module determines molecular symmetry automatically
C
C     vebjorn - 950204 - completely new module
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "mxsymm.h"
      DIMENSION WORK(LWORK), NONT(KATOM)
      CHARACTER KASYM(3,3)*1, CLASS*15

C optinf.h : VRML_SYM
#include "optinf.h"

      CALL QENTER('SYMADD')

C     Set control variables in "mxsymm.h"
      IF (TOLLRN_INPUT .LE. 0.0D0) THEN
         TOLLRN = 5.0D-6
      ELSE
         TOLLRN = TOLLRN_INPUT
      END IF
      TOLLRN_2 = TOLLRN**2
      ZERTOL = 1.0D-12
      MAXAXS = 100
      MAXMIR = 25
      IPRSYM = IPRSYM_INPUT

CSK   only MASTER needs to print the information
      IF( MYTID .eq. MPARID )THEN
        CALL HEADER('SYMADD: Requested addition of symmetry',-1)
        WRITE(LUPRI,'(A,1P,G9.2)') ' Symmetry test threshold: ', TOLLRN
      END IF

      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
      CALL MEMGET('REAL',KCOOR ,6*(MXCENT+1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KROT  ,6*(MAXAXS+1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KIROT ,6*(MAXAXS+1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KMIR  ,6*(MAXMIR+1),WORK,KFREE,LFREE)
      LENTMP = MAX(2*MXCENT,MXCENT+MAXAXS+MAXMIR)
      CALL MEMGET('REAL',KTEMP ,7*(LENTMP+1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP2,6*(MXCENT+1),WORK,KFREE,LFREE)
      CALL FIND_SYMMETRY(WORK(KCOOR),WORK(KROT),WORK(KIROT),WORK(KMIR),
     &            WORK(KTEMP),WORK(KTEMP2),NONTYP,NSYMOP,NONT,
     &            KATOM,KASYM,LENTMP,CLASS,VRML_SYM,WORK,KFREE,LFREE)
      CALL MEMREL('SYMADD',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SYMADD')
      RETURN
      END
C  /* Deck FIND_SYMMETRY */
      SUBROUTINE FIND_SYMMETRY(ATMARR,DRTAXS,DIRAXS,DMRPLN,
     &            TMPARR,TMP2AR,NONTYP,NSYMOP,NONT,
     &            KATOM,KASYM,LENTMP,CLASS,VRML_SYM,WORK,KFREE,LFREE)
#include "implicit.h"
#include "maxaqn.h"
#include "mxsymm.h"
#include "mxcent.h"
#include "maxorb.h"

#include "nuclei.h"
#include "r12int.h"
#include "priunit.h"
#include "cbirea.h"
#include "cbisol.h"
#include "cbihr1.h"
#include "molinp.h"
#include "numder.h"
#include "pcmlog.h"
#include "pcmdef.h"
#include "pcm.h"


! infpar.h: MASTER, MYNUM
#include "infpar.h"
C ***********************************************************
C *   ATMARR - Atomic coordinates
C *   CLASS  - Full point group
C *   DRTAXS - Rotational axes
C *   DIRAXS - Improper rotational axes
C *   DMRPLN - Mirror planes
C *   TMPARR + TMP2AR - Temporary arrays
C ***********************************************************
         DIMENSION ATMARR(6,0:MXCENT), DRTAXS(6,0:MAXAXS)
         DIMENSION DIRAXS(6,0:MAXAXS), DMRPLN(6,0:MAXMIR)
         DIMENSION TMPARR(7,0:LENTMP), WORK(*)
         DIMENSION TMP2AR(6,0:MXCENT), NONT(KATOM), JCO(MXAQN)
         LOGICAL   VRML_SYM, ADDSYM, DOCART, DOOWN, AUTOSY, NOSYM
         LOGICAL   MIDAS_SYM
         PARAMETER (D0 = 0.0D0)
         CHARACTER*1 KASYM(3,3), ID3, CRT
         CHARACTER SYMM*9, TMPLN*80, BSNM*80
         CHARACTER CLASS*15, TEXT*41, NAME*4

C        MIDAS_SYM: write out input coordinates for midas.sym analysis
         MIDAS_SYM = .TRUE.

C ::: IPRTHR - Print level threshold for verbose output. :::
         IPRTHR = 3
C
C     Just to make sure, we zero all temporary arrays
C
         CALL DZERO(ATMARR,6*(MXCENT + 1))
         CALL DZERO(DRTAXS,6*(MAXAXS + 1))
         CALL DZERO(DIRAXS,6*(MAXAXS + 1))
         CALL DZERO(DMRPLN,6*(MAXMIR + 1))
         CALL DZERO(TMPARR,7*(LENTMP + 1))
         CALL DZERO(TMP2AR,6*(MXCENT + 1))
C        IF (SOLVNT .OR. DORLM) THEN
C          N_ATOMS = NUCIND - 1
C        ELSE
           N_ATOMS = NUCIND
C        END IF
         DO 100 I = 1, N_ATOMS
            ATMARR(1,I) = CORD(1,I)
            ATMARR(2,I) = CORD(2,I)
            ATMARR(3,I) = CORD(3,I)
Chjaaj      ATMARR(4,I) = CHARGE(I)
            ATMARR(4,I) = IZATOM(I)
C           ... use true nuclear charge, not effective charge in CHARGE(I)
C               (As different atoms can have same effective charge with
C                ECP, using CHARGE for sorting is not reliable for ECP)
C               /hjaaj Mar 2004
            ATMARR(5,I) = I
            ATMARR(6,I) = ISOTOP(I)
 100     CONTINUE
         ATMARR(1,0) = N_ATOMS
         IF (IPRSYM .GE. IPRTHR) THEN
            TEXT = 'Original Coordinates'
            CALL PRN_ATMARR(TEXT, MXCENT, ATMARR)
         END IF
C
C        MIDAS Interface, write out input coordinates for midas.sym analysis
C
         IF (MIDAS_SYM) THEN
            LUCOOR = -1
            CALL GPOPEN(LUCOOR,'midasifc.cartorig','NEW',
     &                  ' ','FORMATTED',IDUMMY,.FALSE.)
C           CALL PRIGEOLU(LUCOOR,CORD)
            DO ICENT=1,N_ATOMS
               NAME=NAMEX(3*ICENT)(1:4)
               WRITE(LUCOOR,'(2X,A,3(2X,E23.16))')
     &              NAMEX(3*ICENT)(1:4),ATMARR(1,ICENT),
     &              ATMARR(2,ICENT),ATMARR(3,ICENT)
            ENDDO
            CALL GPCLOSE(LUCOOR,'KEEP')
         ENDIF
C
C
C        CALL SRTATM(MXCENT, ATMARR, .TRUE.)
         IF (.NOT. NOMOVE) THEN
C           Place molecular center of mass at the origin and
C           with principal axes of inertia along coordinate axes
            CALL CENTER_MOLECULE(ATMARR)
         END IF
CRF      We have to sort the atoms AFTER we center the molecule
C        Else errors occur when non-standard isotopes are used
         CALL SRTATM(MXCENT, ATMARR, .TRUE.)
C
C        Determine molecular point group
         CALL FIND_PGROUP(MXCENT, ATMARR, DRTAXS, DIRAXS, DMRPLN,
     &                    TMPARR, CLASS)
C ::: If requested, a VRML representation of :::
C ::: the symmetry elements is constructed.  :::
         IF (VRML_SYM) CALL MKVRSY(ATMARR,DRTAXS,MAXAXS,
     &                             DMRPLN,MAXMIR,WORK(KFREE),LFREE)
C
         IF( MYNUM .eq. MASTER )
     &       WRITE(LUPRI,'(/2A)') ' Symmetry class found: ', CLASS
         LENTMP = NINT(ATMARR(1,0) + DRTAXS(1,0) + DMRPLN(1,0))
         IF (.NOT. NOMOVE) THEN
           CALL TURN_MOLECULE(MXCENT,ATMARR,DRTAXS,DIRAXS,DMRPLN,
     &                        TMPARR,LENTMP)
         END IF
         IF (IPRSYM .GE. IPRTHR) THEN
            IF (.NOT. NOMOVE) THEN
              TEXT = 'Centered and Rotated Coordinates'
              CALL PRN_ATMARR(TEXT , MXCENT, ATMARR)
            ELSE
              WRITE (LUPRI,'(/A/)') '   System has not been centered '//
     &           ' or rotated before symmetry analysis.'
            END IF
         END IF
         IF (IPRSYM .GE. IPRTHR) THEN
            IF (NINT(DRTAXS(1,0)) .EQ. 0) THEN
               WRITE(LUPRI,'(/A)') ' No rotational axes were found.'
            ELSE
              TEXT = 'Rotational Axes'
              CALL PRNARR(TEXT, MAXAXS, DRTAXS)
            END IF
            IF (NINT(DIRAXS(1,0)) .LT. 0) THEN
               WRITE(LUPRI,'(/A)') ' Due to number of improper '
     &           // 'rotational axes, they were not determined.'
            ELSE IF (NINT(DIRAXS(1,0)) .EQ. 0) THEN
               WRITE(LUPRI,'(/A)')
     &           ' No unique improper rotational axes were found.'
            ELSE
               TEXT = 'Improper Rotational Axes'
               CALL PRNARR(TEXT, MAXAXS, DIRAXS)
            END IF
            IF (NINT(DMRPLN(1,0)) .EQ. 0) THEN
               WRITE(LUPRI,'(/A)') ' No mirror planes were found.'
            ELSE
               TEXT = 'Mirror Planes'
               CALL PRNARR(TEXT, MAXMIR, DMRPLN)
            END IF
         END IF
         CALL REDUCE(MXCENT, ATMARR,TMPARR, TMP2AR, SYMM)
         IF (IPRSYM .GE. MIN(3,IPRTHR/2)) THEN
            TEXT = 'Symmetry Independent Centres'
            CALL PRN_ATMARR(TEXT, MXCENT, ATMARR)
         END IF
C ::: The new symmetry is added. :::
         NSYMOP = 0
         DO 200 I = 1, 3
            DO 210 J = 1, 3
               KASYM(J,I) = SYMM((I-1)*3+J:(I-1)*3+J)
 210        CONTINUE
            IF (SYMM((I-1)*3+1:(I-1)*3+3) .NE. '   ') NSYMOP=NSYMOP+1
 200     CONTINUE
C ::: The input file is modified to include symmetry. :::
C ::: If there's any Angstrom mark, it's removed.    :::
         TMPLN = MLINE(NMLINE_4)
         CALL UPCASE(TMPLN)
         IF (INDEX(TMPLN,'ATO') .EQ. 0) THEN
            WRITE(MLINE(NMLINE_4)(10:10), '(I1)') NSYMOP
            WRITE(MLINE(NMLINE_4)(11:20), '(A9,A1)') SYMM, ' '
         ELSE
            CALL LINE4(MLINE(NMLINE_4),NONTYP2,NSMOP2,CRT,KCHARG,THRS,
     &                 ADDSYM,KASYM,ID3,DOCART,DOOWN)
            AUTOSY = .FALSE.
            NOSYM  = .FALSE.
            IWHERE = 1
            DO I2 = 1, 3
               DO J2 = 1, 3
                  KASYM(J2,I2) = SYMM(IWHERE:IWHERE)
                  IWHERE = IWHERE + 1
               END DO
            END DO
            ID3 = ' '
            CALL LINE4W(MLINE(NMLINE_4),NONTYP,NSYMOP,KCHARG,THRS,
     &                  AUTOSY,NOSYM,KASYM,ID3,DOCART,DOOWN)
         END IF
         IF (NSYMOP .EQ. 0) THEN
            WRITE(LUPRI, '(/A)') ' No symmetry elements were found.'
         ELSE IF (NSYMOP .EQ. 1) THEN
            WRITE(LUPRI, '(/2A)')
     &         ' The following symmetry element was found:   ',SYMM
         ELSE
            WRITE(LUPRI, '(/2A)')
     &         ' The following symmetry elements were found:   ',SYMM
         END IF
C ::: The reduced coordinates are moved back. :::
         CALL SRTOLD(MXCENT, ATMARR)
         NCTOT = NINT(ATMARR(1,0))
         IF (NCTOT .LE. 0) CALL QUIT ('No atoms defined')

         DO 250 I = 1, NCTOT
            CORD(1:3,I) = ATMARR(1:3,I)
            IOLD = NINT(ATMARR(5,I))
C           ... IOLD .ge. I after CALL SRTOLD
            RSPH(I)   = RSPH(IOLD)
            IDXSPH(I) = IDXSPH(IOLD)
            CHARGE(I) = CHARGE(IOLD)
            IZATOM(I) = IZATOM(IOLD)
C           GNUEXP(I) and MULBSI(I) have been added (WK/UniKA/04-11-2002).
            GNUEXP(I) = GNUEXP(IOLD)
            MULBSI(I) = MULBSI(IOLD)
            NOORBT(I) = NOORBT(IOLD)
            NCLINE(I) = NCLINE(IOLD)
 250     CONTINUE
C ::: Then the different elements are counted, and arrays updated. :::
         NONTYP = 1
         NONT(1:KATOM) = 0
         IELMNT = NINT(ATMARR(4,1))
         DO 275 I = 1, NCTOT
            IF (NINT(ATMARR(4,I)) .NE. IELMNT) THEN
               NONTYP = NONTYP + 1
               IELMNT = NINT(ATMARR(4,I))
            END IF
            NONT(NONTYP) = NONT(NONTYP) + 1
            NAMN(I) = NAMN(NINT(ATMARR(5,I)))
            NAMEX(3*I)   = NAMN(I)//'z'
            NAMEX(3*I-1) = NAMN(I)//'y'
            NAMEX(3*I-2) = NAMN(I)//'x'
            ISOTOP(I) = NINT(ATMARR(6,I))
 275     CONTINUE

         DO 300 I = NCTOT + 1, KATOM
            NAMN(I) = '    '
            NAMEX(3*I)   = '      '
            NAMEX(3*I-1) = '      '
            NAMEX(3*I-2) = '      '
            CORD(1,I) = D0
            CORD(2,I) = D0
            CORD(3,I) = D0
            CHARGE(I) = D0
            IZATOM(I) = 0
            ISOTOP(I) = 1
 300     CONTINUE
C ::: We update the coordinates of the input file :::
         DO 400 I = 1, NCTOT
            NCI   = NCLINE(I)
            TMPLN = MLINE(NCI)
            IF (ABS(CORD(1,I)).ge.10000.0D0 .OR.
     &          ABS(CORD(2,I)).ge.10000.0D0 .OR.
     &          ABS(CORD(3,I)).ge.10000.0D0) THEN
               WRITE(TMPLN(5:65),'(3G20.12,1X)')
     &            CORD(1,I),CORD(2,I),CORD(3,I)
            ELSE
               WRITE(TMPLN(5:65),'(3F20.14,1X)')
     &            CORD(1,I),CORD(2,I),CORD(3,I)
            END IF
C
            IF (ISOTOP(I) .GT. 1) THEN
              QMASS = DISOTP(IZATOM(I),ISOTOP(I),'MASS')
              WRITE(TMPLN(66:76),'(A8,I3)') 'Isotope=', NINT(QMASS)
            END IF
            MLINE(NCI) = TMPLN
 400     CONTINUE
C ::: Finally the number of atoms are updated. :::
         I = 1
         J = NINT(ATMARR(4,I))
         DO K = 1, NONTYP
            NCI   = NCLINE(I) - 1
            TMPLN = MLINE(NCI)
            CALL UPCASE(TMPLN)
            IF (INDEX(TMPLN,'CHA') .EQ. 0) THEN
               WRITE(MLINE(NCI)(12:15),'(I4)') NONT(K)
            ELSE
               CALL LINE5R(MLINE(NCI),Q,NONT2,MBSI,IQM,JCO,MXAQN,
     &                    BASIS,ATOMBA,LMULBS,BSNM,
     &                    RADIUS_PCM, ALPHA_PCM)
!              CALL LINE5W(MLINE(NCI),Q,NONT(K),MBSI,BASIS,ATOMBA,
!    &                    LMULBS,BSNM,IQM,JCO,MXAQN,
!    &                    RADIUS_PCM, ALPHA_PCM)
               CALL LINE5_UPD(MLINE(NCI),NONT(K))
            END IF
 510        CONTINUE
            IF (J .EQ. NINT(ATMARR(4,I))) THEN
               I = I + 1
               GOTO 510
            END IF
            J = NINT(ATMARR(4,I))
         END DO
C        IF (SOLVNT .OR. DORLM) THEN
C           NUCIND = NUCIND + 1
C           NCNTCV = NUCIND
C           NAMN(NUCIND) = 'cav '
C           NCLINE(NUCIND) = 0
C           NAMEX(3*NUCIND)     = NAMN(NUCIND)//' z'
C           NAMEX(3*NUCIND - 1) = NAMN(NUCIND)//' y'
C           NAMEX(3*NUCIND - 2) = NAMN(NUCIND)//' x'
C           CORD(1,NUCIND) = D0
C           CORD(2,NUCIND) = D0
C           CORD(3,NUCIND) = D0
C           CHARGE(NUCIND) = D0
C        END IF
      NUCIND = NCTOT
      RETURN
      END
C  /* Deck prnarr */
      SUBROUTINE PRNARR(TEXT, MXM, ARR)
C ::::::::::::::::::::::::::::::::::::::::::::::
C ::                                          ::
C ::       Prints all elements of array       ::
C ::(VB, part of hersym.F)                    ::
C ::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
         CHARACTER*(*) TEXT
         DIMENSION ARR(6,0:MXM)
         I = LEN(TEXT)
 10      CONTINUE
         IF ((I .GE. 1) .AND. (TEXT(I:I) .EQ. ' ')) THEN
            I = I - 1
            GOTO 10
         END IF
         WRITE(LUPRI,'(/1X,A/1X,80A1)') TEXT, ('-',J=1,I)
         DO 200 I=1,NINT(ARR(1,0))
            WRITE (LUPRI,11) NINT(ARR(4,I)), (ARR(J,I), J=1,3)
 11         FORMAT (I8, ' : ', 3F15.8)
 200     CONTINUE
      RETURN
      END
      SUBROUTINE PRN_ATMARR(TEXT, MXM, ATMARR)
C ::::::::::::::::::::::::::::::::::::::::::::::
C ::                                          ::
C ::       Prints all elements of array       ::
C ::(based on PRNARR, 12-feb-19 hjaaj)        ::
C ::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) TEXT
      INTEGER MXM, I, J
      REAL*8  ATMARR(6,0:MXM)
! Used from
!  nuclei.h : NAMN()
#include "mxcent.h"
#include "nuclei.h"

      I = LNBLNK(TEXT)
      WRITE(LUPRI,'(/1X,A/1X,80A1)') TEXT, ('-',J=1,I)

      DO I = 1,NINT(ATMARR(1,0))
         J = NINT(ATMARR(5,I)) ! old atom number
         WRITE (LUPRI,11) NAMN(J), NINT(ATMARR(4,I)),
     &         (ATMARR(J,I), J=1,3), NINT(ATMARR(6,I))
 11         FORMAT (A4,I4, ' : ', 3F15.8,'  Isotope',I3)
      END DO
      RETURN
      END
C  /* Deck vecang */
      FUNCTION VECANG(V1, V2)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                 ::
C ::       Finds the angle between two vectors       ::
C ::                                                 ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "pi.h"
#include "mxsymm.h"
         DIMENSION V1(3), V2(3)
         TEMP = (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))*
     &          (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3))
         IF (TEMP .GT. ZERTOL) THEN
            TEMP = (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))/SQRT(TEMP)
C ::: Parallel vectors might yield a value here that is slightly greater :::
C ::: than one. ACOS is undefined for these values, so we have to round  :::
C ::: them off, removing the "excess", accumulated numerical error.      :::
C ::: Also, ACOS is extremely sensitive around +/-1, so we round off     :::
C ::: numbers close to these values
            IF ((1.0D0-ABS(TEMP)) .LT. ZERTOL) TEMP = ANINT(TEMP)
            VECANG = ACOS(TEMP)
         ELSE
            VECANG = 0.5D0*PI
         END IF
         RETURN
      END
C  /* Deck add2ar */
      SUBROUTINE ADD2AR(MXM, ARR, IORDER, VEC)
C ::::::::::::::::::::::::::::::::::
C ::                              ::
C ::       Add Vec to array       ::
C ::                              ::
C ::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
      DIMENSION ARR(6,0:MXM), VEC(3)
      IF (ABS(VEC(1)) .LT. TOLLRN) VEC(1) = 0.0D0
      IF (ABS(VEC(2)) .LT. TOLLRN) VEC(2) = 0.0D0
      IF (ABS(VEC(3)) .LT. TOLLRN) VEC(3) = 0.0D0
      VCNORM = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
      IF (VCNORM .LT. TOLLRN) RETURN
      VEC(1) = VEC(1)/VCNORM
      VEC(2) = VEC(2)/VCNORM
      VEC(3) = VEC(3)/VCNORM
      N_ATOMS = NINT(ARR(1,0))
C
Chjaaj Sep08: remove duplicates already here, to include as many axes
Chjaaj        as possible before N_ATOMS .ge. MXM.
      DO I = 1, N_ATOMS
         J = 0
         IF (ABS(ARR(1,I) - VEC(1)) .LT. TOLLRN) J = J + 1
         IF (ABS(ARR(2,I) - VEC(2)) .LT. TOLLRN) J = J + 1
         IF (ABS(ARR(3,I) - VEC(3)) .LT. TOLLRN) J = J + 1
         IF (J .EQ. 3) THEN
Cdbg        WRITE(LUPRI,*) 'ADD2AR info   : duplicate axis',I,N+1
            RETURN
         END IF
      END DO
C
      IF (N_ATOMS .GE. MXM) THEN
Cdbg     WRITE(LUPRI,*) 'ADD2AR warning: axis not added',N_ATOMS+1,MXM
Cdbg     WRITE(LUPRI,*) 'ADD2AR warning: the axis',
Cdbg &      VEC(1),VEC(2),VEC(3)
         RETURN
      ELSE
Cdbg     WRITE(LUPRI,*) 'ADD2AR info   : axis added',N_ATOMS+1,MXM
Cdbg     WRITE(LUPRI,*) 'ADD2AR info   : the axis',
Cdbg &      VEC(1),VEC(2),VEC(3)
      END IF
C ::: The vector is placed in its correct position in the array, :::
C ::: so no sorting is necessary.                                :::
      I = 1
 100  CONTINUE
         IF ((NINT(ARR(4,I)) .GT. IORDER) .AND. (I .LE. N_ATOMS)) THEN
            I = I + 1
            GOTO 100
         END IF
      DO J = N_ATOMS, I, -1
         ARR(1,J+1) = ARR(1,J)
         ARR(2,J+1) = ARR(2,J)
         ARR(3,J+1) = ARR(3,J)
         ARR(4,J+1) = ARR(4,J)
         ARR(5,J+1) = ARR(5,J)
      END DO
      ARR(1,I) = VEC(1)
      ARR(2,I) = VEC(2)
      ARR(3,I) = VEC(3)
      ARR(4,I) = IORDER
      ARR(1,0) = N_ATOMS + 1
C
      RETURN
      END
C  /* Deck dldpax */
      SUBROUTINE DLDPAX(MXM, AXARR)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                             ::
C ::       Deletes all duplicate axes (parallel) in array        ::
C ::                                                             ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
         DIMENSION AXARR(6,0:MXM), V1(3), V2(3)
         N_ATOMS = NINT(AXARR(1,0))
         I = 1
 100     CONTINUE
            V1(1) = AXARR(1,I)
            V1(2) = AXARR(2,I)
            V1(3) = AXARR(3,I)
            J = I + 1
 150        CONTINUE
            IF (J .LE. N_ATOMS) THEN
C ::: Parallel axes yields a zero vector as cross product.  :::
C ::: Since the array is sorted, the last one will have     :::
C ::: equal or lesser order, and should be removed.         :::
               V2(1) =  V1(2)*AXARR(3,J) - AXARR(2,J)*V1(3)
               V2(2) = -V1(1)*AXARR(3,J) + AXARR(1,J)*V1(3)
               V2(3) =  V1(1)*AXARR(2,J) - AXARR(1,J)*V1(2)
               V2NORM2 = V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
               IF (V2NORM2 .LT. TOLLRN*TOLLRN) THEN
                  DO 200 I = 1, N_ATOMS - J
                     AXARR(1,J+I-1) = AXARR(1,J+I)
                     AXARR(2,J+I-1) = AXARR(2,J+I)
                     AXARR(3,J+I-1) = AXARR(3,J+I)
                     AXARR(4,J+I-1) = AXARR(4,J+I)
 200              CONTINUE
                  N_ATOMS = N_ATOMS - 1
               ELSE
                  J = J + 1
               END IF
               GOTO 150
            END IF
         I = I + 1
         IF (I .LT. N_ATOMS) GOTO 100

         AXARR(1,0) = N_ATOMS
         RETURN
      END
C  /* Deck dldpat */
      SUBROUTINE DLDPAT(MXM, ATMARR)
C :::::::::::::::::::::::::::::::::::::::::::::
C ::                                         ::
C ::       Deletes all duplicate atoms       ::
C ::                                         ::
C :::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), V1(3), V2(3)
         N_ATOMS = NINT(ATMARR(1,0))
         CALL SRTATM(MXM, ATMARR, .TRUE.)
         I = 1
 100     CONTINUE
         IF (I .LT. N_ATOMS) THEN
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            J = I + 1
 150        CONTINUE
            IF ((J .LE. N_ATOMS) .AND.
     &        (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J)))) THEN
C ::: Subtracting two position yields a zero vector for duplicate atoms. :::
               V2(1) = ATMARR(1,J) - V1(1)
               V2(2) = ATMARR(2,J) - V1(2)
               V2(3) = ATMARR(3,J) - V1(3)
               IF (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3) .LT.
     &                                           TOLLRN_2) THEN
                  DO 200 K = 1, N_ATOMS - J
                     ATMARR(1,J+K-1) = ATMARR(1,J+K)
                     ATMARR(2,J+K-1) = ATMARR(2,J+K)
                     ATMARR(3,J+K-1) = ATMARR(3,J+K)
                     ATMARR(4,J+K-1) = ATMARR(4,J+K)
                     ATMARR(5,J+K-1) = ATMARR(5,J+K)
                     ATMARR(6,J+K-1) = ATMARR(6,J+K)
 200              CONTINUE
                  N_ATOMS = N_ATOMS - 1
               ELSE
                  J = J + 1
               END IF
               GOTO 150
            END IF
            I = I + 1
            GOTO 100
         END IF
         ATMARR(1,0) = 1.0D0*N_ATOMS
         RETURN
      END
C  /* Deck srtatm */
      SUBROUTINE SRTATM(IUPPER, ARR, USENRM)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                      ::
C ::       Sorts atoms after atom number and distance from origo.         ::
C ::       Rotational axes are sorted after order, the norm does not      ::
C ::       matter for these (USENRM turns on and off the norm criteria).  ::
C ::       The method used is bubble sorting.                             ::
C ::                                                                      ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
         DIMENSION ARR(6,0:IUPPER)
         LOGICAL USENRM, SORTED
C
         N_ATOMS = NINT(ARR(1,0))
         I = 1
         SORTED = .FALSE.
 100     CONTINUE
         IF ((I .LT. N_ATOMS) .AND. (.NOT. SORTED)) THEN
            SORTED = .TRUE.
            DO 200 J = 1, N_ATOMS-I
               IF (USENRM) THEN
C              ... sort atoms
C              test 1) higher charge before lower charge
                  IF ( ARR(4,J+1) .GT. ARR(4,J) )  GOTO 300
                  IF ( ARR(4,J+1) .LT. ARR(4,J) )  GOTO 310
C              test 2) for same charge: higher isotope before lower isotope
                  IF ( ARR(6,J+1) .GT. ARR(6,J) ) GOTO 300
                  IF ( ARR(6,J+1) .LT. ARR(6,J) ) GOTO 310
                  RATM1 = ARR(1,J  )**2 + ARR(2,J  )**2 + ARR(3,J  )**2
                  RATM2 = ARR(1,J+1)**2 + ARR(2,J+1)**2 + ARR(3,J+1)**2
C              test 3) shorter distance to origo before longer distance
                  IF ( RATM1 .GT. RATM2 ) GOTO 300
                  GOTO 310 ! J and J+1 ar in correct order
C
  300             CONTINUE
C                 ... J and J+1 are NOT in correct order, switch them
                     DO K = 1, 6
                        TEMP       = ARR(K,J)
                        ARR(K,J)   = ARR(K,J+1)
                        ARR(K,J+1) = TEMP
                     END DO
                     SORTED = .FALSE.
  310             CONTINUE
               ELSE
C                 ... sort axes
                  IF (ARR(4,J+1) .GT. ARR(4,J)) THEN
C                 ... J and J+1 are NOT in correct order, switch them
                     DO K = 1, 5
                        TEMP       = ARR(K,J)
                        ARR(K,J)   = ARR(K,J+1)
                        ARR(K,J+1) = TEMP
                     END DO
                     SORTED = .FALSE.
                  END IF
               END IF
 200        CONTINUE
            I = I + 1
            GOTO 100
         END IF
         RETURN
      END
C  /* Deck srtold */
      SUBROUTINE SRTOLD(IUPPER, ARR)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                               ::
C ::       Sorts atoms after original order. The method used       ::
C ::       is bubble sorting.                                      ::
C ::                                                               ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
         DIMENSION ARR(6,0:IUPPER)
         DIMENSION ATM1(3), ATM2(3)
         LOGICAL SORTED
         N_ATOMS = NINT(ARR(1,0))
         I = 1
         SORTED = .FALSE.
 100     CONTINUE
         IF ((I .LT. N_ATOMS) .AND. (.NOT. SORTED)) THEN
            SORTED = .TRUE.
            DO 200 J = 1, N_ATOMS-I
               ATM1(1)=ARR(1,J)
               ATM1(2)=ARR(2,J)
               ATM1(3)=ARR(3,J)
               ATM2(1)=ARR(1,J+1)
               ATM2(2)=ARR(2,J+1)
               ATM2(3)=ARR(3,J+1)
               IF (NINT(ARR(5,J)) .GT. NINT(ARR(5,J+1))) THEN
                  DO 300 K = 1, 6
                     TEMP = ARR(K,J)
                     ARR(K,J) = ARR(K,J+1)
                     ARR(K,J+1) = TEMP
 300              CONTINUE
                  SORTED = .FALSE.
               END IF
 200        CONTINUE
            I = I + 1
            GOTO 100
         END IF
         RETURN
      END
C ---
      SUBROUTINE CENTER_MOLECULE(ATMARR)
C :::::::::::::::::::::::::::::::::::::::::::
C ::                                       ::
C ::       Centres molecule in space       ::
C ::                                       ::
C :::::::::::::::::::::::::::::::::::::::::::
C
C   Revision Aug 2011 HJAaJ (merged some of my old code from 2002):
C    - rotate to axes corresponding to
C      principal moment of mass inertia
C    - do it twice to minimize round-off errors
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION ATMARR(6,0:*), SUM(3)
      DIMENSION CINRTP(6),VPCMOM(3,3), TMP(3), ITMP(3)

      N_ATOMS = NINT(ATMARR(1,0))

Cdbg     write (lupri,*) 'CENTER_MOLECULE: initial coord.s etc.'
Cdbg     CALL OUTPUT(ATMARR(1,0),1,6,1,NUCLP+1,6,NUCLP+1,1,LUPRI)
Cdbg     CALL PRN_ATMARR('CENTER_MOLECULE: ATMARR start', NUCLP, ATMARR)

      IF (N_ATOMS .EQ. 1) THEN ! atom, just move to (0,0,0)
         ATMARR(1,1) = 0.0D0
         ATMARR(2,1) = 0.0D0
         ATMARR(3,1) = 0.0D0
         GO TO 100
      END IF

      ITURN = 0
  10  ITURN = ITURN + 1
         SUM(1) = 0.0D0
         SUM(2) = 0.0D0
         SUM(3) = 0.0D0
C ::: Finds the three coordinates of the center of mass,
C ::: then we calculate the inertia tensor around the center of mass,
C ::: then the molecule is moved so that the center of mass coincides with the origin
C ::: and the molecule is rotated such that the principal moment of inertia axes
C ::: are along the coordinate axes.
         TOTMAS = 0.0D0
         DO I = 1,N_ATOMS
            IZ = NINT(ATMARR(4,I))
         IF (IZ .LE. 0) CYCLE
C... hjaaj exclude point charges, MM atoms and floating orbitals
C          in the centering of the molecule.
            AMASS  = DISOTP(IZ,NINT(ATMARR(6,I)),'MASS')
            TOTMAS = TOTMAS + AMASS
            SUM(1) = SUM(1) + ATMARR(1,I)*AMASS
            SUM(2) = SUM(2) + ATMARR(2,I)*AMASS
            SUM(3) = SUM(3) + ATMARR(3,I)*AMASS
         END DO
         SUM(1) = SUM(1)/TOTMAS
         SUM(2) = SUM(2)/TOTMAS
         SUM(3) = SUM(3)/TOTMAS
C        SUM(1:3) now contains center-of-mass coordinates

Cdbg     write(lupri,*) 'Total mass ',TOTMAS
Cdbg     write(lupri,*) 'Center of mass',SUM(1:3)

         CINRTP(1:6) = 0.0D0
         DO I=1,N_ATOMS
            IZ = NINT(ATMARR(4,I))
         IF (IZ .LE. 0) CYCLE
            AMASS  = DISOTP(IZ,NINT(ATMARR(6,I)),'MASS')
            XI = ATMARR(1,I) - SUM(1)
            YI = ATMARR(2,I) - SUM(2)
            ZI = ATMARR(3,I) - SUM(3)
            CINRTP(1) = CINRTP(1) + AMASS*(YI*YI + ZI*ZI)
            CINRTP(2) = CINRTP(2) - AMASS* XI*YI
            CINRTP(3) = CINRTP(3) + AMASS*(XI*XI + ZI*ZI)
            CINRTP(4) = CINRTP(4) - AMASS* XI*ZI
            CINRTP(5) = CINRTP(5) - AMASS* YI*ZI
            CINRTP(6) = CINRTP(6) + AMASS*(XI*XI + YI*YI)
         END DO
Cdbg     write(lupri,*) 'mass inertia matrix'
Cdbg     call outpak(cinrtp,3,1,lupri)
         CALL DUNIT(VPCMOM,3)
         CALL JACO(CINRTP,VPCMOM,3,3,3,TMP,ITMP)
Cdbg     call outpak(cinrtp,3,1,lupri)
         CINRTP(2) = CINRTP(3)
         CINRTP(3) = CINRTP(6)

C        Check if molecule is mostly oblate or mostly prolate

         TMP(1:3) = CINRTP(1:3)
         CALL ORDER(DUMMY,TMP,3,0)
         IF (TMP(1)/TMP(2) .LT. 0.8D0*TMP(2)/TMP(3)) THEN ! prolate type
C        ... definitely more prolate than oblate, maybe (close to) linear,
C            sort after decreasing value to get molecular axis along z-axis
            IF (CINRTP(1) .LT. CINRTP(2)) THEN ! switch 1 and 2
               CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
               CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
               VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
            END IF
            IF (CINRTP(2) .LT. CINRTP(3)) THEN ! switch 2 and 3
               CALL DSWAP(3,VPCMOM(1,2),1,VPCMOM(1,3),1)
               CALL DSWAP(1,CINRTP(2),1,CINRTP(3),1)
               VPCMOM(1:3,2) = -VPCMOM(1:3,2)  ! to not change chirality!!!
               IF (CINRTP(1) .LT. CINRTP(2)) THEN ! switch 1 and 2
                  CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
                  CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
                  VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
               END IF
            END IF
         ELSE ! oblate type
C        ... Sorting after increasing value - will put oblate molecule
C            primarily in the x-y plane.
            IF (CINRTP(1) .GT. CINRTP(2)) THEN ! switch 1 and 2
               CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
               CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
               VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
            END IF
            IF (CINRTP(2) .GT. CINRTP(3)) THEN ! switch 2 and 3
               CALL DSWAP(3,VPCMOM(1,2),1,VPCMOM(1,3),1)
               CALL DSWAP(1,CINRTP(2),1,CINRTP(3),1)
               VPCMOM(1:3,2) = -VPCMOM(1:3,2)  ! to not change chirality!!!
               IF (CINRTP(1) .GT. CINRTP(2)) THEN ! switch 1 and 2
                  CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
                  CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
                  VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
               END IF
            END IF
         END IF

Cdbg     write (lupri,*) 'PCMOM, VPCMOM turn',ITURN
Cdbg     call output(CINRTP,1,1,1,3,1,3,1,lupri)
Cdbg     call output(VPCMOM,1,3,1,3,3,3,1,LUPRI)
C
C        Now we are ready to center and rotate molecule:
C        Note:  VPCMOM describes how to rotate coordinate axes
C               thus molecule should be rotated with VPCMOM(inverse)
C               which is same as VPC(transposed)
C
         DO I=1,N_ATOMS
            XI = ATMARR(1,I) - SUM(1)
            YI = ATMARR(2,I) - SUM(2)
            ZI = ATMARR(3,I) - SUM(3)
            ATMARR(1,I) = XI*VPCMOM(1,1)+YI*VPCMOM(2,1)+ZI*VPCMOM(3,1)
            ATMARR(2,I) = XI*VPCMOM(1,2)+YI*VPCMOM(2,2)+ZI*VPCMOM(3,2)
            ATMARR(3,I) = XI*VPCMOM(1,3)+YI*VPCMOM(2,3)+ZI*VPCMOM(3,3)
         END DO

Cdbg     write (lupri,*) 'CENTER_MOLECULE: transformed coord.s'
Cdbg     CALL OUTPUT(ATMARR(1,0),1,6,1,N_ATOMS+1,6,N_ATOMS+1,1,LUPRI)

      IF (ITURN .LE. 1) GO TO 10

 100  WRITE(LUPRI,'(/A/A)')
     &   '@   The molecule is centered at center of mass and rotated',
     &   '@   so principal axes of inertia are along coordinate axes.'
Cdbg     CALL PRN_ATMARR('CENTER_MOLECULE: coord.s now', N_ATOMS, ATMARR)
      RETURN
      END
C  /* Deck rotmol */
      SUBROUTINE ROTMOL(MXM, ATMARR, AXVEC, DEG)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                         ::
C ::       Rotate molecule in space Deg degrees around       ::
C ::           an arbitrary axis defined by AxVec            ::
C ::                                                         ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), AXVEC(3), V1(3),V2(3),V3(3),V4(3)
         AXNRM=SQRT(AXVEC(1)*AXVEC(1) + AXVEC(2)*AXVEC(2) +
     &                                  AXVEC(3)*AXVEC(3))
         AXVEC(1) = AXVEC(1)/AXNRM
         AXVEC(2) = AXVEC(2)/AXNRM
         AXVEC(3) = AXVEC(3)/AXNRM
         CSD = COS(DEG)
         SND = SIN(DEG)
         DO 100 I=1,NINT(ATMARR(1,0))
            V1(1)=ATMARR(1,I)
            V1(2)=ATMARR(2,I)
            V1(3)=ATMARR(3,I)
            DOT = AXVEC(1)*V1(1)+AXVEC(2)*V1(2)+AXVEC(3)*V1(3)
            VCLEN = ABS(DOT) - SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
C ::: Decomposes position to two vectors, one parallel(V2) and    :::
C ::: one perpendicular(V3) to the rotational axis after checking :::
C ::: wether the rotational axis passes through the atom          :::
            IF (ABS(VCLEN) .GT. ZERTOL) THEN
               V2(1) = (DOT*AXVEC(1))
               V3(1) = V1(1)-(DOT*AXVEC(1))
               V2(2) = (DOT*AXVEC(2))
               V3(2) = V1(2)-(DOT*AXVEC(2))
               V2(3) = (DOT*AXVEC(3))
               V3(3) = V1(3)-(DOT*AXVEC(3))
C ::: Finds another vector(V4) in the plane affected by rotation :::
               V4(1) =  AXVEC(2)*V3(3) - V3(2)*AXVEC(3)
               V4(2) = -AXVEC(1)*V3(3) + V3(1)*AXVEC(3)
               V4(3) =  AXVEC(1)*V3(2) - V3(1)*AXVEC(2)
               VCLEN = SQRT(V4(1)*V4(1)+V4(2)*V4(2)+V4(3)*V4(3))
               V4(1)= V4(1)/VCLEN
               V4(2)= V4(2)/VCLEN
               V4(3)= V4(3)/VCLEN
C ::: Stores the norm of the perpendicular component :::
C ::: (which is to be rotated).                      :::
               VCLEN = SQRT(V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
               V3(1)= V3(1)/VCLEN
               V3(2)= V3(2)/VCLEN
               V3(3)= V3(3)/VCLEN
C ::: V3 and V4 are now normalized and spans a plan perpendicular to the :::
C ::: rotational axis. A normalized vector is then rotated the desired   :::
C ::: amount, before it is multiplied with the stored norm, L.           :::
C ::: Finally the rotated coordinates are given back to the atom. :::
               ATMARR(1,I) = V2(1) + VCLEN*(CSD*V3(1)+SND*V4(1))
               ATMARR(2,I) = V2(2) + VCLEN*(CSD*V3(2)+SND*V4(2))
               ATMARR(3,I) = V2(3) + VCLEN*(CSD*V3(3)+SND*V4(3))
            END IF
 100     CONTINUE
         RETURN
      END
C  /* Deck trnmol */
      SUBROUTINE TURN_MOLECULE(MXM, ATMARR, DRTAXS, DIRAXS, DMRPLN,
     &                   TMPARR, LENTMP)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                    ::
C ::       Turns molecule to sensible orientation       ::
C ::                                                    ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
      DIMENSION ATMARR(6,0:MXM), DRTAXS(6,0:MAXAXS)
      DIMENSION DIRAXS(6,0:MAXAXS), DMRPLN (6,0:MAXMIR)
      DIMENSION V1(3), V2(3), V3(3), TMPARR(6,0:LENTMP)
C
      N_ROTAXS = NINT(DRTAXS(1,0))
      N_MIRROR = NINT(DMRPLN(1,0))
      N_ATOMS  = NINT(ATMARR(1,0))

      IF (N_ROTAXS .LE. 0) GOTO 9000 ! RETURN, no turning if no rotation axes (C1, Cs, or Ci symmetry)
C
      CALL DZERO(TMPARR,6*(LENTMP+1))
C ::: First we make a new array of vectors, consisting of mirror :::
C ::: planes, even-ordered rotational axes and atom positions.   :::
      ITURN = 0
c     NTOT  = N_ATOMS + N_ROTAXS + N_MIRROR
C        ... we count instead
      NTOT  = 0
C 1) rotational axes of even order (C_2, C_4, etc.)
chj Nov 2007: include both even and odd orders
         DO I = 1, N_ROTAXS
ckr
ckr         The following test prevented the molecule from being rotated
ckr         to standard orientation in case of groups in which the highest
ckr         rotation axis was of order 2.
ckr
ckr         IF ((NINT(DRTAXS(4,I)) .GT. 2) .AND.
ckr
chj         IF ((NINT(DRTAXS(4,I)) .GE. 2) .AND.
chj  &         (ABS(MOD(DRTAXS(4,I),2.0D0)) .LT. 1.0D-8)) THEN
            IF (NINT(DRTAXS(4,I)) .GE. 2) THEN
               NTOT = NTOT + 1
               DO K = 1, 4
                  TMPARR(K,NTOT) = DRTAXS(K,I)
               END DO
chj            TMPARR(4,NTOT) = 1.0D0
chj: no, we want to save the order to find highest rotation axis for z-direction
            END IF
         END DO
C 2) mirror planes
         DO I = 1, N_MIRROR
            NTOT = NTOT + 1
            DO 115 K = 1, 3
               TMPARR(K,NTOT) = DMRPLN(K,I)
 115        CONTINUE
            TMPARR(4,NTOT) = 1.0D0
         END DO
C 3) rotational axes of odd order (C_3, C_5, etc.)
chj Nov 2007: now done in "1)" above
chj      DO I = 1, N_ROTAXS
chj         IF ((NINT(DRTAXS(4,I)) .GT. 2) .AND.
chj  &         (ABS(MOD(DRTAXS(4,I),2.0D0)) .GT. 1.0D-2)) THEN
chj            NTOT = NTOT + 1
chj            DO K = 1, 4
chj               TMPARR(K,NTOT) = DRTAXS(K,I)
chj            END DO
chj         END IF
chj      END DO
C 4) and all the nuclei/point charges
         DO 130 I = 1, N_ATOMS
            NTOT = NTOT + 1
            DO 135 K = 1, 3
               TMPARR(K,NTOT) = ATMARR(K,I)
 135        CONTINUE
            TMPARR(4,NTOT) = 0.0D0
            TMPARR(6,NTOT) = 1.0D0
 130     CONTINUE
C
         TMPARR(1,0) = NTOT
C
C        Find first non-zero vector, place it along z-axis
C        (if no symmetry axes, it will become the coordinates of
C        the first of the nuclei with highest charge).
C
         I = 1
         V1(1) = TMPARR(1,1)
         V1(2) = TMPARR(2,1)
         V1(3) = TMPARR(3,1)
 150     CONTINUE
         IF (((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)) .LT. TOLLRN_2)
     &                        .AND. (I .LT. NTOT)) THEN
            I = I + 1
            V1(1) = TMPARR(1,I)
            V1(2) = TMPARR(2,I)
            V1(3) = TMPARR(3,I)
            GOTO 150
         END IF
         IF ((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)).GT.TOLLRN_2)THEN
C           V2 = target direction (along z-axis)
            V2(1) = 0.0D0
            V2(2) = 0.0D0
            V2(3) = 1.0D0
C           V3 = axis to rotate around
            V3(1) =  V1(2)
            V3(2) = -V1(1)
            V3(3) = 0.0D0
            IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                   .GT. TOLLRN_2) THEN
               DEG = VECANG(V1, V2)
               CALL ROTMOL(MXM   , ATMARR, V3, DEG)
               CALL ROTMOL(MAXAXS, DRTAXS, V3, DEG)
               CALL ROTMOL(MAXAXS, DIRAXS, V3, DEG)
               CALL ROTMOL(MAXMIR, DMRPLN, V3, DEG)
               CALL ROTMOL(NTOT  , TMPARR, V3, DEG)
            END IF
            ITURN = 1
         ELSE
            ITURN = 2
         END IF
C
C        Find second non-zero vector, rotate it so y-coordinate zero
C        (place it in xz-plane)
C        First look for a vector with zero z-coordiante
C
         JJ = I + 1
 200     CONTINUE
         IF ((ITURN .LT. 2) .AND. (I .LE. NTOT)) THEN
            IBEF = I
 250        CONTINUE
            IF (I .GT. NTOT) THEN
               I = IBEF
            ELSE IF (ABS(TMPARR(3,I)) .GT. 1.0D-8) THEN
               I = I + 1
               GOTO 250
            END IF
            V1(1) = TMPARR(1,I)
            V1(2) = TMPARR(2,I)
            V1(3) = 0.0D0
            IF ((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .GT. TOLLRN_2)
     &          .AND. (ABS(TMPARR(3,I)) .LT. TOLLRN)) THEN
C              V2 = target direction (along x-axis)
               V2(1) = 1.0D0
               V2(2) = 0.0D0
               V2(3) = 0.0D0
C              V3 = axis to rotate around
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                                        .GT. TOLLRN_2) THEN
                  DEG = VECANG(V1, V2)
                  CALL ROTMOL(MXM   , ATMARR, V3, DEG)
                  CALL ROTMOL(MAXAXS, DRTAXS, V3, DEG)
                  CALL ROTMOL(MAXAXS, DIRAXS, V3, DEG)
                  CALL ROTMOL(MAXMIR, DMRPLN, V3, DEG)
               END IF
               ITURN = 2
            END IF
            I = I + 1
            GOTO 200
         END IF
         IF (ITURN .EQ. 1) THEN
C           did not find a second non-zero vector with zero z-coordinate
C           - zero z-coordinate in temporary array and try again,
C             i.e. do not require true z-coordinate to be zero this time.
            DO 300 K = 1, NTOT
               TMPARR(3,K) = 0.0D0
 300        CONTINUE
            ITURN = -1
            I = JJ
            GOTO 200
         END IF
C
 9000 RETURN
      END
C  /* Deck linmol */
      LOGICAL FUNCTION LINMOL(MXM, ATMARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                              ::
C ::       Checks if the molecule is linear       ::
C ::                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), V1(3), V2(3), V3(3)
         N_ATOMS = NINT(ATMARR(1,0))
         LINMOL = .TRUE.
         IF (N_ATOMS .GT. 2) THEN
C ::: Builds first vector from bond between the two first atoms :::
            V1(1) = ATMARR(1,1)-ATMARR(1,2)
            V1(2) = ATMARR(2,1)-ATMARR(2,2)
            V1(3) = ATMARR(3,1)-ATMARR(3,2)
            I = 2
 100        CONTINUE
            IF((I .LT. N_ATOMS) .AND. LINMOL) THEN
C ::: If molecule is to be linear, all subsequent bonds must be parallel  :::
C ::: with the first one, and hence give a zero vector as cross product.  :::
               V2(1) = ATMARR(1,I)-ATMARR(1,I+1)
               V2(2) = ATMARR(2,I)-ATMARR(2,I+1)
               V2(3) = ATMARR(3,I)-ATMARR(3,I+1)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                   .GT. TOLLRN_2) LINMOL = .FALSE.
               I = I + 1
               GOTO 100
            END IF
         END IF
         RETURN
      END
C  /* Deck plnmol */
      LOGICAL FUNCTION PLNMOL(MXM, ATMARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                              ::
C ::       Checks if the molecule is planar       ::
C ::                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), V1(3), V2(3), V3(3)
         LOGICAL LINMOL
         N_ATOMS = NINT(ATMARR(1,0))
         PLNMOL = .TRUE.
C ::: With three atoms or less, we can be sure it's planar. If it's :::
C ::: linear, it's also planar of course.                           :::
         IF ((N_ATOMS .GT. 3) .AND. (.NOT. LINMOL(MXM, ATMARR))) THEN
C ::: We then have to find two bonds that span the possible molecular :::
C ::: plane. Since we've taken care of linear molecules, two such     :::
C ::: vectors do indeed exist. The cross product, V3, of these two    :::
C ::: vectors will be perpendicular to the plane.                     :::
            V1(1) = ATMARR(1,1) - ATMARR(1,2)
            V2(1) = ATMARR(1,2) - ATMARR(1,3)
            V1(2) = ATMARR(2,1) - ATMARR(2,2)
            V2(2) = ATMARR(2,2) - ATMARR(2,3)
            V1(3) = ATMARR(3,1) - ATMARR(3,2)
            V2(3) = ATMARR(3,2) - ATMARR(3,3)
            V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
            V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
            V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
            I = 3
 100        CONTINUE
            IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3) .LT.
     &                      TOLLRN_2) .AND. (I .LE. N_ATOMS)) THEN
               V2(1) = ATMARR(1,2) - ATMARR(1,I)
               V2(2) = ATMARR(2,2) - ATMARR(2,I)
               V2(3) = ATMARR(3,2) - ATMARR(3,I)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               I = I + 1
               GOTO 100
            END IF
C ::: We the go through the rest of the bonds, verifing that they    :::
C ::: are perpendicular to the perpendicular vector, V3, that is,    :::
C ::: that they lie in the same plane. The test is done by checking  :::
C ::: that the dot product is very close to zero.                    :::
            I = 2
 200        CONTINUE
            IF ((I .LT. N_ATOMS) .AND. PLNMOL) THEN
               V1(1) = ATMARR(1,I) - ATMARR(1,I+1)
               V1(2) = ATMARR(2,I) - ATMARR(2,I+1)
               V1(3) = ATMARR(3,I) - ATMARR(3,I+1)
               PLNMOL = (ABS(V1(1)*V3(1)+V1(2)*V3(2)+V1(3)*V3(3))
     &                                                  .LT. TOLLRN)
               I = I + 1
               GOTO 200
            END IF
         END IF
         RETURN
      END
C  /* Deck invmol */
      LOGICAL FUNCTION INVMOL(MXM, ATMARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                              ::
C ::       Checks if the molecule has a center of inversion       ::
C ::                                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
#include "mxcent.h"
         DIMENSION ATMARR(6,0:MXM), V1(3)
         DIMENSION IUSED(MXCENT)
         N_ATOMS = NINT(ATMARR(1,0))
         CALL IZERO(IUSED,MXCENT)
C ::: The routine depends on the atoms being sorted and centered.    :::
C ::: The matrix 'Used' keeps track of which atoms have been paired. :::
C ::: First, atom(s) lying in origo are dropped.                     :::
         DO 100 I=1,N_ATOMS
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            IF (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .LT. TOLLRN_2) THEN
               IUSED(I) = 1
            ELSE
               IUSED(I) = 0
            END IF
 100     CONTINUE
C ::: We then proceed with the rest of the atoms, trying to pair them. :::
         I = 1
 200     CONTINUE
         IF (I .LT. N_ATOMS) THEN
            IF (IUSED(I) .EQ. 0) THEN
               J = I + 1
 250           CONTINUE
                  IF ((IUSED(I) .EQ. 0) .AND. (J .LE. N_ATOMS)) THEN
                     IF (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J)) .AND.
     &                   NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J))) THEN
                        V1(1) = ATMARR(1,I) + ATMARR(1,J)
                        V1(2) = ATMARR(2,I) + ATMARR(2,J)
                        V1(3) = ATMARR(3,I) + ATMARR(3,J)
C ::: The sum of two anti-parallel vectors of equal length, is the  :::
C ::: zero vector. We use this to find each atoms "partner".        :::
                        IF (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .LT.
     &                       TOLLRN_2 ) THEN
                           IUSED(I) = 1
                           IUSED(J) = 1
                        END IF
                        J = J + 1
                        GOTO 250
                     END IF
                  END IF
               END IF
            I = I + 1
            GOTO 200
         END IF
C ::: Finally, we have to check that all the atoms of the molecule :::
C ::: have been inverted into another atom of the same kind.       :::
         ISUM = 0
         DO 300 I=1,N_ATOMS
            ISUM = ISUM + IUSED(I)
 300     CONTINUE
         INVMOL = (ISUM .EQ. N_ATOMS)
         RETURN
      END
C  /* Deck mirmol */
      LOGICAL FUNCTION MIRMOL(MXM,ORGARR,ATMARR,VEC,N_ATOMS_IN_PLANE)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                   ::
C ::       Checks if the mol. has mirror plane defined by vector       ::
C ::                                                                   ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
Chj Oct 2007: new parameter N_ATOMS_IN_PLANE
C             used for ROTMOL to make sure the convention is followed
C             that as many atoms as possible are placed in the x-z plane.
C             (After z-axis has been chosen as one of the axes of highest
C              order)
C
#include "implicit.h"
#include "mxcent.h"
#include "mxsymm.h"
         DIMENSION ORGARR(6,0:MXM), ATMARR(6,0:MXM), VEC(3)
         DIMENSION VROT(3), V1(3), IUSED(MXCENT)
         N_ATOMS = NINT(ORGARR(1,0))
C ::: Instead of working on the original coordinates, we make a copy.   :::
C ::: By doing so, no numerical error is built up by repeated rotation. :::
         DO 100 I = 0,N_ATOMS
            ATMARR(1,I) = ORGARR(1,I)
            ATMARR(2,I) = ORGARR(2,I)
            ATMARR(3,I) = ORGARR(3,I)
            ATMARR(4,I) = ORGARR(4,I)
            ATMARR(6,I) = ORGARR(6,I)
 100     CONTINUE
C ::: The molecule is rotated so that the possible mirror plane lies :::
C ::: in the XY-plane.                                               :::
         V1(1) = 0.0D0
         V1(2) = 0.0D0
         V1(3) = 1.0D0
         VROT(1) =  V1(2)*VEC(3) - VEC(2)*V1(3)
         VROT(2) = -V1(1)*VEC(3) + VEC(1)*V1(3)
         VROT(3) =  V1(1)*VEC(2) - VEC(1)*V1(2)
         IF (VROT(1)*VROT(1)+VROT(2)*VROT(2)+VROT(3)*VROT(3) .GT.
     &                                            ZERTOL*ZERTOL) THEN
            DEG = VECANG(V1, VEC)
            CALL ROTMOL(MXM, ATMARR, VROT, -DEG)
         END IF
C ::: All atoms lying in the XY-plane are dropped. :::
         N_ATOMS_IN_PLANE = 0
         DO 200 I=1,N_ATOMS
            IF (ATMARR(3,I)*ATMARR(3,I).LT.TOLLRN_2) THEN
               N_ATOMS_IN_PLANE = N_ATOMS_IN_PLANE + 1
               IUSED(I) = 1
            ELSE
               IUSED(I) = 0
            END IF
 200     CONTINUE
C ::: Then we go through the rest, trying to find each atoms mirror image :::
         I = 1
 300     CONTINUE
         IF (I .LT. N_ATOMS) THEN
            IF (IUSED(I) .EQ. 0) THEN
               J = I + 1
 350           CONTINUE
               IF ((IUSED(I) .EQ. 0) .AND. (J .LE. N_ATOMS) .AND.
     &             (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J))) .AND.
     &             (NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J)))) THEN
C ::: The position of the possible mirror atom is investigated.      :::
C ::: First, we subtract the XY-comp. of the first atom. If the two  :::
C ::: project into the same point in the XY-plane, Atom #2 only has  :::
C ::: its Z-comp. left. We add the Z-comp. of #1, and the sum should :::
C ::: be a zero vector if they are each others mirror image.         :::
                  V1(1) = ATMARR(1,J) - ATMARR(1,I)
                  V1(2) = ATMARR(2,J) - ATMARR(2,I)
                  V1(3) = ATMARR(3,J) + ATMARR(3,I)
                  IF (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .LT.
     &                                    10.D0*TOLLRN_2) THEN
                     IUSED(I) = 1
                     IUSED(J) = 1
                  END IF
                  J = J + 1
                  GOTO 350
               END IF
            END IF
            I = I + 1
            GOTO 300
         END IF
C ::: Finally, we have to check that all atoms have been used. :::
         ISUM = 0
         DO 400 I = 1, N_ATOMS
            ISUM = ISUM + IUSED(I)
 400     CONTINUE
         MIRMOL = (ISUM .EQ. N_ATOMS)
         RETURN
      END
C  /* Deck rotaxs */
      LOGICAL FUNCTION ROTAXS(MXM,ORGARR,ATMARR,IORDER,ROTVEC,PROPER)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                   ::
C ::       Checks if the mol. has rotational axis of order IORDER      ::
C ::                     around the given vector                       ::
C ::                                                                   ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "mxsymm.h"
         DIMENSION ORGARR(6,0:MXM), ATMARR(6,0:MXM), ROTVEC(3)
         LOGICAL PROPER,EQL
         DIMENSION VROT(3), V1(3), V2(3)
         DIMENSION IUSED(MXCENT)
         IORDR = IORDER
         IF (.NOT. PROPER) IORDR = 2*IORDER
         CON2PO = 2.0D0*PI/DBLE(IORDR)
C ::: If we are looking for improper rotational axes, we have :::
C ::: to multiply the order by a factor of two.               :::
         N_ATOMS = NINT(ORGARR(1,0))
         ATMARR(1,0) = ORGARR(1,0)
C ::: The original matrix is copied. :::
         DO 100 I=1,N_ATOMS
            ATMARR(1,I) = ORGARR(1,I)
            ATMARR(2,I) = ORGARR(2,I)
            ATMARR(3,I) = ORGARR(3,I)
            ATMARR(4,I) = ORGARR(4,I)
            ATMARR(6,I) = ORGARR(6,I)
 100     CONTINUE
C
         TEMP=SQRT(ROTVEC(1)*ROTVEC(1)+ROTVEC(2)*ROTVEC(2)+
     &             ROTVEC(3)*ROTVEC(3))
         ROTVEC(1)=ROTVEC(1)/TEMP
         ROTVEC(2)=ROTVEC(2)/TEMP
         ROTVEC(3)=ROTVEC(3)/TEMP
C ::: The molcule is rotated so that the possible rotational axis :::
C ::: lies along the Z-axis                                       :::
         V1(1) = 0.0D0
         V1(2) = 0.0D0
         V1(3) = 1.0D0
         VROT(1) =  V1(2)*ROTVEC(3) - ROTVEC(2)*V1(3)
         VROT(2) = -V1(1)*ROTVEC(3) + ROTVEC(1)*V1(3)
         VROT(3) =  V1(1)*ROTVEC(2) - ROTVEC(1)*V1(2)
         IF (VROT(1)*VROT(1)+VROT(2)*VROT(2)+VROT(3)*VROT(3) .GT.
     &                                            TOLLRN_2) THEN
            DEG = VECANG(ROTVEC, V1)
            CALL ROTMOL(MXM, ATMARR, VROT, -DEG)
         END IF
C ::: Initialize IUSED array to all atoms not used :::
         DO I=1,N_ATOMS
            IUSED(I) = 0
         END DO
C ::: Then we check if the (im)proper rotation of atom I hits another atom :::
         TOLROT = 10.0D0*TOLLRN
         I = 1
         ROTAXS = .TRUE.
 300     CONTINUE
         IF (ROTAXS .AND. (I .LE. N_ATOMS)) THEN
            IF (IUSED(I) .EQ. 0) THEN
C ::: We take the XY-component. :::
               V1(1) = ATMARR(1,I)
               V1(2) = ATMARR(2,I)
               V1(3) = 0.0D0
               V1NRM2 = V1(1)*V1(1)+V1(2)*V1(2)
               TOLLDS = MAX( 2.0D0 * TOLLRN * SQRT (V1NRM2), TOLLRN_2)
               IORD = 0
               J = I
 350           CONTINUE
               IF ( (IORD .LT. IORDR) .AND. (J .LE. N_ATOMS) ) THEN
                  EQL = (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J)) .AND.
     &                   NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J)))
               IF ( EQL ) THEN
                  V2(1) = ATMARR(1,J)
                  V2(2) = ATMARR(2,J)
                  IF (PROPER) THEN
C                    Z_I -> +Z_J ?
                     V2(3) = ATMARR(3,J) - ATMARR(3,I)
                  ELSE
C                    Z_I -> -Z_J ?
                     V2(3) = ATMARR(3,J) + ATMARR(3,I)
                  END IF
                  IF (V1NRM2 .GT. TOLLRN_2) THEN
C                    ... atom I is not on z-axis
                     DEG = VECANG(V1, V2)
                     DEG = DEG - ANINT(DEG/CON2PO)*CON2PO
                  ELSE
                     DEG = 0.0D0
                  END IF
C ::: Program checks if the two atoms got the same Z-coordinate,   :::
C ::: that the vectors have the same length, and finally that the  :::
C ::: angle between them is correct based on the rotational order. :::
                  V2NRM2 = V2(1)*V2(1)+V2(2)*V2(2)
                  IF ((ABS(V2(3)) .LT. TOLLRN) .AND.
     &               (ABS(V1NRM2-V2NRM2) .LT. TOLLDS )  ! | r^2 - (r+d)^2 | = | 2*d*r + d^2 |
     &               .AND. (ABS(DEG) .LT. TOLROT)) THEN
                     IUSED(I) = 1
                     IUSED(J) = 1
                     IF (V1NRM2 .GT. TOLLRN_2) THEN
                        IORD = IORD + 1
                     ELSE
C                       hjaaj: any order is OK for atoms on z-axis
                        IORD = IORDR
                     END IF
                  END IF
                  J = J + 1
                  GOTO 350
               END IF
               END IF
C ::: We also have to check that the symmetry is complete :::
               IF (PROPER) THEN
                  ROTAXS = (IORD .GE. IORDR)
               ELSE
                  ROTAXS = (IORD .GE. IORDR/2)
               END IF
            END IF
            I = I + 1
            GOTO 300
         END IF
C ::: Finally, we check that all atoms have been rotated into another atom. :::
         ISUM = 0
         IF (ROTAXS) THEN
            DO 400 I = 1, N_ATOMS
               ISUM = ISUM + IUSED(I)
 400        CONTINUE
         END IF
         ROTAXS = (ISUM .EQ. N_ATOMS)
         RETURN
      END
C  /* Deck findax */
      SUBROUTINE FINDAX(MXM, ATMARR, DRTAXS, TMPARR, PROPER)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                       ::
C ::       Finds all rotational axes in the molecule.      ::
C ::       PROPER determines if we are looking for         ::
C ::       proper or improper rotations.                   ::
C ::                                                       ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
#include "mxcent.h"
         DIMENSION ATMARR(6,0:MXM), DRTAXS(6,0:MAXAXS)
         DIMENSION TMPARR(7,0:MXM)
         DIMENSION V1(3), V2(3), V3(3), VSTORE(3)
         DIMENSION ICNT(2), IUSED(2*MXCENT)
         LOGICAL ROTAXS, PROPER
         DRTAXS(1,0) = 0.0D0
         N_ATOMS = NINT(ATMARR(1,0))
         IF (PROPER) THEN
            MINORD = 2
         ELSE
            MINORD = 1
         END IF
         VSTORE(1) = 0.0D0
         VSTORE(2) = 0.0D0
         VSTORE(3) = 0.0D0
         DO 95 M = 1, MXM
            IUSED(M) = 0
 95      CONTINUE
         I = 1
 100     CONTINUE
         IF (I .LE. N_ATOMS) THEN
C ::: "I" will be the lowest index in a set of equivalent atoms. :::
C ::: We need to find the rest
            J = I+1
            IORD = 1
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            V1NRM2 = V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
            TOLLDS = 2.0D0*TOLLRN*SQRT (V1NRM2)
C ::: The number of equivalent atoms are counted. :::
 110        CONTINUE
            IF (J .LE. N_ATOMS) THEN
               V2(1) = ATMARR(1,J)
               V2(2) = ATMARR(2,J)
               V2(3) = ATMARR(3,J)
               V2NRM2 = V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
               IF ((NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J))) .AND.
     &             (NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J))) .AND.
     &             (ABS( V1NRM2 - V2NRM2 ) .LT. TOLLDS )) THEN  
                 IORD = IORD + 1
                 J = J + 1
                 GOTO 110
               END IF
            END IF
C ::: Order must be two or higher to be of interest. :::
            IF (IORD .GT. 1) THEN
C ::: We "add" all equivalent atoms to yield a vector. If this vector is :::
C ::: different from the zero vector, it is the only possible rotational :::
C ::: axis. If it adds up to the zero vector, the atoms defines a        :::
C ::: linear, planar or highly symmetrical structure (it has to be       :::
C ::: symmetrical around origo).                                         :::
               V1(1) = 0.0D0
               V1(2) = 0.0D0
               V1(3) = 0.0D0
               DO 120 K=I,I+IORD-1
                  V1(1) = V1(1) + ATMARR(1,K)
                  V1(2) = V1(2) + ATMARR(2,K)
                  V1(3) = V1(3) + ATMARR(3,K)
 120           CONTINUE
               V1NRM2 = V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
               IF (V1NRM2 .LT. TOLLRN_2) THEN
C ::: All but the linear case has order higher than two. :::
                  IF (IORD .GT. 2) THEN
                     DO 133 J=I,I+IORD-1
                       V1(1) = ATMARR(1,J)
                       V1(2) = ATMARR(2,J)
                       V1(3) = ATMARR(3,J)
C ::: We let each possible pair of atoms define a plane, and examine :::
C ::: the normal vector of this plane as a possible rotational axis. :::
                       DO 130 K=J,I+IORD-1 ! only check each pair once
                         V2(1)= ATMARR(1,K)
                         V2(2)= ATMARR(2,K)
                         V2(3)= ATMARR(3,K)
                         V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
                         V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
                         V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
                         IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                                       .GT. TOLLRN_2) THEN
                           DO L=IORD,MINORD,-1
                             IF (ROTAXS(MXM,ATMARR,TMPARR,L,V3,PROPER))
     &                         THEN
                               CALL ADD2AR(MAXAXS, DRTAXS, L, V3)
                               CALL DLDPAX(MAXAXS, DRTAXS)
                             END IF
                           END DO
                         END IF
 130                   CONTINUE
C ::: We also examine all axes defined by ONE atom. :::
C     hjaaj: at most of order IORD-1, as one atom is on the axis
C            meaning we have max IORD-1 equivalent atoms
C            not on the axis.
                       IF ( (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
     &                       .GT. TOLLRN_2) THEN
                         DO L = IORD - 1, MINORD, -1
                           IF (ROTAXS(MXM,ATMARR,TMPARR,L,V1,PROPER))
     &                                    THEN
                             CALL ADD2AR(MAXAXS, DRTAXS, L, V1)
                             CALL DLDPAX(MAXAXS, DRTAXS)
                           END IF
                         END DO
                       END IF
 133                 CONTINUE
C ::: If the order is even, and the symmetry high, we might have missed    :::
C ::: axes lying between the atoms. We therefore test all possibilities.   :::
C ::: These tests are so thorough, they only have to be performed once.    :::
                     IF ((MOD(IORD,2).EQ.0).AND.(IUSED(I).EQ.0)) THEN
                        V1(1) = ATMARR(1,I)
                        V1(2) = ATMARR(2,I)
                        V1(3) = ATMARR(3,I)
C ::: We are going to add from 2 to 3 vectors. The matrix ICnt keeps :::
C ::: track of indexes, so that all combinations are tried.          :::
                        DO 160 K = 1, MIN(IORD/2-1,2)
                           ICNT(1) = 1
                           ICNT(2) = 1
C                          DO WHILE (ICNT(K) .LE. IORD)
 166                       CONTINUE
                           IF (ICNT(K) .LE. IORD) THEN
                              V2(1) = V1(1)
                              V2(2) = V1(2)
                              V2(3) = V1(3)
                              DO L = 1, K
                                 M = I-1+ICNT(L)
                                 V2(1) = V2(1) + ATMARR(1,M)
                                 V2(2) = V2(2) + ATMARR(2,M)
                                 V2(3) = V2(3) + ATMARR(3,M)
                              END DO ! L = 1,K
C ::: We check for all possible orders of rotation from Ord/2 and down :::
                              IF((V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3))
     &                                       .GT. TOLLRN_2) THEN
                               DO M = IORD/2, MINORD, -1
                               IF(ROTAXS(MXM,ATMARR,TMPARR,M,V2,PROPER))
     &                            THEN
                                  CALL ADD2AR(MAXAXS, DRTAXS, M, V2)
                                  CALL DLDPAX(MAXAXS, DRTAXS)
                               END IF
                               END DO
                              END IF
                              ICNT(1) = ICNT(1) + 1
                              IF ((K.GT.1).AND.(ICNT(1).GT.IORD)) THEN
                                 ICNT(1) = 1
                                 ICNT(2) = ICNT(2) + 1
                              END IF
                              GOTO 166
                           END IF
C                          ... END DO WHILE
 160                    CONTINUE
C ::: Molecules with high symmetry and an inversion center, might have  :::
C ::: a central axis that is difficult to get at. To find it, we add    :::
C ::: all atoms on one side of the origo. All these have positive dot   :::
C ::: products between themselves.                                      :::
                        V1(1) = 0.0D0
                        V1(2) = 0.0D0
                        V1(3) = 0.0D0
                        DO 169 ICRDAX = 1, 3
                           V1(ICRDAX) = 1.0D0
                           V3(1) = 0.0D0
                           V3(2) = 0.0D0
                           V3(3) = 0.0D0
                           DO L = I, I+IORD-1
                              V2(1) = ATMARR(1,L)
                              V2(2) = ATMARR(2,L)
                              V2(3) = ATMARR(3,L)
                              IF (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
     &                           .GT. 0.0D0) THEN
                                 V3(1) = V3(1) + V2(1)
                                 V3(2) = V3(2) + V2(2)
                                 V3(3) = V3(3) + V2(3)
                              END IF
                           END DO
                           IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &                         .GT. TOLLRN_2) THEN
                             DO M = IORD/2, MINORD, -1
                              IF (ROTAXS(MXM,ATMARR,TMPARR,M,V3,PROPER))
     &                           THEN
                                 CALL ADD2AR(MAXAXS, DRTAXS, M, V3)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                             END DO
                           END IF
                           V1(ICRDAX) = 0.0D0
 169                        CONTINUE
C ::: To be sure to get all the rotational axes in high order symmetry :::
C ::: species, all sums and differences of rotational axes are tested  :::
C ::: as well. Maximum order is the highest order among the axes.      :::
                        NAX = NINT(DRTAXS(1,0))
                        MXORDR = NINT(DRTAXS(4,1))
                        IF (NAX .GE. 2) THEN
                         K = 1
                         ICHANG = 0
 175                     CONTINUE
                         IF (K .LT. NAX) THEN
                           L = K + 1
                           V1(1) = DRTAXS(1,K)
                           V1(2) = DRTAXS(2,K)
                           V1(3) = DRTAXS(3,K)
 177                       CONTINUE
                           IF ((L .LE. NAX) .AND. (ICHANG .EQ. 0))THEN
                             V2(1) = V1(1) + DRTAXS(1,L)
                             V2(2) = V1(2) + DRTAXS(2,L)
                             V2(3) = V1(3) + DRTAXS(3,L)
                             V3(1) = V1(1) - DRTAXS(1,L)
                             V3(2) = V1(2) - DRTAXS(2,L)
                             V3(3) = V1(3) - DRTAXS(3,L)
                             IF (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
     &                                     .LT. ZERTOL) V2(1) = 1.0D0
                             IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &                                     .LT. ZERTOL) V3(1) = 1.0D0
                             DO M = MXORDR, 2, -1
                              IF (ROTAXS(MXM,ATMARR,TMPARR,M,V2,PROPER))
     &                           THEN
                                 CALL ADD2AR(MAXAXS, DRTAXS, M, V2)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                              IF (ROTAXS(MXM,ATMARR,TMPARR,M,V3,PROPER))
     &                           THEN
                                 CALL ADD2AR(MAXAXS, DRTAXS, M, V3)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                             END DO
                             IF(NAX .LT. NINT(DRTAXS(1,0))) ICHANG = 1
                             L = L + 1
                             GOTO 177
                           END IF
                           K = K + 1
                           IF (ICHANG .NE. 0) THEN
                              K = 1
                              NAX = NINT(DRTAXS(1,0))
                              MXORDR = NINT(DRTAXS(4,1))
                              ICHANG = 0
                           END IF
                           GOTO 175
                         END IF
                        END IF
C ::: Finally the atoms are marked as used. :::
                        DO 180 M = I, I+IORD-1
                           IUSED(M) = 1
 180                    CONTINUE
                     END IF
C ::: Then we move on to the linear case :::
                  ELSE
                     V1(1) = ATMARR(1,I)
                     V1(2) = ATMARR(2,I)
                     V1(3) = ATMARR(3,I)
C ::: We check for an axis of infinite (more precisely 99.) order. :::
                     IF (ROTAXS(MXM, ATMARR, TMPARR, 99, V1, PROPER))
     &                  THEN
                        CALL ADD2AR(MAXAXS, DRTAXS, 99, V1)
                        CALL DLDPAX(MAXAXS, DRTAXS)
C ::: Otherwise we look for all other possible orders along the line. :::
                     ELSE
                        DO L = 2, 8
                           IF (ROTAXS(MXM,ATMARR,TMPARR,L,V1,PROPER))
     &                        THEN
                              CALL ADD2AR(MAXAXS, DRTAXS, L, V1)
                              CALL DLDPAX(MAXAXS, DRTAXS)
                           END IF
                        END DO
C ::: Another possibility is a C2 axis perpendicular to the line. :::
C ::: We need another vector to specify the axis, so we store it. :::
C ::: If a vector already is stored, we check if the two span a   :::
C ::: plane. If they do, a C2 axis perp. to the plane is tested.  :::
                        IF (VSTORE(1)*VSTORE(1)+VSTORE(2)*VSTORE(2)+
     &                      VSTORE(3)*VSTORE(3) .LT. TOLLRN_2) THEN
                           VSTORE(1) = V1(1)
                           VSTORE(2) = V1(2)
                           VSTORE(3) = V1(3)
                        ELSE
                           V3(1) =  V1(2)*VSTORE(3) - VSTORE(2)*V1(3)
                           V3(2) = -V1(1)*VSTORE(3) + VSTORE(1)*V1(3)
                           V3(3) =  V1(1)*VSTORE(2) - VSTORE(1)*V1(2)
                           IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                                        .GT. TOLLRN_2) THEN
                              IF (ROTAXS(MXM,ATMARR,TMPARR,2,V3,PROPER))
     &                           THEN
                                 CALL ADD2AR(MAXAXS, DRTAXS, 2, V3)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                           END IF
                        END IF
                     END IF
                  END IF
C ::: The easiest case is when the atoms define the only possible axis. :::
               ELSE
                  IF (ROTAXS(MXM, ATMARR, TMPARR, IORD, V1, PROPER))
     &               THEN
                     CALL ADD2AR(MAXAXS, DRTAXS, IORD, V1)
                     CALL DLDPAX(MAXAXS, DRTAXS)
                  END IF
               END IF
            END IF
CRf   Set I to first atom not in this equivalent set and repeat                           
            I = I + IORD
            GOTO 100
         END IF
         CALL DLDPAX(MAXAXS, DRTAXS)
         CALL SRTATM(MAXAXS, DRTAXS, .FALSE.) ! sort so highest order first
         RETURN
      END
C  /* Deck fndima */
      SUBROUTINE FNDIMA(MXM, ATMARR, AXIMAR, AXARR, TMPARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                ::
C ::       Finds all improper rotational axes in the molecule       ::
C ::                                                                ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), AXIMAR(6,0:MAXAXS)
         DIMENSION AXARR(6,0:MAXAXS), TMPARR(7,0:MXM)
         DIMENSION V1(3), V2(3)
         N_ATOMS = NINT(ATMARR(1,0))
C ::: And we call upon FindAx to find all axes for this molecule. :::
         CALL FINDAX(MXM, ATMARR, AXIMAR, TMPARR, .FALSE.)
C ::: We have to multiply all orders by a factor of two to get S_n not S_2n :::
         DO 100 I = 1, NINT(AXIMAR(1,0))
            AXIMAR(4,I) = 2.0D0*AXIMAR(4,I)
 100     CONTINUE
C ::: Finally we remove all imroper rotational axes with the same order :::
C ::: and parallellity to proper axes.                                   :::
         DO 200 I = 1, NINT(AXARR(1,0))
            V1(1) = AXARR(1,I)
            V1(2) = AXARR(2,I)
            V1(3) = AXARR(3,I)
            ORD = NINT(AXARR(4,I))
            J = 1
 210        CONTINUE
            IF (J .LE. NINT(AXIMAR(1,0))) THEN
               V2(1) =  V1(2)*AXIMAR(3,J) - AXIMAR(2,J)*V1(3)
               V2(2) = -V1(1)*AXIMAR(3,J) + AXIMAR(1,J)*V1(3)
               V2(3) =  V1(1)*AXIMAR(2,J) - AXIMAR(1,J)*V1(2)
               IF ((ORD .GE. NINT(AXIMAR(4,J))) .AND. (V2(1)*V2(1)+
     &               V2(2)*V2(2)+V2(3)*V2(3) .LT. TOLLRN_2)) THEN
                  DO 205 K = 1, NINT(AXIMAR(1,0)) - J
                     AXIMAR(1,J+K-1) = AXIMAR(1,J+K)
                     AXIMAR(2,J+K-1) = AXIMAR(2,J+K)
                     AXIMAR(3,J+K-1) = AXIMAR(3,J+K)
                     AXIMAR(4,J+K-1) = AXIMAR(4,J+K)
 205              CONTINUE
                  AXIMAR(1,0) = AXIMAR(1,0) - 1.0D0
               ELSE
                  J = J + 1
               END IF
            GOTO 210
            END IF
 200     CONTINUE
         RETURN
      END
C  /* Deck fndmir */
      SUBROUTINE FNDMIR(MXM, ATMARR, DMRARR, AXARR, TMPARR)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                     ::
C ::       Finds all mirror planes in the molecule       ::
C ::                                                     ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), DMRARR(6,0:MAXMIR)
         DIMENSION AXARR(6,0:MAXAXS), TMPARR(7,0:MXM)
         DIMENSION V1(3), V2(3), V3(3), VSTORE(3), TMPAX(6)
         LOGICAL   MIRMOL, PLNMOL, LINMOL
         DMRARR(1,0) = 0.0D0
         N_ATOMS = NINT(ATMARR(1,0))
         VSTORE(1) = 0.0D0
         VSTORE(2) = 0.0D0
         VSTORE(3) = 0.0D0
C ::: First we check all planes perpendicular to one atom pos. :::
         DO 100 I = 1, N_ATOMS
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            IF ((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .GT. TOLLRN_2)
     &          .AND. (MIRMOL(MXM, ATMARR, TMPARR, V1,N2))) THEN
               CALL ADD2AR(MAXMIR, DMRARR, N2, V1)
               CALL DLDPAX(MAXMIR, DMRARR)
            END IF
 100     CONTINUE
C ::: Next thing, we go through all pairs of equal atoms, examining      :::
C ::: possible mirror planes through and between them. For the linear    :::
C ::: cases, a vector is stored. If one already is stored, we got enough :::
C ::: information to test a mirror plane.                                :::
         DO 200 I = 1, N_ATOMS-1
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            J = I + 1
            IATMNR = NINT(ATMARR(4,I))
            IATMMS = NINT(ATMARR(6,I))
 210        CONTINUE
            IF ((J .LE. N_ATOMS) .AND. (NINT(ATMARR(4,J)) .EQ. IATMNR)
     &          .AND. (NINT(ATMARR(6,J)) .EQ. IATMMS)) THEN
               V2(1) = ATMARR(1,J)
               V2(2) = ATMARR(2,J)
               V2(3) = ATMARR(3,J)
C ::: Checking the cross product, tells us whether the atoms lie on a line. :::
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)).GT. TOLLRN_2)
     &         THEN
C ::: The cross product defines the plane both atoms lie in. :::
                  IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                     CALL ADD2AR(MAXMIR, DMRARR, N2, V3)
                     CALL DLDPAX(MAXMIR, DMRARR)
                  END IF
C ::: The difference between the two vectors, define the plane between :::
C ::: the atoms.                                                       :::
                  V3(1) = V1(1) - V2(1)
                  V3(2) = V1(2) - V2(2)
                  V3(3) = V1(3) - V2(3)
                  IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                     CALL ADD2AR(MAXMIR, DMRARR, N2, V3)
                     CALL DLDPAX(MAXMIR, DMRARR)
                  END IF
               ELSE
C ::: The linear case calls for another vector to specify a plane. :::
                  IF ((VSTORE(1)*VSTORE(1)+VSTORE(2)*VSTORE(2)+
     &                 VSTORE(3)*VSTORE(3)) .LT. TOLLRN_2) THEN
                     VSTORE(1) = V1(1)
                     VSTORE(2) = V1(2)
                     VSTORE(3) = V1(3)
                  ELSE
                     V3(1) =  V1(2)*VSTORE(3) - VSTORE(2)*V1(3)
                     V3(2) = -V1(1)*VSTORE(3) + VSTORE(1)*V1(3)
                     V3(3) =  V1(1)*VSTORE(2) - VSTORE(1)*V1(2)
                     IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                   .GT. TOLLRN_2) THEN
                        IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                           CALL ADD2AR(MAXMIR, DMRARR, N2, V3)
                           CALL DLDPAX(MAXMIR, DMRARR)
                        END IF
                     END IF
                  END IF
               END IF
               J = J + 1
               GOTO 210
            END IF
 200     CONTINUE
C ::: If the molecule is planar, not linear, and have no other symmetry   :::
C ::: element than this plane, we have not found it. This must be taken   :::
C ::: care of.                                                            :::
         IF ((N_ATOMS .GT. 2) .AND. PLNMOL(MXM, ATMARR) .AND.
     &       (.NOT. LINMOL(MXM, ATMARR))) THEN
            V1(1) = ATMARR(1,1) - ATMARR(1,2)
            V2(1) = ATMARR(1,2) - ATMARR(1,3)
            V1(2) = ATMARR(2,1) - ATMARR(2,2)
            V2(2) = ATMARR(2,2) - ATMARR(2,3)
            V1(3) = ATMARR(3,1) - ATMARR(3,2)
            V2(3) = ATMARR(3,2) - ATMARR(3,3)
C ::: Take cross product of the two bonds. :::
            V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
            V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
            V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
            I = 4
 250        CONTINUE
            IF ((I .LE. N_ATOMS) .AND.
     &         (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3) .LT. TOLLRN_2)) THEN
               V2(1) = ATMARR(1,2) - ATMARR(1,I)
               V2(2) = ATMARR(2,2) - ATMARR(2,I)
               V2(3) = ATMARR(3,2) - ATMARR(3,I)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               I = I + 1
               GOTO 250
            END IF
            IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &                                     .GT. TOLLRN_2) THEN
               IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                  CALL ADD2AR(MAXMIR, DMRARR, N2, V3)
                  CALL DLDPAX(MAXMIR, DMRARR)
               END IF
            END IF
         END IF
C ::: After finding all planes, we would like to classify them, provided :::
C ::: we have found rotational axes. The axis with the highest order is  :::
C ::: used as reference, and the code is as follows:                     :::
C ::: Horizontal - 2   Vertical - 1   Other/Undecided - 0                :::
         IF (NINT(AXARR(1,0)) .GT. 0) THEN
            IFOUND = 0
            J = 1
 275        CONTINUE 
            V1(1) = AXARR(1,J)
            V1(2) = AXARR(2,J)
            V1(3) = AXARR(3,J)
            DO 300 I = 1, NINT(DMRARR(1,0))
               V2(1) = DMRARR(1,I)
               V2(2) = DMRARR(2,I)
               V2(3) = DMRARR(3,I)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
C ::: First the horizontal. :::
               IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &                                    .LT. TOLLRN_2) THEN
                  DMRARR(4,I) = 2.0D6 + DMRARR(4,I) + 0.1D0
                  IFOUND = IFOUND + 1
C ::: Then the vertical. :::
               ELSE IF (ABS(V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))
     &                                     .LT. TOLLRN) THEN
                  DMRARR(4,I) = 1.0D6 + DMRARR(4,I) + 0.1D0
                  IFOUND = IFOUND + 1
C ::: The rest must be diagonal/undecided. :::
               ELSE
                  DMRARR(4,I) =         DMRARR(4,I)
               END IF
C           (Type is recovered later by MIRTYP = 1.0D-6*DMRARR(4,I)
C            --  that is 2, 1, or 0 -- the +0.1D0 is to avoid
C            roundoff errrors when DMRARR(4,I) = 0.0D0)
 300        CONTINUE
CRF   If we have several rotational axes of the highest order, we might
CRF   not have tested the real main axis. If we have been unable to classify 
CRF   any planes thus far, we try the next one.
            IF ( (IFOUND .EQ. 0 )  .AND. (NINT(AXARR(1,0)) .GT. J ) 
     &           .AND. (NINT(AXARR(4,J)) .EQ. NINT(AXARR(4,J+1)) )) 
     &           THEN
               J = J + 1
               GOTO 275  
            END IF
            IF ( (IFOUND .GE. 1) .AND. ( J .GT. 1 ) ) THEN
CRF  We have been able to classify som mirror planes, but not using the first axis
CRF  Swich them around
               TMPAX = AXARR(:,1)
               AXARR(:,1) = AXARR(:,J)
               AXARR(:,J) = TMPAX
            END IF
            
            CALL SRTATM(MAXMIR, DMRARR, .FALSE.)
         END IF
         RETURN
      END
C  /* Deck fndgrp */
      SUBROUTINE FIND_PGROUP(MXM, ATMARR, DRTAXS, DIRAXS, DMRPLN,
     &                       TMPARR, CLASS)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                              ::
C ::       Determines which point group molecule belongs to       ::
C ::                                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
#include "priunit.h"
         DIMENSION ATMARR(6,0:MXM), DRTAXS(6,0:MAXAXS)
         DIMENSION DIRAXS(6,0:MAXAXS), DMRPLN(6,0:MAXMIR)
         DIMENSION TMPARR(7,0:MXM)
         DIMENSION V1(3), V2(3)
         CHARACTER*15 CLASS
         LOGICAL LINMOL, INVMOL
 111     FORMAT(I1)
 222     FORMAT(I2)
C
C ::: Find proper rotation axes
         CALL FINDAX(MXM, ATMARR, DRTAXS, TMPARR, .TRUE.)
C ::: If the number of rotational axes exceeds 10, the number :::
C ::: of improper rotational axes is likely to be very high   :::
C ::: and we choose to skip their determination.              :::
         IF (NINT(DRTAXS(1,0)) .LE. 10) THEN
            CALL FNDIMA(MXM, ATMARR, DIRAXS, DRTAXS, TMPARR)
         ELSE
            DIRAXS(1,0) = -1.0D0
         END IF
         CALL FNDMIR(MXM, ATMARR, DMRPLN, DRTAXS, TMPARR)
         IF (LINMOL(MXM, ATMARR)) THEN
            IF (INVMOL(MXM, ATMARR)) THEN
               CLASS = 'D(oo,h)'
            ELSE
               CLASS = 'C(oo,v)'
            END IF
         ELSE
C ::: We count rotational axes with order > 2 :::
            ITEMP = 0
            DO 100 I = 1, NINT(DRTAXS(1,0))
               IF (NINT(DRTAXS(4,I)) .GT. 2) ITEMP = ITEMP + 1
 100        CONTINUE
            IF (ITEMP .GE. 2) THEN
               IF (.NOT. INVMOL(MXM, ATMARR)) THEN
                  CLASS = 'T(d)'
               ELSE IF (NINT(DRTAXS(4,1)) .EQ. 5) THEN
                  CLASS = 'I(h)'
               ELSE
                  CLASS = 'O(h)'
               END IF
            ELSE IF (NINT(DRTAXS(1,0)) .EQ. 0) THEN
               IF (NINT(DMRPLN(1,0)) .GT. 0) THEN
                  CLASS = 'C(s)'
               ELSE IF (INVMOL(MXM, ATMARR)) THEN
                  CLASS = 'C(i)'
               ELSE
                  CLASS = 'C(1)'
               END IF
            ELSE
               V1(1) = DRTAXS(1,1)
               V1(2) = DRTAXS(2,1)
               V1(3) = DRTAXS(3,1)
C ::: Number of perpendicular C2 axes is counted. :::
               ITEMP = 0
               DO 200 I = 2, NINT(DRTAXS(1,0))
                  V2(1) = DRTAXS(1,I)
                  V2(2) = DRTAXS(2,I)
                  V2(3) = DRTAXS(3,I)
                  IF ((NINT(DRTAXS(4,I)).EQ.2) .AND. (ABS(V1(1)*V2(1)+
     &              V1(2)*V2(2)+V1(3)*V2(3)).LT.TOLLRN)) ITEMP = ITEMP+1
 200           CONTINUE
               IF (ITEMP .GE. NINT(DRTAXS(4,1))) THEN
                  MIRTYP = 1.0D-6*DMRPLN(4,1)
C                 ... any horizontal mirror plane is first in DMRPLN
                  IF ((NINT(DMRPLN(1,0)) .GT. 0) .AND.
     &                (MIRTYP .EQ. 2)) THEN
                     IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                        CLASS = 'D(nh)'
                        WRITE(CLASS(3:3),111) NINT(DRTAXS(4,1))
                     ELSE
                        CLASS = 'D(nnh)'
                        WRITE(CLASS(3:4),222) NINT(DRTAXS(4,1))
                     END IF
                  ELSE
C ::: Number of perpendicular mirror planes is counted. :::
                     ITEMP = 0
                     DO 300 I = 1, NINT(DMRPLN(1,0))
                        MIRTYP = 1.0D-6*DMRPLN(4,I)
                        IF (MIRTYP .EQ. 1) ITEMP = ITEMP + 1
 300                 CONTINUE
                     IF (ITEMP .GE. NINT(DRTAXS(4,1))) THEN
                        IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                           CLASS = 'D(nd)'
                           WRITE(CLASS(3:3),111) NINT(DRTAXS(4,1))
                        ELSE
                           CLASS = 'D(nnd)'
                           WRITE(CLASS(3:4),222) NINT(DRTAXS(4,1))
                        END IF
                     ELSE
                        IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                           CLASS = 'D(n)'
                           WRITE(CLASS(3:3),111) NINT(DRTAXS(4,1))
                        ELSE
                           CLASS = 'D(nn)'
                           WRITE(CLASS(3:4),222) NINT(DRTAXS(4,1))
                        END IF
                     END IF
                  END IF
               ELSE
                MIRTYP = 1.0D-6*DMRPLN(4,1)
                IF (NINT(DMRPLN(1,0)) .GT. 0 .AND.
     &              MIRTYP .EQ. 2) THEN
                  IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                     CLASS = 'C(nh)'
                     WRITE(CLASS(3:3),111) NINT(DRTAXS(4,1))
                  ELSE
                     CLASS = 'C(nnh)'
                     WRITE(CLASS(3:4),222) NINT(DRTAXS(4,1))
                  END IF
                ELSE
C ::: Number of perpendicular mirror planes is counted. :::
                  ITEMP = 0
                  DO 400 I = 1, NINT(DMRPLN(1,0))
                     MIRTYP = 1.0D-6*DMRPLN(4,I)
                     IF (MIRTYP .EQ. 1) ITEMP = ITEMP + 1
 400              CONTINUE
                  IF (ITEMP .GE. NINT(DRTAXS(1,0))) THEN
                     IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                        CLASS = 'C(nv)'
                        WRITE(CLASS(3:3),111) NINT(DRTAXS(4,1))
                     ELSE
                        CLASS = 'C(nnv)'
                        WRITE(CLASS(3:4),222) NINT(DRTAXS(4,1))
                     END IF
                  ELSE
                     IF (NINT(DIRAXS(1,0)) .LT. 0)
     &                  CALL FNDIMA(MXM, ATMARR, DIRAXS, DRTAXS, TMPARR)
                     IF ((NINT(DIRAXS(1,0)).GE.1) .AND.
     &                   (NINT(DIRAXS(4,1)).GE.2*NINT(DRTAXS(4,1))))THEN
                        IF (NINT(DRTAXS(4,1)) .LT. 5) THEN
                           CLASS = 'S(n)'
                           WRITE(CLASS(3:3),111) 2*NINT(DRTAXS(4,1))
                        ELSE
                           CLASS = 'S(nn)'
                           WRITE(CLASS(3:4),222) 2*NINT(DRTAXS(4,1))
                        END IF
                     ELSE
                        IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                           CLASS = 'C(n)'
                           WRITE(CLASS(3:3),111) NINT(DRTAXS(4,1))
                        ELSE
                           CLASS = 'C(nn)'
                           WRITE(CLASS(3:4),222) NINT(DRTAXS(4,1))
                        END IF
                     END IF
                  END IF
                END IF
               END IF
            END IF
         END IF
         RETURN
      END
C  /* Deck tstpar */
      SUBROUTINE TSTPAR(MXM,ARR,TARR,IATNO1,IATNO2,ICODE,VEC,GENELM)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                    ::
C ::       Test pair of atoms for given operation       ::
C ::                                                    ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
         DIMENSION ARR(7, 0:MXM), TARR(6, 0:MXM)
         DIMENSION TMPARR(6, 0:2), VEC(3)
         LOGICAL GENELM, MIRMOL, ROTAXS, INVMOL
         GENELM = .FALSE.
         TMPARR(1,0) = 2.0D0
         DO 100 I = 1,6
            TMPARR(I,1) = ARR(I, IATNO1)
            TMPARR(I,2) = ARR(I, IATNO2)
 100     CONTINUE
C ::: The icode determines what should be tested     :::
C ::: 1 - mirror plane   2 - C2 axis   3 - inversion :::
         IF( ((ICODE .EQ. 1) .AND. (MIRMOL(MXM,TMPARR,TARR,VEC,N2)))
     &   .or.((ICODE .EQ. 2) .AND. (ROTAXS(MXM,TMPARR,TARR,2,
     &                                         VEC,.TRUE.)))
     &   .or.((ICODE .EQ. 3) .AND. (INVMOL(MXM,TMPARR)))) THEN
C             GENELM = ((NINT(ARR(6,IATNO1)) .EQ. 0) .OR.
C     &                          (NINT(ARR(6,IATNO2)) .EQ. 0))
            TEMP = 0.0D0
            DO 110 I = 1, 3
               TEMP = TEMP + TMPARR(I,1) - TMPARR(I,2)
 110        CONTINUE
            IF (TEMP .GT. 0.0D0) THEN
               GENELM = (NINT(ARR(7,IATNO2)) .LT. 1)
               ARR(7, IATNO2) = ANINT(ARR(7,IATNO2)+1.0D0)
            ELSE
               GENELM = (NINT(ARR(7,IATNO1)) .LT. 1)
               ARR(7, IATNO1) = ANINT(ARR(7,IATNO1)+1.0D0)
            END IF
         END IF
         RETURN
      END
C  /* Deck reduce */
      SUBROUTINE REDUCE(MXM, ATMARR, EXTARR, TMPARR, SYMELM)
C :::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                             ::
C ::       Removes symmetry-dependent atoms      ::
C ::                                             ::
C :::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "mxsymm.h"
#include "molinp.h"
         DIMENSION ATMARR(6, 0:MXM), EXTARR(7, 0:MXM)
         DIMENSION TMPARR(6, 0:MXM), V1(3), V2(3), V3(3)
         CHARACTER STRING*9, SYMELM*9
         LOGICAL MIRMOL, ROTAXS, INVMOL, NECESS, GENELM
         N_ATOMS = NINT(ATMARR(1,0))
         EXTARR(1,0) = N_ATOMS
         DO 100 I = 1, N_ATOMS
            DO 110 J = 1, 6
               EXTARR(J,I) = ATMARR(J,I)
 110        CONTINUE
            EXTARR(7,I) = 0.0D0
 100     CONTINUE
         IELMNT = 0
         SYMELM = '         '
         STRING = 'XYZYZXZXY'
         V1(1) = 0.0D0
         V1(2) = 0.0D0
         V1(3) = 0.0D0
         DO 190 I = 1, 3
            V1(I) = 1.0D0
            IF (MIRMOL(MXM, ATMARR, TMPARR, V1, N2)) THEN
               SYMELM(IELMNT*3+1:IELMNT*3+1) = STRING(I:I)
               IELMNT = IELMNT + 1
              J = 1
 130           CONTINUE
               IF (J .LT. N_ATOMS) THEN
                  IF ((ABS(EXTARR(I,J)) .GT. TOLLRN) .AND.
     &                 (NINT(EXTARR(7,J)) .EQ. 0)) THEN
                     K = J + 1
 150                 CONTINUE
                     IF ((K .LE. N_ATOMS) .AND.
     &                (NINT(EXTARR(4,K)) .EQ. NINT(EXTARR(4,J))) .AND.
     &                (NINT(EXTARR(6,K)) .EQ. NINT(EXTARR(6,J)))) THEN
                       CALL TSTPAR(MXM,EXTARR,TMPARR,J,K,1,V1,GENELM)
                        K = K + 1
                        GOTO 150
                     END IF
                  END IF
                  J = J + 1
                  GOTO 130
               END IF
            END IF
            V1(I) = 0.0D0
 190     CONTINUE
         NROT = 0
         DO 290 I = 1, 3
            V1(I) = 1.0D0
            IF( (IELMNT.LT.3) .AND.
     &          ROTAXS(MXM,ATMARR,TMPARR,2,V1,.TRUE.) ) THEN
               SYMELM(IELMNT*3+1:IELMNT*3+2) = STRING(2+2*I:3+2*I)
               IELMNT = IELMNT + 1
               NROT = NROT + 1
               J = 1
               NECESS = .FALSE.
 250           CONTINUE
               IF (J .LT. N_ATOMS) THEN
                  V2(1) = EXTARR(1,J)
                  V2(2) = EXTARR(2,J)
                  V2(3) = EXTARR(3,J)
                  V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
                  V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
                  V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
                  IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &               .GT. TOLLRN_2) THEN
                     K = J + 1
 260                 CONTINUE
                     IF ((K .LE. N_ATOMS) .AND.
     &                (NINT(EXTARR(4,K)) .EQ. NINT(EXTARR(4,J))) .AND.
     &                (NINT(EXTARR(6,K)) .EQ. NINT(EXTARR(6,J)))) THEN
                        CALL TSTPAR(MXM,EXTARR,TMPARR,J,K,2,V1,GENELM)
                        IF (GENELM) NECESS = .TRUE.
                        K = K + 1
                        GOTO 260
                     END IF
                  END IF
                  J = J + 1
                  GOTO 250
               END IF
               IF (.NOT. NECESS) THEN
                  IELMNT = IELMNT - 1
                  SYMELM(IELMNT*3+1:IELMNT*3+2) = '  '
               END IF
C ::: Three rotations (D2) causes a problem, only two rotations are     :::
C ::: necessary. The last is removed here (extremely dirty solution!!!) :::
C ::: 12-1999 Torgeir: Expanded hack to make sure D2 symmetry is ok, and:::
C ::: make use of as much symmetry as possible.                         :::
               IF ((NROT .EQ. 3) .AND. (IELMNT .GT. 1)) THEN
                  CALL D2SYMM(ATMARR,N_ATOMS,MXM,SYMELM)
               END IF
C ::: Some molecules causes a problem, with two mirror planes (A & B) :::
C ::: and  a redundant rotation (AB), the last is removed             :::
C ::: (another dirty solution!!!)                                     :::
               IF ((IELMNT .EQ. 3) .AND. (NROT .EQ. 1)) THEN
                  IF ((SYMELM(1:1)//SYMELM(4:4)).EQ.SYMELM(7:8))THEN
                     IELMNT = IELMNT - 1
                     SYMELM(IELMNT*3+1:IELMNT*3+2) = '  '
                  END IF
               END IF
            END IF
            V1(I) = 0.0D0
 290     CONTINUE
         IF (IELMNT .EQ. 0) THEN
            IF (INVMOL(MXM, ATMARR)) THEN
               SYMELM(IELMNT*3+1:IELMNT*3+3) = 'XYZ'
               IELMNT = IELMNT + 1
               J = 1
 350           CONTINUE
               IF (J .LT. N_ATOMS) THEN
                  V2(1) = EXTARR(1,J)
                  V2(2) = EXTARR(2,J)
                  V2(3) = EXTARR(3,J)
                  IF ( (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3).GT.TOLLRN_2)
     &                .AND. (NINT(EXTARR(7,J)).EQ.0) ) THEN
                     K = J + 1
 360                 CONTINUE
                     IF ((K .LE. N_ATOMS) .AND.
     &                (NINT(EXTARR(4,K)) .EQ. NINT(EXTARR(4,J))) .AND.
     &                (NINT(EXTARR(6,K)) .EQ. NINT(EXTARR(6,J)))) THEN
                        CALL TSTPAR(MXM,EXTARR,TMPARR,J,K,3,V1,GENELM)
                        K = K + 1
                        GOTO 360
                     END IF
                  END IF
                  J = J + 1
                  GOTO 350
               END IF
            END IF
         END IF
C ::: All symmetry dependant atoms are removed. :::
         I = 1
         DO 600 J = 1, N_ATOMS
            IF (NINT(EXTARR(7,J)) .EQ. 0) THEN
               DO 610 K = 1, 6
                  ATMARR(K,I) = EXTARR(K,J)
 610           CONTINUE
               DO 615 K = 1, 3
                  IF (ABS(ATMARR(K,I)).LT.TOLLRN) ATMARR(K,I) = 0.0D0
 615           CONTINUE
               I = I + 1
            ELSE
               NC = NCLINE(NINT(EXTARR(5,J)))
               DO 620 K = NC, NMLINE - 1
                  MLINE(K) = MLINE(K+1)
 620           CONTINUE
               DO 640 K = 1, N_ATOMS
                  IF (NCLINE(NINT(EXTARR(5,K))) .GE. NC)
     &           NCLINE(NINT(EXTARR(5,K)))=NCLINE(NINT(EXTARR(5,K)))-1
 640           CONTINUE
               NMLINE = NMLINE - 1
            END IF
 600     CONTINUE
         ATMARR(1,0) = 1.0D0*(I - 1)
         RETURN
      END
C   /*Deck d2symm*/
      SUBROUTINE D2SYMM(ATMARR,NATOM,MXM,SYMELM)
C
C       ***************************************************************
C       *** This subroutine checks whether there is problems, if we ***
C       *** have D2-symmetry. These problems arise from atoms that  ***
C       *** are degenerate with respect to some rotations, and not  ***
C       *** degenerate with respect to others.                      ***
C       ***************************************************************
C
#include "implicit.h"
#include "priunit.h"
C
        PARAMETER (THRS = 1.0D-10)
        LOGICAL ONAXS, HAVAXS
        CHARACTER*(*) SYMELM
        CHARACTER*2 ROTAXIS
C
        DIMENSION ATMARR(6, 0:MXM), IATAXS(3), ONAXS(3), HAVAXS(3),
     &            ROTAXIS(3)
C
        ROTAXIS(1) = 'YZ'
        ROTAXIS(2) = 'XZ'
        ROTAXIS(3) = 'XY'
C
        DO I = 1, 3
           ONAXS (I) = .FALSE.
           HAVAXS(I) = .FALSE.
           IATAXS(I) = 0
        END DO
C
C       ********************************************************
C       *** Find out if there are pairs of atoms on the axis ***
C       ********************************************************
C
        NPRAXS = 0
        DO IATOM = 1, NATOM
           DIST = SQRT(ATMARR(1,IATOM)**2 + ATMARR(2,IATOM)**2
     &               + ATMARR(3,IATOM)**2)
           DO ICART = 1, 3
              IF ((ABS(DIST-ABS(ATMARR(ICART,IATOM))) .LT. THRS)
     &                                    .AND. (DIST .GE. THRS)) THEN
                 IATAXS(ICART) = IATAXS(ICART) + 1
                 IF (IATAXS(ICART) .EQ. 2) THEN
                    ONAXS(ICART) = .TRUE.
                    NPRAXS = NPRAXS + 1
                 END IF
              END IF
           END DO
        END DO
C
C       ***********************************************************
C       ***We figure out which of the three situations we have. ***
C       ***********************************************************
C
        IF (NPRAXS.EQ.0) THEN
C
C       ************************************************************
C       *** No pairs of atoms on any axis. We can treat it in D2 ***
C       ***             symmetry without modification.           ***
C       ************************************************************
C
        ELSE IF (NPRAXS.EQ.3) THEN
C
C       **********************************************************
C       *** Pairs of atoms on all three axis. We have to treat ***
C       *** it in C2 symmetry.                                 ***
C       **********************************************************
C
           SYMELM(2*3+1:2*3+2) = '  '
           SYMELM(  3+1:  3+2) = '  '
        ELSE
C
C       ****************************************************************
C       *** 1 axis with atom-pairs on it. We have to be careful with ***
C       *** which axes we use to generate the dependent atoms.       ***
C       ****************************************************************
C
C          *** Make sure the axes, that have atoms on it, are used. ***
C
           SYMELM(1:9) = '         '
           IEL = 0
           DO IAXS = 3, 1, -1
              IF (ONAXS(IAXS)) THEN
                 SYMELM(IEL*3+1:IEL*3+2) = ROTAXIS(IAXS)
                 HAVAXS(IAXS) = .TRUE.
                 IEL = IEL + 1
              END IF
           END DO
C
C          *** If necessary the second axis is arbitrarily chosen ***
C
           DO IAXS = 3, 1, -1
              IF ((IEL.LT.2) .AND. (.NOT.HAVAXS(IAXS))) THEN
                 SYMELM(IEL*3+1:IEL*3+2) = ROTAXIS(IAXS)
                 HAVAXS(IAXS) = .TRUE.
                 IEL = IEL + 1
              END IF
           END DO
        END IF
C
        RETURN
        END
! --- end of hersym.F ---
