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
C
C FILE    : graphic.F
C
#include <single.h>
C
C/* Deck MOLPLT */
      SUBROUTINE MOLPLT(KINDS,KOLOR,TITMOL,WORK,LWORK)
C
C Original source code provided by Kim Baldridge at SDSC.
C Punches molecule coordinate information as well as normal coordinate modes.
C Modified for use in DALTON by K.Ruud, September 1995
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
c maxorb for symmet.h
#include "maxorb.h"
#include "dummy.h"
#include "codata.h"
C
      DIMENSION KINDS(MXCENT), KOLOR(MXCENT), WORK(LWORK)
      CHARACTER*(*) TITMOL(2)
C
#include "abainf.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C
      DIMENSION KOLORS(103)
      CHARACTER*2 ATMSYM(103),SKINDS(MXCENT)
C
      DATA KOLORS/5,14,          2*11,10,1,4,2,15,14,
     *    2*11,10,6,9,7,3,14,    2*11,10*8,4*10,9,14,
     *    2*11,10*8,4*10,12,14,  2*11,14*13,10*8,5*10,14,
     *    2*11,14*13,8/
      DATA (ATMSYM(I),I = 1,103)
     1/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     2 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     3 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     4 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     5 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     6 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     7 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     8 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     9 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     O 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     1 'Md', 'No', 'Lr' /
C
C     ----- PUNCH AN INPUT FILE FOR THE -MOLPLT- PROGRAM -----
C
      LUIP = -1
      CALL GPOPEN(LUIP,'DALTON.MOL','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUIP)
      WRITE(LUIP,8000)
C
      CALL GEOANA(CORD,.FALSE.,.FALSE.,NBONDS,-1,WORK,LWORK)
C
      NKINDS    = 1
      NUCZ      = NINT(CHARGE(1))
      KINDS(1)  = NUCZ
      SKINDS(1) = ATMSYM(NUCZ)
      KOLOR(1)  = KOLORS(NUCZ)
      DO 120 IAT = 2, NUCIND
         NUCZ  = NINT(CHARGE(IAT))
         MATCH = 0
         DO 110 I = 1, NKINDS
            IF(NUCZ .EQ. KINDS(I)) MATCH = MATCH + 1
  110    CONTINUE
         IF(MATCH .EQ. 0) THEN
            NKINDS = NKINDS+1
            KINDS(NKINDS)  = NUCZ
            SKINDS(NKINDS) = ATMSYM(NUCZ)
            KOLOR(NKINDS)  = KOLORS(NUCZ)
         END IF
  120 CONTINUE
C
C     ----- OPTIONS CARD -----
C
      WRITE(LUIP,8010) NATOMS,NKINDS,NBONDS,(TITMOL(J),J=1,2)
C
C     ----- PUNCH ATOMIC SYMBOL, KOLOR, SIZE -----
C     THE BALL SIZES ARE DETERMINED BY PLAYING WITH THE RADIAL
C     EXPECTATION VALUES OF THE HIGHEST AO IN C,SI,GE,SN,PB.
C
      CARBON = 0.3D+00
      DO 210 I=1, NKINDS
         NUCZ = KINDS(I)
                        SIZE = 0.75D+00* CARBON
         IF(NUCZ.GT. 2) SIZE =           CARBON
         IF(NUCZ.GT.10) SIZE = 1.6D+00 * CARBON
         IF(NUCZ.GT.18) SIZE = 1.7D+00 * CARBON
         IF(NUCZ.GT.36) SIZE = 1.9D+00 * CARBON
         IF(NUCZ.GT.54) SIZE = 2.0D+00 * CARBON
         IF(NUCZ.GT.86) SIZE = 2.1D+00 * CARBON
         WRITE(LUIP,8020) SKINDS(I),KOLOR(I),SIZE
  210 CONTINUE
C
C     ----- PUNCH ATOMIC COORDINATES IN ANGSTROMS -----
C
      IAT = 0
      DO 310 INUC=1,NUCIND
         DO 311 IA = 0, MAXOPR
            IF (IBTAND(IA,ISTBNU(INUC)) .EQ. 0) THEN
               IAT = IAT + 1
               NUCZ = NINT(CHARGE(INUC))
               X = XTANG * PT(IBTAND(ISYMAX(1,1),IA)) * CORD(1,INUC)
               Y = XTANG * PT(IBTAND(ISYMAX(2,1),IA)) * CORD(2,INUC)
               Z = XTANG * PT(IBTAND(ISYMAX(3,1),IA)) * CORD(3,INUC)
               WRITE(LUIP,8030) ATMSYM(NUCZ),X,Y,Z
            END IF
 311     CONTINUE 
 310  CONTINUE
C
C     ----- PUNCH BONDED ATOM LIST -----
C
      CALL GEOANA(CORD,.FALSE.,.FALSE.,NBONDS,LUIP,WORK,LWORK)
      WRITE (LUIP,8040)
      CALL GPCLOSE(LUIP,'KEEP')
      RETURN
C
 8000 FORMAT('-------- START OF -MOLPLT- INPUT FILE ----------')
 8010 FORMAT('NATOMS=',I4,'   NKINDS=',I4,'    NBONDS=',I4
     &        /T6,A/T6,A)
 8020 FORMAT(A4,I2,F5.2)
 8030 FORMAT(A4,3F12.6)
 8040 FORMAT('-------- END OF -MOLPLT- INPUT FILE ----------')
      END
C/* Deck PLTORB */
      SUBROUTINE PLTORB(WORK,LWORK)
C
C     Original source code provided by Kim Baldridge at SDSC.
C     Punches molecule orbital information
C     Modified for use in DALTON by K.Ruud, September 1995
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "codata.h"
C
      PARAMETER (ONE=1.0D+00, HALF=0.5D+00, PT75 = 0.75D+00, 
     *           PT187 = 1.875D+00,
     *           PT6562 = 6.5625D+00)
C
      CHARACTER*1 TYPE
      CHARACTER*2 ATMSYM(103)
C
      DIMENSION WORK(LWORK)
C
#include "shells.h"
#include "nuclei.h"
#include "symmet.h"
#include "primit.h"
#include "ibtfun.h"
C
      DATA (ATMSYM(I),I = 1,103)
     1/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     2 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     3 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     4 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     5 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     6 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     7 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     8 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     9 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     O 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     1 'Md', 'No', 'Lr' /
C
C     ----- PUNCH AN INPUT FILE FOR THE -PLTORB- PROGRAM -----
C
      LUIP = -1
      CALL GPOPEN(LUIP,'DALTON.ORB','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND (LUIP)
      WRITE(LUIP,8000)
C
      PI32 = PI * SQRTPI
C
C     ----- PUNCH THE OPTION CARDS -----
C
      CALL GEOANA(CORD,.FALSE.,.FALSE.,NBONDS,-1,WORK,LWORK)
C
      NUM = 0
      DO 10 I = 0, MAXREP
         NUM = NUM + NAOS(I + 1)
 10   CONTINUE 
C
      WRITE(LUIP,8010) NATOMS, NBONDS, NUM
C
C     ----- PUNCH ATOMIC COORDINATES IN ANGSTROMS -----
C
      IAT = 0
      DO 310 INUC=1,NUCIND
         DO 311 IA = 0, MAXOPR
            IF (IBTAND(IA,ISTBNU(INUC)) .EQ. 0) THEN
               IAT = IAT + 1
               NUCZ = NINT(CHARGE(INUC))
               X = XTANG * PT(IBTAND(ISYMAX(1,1),IA)) * CORD(1,INUC)
               Y = XTANG * PT(IBTAND(ISYMAX(2,1),IA)) * CORD(2,INUC)
               Z = XTANG * PT(IBTAND(ISYMAX(3,1),IA)) * CORD(3,INUC)
               WRITE(LUIP,8020) ATMSYM(NUCZ),X,Y,Z
            END IF
 311     CONTINUE 
 310  CONTINUE
C
C     ----- WE CANNOT GUESS WHAT PLANE SHOULD BE PLOTTED -----
C
      WRITE(LUIP,8030)
C
C     ----- PUNCH LIST OF BONDED ATOMS -----
C
      CALL GEOANA(CORD,.FALSE.,.FALSE.,NBONDS,-1,WORK,LWORK)
C
C     ----- WE CANNOT GUESS WHAT MO-S SHOULD BE PLOTTED -----
C
      WRITE(LUIP,8060)
C
C     ----- PUNCH THE ATOMIC BASIS SET -----
C
      DO 280 ISHELL = 1, KMAX
         ICENT = NCENT(ISHELL)
         NDEG = NUCDEG(ICENT)
         NUCZ = NINT(CHARGE(ICENT))
         DO 281 IDEG = 1, NDEG
            NSH = 1
            NNUC = NUCO(ISHELL)
            NNRC = NRCO(ISHELL)
            KK = NUMCF(ISHELL)
            IPRIM = JSTRT(ISHELL)
            J = NHKT(ISHELL)
C
C     Write out the number of shells (blocks) on this atom
C
            IF (KK .EQ. 1) WRITE(LUIP,8070) ATMSYM(NUCZ), NSH
C     
            TYPE = ' '
            IF(J.EQ.1) TYPE='S'
            IF(J.EQ.2) TYPE='P'
            IF(J.EQ.3) TYPE='D'
            IF(J.EQ.4) TYPE='F'
            IF(J.EQ.5) TYPE='G'
            IF(TYPE.EQ.' ') THEN
               WRITE(LUPRI,'(/A/A)')
     &            'PLTORB knows S,P,D,F,G shells. The basis set'//
     &            ' contains orbital(s) with higher angular momentum,',
     &            ' and the writing of MOLPLT file is thus abandoned.'
               CALL GPCLOSE(LUIP,'DELETE')
               GOTO 5000
            END IF
C     
C     
            IF (KK .EQ. 1) THEN
               WRITE(LUIP,8080) TYPE,NNUC,NNRC
               WRITE(LUIP,8090) (PRIEXP(IPRIM + NN),NN=1,NNUC)
               DO 260 NN = 1, NNRC
                  WRITE (LUIP,8090) (PRICCF(IPRIM+MM,NN),MM=1,NNUC)
 260           CONTINUE 
            END IF
 281     CONTINUE 
 280  CONTINUE 
C
C     ----- WE CANNOT GUESS WHAT THE ORBITAL TITLE CARDS ARE -----
C
      WRITE(LUIP,8100)
      WRITE(LUIP,8110)
      CALL GPCLOSE(LUIP,'KEEP')
C
 5000 CONTINUE
      RETURN
C
 8000 FORMAT('------ START OF -PLTORB- INPUT FILE -------')
 8010 FORMAT('NATOMS=',I4,'   NBONDS=',I4,'   NAOS=',I4,
     *       '   NMOS=??   NPLOTS=?? '/
     *       'PLANE=????????  KOLOR=1  ANGSTROMS FORMAT (5X,5F15.10)')
 8020 FORMAT(A4,1X,3F20.10)
 8030 FORMAT('?? DESIRED PLOTTING PLANE GOES HERE...'/
     *       '?? PLOTTING PLANE BOUNDARIES GO HERE...')
 8060 FORMAT('PLOTMOS ?? ')
 8070 FORMAT(A4,I5)
 8080 FORMAT(A1,2X,2I4)
 8090 FORMAT(1P,5E14.6,:,' >')
 8100 FORMAT('?? ORBITAL TITLE CARDS GO HERE...')
 8110 FORMAT('------ END OF THE -PLTORB- INPUT FILE ------')
c      END
C/* Deck PUSQL */
c      SUBROUTINE PUSQL(V,M,N,NDIM)
C
C     Original source code provided by Kim Baldridge at SDSC. Punches
C     molecular orbital coefficients
C     Modified for use in DALTON by K.Ruud, September 1995
C
c#include <implicit.h>
C
c      DIMENSION V(NDIM,M)
C
C     ----- PUNCH A RECTANGULAR MATRIX WITH ORDERING LABELS -----
C     -V- IS -N- ROWS BY -M- COLUMNS, WITH TRUE LEAD DIMENSION -NDIM-
C
c      DO 120 J = 1,M
c      IC = 0
c      MAX = 0
c  100 MIN = MAX+1
c      MAX = MAX+5
c      IC = IC+1
c      IF (MAX .GT. N) MAX = N
c      MODJ=MOD(J,100)
c      WRITE (IP,9008) MODJ,IC,(V(I,J),I = MIN,MAX)
c      IF (MAX .LT. N) GO TO 100
c  120 CONTINUE
c      RETURN
C
c 9008 FORMAT(I2,I3,1P,5E15.8)
      END
