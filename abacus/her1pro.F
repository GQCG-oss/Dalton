C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
C
C$Id: her1pro.F,v 1.10 2001-05-31 06:38:05 vebjornb Exp $
C
#include <single.h>
C  /* Deck pr1int */
      SUBROUTINE PR1INT(WORD,WORK,LWORK,IORDER,NPQUAD,
     &                  TRIANG,PROPRI,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxmom.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <iratdef.h>
C
#include <huckel.h>
#include <nuclei.h>
#include <shells.h>
#include <symmet.h>
C
      LOGICAL TRIANG, PROPRI, TOFILE, DOINT(2,2)
      CHARACTER WORD*7, LABINT(3*MXCOOR)*8
      DIMENSION WORK(LWORK)
C
C     MAXTYP is needed for allocation of INTREP AND INTADR
C     The actual number of elements in these arrays is later
C     tested for in NTYPTS.
C
      MAXTYP = 3*MXCOOR
C
C     For DSO integrals the following number of elements are
C     needed. NTYPTS is therefore not called.
C
      IF (WORD .EQ. 'DSO    ') MAXTYP = (3*NUCDEP)**2
C
      NCOMP    = 0
      TOFILE   = .TRUE.
      DOINT(1,1) = .TRUE.
      IF (ADDSTO .AND. WORD .EQ. 'OVERLAP') THEN
         DOINT(2,1) = .TRUE.
         DOINT(1,2) = .TRUE.
         DOINT(2,2) = .TRUE.
         KMAX   = NLRGSH + NSMLSH
         NBASIS = NLARGE + NSMALL
         NPBAS  = NPLRG + NPSML
C
         KKVAL = 1
         KMVAL = KKVAL + MXAQN
         KNVAL = KMVAL + MXAQN
         KIREP = KNVAL + MXAQN
         KLAST = KIREP + MXCORB
         CALL SYMPRO(WORK(KKVAL),WORK(KMVAL),WORK(KNVAL),WORK(KIREP),
     &               .FALSE.)
         CALL ICOPY(8,NCOS(1,2),1,NHUCAO(1),1)
         NHUCBA = 0
         DO 14 I = 1, MAXREP + 1
            NHUCBA = NHUCBA + NHUCAO(I)
 14      CONTINUE 
      ELSE
         DOINT(2,1) = .FALSE.
         DOINT(1,2) = .FALSE.
         DOINT(2,2) = .FALSE.
      END IF
C
      KLBINT = 1
      IF (WORD .EQ. 'DSO    ') THEN
         KINTRP = KLBINT +  (3*NUCDEP)**2
      ELSE
         KINTRP = KLBINT +  3*MXCOOR
      END IF
      KINTAD = KINTRP + (MAXTYP + 1)/IRAT
      IF (WORD .EQ. 'ELFGRDC' .OR. WORD .EQ. 'ELFGRDS') THEN
         KLAST = KINTAD + 9*NUCIND*(MAXREP + 1)
      ELSE
         KLAST = KINTAD + (MAXTYP + 1)/IRAT
      END IF
      IF (KLAST .GT. LWORK) CALL STOPIT('PR1INT',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL PR1IN1(WORK(KLAST),LWRK,WORK(KINTRP),WORK(KINTAD),
     &            WORK(KLBINT),WORD,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,DUMMY,NCOMP,TOFILE,'TRIANG',
     &            DOINT)
      IF (ADDSTO .AND. WORD .EQ. 'OVERLAP') THEN
         KMAX   = NLRGSH
         NBASIS = NLARGE
         NPBAS  = NPLRG
C
         KKVAL = 1
         KMVAL = KKVAL + MXAQN
         KNVAL = KMVAL + MXAQN
         KIREP = KNVAL + MXAQN
         KLAST = KIREP + MXCORB
         CALL SYMPRO(WORK(KKVAL),WORK(KMVAL),WORK(KNVAL),WORK(KIREP),
     &               .FALSE.)
      END IF
      RETURN
      END
C  /* Deck pr1in1 */
      SUBROUTINE PR1IN1(WORK,LWORK,INTREP,INTADR,LABINT,WORD,
     &                  IORDER,NPQUAD,TRIANG,PROPRI,IPRINT,
     &                  SINTMA,NCOMP,TOFILE,MTFORM,DOINT)
C
C     Calculation of one-electron property integrals
C
C     T. Helgaker
C
C     Overlap integrals (28.06.89) (OVERLAP)
C     Dipole integrals (28.06.89) (DIPLEN)
C     Spatial one-electron spin-orbit integrals (23.11.89) (SPNORB)
C     Dipole velocity integrals (17.01.90) (DIPVEL)
C     Quadrupole integrals (17.01.90) (QUADRUP)
C     Cartesian moments integrals (all orders) (28.09.90) (CARMOM)
C     Spherical moments integrals (all orders) (20.10.90) (SPHMOM)
C     One-electron Fermi contact (07.02.91)
C     Paramagnetic spin-orbit integrals (09.02.91)
C     Spin-dipole integrals (10.02.91)
C     Diamagnetic spin-orbit integrals (11.02.91)
C     Half-derivative overlap integrals for 1st-order NACMEs (25.06.91)
C     Cosine and Sine integrals (24.06.93)
C     Mass-velocity and Darwin integrals (23.07.93 ShKi+HJAaJ)
C     Magnetic field derivatives of electric field (280893 KRu)
C
C     Changes for writing symmetry information to property file
C     (OV 08.03.90)
C
C     Solvent flag (22.01.91/HJAAJ+KM) (SOLVENT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
#include <iratdef.h>
#include <efield.h>
      PARAMETER (NTABLE = 72)
      PARAMETER (D2 = 2.D0)
C
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, PROPRI, SQUARE,
     &        TOFILE,DOINT(2,2)
      CHARACTER LABINT(*)*8, WORD*7, TABLE(NTABLE)*7, MTFORM*6,
     &          EFDIR*1
      CHARACTER*8  RTNLBL(2)
      DIMENSION WORK(LWORK), INTREP(*), INTADR(*), SINTMA(*)
C
#include <shells.h>
#include <symmet.h>
#include <pgroup.h>
#include <nuclei.h>
#include <huckel.h>
#include <cbihr1.h>
#include <orgcom.h>
C
      DATA TABLE /'OVERLAP', 'DIPLEN ', 'DIPVEL ', 'QUADRUP',
     &            'SPNORB ', 'SECMOM ', 'THETA  ', 'CARMOM ',
     &            'SPHMOM ', 'SOLVENT', 'FERMI C', 'PSO    ',
     &            'SPIN-DI', 'DSO    ', 'SDFC   ', 'HDO    ',
     &            'S1MAG  ', 'S2MAG  ', 'ANGLON ', 'ANGMOM ',
     &            'LONMOM ', 'MAGMOM ', 'KINENER', 'DSUSNOL',
     &            'DSUSLAN', 'DSUSLH ', 'DIASUS ', 'NUCSNLO',
     &            'NUCSLO ', 'NUCSHI ', 'NEFIELD', 'ELFGRDC',
     &            'ELFGRDS', 'S1MAGL ', 'S1MAGR ', 'HDOBR  ',
     &            'NUCPOT ', 'HBDO   ', 'SQHDO  ', 'DSUSCGO',
     &            'NSTCGO ', 'EXPIKR ', 'MASSVEL', 'DARWIN ',
     &            'CM1    ', 'CM2    ', 'SQHDOR ', 'SQOVLAP',
     &            'LONSOL1', 'LONSOL2', 'NSTCGOS', 'S1ELE  ', 
     &            'S1ELB  ', 'ONEELD ', 'DPLGRA ', 'QUAGRA ',
     &            'OCTGRA ', 'ROTSTR ', 'THRMOM ', 'SOFIELD',
     &            'SOMAGMO', 'DEROVLP', 'DERHAMI', 'ELGDIAN',
     &            'ELGDIAL', 'DPTOVL ', 'DPTPO1 ', 'DPTPO2 ',
     &            'XDDXR3 ', 'PVPINT ', 'POTENER', 'SPHMOML'/
C
#include <ibtfun.h>
C
      KATOM = 1
      KLAST = KATOM + NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('PR1IN1','ATM',KLAST,LWORK)
C
C     **************************
C     ***** Integral types *****
C     **************************
C
      TRASPH = .FALSE.
      SQUARE = .FALSE.
      SOLVNT = .FALSE.
      DO 100 I = 1, NTABLE
         IF (TABLE(I) .EQ. WORD) THEN
            GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &             19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,
     &             34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
     &             49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
     &             64,65,66,67,68,69,70,71,72), I
         END IF
  100 CONTINUE
      WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *    '" not recognized in PR1IN1.'
      CALL PRTAB(NTABLE,TABLE,'PR1IN1 input keywords',LUPRI)
      CALL QUIT('Illegal keyword in PR1IN1.')
C
C     Overlap integrals
C     -----------------
C
    1 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of overlap integrals.'
         INTTYP = 1
         NOPTYP = 1
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'OVERLAP '
         INTREP(1) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Dipole moment (dipole length) integrals
C     ---------------------------------------
C
    2 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of dipole moment (length) integrals.'
         INTTYP = 2
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XDIPLEN '
         LABINT(2) = 'YDIPLEN '
         LABINT(3) = 'ZDIPLEN '
         INTREP(1) = ISYMAX(1,1)
         INTREP(2) = ISYMAX(2,1)
         INTREP(3) = ISYMAX(3,1)
         ANTI = .FALSE.
      GO TO 200
C
C     Dipole velocity integrals
C     -------------------------
C
    3 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of dipole velocity integrals.'
         INTTYP = 3
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XDIPVEL '
         LABINT(2) = 'YDIPVEL '
         LABINT(3) = 'ZDIPVEL '
         INTREP(1) = ISYMAX(1,1)
         INTREP(2) = ISYMAX(2,1)
         INTREP(3) = ISYMAX(3,1)
         ANTI = .TRUE.
      GO TO 200
C
C     Quadrupole moment integrals
C     ---------------------------
C
    4 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of quadrupole moment integrals.'
         INTTYP = 4
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXQUADRU'
         LABINT(2) = 'XYQUADRU'
         LABINT(3) = 'XZQUADRU'
         LABINT(4) = 'YYQUADRU'
         LABINT(5) = 'YZQUADRU'
         LABINT(6) = 'ZZQUADRU'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Spin-orbit integrals
C     --------------------
C
    5 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of spatial one-electron spin-orbit integrals.'
         INTTYP = 5
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'X1SPNORB'
         LABINT(2) = 'Y1SPNORB'
         LABINT(3) = 'Z1SPNORB'
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .TRUE.
      GO TO 200
C
C     Second moments integrals
C     ------------------------
C
    6 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of second moments integrals.'
         INTTYP = 6
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXSECMOM'
         LABINT(2) = 'XYSECMOM'
         LABINT(3) = 'XZSECMOM'
         LABINT(4) = 'YYSECMOM'
         LABINT(5) = 'YZSECMOM'
         LABINT(6) = 'ZZSECMOM'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Traceless theta quadrupole integrals
C     ------------------------------------
C
    7 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of theta quadrupole moments integrals.'
         INTTYP = 7
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXTHETA '
         LABINT(2) = 'XYTHETA '
         LABINT(3) = 'XZTHETA '
         LABINT(4) = 'YYTHETA '
         LABINT(5) = 'YZTHETA '
         LABINT(6) = 'ZZTHETA '
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Cartesian moments integrals
C     ---------------------------
C
    8 CONTINUE
         INTTYP = 8
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I2,A/)')
     &      ' Calculation of Cartesian multipole moment '/
     &      /'integrals of order',IORDER,'.'
         IF (IORDER .GT. MXQNM - 1) THEN
            WRITE (LUPRI,'(2X,A)')
     &         ' Maximum multipole moment order exceeded in PR1IN1.'
            WRITE (LUPRI,'(2X,A,I5,/,A,I5,/,A,I3)')
     &         ' Order requested:',IORDER,
     &         ' Maximum order:  ',MXQNM-1,
     &         ' Increase MXQNM to',IORDER + 1,' and recompile.'
            CALL QUIT('Multipole moment order exceeded in PR1IN1.')
         END IF
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL CMTYP(IORDER,NOPTYP,INTREP,LABINT)
         ANTI = .FALSE.
      GO TO 200
C
C     Spherical moments integrals
C     ---------------------------
C
    9 CONTINUE
         INTTYP = 8
         TRASPH = .TRUE.
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I2,A/)')
     &      ' Calculation of spherical multipole moment '/
     &      /'integrals of order',IORDER,'.'
         IF (IORDER .GT. MXQNM - 1) THEN
            WRITE (LUPRI,'(2X,A)')
     &         ' Maximum multipole moment order exceeded in PR1IN1.'
            WRITE (LUPRI,'(2X,A,I5,/,A,I5,/,A,I3)')
     &         ' Order requested:',IORDER,
     &         ' Maximum order:  ',MXQNM-1,
     &         ' Increase MXQNM to',IORDER + 1,' and recompile.'
            CALL QUIT('Multipole moment order exceeded in PR1IN1.')
         END IF
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL CMTYP(IORDER,NOPTYP,INTREP,LABINT)
         ANTI = .FALSE.
      GO TO 200
C
C     Electronic solvent integrals
C     ----------------------------
C
   10 CONTINUE
         INTTYP = 8
         TRASPH = .TRUE.
         SOLVNT = .TRUE.
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I2,/A/)')
     &      ' Calculation of spherical multipole moment '/
     &      /'integrals of order',IORDER,
     &       ' for cavity solvent model.'
         IF (IORDER .GT. MXQNM - 1) THEN
            WRITE (LUPRI,'(2X,A)')
     &         ' Maximum multipole moment order exceeded in PR1IN1.'
            WRITE (LUPRI,'(2X,A,I5,/,A,I5,/,A,I3)')
     &         ' Order requested:',IORDER,
     &         ' Maximum order:  ',MXQNM-1,
     &         ' Increase MXQNM to',IORDER + 1,' and recompile.'
            CALL QUIT('Multipole moment order exceeded in PR1IN1.')
         END IF
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL CMTYP(IORDER,NOPTYP,INTREP,LABINT)
         ANTI = .FALSE.
      GO TO 200
C
C     One-electron Fermi contact integrals
C     ------------------------------------
C
   11 CONTINUE
         INTTYP = 9
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of one-electron Fermi contact integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL FRMTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM)
         ANTI = .FALSE.
      GO TO 200
C
C     Paramagnetic spin-orbit integrals
C     ---------------------------------
C
   12 CONTINUE
         INTTYP = 10
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of paramagnetic spin-orbit integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL PSOTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),INTADR,NATOM)
         ANTI = .TRUE.
      GO TO 200
C
C     Spin-dipole integrals
C     ---------------------
C
   13 CONTINUE
         INTTYP = 11
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of spin-dipole integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL SDTYP(INTTYP,NOPTYP,INTREP,LABINT,INTADR,WORK(KATOM),
     &              NATOM)
         ANTI = .FALSE.
      GO TO 200
C
C     Diamagnetic spin-orbit integrals
C     --------------------------------
C
   14 CONTINUE
         INTTYP = 12
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of diamagnetic spin-orbit integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL DSOTYP(NOPTYP,INTREP,LABINT,INTADR,WORK(KATOM),NATOM,
     &               TRIANG)
         ANTI = .FALSE.
      GO TO 200
C
C     Spin-dipole + Fermi contact integrals
C     -------------------------------------
C
   15 CONTINUE
         INTTYP = 13
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of spin-dipole + Fermi contact integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL SDTYP(INTTYP,NOPTYP,INTREP,LABINT,INTADR,WORK(KATOM),
     &              NATOM)
         ANTI = .FALSE.
      GO TO 200
C
C     Half-derivative overlap integrals
C     ---------------------------------
C
   16 CONTINUE
         INTTYP = 14
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of half-derivative overlap integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL HDOTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,.TRUE.,
     &               INTTYP)
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GO TO 200
C
C     Contribution from overlap matrix to magnetic properties
C     -------------------------------------------------------
C
 17   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of first magnetic derivative of overlap matrix'
         INTTYP = 15
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'dS/dBX  '
         LABINT(2) = 'dS/dBY  '
         LABINT(3) = 'dS/dBZ  '
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(3,1),ISYMAX(1,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .TRUE.
      GOTO 200
C
C     Second order contribution from overlap matrix to magnetic properties
C     --------------------------------------------------------------------
C
 18   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &    ' Calculation of second magnetic derivative of overlap matrix'
         INTTYP = 16
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'dS/dB2XX'
         LABINT(2) = 'dS/dB2XY'
         LABINT(3) = 'dS/dB2XZ'
         LABINT(4) = 'dS/dB2YY'
         LABINT(5) = 'dS/dB2YZ'
         LABINT(6) = 'dS/dB2ZZ'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
      GOTO 200
C
C     Electronic angular momentum around the nuclei
C     ---------------------------------------------
C
 19   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &      ' Calculation of angular momentum around the nuclei'
         INTTYP = 17
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XANGLON '
         LABINT(2) = 'YANGLON '
         LABINT(3) = 'ZANGLON '
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GO TO 200
C
C     Electronic angular momentum around the origin
C     ---------------------------------------------
C
 20   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &      ' Calculation of angular momentum around the gauge origin'
         INTTYP = 18
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XANGMOM '
         LABINT(2) = 'YANGMOM '
         LABINT(3) = 'ZANGMOM '
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .TRUE.
      GO TO 200
C
C     London orbital contribution to angular momentum
C     -----------------------------------------------
C
 21   CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &  ' Calculation of London orbital contribution to magnetic moment'
         INTTYP = 19
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XLONMOM '
         LABINT(2) = 'YLONMOM '
         LABINT(3) = 'ZLONMOM '
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GO TO 200
C
C     One-electron contribution to magnetic moment
C     --------------------------------------------
C
 22   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &   ' Calculation of one-electron contribution to magnetic moment'
         INTTYP = 20
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'dh/dBX  '
         LABINT(2) = 'dh/dBY  '
         LABINT(3) = 'dh/dBZ  '
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .TRUE.
      GO TO 200
C
C     Electronic kinetic energy
C     -------------------------
C
 23   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &   ' Calculation of electronic kinetic energy'
         INTTYP = 21
         NOPTYP = 1
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'KINENERG'
         INTREP(1) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Diamagnetic susceptiblity without London orbital contribution
C     -------------------------------------------------------------
C
 24   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &     'Calculation of diamagnetic susceptiblity with no London'/
     &     /'orbital contribution'
         INTTYP = 22
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXDSUSNL'
         LABINT(2) = 'XYDSUSNL'
         LABINT(3) = 'XZDSUSNL'
         LABINT(4) = 'YYDSUSNL'
         LABINT(5) = 'YZDSUSNL'
         LABINT(6) = 'ZZDSUSNL'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GO TO 200
C
C     Angular London orbital contribution to diamagnetic susceptibility
C     -----------------------------------------------------------------
C
 25   CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      'Calculation of angular london orbital contribution to'/
     &      /'magnetic susceptibility'
         INTTYP = 23
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXDSUSLL'
         LABINT(2) = 'XYDSUSLL'
         LABINT(3) = 'XZDSUSLL'
         LABINT(4) = 'YYDSUSLL'
         LABINT(5) = 'YZDSUSLL'
         LABINT(6) = 'ZZDSUSLL'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Angular London orbital contribution to diamagnetic susceptibility
C     -----------------------------------------------------------------
C
 26   CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      'Calculation of London orbital contribution to'/
     &      /'magnetic susceptibility'
         INTTYP = 24
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXDSUSLH'
         LABINT(2) = 'XYDSUSLH'
         LABINT(3) = 'XZDSUSLH'
         LABINT(4) = 'YYDSUSLH'
         LABINT(5) = 'YZDSUSLH'
         LABINT(6) = 'ZZDSUSLH'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Angular London orbital contribution to diamagnetic susceptibility
C     -----------------------------------------------------------------
C
 27   CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      'Calculation of diamagnetic susceptibility'
         INTTYP = 25
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXdh/dB2'
         LABINT(2) = 'XYdh/dB2'
         LABINT(3) = 'XZdh/dB2'
         LABINT(4) = 'YYdh/dB2'
         LABINT(5) = 'YZdh/dB2'
         LABINT(6) = 'ZZdh/dB2'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
      GOTO 200
C
C     Nuclear shielding integrals without London orbital contribution
C     ---------------------------------------------------------------
C
 28   CONTINUE
         INTTYP = 26
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &        ' Calculation of nuclear shieldings without London contr.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL NSTTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,'NSNL',
     &               INTADR)
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     London orbital contribution to nuclear shielding tensor integrals
C     -----------------------------------------------------------------
C
 29   CONTINUE
         INTTYP = 27
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &       ' Calculation of London contribution to nuclear shieldings'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL NSTTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,'NSLO',
     &               INTADR)
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Nuclear shielding tensor integrals
C     ----------------------------------
C
 30   CONTINUE
         INTTYP = 28
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &        'Calculation of nuclear shielding tensor integrals'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL NSTTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,' NST',
     &               INTADR)
         ANTI = .FALSE.
      GOTO 200
C
C     Electric field at the individual nuclei
C     --------------------------------------
C
 31   CONTINUE
         INTTYP = 29
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &        'Calculation of electric field strength at the nuclei'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL EFNTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,INTADR,
     &               INTTYP)
         ANTI = .FALSE.
      GOTO 200
C
C     Electric field gradient at the individual nuclei, cartesian
C     -----------------------------------------------------------
C
 32   CONTINUE
         INTTYP = 30
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &        ' Calculation of electric field gradients (cartesian)'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL EFGTYP(NOPTYP,INTREP,INTADR,LABINT,WORK(KATOM),NATOM)
         ANTI = .FALSE.
      GOTO 200
C
C     Electric field gradient at the individual nuclei, spherical
C     -----------------------------------------------------------
C
 33   CONTINUE
         INTTYP = 31
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &        ' Calculation of electric field gradients (spherical)'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL EFGTYP(NOPTYP,INTREP,INTADR,LABINT,WORK(KATOM),NATOM)
         ANTI = .FALSE.
      GOTO 200
C
C     Bra-differentiation of overlap matrix with respect to magnetic field
C     --------------------------------------------------------------------
C
 34   CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     ' Calculation of bra-magnetic derivative of overlap matrix'
         INTTYP = 32
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'd<S|/dBX'
         LABINT(2) = 'd<S|/dBY'
         LABINT(3) = 'd<S|/dBZ'
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(3,1),ISYMAX(1,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Ket-differentiation of overlap matrix with respect to magnetic field
C     --------------------------------------------------------------------
C
 35   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of ket-magnetic derivative of magnetic field'
         INTTYP = 33
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'd|S>/dBX'
         LABINT(2) = 'd|S>/dBY'
         LABINT(3) = 'd|S>/dBZ'
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(3,1),ISYMAX(1,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Ket-differentiation of HDO-integrals with respect to magnetic field
C     -------------------------------------------------------------------
C
 36   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of ket-magnetic derivative of HDO-integrals'
         INTTYP = 34
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL HDBTYP(NOPTYP,INTREP,INTADR,LABINT,WORK(KATOM),NATOM)
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Potential energy at the nuclei
C     -----------------------------
C
 37   CONTINUE
         INTTYP = 35
         IF (MAXREP .GT. 0) THEN
            WRITE (LUPRI,'(/A/)')
     &           ' Program cannot calculate potential energy '//
     &           'with symmetry'
            CALL QUIT('Cannot calculate potential energy with symmetry')
         ELSE
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &           ' Calculation of the potential energy at the nuclei'
            CALL SETATM(WORK(KATOM),NATOM,INTTYP)
            CALL NPETYP(NOPTYP,INTREP,INTADR,LABINT,WORK(KATOM),NATOM)
            ANTI = .FALSE.
         END IF
      GOTO 200
C
C     Half B-differentiated overlap matrix
C     ------------------------------------
C
 38   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of half B-differentiated overlap matrix'
         INTTYP = 36
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = ' HBDO X '
         LABINT(2) = ' HBDO Y '
         LABINT(3) = ' HBDO Z '
         INTREP(1) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(2) = IBTXOR(ISYMAX(3,1),ISYMAX(1,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         ANTI = .FALSE.
      GOTO 200
C
C     Half-derivative overlap integrals not to be antisymmetrized
C     -----------------------------------------------------------
C
 39   CONTINUE
         INTTYP = 14
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of half-derivative overlap integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL HDOTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,.FALSE.,
     &               INTTYP)
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Diamagnetic susceptiblity with common gauge origin
C     --------------------------------------------------
C
 40   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &     'Calculation of diamagnetic susceptiblity with common '/
     &     /'gauge origin'
         INTTYP = 37
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXSUSCGO'
         LABINT(2) = 'XYSUSCGO'
         LABINT(3) = 'XZSUSCGO'
         LABINT(4) = 'YYSUSCGO'
         LABINT(5) = 'YZSUSCGO'
         LABINT(6) = 'ZZSUSCGO'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Nuclear shielding integrals with common gauge origin
C     ----------------------------------------------------
C
 41   CONTINUE
         INTTYP = 38
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     ' Calculation of nuclear shieldings with common gauge origin'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL NSTTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,'NSCO',
     &               INTADR)
         ANTI = .FALSE.
      GOTO 200
C
C     Cosine and Sine integrals
C     -------------------------
C
   42 CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of cos(kr)/k and sin(kr)/k integrals'
         INTTYP = 39
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'COS KX/K'
         LABINT(2) = 'COS KY/K'
         LABINT(3) = 'COS KZ/K'
         LABINT(4) = 'SIN KX/K'
         LABINT(5) = 'SIN KY/K'
         LABINT(6) = 'SIN KZ/K'
         INTREP(1) = 0
         INTREP(2) = 0
         INTREP(3) = 0
         INTREP(4) = ISYMAX(1,1)
         INTREP(5) = ISYMAX(2,1)
         INTREP(6) = ISYMAX(3,1)
         ANTI = .FALSE.
      GOTO 200
C
C     Mass velocity integrals
C     -----------------------
C
  43  CONTINUE
         INTTYP = 40
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of mass velocity integrals'
         NOPTYP = 1
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'MASSVELO'
         INTREP(1) = 0
         ANTI = .FALSE.
      GOTO 200
C
C     Darwin type integrals
C     ---------------------
C
   44 CONTINUE
         INTTYP = 41
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of Darwin type integrals'
         NOPTYP  = 1
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'DARWIN  '
         INTREP(1) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     First order magnetic field derivatives of electric field
C     --------------------------------------------------------
C
 45   CONTINUE
         INTTYP = 42
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     ' Calculation of first magnetic derivative of electric field'
         NOPTYP  = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         EFDIR = FIELD1(1:1)
         LABINT(1) = EFDIR//'-CM1 X'
         LABINT(2) = EFDIR//'-CM1 Y'
         LABINT(3) = EFDIR//'-CM1 Z'
         IF (EFDIR .EQ. 'X') THEN
            IRET = 1
         ELSE IF (EFDIR .EQ. 'Y') THEN
            IRET = 2
         ELSE
            IRET = 3
         END IF
         INTREP(1) = IBTXOR(ISYMAX(1,2),ISYMAX(IRET,1))
         INTREP(2) = IBTXOR(ISYMAX(2,2),ISYMAX(IRET,1))
         INTREP(3) = IBTXOR(ISYMAX(3,2),ISYMAX(IRET,1))
         ANTI      = .TRUE.
      GOTO 200
C
C     Second order magnetic field derivatives of electric field
C     ---------------------------------------------------------
C
 46   CONTINUE
         INTTYP = 43
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of second magnetic derivatives of electric field'
         NOPTYP  = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         EFDIR      = FIELD2(1:1)
         LABINT(1)  = EFDIR//'-CM2XX'
         LABINT(2)  = EFDIR//'-CM2XY'
         LABINT(3)  = EFDIR//'-CM2XZ'
         LABINT(4)  = EFDIR//'-CM2YY'
         LABINT(5)  = EFDIR//'-CM2YZ'
         LABINT(6)  = EFDIR//'-CM2ZZ'
         IF (EFDIR .EQ. 'X') THEN
            IRET = 1
         ELSE IF (EFDIR .EQ. 'Y') THEN
            IRET = 2
         ELSE
            IRET = 3
         END IF
         INTREP(1)  = ISYMAX(IRET,1)
         INTREP(2)  = IBTXOR(IBTXOR(ISYMAX(1,2),ISYMAX(2,2)),
     &                              ISYMAX(IRET,1))
         INTREP(3)  = IBTXOR(IBTXOR(ISYMAX(1,2),ISYMAX(3,2)),
     &                              ISYMAX(IRET,1))
         INTREP(4)  = ISYMAX(IRET,1)
         INTREP(5)  = IBTXOR(IBTXOR(ISYMAX(2,2),ISYMAX(3,2)),
     &                              ISYMAX(IRET,1))
         INTREP(6)  = ISYMAX(IRET,1)
         ANTI       = .FALSE.
      GOTO 200
C
C     Half-derivative overlap integrals not to be antisymmetrized
C     -----------------------------------------------------------
C
 47   CONTINUE
         INTTYP = 44
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of ket-half-derivative overlap integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL HDOTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,.FALSE.,
     &               INTTYP)
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GOTO 200
C
C     Square overlap integrals
C     ------------------------
C
 48   CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of square overlap integrals.'
         INTTYP = 45
         NOPTYP = 1
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'SQOVLAP'
         INTREP(1) = 0
         ANTI = .FALSE.
         SQUARE = .TRUE.
      GO TO 200
C
C     First derivative of solvent Hamiltonian with respect to magnetic field
C     ----------------------------------------------------------------------
C
 49   CONTINUE
         INTTYP = 46
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I2,A/)')
     &      ' Calculation of magnetic differentiated solvent '/
     &      /'integrals of order',IORDER,'.'
         IF (IORDER .GT. MXQNM - 1) THEN
            WRITE (LUPRI,'(2X,A)')
     &         ' Maximum multipole moment order exceeded in PR1IN1.'
            WRITE (LUPRI,'(2X,A,I5,/,A,I5,/,A,I3)')
     &         ' Order requested:',IORDER,
     &         ' Maximum order:  ',MXQNM-1,
     &         ' Increase MXQNM to',IORDER + 1,' and recompile.'
            CALL QUIT('Multipole moment order exceeded in PR1IN1.')
         END IF
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL SL1TYP(IORDER,NOPTYP,INTREP,LABINT)
         ANTI = .TRUE.
      GO TO 200
C
C     Second derivative of solvent Hamiltonian with respect to magnetic field
C     -----------------------------------------------------------------------
C
 50   CONTINUE
         INTTYP = 47
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I2,A/)')
     &      ' Calculation of magnetic differentiated solvent '/
     &      /'integrals of order',IORDER,'.'
         IF (IORDER .GT. MXQNM - 1) THEN
            WRITE (LUPRI,'(2X,A)')
     &         ' Maximum multipole moment order exceeded in PR1IN1.'
            WRITE (LUPRI,'(2X,A,I5,/,A,I5,/,A,I3)')
     &         ' Order requested:',IORDER,
     &         ' Maximum order:  ',MXQNM-1,
     &         ' Increase MXQNM to',IORDER + 1,' and recompile.'
            CALL QUIT('Multipole moment order exceeded in PR1IN1.')
         END IF
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL SL2TYP(IORDER,NOPTYP,INTREP,LABINT)
         ANTI = .FALSE.
      GO TO 200
C
C     Diamagnetic contribution integrals to rotational g factors
C     ----------------------------------------------------------
C
 51   CONTINUE
         INTTYP = 48
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     ' Calculation of nuclear shieldings (gauge at nucleus)'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL NSTTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,'NSTK',
     &               INTADR)
         ANTI = .FALSE.
      GOTO 200
C
C     Contribution from overlap matrix to electric prop. Type A
C     ---------------------------------------------------------
 52   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of first electric derivative of overlap '/
     $     /'matrix. Type A'
         INTTYP = 49
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'dS/dEXa '
         LABINT(2) = 'dS/dEYa '
         LABINT(3) = 'dS/dEZa '
C
         INTREP(1) = ISYMAX(1,1)
         INTREP(2) = ISYMAX(2,1)
         INTREP(3) = ISYMAX(3,1)
         ANTI = .FALSE.
      GOTO 200
C
C     Contribution from overlap matrix to electric prop. Type B
C     ---------------------------------------------------------
 53   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of first electric derivative of overlap '/
     $     / 'matrix. Type B'
         INTTYP = 50
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'dS/dEXb '
         LABINT(2) = 'dS/dEYb '
         LABINT(3) = 'dS/dEZb '
C
         INTREP(1) = ISYMAX(1,1)
         INTREP(2) = ISYMAX(2,1)
         INTREP(3) = ISYMAX(3,1)
         ANTI = .FALSE.
      GOTO 200
C
C     First electric deriv. of one-electron Hamiltonian integrals
C     ___________________________________________________________
C
 54   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     $     'Calculation of electric derivative of 1-electron
     $     Hamiltonian integrals'
         INTTYP = 51
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'dh1/dEX '
         LABINT(2) = 'dh1/dEY '
         LABINT(3) = 'dh1/dEZ '
C
         INTREP(1) = ISYMAX(1,1)
         INTREP(2) = ISYMAX(2,1)
         INTREP(3) = ISYMAX(3,1)
         ANTI = .FALSE.
      GOTO 200
C
C     Dipole gradient integrals
C     -------------------------
C
 55   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of dipole gradient integrals'
         INTTYP = 52
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL DPGTYP(NOPTYP,INTREP,INTADR,LABINT,WORK(KATOM),NATOM)
CCH
Cori     ANTI = .TRUE.
         ANTI = .FALSE.
CCH
         SQUARE = .FALSE.
      GOTO 200                             
C
C     Quadrupole gradient integrals
C     -----------------------------
C
 56   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of quadrupole gradient integrals'
         INTTYP = 53
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL QUGTYP(NOPTYP,INTREP,INTADR,LABINT,WORK(KATOM),NATOM)
CCH
Cori     ANTI = .TRUE.
         ANTI = .FALSE.
CCH
         SQUARE = .FALSE.
      GOTO 200                              
C
C     Octupole gradient integrals
C     ---------------------------
C
 57   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of octupole gradient integrals'
         INTTYP = 54
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL OCGTYP(NOPTYP,INTREP,INTADR,LABINT,WORK(KATOM),NATOM)
CCH
Cori     ANTI = .TRUE.
         ANTI = .FALSE.
CCH
         SQUARE = .FALSE.
      GOTO 200                             
C
C     Rotational strength integrals
C     -----------------------------
C
 58   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &   ' Calculation of rotational strengths integrals'
         INTTYP = 55
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'XXROTSTR'
         LABINT(2) = 'XYROTSTR'
         LABINT(3) = 'XZROTSTR'
         LABINT(4) = 'YYROTSTR'
         LABINT(5) = 'YZROTSTR'
         LABINT(6) = 'ZZROTSTR'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI = .TRUE.
         SQUARE = .FALSE.
      GOTO 200                                   
C
C     Third moments integrals
C     -----------------------
C
 59   CONTINUE
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &      ' Calculation of third moments integrals.'
         INTTYP = 56
         NOPTYP = 10
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1)  = 'XXX 3MOM'
         LABINT(2)  = 'XXY 3MOM'
         LABINT(3)  = 'XXZ 3MOM'
         LABINT(4)  = 'XYY 3MOM'
         LABINT(5)  = 'XYZ 3MOM'
         LABINT(6)  = 'XZZ 3MOM'
         LABINT(7)  = 'YYY 3MOM'
         LABINT(8)  = 'YYZ 3MOM'
         LABINT(9)  = 'YZZ 3MOM'
         LABINT(10) = 'ZZZ 3MOM'
         INTREP(1)  = ISYMAX(1,1)
         INTREP(2)  = ISYMAX(2,1)
         INTREP(3)  = ISYMAX(3,1)
         INTREP(4)  = ISYMAX(1,1)
         INTREP(5)  =IBTXOR(ISYMAX(1,1),IBTXOR(ISYMAX(2,1),ISYMAX(3,1)))
         INTREP(6)  = ISYMAX(1,1)
         INTREP(7)  = ISYMAX(2,1)
         INTREP(8)  = ISYMAX(3,1)
         INTREP(9)  = ISYMAX(2,1)
         INTREP(10) = ISYMAX(3,1)
         ANTI = .FALSE.
         SQUARE = .FALSE.
      GO TO 200                       
C
C     Magnetic-field correction to spin-orbit integrals
C     -------------------------------------------------
C
 60   CONTINUE
         INTTYP = 57
         NOPTYP = 9
         CALL NTYPTS(NOPTYP)
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     ' Magnetic-field correction to spin-orbit integrals'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'SOMF  XX'
         LABINT(2) = 'SOMF  XY'
         LABINT(3) = 'SOMF  XZ'
         LABINT(4) = 'SOMF  YX'
         LABINT(5) = 'SOMF  YY'
         LABINT(6) = 'SOMF  YZ'
         LABINT(7) = 'SOMF  ZX'
         LABINT(8) = 'SOMF  ZY'
         LABINT(9) = 'SOMF  ZZ'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,2),ISYMAX(2,2))
         INTREP(3) = IBTXOR(ISYMAX(1,2),ISYMAX(3,2))
         INTREP(4) = IBTXOR(ISYMAX(1,2),ISYMAX(2,2))
         INTREP(5) = 0
         INTREP(6) = IBTXOR(ISYMAX(2,2),ISYMAX(3,2))
         INTREP(7) = IBTXOR(ISYMAX(1,2),ISYMAX(3,2))
         INTREP(8) = IBTXOR(ISYMAX(2,2),ISYMAX(3,2))
         INTREP(9) = 0
         ANTI = .FALSE.
         SQUARE = .FALSE.
      GOTO 200                                      
C
C     Magnetic-moment correction to spin-orbit integrals
C     --------------------------------------------------
C
 61   CONTINUE
         INTTYP = 58
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     ' Magnetic-moment correction to spin-orbit integrals'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL NSTTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,'SOMM',
     &               INTADR)
         ANTI = .FALSE.
      GOTO 200
C
C     Derivative overlap integrals
C     ----------------------------
C
 62   CONTINUE
         INTTYP = 59
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of derivative overlap integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL HDOTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,.FALSE.,
     &               INTTYP)
         ANTI = .FALSE.
         SQUARE = .FALSE.
      GO TO 200
C
C     Derivative one-electron Hamiltonian integrals
C     ---------------------------------------------
C
 63   CONTINUE
         INTTYP = 70
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of derivative Hamiltonian integrals.'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL EFNTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM,INTADR,
     &               INTTYP)
         ANTI = .FALSE.
         SQUARE = .FALSE.
      GO TO 200                                 
C
C     Diamagnetic one-electron spin-orbit integrals (No-London)
C     ---------------------------------------------------------
C
 64   CONTINUE
         INTTYP = 61
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &        ' Calculation of diamagnetic one-electron SO without '//
     &        'London contribution'
         NOPTYP = 9
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'D1-SO XX'
         LABINT(2) = 'D1-SO XY'
         LABINT(3) = 'D1-SO XZ'
         LABINT(4) = 'D1-SO YX'
         LABINT(5) = 'D1-SO YY'
         LABINT(6) = 'D1-SO YZ'
         LABINT(7) = 'D1-SO ZX'
         LABINT(8) = 'D1-SO ZY'
         LABINT(9) = 'D1-SO ZZ'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(5) = 0
         INTREP(6) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(7) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(8) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(9) = 0
         ANTI      = .FALSE.
      GOTO 200
C
C     Nuclear shielding integrals with London orbital contribution
C     ---------------------------------------------------------------
C
 65   CONTINUE
         INTTYP = 62
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &        ' Calculation of diamagnetic one-electron SO with '//
     &        'London contribution'
         NOPTYP = 9
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'D1-SOLXX'
         LABINT(2) = 'D1-SOLXY'
         LABINT(3) = 'D1-SOLXZ'
         LABINT(4) = 'D1-SOLYX'
         LABINT(5) = 'D1-SOLYY'
         LABINT(6) = 'D1-SOLYZ'
         LABINT(7) = 'D1-SOLZX'
         LABINT(8) = 'D1-SOLZY'
         LABINT(9) = 'D1-SOLZZ'
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(5) = 0
         INTREP(6) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(7) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(8) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(9) = 0
         ANTI = .FALSE.
      GOTO 200
C     Overlap integrals for the small component in dpt ( <opG|opG> )
C     --------------------------------------------------------------
C     ACH     
 66   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     'Calculation of small-component 1 el. overlap '
         INTTYP = 71
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'dd/dxdx '
         LABINT(2) = 'dd/dxdy '
         LABINT(3) = 'dd/dxdz '
         LABINT(4) = 'dd/dydy '
         LABINT(5) = 'dd/dydz '
         LABINT(6) = 'dd/dzdz '        
C
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI=.FALSE.
      GOTO 200
C
C     Pot energy integrals for the small component in dpt ( <opG|V|opG> )
C     ---------------------------------------------------------------.---
C     ACH     
 67   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     'Calculation of small-component 1 el. pot. integrals'
         INTTYP = 72
         NOPTYP = 6
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'DERXXPVP'
         LABINT(2) = 'DERXY+YX'
         LABINT(3) = 'DERXZ+ZX'
         LABINT(4) = 'DERYY   '
         LABINT(5) = 'DERYZ+ZY'
         LABINT(6) = 'DERZZ   '
C
         INTREP(1) = 0
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(3) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(4) = 0
         INTREP(5) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         INTREP(6) = 0
         ANTI=.FALSE.
      GOTO 200
C
C     Pot energy integrals for the small component in dpt ( <opG|V|opG> )
C     ---------------------------------------------------------------.---
C     ACH
 68   CONTINUE
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     'Calculation of small-component 1 el. pot. integrals'
         INTTYP = 73
         NOPTYP = 3
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'DERXY-YX'
         LABINT(2) = 'DERXZ-ZX'
         LABINT(3) = 'DERYZ-ZY'
C
         INTREP(1) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
         INTREP(2) = IBTXOR(ISYMAX(1,1),ISYMAX(3,1))
         INTREP(3) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
         ANTI=.TRUE.
      GOTO 200
C
C     DPT pso-lookalike integrals <(op)(oA)>
C     ---------------------------------
C
 69   CONTINUE
         INTTYP = 74
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A)')
     &      ' Calculation of magnetic-dpt integrals ( <(op)(oA)> )'
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL ALFTYP(INTTYP,NOPTYP,INTREP,LABINT,INTADR,WORK(KATOM),
     &              NATOM)
C         CALL ALFTYP(NOPTYP,INTREP,LABINT,INTADR,WORK(KATOM),NATOM,
C     &        TRIANG)
         SQUARE = .TRUE.
      GO TO 200
C
C     pVp integrals for Douglas-Kroll transformation
C     ----------------------------------------------
C
 70   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &   ' Calculation of Douglas-Kroll pVp integrals'
         INTTYP = 63
         NOPTYP = 1
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'pVpINTEG'
         INTREP(1) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Electronic kinetic energy
C     -------------------------
C
 71   CONTINUE
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A/)')
     &   ' Calculation of electronic one-electron potential energy'
         INTTYP = 64
         NOPTYP = 1
         CALL NTYPTS(NOPTYP)
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         LABINT(1) = 'POTENERG'
         INTREP(1) = 0
         ANTI = .FALSE.
      GO TO 200
C
C     Spherical moments integrals with local origins
C     ----------------------------------------------
C
   72 CONTINUE
         INTTYP = 8
         TRASPH = .TRUE.
         FMMORI = .TRUE.
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I2,A/)')
     &      ' Calculation of spherical multipole moment '/
     &      /'integrals of order',IORDER,'.'
         IF (IORDER .GT. MXQNM - 1) THEN
            WRITE (LUPRI,'(2X,A)')
     &         ' Maximum multipole moment order exceeded in PR1IN1.'
            WRITE (LUPRI,'(2X,A,I5,/,A,I5,/,A,I3)')
     &         ' Order requested:',IORDER,
     &         ' Maximum order:  ',MXQNM-1,
     &         ' Increase MXQNM to',IORDER + 1,' and recompile.'
            CALL QUIT('Multipole moment order exceeded in PR1IN1.')
         END IF
         CALL SETATM(WORK(KATOM),NATOM,INTTYP)
         CALL CMTYP(IORDER,NOPTYP,INTREP,LABINT)
         ANTI = .FALSE.
      GO TO 200
C
C
  200 CONTINUE
C
C     ***** Print section *****
C
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(/,A,I5)') ' Integral type:', INTTYP
         WRITE (LUPRI,'(/,A,L5)') ' TRASPH:', TRASPH
         WRITE (LUPRI,'(/,A,L5)') ' SQUARE:', SQUARE
         WRITE (LUPRI,'(/,A,L5)') ' SOLVNT:', SOLVNT
         WRITE (LUPRI,'(A,I5)') ' Number of operator components:',NOPTYP
         DO 300 I = 1, NOPTYP
            WRITE (LUPRI,'(2A)')   ' Molecule label: ', LABINT(I)
            WRITE (LUPRI,'(A,I5)') ' Symmetry:', INTREP(I)
  300    CONTINUE
      END IF
C
C     ***** Number of basis functions *****
C
      NBAST = 0
      IF (ADDSTO .AND. INTTYP .EQ. 1) THEN
         KMAXT = NSMLSH + NLRGSH
      ELSE
         KMAXT = KMAX
      END IF
      DO 400 ISHELL = 1, KMAXT
         DO 410 KB = 0,MAXREP
            IF (IBTAND(KB,ISTBAO(ISHELL)) .EQ. 0)
     *            NBAST = NBAST + KHKT(ISHELL)
  410    CONTINUE
  400 CONTINUE
      NELMNT = NBAST*(NBAST + 1)/2
      IF (SQUARE) NELMNT = NBAST*NBAST
C
C     *******************************
C     ***** Calculate integrals *****
C     *******************************
C
      KSOINT = KLAST
      KWRK   = KSOINT + NOPTYP*NELMNT
      LWRK   = LWORK - KWRK + 1
      KLAST  = KWRK
      IF (KWRK .GT. LWORK) CALL STOPIT('PR1IN1','NE12',KWRK,LWORK)
      CALL PR1DRV(WORK(KSOINT),NELMNT,WORK(KWRK),LWRK,
     &            NPQUAD,LABINT,INTTYP,INTREP,NOPTYP,NBAST,ANTI,
     &            IORDER,WORK(KATOM),INTADR,TRIANG,NATOM,SQUARE,
     &            IPRINT,DOINT)
C
C     *********************************************************
C     ***** Transform from Cartesian to spherical moments *****
C     *********************************************************
C
      IF (TRASPH) THEN
         CALL SPHTRA(WORK(KSOINT),WORK(KWRK),LWRK,IORDER,NELMNT,
     &               NBAST,NOPTYP,IPRINT)
         CALL SMTYP(IORDER,NOPTYP,INTREP,LABINT)
      END IF
C
C     *******************************************************
C     ***** Transform from Cartesian to Spherical EFG's *****
C     *******************************************************
C
      IF (INTTYP .EQ. 31) THEN
         CALL EFGSPH(WORK(KSOINT),WORK(KWRK),LWRK,NELMNT,NBAST,
     &               NOPTYP,WORK(KATOM),NATOM,IPRINT)
         CALL FGSTYP(NOPTYP,INTREP,LABINT,WORK(KATOM),NATOM)
      END IF
C
C     *************************************************
C     ***** Test diamagnetic spin-orbit integrals *****
C     *************************************************
C
      IF (INTTYP .EQ. 12 .AND. .NOT.TRIANG) THEN
         KDIFF = KLAST
         KLAST = KDIFF + NELMNT
         IF (KLAST .GT. LWORK) CALL STOPIT('PR1IN1','TST',KLAST,LWORK)
         CALL DSOTST(WORK(KSOINT),WORK(KDIFF),NBAST,NELMNT,NOPTYP,
     &               LABINT,WORK(KATOM),NPQUAD,INTADR,IPRINT)
      END IF
C
C     ***********************************
C     ***** Write integrals on file *****
C     ***********************************
C
      IF (TOFILE) THEN
         IF (SOLVNT) THEN
            CALL  WRTSOL(WORK(KSOINT),IORDER,NBAST,NELMNT,NOPTYP,
     &                   INTREP,.NOT.ALLRLM,IPRINT)
C        ... for non-symmtric response (in ABACUS or RESPONS) and
C            Direct Reaction Field ALLRLM must be .true.
         ELSE
            CALL GETDAT(RTNLBL(1),RTNLBL(2))
C           Replace time information with symmetry information
            IF (SQUARE) THEN
               RTNLBL(2)='SQUARE  '
            ELSE
               IF (ANTI) THEN
                  RTNLBL(2)='ANTISYMM'
               ELSE
                  RTNLBL(2)='SYMMETRI'
               END IF
            END IF
C
            IADR = KSOINT
            DO 600 I = 1, NOPTYP
               IF (PROPRI .OR. IPRINT .GT. 4) THEN
                  CALL AROUND('Integrals of operator: '//LABINT(I))
                  WRITE (LUPRI,'(A,2X,A3,A1,I2,A1)')
     &               ' Symmetry of operator:',
     &               REP(INTREP(I)),'(',(INTREP(I) + 1),')'
                  IF (SQUARE) THEN
                     CALL OUTPUT(WORK(IADR),1,NBAST,1,NBAST,NBAST,NBAST,
     &                           1,LUPRI)
                  ELSE
                     CALL OUTPAK(WORK(IADR),NBAST,1,LUPRI)
                  END IF
               END IF
               CALL WRTPRO(WORK(IADR),NELMNT,LABINT(I),RTNLBL)
               IADR = IADR + NELMNT
  600       CONTINUE
C
         END IF
      END IF
C
C     If integrals requested from ABACUS, copy matrices from SOINT to
C     SINTMA for further use in ABACUS. KR, Spring 1992
C
      IF (NCOMP .NE. 0) THEN
         NCOMP = NOPTYP
         IF ((INTTYP .EQ. 15).OR.(INTTYP .EQ. 19).OR.(INTTYP .EQ. 9 )
     &                       .OR.(INTTYP .EQ. 10).OR.(INTTYP .EQ. 11)
     &                       .OR.(INTTYP .EQ. 13).OR.(INTTYP .EQ. 17)
     &                       .OR.(INTTYP .EQ. 18).OR.(INTTYP .EQ. 42)
     &                       .OR.(INTTYP .EQ. 20)) THEN
            IADR = KSOINT
            DO 1010 ICOMP = 1, NCOMP
               IF (TRIANG) THEN
                  IF (MTFORM .EQ. 'TRIANG') THEN
                     CALL DCOPY(NELMNT,WORK(IADR),1,
     &                    SINTMA((ICOMP - 1)*NELMNT + 1),1)
                  ELSE
                     CALL DAPTGE(NBAST,WORK(IADR),
     &                    SINTMA(NBAST*NBAST*(ICOMP - 1) + 1))
                  END IF
               ELSE
                  CALL DCOPY(NELMNT,WORK(IADR),1,
     &                 SINTMA((ICOMP - 1)*NELMNT + 1),1)
               END IF
               IADR = IADR + NELMNT
 1010       CONTINUE
         ELSE
            IADR = KSOINT
            DO 1020 ICOMP = 1, NCOMP
               IF (TRIANG) THEN
                  IF (MTFORM .EQ. 'TRIANG') THEN
                     CALL DCOPY(NELMNT,WORK(IADR),1,
     &                    SINTMA((ICOMP - 1)*NELMNT + 1),1)
                  ELSE
                     CALL DSPTSI(NBAST,WORK(IADR),
     &                    SINTMA(NBAST*NBAST*(ICOMP - 1) + 1))
                  END IF
               ELSE
                  CALL DCOPY(NELMNT,WORK(IADR),1,
     &                 SINTMA((ICOMP - 1)*NELMNT + 1),1)
               END IF
               IADR = IADR + NELMNT
 1020       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck ioden */
      FUNCTION IODEN(I)
#include <implicit.h>
      DIMENSION IOD(14)
      DATA IOD /1,2,3,4,11,12,13,14,5,6,7,8,9,10/
      IODEN = IOD(I)
      RETURN
      END
C  /* Deck setatm */
      SUBROUTINE SETATM(DOATOM,NATOM,INTTYP)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
C
      LOGICAL DOATOM(NUCIND), SPIN
#include <nuclei.h>
#include <symmet.h>
#include <cbiher.h>
C
      SPIN = INTTYP .EQ.  9 .OR. INTTYP .EQ. 10 .OR.
     &       INTTYP .EQ. 11 .OR. INTTYP .EQ. 12 .OR. INTTYP .EQ. 13 .OR.
     &       INTTYP .EQ. 26 .OR. INTTYP .EQ. 27 .OR. INTTYP .EQ. 28 .OR.
     &       INTTYP .EQ. 34 .OR. INTTYP .EQ. 38 .OR. INTTYP .EQ. 48
      IF (.NOT.SPIN .OR. ALLATM) THEN
         NATOM = NUCDEP
         DO 100 I = 1, NUCIND
            DOATOM(I) = .TRUE.
  100    CONTINUE
      ELSE
         DO 200 I = 1, NUCIND
            DOATOM(I) = .FALSE.
  200    CONTINUE
         NATOM = 0
         DO 300 I = 1, NPATOM
            IATOM = IPATOM(I)
            IF (IATOM .GT. NUCIND) THEN
               WRITE (LUPRI,'(A,I3,A)') ' Error in input: atom ',
     &         IATOM,' does not exist.'
               CALL QUIT('Error in SETATM')
            END IF
            NATOM = NATOM + MULT(ISTBNU(IATOM))
            DOATOM(IATOM) = .TRUE.
  300    CONTINUE
      END IF
      RETURN
      END
C  /* Deck ntypts */
      SUBROUTINE NTYPTS(NOPTYP)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      IF (NOPTYP .GT. 9*MXCENT) THEN
         WRITE (LUPRI,'(/A,I5)') ' NOPTYP too large in NTYPTS '
         WRITE (LUPRI,'(/A,I5)') ' MAXTYP: ', 9*MXCENT
         WRITE (LUPRI,'(/A,I5)') ' NOPTYP: ', NOPTYP
         CALL QUIT('ERROR in NTYPTS, NOPTYP out of range')
      END IF
      RETURN
      END
C  /* Deck cmtyp */
      SUBROUTINE CMTYP(IORDER,NOPTYP,INTREP,LABINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP((IORDER + 1)*(IORDER + 2)/2)
      DIMENSION IX(9*MXCENT), IY(9*MXCENT), IZ(9*MXCENT)
      CHARACTER LABINT((IORDER + 1)*(IORDER + 2)/2)*8
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
      NOPTYP = (IORDER + 1)*(IORDER + 2)/2
      CALL NTYPTS(NOPTYP)
      CALL LMNVAL(IORDER+1,NOPTYP,IX,IY,IZ)
      DO 100 I = 1, NOPTYP
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
C
C        Label
C
         LABINT(I) = 'CM'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                   //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                   //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
C
C        Symmetry
C
         ISX = MOD(NX,2)*ISYMAX(1,1)
         ISY = MOD(NY,2)*ISYMAX(2,1)
         ISZ = MOD(NZ,2)*ISYMAX(3,1)
         INTREP(I) = IBTXOR(ISX,IBTXOR(ISY,ISZ))
  100 CONTINUE
      RETURN
      END
C  /* Deck sl1typ */
      SUBROUTINE SL1TYP(IORDER,NOPTYP,INTREP,LABINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP((IORDER + 1)*(IORDER + 2)*3/2)
      DIMENSION IX(9*MXCENT), IY(9*MXCENT), IZ(9*MXCENT)
      CHARACTER LABINT((IORDER + 1)*(IORDER + 2)*3/2)*8
#include <symmet.h>
#include <mgsolt.h>
#include <chrnos.h>
#include <ibtfun.h>
      IF (DOALL) THEN
         NOPTYP = (IORDER + 1)*(IORDER + 2)*3/2
      ELSE
         NOPTYP = (IORDER + 1)*(IORDER + 2)/2
      END IF
      CALL NTYPTS(NOPTYP)
      CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,IX,IY,IZ)
      DO 100 I = 1, (IORDER + 1)*(IORDER + 2)/2
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
C
C        Label
C
         IF (DOALL) THEN
         LABINT(3*(I - 1) + 1) = 'X '//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
         LABINT(3*(I - 1) + 2) = 'Y '//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
         LABINT(3*(I - 1) + 3) = 'Z '//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
C
C        Symmetry
C
         ISX = MOD(NX,2)*ISYMAX(1,1)
         ISY = MOD(NY,2)*ISYMAX(2,1)
         ISZ = MOD(NZ,2)*ISYMAX(3,1)
         INTREP(3*(I - 1) + 1) = IBTXOR(IBTXOR(ISX,
     &                           IBTXOR(ISY,ISZ)),ISYMAX(1,2))
         INTREP(3*(I - 1) + 2) = IBTXOR(IBTXOR(ISX,
     &                           IBTXOR(ISY,ISZ)),ISYMAX(2,2))
         INTREP(3*(I - 1) + 3) = IBTXOR(IBTXOR(ISX,
     &                           IBTXOR(ISY,ISZ)),ISYMAX(3,2))
         ELSE
            ISX = MOD(NX,2)*ISYMAX(1,1)
            ISY = MOD(NY,2)*ISYMAX(2,1)
            ISZ = MOD(NZ,2)*ISYMAX(3,1)
            IF (ICOMP .EQ. 1) THEN
               LABINT(I) = 'X '//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = IBTXOR(IBTXOR(ISX,
     &                           IBTXOR(ISY,ISZ)),ISYMAX(1,2))
            ELSE IF (ICOMP .EQ. 2) THEN
               LABINT(I) = 'Y '//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = IBTXOR(IBTXOR(ISX,
     &                           IBTXOR(ISY,ISZ)),ISYMAX(2,2))
            ELSE
               LABINT(I) = 'Z '//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = IBTXOR(IBTXOR(ISX,
     &                           IBTXOR(ISY,ISZ)),ISYMAX(3,2))
            END IF
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck sl2typ */
      SUBROUTINE SL2TYP(IORDER,NOPTYP,INTREP,LABINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP((IORDER + 1)*(IORDER + 2)*6/2)
      DIMENSION IX(9*MXCENT), IY(9*MXCENT), IZ(9*MXCENT)
      CHARACTER LABINT((IORDER + 1)*(IORDER + 2)*6/2)*8
#include <symmet.h>
#include <mgsolt.h>
#include <chrnos.h>
#include <ibtfun.h>
      IF (DOALL) THEN
         NOPTYP = (IORDER + 1)*(IORDER + 2)*6/2
      ELSE
         NOPTYP = (IORDER + 1)*(IORDER + 2)/2
      END IF
      CALL NTYPTS(NOPTYP)
      CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,IX,IY,IZ)
      DO 100 I = 1, (IORDER + 1)*(IORDER + 2)/2
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
C
C        Label
C
         IF (DOALL) THEN
         LABINT(6*(I - 1) + 1) = 'XX'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
         LABINT(6*(I - 1) + 2) = 'XY'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
         LABINT(6*(I - 1) + 3) = 'XZ'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
         LABINT(6*(I - 1) + 4) = 'YY'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
         LABINT(6*(I - 1) + 5) = 'YZ'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
         LABINT(6*(I - 1) + 6) = 'ZZ'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                               //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                               //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
C
C        Symmetry
C
         ISX = MOD(NX,2)*ISYMAX(1,1)
         ISY = MOD(NY,2)*ISYMAX(2,1)
         ISZ = MOD(NZ,2)*ISYMAX(3,1)
         ISYMC = IBTXOR(ISX,IBTXOR(ISY,ISZ))
         INTREP(6*(I - 1) + 1) = ISYMC
         INTREP(6*(I - 1) + 2) = IBTXOR(ISYMC,
     &                           IBTXOR(ISYMAX(1,2),ISYMAX(2,2)))
         INTREP(6*(I - 1) + 3) = IBTXOR(ISYMC,
     &                           IBTXOR(ISYMAX(1,2),ISYMAX(3,2)))
         INTREP(6*(I - 1) + 4) = ISYMC
         INTREP(6*(I - 1) + 5) = IBTXOR(ISYMC,
     &                           IBTXOR(ISYMAX(2,2),ISYMAX(3,2)))
         INTREP(6*(I - 1) + 6) = ISYMC
         ELSE
            ISX = MOD(NX,2)*ISYMAX(1,1)
            ISY = MOD(NY,2)*ISYMAX(2,1)
            ISZ = MOD(NZ,2)*ISYMAX(3,1)
            ISYMC = IBTXOR(ISX,IBTXOR(ISY,ISZ))
            IF (ICOMP .EQ. 1) THEN
               LABINT(I) = 'XX'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &              //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &              //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = ISYMC
            ELSE IF (ICOMP .EQ. 2) THEN
               LABINT(I) = 'XY'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &              //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &              //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = IBTXOR(ISYMC,IBTXOR(ISYMAX(1,2),ISYMAX(2,2)))
            ELSE IF (ICOMP .EQ. 3) THEN
               LABINT(I) = 'XZ'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &              //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &              //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = IBTXOR(ISYMC,IBTXOR(ISYMAX(1,2),ISYMAX(3,2)))
            ELSE IF (ICOMP .EQ. 4) THEN
               LABINT(I) = 'YY'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &              //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &              //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = ISYMC
            ELSE IF (ICOMP .EQ. 5) THEN
               LABINT(I) = 'YZ'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &              //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &              //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = IBTXOR(ISYMC,IBTXOR(ISYMAX(2,2),ISYMAX(3,2)))
            ELSE 
               LABINT(I) = 'ZZ'//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &              //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &              //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               INTREP(I) = ISYMC
            END IF
C
C        Symmetry
C
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck smtyp */
      SUBROUTINE SMTYP(IORDER,NOPTYP,INTREP,LABINT)
C
C     Calculates labels and symmetries for Spherical Moments
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(2*IORDER + 1)
      DIMENSION IX(9*MXCENT), IY(9*MXCENT), IZ(9*MXCENT)
      CHARACTER LABINT(2*IORDER + 1)*8
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
      NOPTYP = 2*IORDER + 1
      CALL NTYPTS(NOPTYP)
C
      LABINT(1) = 'SM'//CHRNOS(IORDER/10)//CHRNOS(MOD(IORDER,10))//'+00'
      INTREP(1) = IREPLM(IORDER,0)
      DO 100 I = 1, 2*IORDER
         M = (I + 1)/2
         IF (MOD(I,2) .EQ. 1) THEN
            LABINT(I+1) ='SM'//CHRNOS(IORDER/10)//CHRNOS(MOD(IORDER,10))
     &                 //'+'//CHRNOS(M/10)//CHRNOS(MOD(M,10))
            INTREP(I+1) = IREPLM(IORDER,M)
         ELSE
            LABINT(I+1) ='SM'//CHRNOS(IORDER/10)//CHRNOS(MOD(IORDER,10))
     &                 //'-'//CHRNOS(M/10)//CHRNOS(MOD(M,10))
            INTREP(I+1) = IREPLM(IORDER,-M)
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck frmtyp */
      SUBROUTINE FRMTYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
      NOPTYP = NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
                  ITYP = ITYP + 1
                  LABINT(ITYP) = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                                //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &	                              CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
                  INTREP(ITYP) = IREP
               END IF
            END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck nsttyp */
      SUBROUTINE NSTTYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM,LABEL,INTADR)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8, LABEL*4
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
      NOPTYP = 0
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 300 LCOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + LCOOR, IREP,2)
                  IF (ISCOOR .GT. 0) THEN
                     DO 400 ICOOR = 1, 3
                        ITYP = ITYP + 1
                        LABINT(ITYP) = CHRNOS(ISCOOR/10)//
     &                                 CHRNOS(MOD(ISCOOR,10))//LABEL//
     &                                 ' '//CHRXYZ(ICOOR)
                        INTREP(ITYP) = IBTXOR(IREP,ISYMAX(ICOOR,2))
                        LSCOOR = 3*(ISCOOR - 1) + ICOOR
                        INTADR(LSCOOR) = ITYP
                        NOPTYP = NOPTYP + 1
 400                 CONTINUE
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
      CALL NTYPTS(NOPTYP)
      RETURN
      END
C  /* Deck psotyp */
      SUBROUTINE PSOTYP(NOPTYP,INTREP,LABINT,DOATOM,INTADR,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
      NOPTYP = 3*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            DO 300 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,2)
               IF (ISCOOR .GT. 0) THEN
                  ITYP = ITYP + 1
                  LABINT(ITYP) = 'PSO '//CHRNOS(ISCOOR/10)
     &                                 //CHRNOS(MOD(ISCOOR,10))//'  '
                  INTREP(ITYP) = IREP
                  INTADR(ISCOOR) = ITYP
               END IF
  300       CONTINUE
         END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck sdtyp */
      SUBROUTINE SDTYP(INTTYP,NOPTYP,INTREP,LABINT,INTADR,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
      NOPTYP = 9*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            DO 300 ICOOR1 = 1, 3
               ISCOR1 = IPTCNT(3*(IATOM - 1) + ICOOR1,IREP,2)
               IF (ISCOR1 .GT. 0) THEN
                  DO 400 ICOOR2 = 1, 3
                     ITYP = ITYP + 1
                     IF (INTTYP .EQ. 11) THEN
                        LABINT(ITYP) = 'SD  '//CHRNOS(ISCOR1/10)
     &                                       //CHRNOS(MOD(ISCOR1,10))
     &                                       //' '//CHRXYZ(-ICOOR2)
                     ELSE
                        LABINT(ITYP) = 'SDC '//CHRNOS(ISCOR1/10)
     &                                       //CHRNOS(MOD(ISCOR1,10))
     &                                       //' '//CHRXYZ(-ICOOR2)
                     END IF
                     INTREP(ITYP) = IBTXOR(IREP,ISYMAX(ICOOR2,2))
                     ISCOOR = 3*(ISCOR1 - 1) + ICOOR2
                     INTADR(ISCOOR) = ITYP
  400             CONTINUE
               END IF
  300       CONTINUE
         END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck dsotyp */
      SUBROUTINE DSOTYP(NOPTYP,INTREP,LABINT,INTADR,DOATOM,NATOM,TRIANG)
C
C     helgaker
C
C     Note: Integrals are obtained by numerical integration
C           (Gaussian quadrature) of field integrals according to
C           O. Matsuoka and T. Aoyama, JCP 73 , 5718 (1980).
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      PARAMETER (D1 = 1.0D0)
      DIMENSION INTREP(*), INTADR(*)
      LOGICAL DOATOM(NUCIND), TRIANG, SAME
      CHARACTER LABINT(*)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
C
      IF (TRIANG) THEN
         NOPTYP = (3*NATOM*(3*NATOM + 1)/2)
      ELSE
         NOPTYP = 9*NATOM*NATOM
      END IF
C
C     no check since the allocations of INTADR and INTREP are exact!
C
C     CALL NTYPTS(NOPTYP)
C
      ITYP = 0
C
C     Irreps
C
      DO 100 IREPO = 0, MAXREP
C
C        Atoms
C
         DO 200 IATOM1 = 1, NUCIND
         IF (DOATOM(IATOM1)) THEN
            MXATM2 = NUCIND
            IF (TRIANG) MXATM2 = IATOM1
            DO 400 IATOM2 = 1, MXATM2
            IF (DOATOM(IATOM2)) THEN
               SAME = TRIANG .AND. IATOM1.EQ.IATOM2
C
C              Cartesian directions
C
               DO 500 ICOOR1 = 1, 3
                  ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
                  IF (ISCOR1 .GT. 0) THEN
                     MXCR2 = 3
                     IF (SAME) MXCR2 = ICOOR1
                     DO 600 ICOOR2 = 1, MXCR2
                        ISCOR2 = IPTCNT(3*(IATOM2-1)+ICOOR2,IREPO,2)
                        IF (ISCOR2 .GT. 0) THEN
                           ITYP = ITYP + 1
                           IF (TRIANG) THEN
                              MXCOR = MAX(ISCOR1,ISCOR2)
                              MNCOR = MIN(ISCOR1,ISCOR2)
                              ISCOOR = MXCOR*(MXCOR - 1)/2 + MNCOR
                           ELSE
                              ISCOOR = 3*NUCDEP*(ISCOR1 - 1) + ISCOR2
                           END IF
                           LABINT(ITYP) = 'DSO '
     &                       //CHRNOS(ISCOR1/10)//CHRNOS(MOD(ISCOR1,10))
     &                       //CHRNOS(ISCOR2/10)//CHRNOS(MOD(ISCOR2,10))
                           INTREP(ITYP)   = 0
                           INTADR(ISCOOR) = ITYP
                        END IF
  600                CONTINUE
                  END IF
  500          CONTINUE
            END IF
  400       CONTINUE
         END IF
  200    CONTINUE
  100 CONTINUE
      NOPTYP = ITYP
      CALL NTYPTS(NOPTYP)
      RETURN
      END
C  /* Deck hdbtyp */
      SUBROUTINE HDBTYP(NOPTYP,INTREP,INTADR,LABINT,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
      NOPTYP = 9*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 300 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     DO 400 LCOOR = 1, 3
                        ITYP = ITYP + 1
                        LABINT(ITYP) = CHRNOS(ISCOOR/10)//
     &                                 CHRNOS(MOD(ISCOOR,10))//' HDB '//
     &                                 CHRXYZ(LCOOR)
                        INTREP(ITYP) = IBTXOR(IREP,ISYMAX(LCOOR,2))
                        LSCOOR = 3*(ISCOOR - 1) + LCOOR
                        INTADR(LSCOOR) = ITYP
 400                 CONTINUE
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck dpgtyp */
      SUBROUTINE DPGTYP(NOPTYP,INTREP,INTADR,LABINT,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
      NOPTYP = 9*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 300 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     IXY = 0
                     DO 400 LCOOR = 1, 3
                        ITYP = ITYP + 1
                        IXY = IXY + 1
                        LABINT(ITYP) = CHRNOS(ISCOOR/10)//
     &                                 CHRNOS(MOD(ISCOOR,10))//' DPG '//
     &                                 CHRXYZ(LCOOR)
CCH                     wrong copy & paste from QUGTYP ??
CCH                     IAX = IBTXOR(ISYMAX(LCOOR,1),ISYMAX(MCOOR,1))
CCH
                        IAX = ISYMAX(LCOOR,1)
CCH
                        INTREP(ITYP) = IBTXOR(IREP,IAX)
                        LSCOOR = 3*(ISCOOR - 1) + IXY
                        INTADR(LSCOOR) = ITYP
 400                 CONTINUE
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck qugtyp */
      SUBROUTINE QUGTYP(NOPTYP,INTREP,INTADR,LABINT,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(18*MXCENT), INTADR(18*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(18*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
      NOPTYP = 18*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 300 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     IXY = 0
                     DO 400 LCOOR = 1, 3
                     DO 400 MCOOR = LCOOR, 3
                        ITYP = ITYP + 1
                        IXY = IXY + 1
                        LABINT(ITYP) = CHRNOS(ISCOOR/10)//
     &                                 CHRNOS(MOD(ISCOOR,10))//'QDG '//
     &                                 CHRXYZ(LCOOR)//CHRXYZ(MCOOR)
                        IAX = IBTXOR(ISYMAX(LCOOR,1),ISYMAX(MCOOR,1))
                        INTREP(ITYP) = IBTXOR(IREP,IAX)
                        LSCOOR = 6*(ISCOOR - 1) + IXY
                        INTADR(LSCOOR) = ITYP
 400                 CONTINUE
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck ocgtyp */
      SUBROUTINE OCGTYP(NOPTYP,INTREP,INTADR,LABINT,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(30*MXCENT), INTADR(30*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(30*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
      NOPTYP = 30*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 300 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     IXY = 0
                     DO 400 LCOOR = 1, 3
                     DO 400 MCOOR = LCOOR, 3
                     DO 400 NCOOR = MCOOR, 3
                        ITYP = ITYP + 1
                        IXY = IXY + 1
                        LABINT(ITYP) = CHRNOS(ISCOOR/10)//
     &                                 CHRNOS(MOD(ISCOOR,10))//'ODG'//
     &                                 CHRXYZ(LCOOR)//CHRXYZ(MCOOR)//
     &                                 CHRXYZ(NCOOR)
CCH
CCH                     isn't this a wrong copy & paste from QUCTYP ?
CCH
                        IAX = IBTXOR(ISYMAX(LCOOR,1),ISYMAX(MCOOR,1))
CCH
CCH                     guess here it should be introduced an additional
CCH                     line with:
                        IAX = IBTXOR(IAX,ISYMAX(NCOOR,1))
CCH
                        INTREP(ITYP) = IBTXOR(IREP,IAX)
                        LSCOOR = 10*(ISCOOR - 1) + IXY
                        INTADR(LSCOOR) = ITYP
 400                 CONTINUE
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck hdotyp */
      SUBROUTINE HDOTYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM,ANTISY,INTTYP)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(*)
      LOGICAL DOATOM(NUCIND), ANTISY
      CHARACTER LABINT(*)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
      NOPTYP = 3*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 300 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     ITYP = ITYP + 1
                     IF (ANTISY) THEN
                        LABINT(ITYP) = 'HDO '//CHRNOS(ISCOOR/10)
     &                                    //CHRNOS(MOD(ISCOOR,10))//'  '
                     ELSE
                        IF (INTTYP .EQ. 14) THEN
                           LABINT(ITYP) = 'SQHDO'//CHRNOS(ISCOOR/100)
     &                              //CHRNOS(MOD(ISCOOR,100)/10)
     &                              //CHRNOS(MOD((MOD(ISCOOR,100)),10))
                        ELSE IF (INTTYP .EQ. 59) THEN
                           LABINT(ITYP) = '1DOVL'//CHRNOS(ISCOOR/100)
     &                              //CHRNOS(MOD(ISCOOR,100)/10)
     &                              //CHRNOS(MOD((MOD(ISCOOR,100)),10))
                        ELSE IF (INTTYP .EQ. 70) THEN
                           LABINT(ITYP) = '1DHAM'//CHRNOS(ISCOOR/100)
     &                              //CHRNOS(MOD(ISCOOR,100)/10)   
     &                              //CHRNOS(MOD((MOD(ISCOOR,100)),10))
                        ELSE
                           LABINT(ITYP) = 'SQHDR'//CHRNOS(ISCOOR/100)
     &                              //CHRNOS(MOD(ISCOOR,100)/10)
     &                              //CHRNOS(MOD((MOD(ISCOOR,100)),10))
                        END IF
                     END IF
                     INTREP(ITYP) = IREP
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck npetyp */
      SUBROUTINE NPETYP(NOPTYP,INTREP,NPETBL,LABINT,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), NPETBL(NUCIND,0:MAXREP)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8, NAME*2
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
C
      NOPTYP = NATOM
      CALL NTYPTS(NOPTYP)
      CALL DZERO(NPETBL,NUCIND*(MAXREP + 1))
C
      ITYP = 0
      DO 100 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            NAME = NAMEX(3*(IATOM - 1) + 1)(1:2)
            DO 200 IREP = 0, MAXREP
               IF (IBTAND(ISTBNU(IATOM),IREP) .EQ. 0) THEN
                  ITYP = ITYP + 1
                  NPETBL(IATOM,IREP) = ITYP
                  LABINT(ITYP) = 'POT.E '//NAME
                  INTREP(ITYP) = IREP
               END IF
 200        CONTINUE
         END IF
 100  CONTINUE
      RETURN
      END
C  /* Deck efntyp */
      SUBROUTINE EFNTYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM,INTADR,INTTYP)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), INTADR(*)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8, NAME*2
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
      NOPTYP = 3*NATOM
      CALL NTYPTS(NOPTYP)
C
      ITYP = 0
      DO 50 IREP = 0, MAXREP
         DO 100 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 200 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     ITYP = ITYP + 1
                     IF (INTTYP .EQ. 70) THEN
                        LABINT(ITYP) = '1DHAM'//CHRNOS(ISCOOR/100)
     &                           //CHRNOS(MOD(ISCOOR,100)/10)   
     &                           //CHRNOS(MOD((MOD(ISCOOR,100)),10))
                     ELSE
                        LABINT(ITYP) = 'NEF '//CHRNOS(ISCOOR/10)
     &                                       //CHRNOS(MOD(ISCOOR,10))
     &                                       //'  '
                     END IF                                                   
                     INTREP(ITYP) = IREP
                     INTADR(ISCOOR) = ITYP
                  END IF
 200           CONTINUE
            END IF
 100     CONTINUE
 50   CONTINUE 
      RETURN
      END
C  /* Deck efgtyp */
      SUBROUTINE EFGTYP(NOPTYP,INTREP,IFGTBL,LABINT,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(6*MXCENT), IFGTBL(NUCIND,6,0:MAXREP)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(6*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
C
      NOPTYP = 6*NATOM
      CALL NTYPTS(NOPTYP)
      CALL IZERO(IFGTBL,6*NUCIND*(MAXREP + 1))
C
      ITYP = 0
      DO 100 IATOM = 1, NUCIND
      IF (DOATOM(IATOM)) THEN
         IJ = 0
         DO 200 ICOOR1 = 1, 3
         DO 200 ICOOR2 = ICOOR1, 3
            IJ = IJ + 1
            ISYMIJ = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IBTAND(ISTBNU(IATOM),IBTXOR(IREPC,ISYMIJ)).EQ.0) THEN
                  IOFF = IOFF + 1
                  ITYP = ITYP + 1
                  IFGTBL(IATOM,IJ,IREPC) = ITYP
                  LABINT(ITYP) = CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'//
     &               CHRNOS(IATOM/10)//CHRNOS(MOD(IATOM,10))//
     &               CHRNOS(IOFF)
                  INTREP(ITYP) = IREPC
               END IF
 300        CONTINUE
 200     CONTINUE
      END IF
 100  CONTINUE
      RETURN
      END
C  /* Deck efgsph */
      SUBROUTINE EFGSPH(AINT,WORK,LWORK,NELMNT,NBAST,NOPTYP,DOATOM,
     &                  NATOM,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxmom.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
      LOGICAL DOATOM(9*MXCENT)
      DIMENSION AINT(NELMNT,NOPTYP), WORK(LWORK)
#include <nuclei.h>
#include <symmet.h>
C
      IORDER = 2
      NLM    = 5
      NCOMPT = 6
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from EFGSPH','*',103)
         WRITE (LUPRI,'(2X,A,I5)') ' NBAST:  ', NBAST
         WRITE (LUPRI,'(2X,A,I5)') ' NELMNT: ', NELMNT
         WRITE (LUPRI,'(2X,A,I5)') ' NOPTYP: ', NOPTYP
      END IF
      KTRA = 1
      KINT = KTRA + NLM*NOPTYP
      KWRK = KINT + NELMNT*NLM*NATOM
      LWRK = LWORK - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('EFGSPH',' ',KWRK,LWORK)
      CALL EFGSP1(AINT,WORK(KINT),IORDER,WORK(KTRA),WORK(KWRK),LWRK,
     &            NCOMPT,NLM,NELMNT,NATOM,DOATOM,NOPTYP,IPRINT)
      RETURN
      END
C  /* Deck efgsp1 */
      SUBROUTINE EFGSP1(CARINT,SPHINT,IORDER,TRAMAT,WORK,LWORK,NXYZ,
     &                  NLM,NELMNT,NATOM,DOATOM,NOPTYP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxmom.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
      PARAMETER (D0 = 0.D0)
      LOGICAL DOATOM(NUCIND)
      DIMENSION CARINT(NELMNT,NOPTYP), TRAMAT(NXYZ,NLM), JTABLE(6),
     &          SPHINT(NELMNT,NLM*NATOM), ITABLE(6), WORK(LWORK)
#include <nuclei.h>
#include <symmet.h>
      DATA ITABLE /1,2,3,5,6,9/
#include <ibtfun.h>
      CALL SPHCOM(IORDER,TRAMAT,NLM,NXYZ,0,0,WORK,LWORK,IPRINT)
      CALL DZERO(SPHINT,NELMNT*NLM*NATOM)
      ITYP  = 0
      ICOMP = 0
      DO 150 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            NUATOM = 0
            DO 101 IREPAX = 0, MAXREP
               IF (IBTAND(IREPAX,ISTBNU(IATOM)).EQ.0) NUATOM = NUATOM+1
 101        CONTINUE
            DO 102 I = 1, 6
               JTABLE(I) = (I - 1)*NUATOM + 1
 102        CONTINUE
            DO 100 I = 1, NLM
               JATOM = 0
               DO 160 IREP = 0, MAXREP
                  IF (IBTAND(IREP,ISTBNU(IATOM)) .EQ. 0) THEN
                     ITYP   = ITYP + 1
                     ICOMP2 = ICOMP + 6*NUATOM*(I/(5*NUATOM)) + JATOM
                     DO 110 J = 1, NXYZ
                        COEF = TRAMAT(J,I)
                        IF (ABS(COEF) .GT. D0) THEN
                           IF (IPRINT .GT. 5)
     &                          WRITE(LUPRI,'(1X,A,2I5,F12.6)')
     &                          ' I, J, COEF ', I, J, COEF
                           DO 300 K = 1, NELMNT
                              SPHINT(K,ITYP) = SPHINT(K,ITYP)
     &                             + COEF*CARINT(K,ICOMP2+JTABLE(J))
 300                       CONTINUE
                        END IF
 110                 CONTINUE
                     JATOM = JATOM + 1
                  END IF
 160           CONTINUE
 100        CONTINUE
            ICOMP = ICOMP + NUATOM*6
         END IF
 150  CONTINUE
      CALL DCOPY(NLM*NATOM*NELMNT,SPHINT,1,CARINT,1)
      RETURN
      END
C  /* Deck fgstyp */
      SUBROUTINE FGSTYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM)
C
C     K.Ruud, June 1991
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      LOGICAL DOATOM(NUCIND)
      DIMENSION INTREP(5*NATOM)
      CHARACTER LABINT(5*NATOM)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ibtfun.h>
C
      NOPTYP = NATOM*5
      CALL NTYPTS(NOPTYP)
C
      ITYP = 0
      DO 10 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            DO 20 I = 0, 4
               IF (I .EQ. 0) THEN
                  ISYMCO = IREPLM(2,0)
               ELSE
                  M = (I + 1)/2
                  IF (MOD(I,2) .EQ. 1) THEN
                     ISYMCO = IREPLM(2,M)
                  ELSE
                     ISYMCO = IREPLM(2,-M)
                  END IF
               END IF
               DO 30 IREP = 0, MAXREP
               IF (IBTAND(IBTXOR(IREP,ISYMCO),ISTBNU(IATOM)).EQ.0) THEN
                  ITYP = ITYP + 1
                  IF (I .EQ. 0) THEN
                     LABINT(ITYP) = '00'//'EFG '//
     &                  CHRNOS(IATOM/10)//CHRNOS(MOD(IATOM,10))
                  ELSE
                     M = (I + 1)/2
                     IF (MOD(I,2) .EQ. 1) THEN
                        LABINT(ITYP) = '+'//CHRNOS(M)//'EFG '//
     &                     CHRNOS(IATOM/10)//CHRNOS(MOD(IATOM,10))
                     ELSE
                        LABINT(ITYP) = '-'//CHRNOS(M)//'EFG '//
     &                     CHRNOS(IATOM/10)//CHRNOS(MOD(IATOM,10))
                     END IF
                  END IF
                  INTREP(ITYP) = IREP
               END IF
 30            CONTINUE
 20         CONTINUE
         END IF
 10   CONTINUE
      RETURN
      END
C  /* Deck gauher */
      SUBROUTINE GAUHER(X,W,N)
C
C  This routine returns
C  arrays X and W of length N, containing the abscissas and weights of the
C  Gauss-Hermite 2N-points quadrature formula
C
C
#include <implicit.h>
      DIMENSION X(N),W(N)
#include <pi.h>
C
C  Roots and abscissas from "Handbook of Mathematical Functions"
C  ed. M.Abramowitz and I.A.Stegun (Dover)
C
      DIMENSION W5(5),X5(5)
      DATA X5/.34290 13272 23705 D0,
     *       1.03661 08297 89514 D0,
     *       1.75668 36492 99882 D0,
     *       2.53273 16742 32790 D0,
     *       3.43615 91188 37738 D0/
      DATA W5/.68708 18539 513 D0,
     *        .70329 63231 049 D0,
     *        .74144 19319 436 D0,
     *        .82066 61264 048 D0,
     *       1.02545 16913 657 D0/
      DIMENSION W6(6),X6(6)
      DATA X6/.31424 03762 54359 D0,
     *        .94778 83912 40164 D0,
     *       1.59768 26351 52605 D0,
     *       2.27950 70805 01060 D0,
     *       3.02063 70251 20890 D0,
     *       3.88972 48978 69782 D0/
      DATA W6/.62930 78743 695 D0,
     *        .63962 12320 203 D0,
     *        .66266 27732 669 D0,
     *        .70522 03661 122 D0,
     *        .78664 39394 633 D0,
     *        .98969 90470 923 D0/
      DIMENSION W8(8),X8(8)
      DATA X8/.27348 10461 3815 D0,
     *        .82295 14491 4466 D0,
     *       1.38025 85391 9888 D0,
     *       1.95178 79909 1625 D0,
     *       2.54620 21578 4748 D0,
     *       3.17699 91619 7996 D0,
     *       3.86944 79048 6012 D0,
     *       4.68873 89393 0582 D0/
      DATA W8/.54737 52050 378 D0,
     *        .55244 19573 675 D0,
     *        .56321 78290 882 D0,
     *        .58124 72754 009 D0,
     *        .60973 69582 560 D0,
     *        .65575 56728 761 D0,
     *        .73824 56222 777 D0,
     *        .93687 44928 841 D0/
      DIMENSION W10(10),X10(10)
      DATA X10/.24534 07083 009 D0,
     *         .73747 37285 454 D0,
     *        1.23407 62153 953 D0,
     *        1.73853 77121 166 D0,
     *        2.25497 40020 893 D0,
     *        2.78880 60584 281 D0,
     *        3.34785 45673 832 D0,
     *        3.94476 40401 156 D0,
     *        4.60368 24495 507 D0,
     *        5.38748 08900 112 D0/
      DATA W10/.49092 15006 667 D0,
     *         .49384 33852 721 D0,
     *         .49992 08713 363 D0,
     *         .50967 90271 175 D0,
     *         .52408 03509 486 D0,
     *         .54485 17423 644 D0,
     *         .57526 24428 525 D0,
     *         .62227 86961 914 D0,
     *         .70433 29611 769 D0,
     *         .89859 19614 532 D0/
      IF (N.EQ.5) THEN
         DO 5 I=1,5
         X(I) = X5(I)
         W(I) = W5(I)
 5       CONTINUE
      ELSE IF (N.EQ.6) THEN
         DO 6 I=1,6
         X(I) = X6(I)
         W(I) = W6(I)
 6       CONTINUE
      ELSE IF (N.EQ.8) THEN
         DO 8 I=1,8
         X(I) = X8(I)
         W(I) = W8(I)
 8       CONTINUE
      ELSE IF (N.EQ.10) THEN
         DO 10 I=1,10
         X(I) = X10(I)
         W(I) = W10(I)
 10      CONTINUE
      ELSE
         PRINT '(A)', 'WRONG ORDER FOR GAUSS-HERMITE QUADRATURE'
         PRINT '(A)', 'YOUR VALUE: ',N
         PRINT '(A,4I4)', 'ALLOWED VALUES: ',5,6,8,10
      END IF
      RETURN
      END
C  /* Deck gauleg */
      SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C  Given lower and uper limits of integration X1 and X2, this routine
C  returns arrays X and W of length N, containing the abscissas and
C  weights of the Gauss-Legendre N-points quadrature formula
C
C  Written by G. Rybicki
C  Copied from "Numerical Recipes" by W.H. Press et.al.
C  by Olav Vahtras (910208)
C
#include <implicit.h>
#include <pi.h>
      PARAMETER (EPS=1.D-15)
      DIMENSION X(N),W(N)
C
C  Roots are symmetric in interval
C  Only necessary to find N+1/2 roots
C
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
C
C  Start guess of i:th zero
C
         Z=COS(PI*(I-.25D0)/(N+.5D0))
1        CONTINUE
            P1=1D0
            P2=0D0
            DO 11 J=1,N
               P3=P2
               P2=P1
               P1=((2D0*J-1D0)*Z*P2-(J-1D0)*P3)/J
11          CONTINUE
C
C  P1 is now the desired Legendre polynomial. We next compute PP, its
C  derivative, by a standard relation involving also P2, the polynomial of one
C  lower order
C
            PP=N*(Z*P1-P2)/(Z*Z-1D0)
            Z1=Z
C
C  Newton's Method
C
            Z=Z1-P1/PP
         IF (ABS(Z-Z1).GT.EPS) GO TO 1
C
C  Scale the root to the desired interval and put in its symmetric counterpart.
C
         X(I) = XM-XL*Z
         X(N+1-I) = XM+XL*Z
C
C  Compute the weight and its symmetric counterpart.
C

         W(I) = 2D0*XL/((1D0-Z*Z)*PP*PP)
         W(N+1-I) = W(I)
12    CONTINUE
      RETURN
      END
C  /* Deck sphtra */
      SUBROUTINE SPHTRA(AINT,WORK,LWORK,IORDER,NELMNT,NBAST,NOPTYP,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
      DIMENSION AINT(NELMNT,NOPTYP), WORK(LWORK)
C
      NLM  = 2*IORDER + 1
      NXYZ = (IORDER + 1)*(IORDER + 2)/2
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from SPHTRA','*',103)
         WRITE (LUPRI,'(2X,A,I5)') ' Multipole order:     ',IORDER
         WRITE (LUPRI,'(2X,A,I5)') ' Spherical components:',NLM
         WRITE (LUPRI,'(2X,A,I5)') ' Cartesian components:',NXYZ
         WRITE (LUPRI,'(2X,A,I5)') ' NBAST:  ',NBAST
         WRITE (LUPRI,'(2X,A,I5)') ' NELMNT: ',NELMNT
         WRITE (LUPRI,'(2X,A,I5)') ' NOPTYP: ',NOPTYP
         CALL AROUND('Cartesian integrals in SPHTRA.')
         DO 100 I = 1, NXYZ
            WRITE (LUPRI,'(//,2X,A,I5)') ' Cartesian component:',I
            CALL OUTPAK(AINT(1,I),NBAST,1,LUPRI)
  100    CONTINUE
      END IF
      KTRA = 1
      KINT = KTRA  + NLM*NXYZ
      KWRK = KINT  + NELMNT*NLM
      LWRK = LWORK - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('SPHTRA',' ',KWRK,LWORK)
      CALL SPHTR1(AINT(1,1),WORK(KINT),IORDER,WORK(KTRA),WORK(KWRK),
     &            LWRK,NXYZ,NLM,NELMNT,IPRINT)
      IF (IPRINT .GT. 5) THEN
         CALL AROUND('Spherical integrals in SPHTRA.')
         DO 200 I = 1, NLM
            WRITE (LUPRI,'(//,2X,A,I5)') ' Spherical component:',I
            CALL OUTPAK(AINT(1,I),NBAST,1,LUPRI)
  200    CONTINUE
      END IF
      RETURN
      END
C  /* Deck sphtr1 */
      SUBROUTINE SPHTR1(CARINT,SPHINT,IORDER,TRAMAT,WORK,LWORK,
     &                  NXYZ,NLM,NELMNT,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.D0)
      DIMENSION TRAMAT(NXYZ,NLM), WORK(LWORK),
     &          CARINT(NELMNT,NXYZ), SPHINT(NELMNT,NLM)
      CALL SPHCOM(IORDER,TRAMAT,NLM,NXYZ,0,0,WORK,LWORK,IPRINT)
      CALL DZERO(SPHINT,NELMNT*NLM)
      DO 100 I = 1, NLM
         DO 200 J = 1, NXYZ
            COEF = TRAMAT(J,I)
            IF (ABS(COEF) .GT. D0) THEN
               IF (IPRINT .GT. 5) WRITE (LUPRI,'(1X,A,2I5,F12.6)')
     &                               ' I, J , COEF ', I, J, COEF
               DO 300 K = 1, NELMNT
                  SPHINT(K,I) = SPHINT(K,I) + COEF*CARINT(K,J)
  300          CONTINUE
            END IF
  200    CONTINUE
  100 CONTINUE
      CALL DCOPY(NLM*NELMNT,SPHINT,1,CARINT,1)
      RETURN
      END
C  /* Deck sphcom */
      SUBROUTINE SPHCOM(LVAL,TRAMAT,NLM,NXYZ,MORDER,MINTEG,WORK,LWORK,
     &                  IPRINT)
C
C     This routine generates coefficients for transforming from
C     Cartesian to spherical components. It is based on notes by
C     P. Wormer, September 90.
C
C     MORDER eq 0 : M_l order 0, +1, -1, +2, -2, ..., +LVAL, -LVAL
C     MORDER ne 0 : M_l order -LVAL, ..., -1, 0, 1, ..., +LVAL
C
C     MINTEG eq 0 : old Hermit normalization
C     MINTEG eq 1 : the smallest coefficient is one (abs.val.).
C     MINTEG eq 2 : the spherical components are normalized to unity 
C
C     tuh 11.10.90
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION TRAMAT(NXYZ,NLM), WORK(LWORK)
C
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from SPHCOM','*',103)
         WRITE (LUPRI,'(5X,A,I5)') ' LVAL  ', LVAL
         WRITE (LUPRI,'(5X,A,I5)') ' NLM   ', NLM
         WRITE (LUPRI,'(5X,A,I5)') ' NXYZ  ', NXYZ
         WRITE (LUPRI,'(5X,A,I5)') ' MINTEG', MINTEG
      END IF
      CALL FLSHFO(LUPRI)
C
      LCS   = 1
      LPC   = LCS + (LVAL + 1)*(LVAL + 1)
      KWORK = LPC + LVAL + 1
      IF (KWORK .GT. LWORK) CALL STOPIT('SPHCOM',' ',KWORK,LWORK)
      CALL SPHCO1(LVAL,TRAMAT,NLM,NXYZ,WORK(LCS),WORK(LPC),
     &            MORDER,MINTEG,IPRINT)
      RETURN
      END
C  /* Deck sphco1 */
      SUBROUTINE SPHCO1(LVAL,TRAMAT,NLM,NXYZ,COSSIN,PL,MORDER,MINTEG,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, D3 = 3.0D0)
      DIMENSION TRAMAT(NXYZ,NLM), COSSIN(0:LVAL,0:LVAL), PL(0:LVAL)
C
#include <sdpre.h>
C
      LVAL1 = LVAL + 1
C
C     Legendre coefficients
C     ---------------------
C
      CALL DZERO(PL(0),LVAL1)
      DO 100 K = 0, LVAL/2
         PL(LVAL-2*K) = (SDPRE((-1)**K)/SDPRE(2**LVAL))
     &                    *BINOM(LVAL,K)*BINOM(2*(LVAL-K),LVAL)
  100 CONTINUE
C
C     Cosine and sine coefficients
C     ----------------------------
C
      CALL DZERO(COSSIN(0,0),LVAL1*LVAL1)
      DO 200 M = 0, LVAL
         COSSIN(M,0) = D1
         DO 210 K = 1, M
            COSSIN(M,K) = COSSIN(M-1,K-1)*SDPRE((-1)**(K-1))
            IF (M .GT. K) COSSIN(M,K) = COSSIN(M,K) + COSSIN(M-1,K)
  210    CONTINUE
  200 CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND('Legendre polynomial')
         CALL OUTPUT(PL(0),1,1,1,LVAL1,1,LVAL1,1,LUPRI)
         CALL AROUND('Cosine and sine factors')
         CALL OUTPUT(COSSIN(0,0),1,LVAL1,1,LVAL1,LVAL1,LVAL1,1,LUPRI)
      END IF
      CALL FLSHFO(LUPRI)
C
C     Transformation coefficients
C     ---------------------------
C
      CALL DZERO(TRAMAT,NXYZ*NLM)
      DO 300 M = 0, LVAL
         CM = SQRT(D2*FACULT(LVAL-M)/FACULT(LVAL+M))
         IF (M .EQ. 0) CM = D1
         IF (MINTEG.EQ.2) CM = CM/SQRT(FACUL2(2*LVAL-1))
         DO 400 K = MOD(LVAL - M,2), LVAL - M, 2
            IF (M .GT. 0) PL(K) = SDPRE((K+1))*PL(K+1)
            CMK = CM*PL(K)
            DO 500 I = 0, (LVAL - K - M)/2
               CMKI = CMK*BINOM((LVAL - K - M)/2,I)
               DO 600 J = 0, I
                  CMKIJ = CMKI*BINOM(I,J)
                  DO 700 N = 0, M
                    IX = LVAL - 2*J - M + N
                    IX = IX*(IX + 1)/2 + LVAL1 - M - 2*I
                    IF (MORDER .EQ. 0) THEN
                       ILM = MAX(1,2*M + MOD(N,2))
                    ELSE
                       IF (MOD(N,2) .EQ. 1) THEN
                          ILM = 1 + LVAL - M
                       ELSE
                          ILM = 1 + LVAL + M
                       END IF
                    END IF
                    TRAMAT(IX,ILM) = TRAMAT(IX,ILM) + CMKIJ*COSSIN(M,N)
  700             CONTINUE
  600          CONTINUE
  500       CONTINUE
  400    CONTINUE
  300 CONTINUE
C
C     Renormalize if requested with MINTEG
C
      IF (MINTEG.EQ.1) THEN
         IOFF = LVAL*LVAL
         DO 800 I = 1, 2*LVAL + 1
            TMIN = TRAMAT(IDAMAX(NXYZ,TRAMAT(1,I),1),I)
            DO 810 J = 1, NXYZ
               TJI = ABS(TRAMAT(J,I))
               IF ((TJI.GT.D0).AND.(TJI.LT.TMIN)) TMIN = TJI
  810       CONTINUE
            TMIN = D1 / TMIN
            CALL DSCAL(NXYZ,TMIN,TRAMAT(1,I),1)
  800    CONTINUE
      END IF
      IF (IPRINT .GT. 4) THEN
         CALL AROUND('Cartesian to spherical transformation matrix')
         WRITE (LUPRI,'(29X,A,I2)') ' Moment order:',LVAL
         IXYZ = (LVAL+1)*(LVAL+2)/2
         ILM  = 2*LVAL + 1
         CALL OUTPUT(TRAMAT,1,IXYZ,1,ILM,NXYZ,NLM,1,LUPRI)
      END IF
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck facult */
      FUNCTION FACULT(N)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1=1.D0)
#include <sdpre.h>
      IF (N .LT. 0) THEN
         WRITE (LUPRI,'(/,A,I10,/A)')
     &         ' Argument less than zero in FACULT:',N,
     &         ' Program cannot continue.'
         CALL QUIT('Illegal argument in FACULT')
      ELSE
        FACULT = D1
        DO 100 I = 1, N
           FACULT = FACULT*SDPRE(I)
  100   CONTINUE
      END IF
      RETURN
      END
C  /* Deck facul2 */
      FUNCTION FACUL2(N)
C
C     tuh
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1=1.D0)
#include <sdpre.h>
C
C     N < 0
C
      IF (N .LT. 0) THEN
         FACUL2 = SDPRE(N + 2)
         DO I = N + 4, 1, 2
            FACUL2 = FACUL2*SDPRE(I)
         END DO
         IF (FACUL2 .EQ. SDPRE(0)) THEN
            WRITE (LUPRI,'(/,A,I10,/A)')
     &            ' Double factorial undefined for ',N, 
     &            ' Program cannot continue.'
            CALL QUIT('Illegal argument in FACUL2')
         ELSE
            FACUL2 = D1/FACUL2
         END IF
C
C     N = 0
C
      ELSE IF (N.EQ.0) THEN
         FACUL2 = D1
C
C     N > 0
C
      ELSE
        FACUL2 = SDPRE(N) 
        DO I = N - 2, 1, -2 
           FACUL2 = FACUL2*SDPRE(I)
        END DO
      END IF
      RETURN
      END
C  /* Deck binom */
      FUNCTION BINOM(I,J)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1=1.D0)
      IF (I .LT. J) THEN
         WRITE (LUPRI,'(/,A,2I5,/A)')
     &         ' Second argument larger than first argument in BINOM:',
     &         I,J,' Program cannot continue.'
         CALL QUIT('Illegal arguments in BINOM')
      ELSE
        BINOM = FACULT(I)/(FACULT(I-J)*FACULT(J))
      END IF
      RETURN
      END
C  /* Deck wrtpro */
      SUBROUTINE WRTPRO(AINT,LENGTH,LABEL,RTNLBL)
C
C     290689 Henrik Koch
C     241189 tuh
C     080390 OV
C     250691 tuh - square matrices
C
C     Purpose: Write integrals on property file.
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <inftap.h>
#include <gnrinf.h>
C
      CHARACTER*8 LABEL,RTNLBL(2)
      LOGICAL FIRST
      DIMENSION AINT(LENGTH)
C
C     Decide whether we start to write at the beginning or the end of the
C     file. Please note that we now require to be open when reaching this
C     routine.
C
      IF (NEWPRP) THEN
         NEWPRP = .FALSE.
         REWIND LUPROP
      ELSE
         REWIND (LUPROP)
         CALL FNDLAB('EOFLABEL',LUPROP)
         BACKSPACE LUPROP
      END IF
C
C     Write label
C
C     CALL NEWLAB(LABEL,LUPROP,LUPRI)
      CALL NEWLB2(LABEL,RTNLBL,LUPROP,LUPRI)
C
C     Write integrals
C
      LEN = MAX(4,LENGTH)
      CALL WRITI(LUPROP,IRAT*LEN,AINT)
C
C     We add an extra label to signify EOF
C
      CALL NEWLB2('EOFLABEL',RTNLBL,LUPROP,LUPRI)
C      
      RETURN
      END
C  /* Deck dsotst */
      SUBROUTINE DSOTST(SOINT,DIFFER,NBAST,NELMNT,NOPTYP,LABINT,DOATOM,
     &                  NPQUAD,INTADR,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxmom.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
      PARAMETER (D0 = 0.D0)
      DIMENSION SOINT(NELMNT,NOPTYP), DIFFER(NELMNT), INTADR(NOPTYP)
      LOGICAL DOATOM(NUCIND), SAME
      CHARACTER LABINT(9*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
C
      DIFMAX = D0
C
C     First atom
C
      DO 100 IREPO = 0, MAXREP
         DO 200 IATOM1 = 1, NUCIND
         IF (DOATOM(IATOM1)) THEN
C
C           Second atom
C
            DO 400 IATOM2 = 1, IATOM1
            IF (DOATOM(IATOM2)) THEN
               SAME = IATOM1 .EQ. IATOM2
C
C              Cartesian directions
C
               DO 500 ICOOR1 = 1, 3
                  ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
                  IF (ISCOR1 .GT. 0) THEN
                     MXCR2 = 3
                     IF (SAME) MXCR2 = ICOOR1
                     DO 600 ICOOR2 = 1, MXCR2
                        ISCOR2 = IPTCNT(3*(IATOM2-1)+ICOOR2,IREPO,2)
                        IF (ISCOR2 .GT. 0) THEN
                          IJ = INTADR(3*NUCDEP*(ISCOR1 - 1) + ISCOR2)
                          JI = INTADR(3*NUCDEP*(ISCOR2 - 1) + ISCOR1)
                          DFMAX = D0
                          DO 700 I = 1, NELMNT
                            DIFFER(I) = ABS(SOINT(I,IJ)-SOINT(I,JI))
                            DFMAX = MAX(DFMAX,DIFFER(I))
  700                     CONTINUE
                          DIFMAX = MAX(DIFMAX,DFMAX)
                          IF (IPRINT .GT. 4) THEN
                            CALL AROUND('Difference between '
     &                         //LABINT(IJ)//' and '//LABINT(JI))
                            WRITE (LUPRI,'(A,E12.6)')
     &                      ' Largest difference for this matrix: ',
     &                       DFMAX
                            CALL OUTPAK(DIFFER,NBAST,1,LUPRI)
                          END IF
                        END IF
  600                CONTINUE
                  END IF
  500          CONTINUE
            END IF
  400       CONTINUE
         END IF
  200    CONTINUE
  100 CONTINUE
      WRITE (LUPRI,'(/A,I4,A,E12.6)')
     &   ' Largest difference found for',NPQUAD,
     &   ' quadrature points: ',DIFMAX
      RETURN
      END
C  /* Deck alftyp */
      SUBROUTINE ALFTYP(INTTYP,NOPTYP,INTREP,LABINT,INTADR,DOATOM,NATOM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
      CHARACTER LABINT(9*MXCENT)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
      NOPTYP = 9*NATOM
      CALL NTYPTS(NOPTYP)
      ITYP = 0
      DO 100 IREP = 0, MAXREP
         DO 200 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
               DO 300 ICOOR1 = 1, 3               
                  ISCOR1 = IPTCNT(3*(IATOM - 1) + ICOOR1,IREP,2)
                  IF (ISCOR1 .GT. 0) THEN
                     DO 400 ICOOR2 = 1, 3
                        ITYP = ITYP + 1                        
                        LABINT(ITYP) = 'ALF '//CHRXYZ(ICOOR1)//
     &                       CHRXYZ(ICOOR2)//
     &                       CHRNOS(ITYP/10)//CHRNOS(MOD(ITYP,10))
                        INTREP(ITYP) = IBTXOR(IREP,ISYMAX(ICOOR2,2))
                        ISCOOR = 3*(ISCOR1 - 1) + ICOOR2
                        INTADR(ISCOOR) = ITYP
 400                 CONTINUE
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
