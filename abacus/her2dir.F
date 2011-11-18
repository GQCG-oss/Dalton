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
C  /* Deck direct */
      SUBROUTINE DIRECT(BIGVEC,ICENTA,ICENTB,ICENTC,ICENTD,
     &                  NCENTA,NCENTB,NCENTC,NCENTD,
     &                  ISOA,ISOB,ISOC,ISOD,
     &                  SIGNAX,SIGNAY,SIGNAZ,SIGNBX,SIGNBY,SIGNBZ,
     &                  SIGNCX,SIGNCY,SIGNCZ,SIGNDX,SIGNDY,SIGNDZ,
     &                  NCCINT,MAXDER,EXPECT,LONDON,SPNORB,DIA2SO,
     &                  ZFS2EL,EXPGRA,IATOM,MULATM,NINTYP,IPRINT)
C
C     This subroutine determines type of two-electron integral
C     and sets directives for integral calculation.
C
C     tuh
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      IMPLICIT REAL (F,S), LOGICAL(D)
      REAL ONE, TWO
      PARAMETER (ONE = 1.00 E00, TWO = 2.00 E00)
C
#endif
#if !defined (SYS_CRAY) && !defined (SYS_T3D) && !defined (VAR_REAL) && !defined (SYS_T90)
      IMPLICIT DOUBLE PRECISION (F,S), LOGICAL(D)
      DOUBLE PRECISION ONE, TWO
      PARAMETER (ONE = 1.00 D00, TWO = 2.00 D00)
C
#endif
#include <priunit.h>
      LOGICAL ONECEN, DERECC, DERECD, EXPECT, BIGVEC, SPNORB, DIA2SO,
     &        LONDON, ZFS2EL, EXPGRA
      CHARACTER*4 ATYPE
#include <r12int.h>
#include <doxyz.h>
#include <crossd.h>
#include <subdir.h>
#include <crsdir.h>
#include <expcom.h>
#include <dirprt.h>
#include <derzer.h>
#include <drw2el.h>
      COMMON /DEODIR/ DEOVEC(54)
      COMMON /DETDIR/ DETVEC(36)
      COMMON /DHODIR/ DHOVEC(18)
      COMMON /DEOADR/ IEOVEC(54)
      COMMON /DETADR/ IETVEC(36)
      COMMON /DHOADR/ IHOVEC(18)
      COMMON /DEOFAC/ FEOVEC(54)
      COMMON /DETFAC/ FETVEC(36)
      COMMON /DHOFAC/ FHOVEC(18)
C
C     **********************************************
C     *************** Initialization ***************
C     **********************************************
C
C     Number of centers
C
      ONECEN = .FALSE.
      TWOCEN = .FALSE.
      THRCEN = .FALSE.
      FOUCEN = .FALSE.
C
C     Expansion coefficient directives
C
      DERECC = .FALSE.
      DERECD = .FALSE.
C
C     Subroutine directives
C
      DZER  = .FALSE.
      DC101 = .FALSE.
      DC1H1 = .FALSE.
      DC1E1 = .FALSE.
      DC2H1 = .FALSE.
      DC2E1 = .FALSE.
      DC102 = .FALSE.
      DC1H2 = .FALSE.
      DC1E2 = .FALSE.
      DC2H2 = .FALSE.
      DC2E2 = .FALSE.
C
C     Directives for cross differentiation
C
      CROSS1 = .FALSE.
      CROSS2 = .FALSE.
      DHCHX  = .FALSE.
      DHCHY  = .FALSE.
      DHCHZ  = .FALSE.
      DHCEX  = .FALSE.
      DHCEY  = .FALSE.
      DHCEZ  = .FALSE.
      DHCEX1 = .FALSE.
      DHCEY1 = .FALSE.
      DHCEZ1 = .FALSE.
      DHCEX2 = .FALSE.
      DHCEY2 = .FALSE.
      DHCEZ2 = .FALSE.
C
C     Calculation paths
C
      DPATH1 = .FALSE.
      DPATH2 = .FALSE.
C
C     Derivative types DEOVEC, DETVEC and DHOVEC
C
      DO 100 I = 1, 54
         DEOVEC(I) = .FALSE.
  100 CONTINUE
      DO 110 I = 1, 36
         DETVEC(I) = .FALSE.
  110 CONTINUE
      DO 120 I = 1, 18
         DHOVEC(I) = .FALSE.
  120 CONTINUE
C
C     Orbital-exponent derivatives
C
      IF (EXPGRA) THEN
         ATYPE  = 'EXPG'
         DPATH1 = .TRUE.
         DPATH2 = .TRUE.
         DC101  = .TRUE.
         DC2H1  = .TRUE.
         DC102  = .TRUE.
         DC2H2  = .TRUE.
         NINTYP = 4
         NTOTAL = NINTYP*NCCINT
         GO TO 900
      ELSE IF (BPH2OO) THEN
         ATYPE  = 'BPOO'
         DPATH1 = .TRUE.
         DC2H1  = .TRUE.
         NINTYP = 1
         NTOTAL = NINTYP*NCCINT
         GO TO 900
C
C     Integrals for r12 methods (WK/UniKA/15-11-2002)
C
      ELSE IF (R12EIN .OR. R12INT .OR. U12INT) THEN
         ATYPE  = 'R12 '
         DPATH1 = .TRUE.
         DZER   = .TRUE.
         FZERO  = ONE
         IZERO  = 1
         DC101  = .TRUE.
         DC2H1  = .TRUE.
         NINTYP = 0
         IF (R12EIN) THEN
            NINTYP = NINTYP + 1
         ELSE
            IF (V12INT) NINTYP = NINTYP + 1
            IF (R12INT) NINTYP = NINTYP + 1
            IF (U21INT) NINTYP = NINTYP + 1
            IF (U12INT) NINTYP = NINTYP + 1
         END IF
         IF (U12INT) THEN
            DPATH2 = .TRUE.
            DC102  = .TRUE.
         END IF
         NTOTAL = NCCINT
         GO TO 900
C
C     Spin-orbit integrals
C
      ELSE IF (SPNORB) THEN
         ATYPE  = 'SO  '
         DPATH1 = .TRUE.
         DC101  = .TRUE.
         DC2E1  = .TRUE.
         NINTYP = 3
         NTOTAL = NINTYP*NCCINT
         CALL DIRESO(1,2,3,ONE,ONE,ONE)
         GO TO 900
C
C     Diamagnetic spin-orbit integrals
C
      ELSE IF (DIA2SO) THEN
         ATYPE  = 'DSO2'
         DPATH1 = .TRUE.
         CROSS1 = .TRUE.
         DC101  = .TRUE.
         DC1E1  = .TRUE.
         DC2E1  = .TRUE.
C
         DHCEX  = .TRUE.
         DHCEY  = .TRUE.
         DHCEZ  = .TRUE.
         DHCEX1 = .TRUE.
         DHCEY1 = .TRUE.
         DHCEZ1 = .TRUE.
C
         NCENT1 = NCENTA
         NCENT2 = NCENTB
         NCENT3 = NCENTC
         NCENT4 = NCENTD
         ICENT1 = ICENTA
         ICENT2 = ICENTB
         ICENT3 = ICENTC
         ICENT4 = ICENTD
         ISO1   = ISOA
         ISO2   = ISOB
         ISO3   = ISOC
         ISO4   = ISOD
         SIGN1X = SIGNAX
         SIGN2X = SIGNBX
         SIGN3X = SIGNCX
         SIGN4X = SIGNDX
         SIGN1Y = SIGNAY
         SIGN2Y = SIGNBY
         SIGN3Y = SIGNCY
         SIGN4Y = SIGNDY
         SIGN1Z = SIGNAZ
         SIGN2Z = SIGNBZ
         SIGN3Z = SIGNCZ
         SIGN4Z = SIGNDZ
C
         CALL DIREAC(1,2,3,4,5,6,7,8,9,
     &               ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
         NINTYP = 9
         NTOTAL = NINTYP*NCCINT
         GO TO 900
C
C     Zero-field electronic spin-spin interaction integrals
C
      ELSE IF (ZFS2EL) THEN
         DPATH2 = .TRUE.
         DC102 = .TRUE.
         DC2E2 = .TRUE.
C
         NCENT1 = NCENTA
         NCENT2 = NCENTB
         NCENT3 = NCENTC
         NCENT4 = NCENTD
         ICENT1 = ICENTA
         ICENT2 = ICENTB
         ICENT3 = ICENTC
         ICENT4 = ICENTD
C
         ISO1   = ISOA
         ISO2   = ISOB
         ISO3   = ISOC
         ISO4   = ISOD
C
         SIGN1X = SIGNAX
         SIGN2X = SIGNBX
         SIGN3X = SIGNCX
         SIGN4X = SIGNDX
         SIGN1Y = SIGNAY
         SIGN2Y = SIGNBY
         SIGN3Y = SIGNCY
         SIGN4Y = SIGNDY
         SIGN1Z = SIGNAZ
         SIGN2Z = SIGNBZ
         SIGN3Z = SIGNCZ
         SIGN4Z = SIGNDZ
C
         CALL DIREAA(1,2,3,4,5,6,ONE,ONE,ONE,ONE,ONE,ONE)
         CALL DIREBB(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
         CALL DIREAB(13,14,15,16,17,18,19,20,21,
     *               ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
         NINTYP = 21
C
C     London integrals
C
      ELSE IF (LONDON) THEN
         ATYPE  = 'LNDN'
         DPATH1 = .TRUE.
         DPATH2 = .TRUE.
         DC101  = .TRUE.
         DC102  = .TRUE.
         DC2E1  = .TRUE.
         DC2E2  = .TRUE.
         IF (MAXDER .EQ. 1) THEN
            CALL DIREA0(1,2,3,ONE,ONE,ONE)
            CALL DIREC0(4,5,6,ONE,ONE,ONE)
            NINTYP = 6
         ELSE
            CROSS1 = .TRUE.
            DC1E1  = CROSS1
            DC1E2  = CROSS2
            DHCEX  = .TRUE.
            DHCEY  = .TRUE.
            DHCEZ  = .TRUE.
            DHCEX1 = .TRUE.
            DHCEY1 = .TRUE.
            DHCEZ1 = .TRUE.
            NCENT1 = NCENTA
            NCENT2 = NCENTB
            NCENT3 = NCENTC
            NCENT4 = NCENTD
            ICENT1 = ICENTA
            ICENT2 = ICENTB
            ICENT3 = ICENTC
            ICENT4 = ICENTD
            ISO1   = ISOA
            ISO2   = ISOB
            ISO3   = ISOC
            ISO4   = ISOD
            SIGN1X = SIGNAX
            SIGN2X = SIGNBX
            SIGN3X = SIGNCX
            SIGN4X = SIGNDX
            SIGN1Y = SIGNAY
            SIGN2Y = SIGNBY
            SIGN3Y = SIGNCY
            SIGN4Y = SIGNDY
            SIGN1Z = SIGNAZ
            SIGN2Z = SIGNBZ
            SIGN3Z = SIGNCZ
            SIGN4Z = SIGNDZ
            CALL DIREA0(1,2,3,ONE,ONE,ONE)
            CALL DIREC0(4,5,6,ONE,ONE,ONE)
            CALL DIREAA( 7, 8, 9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
            CALL DIRECC(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
            CALL DIREAC(19,20,21,22,23,24,25,26,27,
     &                  ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
            NINTYP = 27
         END IF
         NTOTAL = NINTYP*NCCINT
         GO TO 900
      ELSE IF (MAXDER .EQ. 0) THEN
         ATYPE  = 'UNDF'
         DPATH1 = .TRUE.
         DZER   = .TRUE.
         FZERO  = ONE
         IZERO  = 1
         DC101  = .TRUE.
         DC2H1  = .TRUE.
         NINTYP = 1
         NTOTAL = NCCINT
         RETURN
      ELSE IF (MAXDER .EQ. 1) THEN
         DERONE = .TRUE.
         DERTWO = .FALSE.
      ELSE
         DERONE = .TRUE.
         DERTWO = .TRUE.
      END IF
C     If only derivatives with respect to IATOM are needed,
C     check to see if this atom contributes to the integral.
C
      IF (.NOT. EXPECT) THEN
         IF (NCENTA .NE. IATOM .AND. NCENTB .NE. IATOM .AND.
     *       NCENTC .NE. IATOM .AND. NCENTD .NE. IATOM) RETURN
      END IF
C
C     Special case: PERTURB and differentiation on symmetry
C     related centers. All atoms are differentiated separately.
C
      IF (.NOT. EXPECT) THEN
         ATYPE  = '1STD'
         IF (MULATM .EQ. 1) THEN
            NOATMS = 1
            IATOMS(1) = IATOM
            TRANS = .FALSE.
            IF (IATOM .EQ. NCENTA) THEN
               ISOPDR(1)   = ISOA
               ICNTDR(1)   = ICENTA
               NCNTDR(1)   = NCENTA
               SIGNDR(1,1) = SIGNAX
               SIGNDR(1,2) = SIGNAY
               SIGNDR(1,3) = SIGNAZ
            ELSE IF (IATOM .EQ. NCENTB) THEN
               ISOPDR(1)   = ISOB
               ICNTDR(1)   = ICENTB
               NCNTDR(1)   = NCENTB
               SIGNDR(1,1) = SIGNBX
               SIGNDR(1,2) = SIGNBY
               SIGNDR(1,3) = SIGNBZ
            ELSE IF (IATOM .EQ. NCENTC) THEN
               ISOPDR(1)   = ISOC
               ICNTDR(1)   = ICENTC
               NCNTDR(1)   = NCENTC
               SIGNDR(1,1) = SIGNCX
               SIGNDR(1,2) = SIGNCY
               SIGNDR(1,3) = SIGNCZ
            ELSE
               ISOPDR(1)   = ISOD
               ICNTDR(1)   = ICENTD
               NCNTDR(1)   = NCENTD
               SIGNDR(1,1) = SIGNDX
               SIGNDR(1,2) = SIGNDY
               SIGNDR(1,3) = SIGNDZ
            END IF
         ELSE
            DOA = IATOM .EQ. NCENTA
            DOB = IATOM .EQ. NCENTB
            DOC = IATOM .EQ. NCENTC
            DOD = IATOM .EQ. NCENTD
            TRANS = .FALSE.
C           DOD = IATOM .EQ. NCENTD .AND. .NOT.(DOA.AND.DOB.AND.DOC)
C           TRANS = IATOM .EQ. NCENTD .AND. (DOA.AND.DOB.AND.DOC)
            INT = 0
            IF (DOA) THEN
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOA
               ICNTDR(INT)   = ICENTA
               NCNTDR(INT)   = NCENTA
               SIGNDR(INT,1) = SIGNAX
               SIGNDR(INT,2) = SIGNAY
               SIGNDR(INT,3) = SIGNAZ
            END IF
            IF (DOB) THEN
               CALL DIREB0(3*INT+1,3*INT+2,3*INT+3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOB
               ICNTDR(INT)   = ICENTB
               NCNTDR(INT)   = NCENTB
               SIGNDR(INT,1) = SIGNBX
               SIGNDR(INT,2) = SIGNBY
               SIGNDR(INT,3) = SIGNBZ
            END IF
            IF (DOC) THEN
               CALL DIREC0(3*INT+1,3*INT+2,3*INT+3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOC
               ICNTDR(INT)   = ICENTC
               NCNTDR(INT)   = NCENTC
               SIGNDR(INT,1) = SIGNCX
               SIGNDR(INT,2) = SIGNCY
               SIGNDR(INT,3) = SIGNCZ
            END IF
            IF (DOD) THEN
               CALL DIRED0(3*INT+1,3*INT+2,3*INT+3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOD
               ICNTDR(INT)   = ICENTD
               NCNTDR(INT)   = NCENTD
               SIGNDR(INT,1) = SIGNDX
               SIGNDR(INT,2) = SIGNDY
               SIGNDR(INT,3) = SIGNDZ
            END IF
            NOATMS = INT
            NINTYP = 3*INT
            IF (DOA .OR. DOB) THEN
               DPATH2 = .TRUE.
               DC102  = .TRUE.
               DC2E2  = .TRUE.
            END IF
            IF (DOC .OR. DOD) THEN
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
            END IF
            GO TO 200
         END IF
      END IF
C
C
C     ***********************************************************
C     ********* Determine type of multicenter integral **********
C     ***********************************************************
C
      IF(ICENTA.EQ.ICENTB) THEN
         IF(ICENTA.EQ.ICENTC) THEN
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1111'
               GO TO 111
            ELSE
               ATYPE = '1112'
               GO TO 212
            END IF
         ELSE
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1121'
               GO TO 212
            ELSE IF(ICENTC.EQ.ICENTD) THEN
               ATYPE = '1122'
               GO TO 211
            ELSE
               ATYPE = '1123'
               GO TO 312
            END IF
         END IF
      ELSE
         IF(ICENTA.EQ.ICENTC) THEN
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1211'
               GO TO 221
            ELSE IF(ICENTB.EQ.ICENTD) THEN
               ATYPE = '1212'
               GO TO 222
            ELSE
               ATYPE = '1213'
               GO TO 322
            END IF
         ELSE IF(ICENTB.EQ.ICENTC) THEN
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1221'
               GO TO 222
            ELSE IF(ICENTB.EQ.ICENTD) THEN
               ATYPE = '1222'
               GO TO 221
            ELSE
               ATYPE = '1223'
               GO TO 322
            END IF
         ELSE
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1231'
               GO TO 322
            ELSE IF(ICENTB.EQ.ICENTD) THEN
               ATYPE = '1232'
               GO TO 322
            ELSE IF(ICENTC.EQ.ICENTD) THEN
               ATYPE = '1233'
               GO TO 321
            ELSE
               ATYPE = '1234'
               GO TO 422
            END IF
         END IF
      END IF
C
C     ************************************
C     **************** 111 ***************
C     ************************************
C
  111 CONTINUE
         ONECEN = .TRUE.
      RETURN
C
C     ************************************
C     *************** 212 ****************
C     ************************************
C
  212 CONTINUE
         TWOCEN = .TRUE.
         DPATH1 = .TRUE.
         DC101  = .TRUE.
         DC2E1  = .TRUE.
         NINTYP = 0
C
C        ***** 1112 *****
C
         IF (ATYPE .EQ. '1112') THEN
            DERECD = .TRUE.
            IF (EXPECT) THEN
               NCENT1 = NCENTD
               NCENT2 = NCENTA
               ICENT1 = ICENTD
               ICENT2 = ICENTA
               ISO1   = ISOD
               ISO2   = ISOA
               SIGN1X = SIGNDX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNDY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNDZ
               SIGN2Z = SIGNAZ
               IF (DERONE) THEN
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREDD(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  CALL DIRED0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1121 *****
C
         ELSE
            DERECC = .TRUE.
            IF (EXPECT) THEN
               NCENT1 = NCENTC
               NCENT2 = NCENTA
               ICENT1 = ICENTC
               ICENT2 = ICENTA
               ISO1   = ISOC
               ISO2   = ISOA
               SIGN1X = SIGNCX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNCY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNCZ
               SIGN2Z = SIGNAZ
               IF (DERONE) THEN
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIRECC(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  CALL DIREC0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 221 ***************
C     ***********************************
C
  221 CONTINUE
         TWOCEN = .TRUE.
         DPATH2 = .TRUE.
         DC102 = .TRUE.
         DC2E2 = .TRUE.
         NINTYP = 0
C
C        ***** 1211 *****
C
         IF (ATYPE .EQ. '1211') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTA
               ICENT1 = ICENTB
               ICENT2 = ICENTA
               ISO1   = ISOB
               ISO2   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNAZ
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREBB(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  CALL DIREB0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1222 *****
C
         ELSE
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               ICENT1 = ICENTA
               ICENT2 = ICENTB
               ISO1   = ISOA
               ISO2   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREAA(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
               ELSE
                  CALL DIREA0(1,2,3,-ONE,-ONE,-ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 211 ***************
C     ***********************************
C
  211 CONTINUE
         TWOCEN = .TRUE.
         CALL DP211(MAXDER,DPATH1,DPATH2)
         NINTYP = 0
         IF (EXPECT) THEN
            IF (DPATH1) THEN
               NCENT1 = NCENTC
               NCENT2 = NCENTA
               ICENT1 = ICENTC
               ICENT2 = ICENTA
               ISO1   = ISOC
               ISO2   = ISOA
               SIGN1X = SIGNCX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNCY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNCZ
               SIGN2Z = SIGNAZ
               DC101  = .TRUE.
               DC2H1  = .TRUE.
               IF (DERONE) THEN
                  CALL DIREQ0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREQQ(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               NCENT1 = NCENTA
               NCENT2 = NCENTC
               ICENT1 = ICENTA
               ICENT2 = ICENTC
               ISO1   = ISOA
               ISO2   = ISOC
               SIGN1X = SIGNAX
               SIGN2X = SIGNCX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNCY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNCZ
               DC102 = .TRUE.
               DC2H2 = .TRUE.
               IF (DERONE) THEN
                  CALL DIREP0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREPP(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            END IF
         ELSE
            IF (DPATH1) THEN
               DC101 = .TRUE.
               DC2H1 = .TRUE.
               IF (IATOM .EQ. NCENTA) THEN
                  CALL DIREQ0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIREQ0(1,2,3,ONE,ONE,ONE)
               END IF
            ELSE
               DC102 = .TRUE.
               DC2H2 = .TRUE.
               IF (IATOM .EQ. NCENTA) THEN
                  CALL DIREP0(1,2,3,ONE,ONE,ONE)
               ELSE
                  CALL DIREP0(1,2,3,-ONE,-ONE,-ONE)
               END IF
            END IF
            NINTYP = NINTYP + 3
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 222 ***************
C     ***********************************
C
  222 CONTINUE
         TWOCEN = .TRUE.
         DPATH1 = .TRUE.
         DPATH2 = .TRUE.
         DC101  = .TRUE.
         DC2E1  = .TRUE.
         DC102  = .TRUE.
         DC2E2  = .TRUE.
         IF (DERTWO) THEN
            CALL DCROSS(CROSS1,CROSS2)
            DC1E1 = CROSS1
            DC1E2 = CROSS2
         END IF
         NINTYP = 0
C
C        ***** 1212 *****
C
         IF (ATYPE .EQ. '1212') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               ICENT1 = ICENTA
               ICENT2 = ICENTB
               ISO1   = ISOA
               ISO2   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  DHCEX1 = DOX
                  DHCEY1 = DOY
                  DHCEZ1 = DOZ
                  CALL DIREAA(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAC(4,5,6,5,7,8,6,8,9,
     *                        TWO,ONE,ONE,ONE,TWO,ONE,ONE,ONE,TWO)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  DERECC = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DERECD = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1221 *****
C
         ELSE
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               ICENT1 = ICENTA
               ICENT2 = ICENTB
               ISO1   = ISOA
               ISO2   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  IF (CROSS1) THEN
                     DHCEX1 = DOX
                     DHCEY1 = DOY
                     DHCEZ1 = DOZ
                  ELSE
                     DHCEX2 = DOX
                     DHCEY2 = DOY
                     DHCEZ2 = DOZ
                  END IF
                  CALL DIREAA(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAD(4,5,6,5,7,8,6,8,9,
     *                        TWO,ONE,ONE,ONE,TWO,ONE,ONE,ONE,TWO)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  DERECD = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DERECC = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 312 ***************
C     ***********************************
C
  312 CONTINUE
         THRCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            NCENT1 = NCENTC
            NCENT2 = NCENTD
            NCENT3 = NCENTA
            ICENT1 = ICENTC
            ICENT2 = ICENTD
            ICENT3 = ICENTA
            ISO1   = ISOC
            ISO2   = ISOD
            ISO3   = ISOA
            SIGN1X = SIGNCX
            SIGN2X = SIGNDX
            SIGN3X = SIGNAX
            SIGN1Y = SIGNCY
            SIGN2Y = SIGNDY
            SIGN3Y = SIGNAY
            SIGN1Z = SIGNCZ
            SIGN2Z = SIGNDZ
            SIGN3Z = SIGNAZ
            DPATH1 = .TRUE.
            DC101  = .TRUE.
            DC2E1  = .TRUE.
            DERECC = .TRUE.
            DERECD = .TRUE.
            IF (DERONE) THEN
               CALL DIREC0(1,2,3,ONE,ONE,ONE)
               CALL DIRED0(4,5,6,ONE,ONE,ONE)
               NINTYP = NINTYP + 6
            END IF
            IF (DERTWO) THEN
               CALL DIRECC(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIREDD(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIRECD(19,20,21,22,23,24,25,26,27,
     *                     ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
               NINTYP = NINTYP + 21
            END IF
         ELSE
            IF (IATOM .EQ. NCENTA) THEN
               DPATH2 = .TRUE.
               DC102  = .TRUE.
               DC2H2  = .TRUE.
               CALL DIREP0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. NCENTC) THEN
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECC = .TRUE.
               CALL DIREC0(1,2,3,ONE,ONE,ONE)
            ELSE
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECD = .TRUE.
               CALL DIRED0(1,2,3,ONE,ONE,ONE)
            END IF
            NINTYP = NINTYP + 3
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 322 ***************
C     ***********************************
C
  322 CONTINUE
         THRCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            DPATH1 = .TRUE.
            DPATH2 = .TRUE.
            DC101  = .TRUE.
            DC2E1  = .TRUE.
            DC102  = .TRUE.
            DC2E2  = .TRUE.
            IF (DERTWO) THEN
               CALL DCROSS(CROSS1,CROSS2)
               DC1E1 = CROSS1
               DC1E2 = CROSS2
            END IF
         END IF
C
C        ***** 1213 *****
C
         IF (ATYPE .EQ. '1213') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTD
               NCENT3 = NCENTA
               ICENT1 = ICENTB
               ICENT2 = ICENTD
               ICENT3 = ICENTA
               ISO1   = ISOB
               ISO2   = ISOD
               ISO3   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNDX
               SIGN3X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNDY
               SIGN3Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNDZ
               SIGN3Z = SIGNAZ
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  DHCEX2 = DOX
                  DHCEY2 = DOY
                  DHCEZ2 = DOZ
                  CALL DIREBB(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBD(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. NCENTB) THEN
                  DPATH2 = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1223 *****
C
         ELSE IF (ATYPE .EQ. '1223') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTD
               NCENT3 = NCENTB
               ICENT1 = ICENTA
               ICENT2 = ICENTD
               ICENT3 = ICENTB
               ISO1   = ISOA
               ISO2   = ISOD
               ISO3   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNDX
               SIGN3X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNDY
               SIGN3Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNDZ
               SIGN3Z = SIGNBZ
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  IF (CROSS1) THEN
                     DHCEX1 = DOX
                     DHCEY1 = DOY
                     DHCEZ1 = DOZ
                  ELSE
                     DHCEX2 = DOX
                     DHCEY2 = DOY
                     DHCEZ2 = DOZ
                  END IF
                  CALL DIREAA(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAD(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  DPATH2 = .TRUE.
                  DC102 = .TRUE.
                  DC2E2 = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. NCENTB) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1231 *****
C
         ELSE IF (ATYPE .EQ. '1231') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTC
               NCENT3 = NCENTA
               ICENT1 = ICENTB
               ICENT2 = ICENTC
               ICENT3 = ICENTA
               ISO1   = ISOB
               ISO2   = ISOC
               ISO3   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNCX
               SIGN3X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNCY
               SIGN3Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNCZ
               SIGN3Z = SIGNAZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  IF (CROSS1) THEN
                     DHCEX2 = DOX
                     DHCEY2 = DOY
                     DHCEZ2 = DOZ
                  ELSE
                     DHCEX1 = DOX
                     DHCEY1 = DOY
                     DHCEZ1 = DOZ
                  END IF
                  CALL DIREBB(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBC(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC102  = .TRUE.
                  DC2E1  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. NCENTB) THEN
                  DPATH2 = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1232 *****
C
         ELSE
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTC
               NCENT3 = NCENTB
               ICENT1 = ICENTA
               ICENT2 = ICENTC
               ICENT3 = ICENTB
               ISO1   = ISOA
               ISO2   = ISOC
               ISO3   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNCX
               SIGN3X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNCY
               SIGN3Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNCZ
               SIGN3Z = SIGNBZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  DHCEX1 = DOX
                  DHCEY1 = DOY
                  DHCEZ1 = DOZ
                  CALL DIREAA(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAC(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. NCENTA) THEN
                  DPATH2 = .TRUE.
                  DC102 = .TRUE.
                  DC2E2 = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. NCENTB) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC102  = .TRUE.
                  DC2E1  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 321 ***************
C     ***********************************
C
  321 CONTINUE
         THRCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            NCENT1 = NCENTA
            NCENT2 = NCENTB
            NCENT3 = NCENTC
            ICENT1 = ICENTA
            ICENT2 = ICENTB
            ICENT3 = ICENTC
            ISO1   = ISOA
            ISO2   = ISOB
            ISO3   = ISOC
            SIGN1X = SIGNAX
            SIGN2X = SIGNBX
            SIGN3X = SIGNCX
            SIGN1Y = SIGNAY
            SIGN2Y = SIGNBY
            SIGN3Y = SIGNCY
            SIGN1Z = SIGNAZ
            SIGN2Z = SIGNBZ
            SIGN3Z = SIGNCZ
            DPATH2 = .TRUE.
            DC102 = .TRUE.
            DC2E2 = .TRUE.
            IF (DERONE) THEN
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
               CALL DIREB0(4,5,6,ONE,ONE,ONE)
               NINTYP = NINTYP + 6
            END IF
            IF (DERTWO) THEN
               CALL DIREAA(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIREBB(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIREAB(19,20,21,22,23,24,25,26,27,
     *                     ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
               NINTYP = NINTYP + 21
            END IF
         ELSE
            IF (IATOM .EQ. NCENTA) THEN
               DPATH2 = .TRUE.
               DC102 = .TRUE.
               DC2E2 = .TRUE.
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. NCENTB) THEN
               DPATH2 = .TRUE.
               DC102 = .TRUE.
               DC2E2 = .TRUE.
               CALL DIREB0(1,2,3,ONE,ONE,ONE)
            ELSE
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2H1  = .TRUE.
               CALL DIREQ0(1,2,3,ONE,ONE,ONE)
            END IF
            NINTYP = NINTYP + 3
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 422 ***************
C     ***********************************
C
  422 CONTINUE
         FOUCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            DPATH1 = .TRUE.
            DPATH2 = .TRUE.
            DC101  = .TRUE.
            DC2E1  = .TRUE.
            DC102  = .TRUE.
            DC2E2  = .TRUE.
            CALL DP422(MAXDER,CROSS1,CROSS2)
            IF (DERTWO) THEN
               DC1E1 = CROSS1
               DC1E2 = CROSS2
            END IF
            IF (CROSS1) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTC
               NCENT3 = NCENTD
               NCENT4 = NCENTA
               ICENT1 = ICENTB
               ICENT2 = ICENTC
               ICENT3 = ICENTD
               ICENT4 = ICENTA
               ISO1   = ISOB
               ISO2   = ISOC
               ISO3   = ISOD
               ISO4   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNCX
               SIGN3X = SIGNDX
               SIGN4X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNCY
               SIGN3Y = SIGNDY
               SIGN4Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNCZ
               SIGN3Z = SIGNDZ
               SIGN4Z = SIGNAZ
               DERECC = .TRUE.
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(4,5,6,ONE,ONE,ONE)
                  CALL DIRED0(7,8,9,ONE,ONE,ONE)
                  NINTYP = NINTYP + 9
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  DHCEX2 = DOX
                  DHCEY2 = DOY
                  DHCEZ2 = DOZ
                  CALL DIREBB(10,11,12,13,14,15,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(16,17,18,19,20,21,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(22,23,24,25,26,27,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBC(28,29,30,31,32,33,34,35,36,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBD(37,38,39,40,41,42,43,44,45,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECD(46,47,48,49,50,51,52,53,54,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 45
               END IF
            ELSE
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               NCENT3 = NCENTC
               NCENT4 = NCENTD
               ICENT1 = ICENTA
               ICENT2 = ICENTB
               ICENT3 = ICENTC
               ICENT4 = ICENTD
               ISO1   = ISOA
               ISO2   = ISOB
               ISO3   = ISOC
               ISO4   = ISOD
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN3X = SIGNCX
               SIGN4X = SIGNDX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN3Y = SIGNCY
               SIGN4Y = SIGNDY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               SIGN3Z = SIGNCZ
               SIGN4Z = SIGNDZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREB0(4,5,6,ONE,ONE,ONE)
                  CALL DIREC0(7,8,9,ONE,ONE,ONE)
                  NINTYP = NINTYP + 9
               END IF
               IF (DERTWO) THEN
                  DHCEX = DOX
                  DHCEY = DOY
                  DHCEZ = DOZ
                  DHCEX1 = DOX
                  DHCEY1 = DOY
                  DHCEZ1 = DOZ
                  CALL DIREAA(10,11,12,13,14,15,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBB(16,17,18,19,20,21,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(22,23,24,25,26,27,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAB(28,29,30,31,32,33,34,35,36,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAC(37,38,39,40,41,42,43,44,45,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBC(46,47,48,49,50,51,52,53,54,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 45
               END IF
            END IF
         ELSE
            IF (IATOM .EQ. NCENTA) THEN
               DPATH2 = .TRUE.
               DC102 = .TRUE.
               DC2E2 = .TRUE.
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. NCENTB) THEN
               DPATH2 = .TRUE.
               DC102  = .TRUE.
               DC2E2  = .TRUE.
               CALL DIREB0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. NCENTC) THEN
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECC = .TRUE.
               CALL DIREC0(1,2,3,ONE,ONE,ONE)
            ELSE
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECD = .TRUE.
               CALL DIRED0(1,2,3,ONE,ONE,ONE)
            END IF
            NINTYP = NINTYP + 3
         END IF
  200 CONTINUE
      NTOTAL = NINTYP*NCCINT
  900 CONTINUE
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1005) EXPECT, IATOM
         WRITE (LUPRI, 1010) ONECEN, TWOCEN, THRCEN, FOUCEN
         WRITE (LUPRI, 1020) ATYPE
         WRITE (LUPRI, 1030) DPATH1, DPATH2
         WRITE (LUPRI, 1035) NCCINT
         WRITE (LUPRI, 1040) NTOTAL
         WRITE (LUPRI, 1045) CROSS1, CROSS2
         WRITE (LUPRI, 1050) DC101, DC102
         WRITE (LUPRI, 1060) DC1H1, DC1H2
         WRITE (LUPRI, 1070) DC1E1, DC1E2
         WRITE (LUPRI, 1080) DC2H1, DC2H2
         WRITE (LUPRI, 1090) DC2E1, DC2E2
         WRITE (LUPRI, 1095) DERECC, DERECD
         WRITE (LUPRI, 1100) (DEOVEC(I),I=1,54)
         WRITE (LUPRI, 1110) (DETVEC(I),I=1,36)
         WRITE (LUPRI, 1120) (DHOVEC(I),I=1,18)
         WRITE (LUPRI, 1130) (IEOVEC(I),I=1,54)
         WRITE (LUPRI, 1140) (IETVEC(I),I=1,36)
         WRITE (LUPRI, 1150) (IHOVEC(I),I=1,18)
         WRITE (LUPRI, 1160) (NINT(FEOVEC(I)),I=1,54)
         WRITE (LUPRI, 1170) (NINT(FETVEC(I)),I=1,36)
         WRITE (LUPRI, 1180) (NINT(FHOVEC(I)),I=1,18)
      END IF
 1000 FORMAT (/,1X,'<<<<<<<<<< SUBROUTINE DIRECT >>>>>>>>>>',/)
 1005 FORMAT (1X,'EXPECT/IATOM     ',L5,I5)
 1010 FORMAT (1X,'ONECEN,... FOUCEN',4L5)
 1020 FORMAT (1X,'INTEGRAL TYPE    ',A4)
 1030 FORMAT (1X,'DPATH1/2         ',2L5)
 1035 FORMAT (1X,'NCCINT           ',I5)
 1040 FORMAT (1X,'NTOTAL           ',I5)
 1045 FORMAT (1X,'CROSS1/2         ',2L5)
 1050 FORMAT (1X,'DC101/2          ',2L5)
 1060 FORMAT (1X,'DC1H1/2          ',2L5)
 1070 FORMAT (1X,'DC1E1/2          ',2L5)
 1080 FORMAT (1X,'DC2H1/2          ',2L5)
 1090 FORMAT (1X,'DC2E1/2          ',2L5)
 1095 FORMAT (1X,'DERECC/D         ',2L5)
 1100 FORMAT (1X,'COMMON /DEODIR/',
     *        1X,4L5,/,2(17X,4L5,/),2(17X,3L5,/),4(17X,4L5,/),
     *        2(17X,3L5,/),2(17X,4L5,/),2(17X,3L5,/))
 1110 FORMAT (1X,'COMMON /DETDIR/',1X,3L5,/,11(17X,3L5,/))
 1120 FORMAT (1X,'COMMON /DHODIR/',1X,3L5,/,17X,3L5,/,2(17X,6L5,/))
 1130 FORMAT (1X,'COMMON /DEOADR/',
     *        1X,4I5,/,2(17X,4I5,/),2(17X,3I5,/),4(17X,4I5,/),
     *        2(17X,3I5,/),2(17X,4I5,/),2(17X,3I5,/))
 1140 FORMAT (1X,'COMMON /DETADR/',1X,3I5,/,11(17X,3I5,/))
 1150 FORMAT (1X,'COMMON /DHOADR/',1X,3I5,/,17X,3I5,/,2(17X,6I5,/))
 1160 FORMAT (1X,'COMMON /DEOFAC/',
     *        1X,4I5,/,2(17X,4I5,/),2(17X,3I5,/),4(17X,4I5,/),
     *        2(17X,3I5,/),2(17X,4I5,/),2(17X,3I5,/))
 1170 FORMAT (1X,'COMMON /DETFAC/',1X,3I5,/,11(17X,3I5,/))
 1180 FORMAT (1X,'COMMON /DHOFAC/',1X,3I5,/,17X,3I5,/,2(17X,6I5,/))
      RETURN
      END
C  /* Deck dire00 */
      SUBROUTINE DIRE00
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      IMPLICIT REAL (F,S), LOGICAL (D)
#endif
#if !defined (SYS_CRAY) && !defined (SYS_T3D) && !defined (VAR_REAL) && !defined (SYS_T90)
      IMPLICIT DOUBLE PRECISION (F,S), LOGICAL (D)
#endif
#include <expcom.h>
      COMMON /DEODIR/ DEXA00, DEXB00, DEXC00, DEXD00,
     *                DEYA00, DEYB00, DEYC00, DEYD00,
     *                DEZA00, DEZB00, DEZC00, DEZD00,
     *                DEXXAA, DEXXAB, DEXXBB,
     *                DEXXCC, DEXXCD, DEXXDD,
     *                DEXYAA, DEXYAB, DEXYBA, DEXYBB,
     *                DEXYCC, DEXYCD, DEXYDC, DEXYDD,
     *                DEXZAA, DEXZAB, DEXZBA, DEXZBB,
     *                DEXZCC, DEXZCD, DEXZDC, DEXZDD,
     *                DEYYAA, DEYYAB, DEYYBB,
     *                DEYYCC, DEYYCD, DEYYDD,
     *                DEYZAA, DEYZAB, DEYZBA, DEYZBB,
     *                DEYZCC, DEYZCD, DEYZDC, DEYZDD,
     *                DEZZAA, DEZZAB, DEZZBB,
     *                DEZZCC, DEZZCD, DEZZDD
      COMMON /DETDIR/ DEXXAC, DEXYAC, DEXZAC,
     *                DEXYCA, DEYYAC, DEYZAC,
     *                DEXZCA, DEYZCA, DEZZAC,
     *                DEXXAD, DEXYAD, DEXZAD,
     *                DEXYDA, DEYYAD, DEYZAD,
     *                DEXZDA, DEYZDA, DEZZAD,
     *                DEXXBC, DEXYBC, DEXZBC,
     *                DEXYCB, DEYYBC, DEYZBC,
     *                DEXZCB, DEYZCB, DEZZBC,
     *                DEXXBD, DEXYBD, DEXZBD,
     *                DEXYDB, DEYYBD, DEYZBD,
     *                DEXZDB, DEYZDB, DEZZBD
      COMMON /DHODIR/ DHXP00, DHYP00, DHZP00,
     *                DHXQ00, DHYQ00, DHZQ00,
     *                DHXXPP, DHXYPP, DHXZPP, DHYYPP, DHYZPP, DHZZPP,
     *                DHXXQQ, DHXYQQ, DHXZQQ, DHYYQQ, DHYZQQ, DHZZQQ
      COMMON /DEOADR/ IEXA00, IEXB00, IEXC00, IEXD00,
     *                IEYA00, IEYB00, IEYC00, IEYD00,
     *                IEZA00, IEZB00, IEZC00, IEZD00,
     *                IEXXAA, IEXXAB, IEXXBB,
     *                IEXXCC, IEXXCD, IEXXDD,
     *                IEXYAA, IEXYAB, IEXYBA, IEXYBB,
     *                IEXYCC, IEXYCD, IEXYDC, IEXYDD,
     *                IEXZAA, IEXZAB, IEXZBA, IEXZBB,
     *                IEXZCC, IEXZCD, IEXZDC, IEXZDD,
     *                IEYYAA, IEYYAB, IEYYBB,
     *                IEYYCC, IEYYCD, IEYYDD,
     *                IEYZAA, IEYZAB, IEYZBA, IEYZBB,
     *                IEYZCC, IEYZCD, IEYZDC, IEYZDD,
     *                IEZZAA, IEZZAB, IEZZBB,
     *                IEZZCC, IEZZCD, IEZZDD
      COMMON /DETADR/ IEXXAC, IEXYAC, IEXZAC,
     *                IEXYCA, IEYYAC, IEYZAC,
     *                IEXZCA, IEYZCA, IEZZAC,
     *                IEXXAD, IEXYAD, IEXZAD,
     *                IEXYDA, IEYYAD, IEYZAD,
     *                IEXZDA, IEYZDA, IEZZAD,
     *                IEXXBC, IEXYBC, IEXZBC,
     *                IEXYCB, IEYYBC, IEYZBC,
     *                IEXZCB, IEYZCB, IEZZBC,
     *                IEXXBD, IEXYBD, IEXZBD,
     *                IEXYDB, IEYYBD, IEYZBD,
     *                IEXZDB, IEYZDB, IEZZBD
      COMMON /DHOADR/ IHXP00, IHYP00, IHZP00,
     *                IHXQ00, IHYQ00, IHZQ00,
     *                IHXXPP, IHXYPP, IHXZPP, IHYYPP, IHYZPP, IHZZPP,
     *                IHXXQQ, IHXYQQ, IHXZQQ, IHYYQQ, IHYZQQ, IHZZQQ
      COMMON /DEOFAC/ FEXA00, FEXB00, FEXC00, FEXD00,
     *                FEYA00, FEYB00, FEYC00, FEYD00,
     *                FEZA00, FEZB00, FEZC00, FEZD00,
     *                FEXXAA, FEXXAB, FEXXBB,
     *                FEXXCC, FEXXCD, FEXXDD,
     *                FEXYAA, FEXYAB, FEXYBA, FEXYBB,
     *                FEXYCC, FEXYCD, FEXYDC, FEXYDD,
     *                FEXZAA, FEXZAB, FEXZBA, FEXZBB,
     *                FEXZCC, FEXZCD, FEXZDC, FEXZDD,
     *                FEYYAA, FEYYAB, FEYYBB,
     *                FEYYCC, FEYYCD, FEYYDD,
     *                FEYZAA, FEYZAB, FEYZBA, FEYZBB,
     *                FEYZCC, FEYZCD, FEYZDC, FEYZDD,
     *                FEZZAA, FEZZAB, FEZZBB,
     *                FEZZCC, FEZZCD, FEZZDD
      COMMON /DETFAC/ FEXXAC, FEXYAC, FEXZAC,
     *                FEXYCA, FEYYAC, FEYZAC,
     *                FEXZCA, FEYZCA, FEZZAC,
     *                FEXXAD, FEXYAD, FEXZAD,
     *                FEXYDA, FEYYAD, FEYZAD,
     *                FEXZDA, FEYZDA, FEZZAD,
     *                FEXXBC, FEXYBC, FEXZBC,
     *                FEXYCB, FEYYBC, FEYZBC,
     *                FEXZCB, FEYZCB, FEZZBC,
     *                FEXXBD, FEXYBD, FEXZBD,
     *                FEXYDB, FEYYBD, FEYZBD,
     *                FEXZDB, FEYZDB, FEZZBD
      COMMON /DHOFAC/ FHXP00, FHYP00, FHZP00,
     *                FHXQ00, FHYQ00, FHZQ00,
     *                FHXXPP, FHXYPP, FHXZPP, FHYYPP, FHYZPP, FHZZPP,
     *                FHXXQQ, FHXYQQ, FHXZQQ, FHYYQQ, FHYZQQ, FHZZQQ
#include <doxyz.h>
C
      ENTRY DIREA0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXA00 = DOX
         DEYA00 = DOY
         DEZA00 = DOZ
         IEXA00 = INT1
         IEYA00 = INT2
         IEZA00 = INT3
         FEXA00 = FAC1
         FEYA00 = FAC2
         FEZA00 = FAC3
      RETURN
      ENTRY DIREB0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXB00 = DOX
         DEYB00 = DOY
         DEZB00 = DOZ
         IEXB00 = INT1
         IEYB00 = INT2
         IEZB00 = INT3
         FEXB00 = FAC1
         FEYB00 = FAC2
         FEZB00 = FAC3
      RETURN
      ENTRY DIREC0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXC00 = DOX
         DEYC00 = DOY
         DEZC00 = DOZ
         IEXC00 = INT1
         IEYC00 = INT2
         IEZC00 = INT3
         FEXC00 = FAC1
         FEYC00 = FAC2
         FEZC00 = FAC3
      RETURN
      ENTRY DIRED0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXD00 = DOX
         DEYD00 = DOY
         DEZD00 = DOZ
         IEXD00 = INT1
         IEYD00 = INT2
         IEZD00 = INT3
         FEXD00 = FAC1
         FEYD00 = FAC2
         FEZD00 = FAC3
      RETURN
      ENTRY DIREP0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DHXP00 = DOX
         DHYP00 = DOY
         DHZP00 = DOZ
         IHXP00 = INT1
         IHYP00 = INT2
         IHZP00 = INT3
         FHXP00 = FAC1
         FHYP00 = FAC2
         FHZP00 = FAC3
      RETURN
      ENTRY DIREQ0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DHXQ00 = DOX
         DHYQ00 = DOY
         DHZQ00 = DOZ
         IHXQ00 = INT1
         IHYQ00 = INT2
         IHZQ00 = INT3
         FHXQ00 = FAC1
         FHYQ00 = FAC2
         FHZQ00 = FAC3
      RETURN
      ENTRY DIREAA(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXAA = DOX
         DEXYAA = DOX .AND. DOY
         DEXZAA = DOX .AND. DOZ
         DEYYAA = DOY
         DEYZAA = DOY .AND. DOZ
         DEZZAA = DOZ
         IEXXAA = INT1
         IEXYAA = INT2
         IEXZAA = INT3
         IEYYAA = INT4
         IEYZAA = INT5
         IEZZAA = INT6
         FEXXAA = FAC1
         FEXYAA = FAC2
         FEXZAA = FAC3
         FEYYAA = FAC4
         FEYZAA = FAC5
         FEZZAA = FAC6
      RETURN
      ENTRY DIREBB(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXBB = DOX
         DEXYBB = DOX .AND. DOY
         DEXZBB = DOX .AND. DOZ
         DEYYBB = DOY
         DEYZBB = DOY .AND. DOZ
         DEZZBB = DOZ
         IEXXBB = INT1
         IEXYBB = INT2
         IEXZBB = INT3
         IEYYBB = INT4
         IEYZBB = INT5
         IEZZBB = INT6
         FEXXBB = FAC1
         FEXYBB = FAC2
         FEXZBB = FAC3
         FEYYBB = FAC4
         FEYZBB = FAC5
         FEZZBB = FAC6
      RETURN
      ENTRY DIRECC(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXCC = DOX
         DEXYCC = DOX .AND. DOY
         DEXZCC = DOX .AND. DOZ
         DEYYCC = DOY
         DEYZCC = DOY .AND. DOZ
         DEZZCC = DOZ
         IEXXCC = INT1
         IEXYCC = INT2
         IEXZCC = INT3
         IEYYCC = INT4
         IEYZCC = INT5
         IEZZCC = INT6
         FEXXCC = FAC1
         FEXYCC = FAC2
         FEXZCC = FAC3
         FEYYCC = FAC4
         FEYZCC = FAC5
         FEZZCC = FAC6
      RETURN
      ENTRY DIREDD(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXDD = DOX
         DEXYDD = DOX .AND. DOY
         DEXZDD = DOX .AND. DOZ
         DEYYDD = DOY
         DEYZDD = DOY .AND. DOZ
         DEZZDD = DOZ
         IEXXDD = INT1
         IEXYDD = INT2
         IEXZDD = INT3
         IEYYDD = INT4
         IEYZDD = INT5
         IEZZDD = INT6
         FEXXDD = FAC1
         FEXYDD = FAC2
         FEXZDD = FAC3
         FEYYDD = FAC4
         FEYZDD = FAC5
         FEZZDD = FAC6
      RETURN
      ENTRY DIREPP(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DHXXPP = DOX
         DHXYPP = DOX .AND. DOY
         DHXZPP = DOX .AND. DOZ
         DHYYPP = DOY
         DHYZPP = DOY .AND. DOZ
         DHZZPP = DOZ
         IHXXPP = INT1
         IHXYPP = INT2
         IHXZPP = INT3
         IHYYPP = INT4
         IHYZPP = INT5
         IHZZPP = INT6
         FHXXPP = FAC1
         FHXYPP = FAC2
         FHXZPP = FAC3
         FHYYPP = FAC4
         FHYZPP = FAC5
         FHZZPP = FAC6
      RETURN
      ENTRY DIREQQ(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DHXXQQ = DOX
         DHXYQQ = DOX .AND. DOY
         DHXZQQ = DOX .AND. DOZ
         DHYYQQ = DOY
         DHYZQQ = DOY .AND. DOZ
         DHZZQQ = DOZ
         IHXXQQ = INT1
         IHXYQQ = INT2
         IHXZQQ = INT3
         IHYYQQ = INT4
         IHYZQQ = INT5
         IHZZQQ = INT6
         FHXXQQ = FAC1
         FHXYQQ = FAC2
         FHXZQQ = FAC3
         FHYYQQ = FAC4
         FHYZQQ = FAC5
         FHZZQQ = FAC6
      RETURN
      ENTRY DIREAB(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXAB = DOX
         DEXYAB = DOX .AND. DOY
         DEXZAB = DOX .AND. DOZ
         DEXYBA = DOY .AND. DOX
         DEYYAB = DOY
         DEYZAB = DOY .AND. DOZ
         DEXZBA = DOZ .AND. DOX
         DEYZBA = DOZ .AND. DOY
         DEZZAB = DOZ
         IEXXAB = INT1
         IEXYAB = INT2
         IEXZAB = INT3
         IEXYBA = INT4
         IEYYAB = INT5
         IEYZAB = INT6
         IEXZBA = INT7
         IEYZBA = INT8
         IEZZAB = INT9
         FEXXAB = FAC1
         FEXYAB = FAC2
         FEXZAB = FAC3
         FEXYBA = FAC4
         FEYYAB = FAC5
         FEYZAB = FAC6
         FEXZBA = FAC7
         FEYZBA = FAC8
         FEZZAB = FAC9
      RETURN
      ENTRY DIREAC(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXAC = DOX
         DEXYAC = DOX .AND. DOY
         DEXZAC = DOX .AND. DOZ
         DEXYCA = DOY .AND. DOX
         DEYYAC = DOY
         DEYZAC = DOY .AND. DOZ
         DEXZCA = DOZ .AND. DOX
         DEYZCA = DOZ .AND. DOY
         DEZZAC = DOZ
         IEXXAC = INT1
         IEXYAC = INT2
         IEXZAC = INT3
         IEXYCA = INT4
         IEYYAC = INT5
         IEYZAC = INT6
         IEXZCA = INT7
         IEYZCA = INT8
         IEZZAC = INT9
         FEXXAC = FAC1
         FEXYAC = FAC2
         FEXZAC = FAC3
         FEXYCA = FAC4
         FEYYAC = FAC5
         FEYZAC = FAC6
         FEXZCA = FAC7
         FEYZCA = FAC8
         FEZZAC = FAC9
      RETURN
      ENTRY DIREAD(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXAD = DOX
         DEXYAD = DOX .AND. DOY
         DEXZAD = DOX .AND. DOZ
         DEXYDA = DOY .AND. DOX
         DEYYAD = DOY
         DEYZAD = DOY .AND. DOZ
         DEXZDA = DOZ .AND. DOX
         DEYZDA = DOZ .AND. DOY
         DEZZAD = DOZ
         IEXXAD = INT1
         IEXYAD = INT2
         IEXZAD = INT3
         IEXYDA = INT4
         IEYYAD = INT5
         IEYZAD = INT6
         IEXZDA = INT7
         IEYZDA = INT8
         IEZZAD = INT9
         FEXXAD = FAC1
         FEXYAD = FAC2
         FEXZAD = FAC3
         FEXYDA = FAC4
         FEYYAD = FAC5
         FEYZAD = FAC6
         FEXZDA = FAC7
         FEYZDA = FAC8
         FEZZAD = FAC9
      RETURN
      ENTRY DIREBC(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXBC = DOX
         DEXYBC = DOX .AND. DOY
         DEXZBC = DOX .AND. DOZ
         DEXYCB = DOY .AND. DOX
         DEYYBC = DOY
         DEYZBC = DOY .AND. DOZ
         DEXZCB = DOZ .AND. DOX
         DEYZCB = DOZ .AND. DOY
         DEZZBC = DOZ
         IEXXBC = INT1
         IEXYBC = INT2
         IEXZBC = INT3
         IEXYCB = INT4
         IEYYBC = INT5
         IEYZBC = INT6
         IEXZCB = INT7
         IEYZCB = INT8
         IEZZBC = INT9
         FEXXBC = FAC1
         FEXYBC = FAC2
         FEXZBC = FAC3
         FEXYCB = FAC4
         FEYYBC = FAC5
         FEYZBC = FAC6
         FEXZCB = FAC7
         FEYZCB = FAC8
         FEZZBC = FAC9
      RETURN
      ENTRY DIREBD(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXBD = DOX
         DEXYBD = DOX .AND. DOY
         DEXZBD = DOX .AND. DOZ
         DEXYDB = DOY .AND. DOX
         DEYYBD = DOY
         DEYZBD = DOY .AND. DOZ
         DEXZDB = DOZ .AND. DOX
         DEYZDB = DOZ .AND. DOY
         DEZZBD = DOZ
         IEXXBD = INT1
         IEXYBD = INT2
         IEXZBD = INT3
         IEXYDB = INT4
         IEYYBD = INT5
         IEYZBD = INT6
         IEXZDB = INT7
         IEYZDB = INT8
         IEZZBD = INT9
         FEXXBD = FAC1
         FEXYBD = FAC2
         FEXZBD = FAC3
         FEXYDB = FAC4
         FEYYBD = FAC5
         FEYZBD = FAC6
         FEXZDB = FAC7
         FEYZDB = FAC8
         FEZZBD = FAC9
      RETURN
      ENTRY DIRECD(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXCD = DOX
         DEXYCD = DOX .AND. DOY
         DEXZCD = DOX .AND. DOZ
         DEXYDC = DOY .AND. DOX
         DEYYCD = DOY
         DEYZCD = DOY .AND. DOZ
         DEXZDC = DOZ .AND. DOX
         DEYZDC = DOZ .AND. DOY
         DEZZCD = DOZ
         IEXXCD = INT1
         IEXYCD = INT2
         IEXZCD = INT3
         IEXYDC = INT4
         IEYYCD = INT5
         IEYZCD = INT6
         IEXZDC = INT7
         IEYZDC = INT8
         IEZZCD = INT9
         FEXXCD = FAC1
         FEXYCD = FAC2
         FEXZCD = FAC3
         FEXYDC = FAC4
         FEYYCD = FAC5
         FEYZCD = FAC6
         FEXZDC = FAC7
         FEYZDC = FAC8
         FEZZCD = FAC9
      RETURN
      ENTRY DIRESO(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXYCD = .TRUE.
         DEXZCD = .TRUE.
         DEYZCD = .TRUE.
         IEXYCD = INT3
         IEXZCD = INT2
         IEYZCD = INT1
         FEXYCD =   FAC3
         FEXZCD = - FAC2
         FEYZCD =   FAC1
      RETURN
      END
C  /* Deck dp211 */
      SUBROUTINE DP211(MAXDER,DALL2,DALL1)
#include <implicit.h>
#include <twocom.h>
      LOGICAL DALL1, DALL2
      NR(I) = (I + 1)*(I + 2)*(I + 3)/6
      LC101(LOOP,IDER) = LOOP*(NUCAB*(11 + 4*NUCCD)
     *                 + NR(MAXCD + IDER)*(10 + 12*NUCAB*NUCCD))
      LC102(LOOP,IDER) = LOOP*(NUCCD*(11 + 4*NUCAB)
     *                 + NR(MAXAB + IDER)*(10 + 12*NUCAB*NUCCD))
      LC2H1(NTYPE,LOOP) = LOOP*(NUCCD*(12 + NORBAB*4)
     *                  + NTYPE*(20 + KCKTAB*(4 + 12*NORBAB*NUCCD)))
     *                  + NTYPE*KCKTAB*KCKTCD*NORBAB*NUCCD*6
      LC2H2(NTYPE,LOOP) = LOOP*(NUCAB*(12 + NORBCD*4)
     *                  + NTYPE*(20 + KCKTCD*(4 + 12*NORBCD*NUCAB)))
     *                  + NTYPE*KCKTCD*KCKTAB*NORBCD*NUCAB*6
C
      CALL DLOOP(NHKTA,NHKTB,KCKTA,KCKTB,DIAGAB,IAB0X,IAB0Y,IAB0Z,
     *           .TRUE.,.FALSE.,.FALSE.,LOOP01,LOOP11,LOOP21)
      CALL DLOOP(NHKTC,NHKTD,KCKTC,KCKTD,DIAGCD,ICD0X,ICD0Y,ICD0Z,
     *           .TRUE.,.FALSE.,.FALSE.,LOOP02,LOOP12,LOOP22)
C
      IF (MAXDER .EQ. 1) THEN
         LALL1 = LC102(LOOP02,1) + LC2H2(3,LOOP01)
         LALL2 = LC101(LOOP01,1) + LC2H1(3,LOOP02)
      ELSE
         LALL1 = LC102(LOOP02,2) + LC2H2(9,LOOP01)
         LALL2 = LC101(LOOP01,2) + LC2H1(9,LOOP02)
      END IF
      DALL1 = LALL2 .GE. LALL1
      DALL2 = .NOT.DALL1
      RETURN
      END
C  /* Deck dp422 */
      SUBROUTINE DP422(MAXDER,DALL2,DALL1)
#include <implicit.h>
#include <twocom.h>
      LOGICAL DALL1, DALL2
      NR(I) = (I + 1)*(I + 2)*(I + 3)/6
      LC1E1(NTYPE,LOOP,IDER) = NTYPE*(LOOP*(NUCAB*(11 + 4*NUCCD)
     *                       + NR(MAXCD + IDER)*(10 + 12*NUCAB*NUCCD)))
      LC1E2(NTYPE,LOOP,IDER) = NTYPE*(LOOP*(NUCCD*(11 + 4*NUCAB)
     *                       + NR(MAXAB + IDER)*(10 + 12*NUCAB*NUCCD)))
      LC2E1(NTYPE,LOOP) = NTYPE*LOOP*(12*NUCCD + 5*NORBAB*NUCCD
     *                  + KCKTAB*(5 + 12*NORBAB*NUCCD))
     *                  + NTYPE*KCKTAB*KCKTCD*NORBAB*NUCCD*5
      LC2E2(NTYPE,LOOP) = NTYPE*LOOP*(12*NUCAB + 5*NORBCD*NUCAB
     *                  + KCKTCD*(5 + 12*NORBCD*NUCAB))
     *                  + NTYPE*KCKTCD*KCKTAB*NORBCD*NUCAB*5
C
      IF (MAXDER .EQ. 1) THEN
         CALL DLOOP(NHKTA,NHKTB,KCKTA,KCKTB,DIAGAB,IAB0X,IAB0Y,IAB0Z,
     *              .FALSE.,.TRUE.,.FALSE.,LOOP01,LOOP11,LOOP21)
         CALL DLOOP(NHKTC,NHKTD,KCKTC,KCKTD,DIAGCD,ICD0X,ICD0Y,ICD0Z,
     *              .FALSE.,.TRUE.,.FALSE.,LOOP02,LOOP12,LOOP22)
         LALL1 = LC2E2(6,LOOP11) + LC2E1(3,LOOP12)
         LALL2 = LC2E1(6,LOOP12) + LC2E2(3,LOOP11)
      ELSE
         CALL DLOOP(NHKTA,NHKTB,KCKTA,KCKTB,DIAGAB,IAB0X,IAB0Y,IAB0Z,
     *              .FALSE.,.TRUE.,.TRUE.,LOOP01,LOOP11,LOOP21)
         CALL DLOOP(NHKTC,NHKTD,KCKTC,KCKTD,DIAGCD,ICD0X,ICD0Y,ICD0Z,
     *              .FALSE.,.TRUE.,.TRUE.,LOOP02,LOOP12,LOOP22)
         LALL1 = LC2E1(6,LOOP22) + LC2E2(21,LOOP21)
     *         + LC1E2(3,LOOP12,1) + LC2E2(18,LOOP11)
         LALL2 = LC2E1(21,LOOP22) + LC2E2(6,LOOP21)
     *         + LC1E1(3,LOOP11,1) + LC2E1(18,LOOP12)
      END IF
      DALL1 = LALL2 .GE. LALL1
      DALL2 = .NOT.DALL1
      RETURN
      END
C  /* Deck dcross */
      SUBROUTINE DCROSS(CROSS1,CROSS2)
#include <implicit.h>
#include <twocom.h>
      LOGICAL CROSS1, CROSS2
      NR(I) = (I + 1)*(I + 2)*(I + 3)/6
      LC1E1(NTYPE,LOOP,IDER) = NTYPE*(LOOP*(NUCAB*(11 + 4*NUCCD)
     *                       + NR(MAXCD + IDER)*(10 + 12*NUCAB*NUCCD)))
      LC1E2(NTYPE,LOOP,IDER) = NTYPE*(LOOP*(NUCCD*(11 + 4*NUCAB)
     *                       + NR(MAXAB + IDER)*(10 + 12*NUCAB*NUCCD)))
      LC2E1(NTYPE,LOOP) = NTYPE*LOOP*(12*NUCCD + 5*NORBAB*NUCCD
     *                  + KCKTAB*(5 + 12*NORBAB*NUCCD))
     *                  + NTYPE*KCKTAB*KCKTCD*NORBAB*NUCCD*5
      LC2E2(NTYPE,LOOP) = NTYPE*LOOP*(12*NUCAB + 5*NORBCD*NUCAB
     *                  + KCKTCD*(5 + 12*NORBCD*NUCAB))
     *                  + NTYPE*KCKTCD*KCKTAB*NORBCD*NUCAB*5
C
      CALL DLOOP(NHKTA,NHKTB,KCKTA,KCKTB,DIAGAB,IAB0X,IAB0Y,IAB0Z,
     *           .FALSE.,.TRUE.,.FALSE.,LOOP01,LOOP11,LOOP21)
      CALL DLOOP(NHKTC,NHKTD,KCKTC,KCKTD,DIAGCD,ICD0X,ICD0Y,ICD0Z,
     *           .FALSE.,.TRUE.,.FALSE.,LOOP02,LOOP12,LOOP22)

      LCROS1 = LC1E1(3,LOOP11,1) + LC2E1(9,LOOP12)
      LCROS2 = LC1E2(3,LOOP12,1) + LC2E2(9,LOOP11)
      CROSS1 = LCROS2 .GE. LCROS1
      CROSS2 = .NOT.CROSS1
      RETURN
      END
C  /* Deck dloop */
      SUBROUTINE DLOOP(NHKT1,NHKT2,KCKT1,KCKT2,DIAG12,IAB0X,IAB0Y,IAB0Z,
     *                 DER0,DER1,DER2,LOOP0,LOOP1,LOOP2)
#include <implicit.h>
      LOGICAL DIAG12, DER0, DER1, DER2
#include <maxmom.h>
#include <xyzpow.h>
      IAB1X = IAB0X + 1
      IAB1Y = IAB0Y + 1
      IAB1Z = IAB0Z + 1
      DO 100 ICOMP1 = 1,KCKT1
         LVAL1 = NHKT1 - ISTEP(ICOMP1) + IAB1X
         MVAL1 = MVAL(ICOMP1) + IAB1Y
         NVAL1 = NVAL(ICOMP1) + IAB1Z
         IF (DIAG12) THEN
            MAX2 = ICOMP1
         ELSE
            MAX2 = KCKT2
         END IF
         DO 200 ICOMP2 = 1,MAX2
            LVAL12 = LVAL1 + NHKT2 - ISTEP(ICOMP2)
            MVAL12 = MVAL1 + MVAL(ICOMP2)
            NVAL12 = NVAL1 + NVAL(ICOMP2)
            NUMT0 = LVAL12/IAB1X
            NUMU0 = MVAL12/IAB1Y
            NUMV0 = NVAL12/IAB1Z
            IF (DER0) LOOP0 = NUMT0*NUMU0*NUMV0
            IF (DER1) THEN
               NUMT1 = (LVAL12 + 1)/IAB1X
               NUMU1 = (MVAL12 + 1)/IAB1Y
               NUMV1 = (NVAL12 + 1)/IAB1Z
               LOOP1 = NUMT1*NUMU0*NUMV0 + NUMT0*NUMU1*NUMV0
     *               + NUMT0*NUMU0*NUMV1
            END IF
            IF (DER2) THEN
               NUMT1 = (LVAL12 + 1)/IAB1X
               NUMU1 = (MVAL12 + 1)/IAB1Y
               NUMV1 = (NVAL12 + 1)/IAB1Z
               NUMT2 = (LVAL12 + 2)/IAB1X
               NUMU2 = (MVAL12 + 2)/IAB1Y
               NUMV2 = (NVAL12 + 2)/IAB1Z
               LOOP2 = NUMT2*NUMU0*NUMV0 + NUMT1*NUMU1*NUMV0
     *               + NUMT1*NUMU0*NUMV1 + NUMT0*NUMU2*NUMV0
     *               + NUMT0*NUMU1*NUMV1 + NUMT0*NUMU0*NUMV2
            END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
