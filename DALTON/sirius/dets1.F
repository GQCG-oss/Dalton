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
C
#ifdef UNDEF
/* Comdeck last_change */
C  CHANGES INCLUDED STARTING MAY 9, 1987 JO
C  14-May-1987 hjaaj: NTEST as parameter
C  23-JUNE-87  jo: ADDITIONS DUE TO RAS MODIFICATIONS
C  27-Jun-1987 hjaaj: IDMTYP in DENSI1
C   1-Jul-1987 hjaaj: SNGLET in CISIG1
#endif
C  /* Deck todsc */
      SUBROUTINE TODSC(A,NDIM,NBLOCK,IFIL)
C TRANSFER ARRAY REAL*8 A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
#include "implicit.h"
      DIMENSION A(NDIM)
      INTEGER START,STOP
      IF( NBLOCK .GT.0 ) THEN
C
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (A(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
      END IF
C
CIBM      IF( NBLOCK .LT. 0 ) THEN
CIBM   CALL SQFILE(IFIL,6,A,2*NDIM)
CIBM    END IF
C
      RETURN
      END
C  /* Deck cidia4 */
      SUBROUTINE CIDIA4(NAEL,NASTR,IASTR,
     &                  NBEL,NBSTR,IBSTR,
     &                  NORB,DIAG,NDET,NSTASA,NSTASB,
     &                  MAXSYM,H,
     &                  ISTBAA,ISTBAB,
     &                  ICSYM,XA,XB,SCR,RJ,RK,
     &                  SYMPRO,NSSOA,NSSOB,IOCOC,NOCTPA,NOCTPB,
     &                  ISSOA,ISSOB,ECORE,ICOOS,NTEST)
C
C CALCULATE DIAGONAL OF CIMATRIX IN DETERMINANT BASIS
C RAS VERSION
C TURBO VERSION
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION NSSOA(NOCTPA,MAXSYM),NSSOB(NOCTPB,MAXSYM)
      DIMENSION ISSOA(NOCTPA,MAXSYM),ISSOB(NOCTPB,MAXSYM)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION IASTR(NAEL,NASTR),IBSTR(NBEL,NBSTR)
      DIMENSION NSTASA(MAXSYM),NSTASB(MAXSYM)
      DIMENSION ISTBAB(MAXSYM),ISTBAA(MAXSYM)
      DIMENSION H(NORB),DIAG(NDET)
      DIMENSION RJ(NORB,NORB),RK(NORB,NORB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION ICOOS(NOCTPB,NOCTPA,MAXSYM)
      LOGICAL EQUAL
      INTEGER SYMPRO(8,8)
C
      IF (NTEST .GT. 0) WRITE(LUPRI,*) ' CIDIA4 IN USE '
      IF ( NAEL .EQ. NBEL .AND. ICSYM.EQ. 1 ) THEN
        EQUAL = .TRUE.
      ELSE
        EQUAL = .FALSE.
      END IF
C
C** 1 : OBTAIN EXCHANGE INTEGRALS (IJ!IJ) IN RK(I,J)
C              COULOMB  INTEGRALS (II!JJ) IN RJ(I,J)
C ALREADY IN IN  THIS VERSION
C
C** 2 : MODIFY ONE ELECTRON INTEGRALS ACCORDING TO
C       H(I) -1/2 * SUM(K)K(I,K)
C
      DO 100 K = 1, NORB
        CALL VECSUM(H,H,RK(1,K),1.0D0,-0.5D0,NORB)
  100 CONTINUE
C
C MODIFY BACK ( FOR THE FUN OF IT )
      DO 101 K = 1, NORB
        CALL VECSUM(H,H,RK(1,K),1.0D0,+0.5D0,NORB)
  101 CONTINUE
C
C K GOES TO J - K
      CALL VECSUM(RK,RK,RJ,-1.0D0,+1.0D0,NORB**2)
C
C
C**3 DIAGONAL ELEMENTS ACCORDING TO HANDYS FORMULAE ( CORRECTED FOR ERRO
C
C   DIAG(IDET) = HII*(NIA+NIB)
C              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
C              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
C              +         J(I,J) * NIA*NJB
C
      IDET = 0
      DO 1000 IASM = 1,MAXSYM
      IBSM = SYMPRO(IASM,ICSYM)
      DO 999  IATP = 1,NOCTPA
        DO 900 IBTP = 1,NOCTPB
        IF(IOCOC(IATP,IBTP) .NE. 1 ) GOTO 900
        IBSTRT = ISSOB(IBTP,IBSM)
        IBSTOP = IBSTRT + NSSOB(IBTP,IBSM)-1
        DO 899 IB = IBSTRT,IBSTOP
C
C TERMS DEPENDING ONLY ON IB
C
          CALL SETVEC(XB,0.0D0,NORB)
       HB = 0.0D0
       RJBB = 0.0D0
       CALL SETVEC(XB,0.0D0,NORB)
C
       DO 990 IEL = 1, NBEL
         IBEL = IBSTR(IEL,IB)
         HB = HB + H(IBEL )
C
         DO 980 JEL = 1, NBEL
           RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
  980    CONTINUE
C
         DO 970 IORB = 1, NORB
           XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
  970    CONTINUE
C
  990 CONTINUE
C
          EB = HB + 0.5D0*RJBB + ECORE
C
C
          IASTRT = ISSOA(IATP,IASM)
          IASTOP = ISSOA(IATP,IASM) + NSSOA(IATP,IASM) - 1
          DO 800 IA = IASTRT,IASTOP
            IDET = IDET + 1
C
            X1 = EB
            X2 = 0.0D0
            DO 890 IEL = 1, NAEL
              IAEL = IASTR(IEL,IA)
              X1 = X1 + ( H(IAEL )+XB(IAEL) )
C
              DO 880 JEL = 1, NAEL
               X2 = X2 + RK(IASTR(JEL,IA),IAEL )
  880         CONTINUE
C
  890 CONTINUE
C?          WRITE(LUPRI,*) ' X1 X2 ',X1,X2
            DIAG(IDET) = X1 + 0.5D0*X2
  800     CONTINUE
  899   CONTINUE
  900   CONTINUE
C
  999   CONTINUE
        IF ( EQUAL ) THEN
         DO 2000 IATP = 1, NOCTPA
           NA = NSSOA(IATP,IASM)
           IF ( NA .EQ. 0 ) GO TO 2000
           DO 1900 IBTP = 1, NOCTPB
             IF( IOCOC(IATP,IBTP) .NE. 1 ) GOTO 1900
             NB = NSSOB(IBTP,IBSM)
             IF ( NB .EQ. 0 ) GO TO 1900
             IOFF1 = ICOOS(IBTP,IATP,IASM)
             IOFF2 = ICOOS(IATP,IBTP,IBSM)
             IF(IATP .EQ. IBTP ) THEN
               CALL TRPAD(DIAG(IOFF1),1.0D0,NA)
               CALL SCALVE(DIAG(IOFF1),0.5D0,NA**2)
             ELSE
               CALL TRPMAT(DIAG(IOFF1),NA,NB,DIAG(IOFF2))
             END IF
 1900     CONTINUE
 2000    CONTINUE
       END IF
C
 1000 CONTINUE
C
      IF(NTEST .GE. 50 ) THEN
        WRITE(LUPRI,*) ' CIDIAGONAL '
        CALL WRTMAT(DIAG(1),1,NDET,1,NDET,0)
      END IF
C
      RETURN
      END
C  /* Deck rewino */
      SUBROUTINE REWINO( LU )
C
C REWIND SEQ FILE LU WITH FASTIO ROUTINES
C
C     CALL SQFILE(LU,10,IDUM,IDUM)
      REWIND LU
C
      RETURN
      END
C  /* Deck tripak */
      SUBROUTINE TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
C
C ( NOT A SIMPLIFIED VERSION OF TETRAPAK )
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C   IWAY = 2 : PACKED TO FULL FORM
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 I = 1,NDIM
          DO 50  J = 1, I
           APAK(IJ+J) = AUTPAK(J,I)
   50     CONTINUE
          IJ = IJ + I
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 I = 1,NDIM
          DO 150  J = 1, I
           AUTPAK(I,J) = APAK(IJ+J)
           AUTPAK(J,I) = APAK(IJ+J)
  150     CONTINUE
          IJ = IJ + I
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' AUTPAK AND APAK FROM TRIPAK '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM,0)
        CALL PRSYM(APAK,NDIM)
      END IF
C
      RETURN
      END
C  /* Deck frmdsc */
      SUBROUTINE FRMDSC(ARRAY,NDIM,NBLOCK,IFILE)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
#include "implicit.h"
      DIMENSION ARRAY(NDIM)
C
C
      IF( NBLOCK .GT. 0 ) THEN
      IREST=NDIM
      IBASE=0
  100 CONTINUE
       IF(IREST.GT.NBLOCK) THEN
        READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
        IBASE=IBASE+NBLOCK
        IREST=IREST-NBLOCK
       ELSE
        READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
        IREST=0
       END IF
      IF( IREST .GT. 0 ) GOTO 100
      END IF
C
C     IF( NBLOCK .LT. 0 ) THEN
C      CALL SQFILE(IFILE,7,ARRAY,2*NDIM)
C     END IF
      RETURN
      END
