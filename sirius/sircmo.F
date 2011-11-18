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
C  /* Deck cmodel */
      SUBROUTINE CMODEL(CMO,GETCMO)
C
C     Based on TRGEN in Per-Olof Widmark's SCF program.
C
C     Purpose: delete 3s component in d; 3s,4p in f; etc.
C     If (not GETCMO) then calculate NORB(8) but do not construct CMO.
C
C     Written 870710-hjaaj
C     revised 89050 -hjaaj: + GETCMO
C             890703-hjaaj: MOLFTUC4
C
#include <implicit.h>
C
      DIMENSION CMO(*)
      LOGICAL   GETCMO
C
      PARAMETER ( D0 = 0.0D0, DP5= 0.5D0, D1 = 1.0D0, D2   =   2.0D0)
      PARAMETER ( D3 = 3.0D0, D4 = 4.0D0, D6 = 6.0D0, D8   =   8.0D0)
      PARAMETER ( D10=10.0D0, D24=24.0D0, D34=34.0D0, D38  =  38.0D0)
      PARAMETER ( D40=40.0D0, D60=60.0D0, D74=74.0D0, D1270=1270.0D0)
C
C Used from common blocks:
C
C     INFINP : KDEL, TYPE(*)
C     INFORB : NSYM, NORB(8), NBAS(8)
C     PRIUNIT : IPRSTAT, LUSTAT
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
      CHARACTER*4 MOLFTUC4, CTMPMO(MXCORB)
C
      CALL QENTER('CMODEL')
      IF (IPRSTAT .GE. 99) THEN
         WRITE (LUSTAT,'(//A/)') ' >>> Test output from CMODEL'
         WRITE (LUSTAT,*) 'KDEL,GETCMO  ',KDEL,GETCMO
      END IF
C
      IF (KDEL.EQ.0) THEN
         IND = 1
         DO 100 ISYM = 1,NSYM
            NORB(ISYM) = NBAS(ISYM)
            IF (GETCMO) CALL DUNIT(CMO(IND),NBAS(ISYM))
            IND = IND + NORB(ISYM)*NBAS(ISYM)
  100    CONTINUE
      ELSE
         IF (IPRSIR .GT. 0 .AND. GETCMO) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 's-component of d-orbitals deleted'
            WRITE(LUPRI,*) 'p-component of f-orbitals deleted'
            WRITE(LUPRI,*) 's and d-components of g-orbitals deleted'
         END IF
         DO 150 JBAS = 1,NBAST
            CTMPMO(JBAS) = TYPE(JBAS)
  150    CONTINUE
         ITR0  = 0
         JBAS0 = 0
         DO 200 ISYM = 1,NSYM
            NORB(ISYM) = 0
            DO 210 JBAS = 1,NBAS(ISYM)
C
C           ... CONSTRUCT PROPER D ORBITALS
C
               IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XX  ') THEN
                  NORB(ISYM)=NORB(ISYM)+2
                  JBAS1=JBAS+1
211               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'YY  ') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 211
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
212               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'ZZ  ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 212
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 219 I=1,2*NBAS(ISYM)
                     CMO(I+ITR0) = D0
219               CONTINUE
                  CMO(JBAS +ITR0)= SQRT(DP5)
                  CMO(JBAS1+ITR0)=-SQRT(DP5)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)=-D1/SQRT(D6)
                  CMO(JBAS1+ITR0)=-D1/SQRT(D6)
                  CMO(JBAS2+ITR0)= D2/SQRT(D6)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
C
C              ... CONSTRUCT PROPER F ORBITALS
C
               ELSE IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XXX ') THEN
                  NORB(ISYM)=NORB(ISYM)+2
                  JBAS1=JBAS+1
221               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'XYY ') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 221
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
222               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'XZZ ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 222
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 229 I=1,2*NBAS(ISYM)
                     CMO(I+ITR0)=D0
229               CONTINUE
                  CMO(JBAS +ITR0)= D1/SQRT(D24)
                  CMO(JBAS1+ITR0)=-D3/SQRT(D24)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)=-D1/SQRT(D40)
                  CMO(JBAS1+ITR0)=-D1/SQRT(D40)
                  CMO(JBAS2+ITR0)= D4/SQRT(D40)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
               ELSE IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XXY ') THEN
                  NORB(ISYM)=NORB(ISYM)+2
                  JBAS1=JBAS+1
231               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'YYY ') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 231
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
232               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'YZZ ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 232
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 239 I=1,2*NBAS(ISYM)
                     CMO(I+ITR0)=D0
239               CONTINUE
                  CMO(JBAS +ITR0)=-D3/SQRT(D24)
                  CMO(JBAS1+ITR0)= D1/SQRT(D24)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)=-D1/SQRT(D40)
                  CMO(JBAS1+ITR0)=-D1/SQRT(D40)
                  CMO(JBAS2+ITR0)= D4/SQRT(D40)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
               ELSE IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XXZ ') THEN
                  NORB(ISYM)=NORB(ISYM)+2
                  JBAS1=JBAS+1
241               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'YYZ ') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 241
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
242               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'ZZZ ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 242
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 249 I=1,2*NBAS(ISYM)
                     CMO(I+ITR0)=D0
249               CONTINUE
                  CMO(JBAS +ITR0)= DP5
                  CMO(JBAS1+ITR0)=-DP5
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)=-D3/SQRT(D60)
                  CMO(JBAS1+ITR0)=-D3/SQRT(D60)
                  CMO(JBAS2+ITR0)= D2/SQRT(D60)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
C
C              ... CONSTRUCT PROPER G ORBITALS
C
               ELSE IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XXXX') THEN
                  NORB(ISYM)=NORB(ISYM)+3
                  JBAS1=JBAS+1
251               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'XXYY') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 251
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
252               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'XXZZ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 252
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
                  JBAS3=JBAS2+1
253               IF(MOLFTUC4(CTMPMO(JBAS3+JBAS0)).NE.'YYYY') THEN
                     JBAS3=JBAS3+1
                     IF(JBAS3.GT.NBAS(ISYM)) GOTO 901
                     GOTO 253
                  END IF
                  CTMPMO(JBAS3+JBAS0)='KURT'
                  JBAS4=JBAS3+1
254               IF(MOLFTUC4(CTMPMO(JBAS4+JBAS0)).NE.'YYZZ') THEN
                     JBAS4=JBAS4+1
                     IF(JBAS4.GT.NBAS(ISYM)) GOTO 901
                     GOTO 254
                  END IF
                  CTMPMO(JBAS4+JBAS0)='KURT'
                  JBAS5=JBAS4+1
255               IF(MOLFTUC4(CTMPMO(JBAS5+JBAS0)).NE.'ZZZZ') THEN
                     JBAS5=JBAS5+1
                     IF(JBAS5.GT.NBAS(ISYM)) GOTO 901
                     GOTO 255
                  END IF
                  CTMPMO(JBAS5+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 259 I=1,3*NBAS(ISYM)
                     CMO(I+ITR0)=D0
259               CONTINUE
                  CMO(JBAS +ITR0)= -D3/SQRT(D1270)
                  CMO(JBAS1+ITR0)= -D6/SQRT(D1270)
                  CMO(JBAS2+ITR0)= D24/SQRT(D1270)
                  CMO(JBAS3+ITR0)= -D3/SQRT(D1270)
                  CMO(JBAS4+ITR0)= D24/SQRT(D1270)
                  CMO(JBAS5+ITR0)= -D8/SQRT(D1270)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)= -D1/SQRT(D74)
                  CMO(JBAS2+ITR0)=  D6/SQRT(D74)
                  CMO(JBAS3+ITR0)=  D1/SQRT(D74)
                  CMO(JBAS4+ITR0)= -D6/SQRT(D74)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)= -D1/SQRT(D38)
                  CMO(JBAS1+ITR0)=  D6/SQRT(D38)
                  CMO(JBAS3+ITR0)= -D1/SQRT(D38)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
               ELSE IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XXXY') THEN
                  NORB(ISYM)=NORB(ISYM)+2
                  JBAS1=JBAS+1
261               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'XYYY') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 261
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
262               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'XYZZ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 262
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 269 I=1,2*NBAS(ISYM)
                     CMO(I+ITR0)=D0
269               CONTINUE
                  CMO(JBAS +ITR0)= -D1/SQRT(D8)
                  CMO(JBAS1+ITR0)= -D1/SQRT(D8)
                  CMO(JBAS2+ITR0)=  D6/SQRT(D8)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)= -D1/SQRT(D2)
                  CMO(JBAS1+ITR0)=  D1/SQRT(D2)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
               ELSE IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XXXZ') THEN
                  NORB(ISYM)=NORB(ISYM)+2
                  JBAS1=JBAS+1
271               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'XYYZ') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 271
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
272               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'XZZZ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 272
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 279 I=1,2*NBAS(ISYM)
                     CMO(I+ITR0)=D0
279               CONTINUE
                  CMO(JBAS +ITR0)= -D3/SQRT(D34)
                  CMO(JBAS1+ITR0)= -D3/SQRT(D34)
                  CMO(JBAS2+ITR0)=  D4/SQRT(D34)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)= -D1/SQRT(D10)
                  CMO(JBAS1+ITR0)=  D3/SQRT(D10)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
               ELSE IF(MOLFTUC4(CTMPMO(JBAS+JBAS0)).EQ.'XXYZ') THEN
                  NORB(ISYM)=NORB(ISYM)+2
                  JBAS1=JBAS+1
281               IF(MOLFTUC4(CTMPMO(JBAS1+JBAS0)).NE.'YYYZ') THEN
                     JBAS1=JBAS1+1
                     IF(JBAS1.GT.NBAS(ISYM)) GOTO 901
                     GOTO 281
                  END IF
                  CTMPMO(JBAS1+JBAS0)='KURT'
                  JBAS2=JBAS1+1
282               IF(MOLFTUC4(CTMPMO(JBAS2+JBAS0)).NE.'YZZZ') THEN
                     JBAS2=JBAS2+1
                     IF(JBAS2.GT.NBAS(ISYM)) GOTO 901
                     GOTO 282
                  END IF
                  CTMPMO(JBAS2+JBAS0)='KURT'
               IF(GETCMO) THEN
                  DO 289 I=1,2*NBAS(ISYM)
                     CMO(I+ITR0)=D0
289               CONTINUE
                  CMO(JBAS +ITR0)= -D3/SQRT(D34)
                  CMO(JBAS1+ITR0)= -D3/SQRT(D34)
                  CMO(JBAS2+ITR0)=  D4/SQRT(D34)
                  ITR0=ITR0+NBAS(ISYM)
                  CMO(JBAS +ITR0)=  D3/SQRT(D10)
                  CMO(JBAS1+ITR0)= -D1/SQRT(D10)
                  ITR0=ITR0+NBAS(ISYM)
               END IF
C
C              ... REST OF BASIS FUNCTIONS
C              ... IGNORE DELETED TYPES
C
               ELSE IF(.NOT.(CTMPMO(JBAS+JBAS0).EQ.'KURT')) THEN
                  NORB(ISYM)=NORB(ISYM)+1
               IF(GETCMO) THEN
                  DO 295 I=1,NBAS(ISYM)
                     CMO(I+ITR0)=D0
295               CONTINUE
                  CMO(JBAS+ITR0)=D1
                  ITR0=ITR0+NBAS(ISYM)
               END IF
               END IF
210         CONTINUE
            JBAS0=JBAS0+NBAS(ISYM)
200      CONTINUE
      END IF
      IF(IPRSTAT.GE.9) THEN
         WRITE(LUSTAT,'(A,8I4)') ' NBAS:',(NBAS(I),I=1,NSYM)
         WRITE(LUSTAT,'(A,8I4)') ' NORB:',(NORB(I),I=1,NSYM)
         IF (GETCMO .AND. IPRSTAT.GE.20) THEN
            WRITE(LUSTAT,*) 'Final CMO in CMODEL'
            CALL PRORB(CMO,.FALSE.,LUSTAT)
         END IF
      END IF
      CALL QEXIT('CMODEL')
      RETURN
901   CONTINUE
      WRITE(LUERR,*) 'INFINITE LOOP IN CMODEL DETECTED.'
      WRITE(LUERR,*) 'CONSTRUCTING FOR ',CTMPMO(JBAS+JBAS0)
      CALL QTRACE(LUERR)
      CALL QUIT('ERROR CMODEL, INFINITE LOOP (BASIS TYPE NOT FOUND)')
C     ... end of cmodel.
      END
C  /* Deck molftuc4 */
      CHARACTER*4 FUNCTION MOLFTUC4( TEXT )
C     ( LEFT-ADJUST IN UPPER CASE )
C
C  5-May-1989 -- hjaaj -- force upper case
C  3-Jul-1989 -- hjaaj -- + left adjust and remove blanks
C 12-Mar-1991 -- hjaaj -- removed leading p,d,f,g,... used by HERMIT
C 18-Mar-1993 -- hjaaj -- translate e.g. 'g211' from HERMIT to 'XXYZ'
C
      CHARACTER*4 TEXT, TEXTUC
      CHARACTER*1 CHRUC
      INTEGER     ILCA, ILCZ, UCSHFT, ITEXT
      LOGICAL     FIRST
C
      SAVE        FIRST, ILCA, ILCZ, UCSHFT
      DATA        FIRST /.TRUE./
C
      IF (FIRST) THEN
         ILCA   = ICHAR('a')
         ILCZ   = ICHAR('z')
         UCSHFT = ICHAR('A') - ILCA
         FIRST  = .FALSE.
      END IF
      TEXTUC = '    '
      J = 0
      IF (TEXT(1:1) .EQ. 'g') THEN
C     ... handle Cartesian g-orbtals from HERMIT
C         (spherical g-orbitals are named '5g-4',...,'5g+4')
         READ (TEXT,'(1X,3I1)') L,M,N
         DO 21 I = 1,L
            J = J + 1
            TEXTUC(J:J) = 'X'
   21    CONTINUE
         DO 22 I = 1,M
            J = J + 1
            TEXTUC(J:J) = 'Y'
   22    CONTINUE
         DO 23 I = 1,N
            J = J + 1
            TEXTUC(J:J) = 'Z'
   23    CONTINUE
      ELSE
         DO 100 I = 1,4
            IF (TEXT(I:I) .NE. ' ') THEN
               ITEXT = ICHAR(TEXT(I:I))
               IF (ITEXT .GE. ILCA .AND. ITEXT .LE. ILCZ) THEN
                  CHRUC = CHAR( ITEXT + UCSHFT )
               ELSE
                  CHRUC = TEXT(I:I)
               END IF
C           include only 'X', 'Y', or 'Z' in TEXTUC
               IF (INDEX('XYZ',CHRUC) .NE. 0) THEN
                  J = J + 1
                  TEXTUC(J:J) = CHRUC
               END IF
            END IF
  100    CONTINUE
      END IF
      MOLFTUC4 = TEXTUC
      RETURN
      END
C  /* Deck delmo */
      SUBROUTINE DELMO(CMO,SCRA,LSCRA,THROVL,CMAXMO,DELEMO)
C
C Written 18-Jan-198* by Hans Jorgen Aa. Jensen and Hans Agren
C Revised 26-Aug-1992 by OV+HJAaJ
C Revised 2004-Feb.2005, Jan.2007 by HJAaJ
C
C Purpose:
C  Obtain initial guess for molecular orbitals by
C  diagonalizing the overlap matrix and discarding
C  numerically linar dependent orbitals
C  (defined by eigenvalue of overlap matrix .lt. THROVL).
C
C Input:
C  CMO contains initial molecular orbital coefficients
C      (unit matrix, one where 3s in d etc. are deleted, or ?).
C  THROVL: min overlap eigenvalue
C  CMAXMO: max MO coefficient
C          (Generally CMAXMO \approx sqrt(d1/THROVL)
C           but they can be set independently by user.
C           Therefore we require both criteria fulfilled. /hjaaj jan 07)
C Output:
C  CMO; molecular orbitals
C  DELEMO set true if orbitals have been deleted
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION CMO(*),SCRA(*)
      LOGICAL   DELEMO
      PARAMETER (D1=1.0D0)
C
C Used from common blocks:
C   INFORB : NNBAST,...
C   INFDIM : NNBASM, NBASMA
C   INFPRI : LUW4
C   CBIREA : LMULBS
C   R12INT : LAUXBS
C
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <cbirea.h>
#include <r12int.h>
C
      CALL QENTER('DELMO ')
C
C Core allocation
C
      KOVLP = 1
      KS1T  = KOVLP + NNBAST
      KSCR1 = KS1T  + NNBASM
      KSCR2 = KSCR1 + NBASMA
      IF (KSCR2+NBASMA .GT. LSCRA)
     *   CALL ERRWRK('DELMO',KSCR2+NBASMA,LSCRA)
C
C Read the overlap matrix in AO-basis.
C
      CALL RDONEL('OVERLAP ',.TRUE.,SCRA(KOVLP),NNBAST)
C
      DELEMO = .FALSE.
C
      ICSYM = 1
      DO 200 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         ISSYM = KOVLP + IIBAS(ISYM)
         IF (LAUXBS) THEN
            NORBI = NORB2(ISYM)
         ELSE
            NORB1(ISYM) = NORB(ISYM)
            NORBI = NORB(ISYM)
            ICSYM = ICMO(ISYM) + 1
         END IF
      IF (NORBI.EQ.0) GO TO 200
         IF (LMULBS) THEN
            IF (LAUXBS) THEN
              IF (.NOT. (R12ECO .OR. R12CBS)) THEN
C              Overlap matrix elements that belong to the 
C              primary basis are zeroed (WK/UniKA/04-11-2002).
C              This routine will then delete the corresponding vectors from CMO.
               MBSYM = ISSYM 
               MBNUL = MBAS1(ISYM) * (MBAS1(ISYM) + 1) / 2
               CALL DZERO(SCRA(MBSYM),MBNUL)
               MBSYM = MBSYM + MBNUL 
               DO 202 K = 1, MBAS2(ISYM)
                  CALL DZERO(SCRA(MBSYM),MBAS1(ISYM))
                  MBSYM = MBSYM + MBAS1(ISYM) + K
  202          CONTINUE
              END IF
            ELSE
C              Overlap matrix elements that belong to the 
C              secondary basis are zeroed (WK/UniKA/04-11-2002).
C              This routine will then delete the corresponding vectors from CMO.
               MBSYM = ISSYM + MBAS1(ISYM) * (MBAS1(ISYM) + 1) / 2
               MBNUL = NNBAS(ISYM) - MBSYM + ISSYM
               CALL DZERO(SCRA(MBSYM),MBNUL)
            ENDIF
         END IF
C
C        Update SCR(KS1T) to the new CMO vectors:
C
         CALL UTHU(SCRA(ISSYM),SCRA(KS1T),CMO(ICSYM),SCRA(KSCR1),
     *             NBASI,NORBI)
C
C        Diagonalize overlap matrix ( S * U = U * Seig )
C
         CALL JACO(SCRA(KS1T),CMO(ICSYM),NORBI,NORBI,NBASI,
     *             SCRA(KSCR1),SCRA(KSCR2))
C        CALL JACO (F,V,NB,NMAX,NROWV,BIG,JBIG)
         II = 0
         DO 175 I = 1,NORBI
            II = II + I
            SCRA(KS1T-1+I)=SCRA(KS1T-1+II)
  175    CONTINUE
         CALL ORDER2(CMO(ICSYM),SCRA(KS1T),NORBI,NBASI)
         IF (IPRI6 .GT. 10 .OR. P6FLAG(38)) THEN
            IF (LMULBS .AND. LAUXBS) THEN
              WRITE (LUPRI,'(/A,I5/)')
     &        ' DELMO: eigenvalues of overlap matrix for '//
     &        'auxiliary basis for symmetry',ISYM
            ELSE
              WRITE (LUPRI,'(/A,I5/)')
     &        ' DELMO: eigenvalues of overlap matrix for symmetry',ISYM
            END IF
            WRITE (LUPRI,'(1P,5D15.5)') (SCRA(KS1T-1+I),I=1,NORBI)
         END IF
C
C Delete orbitals with small ( THROVL ) eigenvalues
C
         IDEL = 0
         DO 275 J = 1, NORBI
           EIGVAL = SCRA(KS1T + J - 1)
           JCMO = ICSYM + (J-1)*NBASI
           ICMOMAX = IDAMAX(NBASI,CMO(JCMO),1)
           CMOMAX  = ABS( CMO(JCMO-1+ICMOMAX) ) / SQRT(EIGVAL)
           IF (EIGVAL .LT. THROVL .OR. CMOMAX .GT. CMAXMO) THEN
              IDEL  = (NORBI + 1 - J)
C
              IF (LAUXBS) THEN
C             ... secondary basis
                 JDEL = IDEL
                 IF (.NOT.DELEMO) WRITE (LUW4,'(//A/A/A,1P,D10.2)')
     *           '@ Auxiliary orbitals are deleted during canonical'//
     *           ' orthonormalization',
     *           '@   because of detected numerical linear dependence.',
     *           '@ Eigenvalue threshold for num. lin. dep. =',THROVL
                 WRITE (LUW4,'(/A,I4,A,I5/A/,(1P,5D10.2))')
     *           '@',IDEL,' MO components are deleted in symmetry',ISYM,
     *           ' Overlap eigenvalues of the deleted orbitals:',
     *          (SCRA(KS1T + J - 2 + I),I = 1,IDEL)
                 DELEMO = .TRUE.
                IF (ISYM .LT. NSYM) THEN
                 ICTO  = ICSYM + (NORB2(ISYM)-IDEL)*NBASI
                 ICFROM= ICSYM + NORB2(ISYM)*NBASI
                 DO 260 JSYM = ISYM+1, NSYM
                   DO 250 I = 0, NORB2(JSYM)*NBAS(JSYM) - 1
                     CMO(ICTO + I) = CMO(ICFROM + I)
  250              CONTINUE
                   ICTO = ICTO + NORB2(JSYM)*NBAS(JSYM)     
                   ICFROM = ICFROM + NORB2(JSYM)*NBAS(JSYM)
  260            CONTINUE
                END IF
                NORB2(ISYM) = NORB2(ISYM) - IDEL
              ELSE
C             ... primary basis
                JDEL = IDEL - MBAS2(ISYM)
                IF (JDEL .GT. 0) THEN
                   IF (.NOT.DELEMO) WRITE (LUW4,'(//A/A/A,1P,D10.2)')
     *           '@ Primary orbitals are deleted during canonical'//
     *            ' orthonormalization',
     *           '@   because of detected numerical linear dependence'//
     &            ' or too large coefficients.',
     *           '@ Eigenvalue threshold for num. lin. dep. =',THROVL
                   WRITE (LUW4,'(/A,I4,A,I5/A/,(1P,5D10.2))')
     *           '@',JDEL,' MO components are deleted in symmetry',ISYM,
     *           ' Overlap eigenvalues of the deleted orbitals:',
     *            (SCRA(KS1T + J - 2 + I),I = 1,JDEL)
                   WRITE (LUW4,'(/A,1P,2D10.2)')
     &           ' Max MO coeff. of last deleted MO & limit:',
     &             CMOMAX,CMAXMO

                 DELEMO = .TRUE.
                END IF
                IF (ISYM .LT. NSYM) THEN
                 ICTO  = ICMO(ISYM) + (NORB(ISYM)-IDEL)*NBASI
                 ICFROM= ICMO(ISYM+1)
                 NCMOVE= NCMOT - ICFROM
                 DO 251 I = 1,NCMOVE
                    CMO(ICTO + I) = CMO(ICFROM + I)
  251            CONTINUE
                 DO 261 JSYM=ISYM+1,NSYM
                    ICMO(JSYM) = ICMO(JSYM) - IDEL*NBASI
  261            CONTINUE
                END IF
                NORB(ISYM) = NORB(ISYM) - IDEL
                NCMOT = NCMOT - IDEL*NBASI
              END IF
              GOTO 280
           END IF
  275    CONTINUE
  280    CONTINUE
C
C        Finish canonical orthonormalization by 
C        normalization (C = S**(-0.5)*U = U * Seig**(-0.5) )
C
         IF (LAUXBS) THEN
            NORBI = NORB2(ISYM)
         ELSE
            NORBI = NORB(ISYM)
            NORB1(ISYM) = NORB(ISYM)
         END IF
         DO J = 1, NORBI
           EIGVAL = SCRA(KS1T + J -1)
           EIGVAL = D1/SQRT(EIGVAL)
           CALL DSCAL(NBASI,EIGVAL,CMO(ICSYM+(J-1)*NBASI),1)
         END DO
C
         ICSYM = ICSYM + NORB2(ISYM)*NBASI       
  200 CONTINUE
C
C *** Make sure primary orbitals are orthonormal /hjaaj-dec-04
C     (there may be numerical round-off errors ...)
C     (ORTHO cannot be used for secondary orbitals as it uses
C      NORB(i) and not NORB2(i). )
C
      IF (.NOT. LAUXBS) THEN
         KSMAT = 1
         KWRK  = KSMAT + N2BASX
         LWRK  = LSCRA - KWRK
         CALL ORTHO(CMO,SCRA(KSMAT),SCRA(KWRK),LWRK)
      END IF
C
C *** OUTPUT SECTION
C
      IF (LAUXBS) THEN
        WRITE(LUPRI,'(/A/A)')
     *      '  Canonical primary and secondary basis sets',
     *      '  ISYM       NORB1 MBAS1       NORB2 MBAS2'
        DO ISYM = 1, NSYM
           WRITE(LUPRI,'(I6,6X,2I6,6X,2I6)') 
     *     ISYM, NORB1(ISYM), MBAS1(ISYM), NORB2(ISYM), MBAS2(ISYM)
        END DO
      END IF
C
C *** end of subroutine DELMO
C
      CALL QEXIT('DELMO ')
      RETURN
      END
C  /* Deck reord */
      SUBROUTINE REORD(CMO,CSCR,IORD)
C
C H.AGREN 19.NOV 84
C
C Purpose: Reorder mo:s according to IORD(*) array
C          so new_mo(i) = old_mo(iord(i)).
C
C Input: CMO : MO:s in normal order
C
C Output: CMO : MO:s in IORD order
C
#include <implicit.h>
      DIMENSION CMO(*),CSCR(*),IORD(*)
C
C  INFORB : NCMOT, ICMO(8), NBAS(8), NORB(8)
C
#include <inforb.h>
C
      INEW1 = 1
      DO 10 ISYM = 1,NSYM
         ICMO1 = ICMO(ISYM)+1
         NBASI = NBAS(ISYM)
         IORBI = IORB(ISYM)
         DO 20 I = 1,NORB(ISYM)
            IOLD1 = ICMO1 + NBASI*(IORD(IORBI+I) - (IORBI+1))
            CALL DCOPY(NBASI,CMO(IOLD1),1,CSCR(INEW1),1)
            INEW1 = INEW1 + NBASI
   20    CONTINUE
   10 CONTINUE
C
      CALL DCOPY(NCMOT,CSCR,1,CMO,1)
C
      RETURN
      END
C  /* Deck ortho */
      SUBROUTINE ORTHO(CMO,S,SIN,LSIN)
C
C Original: CASSCF RELEASE 79 11 23
C Revisions:
C  4-Apr-1984 hjaaj
C    Apr-1985 hjaaj (symort, and flag(32) for control)
C  5-Jul-1989 hjaaj (use PRORB to print orbitals)
C
C     OBJECTIVE :
C         ORTHOGONALIZES TRIAL MOLECULAR ORBITALS
C         TRANSFERRED FROM SIRINP VIA THE VECTOR CMO
C
C     SUBROUTINES CALLED:
C         NORM   (GRAM-SCHMIDT ORTHONORMALIZATION)
C         SYMORT (SYMMETRICAL ORTHONORMALIZATION)
C         MOLLAB (OVERLAP MATRIX ON LUONEL)
C
#include <implicit.h>
      DIMENSION CMO(*),S(*),SIN(LSIN)
C
      PARAMETER (D0 = 0.0D0)
C
C Used from common blocks:
C   INFINP : NWARN,CMAXMO,THROVL,?
C   INFORB : NNBAST,NCMOT
C   INFPRI : P4FLAG(*),P6FLAG(*)
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
      LOGICAL PROVLP
      SAVE    PROVLP
      DATA    PROVLP /.TRUE./
C
      CALL QENTER('ORTHO ')
C
C     ***** READ OVERLAP MATRIX S FROM LUONEL *****
C
      CALL RDONEL('OVERLAP ',.TRUE.,S,NNBAST)
C
      IF (PROVLP .AND. P6FLAG(38)) THEN
         PROVLP = .FALSE.
         WRITE(LUPRI,'(/A)')
     *     ' (ORTHO) Overlap between AO basis functions :'
         CALL OUTPKB(S,NBAS,NSYM,-1,LUPRI)
      END IF
C
C     ***** ORTHOGONALIZE ORBITALS SYMMETRY BY SYMMETRY *****
C
      ISTBAS = 0
      ISTS   = 1
      ISTC   = 1
      DO 100 ISYM=1,NSYM
         NORBI=NORB(ISYM)
         NBASI=NBAS(ISYM)
         IF(NORBI.EQ.0) GO TO 101
C
C
         IF (.NOT.FLAG(32)) THEN
            CALL NORM(S(ISTS),CMO(ISTC),NBASI,NORBI,SIN,THROVL,IRETUR)
C           ... error exit with IRETUR.ne.0 if norm**2 < THROVL
C               after Gram-Schmidt orthogonalization to prev. vectors
C           hjaaj may 2000: use THROVL instead of fixed THNORM in CALL NORM
         ELSE
            CALL SYMORT(CMO(ISTC),S(ISTS),NBASI,NORBI,SIN,LSIN,IRETUR)
            IF (IRETUR .NE. 0 .AND. IRETUR .NE. 8888) THEN
C           ... if (not converged and not numerical round-off
C               errors) then ...
              NWARN = NWARN + 1
              WRITE (LUPRI,4020) ISYM,IRETUR
            END IF
            CALL NORM(S(ISTS),CMO(ISTC),NBASI,NORBI,SIN,THROVL,IRETUR)
C           ... 951201: now always call NORM to ensure orthonormality;
C           Cases have been seen where SYMORT have deviation from
C           orthonormality of the order 1.0D-8 because of numerical
C           round-off errors (IRETUR = 8888)
         END IF
         IF (IRETUR.NE.0) THEN
            WRITE (LUPRI,1000)
            WRITE(LUW4,4010) ISYM, IRETUR
            IF (LUPRI.NE.LUW4) WRITE(LUPRI,4010) ISYM, IRETUR
            WRITE(LUERR,4010) ISYM, IRETUR
            CALL PRORB(CMO,.FALSE.,LUPRI)
            GO TO 5000
         END IF
C
  101    ISTBAS = ISTBAS + NBASI
         ISTS   = ISTS   + NBASI*(NBASI+1)/2
         ISTC   = ISTC   + NORBI*NBASI
  100 CONTINUE
C
 4010 FORMAT(/'@*** ORTHO-FATAL ERROR *** Linear dependency in',
     &        ' symmetry =',I3,', CODE =',I4)
 4020 FORMAT(/'@*** ORTHO-WARNING *** Symmetric orthonormalization'//
     &        ' failed for symmetry',I2,
     &       /'@    Will attempt Gram-Schmidt orthonormalization.')
 4030 FORMAT(/'@(ORTHO) This error message will now be suppressed',
     &   ' because it has been given max. no. of times (= 3 times).')
C
      IMAX = IDAMAX(NCMOT,CMO,1)
      CMAX = ABS( CMO(IMAX) )
c     IF (CMAX .GE. CMAXMO) P6FLAG(6) = .TRUE.
C
C     ***** PRINT MOLECULAR ORBITALS ON UNIT LUPRI *****
C
      IF (P6FLAG(6)) THEN
         IF (.NOT.FLAG(32)) THEN
            WRITE (LUPRI,1000)
         ELSE
            WRITE (LUPRI,1010)
         END IF
         CALL PRORB(CMO,.FALSE.,LUPRI)
      END IF
 1000 FORMAT(/' Trial molecular orbitals Gram-Schmidt orthogonalized.')
 1010 FORMAT(/' Trial molecular orbitals symmetrically orthogonalized.')
C
C     ***** PRINT MOLECULAR ORBITALS ON UNIT LUW4 *****
C
      IF (P4FLAG(12) .AND. ( LUW4.NE.LUPRI .OR. .NOT.P6FLAG(6) )) THEN
         IF (.NOT.FLAG(32)) THEN
            WRITE(LUW4,1000)
         ELSE
            WRITE(LUW4,1010)
         END IF
         CALL PRORB(CMO,.TRUE.,LUW4)
      END IF
C
      IF (CMAX.GE.CMAXMO) GO TO 104
C
      CALL QEXIT('ORTHO ')
      RETURN
C
C     ***** ERROR BRANCHES
C     ***** (LINEAR DEPENDENCIES IN ATOMIC BASIS SET)
C
C
  104 CONTINUE
      DO ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         NBASI = NBAS(ISYM)
         NCMOI = NORBI*NBASI
         IMAX = IDAMAX(NCMOI,CMO(ICMO(ISYM)+1),1)
         CMAX = ABS( CMO(ICMO(ISYM)+IMAX) )
      IF (CMAX .GT. CMAXMO) THEN
         IORBI = (IMAX-1)/NBASI + 1
         WRITE(LUERR,3010) CMAX,ISYM,IORBI,NORBI,CMAXMO
         WRITE(LUW4 ,3010) CMAX,ISYM,IORBI,NORBI,CMAXMO
         IF (LUPRI.NE.LUW4)
     &   WRITE(LUPRI,3010) CMAX,ISYM,IORBI,NORBI,CMAXMO
      END IF
      END DO
      WRITE (LUW4,3020)
      IF (LUPRI.NE.LUW4) WRITE(LUPRI,3020)
 3010 FORMAT(/' *** ORTHO-FATAL ERROR ***',
     &    /' Largest molecular orbital coefficient/sym',1P,D20.6,3I5,
     &    /' This number is larger than accepted limit',D20.6)
 3020 FORMAT(
     &    /' Significant loss of accuaracy is probable ',
     &    /' in transformation of two-electron integrals,'
     &   //'*** Program stops here. ***   Your options:',
     &    /' 1) modify basis set,',
     &    /' 2) decrease numerical linear dependency by increasing',
     &     ' ".AO DELETE" limit, or'
     &    /' 3) take the chance and increase',
     &     ' limit with ".CMOMAX" in "*ORBITAL INPUT".')
      GO TO 5000
C
C
 5000 CALL QTRACE(LUERR)
      CALL QUIT('*** ERROR *** FATAL ERROR IN ORTHO')
C
C     end of ORTHO
C
      END
C  /* Deck h1mo */
      SUBROUTINE H1MO(CMO,SH1,SCRA,LSCRA)
C
C Written 15-Apr-1984 by Hans Jorgen Aa. Jensen and Hans Agren
C Last revision 8-Oct-1984 hjaaj
C
C Purpose:
C  Obtain initial guess for molecular orbitals by
C  diagonalizing the one-electron Hamiltonian H1.
C
C Input:
C  H1; the one-electron Hamiltonian in AO-basis.
C
C Output:
C  CMO; molecular orbitals
C
#include <implicit.h>
      DIMENSION CMO(*),SH1(*),SCRA(*)
C
C Used from common blocks:
C   INFINP : FLAG(32)
C
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <infdim.h>
C
      LOGICAL LSAVE1,LSAVE2,LSAVE3
C
      CALL QENTER('H1MO  ')
C
C Step 1: Schmidt orthogonalize atomic orbitals
C         (code: flag(32) false)
C
C         CMO contains initial matrix
C         - either unit matrix or one where some orbitals are deleted
C           (because of e.g. num. lin. dep. or removal of 3s in d).
C
      LSAVE1    = P4FLAG(12)
      LSAVE2    = P6FLAG(6)
      LSAVE3    = FLAG(32)
      P4FLAG(12)= .FALSE.
      P6FLAG(6) = .FALSE.
      FLAG(32)  = .FALSE.
      CALL ORTHO(CMO,SH1,SCRA,LSCRA)
C
C Step 2: Diagonalize one-electron Hamiltonian
C
C
C  Get one-electron Hamiltonian
C
      CALL SIRH1(SH1,SCRA,LSCRA)
C
      KH1T  = 1
      KSCR1 = KH1T + IROW(NBASMA+1)
      KSCR2 = KSCR1 + NBASMA
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         NBASI = NBAS(ISYM)
         ISSYM = IIBAS(ISYM) + 1
         ICSYM = ICMO(ISYM) + 1
C
         CALL UTHU(SH1(ISSYM),SCRA(KH1T),CMO(ICSYM),SCRA(KSCR1),
     *             NBASI,NORBI)
C        CALL UTHU(H,HT,U,WRK,NAO,NMO)
C
         CALL JACO(SCRA(KH1T),CMO(ICSYM),NORBI,NORBI,NBASI,
     *             SCRA(KSCR1),SCRA(KSCR2))
C        CALL JACO (F,V,NB,NMAX,NROWV,BIG,JBIG)
         II = 0
         DO 175 I = 1,NORBI
            II = II + I
            SCRA(KH1T-1+I)=SCRA(KH1T-1+II)
  175    CONTINUE
         CALL ORDER (CMO(ICSYM),SCRA(KH1T),NORBI,NBASI)
C
  200 CONTINUE
C
C Step 3: Reorthogonalize new mo's using Gram-Schmidt
C
      CALL ORTHO(CMO,SH1,SCRA(1),LSCRA)
      P4FLAG(12)= LSAVE1
      P6FLAG(6) = LSAVE2
      FLAG(32)  = LSAVE3
C
C *** end of subroutine H1MO
C
      CALL QEXIT('H1MO  ')
      RETURN
      END
C  /* Deck h1occ */
      SUBROUTINE H1OCC(CMO,WRK,KFRSAV,LFRSAV)
C
C Automatic determination of initial HF-occupation
C from diagonal elements of one-electron Hamiltonian H1.
C
C Written 25-Aug-1995 by Hans Jorgen Aa. Jensen
C
C Based in part of modifications originally made in H1MO
C by K.Ruud-May 1995 (H1MO now restored to previous content).
C
C Input:
C  CMO; molecular orbitals
C
#include <implicit.h>
      DIMENSION CMO(*),WRK(*)
C
C Used from common blocks:
C  SCBRHF : IOPRHF
C  INFORB : NRHF(), NNBAST,NNORBT,NBAST,...
C  INFIND : IROW()
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <scbrhf.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      CALL QENTER('H1OCC ')
C
      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET('REAL',KHMO,NNORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KHAO,NNBAST,WRK,KFREE,LFREE)
C
C ***** retrieve atomic ONE ELECTRON HAMILTONIAN matrix
C
      CALL SIRH1(WRK(KHAO),WRK(KFREE),LFREE)
C
C     Transform ONE ELECTRON HAMILTONIAN to MO basis
C
      CALL MEMGET('REAL',KWRK,N2BASX,WRK,KFREE,LFREE)
      CALL UTHUB(WRK(KHAO),WRK(KHMO),CMO,WRK(KWRK),NSYM,NBAS,NORB)
C
C     ***** Test output of ONE ELECTRON HAMILTONIAN  matrices *****
C
      IF (IPRI6 .GE. 15) THEN
         WRITE(LUPRI,1000)
         CALL OUTPKB(WRK(KHAO),NBAS,NSYM,-1,LUPRI)
         WRITE(LUPRI,1200)
         CALL OUTPKB(WRK(KHMO),NORB,NSYM,-1,LUPRI)
      END IF
      CALL MEMREL('H1OCC.1',WRK,KFRSAV,KHAO,KFREE,LFREE)
C
 1000 FORMAT(/' H1OCC: TEST OF ONE ELECTRON HAMILTONIAN (AO-basis)')
 1200 FORMAT(/' H1OCC: TEST OF ONE ELECTRON HAMILTONIAN (MO-basis)')
C
C     Extract H1 diagonal in WRK(KH1D) and
C     orbital symmetries in WRK(KSMO)
C
      CALL MEMGET('REAL',KH1D,NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMO,NORBT,WRK,KFREE,LFREE)
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         JHMO = KHMO-1 + IIORB(ISYM)
         JH1D = KH1D-1 + IORB(ISYM)
         JSMO = KSMO-1 + IORB(ISYM)
         DO 175 I = 1,NORBI
            WRK(JH1D+I) = WRK(JHMO+IROW(I+1))
            WRK(JSMO+I) = ISYM
  175    CONTINUE
  200 CONTINUE
C
C     Sort according to diagonal element and
C     determine HF occupation.
C
      CALL ORDER(WRK(KSMO),WRK(KH1D),NORBT,1)
      CALL IZERO(NRHF,8)
      MOCC = NRHFEL/2
      DO 20 I = 1, MOCC
         ISYM = NINT(WRK(KSMO-1+I))
         NRHF(ISYM) = NRHF(ISYM) + 1
   20 CONTINUE
      IF (2*MOCC .NE. NRHFEL) IOPRHF = NINT(WRK(KSMO+MOCC))
C
C *** end of subroutine H1OCC
C
      CALL QEXIT('H1OCC ')
      RETURN
      END
C  /* Deck fcmo */
      SUBROUTINE FCMO(MXFCK,CMO,FC,SCRA,LSCRA)
C
C Written 4-May-1984 by Hans Jorgen Aa. Jensen
C Revisions:
C   8-Oct-1984 hjaaj
C   1-Nov-1984 hjaaj (use NORB(*), not NBAS(*) for FC)
C         1985 hjaaj (use NRHF(*) for RHF occupation)
C   4-Mar-1997 tsaue include screening
C
C Purpose:
C  Do MXFCK restricted Roothaan-Hartree-Fock iterations.
C
C  -- idea: grand-canonical Hartree-Fock can be specified
C           by NASH(*) = NRHF(*), NASHT = NRHFT,
C           DV(ij) = delta(ij) occ(i); but NISHT = 0 and NISH(*) = 0.
C           Then GC Fock matrix is FC + FV = h1 + FV. (860117)
C
C Input:
C  CMO; initial molecular orbitals used to build Fock matrix,
C       assumed to be orthonormal.
C  MXFCK; maximum number of Fock iterations (always one iteration)
C
C Output:
C  CMO; molecular orbitals diagonalizing Fock matrix
C
C Scratch:
C  FC; the inactive Fock matrix and scratch area for overlap matrix
C  SCRA; general scratch area
C
#include <implicit.h>
      DIMENSION CMO(*),FC(*),SCRA(*)
C
C
      PARAMETER (D0=0.0D0, EMYCNV = 1.D-4)
#include <dummy.h>
C
C Used from common blocks:
C  INFINP : 
C  INFOPT : EPOT,EMY,EMCSCF
C  SCBRHF : NFRRHF(*)
C  INFIND : ...,ISSMO(*),?
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infopt.h>
#include <inforb.h>
#include <scbrhf.h>
#include <infind.h>
#include <infpri.h>
#include <gnrinf.h>
C
      LOGICAL LSAVE4,LSAVE6
C
      CALL QENTER('FCMO  ')
      WRITE (LUPRI,'(A//A/A//A,I5/A,8I5)') '1',
     *   ' Restricted Roothaan-Hartree-Fock iterations',
     *   ' -------------------------------------------',
     *   ' Number of electrons :',2*NRHFT,
     *   ' Orbital occupations :',(NRHF(I),I=1,NSYM)
C
      NFRHFT = ISUM(NSYM,NFRRHF,1)
      IF (NFRHFT .GT. 0) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//A/A/)')
     *     '@ WARNING from FCMO: ".FROZEN" is not implemented'//
     *     ' for Fock iterations,',
     *     '@ Fock iterations abandoned.'
         GO TO 9999
      END IF
C
      LSAVE4 = P4FLAG(12)
      LSAVE6 = P6FLAG(6)
      P4FLAG(12)= .FALSE.
      P6FLAG(6) = .FALSE.
      DO 5 ISYM = 1,NSYM
         ISWAP      = NRHF(ISYM)
         NRHF(ISYM) = NISH(ISYM)
         NISH(ISYM) = ISWAP
    5 CONTINUE
      ISWAP = NRHFT
      NRHFT = NISHT
      NISHT = ISWAP
      IEXIT  = 0
      ITFCK  = 0
      EMY    = D0
  100 CONTINUE
         ITFCK  = ITFCK + 1
         EMYOLD = EMY
C
C        Step 1: Construct inactive Fock matrix
C
         CALL FCKMAT(.TRUE.,DUMMY,CMO,EMY,FC,DUMMY,SCRA,LSCRA)
C        CALL FCKMAT(ONLYFC,DV,CMO,EMY,FC,FV,WRK,LFREE)
         IF (SUPSYM) CALL AVEFCK(FC)
         EMCSCF = EPOT + EMY
      IF (IEXIT .EQ. 1) GO TO 500
         IF (IPRSIR .GT. 0) WRITE (LUPRI,1730) ITFCK,EMY,EMCSCF
C
C        Step 2: Diagonalize inactive Fock-matrix:
C
         DO 200 ISYM = 1,NSYM
            NORBI = NORB(ISYM)
         IF (NORBI.EQ.0) GO TO 200
            IORBI = IORB(ISYM)
            NBASI = NBAS(ISYM)
            ISSYM = IIORB(ISYM) + 1
            ICSYM = ICMO(ISYM) + 1
C
            CALL JACO(FC(ISSYM),CMO(ICSYM),NORBI,NORBI,NBASI,
     *                SCRA,SCRA(1+NORBI))
C           CALL JACO (F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
            II = ISSYM - 1
            DO 175 I=1,NORBI
               II = II + I
               SCRA(I)=FC(II)
  175       CONTINUE
            CALL ORDRSS(CMO(ICSYM),SCRA,ISSMO(IORBI+1),NORBI,NBASI)
            IF (IPRSIR .GE. 5) THEN
               IF (IPRSIR .GE. 12) THEN
                  IEND = NORBI
               ELSE
                  IEND = MIN(NORBI,NISH(ISYM)+2)
               END IF
               WRITE (LUPRI,1740) ISYM
               WRITE (LUPRI,1750) (SCRA(I),I=1,IEND)
            END IF
  200    CONTINUE
         IF (SUPSYM) CALL AVEORD()
C        ... remake ISSORD() as ISSMO() may have changed in ORDRSS
C
 1730    FORMAT(//' Fock iteration number',I3,
     *             '; inactive energy:',F25.12,
     *          /T28,'total    energy:',F25.12)
 1740    FORMAT(/' Fock eigenvalues symmetry',I2/)
 1750    FORMAT(4X,5F15.8)
C
C        Step 3: Reorthogonalize new mo's
C
         CALL ORTHO(CMO,FC,SCRA(1),LSCRA)
C
C        Another Fock iteration?
C
ckr         IF (ITFCK .GT. 1 .AND. EMY .GT. EMYOLD) THEN
C
C           Stop because of oscillation, i.e. energy has gone up.
C           One more iteration has been taken. This should lead to
C           lower energy if the oscillation is typical.
C
ckr            WRITE (LUPRI,'(//2A/)') ' *** Energy is oscillating,',
ckr     *         ' exit from Roothaan-Hartree-Fock iterations.'
ckr         ELSE
            IF (ITFCK .GT. 1 .AND. abs(EMYOLD-EMY) .LE. EMYCNV) THEN
               WRITE (LUPRI,'(//2A,1P,D10.2/)')
     *            ' *** Roothaan-Hartree-Fock energy difference',
     *            ' converged to',(EMYOLD-EMY)
            ELSE IF (ITFCK .LT. MXFCK) THEN
               GO TO 100
C     ^-----------------
            ELSE
               WRITE (LUPRI,'(//A/)')
     &            ' *** Max. number of iterations reached.'
            END IF
ckr         END IF
         IEXIT = 1
         GO TO 100
C        ... go calculate final energy
C
C
C
  500 CONTINUE
      WRITE (LUPRI,1730) ITFCK,EMY,EMCSCF
      P4FLAG(12)= LSAVE4
      P6FLAG(6) = LSAVE6
      DO 800 ISYM = 1,NSYM
         ISWAP      = NRHF(ISYM)
         NRHF(ISYM) = NISH(ISYM)
         NISH(ISYM) = ISWAP
  800 CONTINUE
      ISWAP = NRHFT
      NRHFT = NISHT
      NISHT = ISWAP
C
C *** end of subroutine FCMO
C
 9999 CALL QEXIT('FCMO  ')
      RETURN
      END
C  /* Deck fcvirt */
      SUBROUTINE FCVIRT(CMO,WRK,LFREE)
C
C  2-Oct-1986 Poul Joergensen
C  Revised 28-Aug-1995 hjaaj
C
C  Purpose : To use one electron hamiltonian or modified FC
C            to modify virtual Hartree-Fock orbitals for CI/MC
C
C  Reference for .FC MVO: C.W. Bauschlicher, JCP 72 (1980) 880.
C
#include <implicit.h>
C
      DIMENSION CMO(*), WRK(LFREE)
C
      PARAMETER ( D0 = 0.0D0 )
#include <dummy.h>
C
C   INFORB : NSYM, NNBAST, NNORBT, ...
C   SCBRHF : IOPRHF, NMVO(), NMVOT
C   INFIND : IROW(*)
C   INFPRI : IPRI6
C   INFIND : SUPSYM
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <scbrhf.h>
#include <infind.h>
#include <infpri.h>
#include <infinp.h>
C
      DIMENSION MISH(8), MRHF(8)
C
      CALL QENTER('FCVIRT')
C
      CALL ICOPY(8,NRHF,1,MRHF,1)
      IF (IOPRHF .GT. 0) MRHF(IOPRHF) = MRHF(IOPRHF) + 1
C        ... if NRHF() was not defined in *RHF CALC
C            it will contain NISH() from *WAVE FUNC
C
      IF (NMVOT .EQ. 0) THEN
         WRITE (LUPRI,'(//A/A)')
     &   ' >>> Modified virtual orbitals generated by diagonalization',
     &   ' >>> of virtual block of one-electron Hamiltonian.'
      ELSE
         WRITE (LUPRI,'(//A/A/A,8I4)')
     &   ' >>> Modified virtual orbitals generated by diagonalization',
     &   ' >>> of virtual block of core Fock matrix.',
     &   ' Number of core orbitals in each symmetry    : ',
     &   (NMVO(I),I=1,NSYM)
      END IF
      WRITE (LUPRI,'(A,8I4/)')
     &   ' Number of occupied orbitals in each symmetry: ',
     &   (MRHF(I),I=1,NSYM)
C
      KFC_MVO  = 1
      KWRK = KFC_MVO  + NNORBT
      LNEED= KWRK + 2*NBAST
      LWRK = LFREE - KWRK
      IF (LNEED .GT. LFREE) CALL ERRWRK('FCVIRT',LNEED,LFREE)
C
C ***** Calculate a Fock matrix "FC_MVO" based on NMVO(*) occupied orbitals
C       instead fo NISH(*).
C
      CALL ICOPY(8,NISH,1,MISH,1)
      CALL ICOPY(8,NMVO,1,NISH,1)
      MISHT = NISHT
      NISHT = NMVOT
      CALL FCKMAT(.TRUE.,DUMMY,CMO,EMCMY,WRK(KFC_MVO),DUMMY,
     &             WRK(KWRK),LWRK)
C     CALL FCKMAT(ONLYFC,DV,CMO,EMCMY,FC,FV,WRK,LFRSAV)
      CALL ICOPY(8,MISH,1,NISH,1)
      NISHT = MISHT
C
C ***** Zero all elements in FC_MVO matrix except virtual Hartree-Fock block
C
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         IORBI = IORB(ISYM)
         NBASI = NBAS(ISYM)
         NOCCI = MRHF(ISYM)
         JXFC  = KFC_MVO + IIORB(ISYM)
         JCMO =  1 + ICMO(ISYM)
C
C
         DO 167 I = 1 , NORBI
            JSTA = JXFC + IROW(I)
            JEND = JSTA - 1 + MIN(I,NOCCI)
            DO 164 J = JSTA,JEND
               WRK(J) = D0
  164       CONTINUE
  167    CONTINUE
C
C
         KBIG  = KWRK
         KJBIG = KBIG + NORBI
         CALL JACO(WRK(JXFC),CMO(JCMO),NORBI,NORBI,NBASI,
     *             WRK(KBIG),WRK(KJBIG))
C        CALL JACO (F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
         DO 175 I = 1,NORBI
            WRK(I) = WRK(JXFC-1 + IROW(I+1))
  175    CONTINUE
C
         NSSHI = NORBI - NOCCI
         IF (NSSHI .GT. 0) THEN
            JCMO = JCMO + NOCCI*NBASI
C           order virtual HARTREE-FOCK orbitals
            CALL ORDRSS(CMO(JCMO),WRK(1+NOCCI),
     &                  ISSMO(IORBI+NOCCI+1),NSSHI,NBASI)
         END IF
C
         IF (IPRI6 .GE. 10)
     *   WRITE (LUPRI,'(/A/A,I3,//,(5(I3,F12.6)))')
     *   ' Super symmetry and eigenvalues of virtual one-electron'//
     *   ' Hamiltonian',' Symmetry',ISYM,
     *   (ISSMO(IORBI+I),WRK(I),I=(NOCCI+1),NORBI)
  200 CONTINUE
      IF (SUPSYM) CALL AVEORD()
C     ... remake ISSORD() as ISSMO() may have changed in ORDRSS
C
      CALL QEXIT('FCVIRT')
      RETURN
      END
C  /* Deck prorb */
      SUBROUTINE PRORB(CMO,PROCC,IOUT)
C
C Written 11-Dec-1983 by ha/hjaaj
C Revised  8-Jul-1992 hjaaj
C
C Purpose:
C  Print molecular orbital coefficients on unit IOUT.
C
C Input:
C  CMO:    MO orbital coefficients (symmetry blocked)
C  PROCC:  print only occupied orbitals
C  IOUT:   output file unit
C
#include <implicit.h>
      DIMENSION CMO(*)
      LOGICAL   PROCC
C
C Used from common blocks:
C   INFINP : CENT,TYPE,SUPSYM,?
C   INFORB : NSYM,...
C   INFIND : ISSMO(),
C   CBIREA : LMULBS
C
#include <maxorb.h>
#include <maxash.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <cbirea.h>
#include <r12int.h>
C
C *** FORMAT statements
C
 2400 FORMAT(/' Molecular orbitals for symmetry species',I3,
     &       /' ------------------------------------------')
 2600 FORMAT(/'    Orbital ',7I9)
 2620 FORMAT( '    Sup.sym.',7I9)
 2800 FORMAT(I4,1X,A4,':',A4,7F9.4)
C
C *** Print molecular orbitals
C
      IF (PROCC) THEN
         THRPRI = 1.0D-2
         WRITE(IOUT,'(/A,F6.4,A)')
     &   ' (Only coefficients >',THRPRI,' are printed.)'
      ELSE
         THRPRI = 1.0D-4
      END IF
C
      ISTBAS = 0
      DO 400 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         MBASI = NBAS(ISYM)
         IF (LMULBS) MBASI = MBAS1(ISYM)
C
         IF (PROCC) THEN
            IF (NBAST .LE. 25) THEN
C           print all occupied and 2 lowest unoccupied this symmetry
               IEND = 0
            ELSE
C           print 5 highest occupied and 2 lowest unoccupied this symmetry
C           (output becomes too big otherwise ... /hjaaj Dec08)
               IEND  = MAX(0,NOCC(ISYM)-5)
            END IF
            NENDI = NOCC(ISYM)
     &            + MIN(2,NORB(ISYM)-NOCC(ISYM) - NRXR12(ISYM))
         ELSE
            IEND  = 0
            NENDI = NORB(ISYM)
         END IF
      IF (NENDI.EQ.0) GO TO 300
         WRITE(IOUT,2400) ISYM
C
         ICMOI = ICMO(ISYM)
         ISTORB = IORB(ISYM)
  100       IST =IEND+1
            ISTMO=IEND*NBASI+ICMOI
            IEND=IEND+7
            IF(IEND.GT.NENDI) IEND=NENDI
            IEMO=NBASI*(IEND-1)+ICMOI
            WRITE(IOUT,2600) (I,I=IST,IEND)
            IF (SUPSYM) WRITE(IOUT,2620) (ISSMO(ISTORB+I),I=IST,IEND)
            DO 200 I=1,MBASI
               JSMO=ISTMO+I
               JEMO=IEMO+I
               JJ = 0
               DO J=JSMO,JEMO,NBASI
                  IF ( ABS(CMO(J)) .GE. THRPRI ) JJ = 1
               END DO
               IF (JJ .EQ. 1) WRITE(IOUT,2800) I,CENT(I+ISTBAS),
     *            TYPE(I+ISTBAS),(CMO(J),J=JSMO,JEMO,NBASI)
  200       CONTINUE
         IF (IEND.NE.NENDI) GO TO 100
C
  300 CONTINUE
        ISTBAS = ISTBAS + NBASI
  400 CONTINUE
C
C *** End of subroutine PRORB
C
      RETURN
      END
C  /* Deck punmo */
        SUBROUTINE PUNMO(IPCTL,CMO,OCC)
C
C 25-May-1985 hjaaj
C l.r. 910715-hjaaj (added formats 7010-7050)
C
C Purpose:
C   punch mo coefficients and, conditionally, occupation numbers
C   to LUPNCH
C
C Control:
C   IPCTL .eq. 1: also punch occupation numbers
C   else        : do not punch occupation numbers
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*),OCC(*)
      PARAMETER (D0 = 0.0D0)
C
C Used from common blocks:
C   INFORB : NSYM,NORB(*),NBAS(*)
C
#include <inforb.h>
C
      CHARACTER*8 LBLDAT(2)
C
      CALL QENTER('PUNMO')
      CALL GETDAT(LBLDAT(1),LBLDAT(2))
      LUPNCH = -1
      CALL GPOPEN(LUPNCH,'DALTON.MOPUN','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND LUPNCH
      IF (IPCTL.EQ.1) THEN
         WRITE (LUPNCH,'(A,A8,2X,A8,A)')
     *      '**NATORB   ( punched by SIRIUS ',LBLDAT,' )'
      ELSE
         WRITE (LUPNCH,'(A,A8,2X,A8,A)')
     *      '**MOLORB   ( punched by SIRIUS ',LBLDAT,' )'
      END IF
C
      IEND =0
      DO 700 ISYM=1,NSYM
         NORBI=NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 700
         NBASI=NBAS(ISYM)
         DO 705 NI=1,NORBI
            IST=IEND+1
            IEND=IEND+NBASI
            IMOMXV = IDAMAX(NBASI,CMO(IST),1)
            CMOMXV = ABS(CMO(IST-1+IMOMXV))
            IF (CMOMXV .LT. 1.0D2) THEN
               WRITE (LUPNCH,7000) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D3) THEN
               WRITE (LUPNCH,7010) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D4) THEN
               WRITE (LUPNCH,7020) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D5) THEN
               WRITE (LUPNCH,7030) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D6) THEN
               WRITE (LUPNCH,7040) (CMO(NP),NP=IST,IEND)
            ELSE
               WRITE (LUPNCH,7050) (CMO(NP),NP=IST,IEND)
            END IF
  705    CONTINUE
         NDELI = NBASI-NORBI
      IF (NDELI.EQ.0) GO TO 700
         DO 706 NI=1,NDELI
            WRITE (LUPNCH,7000) (D0,NP=1,NBASI)
  706    CONTINUE
  700 CONTINUE
      CALL GPCLOSE(LUPNCH,'KEEP')
C
      IF (IPCTL.EQ.1) THEN
         WRITE (LUPNCH,'(A,A8,2X,A8,A)')
     *   '**NATOCC   ( punched by SIRIUS ',LBLDAT,' )'
         WRITE(LUPNCH,7000) (OCC(NP),NP=1,NORBT)
      END IF
C
      CALL QEXIT('PUNMO')
      RETURN
C
 7000 FORMAT(4F18.14)
 7010 FORMAT(4F18.13)
 7020 FORMAT(4F18.12)
 7030 FORMAT(4F18.11)
 7040 FORMAT(4F18.10)
 7050 FORMAT(4G18.11)
C
C end of PUNMO
C
      END
C --- end of sirius/sircmo.F ---
