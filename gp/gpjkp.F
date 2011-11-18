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
#include <single.h>
C /* Deck erf */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION ERF(X)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Mar. 2003
C
C     Purpose : Calculate value of error-function in point X
C               Based on code from numerical recipies.
C
C*****************************************************************************
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
      IF(X.LT.D0)THEN
        ERF = -GAMMP2(DP5,X**2)
      ELSE
        ERF =  GAMMP2(DP5,X**2)
      ENDIF
      RETURN
      END
C
C*****************************************************************************
C
      FUNCTION GAMMP2(A,X)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      IF (X .LT. D0 .OR. A.LE. D0) CALL QUIT('BAD ARGUMENTS IN GAMMP2')
      IF (X.LT.(A+D1)) THEN
        CALL GSER2(GAMSER,A,X)
        GAMMP2=GAMSER
      ELSE
        CALL GCF2(GAMMCF,A,X)
        GAMMP2=D1-GAMMCF
      ENDIF
      RETURN
      END
C
C*****************************************************************************
C
      SUBROUTINE GSER2(GAMSER,A,X)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ITMAX=100, EPS=3.0D-7, D0=0.0D0, D1 = 1.0D0)
      GLN=GAMMLN2(A)
      IF(X.LE.D0)THEN
        IF(X.LT.D0) 
     >    WRITE(LUPRI,'(//3X//,A)') 'WARNING IN ERF : X < 0 IN GSER'
        GAMSER=D0
        RETURN
      ENDIF
      AP=A
      SUM=D1/A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+D1
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GOTO 1
11    CONTINUE
      IF(X.LT.D0) 
     >  WRITE(LUPRI,'(//3X//,A)') 'A TOO LARGE, ITMAX TOO SMALL IN GSER'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C
C*****************************************************************************
C
      SUBROUTINE GCF2(GAMMCF,A,X)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ITMAX=100,EPS=3.0D-7,FPMIN=1.0D-30, D1 = 1.0D0,
     >           D2 = 2.0D0)
      GLN=GAMMLN(A)
      B=X+D1-A
      C=D1/FPMIN
      D=D1/B
      H=D
      DO 11 I=1,ITMAX
        AN=-I*(I-A)
        B=B+D2
        D=AN*D+B
        IF(ABS(D).LT.FPMIN)D=FPMIN
        C=B+AN/C
        IF(ABS(C).LT.FPMIN)C=FPMIN
        D=D1/D
        DEL=D*C
        H=H*DEL
        IF(ABS(DEL-D1).LT.EPS)GOTO 1
11    CONTINUE
      WRITE(LUPRI,'(//3X//,A)') 'A TOO LARGE, ITMAX TOO SMALL IN GCF2'
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*H
      RETURN
      END
C
C*****************************************************************************
C
      FUNCTION GAMMLN2(XX)
#include <implicit.h>
      DOUBLE PRECISION STP,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,
     *24.01409824083091D0,-1.231739572450155D0,.1208650973866179D-2,
     *-.5395239384953D-5,2.5066282746310005D0/
      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO 11 J=1,6
        Y=Y+1.0D0
        SER=SER+COF(J)/Y
11    CONTINUE
      GAMMLN2=TMP+LOG(STP*SER/X)
      RETURN
      END
C*****************************************************************************
      FUNCTION DERF(X)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Mar. 2003
C
C     Purpose : Calculate value of error-function in point X
C
C               Based on f90-code from : 
C               Naval Surface Warfare Center Mathematical Library
C               (http://www.math.iastate.edu/burkardt/f_src/nswc/nswc.html)
C
C*****************************************************************************
#include <implicit.h>
      PARAMETER (C = .564189583547756D0, ONE = 1.0D0, HALF = 0.5D0,
     >           ZERO = 0.0D0)
      DOUBLE PRECISION A(5), B(3), P(8), Q(8), R(5), S(4)
      SAVE A,B,P,Q,R,S
      DATA A / .771058495001320D-04, -.133733772997339D-02,
     >         .323076579225834D-01,  .479137145607681D-01,
     >         .128379167095513D+00 /
      DATA B / .301048631703895D-02,  .538971687740286D-01,
     >         .375795757275549D+00 /
      DATA P / -1.36864857382717D-07, 5.64195517478974D-01,
     >          7.21175825088309D+00, 4.31622272220567D+01,
     >          1.52989285046940D+02, 3.39320816734344D+02,
     >          4.51918953711873D+02, 3.00459261020162D+02 /
      DATA Q /  1.00000000000000D+00, 1.27827273196294D+01,
     >          7.70001529352295D+01, 2.77585444743988D+02,
     >          6.38980264465631D+02, 9.31354094850610D+02,
     >          7.90950925327898D+02, 3.00459260956983D+02 /
      DATA R /  2.10144126479064D+00, 2.62370141675169D+01,
     >          2.13688200555087D+01, 4.65807828718470D+00,
     >          2.82094791773523D-01 /
      DATA S /  9.41537750555460D+01, 1.87114811799590D+02,
     >          9.90191814623914D+01, 1.80124575948747D+01 /
C
C
      AX = ABS(X)
C
      IF (AX .LE. HALF) THEN
        T = X*X
        TOP = ((((A(1)*T + A(2))*T + A(3))*T + A(4))*T + A(5)) + ONE
        BOT = ((B(1)*T + B(2))*T + B(3))*T + ONE
        FN_VAL = X*(TOP/BOT)
      ELSE IF (AX .LE. 4.0D0) THEN
        TOP = ((((((P(1)*AX + P(2))*AX + P(3))*AX + P(4))*AX + P(5))*AX
     >        + P(6))*AX + P(7))*AX + P(8)
        BOT = ((((((Q(1)*AX + Q(2))*AX + Q(3))*AX + Q(4))*AX + Q(5))*AX
     >        + Q(6))*AX + Q(7))*AX + Q(8)
        FN_VAL = HALF + (HALF - DEXP(-X*X)*TOP/BOT)
        IF (X .LT. ZERO) FN_VAL = -FN_VAL
      ELSE IF (AX .LT. 5.8D0) THEN
        X2 = X*X
        T = ONE / X2
        TOP = (((R(1)*T + R(2))*T + R(3))*T + R(4))*T + R(5)
        BOT = (((S(1)*T + S(2))*T + S(3))*T + S(4))*T + ONE
        FN_VAL = (C - TOP/(X2*BOT)) / AX
        FN_VAL = HALF + (HALF - DEXP(-X2)*FN_VAL)
        IF (X .LT. ZERO) FN_VAL = -FN_VAL
      ELSE
        FN_VAL = SIGN(ONE, X)
      END IF
C      
      DERF = FN_VAL
      RETURN
      END 
C*****************************************************************************
      FUNCTION DERFC(X)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Mar. 2003
C
C     Purpose : Calculate value of complimentary error-function in point X
C
C               Based on f90-code from : 
C               Naval Surface Warfare Center Mathematical Library
C               (http://www.math.iastate.edu/burkardt/f_src/nswc/nswc.html)
C
C*****************************************************************************
#include <implicit.h>
      DERFC = 1.0D0 - DERF(X)
      RETURN
      END 
C*****************************************************************************
