C
C$Id: gamcom.h,v 1.1.1.1 2001-02-08 13:33:25 hjj Exp $
C
      PARAMETER (MAXJ = 4*(MXQN - 1) + 2)
      COMMON /GAMCOM/ WVAL, FJW(0:MAXJ), TABFJW(121*(MAXJ + 7)), JMAX0
