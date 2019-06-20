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
      SUBROUTINE MXMA(A,ICA,IRA,B,ICB,IRB,C,ICC,IRC,NROW,NSUM,NCOL)
C
C    This routine does the matrix multiply C = A.B. We replace all
C    call to DGEMM which do only that by this routines, since
C    this one is cheaper.
C
#include "implicit.h"
      PARAMETER(D0=0.0D0)
#include "priunit.h"
      DIMENSION A(*),B(*),C(*)
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEX=1+ICC*IROW+IRC*ICOL
         C(INDEX)=SUM
100   CONTINUE
      RETURN
      END
      SUBROUTINE MXMAA(A,ICA,IRA, B,ICB,IRB, C,ICC,IRC, NROW,NSUM,NCOL)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
 
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEX=1+ICC*IROW+IRC*ICOL
         C(INDEX)=C(INDEX)+SUM
100   CONTINUE
      RETURN
      END
      SUBROUTINE MXMAM(A,ICA,IRA, B,ICB,IRB, C,ICC,IRC, NROW,NSUM,NCOL)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
 
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEX=1+ICC*IROW+IRC*ICOL
         C(INDEX)=C(INDEX)-SUM
100   CONTINUE
      RETURN
      END
      SUBROUTINE MXMT(A,ICA,IRA, B,ICB,IRB, C, NROW,NSUM)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
      IND=0
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,IROW
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         IND=IND+1
         C(IND)=SUM
100   CONTINUE
      RETURN
      END
