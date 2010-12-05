C  /* Deck so_secgp */
      SUBROUTINE SO_SECGP(GPVC1,LGPVC1,PRP1,LPRP1,ISYMTR,DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C
C     PURPOSE: Calculate second order contribution to 1p-1h part of
C              gradient property vectors.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <soppinf.h>
C
      PARAMETER   (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION   GPVC1(LGPVC1), PRP1(LPRP1)
      DIMENSION   DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_SECGP')
C
      SQ2 = DSQRT(TWO)
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMTR)
         ISYMJ = ISYMI
C
         KOFF1 = IFCRHF(ISYMA,ISYMJ) + NRHF(ISYMA) + 1
         KOFF2 = IIJDEN(ISYMJ,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1 
C
         NVIRA = MAX(NVIR(ISYMA),1)
         NORBA = MAX(NORB(ISYMA),1)
         NRHFJ = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              SQ2,PRP1(KOFF1),NORBA,DENSIJ(KOFF2),NRHFJ,
     &              ONE,GPVC1(KOFF3),NVIRA)
C
         ISYMJ = ISYMA
C
         KOFF4 = IAIDEN(ISYMA,ISYMJ) + 1
         KOFF5 = IFCRHF(ISYMJ,ISYMI) + 1
C
         NORBJ = MAX(NORB(ISYMJ),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              -SQ2,DENSAI(KOFF4),NVIRA,PRP1(KOFF5),NORBJ,
     &              ONE,GPVC1(KOFF3),NVIRA)
C
         ISYMB = ISYMI
C
         KOFF6 = IFCVIR(ISYMA,ISYMB) + NRHF(ISYMA) + 1
         KOFF7 = IAIDEN(ISYMB,ISYMI) + 1
C
         NVIRB = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              SQ2,PRP1(KOFF6),NORBA,DENSAI(KOFF7),NVIRB,
     &              ONE,GPVC1(KOFF3),NVIRA)
C
         ISYMB = ISYMA
C
         KOFF8 = IABDEN(ISYMA,ISYMB) + 1
         KOFF9 = IFCRHF(ISYMB,ISYMI) + NRHF(ISYMB) + 1
C
         NORBB = MAX(NORB(ISYMB),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              -SQ2,DENSAB(KOFF8),NVIRA,PRP1(KOFF9),NORBB,
     &              ONE,GPVC1(KOFF3),NVIRA)
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_SECGP')
C   
      RETURN
      END
