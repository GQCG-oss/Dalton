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
      SUBROUTINE CDENS1(ISYMB,CMO,ZYMB,DXCAO,WRK,LWRK)
C
C FEB-26 1995
C
C Purpose:
C To compute and add together 2 density matrices to give
C D-m = DXCAO = (see notes)
C
C needed in DDCRPA (closed shell) for the construction of
C Fock matrix
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   ZYMB(*)  kappa matrix 1  unpacked
C
#include <implicit.h>
      DIMENSION ZYMB(NORBT,NORBT)
      DIMENSION CMO(*),DXCAO(*),WRK(*)
C
      PARAMETER (D1 = 1.D0, DM1 = -1.D0)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
C
C
      CALL DZERO(DXCAO,N2BASX)
C
C
      CALL CDEN1(ISYMB,CMO,1,D1 ,ZYMB,DXCAO,WRK,LWRK)
      CALL CDEN1(ISYMB,CMO,2,DM1,ZYMB,DXCAO,WRK,LWRK)
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in CDENS1'
         CALL OUTPUT(DXCAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE CDEN1(ISYMB,CMO,JINDX,SIGN,
     &                  ZYMB,DXCAO,WRK,LWRK)
C
C
C Purpose: Computes CMO*ZYMB*transpose(CMO)
C          and adds the result multiplied with SIGN to DXCAO
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   ZYMB(*)  kappa matrix 1  unpacked
C   JINDX    indicates the inactive sumation index
C
#include <implicit.h>
      DIMENSION ZYMB(NORBT,NORBT)
      DIMENSION CMO(*),DXCAO(*),WRK(*)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
C
C
      KDXAO1 = 1
      KDXAO2 = KDXAO1 + N2BASX
      KFREE  = KDXAO2 + N2ORBX
      LFREE  = LWRK - KFREE
      IF ( LFREE .LT. 0 ) CALL ERRWRK('CDEN1',LFREE,LWRK)
C
      CALL DZERO(WRK(KDXAO2),N2BASX)
C
C loop over symmetries then put result in DXCAO
C
      DO 1000 ISYM1 = 1,NSYM
         ISYM2  = MULD2H(ISYM1,ISYMB)
         IORB1  = IORB(ISYM1)
         IORB2  = IORB(ISYM2)
         NORB1  = NORB(ISYM1)
         NORB2  = NORB(ISYM2)
         ICMO1  = ICMO(ISYM1)
         ICMO2  = ICMO(ISYM2)
         NBAS1  = NBAS(ISYM1)
         NBAS2  = NBAS(ISYM2)
         IF (JINDX.EQ.1) NORB1 = NISH(ISYM1)
         IF (NORB1.EQ.0) GOTO 1000
         IF (JINDX.EQ.2) NORB2 = NISH(ISYM2)
         IF (NORB2.EQ.0) GOTO 1000
C
         CALL DGEMM('N','N',NBAS1,NORB2,NORB1,1.D0,
     &              CMO(ICMO1+1),NBAS1,
     &              ZYMB(IORB1+1,IORB2+1),NORBT,0.D0,
     &              WRK(KDXAO1),NBAS1)
         IDXAO2 = KDXAO2 + IBAS(ISYM2)*NBAST + IBAS(ISYM1)
         CALL DGEMM('N','T',NBAS1,NBAS2,NORB2,1.D0,
     &              WRK(KDXAO1),NBAS1,
     &              CMO(ICMO2+1),NBAS2,0.D0,
     &              WRK(IDXAO2),NBAST)
C
 1000 CONTINUE
C
      CALL DAXPY(N2BASX,SIGN,WRK(KDXAO2),1,DXCAO,1)
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in CDEN1'
         CALL OUTPUT(WRK(KDXAO2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE CDENS2(ISYMB,ISYMC,CMO,
     &                  ZYMB,ZYMC,DXCAO,DXAO1,DXAO2,XKAP)
C
C FEB-26 1995
C
C Purpose:
C To compute and add together 4 density matrices to give
C D-mn = DXCAO = (see notes)
C
C needed in DDCRPA (closed shell) for the construction of
C Fock matrix
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   ZYMB(*)  kappa matrix 1  unpacked
C   ZYMC(*)  kappa matrix 2  unpacked
C
#include <implicit.h>
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT)
      DIMENSION CMO(*),DXCAO(*),DXAO1(*),DXAO2(*),XKAP(*)
C
      PARAMETER (D1 = 1.D0, DM1 = -1.D0)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
C
      CALL CDEN2(ISYMB,ISYMC,CMO,1,D1 ,ZYMB,ZYMC,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN2(ISYMB,ISYMC,CMO,2,DM1,ZYMB,ZYMC,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN2(ISYMC,ISYMB,CMO,2,DM1,ZYMC,ZYMB,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN2(ISYMC,ISYMB,CMO,3,D1 ,ZYMC,ZYMB,DXCAO,DXAO1,DXAO2,XKAP)
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in CDENS2'
         CALL OUTPUT(DXCAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
C
C *** end of subroutine CDENS2
C
      END
      SUBROUTINE CDEN2(ISYMB,ISYMC,CMO,JINDX,SIGN,
     &                  ZYMB,ZYMC,DXCAO,DXAO1,DXAO2,XKAP)
C
C
C Purpose: Computes CMO*ZYMB*ZYMC*transpose(CMO)
C          and adds the result multiplied with SIGN to DXCAO
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   ZYMB(*)  kappa matrix 1  unpacked
C   ZYMC(*)  kappa matrix 2  unpacked
C   JINDX    indicates the inactive sumation index
C
#include <implicit.h>
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT)
      DIMENSION CMO(*),DXCAO(*),DXAO1(N2BASX),DXAO2(N2BASX),XKAP(N2ORBX)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
C
C
      CALL DZERO(DXAO2,N2BASX)
C
C loop over symmetries then put result in DXCAO
C
      DO 1000 ISYM1 = 1,NSYM
         ISYM2  = MULD2H(ISYM1,ISYMB)
         ISYM3  = MULD2H(ISYM2,ISYMC)
         IORB1  = IORB(ISYM1)
         IORB2  = IORB(ISYM2)
         IORB3  = IORB(ISYM3)
         NORB1  = NORB(ISYM1)
         NORB2  = NORB(ISYM2)
         NORB3  = NORB(ISYM3)
         ICMO1  = ICMO(ISYM1)
         ICMO3  = ICMO(ISYM3)
         NBAS1  = NBAS(ISYM1)
         NBAS3  = NBAS(ISYM3)
         IF (JINDX.EQ.1) NORB1 = NISH(ISYM1)
         IF (NORB1.EQ.0) GOTO 1000
         IF (JINDX.EQ.2) NORB2 = NISH(ISYM2)
         IF (NORB2.EQ.0) GOTO 1000
         IF (JINDX.EQ.3) NORB3 = NISH(ISYM3)
         IF (NORB3.EQ.0) GOTO 1000
C
         CALL DGEMM('N','N',NORB1,NORB3,NORB2,1.D0,
     &              ZYMB(IORB1+1,IORB2+1),NORBT,
     &              ZYMC(IORB2+1,IORB3+1),NORBT,0.D0,
     &              XKAP,NORB1)
         CALL DGEMM('N','N',NBAS1,NORB3,NORB1,1.D0,
     &              CMO(ICMO1+1),NBAS1,
     &              XKAP,NORB1,0.D0,
     &              DXAO1,NBAS1)
         IDXAO2 = 1 + IBAS(ISYM3)*NBAST + IBAS(ISYM1)
         CALL DGEMM('N','T',NBAS1,NBAS3,NORB3,1.D0,
     &              DXAO1,NBAS1,
     &              CMO(ICMO3+1),NBAS3,0.D0,
     &              DXAO2(IDXAO2),NBAST)
 1000 CONTINUE
C
      CALL DAXPY(N2BASX,SIGN,DXAO2,1,DXCAO,1)
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in CDEN2'
         CALL OUTPUT(DXAO2,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE CDENS3(ISYMB,ISYMC,ISYMD,CMO,
     &                  ZYMB,ZYMC,ZYMD,DXCAO,DXAO1,DXAO2,XKAP)
C
C FEB-26 1995
C
C Purpose:
C To compute and add together 8 density matrices to give
C D-123 = DXCAO = (see notes)
C
C needed in DDCRPA (closed shell) for the construction of
C Fock matrix
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   ZYMB(*)  kappa matrix 1  unpacked
C   ZYMC(*)  kappa matrix 2  unpacked
C   ZYMD(*)  kappa matrix 3  unpacked
C
#include <implicit.h>
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT),ZYMD(NORBT,NORBT)
      DIMENSION CMO(*),DXCAO(*),DXAO1(*),DXAO2(*),XKAP(*)
C
      PARAMETER (D1 = 1.D0, DM1 = -1.D0)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
C
      CALL CDEN3(ISYMB,ISYMC,ISYMD,CMO,1,D1,
     &                  ZYMB,ZYMC,ZYMD,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN3(ISYMB,ISYMC,ISYMD,CMO,2,DM1,
     &                  ZYMB,ZYMC,ZYMD,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN3(ISYMC,ISYMB,ISYMD,CMO,2,DM1,
     &                  ZYMC,ZYMB,ZYMD,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN3(ISYMC,ISYMB,ISYMD,CMO,3,D1,
     &                  ZYMC,ZYMB,ZYMD,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN3(ISYMD,ISYMB,ISYMC,CMO,2,DM1,
     &                  ZYMD,ZYMB,ZYMC,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN3(ISYMD,ISYMB,ISYMC,CMO,3,D1,
     &                  ZYMD,ZYMB,ZYMC,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN3(ISYMD,ISYMC,ISYMB,CMO,3,D1,
     &                  ZYMD,ZYMC,ZYMB,DXCAO,DXAO1,DXAO2,XKAP)
      CALL CDEN3(ISYMD,ISYMC,ISYMB,CMO,4,DM1,
     &                  ZYMD,ZYMC,ZYMB,DXCAO,DXAO1,DXAO2,XKAP)
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in CDENS3'
         CALL OUTPUT(DXCAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
C
C *** end of subroutine CDENS3
C
      END
      SUBROUTINE CDEN3(ISYMB,ISYMC,ISYMD,CMO,JINDX,SIGN,
     &                  ZYMB,ZYMC,ZYMD,DXCAO,DXAO1,DXAO2,XKAP)
C
C
C Purpose: Computes CMO*ZYMB*ZYMC*ZYMD*transpose(CMO)
C          and adds the result multiplied with SIGN to DXCAO
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   ZYMB(*)  kappa matrix 1  unpacked
C   ZYMC(*)  kappa matrix 2  unpacked
C   ZYMD(*)  kappa matrix 3  unpacked
C   JINDX    indicates the inactive sumation index
C
#include <implicit.h>
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT),ZYMD(NORBT,NORBT)
      DIMENSION CMO(*),DXCAO(*),DXAO1(N2BASX),DXAO2(N2BASX),XKAP(N2ORBX)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
C
C
      CALL DZERO(DXAO2,N2BASX)
C
C loop over symmetries then put result in DXCAO
C
      DO 1000 ISYM1 = 1,NSYM
         ISYM2  = MULD2H(ISYM1,ISYMB)
         ISYM3  = MULD2H(ISYM2,ISYMC)
         ISYM4  = MULD2H(ISYM3,ISYMD)
         IORB1  = IORB(ISYM1)
         IORB2  = IORB(ISYM2)
         IORB3  = IORB(ISYM3)
         IORB4  = IORB(ISYM4)
         NORB1  = NORB(ISYM1)
         NORB2  = NORB(ISYM2)
         NORB3  = NORB(ISYM3)
         NORB4  = NORB(ISYM4)
         ICMO1  = ICMO(ISYM1)
         ICMO4  = ICMO(ISYM4)
         NBAS1  = NBAS(ISYM1)
         NBAS4  = NBAS(ISYM4)
         IF (JINDX.EQ.1) NORB1 = NISH(ISYM1)
         IF (NORB1.EQ.0) GOTO 1000
         IF (JINDX.EQ.2) NORB2 = NISH(ISYM2)
         IF (NORB2.EQ.0) GOTO 1000
         IF (JINDX.EQ.3) NORB3 = NISH(ISYM3)
         IF (NORB3.EQ.0) GOTO 1000
         IF (JINDX.EQ.4) NORB4 = NISH(ISYM4)
         IF (NORB4.EQ.0) GOTO 1000
C
         CALL DGEMM('N','N',NORB1,NORB3,NORB2,1.D0,
     &              ZYMB(IORB1+1,IORB2+1),NORBT,
     &              ZYMC(IORB2+1,IORB3+1),NORBT,0.D0,
     &              DXAO1,NORB1)
         CALL DGEMM('N','N',NORB1,NORB4,NORB3,1.D0,
     &              DXAO1,NORB1,
     &              ZYMD(IORB3+1,IORB4+1),NORBT,0.D0,
     &              XKAP,NORB1)
         CALL DGEMM('N','N',NBAS1,NORB4,NORB1,1.D0,
     &              CMO(ICMO1+1),NBAS1,
     &              XKAP,NORB1,0.D0,
     &              DXAO1,NBAS1)
         IDXAO2 = 1 + IBAS(ISYM4)*NBAST + IBAS(ISYM1)
         CALL DGEMM('N','T',NBAS1,NBAS4,NORB4,1.D0,
     &              DXAO1,NBAS1,
     &              CMO(ICMO4+1),NBAS4,0.D0,
     &              DXAO2(IDXAO2),NBAST)
C
 1000 CONTINUE
C
      CALL DAXPY(N2BASX,SIGN,DXAO2,1,DXCAO,1)
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in CDEN3'
         CALL OUTPUT(DXAO2,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE ZYMUL2(ISYMB,ISYMC,ZYMB,ZYMC,ZYMBC)
C
C
C Purpose: Computes ZYMBC = ZYMBC + ZYMB*ZYMC
C
C Input:
C   ZYMB(*)  kappa matrix 1  unpacked
C   ZYMC(*)  kappa matrix 2  unpacked
C
#include <implicit.h>
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT)
      DIMENSION ZYMBC(NORBT,NORBT)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
C
C loop over symmetries then put result in ZYMBC
C
      DO 1000 ISYM1 = 1,NSYM
         ISYM2  = MULD2H(ISYM1,ISYMB)
         ISYM3  = MULD2H(ISYM2,ISYMC)
         IOCC1  = IORB(ISYM1)
         IOCC2  = IORB(ISYM2)
         IOCC3  = IORB(ISYM3)
         NOCC1  = NOCC(ISYM1)
         NOCC2  = NOCC(ISYM2)
         NOCC3  = NOCC(ISYM3)
         IVIR1  = IOCC1 + NOCC1
         IVIR2  = IOCC2 + NOCC2
         IVIR3  = IOCC3 + NOCC3
         NVIR1  = NVIR(ISYM1)
         NVIR2  = NVIR(ISYM2)
         NVIR3  = NVIR(ISYM3)
C
         IF (NOCC1*NOCC3*NVIR2 .NE. 0)
     &   CALL DGEMM('N','N',NOCC1,NOCC3,NVIR2,1.D0,
     &              ZYMB(IOCC1+1,IVIR2+1),NORBT,
     &              ZYMC(IVIR2+1,IOCC3+1),NORBT,1.D0,
     &              ZYMBC(IOCC1+1,IOCC3+1),NORBT)
         IF (NVIR1*NVIR3*NOCC2 .NE. 0)
     &   CALL DGEMM('N','N',NVIR1,NVIR3,NOCC2,1.D0,
     &              ZYMB(IVIR1+1,IOCC2+1),NORBT,
     &              ZYMC(IOCC2+1,IVIR3+1),NORBT,1.D0,
     &              ZYMBC(IVIR1+1,IVIR3+1),NORBT)
C
 1000 CONTINUE
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in ZYMUL2'
         CALL OUTPUT(ZYMBC,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE ZYMUL3(ISYMBC,ISYMD,ZYMBC,ZYMD,ZYMBCD)
C
C
C Purpose: Computes ZYMBCD = ZYMBCD + ZYMBC*ZYMD
C
C Input:
C   ZYMBC(*)  kappa matrix product
C   ZYMD (*)  kappa matrix unpacked
C
#include <implicit.h>
      DIMENSION ZYMBC(NORBT,NORBT),ZYMD(NORBT,NORBT)
      DIMENSION ZYMBCD(NORBT,NORBT)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
C
C loop over symmetries then put result in ZYMBCD
C
      DO 1000 ISYM1 = 1,NSYM
         ISYM2  = MULD2H(ISYM1,ISYMBC)
         ISYM3  = MULD2H(ISYM2,ISYMD)
         IOCC1  = IORB(ISYM1)
         IOCC2  = IORB(ISYM2)
         IOCC3  = IORB(ISYM3)
         NOCC1  = NOCC(ISYM1)
         NOCC2  = NOCC(ISYM2)
         NOCC3  = NOCC(ISYM3)
         IVIR1  = IOCC1 + NOCC1
         IVIR2  = IOCC2 + NOCC2
         IVIR3  = IOCC3 + NOCC3
         NVIR1  = NVIR(ISYM1)
         NVIR2  = NVIR(ISYM2)
         NVIR3  = NVIR(ISYM3)
C
         IF (NOCC1*NVIR3*NOCC2 .NE. 0)
     &   CALL DGEMM('N','N',NOCC1,NVIR3,NOCC2,1.D0,
     &              ZYMBC(IOCC1+1,IOCC2+1),NORBT,
     &              ZYMD(IOCC2+1,IVIR3+1),NORBT,1.D0,
     &              ZYMBCD(IOCC1+1,IVIR3+1),NORBT)
         IF (NVIR1*NOCC3*NVIR2 .NE. 0)
     &   CALL DGEMM('N','N',NVIR1,NOCC3,NVIR2,1.D0,
     &              ZYMBC(IVIR1+1,IVIR2+1),NORBT,
     &              ZYMD(IVIR2+1,IOCC3+1),NORBT,1.D0,
     &              ZYMBCD(IVIR1+1,IOCC3+1),NORBT)
C
 1000 CONTINUE
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in ZYMUL3'
         CALL OUTPUT(ZYMBCD,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE ZYTRA1(ISYM,ZYM,CMO,TMP,DAO)
C
C
C Purpose: Computes DAO = CMO * ZYM[0,+;-,0] * transp(CMO)
C                         = CMO * transp( CMO * transp(ZYM[]) )
C          
C          ZYM[0,+;-,0](M,N) = +ZYM(M,N), M virtual, N occupied
C                            = -ZYM(M,N), N virtual, M occupied
C                            = 0, otherwhise
C
C Input:
C   ZYM(*)  kappa matrix
C   CMO(*)  MO coefficients
C   TMP(*)  temporary matrix
C
#include <implicit.h>
      DIMENSION ZYM(NORBT,NORBT),TMP(NORBT,NORBT)
      DIMENSION CMO(*),DAO(NBAST,NBAST)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
C
C loop over symmetries
C
      DO 1000 ISYM1 = 1,NSYM
         ISYM2  = MULD2H(ISYM1,ISYM)
         IOCC1  = IORB(ISYM1)
         IOCC2  = IORB(ISYM2)
         NOCC1  = NOCC(ISYM1)
         NOCC2  = NOCC(ISYM2)
         IVIR1  = IOCC1 + NOCC1
         IVIR2  = IOCC2 + NOCC2
         NVIR1  = NVIR(ISYM1)
         NVIR2  = NVIR(ISYM2)
         IBAS1  = IBAS(ISYM1)
         IBAS2  = IBAS(ISYM2)
         NBAS1  = NBAS(ISYM1)
         NBAS2  = NBAS(ISYM2)
         NORB1  = NORB(ISYM1)
         NORB2  = NORB(ISYM2)
         ICMO1  = ICMO(ISYM1)
         ICMO2  = ICMO(ISYM2)
C
         IF (NOCC1*NVIR2 .NE. 0) THEN
            CALL DGEMM('N','T',NBAS2,NOCC1,NVIR2,1.D0,
     &                 CMO(ICMO2+NOCC2*NBAS2+1),NBAS2,
     &                 ZYM(IOCC1+1,IVIR2+1),NORBT,0.D0,
     &                 TMP(IBAS2+1,IOCC1+1),NBAST)
         END IF
         NTMP = NBAS1*NBAS2*NOCC1
         IF (NTMP .NE. 0) THEN
            CALL DGEMM('N','T',NBAS1,NBAS2,NOCC1,1.D0,
     &                 CMO(ICMO1+1),NBAS1,
     &                 TMP(IBAS2+1,IOCC1+1),NBAST,1.D0,
     &                 DAO(IBAS1+1,IBAS2+1),NBAST)
         END IF
         IF (NOCC2*NVIR1 .NE. 0) THEN
            CALL DGEMM('N','T',NBAS2,NVIR1,NOCC2,1.D0,
     &                 CMO(ICMO2+1),NBAS2,
     &                 ZYM(IVIR1+1,IOCC2+1),NORBT,0.D0,
     &                 TMP(IBAS2+1,IVIR1+1),NBAST)
         END IF
         NTMP = NBAS1*NBAS2*NVIR1
         IF (NTMP .NE. 0) THEN
            CALL DGEMM('N','T',NBAS1,NBAS2,NVIR1,-1.D0,
     &                 CMO(ICMO1+NOCC1*NBAS1+1),NBAS1,
     &                 TMP(IBAS2+1,IVIR1+1),NBAST,1.D0,
     &                 DAO(IBAS1+1,IBAS2+1),NBAST)
         END IF
C
 1000 CONTINUE
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in ZYTRA1'
         CALL OUTPUT(DAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE ZYTRA2(ISYM,ZYM,CMO,TMP,DAO)
C
C
C Purpose: Computes DAO = CMO * ZYM[+,0;0,-] * transp(CMO)
C                       = CMO * transp( CMO * transp(ZYM[]) )
C          
C          ZYM[+,0;0,-](M,N) = +ZYM(M,N), M,N occupied
C                            = -ZYM(M,N), M,N virtual
C                            = 0, otherwhise
C
C Input:
C   ZYM(*)  kappa matrix
C   CMO(*)  MO coefficients
C   TMP(*)  temporary matrix
C
#include <implicit.h>
      DIMENSION ZYM(NORBT,NORBT),TMP(NORBT,NORBT)
      DIMENSION CMO(*),DAO(NBAST,NBAST)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
C
C loop over symmetries
C
      DO 1000 ISYM1 = 1,NSYM
         ISYM2  = MULD2H(ISYM1,ISYM)
         IOCC1  = IORB(ISYM1)
         IOCC2  = IORB(ISYM2)
         NOCC1  = NOCC(ISYM1)
         NOCC2  = NOCC(ISYM2)
         IVIR1  = IOCC1 + NOCC1
         IVIR2  = IOCC2 + NOCC2
         NVIR1  = NVIR(ISYM1)
         NVIR2  = NVIR(ISYM2)
         IBAS1  = IBAS(ISYM1)
         IBAS2  = IBAS(ISYM2)
         NBAS1  = NBAS(ISYM1)
         NBAS2  = NBAS(ISYM2)
         NORB1  = NORB(ISYM1)
         NORB2  = NORB(ISYM2)
         ICMO1  = ICMO(ISYM1)
         ICMO2  = ICMO(ISYM2)
C
         IF (NOCC1*NOCC2 .NE. 0) THEN
            CALL DGEMM('N','T',NBAS2,NOCC1,NOCC2,1.D0,
     &                 CMO(ICMO2+1),NBAS2,
     &                 ZYM(IOCC1+1,IOCC2+1),NORBT,0.D0,
     &                 TMP(IBAS2+1,IOCC1+1),NBAST)
         END IF
         NTMP = NBAS1*NBAS2*NOCC1
         IF (NTMP .NE.0) THEN
            CALL DGEMM('N','T',NBAS1,NBAS2,NOCC1,1.D0,
     &                 CMO(ICMO1+1),NBAS1,
     &                 TMP(IBAS2+1,IOCC1+1),NBAST,1.D0,
     &                 DAO(IBAS1+1,IBAS2+1),NBAST)
         END IF
         IF (NVIR1*NVIR2 .NE. 0) THEN
            CALL DGEMM('N','T',NBAS2,NVIR1,NVIR2,1.D0,
     &                 CMO(ICMO2+NOCC2*NBAS2+1),NBAS2,
     &                 ZYM(IVIR1+1,IVIR2+1),NORBT,0.D0,
     &                 TMP(IBAS2+1,IVIR1+1),NBAST)
         END IF
         NTMP = NBAS1*NBAS2*NVIR1
         IF (NTMP .NE.0) THEN
            CALL DGEMM('N','T',NBAS1,NBAS2,NVIR1,-1.D0,
     &                 CMO(ICMO1+NOCC1*NBAS1+1),NBAS1,
     &                 TMP(IBAS2+1,IVIR1+1),NBAST,1.D0,
     &                 DAO(IBAS1+1,IBAS2+1),NBAST)
         END IF
C
 1000 CONTINUE
C
      IF ( IPRRSP.GT.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in ZYTRA2'
         CALL OUTPUT(DAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
      END
