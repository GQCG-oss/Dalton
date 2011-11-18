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
C  /* Deck cc_inifro */
      SUBROUTINE CC_INIFRO(WORK,LWORK)
C
C     Written by Asger Halkier 22/5 - 1998.
C
C     Set up index arrays needed for frozen core gradients.
C
#include <implicit.h>
C
      DIMENSION WORK(LWORK)
C
      EXTERNAL INIDAT
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <symsq.h>
#include <ccisao.h>
#include <r12int.h>
C
C----------------------------------------
C     Calculate index arrays for lengths.
C----------------------------------------
C
      DO 100 ISYMAI = 1,NSYM
         NT1FRO(ISYMAI) = 0
         NF1FRO(ISYMAI) = 0
         NDSFRO(ISYMAI) = 0
         NCOFRO(ISYMAI) = 0
         NFROIJ(ISYMAI) = 0
         NFROFR(ISYMAI) = 0
         NV1FRO(ISYMAI) = 0
         DO 110 ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            NT1FRO(ISYMAI) = NT1FRO(ISYMAI) + NVIR(ISYMA)*NRHFFR(ISYMI)
            NF1FRO(ISYMAI) = NF1FRO(ISYMAI) + NRHF(ISYMA)*NRHFFR(ISYMI)
            NDSFRO(ISYMAI) = NDSFRO(ISYMAI) + NNBST(ISYMA)*NRHFFR(ISYMI)
            NCOFRO(ISYMAI) = NCOFRO(ISYMAI) + NRHF(ISYMA)*NRHFFR(ISYMI)
            NFROIJ(ISYMAI) = NFROIJ(ISYMAI) + NRHFS(ISYMA)*NRHFS(ISYMI)
            NFROFR(ISYMAI) = NFROFR(ISYMAI)
     *                     + NRHFFR(ISYMA)*NRHFFR(ISYMI)
celena
            NV1FRO(ISYMAI) = NV1FRO(ISYMAI) + NRHFFR(ISYMI)*
     &                          (NORB1(ISYMA)-NRHFFR(ISYMA))
celena

  110    CONTINUE
  100 CONTINUE
C
      DO 120 ISYMAI = 1,NSYM
         NT2FRO(ISYMAI) = 0
         NF2FRO(ISYMAI) = 0
         NALLAI(ISYMAI) = 0
         NOFROO(ISYMAI) = 0
         NF2IJG(ISYMAI) = 0
         NFRIJK(ISYMAI) = 0
         NCOFRF(ISYMAI) = 0
celena
         NFROVF(ISYMAI) = 0
         NFROVR(ISYMAI) = 0
celena
         DO 130 ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            NT2FRO(ISYMAI) = NT2FRO(ISYMAI) + NT1FRO(ISYMA)*NT1AM(ISYMI)
            NF2FRO(ISYMAI) = NF2FRO(ISYMAI) 
     *                     + NT1FRO(ISYMA)*NF1FRO(ISYMI)
            NALLAI(ISYMAI) = NALLAI(ISYMAI) + NVIRS(ISYMA)*NRHFS(ISYMI)
            NOFROO(ISYMAI) = NOFROO(ISYMAI) + NCOFRO(ISYMA)*NRHF(ISYMI)
            NF2IJG(ISYMAI) = NF2IJG(ISYMAI) + NFROIJ(ISYMA)*NBAS(ISYMI)
            NFRIJK(ISYMAI) = NFRIJK(ISYMAI) + NFROIJ(ISYMA)*NRHFS(ISYMI)
            NCOFRF(ISYMAI) = NCOFRF(ISYMAI)
     *                     + NCOFRO(ISYMA)*NRHFFR(ISYMI)
celena
            NFROVF(ISYMAI) =  NFROVF(ISYMAI)+NFROFR(ISYMI)*NT1FRO(ISYMA)
            NFROVR(ISYMAI) =  NFROVR(ISYMAI)+NV1FRO(ISYMI)*NT1FRO(ISYMA)
celena
  130    CONTINUE
  120 CONTINUE
C
C--------------------------------------------
C     Calculate symmetry offset index arrays.
C--------------------------------------------
C
      DO 140 ISYMAI = 1,NSYM
         ICOUN1  = 0
         ICOUN2  = 0
         ICOUN3  = 0
         ICOUN4  = 0
         ICOUN5  = 0
         ICOUN6  = 0
         ICOUN7  = 0
         ICOUN8  = 0
         ICOUN9  = 0
         ICOUN10 = 0
         ICOUN11 = 0
         ICOUN12 = 0
         ICOUN13 = 0
!sonia
         ICOUN14 = 0
!wim
         ICOUN15 = 0
         ICOUN16 = 0

         DO 150 ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            IT1FRO(ISYMA,ISYMI) = ICOUN1
            IT2FRO(ISYMA,ISYMI) = ICOUN2
            IDSFRO(ISYMA,ISYMI) = ICOUN3
            ICOFRO(ISYMA,ISYMI) = ICOUN4
            IALLAI(ISYMA,ISYMI) = ICOUN5
            IOFROO(ISYMA,ISYMI) = ICOUN6
            IF2IJG(ISYMA,ISYMI) = ICOUN7
            IFROIJ(ISYMA,ISYMI) = ICOUN8
            IFRIJK(ISYMA,ISYMI) = ICOUN9
            IFROFR(ISYMA,ISYMI) = ICOUN10
            ICOFRF(ISYMA,ISYMI) = ICOUN11
            IOFOAO(ISYMA,ISYMI) = ICOUN12
            IOFFAO(ISYMA,ISYMI) = ICOUN13
!sonia
            ICDKFR(ISYMA,ISYMI) = ICOUN14
!wim
            IF1FRO(ISYMA,ISYMI) = ICOUN15
            IF2FRO(ISYMA,ISYMI) = ICOUN16

            ICOUN1  = ICOUN1  + NVIR(ISYMA)*NRHFFR(ISYMI)
            ICOUN2  = ICOUN2  + NT1FRO(ISYMA)*NT1AM(ISYMI)
            ICOUN3  = ICOUN3  + NNBST(ISYMA)*NRHFFR(ISYMI)
            ICOUN4  = ICOUN4  + NRHF(ISYMA)*NRHFFR(ISYMI)
            ICOUN5  = ICOUN5  + NVIRS(ISYMA)*NRHFS(ISYMI)
            ICOUN6  = ICOUN6  + NCOFRO(ISYMA)*NRHF(ISYMI)
            ICOUN7  = ICOUN7  + NFROIJ(ISYMA)*NBAS(ISYMI)
            ICOUN8  = ICOUN8  + NRHFS(ISYMA)*NRHFS(ISYMI)
            ICOUN9  = ICOUN9  + NFROIJ(ISYMA)*NRHFS(ISYMI)
            ICOUN10 = ICOUN10 + NRHFFR(ISYMA)*NRHFFR(ISYMI)
            ICOUN11 = ICOUN11 + NCOFRO(ISYMA)*NRHFFR(ISYMI)
            ICOUN12 = ICOUN12 + NOFROO(ISYMA)*NBAS(ISYMI)
            ICOUN13 = ICOUN13 + NCOFRF(ISYMA)*NBAS(ISYMI)
!sonia
            ICOUN14 = ICOUN14 + NMATAB(ISYMA)*NRHF(ISYMI)
!wim
            ICOUN15 = ICOUN15 + NRHFFR(ISYMA)*NRHF(ISYMI)
            ICOUN16 = ICOUN16 + NT1FRO(ISYMA)*NF1FRO(ISYMI)

  150    CONTINUE
  140 CONTINUE
C
      ICOUN1 = 0
      ICOUN2 = NRHFTS
      DO 160 ISYM = 1,NSYM
         IRHFS(ISYM) = ICOUN1
         IVIRS(ISYM) = ICOUN2
         ICOUN1 = ICOUN1 + NRHFS(ISYM)
         ICOUN2 = ICOUN2 + NVIRS(ISYM)
  160 CONTINUE
C
C----------------------
C     Print if desired.
C----------------------
C
      IF (IPRINT .GT. 6) THEN
         CALL AROUND('Information from CC_INIFRO')
         WRITE(LUPRI,1) 'NT1FRO :',(NT1FRO(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NT2FRO :',(NT2FRO(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NDSFRO :',(NDSFRO(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NCOFRO :',(NCOFRO(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NALLAI :',(NALLAI(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NOFROO :',(NOFROO(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NF2IJG :',(NF2IJG(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NFROIJ :',(NFROIJ(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NFRIJK :',(NFRIJK(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NFROFR :',(NFROFR(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NCOFRF :',(NCOFRF(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'IRHFS  :',(IRHFS(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'IVIRS  :',(IVIRS(I),   I=1,NSYM)
      ENDIF
      IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI,*)
         DO 11 I = 1,NSYM
            WRITE(LUPRI,1) 'IT1FRO :',(IT1FRO(I,J), J=1,NSYM)
  11     CONTINUE
         DO 12 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2FRO :',(IT2FRO(I,J), J=1,NSYM)
  12     CONTINUE
         DO 13 I = 1,NSYM
            WRITE(LUPRI,1) 'IDSFRO :',(IDSFRO(I,J), J=1,NSYM)
  13     CONTINUE
         DO 14 I = 1,NSYM
            WRITE(LUPRI,1) 'ICOFRO :',(ICOFRO(I,J), J=1,NSYM)
  14     CONTINUE
         DO 15 I = 1,NSYM
            WRITE(LUPRI,1) 'IALLAI :',(IALLAI(I,J), J=1,NSYM)
  15     CONTINUE
         DO 16 I = 1,NSYM
            WRITE(LUPRI,1) 'IOFROO :',(IOFROO(I,J), J=1,NSYM)
  16     CONTINUE
         DO 17 I = 1,NSYM
            WRITE(LUPRI,1) 'IF2IJG :',(IF2IJG(I,J), J=1,NSYM)
  17     CONTINUE
         DO 18 I = 1,NSYM
            WRITE(LUPRI,1) 'IFROIJ :',(IFROIJ(I,J), J=1,NSYM)
  18     CONTINUE
         DO 19 I = 1,NSYM
            WRITE(LUPRI,1) 'IFRIJK :',(IFRIJK(I,J), J=1,NSYM)
  19     CONTINUE
         DO 20 I = 1,NSYM
            WRITE(LUPRI,1) 'IFROFR :',(IFROFR(I,J), J=1,NSYM)
  20     CONTINUE
         DO 21 I = 1,NSYM
            WRITE(LUPRI,1) 'ICOFRF :',(ICOFRF(I,J), J=1,NSYM)
  21     CONTINUE
         DO 22 I = 1,NSYM
            WRITE(LUPRI,1) 'IOFOAO :',(IOFOAO(I,J), J=1,NSYM)
  22     CONTINUE
         DO 23 I = 1,NSYM
            WRITE(LUPRI,1) 'IOFFAO :',(IOFFAO(I,J), J=1,NSYM)
  23     CONTINUE
         WRITE(LUPRI,*)
      END IF
C
      RETURN
C
  1   FORMAT(3X,A8,8I8)
C
      END
C  /* Deck cc_kabre */
      SUBROUTINE CC_KABRE(SOLU,ZKAM,WORK,LWORK)
C
C     Written by Asger Halkier 22/5 - 1998.
C
C     To reorder the solution vector from the abacus-response
C     solver to cc-standards.
C
#include <implicit.h>
C
      DIMENSION SOLU(*), ZKAM(*), WORK(LWORK)
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      IF (FROIMP) THEN
C
         KOFF1 = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
         CALL DZERO(ZKAM(KOFF1),2*NT1FRO(1))
C
         DO 100 ISYMI = 1,NSYM
            ISYMA = ISYMI
C
            IF (NRHFFR(ISYMI) .GT. 0) THEN
               LEN1  = NVIR(ISYMA)*NRHFFR(ISYMI)
               KOFF2 = IALLAI(ISYMA,ISYMI) + 1
               KOFF3 = KOFF1 + IT1FRO(ISYMA,ISYMI)
               KOFF4 = KOFF3 + NT1FRO(1)
               CALL DCOPY(LEN1,SOLU(KOFF2),1,ZKAM(KOFF3),1)
               CALL DCOPY(LEN1,SOLU(KOFF2),1,ZKAM(KOFF4),1)
            ENDIF
C
            LEN2  = NVIR(ISYMA)*NRHF(ISYMI)
            KOFF5 = IALLAI(ISYMA,ISYMI) + NVIR(ISYMA)*NRHFFR(ISYMI) + 1
            KOFF6 = NMATAB(1) + NMATIJ(1) + IT1AM(ISYMA,ISYMI) + 1
            KOFF7 = KOFF6 + NT1AMX
            CALL DCOPY(LEN2,SOLU(KOFF5),1,ZKAM(KOFF6),1)
            CALL DCOPY(LEN2,SOLU(KOFF5),1,ZKAM(KOFF7),1)
C
            IF (IPRINT .GT. 20) THEN
C
               WRITE(LUPRI,*) ' '
               WRITE(LUPRI,1) 'Symmetry block number:', ISYMI
  1            FORMAT(3X,A22,2X,I1)
C
               KOFFTO = IALLAI(ISYMA,ISYMI) + 1
               KOFFCP = KOFF6
               KOFFFP = KOFF1 + IT1FRO(ISYMA,ISYMI)
C
               CALL AROUND('Total solution vector')
               CALL OUTPUT(SOLU(KOFFTO),1,NVIR(ISYMA),1,NRHFS(ISYMI),
     *                     NVIR(ISYMA),NRHFS(ISYMI),1,LUPRI)
               CALL AROUND('Correlated part')
               CALL OUTPUT(ZKAM(KOFFCP),1,NVIR(ISYMA),1,NRHF(ISYMI),
     *                     NVIR(ISYMA),NRHF(ISYMI),1,LUPRI)
               CALL AROUND('Frozen part')
               CALL OUTPUT(ZKAM(KOFFFP),1,NVIR(ISYMA),1,NRHFFR(ISYMI),
     *                     NVIR(ISYMA),NRHFFR(ISYMI),1,LUPRI)
C
            ENDIF
C
  100    CONTINUE
      ELSE
C
         KOFF8 = NMATIJ(1) + NMATAB(1) + 1
         KOFF9 = KOFF8 + NT1AMX
         CALL DCOPY(NT1AMX,SOLU(1),1,ZKAM(KOFF8),1)
         CALL DCOPY(NT1AMX,SOLU(1),1,ZKAM(KOFF9),1)
C
      ENDIF
C
      RETURN
      END
C  /* Deck cc_etare */
      SUBROUTINE CC_ETARE(ETAAI,AFROI,WORK,LWORK)
C
C     Written by Asger Halkier 22/5 - 1998.
C
C     To reorder the eta-kappa-0 vector for the abacus-response
C     solver.
C
#include <implicit.h>
C
      DIMENSION ETAAI(*), AFROI(*), WORK(LWORK)
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      IF (FROIMP) THEN
C
         IF (LWORK .LT. NT1AMX) THEN
            WRITE(LUPRI,*) 'Needed:', NT1AMX, 'Available:', LWORK
            CALL QUIT('Insufficient memory in cc_etare')
         ENDIF
C
         CALL DCOPY(NT1AMX,ETAAI,1,WORK,1)
         CALL DZERO(ETAAI,NT1AMX)
C
         DO 100 ISYMI = 1,NSYM
            ISYMA = ISYMI
C
            IF (NRHFFR(ISYMI) .GT. 0) THEN
               LEN1  = NVIR(ISYMA)*NRHFFR(ISYMI)
               KOFF1 = IT1FRO(ISYMA,ISYMI) + 1
               KOFF2 = IALLAI(ISYMA,ISYMI) + 1
               CALL DCOPY(LEN1,AFROI(KOFF1),1,ETAAI(KOFF2),1)
            ENDIF
C
            LEN2  = NVIR(ISYMA)*NRHF(ISYMI)
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
            KOFF4 = IALLAI(ISYMA,ISYMI) + NVIR(ISYMA)*NRHFFR(ISYMI) + 1
            CALL DCOPY(LEN2,WORK(KOFF3),1,ETAAI(KOFF4),1)
C
            IF (IPRINT .GT. 20) THEN
C
               WRITE(LUPRI,*) ' '
               WRITE(LUPRI,1) 'Symmetry block number:', ISYMI
  1            FORMAT(3X,A22,2X,I1)
C
               KOFFTO = IALLAI(ISYMA,ISYMI) + 1
               KOFFCP = KOFF3
               KOFFFP = IT1FRO(ISYMA,ISYMI) + 1
C
               CALL AROUND('Total RHS vector')
               CALL OUTPUT(ETAAI(KOFFTO),1,NVIR(ISYMA),1,NRHFS(ISYMI),
     *                     NVIR(ISYMA),NRHFS(ISYMI),1,LUPRI)
               CALL AROUND('Correlated part')
               CALL OUTPUT(WORK(KOFFCP),1,NVIR(ISYMA),1,NRHF(ISYMI),
     *                     NVIR(ISYMA),NRHF(ISYMI),1,LUPRI)
               CALL AROUND('Frozen part')
               CALL OUTPUT(AFROI(KOFFFP),1,NVIR(ISYMA),1,NRHFFR(ISYMI),
     *                     NVIR(ISYMA),NRHFFR(ISYMI),1,LUPRI)
C
            ENDIF
C
  100    CONTINUE
      ENDIF
C
      RETURN
      END
C  /* Deck mp2_zkfcb */
      SUBROUTINE MP2_ZKFCB(IPDD,R12PRP,RES,T2MAM,WORK,LWORK)
C
C     Written by Asger Halkier 25/5 - 1998.
C
C     To calculate the frozen-occupied block (iJ) of kappa-bar-0.
C
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
      DIMENSION RES(*), T2MAM(*), WORK(LWORK)
      INTEGER YIJIJ, YIJIJT
      LOGICAL   R12PRP
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KINTFR = 1
      KPIMIJ = KINTFR + NT2FRO(1)
      KFOCKD = KPIMIJ + NCOFRO(1)
      KEND1  = KFOCKD + NORBTS
      LWRK1  = LWORK  - KEND1

celena
C-------------------------------------
C     Read R12-contributions
C-------------------------------------
 
      IF (R12PRP .AND. (IPDD .EQ. 2 .OR. IPDD .EQ. 3 
     &    .OR. IPDD .EQ. 5)) THEN
         DO ISYMAJ = 1,NSYM
            ISYMIJ = ISYMAJ
            NCVIJ = 0
            DO ISYMA = 1,NSYM
               ISYMJ = MULD2H(ISYMAJ,ISYMA)
               ISYMI = MULD2H(ISYMIJ,ISYMJ)
               NCVIJ = NCVIJ + NRHFFR(ISYMA)*NRHF(ISYMI)
            ENDDO
         ENDDO
         YIJIJ  = KEND1
         YIJIJT = YIJIJ + NCVIJ  
         KEND1  = YIJIJT+ NCVIJ
         LWRK1  = LWORK  - KEND1
         LUVAJKL = -1
         IF (IPDD .EQ. 2) THEN
            CALL GPOPEN(LUVAJKL,'CCR12YIJIJ','UNKNOWN',' ',
     *                  'UNFORMATTED',IDUMMY,.FALSE.)
         ELSEIF (IPDD .EQ. 3) THEN
            CALL GPOPEN(LUVAJKL,'CCR12ZIJIJ','UNKNOWN',' ',
     *                  'UNFORMATTED',IDUMMY,.FALSE.)
         ELSEIF (IPDD .EQ. 5) THEN
            CALL GPOPEN(LUVAJKL,'CCR12XIJIJ','UNKNOWN',' ',
     *                  'UNFORMATTED',IDUMMY,.FALSE.)
         ENDIF
         REWIND(LUVAJKL)
         READ(LUVAJKL) (WORK(YIJIJ+I-1),I=1,NCVIJ)
         CALL GPCLOSE(LUVAJKL,'DELETE')
         DO ISYMAJ = 1,NSYM
            ISYMIJ = ISYMAJ
            DO ISYMA = 1,NSYM
               ISYMJ = MULD2H(ISYMA,ISYMAJ)
               ISYMI = MULD2H(ISYMA,ISYMIJ)
               do i =1,nrhf(isymi)
                  do j =1,nrhffr(isymj)
                     idxij = icofro(isymj,isymi)+
     &                          nrhffr(isymj)*(i-1) + j
                     idxji = icofro(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                     work(YIJIJT+idxji-1) = work(YIJIJ+idxij-1)
                  enddo
               enddo
            ENDDO
         ENDDO
      ENDIF
celena


C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_ZKFCB')
      ENDIF
C
      CALL DZERO(WORK(KPIMIJ),NCOFRO(1))
      CALL DZERO(WORK(KFOCKD),NORBTS)
C
C--------------------------------------
C     Read integrals (cJ|dk) from disk.
C--------------------------------------
C
      LUCJDK = -1
      CALL GPOPEN(LUCJDK,'INCJDK','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCJDK)
      READ(LUCJDK) (WORK(KINTFR+I-1), I = 1,NT2FRO(1))
      CALL GPCLOSE(LUCJDK,'KEEP')

C
      IF (IPRINT. GT. 20) THEN
         XFRNOR = DDOT(NT2FRO(1),WORK(KINTFR),1,WORK(KINTFR),1)
         WRITE(LUPRI,*) 'Norm of integrals (cJdk):', XFRNOR
      ENDIF
C
C---------------------------------------
C     Contract integrals with amplitude.
C---------------------------------------
C
      DO 100 ISYMDK = 1,NSYM
         ISYMCI = ISYMDK
         ISYMCJ = ISYMDK
         DO 110 NDK = 1,NT1AM(ISYMDK)
            DO 120 ISYMC = 1,NSYM
               ISYMI = MULD2H(ISYMC,ISYMCI)
               ISYMJ = MULD2H(ISYMC,ISYMCJ)
C
               KOFF1 = IT2SQ(ISYMCI,ISYMDK)
     *               + NT1AM(ISYMCI)*(NDK - 1)
     *               + IT1AM(ISYMC,ISYMI) + 1
               KOFF2 = KINTFR + IT2FRO(ISYMCJ,ISYMDK)
     *               + NT1FRO(ISYMCJ)*(NDK - 1)
     *               + IT1FRO(ISYMC,ISYMJ)
               KOFF3 = KPIMIJ + ICOFRO(ISYMI,ISYMJ)
C
               NTOTC = MAX(NVIR(ISYMC),1)
               NTOTI = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),
     *                    NVIR(ISYMC),ONE,T2MAM(KOFF1),NTOTC,
     *                    WORK(KOFF2),NTOTC,ONE,WORK(KOFF3),NTOTI)
cC
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
celena
      IF (R12PRP) THEN
         CALL DAXPY(NCVIJ,ONE,WORK(YIJIJT),1,WORK(KPIMIJ),1)
      END IF
celena
C
C--------------------------
C     Get orbital energies.
C--------------------------
C
      CALL FOCK_ALL(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C---------------------------
C     Calculate the results.
C---------------------------
C
      DO 130 ISYMI = 1,NSYM
         ISYMJ = ISYMI
         DO 140 J = 1,NRHFFR(ISYMJ)
            DO 150 I = 1,NRHF(ISYMI)
               NIJ = ICOFRO(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
               KOI = KFOCKD + IRHFS(ISYMI) + NRHFFR(ISYMI) + I - 1
               KOJ = KFOCKD + IRHFS(ISYMJ) + J - 1
C
               DENOM    = WORK(KOI) - WORK(KOJ)
               RES(NIJ) = HALF*WORK(KPIMIJ+NIJ-1)/DENOM
C
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      CALL DCOPY(NCOFRO(1),RES(1),1,RES(1+NCOFRO(1)),1)
C
      RETURN
      END
C  /* Deck cmo_all */
      SUBROUTINE CMO_ALL(CMO,WORK,LWORK)
C
C     Written by Asger Halkier 25/5 - 1998
C
C     To read in and change the symmetry ordering of the FULL
C     MO coefficient matrix from Sirius to CC-ordering.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      DIMENSION CMO(*),WORK(LWORK)
#include <inftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <r12int.h>
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
C     Use LABEL instead of 'TRCCINT ' (WK/UniKA/04-11-2002).
      CALL MOLLAB(LABEL,LUSIFC,LUPRI)
      READ (LUSIFC)
C
      READ (LUSIFC)
      READ (LUSIFC) (CMO(I), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KSCR1 = 1
      KEND1 = KSCR1 + NLAMDS
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in CMO_ALL')
      ENDIF
C
C----------------------------------
C     Reorder all orbitals in work.
C----------------------------------
C
      ICRHF  = KSCR1
      ICVIR  = KSCR1 + NLRHSI
C
      ICOUNT = 1
      DO 100 ISYM = 1,NSYM
C
         CALL DCOPY(NBAS(ISYM)*NRHFS(ISYM),CMO(ICOUNT),1,WORK(ICRHF),1)
         ICRHF  = ICRHF  + NBAS(ISYM)*NRHFS(ISYM)
         ICOUNT = ICOUNT + NBAS(ISYM)*NRHFS(ISYM)
C
         CALL DCOPY(NBAS(ISYM)*NVIRS(ISYM),CMO(ICOUNT),1,WORK(ICVIR),1)
         ICVIR  = ICVIR  + NBAS(ISYM)*NVIRS(ISYM)
         ICOUNT = ICOUNT + NBAS(ISYM)*NVIRS(ISYM)
C
  100 CONTINUE
C
C-----------------------------------------
C     Copy reordered result back into CMO.
C-----------------------------------------
C
      CALL DCOPY(NLAMDS,WORK(KSCR1),1,CMO,1)
C
      RETURN
      END
C  /* Deck cc_frcoin */
      SUBROUTINE CC_FRCOIN(XFRIN,XINT,CMO,WORK,LWORK,IDEL,ISYDEL)
C
C     Written by Asger Halkier 25/5 - 1998.
C
C     To calculate the integrals (cJ|dk) needed for frozen-core gradients.
C     AO integrals are assumed totally symmetric.
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XFRIN(*), XINT(*), CMO(*), WORK(LWORK)
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <r12int.h>
C
C---------------------------------------
C     Initial symmetries and allocation.
C---------------------------------------
C
      ISYDIS = ISYDEL
      ISYMD  = ISYDEL
C
      KDVEC  = 1
      KEND1  = KDVEC + NVIR(ISYMD)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CC_FRCOIN')
      ENDIF
C
      CALL DZERO(WORK(KDVEC),NVIR(ISYMD))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYMD) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMD),CMO(KOFF1),NBAS(ISYDEL),WORK(KDVEC),1)
C
C--------------------------------------------------------
C     Outer symmetry loop and next work space allocation.
C--------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMG  = ISYMK
         ISYMDK = MULD2H(ISYMD,ISYMK)
         ISYMCJ = ISYMDK
         ISALBE = ISYMCJ 
C
         KINAOK = KEND1
         KSCR1  = KINAOK + NNBST(ISALBE)*NRHF(ISYMK)
         KSCR2  = KSCR1  + N2BST(ISALBE)
         KEND2  = KSCR2  + NT1FRO(ISYMCJ)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient space for allocation in CC_FRCOIN')
         ENDIF
C
         CALL DZERO(WORK(KINAOK),NNBST(ISALBE)*NRHF(ISYMK))
C
C--------------------------------------------
C        Transform gamma-index to correlated.
C--------------------------------------------
C
         KOFF2  = IDSAOG(ISYMG,ISYDIS) + 1
         KOFF3  = ILRHSI(ISYMG) + NBAS(ISYMG)*NRHFFR(ISYMK) + 1
C
         NTOTAB = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHF(ISYMK),NBAS(ISYMG),ONE,
     *              XINT(KOFF2),NTOTAB,CMO(KOFF3),NTOTG,ZERO,
     *              WORK(KINAOK),NTOTAB)


         DO 110 K = 1,NRHF(ISYMK)
C
            CALL DZERO(WORK(KSCR1),N2BST(ISALBE))
            CALL DZERO(WORK(KSCR2),NT1FRO(ISYMCJ))
C
C-----------------------------------
C           Square up the integrals.
C-----------------------------------
C
            KOFF4 = KINAOK + NNBST(ISALBE)*(K - 1)
            CALL CCSD_SYMSQ(WORK(KOFF4),ISALBE,WORK(KSCR1))
C
C---------------------------------------------------------------
C           Inner symmetry loop and final work space allocation.
C---------------------------------------------------------------
C
            DO 120 ISYMAL = 1,NSYM
C
               ISYMC  = ISYMAL
               ISYMBE = MULD2H(ISALBE,ISYMAL)
               ISYMJ  = ISYMBE
C
               KSCR3 = KEND2
               KEND3 = KSCR3 + NBAS(ISYMAL)*NRHFFR(ISYMJ)
               LWRK3 = LWORK - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
                  CALL QUIT('Insufficient space for allocation '//
     &                      'in CC_FRCOIN')
               ENDIF
C
C------------------------------------------------
C              Construct the integrals (cJ|kdel).
C------------------------------------------------
C
               KOFF5  = KSCR1 + IAODIS(ISYMAL,ISYMBE)
               KOFF6  = ILRHSI(ISYMBE) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF5),NTOTAL,
     *                    CMO(KOFF6),NTOTBE,ZERO,WORK(KSCR3),NTOTAL)
C
               KOFF7  = ILVISI(ISYMAL) + 1
               KOFF8  = KSCR2 + IT1FRO(ISYMC,ISYMJ)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHFFR(ISYMJ),
     *                    NBAS(ISYMAL),ONE,CMO(KOFF7),NTOTAL,
     *                    WORK(KSCR3),NTOTAL,ZERO,WORK(KOFF8),NTOTC)
C
  120       CONTINUE
C
C----------------------------------------------------------------
C           Final scaling with CMO element and storage in result.
C----------------------------------------------------------------
C
            IF (R12TRA .AND. .NOT. R12PRP) THEN
               DO 131 D = 1,NRHFFR(ISYMD)
C
                  NDK = IF1FRO(ISYMD,ISYMK) + NRHFFR(ISYMD)*(K - 1) + D
                  KOFF9 = KDVEC + D - 1
                  KOFF10 = IF2FRO(ISYMCJ,ISYMDK)
     *                   + NT1FRO(ISYMCJ)*(NDK - 1) + 1
C
                  CALL DAXPY(NT1FRO(ISYMCJ),WORK(KOFF9),WORK(KSCR2),1,
     *                       XFRIN(KOFF10),1)
C
  131          CONTINUE
C
            ELSE
               DO 130 D = 1,NVIR(ISYMD)
C
                  NDK = IT1AM(ISYMD,ISYMK) + NVIR(ISYMD)*(K - 1) + D
                  KOFF9 = KDVEC + D - 1
                  KOFF10 = IT2FRO(ISYMCJ,ISYMDK)
     *                   + NT1FRO(ISYMCJ)*(NDK - 1) + 1
C
                  CALL DAXPY(NT1FRO(ISYMCJ),WORK(KOFF9),WORK(KSCR2),1,
     *                       XFRIN(KOFF10),1)
  130          CONTINUE
            END IF
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcogr */
      SUBROUTINE CC_FRCOGR(XFRIN,XINT,CMO,WORK,LWORK,IDEL,ISYDEL)
C
C
C     To calculate the integrals (pJ|Ik) needed for frozen-core gradients.
C     AO integrals are assumed totally symmetric.
c     based on cc_frcoin (Elena Vollmer 29/09/04) 
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XFRIN(*), XINT(*), CMO(*), WORK(LWORK)
      integer  ilorbi(8),iv2fro(8,8),iv1fro(8,8)
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <r12int.h>
C
C---------------------------------------
C     Initial symmetries and allocation.
C---------------------------------------
C
celena
      ISYDIS = ISYDEL
      ISYMD  = ISYDEL


      noffset = 0
      icoun1  = 0
      do isym = 1,nsym
         ilorbi(isym) = icoun1 
         icoun1  =  icoun1 + nbas(isym)*(nvir(isym)-nrhffr(isym))
         noffset = noffset + nbas(isym)*(nrhffr(isym)) 
      enddo

      do isymk = 1, nsym
         icoun3 = 0
         icoun4 = 0
         do isymj = 1,nsym
            isymi  = muld2h(isymk,isymj)
            iv2fro(isymi,isymj) = icoun3
            iv1fro(isymi,isymj) = icoun4 
            icoun3 = icoun3 + nt1fro(isymi)*nv1fro(isymj)
            icoun4 = icoun4 + (nrhffr(isymi))
     &                   *(norb1(isymj)-nrhffr(isymj))
         enddo
      enddo

celena

C
      KDVEC  = 1
      KEND1  = KDVEC + NVIR(ISYMD)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CC_FRCOIN')
      ENDIF
C
      CALL DZERO(WORK(KDVEC),NVIR(ISYMD))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYMD) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMD),CMO(KOFF1),NBAS(ISYDEL),WORK(KDVEC),1)
C
C--------------------------------------------------------
C     Outer symmetry loop and next work space allocation.
C--------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMG  = ISYMK
         ISYMDK = MULD2H(ISYMD,ISYMK)
         ISYMCJ = ISYMDK
         ISALBE = ISYMCJ 
C
         KINAOK = KEND1
         KSCR1  = KINAOK + NNBST(ISALBE)*(NORB1(ISYMK)-NRHFFR(ISYMK))
         KSCR2  = KSCR1  + N2BST(ISALBE)
         KEND2  = KSCR2  + NT1FRO(ISYMCJ)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient space for allocation in CC_FRCOGR')
         ENDIF
C
         CALL DZERO(WORK(KINAOK),NNBST(ISALBE)*(NORB1(ISYMK)
     &              -NRHFFR(ISYMK)))
C
C--------------------------------------------
C        Transform gamma-index to correlated.
C--------------------------------------------
C
         KOFF2  = IDSAOG(ISYMG,ISYDIS) + 1
         KOFF3  = NBAS(ISYMG)*NRHFFR(ISYMK)+ILVISI(ISYMG) + 1
C
         NTOTAB = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),(NORB1(ISYMK)-NRHFFR(ISYMK)),
     &              NBAS(ISYMG),ONE,
     *              XINT(KOFF2),NTOTAB,CMO(KOFF3),NTOTG,ZERO,
     *              WORK(KINAOK),NTOTAB)
         


         DO 110 K = 1,(NORB1(ISYMK)-NRHFFR(ISYMK))
C
            CALL DZERO(WORK(KSCR1),N2BST(ISALBE))
            CALL DZERO(WORK(KSCR2),NT1FRO(ISYMCJ))
C
C-----------------------------------
C           Square up the integrals.
C-----------------------------------
C
            KOFF4 = KINAOK + NNBST(ISALBE)*(K - 1)
            CALL CCSD_SYMSQ(WORK(KOFF4),ISALBE,WORK(KSCR1))
C
C---------------------------------------------------------------
C           Inner symmetry loop and final work space allocation.
C---------------------------------------------------------------
C
            DO 120 ISYMAL = 1,NSYM
C
               ISYMC  = ISYMAL
               ISYMBE = MULD2H(ISALBE,ISYMAL)
               ISYMJ  = ISYMBE
C
               KSCR3 = KEND2
               KEND3 = KSCR3 + NBAS(ISYMAL)*NRHFFR(ISYMJ)
               LWRK3 = LWORK - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
                  CALL QUIT('Insufficient space for allocation '//
     &                      'in CC_FRCOIN')
               ENDIF
C
C------------------------------------------------
C              Construct the integrals (cJ|kdel).
C------------------------------------------------
C
               KOFF5  = KSCR1 + IAODIS(ISYMAL,ISYMBE)
               KOFF6  = ILRHSI(ISYMBE) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF5),NTOTAL,
     *                    CMO(KOFF6),NTOTBE,ZERO,WORK(KSCR3),NTOTAL)
C

               KOFF7  = ILVISI(ISYMAL) + 1
               KOFF8  = KSCR2 + IT1FRO(ISYMC,ISYMJ)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHFFR(ISYMJ),
     *                    NBAS(ISYMAL),ONE,CMO(KOFF7),NTOTAL,
     *                    WORK(KSCR3),NTOTAL,ZERO,WORK(KOFF8),NTOTC)
C
  120       CONTINUE
C
C----------------------------------------------------------------
C           Final scaling with CMO element and storage in result.
C----------------------------------------------------------------
C
               DO 131 D = 1,NRHFFR(ISYMD)
C
                  NDK = IV1FRO(ISYMD,ISYMK) + NRHFFR(ISYMD)*(K - 1) + D
                  KOFF9 = KDVEC + D - 1
                  KOFF10 = IV2FRO(ISYMCJ,ISYMDK)
     *                   + NT1FRO(ISYMCJ)*(NDK - 1) + 1
C
                  CALL DAXPY(NT1FRO(ISYMCJ),WORK(KOFF9),WORK(KSCR2),1,
     *                       XFRIN(KOFF10),1)
C
  131          CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcogr1 */
      SUBROUTINE CC_FRCOGR1(XFRIN,XINT,CMO,WORK,LWORK,IDEL,ISYDEL)
C
C
C     To calculate the integrals (pJ|IK) needed for frozen-core gradients.
C     AO integrals are assumed totally symmetric.
c     based on cc_frcoin (Elena Vollmer 29/09/04) 
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XFRIN(*), XINT(*), CMO(*), WORK(LWORK)
      integer  iglmrf(8,8),ilorbi(8),iv2fro(8,8),iv1fro(8,8)
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <r12int.h>
C
C---------------------------------------
C     Initial symmetries and allocation.
C---------------------------------------
C
celena
      ISYDIS = ISYDEL
      ISYMD  = ISYDEL

      do isymmua = 1,nsym
         icou2 = 0
         do isyma = 1,nsym
            isymmu = muld2h(isymmua,isyma)
            iglmrf(isymmu,isyma) = icou2
            icou2 = icou2 + nbas(isymmu)*norb1(isyma)
         enddo
      enddo


      noffset = 0
      icoun1  = 0
      do isym = 1,nsym
         ilorbi(isym) = icoun1 
         icoun1  =  icoun1 + nbas(isym)*(nrhfs(isym))
         noffset = noffset + nbas(isym)*(nrhffr(isym)) 
      enddo

      do isymk = 1, nsym
         icoun3 = 0
         icoun4 = 0
         do isymj = 1,nsym
            isymi  = muld2h(isymk,isymj)
            iv2fro(isymi,isymj) = icoun3
            iv1fro(isymi,isymj) = icoun4 
            icoun3 = icoun3 + nt1fro(isymi)*nfrofr(isymj)
            icoun4 = icoun4 + (nrhffr(isymi))
     &                   *(norb1(isymj)-nrhffr(isymj))
         enddo
      enddo

celena

C
      KDVEC  = 1
      KEND1  = KDVEC + NVIR(ISYMD)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CC_FRCOIN')
      ENDIF
C
      CALL DZERO(WORK(KDVEC),NVIR(ISYMD))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYMD) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMD),CMO(KOFF1),NBAS(ISYDEL),WORK(KDVEC),1)
C
C--------------------------------------------------------
C     Outer symmetry loop and next work space allocation.
C--------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMG  = ISYMK
         ISYMDK = MULD2H(ISYMD,ISYMK)
         ISYMCJ = ISYMDK
         ISALBE = ISYMCJ 
C
         KINAOK = KEND1
         KSCR1  = KINAOK + NNBST(ISALBE)*NRHFFR(ISYMK)
         KSCR2  = KSCR1  + N2BST(ISALBE)
         KEND2  = KSCR2  + NT1FRO(ISYMCJ)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient space for allocation in CC_FRCOGR')
         ENDIF
C
         CALL DZERO(WORK(KINAOK),NNBST(ISALBE)*NRHFFR(ISYMK))
C
C--------------------------------------------
C        Transform gamma-index to correlated.
C--------------------------------------------
C
         KOFF2  = IDSAOG(ISYMG,ISYDIS) + 1
c         KOFF3  = 1 + iglmrf(ISYMg,ISYMk)
         KOFF3  = 1 + ilorbi(ISYMg)
C
         NTOTAB = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHFFR(ISYMK),
     &              NBAS(ISYMG),ONE,
     *              XINT(KOFF2),NTOTAB,CMO(KOFF3),NTOTG,ZERO,
     *              WORK(KINAOK),NTOTAB)


         DO 110 K = 1,NRHFFR(ISYMK)
C
            CALL DZERO(WORK(KSCR1),N2BST(ISALBE))
            CALL DZERO(WORK(KSCR2),NT1FRO(ISYMCJ))
C
C-----------------------------------
C           Square up the integrals.
C-----------------------------------
C
            KOFF4 = KINAOK + NNBST(ISALBE)*(K - 1)
            CALL CCSD_SYMSQ(WORK(KOFF4),ISALBE,WORK(KSCR1))
C
C---------------------------------------------------------------
C           Inner symmetry loop and final work space allocation.
C---------------------------------------------------------------
C
            DO 120 ISYMAL = 1,NSYM
C
               ISYMC  = ISYMAL
               ISYMBE = MULD2H(ISALBE,ISYMAL)
               ISYMJ  = ISYMBE
C
               KSCR3 = KEND2
               KEND3 = KSCR3 + NBAS(ISYMAL)*NRHFFR(ISYMJ)
               LWRK3 = LWORK - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
                  CALL QUIT('Insufficient space for allocation '//
     &                      'in CC_FRCOIN')
               ENDIF
C
C------------------------------------------------
C              Construct the integrals (cJ|kdel).
C------------------------------------------------
C
               KOFF5  = KSCR1 + IAODIS(ISYMAL,ISYMBE)
               KOFF6  = ILRHSI(ISYMBE) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF5),NTOTAL,
     *                    CMO(KOFF6),NTOTBE,ZERO,WORK(KSCR3),NTOTAL)
C

               KOFF7  = ILVISI(ISYMAL) + 1
               KOFF8  = KSCR2 + IT1FRO(ISYMC,ISYMJ)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHFFR(ISYMJ),
     *                    NBAS(ISYMAL),ONE,CMO(KOFF7),NTOTAL,
     *                    WORK(KSCR3),NTOTAL,ZERO,WORK(KOFF8),NTOTC)
C
  120       CONTINUE
C
C----------------------------------------------------------------
C           Final scaling with CMO element and storage in result.
C----------------------------------------------------------------
C
               DO 131 D = 1,NRHFFR(ISYMD)
C
                  NDK = IFROFR(ISYMD,ISYMK) + NRHFFR(ISYMD)*(K - 1) + D
                  KOFF9 = KDVEC + D - 1
                  KOFF10 = Iv2FRO(ISYMCJ,ISYMDK)
     *                   + NT1FRO(ISYMCJ)*(NDK - 1) + 1
C
                  CALL DAXPY(NT1FRO(ISYMCJ),WORK(KOFF9),WORK(KSCR2),1,
     *                       XFRIN(KOFF10),1)
C
  131          CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcr12 */
      SUBROUTINE CC_FRCR12(XFRIN,XINT,CMO,WORK,LWORK,IDEL,ISYDEL)
C
C     R12 version of CC_FRCOIN (WK/UniKA/20-08-2003)
C
C     Transform (ab|cd) into (p'J|kL), with J and L frozen.
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XFRIN(*), XINT(*), CMO(*), WORK(LWORK)
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <r12int.h>
C
C---------------------------------------
C     Initial symmetries and allocation.
C---------------------------------------
C
      ISYDIS = ISYDEL
      ISYMD  = ISYDEL
C
      KDVEC  = 1
      KEND1  = KDVEC + NVIR(ISYMD)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CC_FRCR12')
      ENDIF
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYMD) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMD),CMO(KOFF1),NBAS(ISYDEL),WORK(KDVEC),1)
C
C--------------------------------------------------------
C     Outer symmetry loop and next work space allocation.
C--------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMG  = ISYMK
         ISYMDK = MULD2H(ISYMD,ISYMK)
         ISYMCJ = ISYMDK
         ISALBE = ISYMCJ 
C
         KINAOK = KEND1
         KSCR1  = KINAOK + NNBST(ISALBE)*NRHFFR(ISYMK)
         KSCR2  = KSCR1  + N2BST(ISALBE)
         KEND2  = KSCR2  + NF1FRO(ISYMCJ)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient space for allocation in CC_FRCR12')
         ENDIF
C
         CALL DZERO(WORK(KINAOK),NNBST(ISALBE)*NRHFFR(ISYMK))
C
C--------------------------------------------
C        Transform gamma-index to frozen.
C--------------------------------------------
C
         KOFF2  = IDSAOG(ISYMG,ISYDIS) + 1
         KOFF3  = ILRHSI(ISYMG) + 1
C
         NTOTAB = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHFFR(ISYMK),NBAS(ISYMG),ONE,
     *              XINT(KOFF2),NTOTAB,CMO(KOFF3),NTOTG,ZERO,
     *              WORK(KINAOK),NTOTAB)
C
         DO 110 K = 1,NRHFFR(ISYMK)
C
            CALL DZERO(WORK(KSCR1),N2BST(ISALBE))
            CALL DZERO(WORK(KSCR2),NF1FRO(ISYMCJ))
C
C-----------------------------------
C           Square up the integrals.
C-----------------------------------
C
            KOFF4 = KINAOK + NNBST(ISALBE)*(K - 1)
            CALL CCSD_SYMSQ(WORK(KOFF4),ISALBE,WORK(KSCR1))
C
C---------------------------------------------------------------
C           Inner symmetry loop and final work space allocation.
C---------------------------------------------------------------
C
            DO 120 ISYMAL = 1,NSYM
C
               ISYMI  = ISYMAL
               ISYMBE = MULD2H(ISALBE,ISYMAL)
               ISYMJ  = ISYMBE
C
               KSCR3 = KEND2
               KEND3 = KSCR3 + NBAS(ISYMAL)*NRHF(ISYMJ)
               LWRK3 = LWORK - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
                  CALL QUIT('Insufficient space for allocation '//
     &                      'in CC_FRCR12')
               ENDIF
C
C------------------------------------------------
C              Construct the integrals (iJ|Kd).
C------------------------------------------------
C
               KOFF5  = KSCR1 + IAODIS(ISYMAL,ISYMBE)
               KOFF6  = ILRHSI(ISYMBE) + NBAS(ISYMBE)*NRHFFR(ISYMJ) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF5),NTOTAL,
     *                    CMO(KOFF6),NTOTBE,ZERO,WORK(KSCR3),NTOTAL)
C
               KOFF7  = ILRHSI(ISYMAL) + 1
               KOFF8  = KSCR2 + IF1FRO(ISYMI,ISYMJ)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTI  = MAX(NRHFFR(ISYMI),1)
C
               CALL DGEMM('T','N',NRHFFR(ISYMI),NRHF(ISYMJ),
     *                    NBAS(ISYMAL),ONE,CMO(KOFF7),NTOTAL,
     *                    WORK(KSCR3),NTOTAL,ZERO,WORK(KOFF8),NTOTI)
C
  120       CONTINUE
C
C----------------------------------------------------------------
C           Final scaling with CMO element and storage in result.
C----------------------------------------------------------------
C
            DO 130 D = 1,NVIR(ISYMD)
               NDK = IT1FRO(ISYMD,ISYMK) + NVIR(ISYMD)*(K - 1) + D
               KOFF9 = KDVEC + D - 1
               KOFF10 = IF2FRO(ISYMCJ,ISYMDK) + NDK
               CALL DAXPY(NF1FRO(ISYMCJ),WORK(KOFF9),WORK(KSCR2),1,
     *                    XFRIN(KOFF10),NT1FRO(ISYMDK))
  130       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck fock_all */
      SUBROUTINE FOCK_ALL(FOCKD,WORK,LWORK)
C
C     Written by Asger Halkier 25/5 - 1998
C
C     To read in and change the symmetry ordering of the FULL
C     Fock matrix diagonal from Sirius to CC-ordering.
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION FOCKD(*),WORK(LWORK)
#include <priunit.h>
#include <inftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      IF (LWORK .LT. NORBTS) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NORBTS
         CALL QUIT('Insufficient work space for allocation in '//
     &             'FOCK_ALL')
      ENDIF
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUSIFC)
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(I), I = 1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C------------------------------
C     Do the actual reordering.
C------------------------------
C
      IDRHF  = 0
      IDVIR  = NRHFTS
      ICOUNT = 0
C
      DO 100 ISYM = 1,NSYM
         DO 110 I = 1,NRHFS(ISYM)
            IDRHF  = IDRHF  + 1
            ICOUNT = ICOUNT + 1
            FOCKD(IDRHF) = WORK(ICOUNT)
  110    CONTINUE
C
         DO 120 A = 1,NVIRS(ISYM)
            IDVIR  = IDVIR  + 1
            ICOUNT = ICOUNT + 1
            FOCKD(IDVIR) = WORK(ICOUNT)
  120    CONTINUE
  100 CONTINUE
C
      IF (IPRINT .GT. 20) THEN
         CALL AROUND('Fock matrix diagonal in FOCK_ALL')
         WRITE(LUPRI,1)
         DO 200 I = 1,NORBTS
            WRITE(LUPRI,2) WORK(I), FOCKD(I)
  200    CONTINUE
      END IF
C
      RETURN
C
    1 FORMAT(7X,'Sirius order',5X,'CC order')
    2 FORMAT(6X,F14.10,3X,F14.10)
C
      END
C  /* Deck cc_d1fcb */
      SUBROUTINE CC_D1FCB(AODEN,DIJ,DJI,DAI,DIA,WORK,LWORK,ISDEN,ICON)
C
C     Written by Asger Halkier 27/5 - 1998.
C
C     To calculate the contributions to the AO one electron
C     density AODEN from the frozen core blocks.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION AODEN(*), DIJ(*), DJI(*), DAI(*), DIA(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CHARACTER MODEL*5
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KCMO  = 1
      KEND1 = KCMO  + NLAMDS
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CC_D1FCB')
      ENDIF
C
C----------------------------------------
C     Get the FULL MO coefficient matrix.
C----------------------------------------
C
      CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
C-----------------------------------------------
C     First symmetry loop and second allocation.
C-----------------------------------------------
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISDEN,ISYMAL)
         ISYMI1 = ISYMAL
         ISYMJ1 = ISYMBE
         ISYMI2 = ISYMBE
         ISYMJ2 = ISYMAL
C
         KSCR1  = KEND1
         KSCR2  = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMJ1)
         KEND2  = KSCR2 + NBAS(ISYMBE)*NRHFFR(ISYMJ2)
         LWRK2  = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for allocation '//
     &                'in CC_D1FCB')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND2-KEND1)
C
C-----------------------------------------------------
C        Calculate the contribution from the iJ block.
C-----------------------------------------------------
C
         KOFF1  = KCMO + ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMI1)
         KOFF2  = ICOFRO(ISYMI1,ISYMJ1) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTI  = MAX(NRHF(ISYMI1),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ1),NRHF(ISYMI1),
     *              ONE,WORK(KOFF1),NTOTAL,DIJ(KOFF2),NTOTI,ZERO,
     *              WORK(KSCR1),NTOTAL)
C
         KOFF3  = KCMO + ILRHSI(ISYMBE)
         KOFF4  = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMJ1),
     *              ONE,WORK(KSCR1),NTOTAL,WORK(KOFF3),NTOTBE,ONE,
     *              AODEN(KOFF4),NTOTAL)
C
C-----------------------------------------------------
C        Calculate the contribution from the Ji block.
C-----------------------------------------------------
C
         KOFF5  = KCMO + ILRHSI(ISYMBE) + NBAS(ISYMBE)*NRHFFR(ISYMI2)
         KOFF6  = ICOFRO(ISYMI2,ISYMJ2) + 1
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTI  = MAX(NRHF(ISYMI2),1)
C
         CALL DGEMM('N','N',NBAS(ISYMBE),NRHFFR(ISYMJ2),NRHF(ISYMI2),
     *              ONE,WORK(KOFF5),NTOTBE,DJI(KOFF6),NTOTI,ZERO,
     *              WORK(KSCR2),NTOTBE)
C
         KOFF7  = KCMO + ILRHSI(ISYMAL)
         KOFF8  = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMJ2),
     *              ONE,WORK(KOFF7),NTOTAL,WORK(KSCR2),NTOTBE,ONE,
     *              AODEN(KOFF8),NTOTAL)
C
  100 CONTINUE
C
C-----------------------------------------------
C     Second symmetry loop and final allocation.
C-----------------------------------------------
C
      DO 110 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISDEN,ISYMAL)
         ISYMA1 = ISYMAL
         ISYMI1 = ISYMBE
         ISYMA2 = ISYMBE
         ISYMI2 = ISYMAL
C
         KSCR1  = KEND1
         KSCR2  = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMI1)
         KEND3  = KSCR2 + NBAS(ISYMBE)*NRHFFR(ISYMI2)
         LWRK3  = LWORK - KEND3
C
         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
            CALL QUIT('Insufficient work space for allocation '//
     &                'in CC_D1FCB')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND3-KEND1)
C
C-----------------------------------------------------
C        Calculate the contribution from the aI block.
C-----------------------------------------------------
C
         KOFF9  = KCMO + ILVISI(ISYMAL)
         KOFF10 = IT1FRO(ISYMA1,ISYMI1) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA1),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMI1),NVIR(ISYMA1),
     *              ONE,WORK(KOFF9),NTOTAL,DAI(KOFF10),NTOTA,ZERO,
     *              WORK(KSCR1),NTOTAL)
C
         KOFF11 = KCMO + ILRHSI(ISYMBE)
         KOFF12 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMI1),
     *              ONE,WORK(KSCR1),NTOTAL,WORK(KOFF11),NTOTBE,ONE,
     *              AODEN(KOFF12),NTOTAL)
C
C-----------------------------------------------------
C        Calculate the contribution from the Ia block.
C-----------------------------------------------------
C
         KOFF13 = KCMO + ILVISI(ISYMBE)
         KOFF14 = IT1FRO(ISYMA2,ISYMI2) + 1
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTA  = MAX(NVIR(ISYMA2),1)
C
         CALL DGEMM('N','N',NBAS(ISYMBE),NRHFFR(ISYMI2),NVIR(ISYMA2),
     *              ONE,WORK(KOFF13),NTOTBE,DIA(KOFF14),NTOTA,ZERO,
     *              WORK(KSCR2),NTOTBE)
C
         KOFF15 = KCMO + ILRHSI(ISYMAL)
         KOFF16 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMI2),
     *              ONE,WORK(KOFF15),NTOTAL,WORK(KSCR2),NTOTBE,ONE,
     *              AODEN(KOFF16),NTOTAL)
C
  110 CONTINUE
C
C----------------------------------------
C     Add terms from frozen-frozen block.
C----------------------------------------
C
      IF (ICON .EQ. 1) THEN
         MODEL = 'DUMMY'
         CALL CC_FCD1AO(AODEN,WORK(KEND1),LWRK1,MODEL)
      ENDIF
C
      RETURN
      END
C  /* Deck cc_gtofro */
      SUBROUTINE CC_GTOFRO(XINT,DSFRO,CMO,WORK,LWORK,ISYDIS)
C
C     Written by Asger Halkier 28/5 - 1998.
C
C     To calculate the integral batch (al be | fro del).
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*), DSFRO(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      DO 100 ISYMG = 1,NSYM
C
         ISYMI  = ISYMG
         ISALBE = MULD2H(ISYMG,ISYDIS)
C
         KOFF1  = IDSAOG(ISYMG,ISYDIS) + 1
         KOFF2  = ILRHSI(ISYMG) + 1
         KOFF3  = IDSFRO(ISALBE,ISYMI) + 1
C
         NTOTAB = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHFFR(ISYMI),NBAS(ISYMG),
     *              ONE,XINT(KOFF1),NTOTAB,CMO(KOFF2),NTOTG,ZERO,
     *              DSFRO(KOFF3),NTOTAB)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_ofroin */
      SUBROUTINE CC_OFROIN(DSRHF,DSRES,CMO,WORK,LWORK,ISYDIS)
C
C     Written by Asger Halkier 28/5 - 1998.
C
C     To calculate the integral batch (cor fro | cor del).
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DSRHF(*), DSRES(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C-------------------------------------------------------
C     Outer symmetry loop and work space allocation one.
C-------------------------------------------------------
C
      DO 100 ISYML = 1,NSYM
C
         ISALBE = MULD2H(ISYDIS,ISYML)
         ISYMKI = MULD2H(ISYDIS,ISYML)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KAOINT = 1
         KEND1  = KAOINT + N2BST(ISALBE)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for allocation '//
     &                'in CC_OFROIN')
         ENDIF
C
         CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
C
         DO 110 L = 1,NRHF(ISYML)
C
C----------------------------------------
C           Unpack integral distribution.
C----------------------------------------
C
            KOFF1 = IDSRHF(ISALBE,ISYML) + NNBST(ISALBE)*(L - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KAOINT))
C
C-------------------------------------------------------------
C           Inner symmetry loop and work space allocation two.
C-------------------------------------------------------------
C
            DO 120 ISYMAL = 1,NSYM
C
               ISYMBE = MULD2H(ISYMAL,ISALBE)
               ISYMK  = ISYMAL
               ISYMI  = ISYMBE
C
               KSCR1 = KEND1
               KEND2 = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMI)
               LWRK2 = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
                  CALL QUIT('Insufficient space for allocation '//
     &                      'in CC_OFROIN')
               ENDIF
C
               CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NRHFFR(ISYMI))
C
C----------------------------------------------------------------
C              Perform the contractions for obtaining the result.
C----------------------------------------------------------------
C
               KOFF2  = KAOINT + IAODIS(ISYMAL,ISYMBE)
               KOFF3  = ILRHSI(ISYMBE) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMI),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     *                    CMO(KOFF3),NTOTBE,ZERO,WORK(KSCR1),NTOTAL)
C
               KOFF4  = ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMK) + 1
               KOFF5  = IOFROO(ISYMKI,ISYML) + NCOFRO(ISYMKI)*(L - 1)
     *                + ICOFRO(ISYMK,ISYMI)  + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NRHFFR(ISYMI),
     *                    NBAS(ISYMAL),ONE,CMO(KOFF4),NTOTAL,
     *                    WORK(KSCR1),NTOTAL,ZERO,DSRES(KOFF5),NTOTK)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_etfrdi */
      SUBROUTINE MP2_ETFRD(ETAFRO,DSINT,TINT,ISYDIS)
C
C     Written by Asger Halkier 28/5 - 1998.
C
C     To calculate the direct contribution to the frozen part of
C     eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAFRO(*), DSINT(*), TINT(*)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 ISYML = 1,NSYM
C
         ISYMKI = MULD2H(ISYDIS,ISYML)
         ISYMAK = MULD2H(ISYDIS,ISYML)
C
         DO 110 L = 1,NRHF(ISYML)
            DO 120 ISYMA = 1,NSYM
C
               ISYMK = MULD2H(ISYMAK,ISYMA)
               ISYMI = ISYMA
C
               KOFF1 = IT2BCD(ISYMAK,ISYML) + NT1AM(ISYMAK)*(L - 1)
     *               + IT1AM(ISYMA,ISYMK)   + 1
               KOFF2 = IOFROO(ISYMKI,ISYML) + NCOFRO(ISYMKI)*(L - 1)
     *               + ICOFRO(ISYMK,ISYMI)  + 1
               KOFF3 = IT1FRO(ISYMA,ISYMI)  + 1
C
               NTOTA = MAX(NVIR(ISYMA),1)
               NTOTK = MAX(NRHF(ISYMK),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),
     *                    NRHF(ISYMK),-ONE,TINT(KOFF1),NTOTA,
     *                    DSINT(KOFF2),NTOTK,ONE,ETAFRO(KOFF3),NTOTA)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidv1 */
      SUBROUTINE MP2_EIDV1(ETAFRO,XINTFR,ZKAB,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 29/5 - 1998.
C
C     To calculate the coulomb part of the indirect virtual
C     contribution to the frozen part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      DIMENSION ETAFRO(*), XINTFR(*), ZKAB(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMBC = 1
      ISALBE = ISYMBC
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
      KAVEC  = 1
      KAOINT = KAVEC  + NVIR(ISYMA)
      KMOINT = KAOINT + N2BST(ISALBE)
      KEND1  = KMOINT + NMATAB(ISYMBC)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDV1')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
      DO 100 I = 1,NRHFFR(ISYMI)
C
         CALL DZERO(WORK(KMOINT),NMATAB(ISYMBC))
C
C--------------------------------------------------------
C        Square up integral distribution (al be | I del).
C--------------------------------------------------------
C
         KOFF2 = IDSFRO(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
         CALL CCSD_SYMSQ(XINTFR(KOFF2),ISALBE,WORK(KAOINT))
C
C------------------------------------------------------
C        Inner symmetry loop and final work allocation.
C------------------------------------------------------
C
         DO 110 ISYMB = 1,NSYM
C
            ISYMC  = MULD2H(ISYMBC,ISYMB)
            ISYMAL = ISYMB
            ISYMBE = ISYMC
C
            KSCR1 = KEND1
            KEND2 = KSCR1 + NBAS(ISYMAL)*NVIR(ISYMC)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT('Insufficient memory for allocation '//
     &                   'in MP2_EIDV1')
            ENDIF
C
            CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NVIR(ISYMC))
C
C---------------------------------------------
C           Calculate integrals (b c | I del).
C---------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILVISI(ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMC),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,CMO(KOFF4),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF5  = ILVISI(ISYMAL) + 1
            KOFF6  = KMOINT + IMATAB(ISYMB,ISYMC)
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTB  = MAX(NVIR(ISYMB),1)
C
            CALL DGEMM('T','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF5),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KOFF6),NTOTB)
C
  110    CONTINUE
C
C-----------------------------------
C        Calculate the contribution.
C-----------------------------------
C
         FACT  = FOUR*DDOT(NMATAB(ISYMBC),WORK(KMOINT),1,ZKAB(1),1)
         KOFF7 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
         CALL DAXPY(NVIR(ISYMA),FACT,WORK(KAVEC),1,ETAFRO(KOFF7),1)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidv2 */
      SUBROUTINE MP2_EIDV2(ETAFRO,XINTFR,ZKAB,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 29/5 - 1998.
C
C     To calculate the exchange part of the indirect virtual
C     contribution to the frozen part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAFRO(*), XINTFR(*), ZKAB(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C------------------------------------------------
C     Initial symmetries & work space allocation.
C------------------------------------------------
C
      ISYMAI = 1
      ISYMBC = 1
      ISYMC  = ISYDEL
      ISYMB  = MULD2H(ISYMBC,ISYMC)
C
      KCVEC  = 1
      KBVEC  = KCVEC + NVIR(ISYMC)
      KEND1  = KBVEC + NVIR(ISYMB)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDV2')
      ENDIF
C
      CALL DZERO(WORK(KCVEC),NVIR(ISYMC))
      CALL DZERO(WORK(KBVEC),NVIR(ISYMB))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFFC = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMC),CMO(KOFFC),NBAS(ISYDEL),WORK(KCVEC),1)
C
C-----------------------------------
C     Contract ZKAB with CMO vector.
C-----------------------------------
C
      KOFF1 = IMATAB(ISYMB,ISYMC) + 1 
C
      NTOTB = MAX(NVIR(ISYMB),1)
C
      CALL DGEMV('N',NVIR(ISYMB),NVIR(ISYMC),ONE,ZKAB(KOFF1),NTOTB,
     *           WORK(KCVEC),1,ZERO,WORK(KBVEC),1)
C
C----------------------------------------------------
C     Symmetry loop and second work space allocation.
C----------------------------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         IF (NRHFFR(ISYMI) .EQ. 0) GOTO 100
C
         ISYMA  = MULD2H(ISYMAI,ISYMI)
         ISYMAL = ISYMA
         ISYMBE = ISYMB
         ISALBE = MULD2H(ISYMAL,ISYMBE)
C
         KAOINT = KEND1
         KSCR1  = KAOINT + N2BST(ISALBE)
         KSCR2  = KSCR1  + NBAS(ISYMAL)*NVIR(ISYMB)
         KEND2  = KSCR2  + NVIR(ISYMA)*NVIR(ISYMB)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for allocation '//
     &                'in MP2_EIDV2')
         ENDIF
C
         CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
         CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NVIR(ISYMB))
         CALL DZERO(WORK(KSCR2),NVIR(ISYMA)*NVIR(ISYMB))
C
         DO 110 I = 1,NRHFFR(ISYMI)
C
C----------------------------------------
C           Unpack integral distribution.
C----------------------------------------
C
            KOFF2 = IDSFRO(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
            CALL CCSD_SYMSQ(XINTFR(KOFF2),ISALBE,WORK(KAOINT))
C
C-------------------------------------------
C           Transform integrals to MO basis.
C-------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILVISI(ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMB),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,CMO(KOFF4),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF5  = ILVISI(ISYMAL) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTA  = MAX(NVIR(ISYMA),1)
C
            CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF5),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KSCR2),NTOTA)
C
C--------------------------------------------------
C           Final matrix vector product for result.
C--------------------------------------------------
C
            KOFF6  = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
            NTOTA  = MAX(NVIR(ISYMA),1)
C
            CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMB),-TWO,WORK(KSCR2),
     *                 NTOTA,WORK(KBVEC),1,ONE,ETAFRO(KOFF6),1)
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidc1 */
      SUBROUTINE MP2_EIDC1(ETAFRO,XINTFR,ZKIJ,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 30/5 - 1998.
C
C     To calculate the coulomb part of the indirect correlated
C     contribution to the frozen part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      DIMENSION ETAFRO(*), XINTFR(*), ZKIJ(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
      ISALBE = ISYMJK
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
      KAVEC  = 1
      KAOINT = KAVEC  + NVIR(ISYMA)
      KMOINT = KAOINT + N2BST(ISALBE)
      KEND1  = KMOINT + NMATIJ(ISYMJK)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDC1')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
      DO 100 I = 1,NRHFFR(ISYMI)
C
         CALL DZERO(WORK(KMOINT),NMATIJ(ISYMJK))
C
C--------------------------------------------------------
C        Square up integral distribution (al be | I del).
C--------------------------------------------------------
C
         KOFF2 = IDSFRO(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
         CALL CCSD_SYMSQ(XINTFR(KOFF2),ISALBE,WORK(KAOINT))
C
C------------------------------------------------------
C        Inner symmetry loop and final work allocation.
C------------------------------------------------------
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMJK,ISYMJ)
            ISYMAL = ISYMJ
            ISYMBE = ISYMK
C
            KSCR1 = KEND1
            KEND2 = KSCR1 + NBAS(ISYMAL)*NRHF(ISYMK)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT('Insufficient memory for allocation '//
     &                   'in MP2_EIDC1')
            ENDIF
C
            CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NRHF(ISYMK))
C
C---------------------------------------------
C           Calculate integrals (j k | I del).
C---------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILRHSI(ISYMBE) + NBAS(ISYMBE)*NRHFFR(ISYMK) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,CMO(KOFF4),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF5  = ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMJ) + 1
            KOFF6  = KMOINT + IMATIJ(ISYMJ,ISYMK)
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTJ  = MAX(NRHF(ISYMJ),1)
C
            CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMK),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF5),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KOFF6),NTOTJ)
C
  110    CONTINUE
C
C-----------------------------------
C        Calculate the contribution.
C-----------------------------------
C
         FACT  = FOUR*DDOT(NMATIJ(ISYMJK),WORK(KMOINT),1,ZKIJ(1),1)
         KOFF7 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
         CALL DAXPY(NVIR(ISYMA),FACT,WORK(KAVEC),1,ETAFRO(KOFF7),1)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidc2 */
      SUBROUTINE MP2_EIDC2(ETAFRO,XINOFO,ZKIJ,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 30/5 - 1998.
C
C     To calculate the exchange part of the indirect correlated
C     contribution to the frozen part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAFRO(*), XINOFO(*), ZKIJ(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
      KAVEC = 1
      KIVEC = KAVEC  + NVIR(ISYMA)
      KEND1 = KIVEC  + NRHFFR(ISYMI)
      LWRK1 = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDC2')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KIVEC),NRHFFR(ISYMI))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C-------------------------------------------------
C     Outer loops over third integral batch index.
C-------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMJK,ISYMK)
         ISYMJI = MULD2H(ISYMJ,ISYMI)
C
         DO 110 K = 1,NRHF(ISYMK)
C
C-------------------------------------------------------
C           Inner loop over frozen index and contraction
C           of integrals and ZKIJ-vector.
C-------------------------------------------------------
C
            DO 120 I = 1,NRHFFR(ISYMI)
C
               NJK  = IMATIJ(ISYMJ,ISYMK)  + NRHF(ISYMJ)*(K - 1) + 1
               NJIK = IOFROO(ISYMJI,ISYMK) + NCOFRO(ISYMJI)*(K - 1)
     *              + ICOFRO(ISYMJ,ISYMI)  + NRHF(ISYMJ)*(I - 1) + 1
C
               FACT = DDOT(NRHF(ISYMJ),ZKIJ(NJK),1,XINOFO(NJIK),1)
               WORK(KIVEC+I-1) = WORK(KIVEC+I-1) + FACT
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------
C     Final storage in result.
C-----------------------------
C
      DO 130 I = 1,NRHFFR(ISYMI)
C
         KOFF2 = KIVEC + I - 1
         KOFF3 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),-TWO*WORK(KOFF2),WORK(KAVEC),1,
     *              ETAFRO(KOFF3),1)
C
  130 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidf1 */
      SUBROUTINE MP2_EIDF1(ETAAI,XINOFO,ZKJK,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 30/5 - 1998.
C
C     To calculate the coulomb part of the indirect frozen
C     contribution to the correlated part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, EIGHT = 8.0D0)
      DIMENSION ETAAI(*), XINOFO(*), ZKJK(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHF(ISYMI) .EQ. 0) RETURN
C
      KAVEC = 1
      KEND1 = KAVEC  + NVIR(ISYMA)
      LWRK1 = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDF1')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C----------------------
C     Calculate result.
C----------------------
C
      DO 100 I = 1,NRHF(ISYMI)
C
         KOFF1 = IOFROO(ISYMJK,ISYMI) + NCOFRO(ISYMJK)*(I - 1) + 1
         KOFF2 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
         FACT  = EIGHT*DDOT(NCOFRO(ISYMJK),ZKJK(1),1,XINOFO(KOFF1),1)
         CALL DAXPY(NVIR(ISYMA),FACT,WORK(KAVEC),1,ETAAI(KOFF2),1)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidf2 */
      SUBROUTINE MP2_EIDF2(ETAFRO,XINTFR,ZKJK,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 30/5 - 1998.
C
C     To calculate the coulomb part of the indirect frozen
C     contribution to the frozen part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, EIGHT = 8.0D0)
      DIMENSION ETAFRO(*), XINTFR(*), ZKJK(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
      ISALBE = ISYMJK
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
      KAVEC  = 1
      KAOINT = KAVEC  + NVIR(ISYMA)
      KMOINT = KAOINT + N2BST(ISALBE)
      KEND1  = KMOINT + NCOFRO(ISYMJK)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDF2')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
      DO 100 I = 1,NRHFFR(ISYMI)
C
         CALL DZERO(WORK(KMOINT),NCOFRO(ISYMJK))
C
C--------------------------------------------------------
C        Square up integral distribution (al be | I del).
C--------------------------------------------------------
C
         KOFF2 = IDSFRO(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
         CALL CCSD_SYMSQ(XINTFR(KOFF2),ISALBE,WORK(KAOINT))
C
C------------------------------------------------------
C        Inner symmetry loop and final work allocation.
C------------------------------------------------------
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMJK,ISYMJ)
            ISYMAL = ISYMJ
            ISYMBE = ISYMK
C
            KSCR1 = KEND1
            KEND2 = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMK)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT('Insufficient memory for allocation '//
     &                   'in MP2_EIDF2')
            ENDIF
C
            CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NRHFFR(ISYMK))
C
C---------------------------------------------
C           Calculate integrals (j K | I del).
C---------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILRHSI(ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMK),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,CMO(KOFF4),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF5  = ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMJ) + 1
            KOFF6  = KMOINT + ICOFRO(ISYMJ,ISYMK)
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTJ  = MAX(NRHF(ISYMJ),1)
C
            CALL DGEMM('T','N',NRHF(ISYMJ),NRHFFR(ISYMK),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF5),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KOFF6),NTOTJ)
C
  110    CONTINUE
C
C-----------------------------------
C        Calculate the contribution.
C-----------------------------------
C
         FACT  = EIGHT*DDOT(NCOFRO(ISYMJK),WORK(KMOINT),1,ZKJK(1),1)
         KOFF7 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
         CALL DAXPY(NVIR(ISYMA),FACT,WORK(KAVEC),1,ETAFRO(KOFF7),1)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidf3 */
      SUBROUTINE MP2_EIDF3(ETAAI,XINOFO,ZKJK,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 30/5 - 1998.
C
C     To calculate the first exchange part of the indirect frozen
C     contribution to the correlated part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAAI(*), XINOFO(*), ZKJK(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHF(ISYMI) .EQ. 0) RETURN
C
      KAVEC = 1
      KIVEC = KAVEC + NVIR(ISYMA)
      KEND1 = KIVEC + NRHF(ISYMI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDF3')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KIVEC),NRHF(ISYMI))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C----------------------------------------------
C     Loops over third index in integral batch.
C----------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMK  = MULD2H(ISYMJK,ISYMJ)
         ISYMIK = MULD2H(ISYMI,ISYMK)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            KOFF2 = IOFROO(ISYMIK,ISYMJ) + NCOFRO(ISYMIK)*(J - 1)
     *            + ICOFRO(ISYMI,ISYMK)  + 1
            KOFF3 = ICOFRO(ISYMJ,ISYMK)  + J
C
            NTOTI = MAX(NRHF(ISYMI),1)
C
            CALL DGEMV('N',NRHF(ISYMI),NRHFFR(ISYMK),ONE,
     *                 XINOFO(KOFF2),NTOTI,ZKJK(KOFF3),NRHF(ISYMJ),
     *                 ONE,WORK(KIVEC),1)
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------
C     Final storage in result.
C-----------------------------
C
      DO 120 I = 1,NRHF(ISYMI)
C
         KOFF4 = KIVEC + I - 1
         KOFF5 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),-TWO*WORK(KOFF4),WORK(KAVEC),1,
     *              ETAAI(KOFF5),1)
C
  120 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidf4 */
      SUBROUTINE MP2_EIDF4(ETAAI,DSFRO,ZKJK,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 30/5 - 1998.
C
C     To calculate the second exchange part of the indirect frozen
C     contribution to the correlated part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAAI(*), DSFRO(*), ZKJK(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHF(ISYMI) .EQ. 0) RETURN
C
      KAVEC = 1
      KIVEC = KAVEC + NVIR(ISYMA)
      KEND1 = KIVEC + NRHF(ISYMI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDF4')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KIVEC),NRHF(ISYMI))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C---------------------------------------------------------------------------
C     Loops over third index in integral batch & final wok space allocation.
C---------------------------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMJK,ISYMK)
         ISYMIJ = MULD2H(ISYMI,ISYMJ)
         ISALBE = ISYMIJ
         ISYMAL = ISYMI
         ISYMBE = ISYMJ
C
         KAOINT = KEND1
         KSCR1  = KAOINT + N2BST(ISALBE)
         KMOINT = KSCR1  + NBAS(ISYMAL)*NRHF(ISYMJ)
         KEND2  = KMOINT + NRHF(ISYMI)*NRHF(ISYMJ)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for allocation '//
     &                'in MP2_EIDF4')
         ENDIF
C
         CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
C
         DO 110 K = 1,NRHFFR(ISYMK)
C
C-----------------------------------------------------------
C           Square up integral distribution (al be | K del).
C-----------------------------------------------------------
C
            KOFF2 = IDSFRO(ISALBE,ISYMK) + NNBST(ISALBE)*(K - 1) + 1
C
            CALL CCSD_SYMSQ(DSFRO(KOFF2),ISALBE,WORK(KAOINT))
C
C---------------------------------------------
C           Calculate integrals (i j | K del).
C---------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILRHSI(ISYMBE) + NBAS(ISYMBE)*NRHFFR(ISYMJ) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,CMO(KOFF4),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF5  = ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMI) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTI  = MAX(NRHF(ISYMI),1)
C
            CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF5),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KMOINT),NTOTI)
C
C-------------------------------------------
C           Contract MO integrals with ZKJK.
C-------------------------------------------
C
            KOFF6  = ICOFRO(ISYMJ,ISYMK) + NRHF(ISYMJ)*(K - 1) + 1
C
            NTOTI  = MAX(NRHF(ISYMI),1)
C
            CALL DGEMV('N',NRHF(ISYMI),NRHF(ISYMJ),ONE,WORK(KMOINT),
     *                 NTOTI,ZKJK(KOFF6),1,ONE,WORK(KIVEC),1)
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------
C     Final storage in result.
C-----------------------------
C
      DO 120 I = 1,NRHF(ISYMI)
C
         KOFF7 = KIVEC + I - 1
         KOFF8 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),-TWO*WORK(KOFF7),WORK(KAVEC),1,
     *              ETAAI(KOFF8),1)
C
  120 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidf6 */
      SUBROUTINE MP2_EIDF6(ETAFRO,DSFRO,ZKJK,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 31/5 - 1998.
C
C     To calculate the second exchange part of the indirect frozen
C     contribution to the frozen part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAFRO(*), DSFRO(*), ZKJK(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
      KAVEC = 1
      KIVEC = KAVEC + NVIR(ISYMA)
      KEND1 = KIVEC + NRHFFR(ISYMI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDF6')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KIVEC),NRHFFR(ISYMI))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C---------------------------------------------------------------------------
C     Loops over third index in integral batch & final wok space allocation.
C---------------------------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMJK,ISYMK)
         ISYMIJ = MULD2H(ISYMI,ISYMJ)
         ISALBE = ISYMIJ
         ISYMAL = ISYMI
         ISYMBE = ISYMJ
C
         KAOINT = KEND1
         KSCR1  = KAOINT + N2BST(ISALBE)
         KMOINT = KSCR1  + NBAS(ISYMAL)*NRHF(ISYMJ)
         KEND2  = KMOINT + NRHFFR(ISYMI)*NRHF(ISYMJ)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for allocation '//
     &                'in MP2_EIDF6')
         ENDIF
C
         CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
C
         DO 110 K = 1,NRHFFR(ISYMK)
C
C-----------------------------------------------------------
C           Square up integral distribution (al be | K del).
C-----------------------------------------------------------
C
            KOFF2 = IDSFRO(ISALBE,ISYMK) + NNBST(ISALBE)*(K - 1) + 1
C
            CALL CCSD_SYMSQ(DSFRO(KOFF2),ISALBE,WORK(KAOINT))
C
C---------------------------------------------
C           Calculate integrals (I j | K del).
C---------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILRHSI(ISYMBE) + NBAS(ISYMBE)*NRHFFR(ISYMJ) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,CMO(KOFF4),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF5  = ILRHSI(ISYMAL) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTI  = MAX(NRHFFR(ISYMI),1)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHF(ISYMJ),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF5),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KMOINT),NTOTI)
C
C-------------------------------------------
C           Contract MO integrals with ZKJK.
C-------------------------------------------
C
            KOFF6  = ICOFRO(ISYMJ,ISYMK) + NRHF(ISYMJ)*(K - 1) + 1
C
            NTOTI  = MAX(NRHFFR(ISYMI),1)
C
            CALL DGEMV('N',NRHFFR(ISYMI),NRHF(ISYMJ),ONE,WORK(KMOINT),
     *                 NTOTI,ZKJK(KOFF6),1,ONE,WORK(KIVEC),1)
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------
C     Final storage in result.
C-----------------------------
C
      DO 120 I = 1,NRHFFR(ISYMI)
C
         KOFF7 = KIVEC + I - 1
         KOFF8 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),-TWO*WORK(KOFF7),WORK(KAVEC),1,
     *              ETAFRO(KOFF8),1)
C
  120 CONTINUE
C
      RETURN
      END
C  /* Deck mp2_eidf5 */
      SUBROUTINE MP2_EIDF5(ETAFRO,DSRHF,ZKJK,CMO,WORK,LWORK,IDEL,
     *                     ISYDEL)
C
C     Written by Asger Halkier 31/5 - 1998.
C
C     To calculate the first exchange part of the indirect frozen
C     contribution to the frozen part of eta in MP2 calculations.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAFRO(*), DSRHF(*), ZKJK(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
      KAVEC = 1
      KIVEC = KAVEC + NVIR(ISYMA)
      KEND1 = KIVEC + NRHFFR(ISYMI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in MP2_EIDF5')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KIVEC),NRHFFR(ISYMI))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILVISI(ISYDEL) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C---------------------------------------------------------------------------
C     Loops over third index in integral batch & final wok space allocation.
C---------------------------------------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMK  = MULD2H(ISYMJK,ISYMJ)
         ISYMIK = MULD2H(ISYMI,ISYMK)
         ISALBE = ISYMIK
         ISYMAL = ISYMI
         ISYMBE = ISYMK
C
         IF (NRHFFR(ISYMK) .EQ. 0) GOTO 100
C
         KAOINT = KEND1
         KSCR1  = KAOINT + N2BST(ISALBE)
         KMOINT = KSCR1  + NBAS(ISYMAL)*NRHFFR(ISYMK)
         KEND2  = KMOINT + NRHFFR(ISYMI)*NRHFFR(ISYMK)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for allocation '//
     &                'in MP2_EIDF5')
         ENDIF
C
         CALL DZERO(WORK(KAOINT),KEND2-KEND1)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
C-----------------------------------------------------------
C           Square up integral distribution (al be | j del).
C-----------------------------------------------------------
C
            KOFF2 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE)*(J - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF2),ISALBE,WORK(KAOINT))
C
C---------------------------------------------
C           Calculate integrals (I K | j del).
C---------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILRHSI(ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMK),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,CMO(KOFF4),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF5  = ILRHSI(ISYMAL) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTI  = MAX(NRHFFR(ISYMI),1)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMK),
     *                 NBAS(ISYMAL),ONE,CMO(KOFF5),NTOTAL,WORK(KSCR1),
     *                 NTOTAL,ZERO,WORK(KMOINT),NTOTI)
C
C-------------------------------------------
C           Contract MO integrals with ZKJK.
C-------------------------------------------
C
            KOFF6  = ICOFRO(ISYMJ,ISYMK) + J
C
            NTOTI  = MAX(NRHFFR(ISYMI),1)
C
            CALL DGEMV('N',NRHFFR(ISYMI),NRHFFR(ISYMK),ONE,
     *                 WORK(KMOINT),NTOTI,ZKJK(KOFF6),NRHF(ISYMJ),
     *                 ONE,WORK(KIVEC),1)
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------
C     Final storage in result.
C-----------------------------
C
      DO 120 I = 1,NRHFFR(ISYMI)
C
         KOFF7 = KIVEC + I - 1
         KOFF8 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),-TWO*WORK(KOFF7),WORK(KAVEC),1,
     *              ETAFRO(KOFF8),1)
C
  120 CONTINUE
C
      RETURN
      END
C  /* Deck ccifromo */
      SUBROUTINE CCIFROMO(XFIJ,XFJI,XFAI,XFIA,XFII,XAOINT,
     *                    XLAMDP,XLAMDH,CMO,WORK,LWORK,ISYM)
C
C     Written by Asger Halkier 13/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To transform the AO integrals in (XAOINT) to
C              frozen MO-basis (stored in XINTIJ through XINTAI)!
C              Note that the AO integrals either can be the one
C              electron integrals or a submatrix alfa,beta of the two
C              electron integrals (alfa beta|gamma delta)!
C              ISYM is the integral symmetry!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XFIJ(*), XFJI(*), XFAI(*), XFIA(*), XFII(*)
      DIMENSION XAOINT(*), XLAMDP(*), XLAMDH(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C------------------------------
C     Initialize result arrays.
C------------------------------
C
      CALL DZERO(XFIJ,NCOFRO(ISYM))
      CALL DZERO(XFJI,NCOFRO(ISYM))
      CALL DZERO(XFAI,NT1FRO(ISYM))
      CALL DZERO(XFIA,NT1FRO(ISYM))
      CALL DZERO(XFII,NFROFR(ISYM))
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMAL,ISYM)
         ISYMJ  = ISYMBE
         ISYMI  = ISYMAL
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         LSCR = MAX(NBAS(ISYMAL)*NVIR(ISYMBE),
     *              NBAS(ISYMAL)*NRHFFR(ISYMBE))
C
         KSCR1 = 1
         KEND1 = KSCR1 + LSCR
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory in CCIFROMO')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
C--------------------------------------------------
C        Transform second integral index to frozen.
C--------------------------------------------------
C
         KOFF1  = IAODIS(ISYMAL,ISYMBE) + 1
         KOFF2  = ILRHSI(ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ),NBAS(ISYMBE),
     *              ONE,XAOINT(KOFF1),NTOTAL,CMO(KOFF2),NTOTBE,ZERO,
     *              WORK(KSCR1),NTOTAL)
C
C----------------------------------------------
C        Calculate the frozen-frozen integrals.
C----------------------------------------------
C
         KOFF3  = ILRHSI(ISYMAL) + 1
         KOFF4  = IFROFR(ISYMI,ISYMJ) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTI  = MAX(NRHFFR(ISYMI),1)
C
         CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NBAS(ISYMAL),
     *              ONE,CMO(KOFF3),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *              XFII(KOFF4),NTOTI)
C
C--------------------------------------------------
C        Calculate the correlated-frozen integrals.
C--------------------------------------------------
C
         KOFF5  = ILMRHF(ISYMAL) + 1
         KOFF6  = ICOFRO(ISYMI,ISYMJ) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NBAS(ISYMAL),
     *              ONE,XLAMDP(KOFF5),NTOTAL,WORK(KSCR1),NTOTAL,
     *              ZERO,XFIJ(KOFF6),NTOTI)
C
C-----------------------------------------------
C        Calculate the virtual-frozen integrals.
C-----------------------------------------------
C
         ISYMA  = ISYMAL
C
         KOFF7  = ILMVIR(ISYMAL) + 1
         KOFF8  = IT1FRO(ISYMA,ISYMJ) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHFFR(ISYMJ),NBAS(ISYMAL),
     *              ONE,XLAMDP(KOFF7),NTOTAL,WORK(KSCR1),NTOTAL,
     *              ZERO,XFAI(KOFF8),NTOTA)
C
C------------------------------------------------------
C        Transform second integral index to correlated.
C------------------------------------------------------
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
         KOFF9  = IAODIS(ISYMAL,ISYMBE) + 1
         KOFF10 = ILMRHF(ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NBAS(ISYMBE),
     *              ONE,XAOINT(KOFF9),NTOTAL,XLAMDH(KOFF10),NTOTBE,
     *              ZERO,WORK(KSCR1),NTOTAL)
C
C--------------------------------------------------
C        Calculate the frozen-correlated integrals.
C--------------------------------------------------
C
         KOFF11 = ILRHSI(ISYMAL) + 1
         KOFF12 = ICOFRO(ISYMJ,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('T','N',NRHF(ISYMJ),NRHFFR(ISYMI),NBAS(ISYMAL),
     *              ONE,WORK(KSCR1),NTOTAL,CMO(KOFF11),NTOTAL,ZERO,
     *              XFJI(KOFF12),NTOTJ)
C
C---------------------------------------------------
C        Transform second integral index to virtual.
C---------------------------------------------------
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
         ISYMA  = ISYMBE
C
         KOFF13 = IAODIS(ISYMAL,ISYMBE) + 1
         KOFF14 = ILMVIR(ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMA),NBAS(ISYMBE),
     *              ONE,XAOINT(KOFF13),NTOTAL,XLAMDH(KOFF14),NTOTBE,
     *              ZERO,WORK(KSCR1),NTOTAL)
C
C-----------------------------------------------
C        Calculate the frozen-virtual integrals.
C-----------------------------------------------
C
         KOFF13 = ILRHSI(ISYMAL) + 1
         KOFF14 = IT1FRO(ISYMA,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHFFR(ISYMI),NBAS(ISYMAL),
     *              ONE,WORK(KSCR1),NTOTAL,CMO(KOFF13),NTOTAL,ZERO,
     *              XFIA(KOFF14),NTOTA)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccfretai */
      SUBROUTINE CCFRETAI(ETAAI,XFIJ,XFJI,XFAI,XFIA,DFIJ,DFJI,
     *                    DFAI,DFIA,T1AM,WORK,LWORK,ISYM)
C
C     Written by Asger Halkier 13/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to eta-ai-0 from
C              intermediate looping over frozen.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAI(*), XFIJ(*), XFJI(*), XFAI(*), XFIA(*), DFIJ(*)
      DIMENSION DFJI(*), DFAI(*), DFIA(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C------------------------------
C     Do the four direct terms.
C------------------------------
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMJ = MULD2H(ISYMA,ISYM)
C
         IF (NRHFFR(ISYMJ) .EQ. 0) GOTO 100
C
         KOFFR = IT1AM(ISYMA,ISYMI)  + 1
         KOFF1 = IT1FRO(ISYMA,ISYMJ) + 1
         KOFF2 = ICOFRO(ISYMI,ISYMJ) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTI = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHFFR(ISYMJ),
     *              ONE,XFIA(KOFF1),NTOTA,DFJI(KOFF2),NTOTI,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHFFR(ISYMJ),
     *              -ONE,DFAI(KOFF1),NTOTA,XFIJ(KOFF2),NTOTI,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHFFR(ISYMJ),
     *              -ONE,DFIA(KOFF1),NTOTA,XFJI(KOFF2),NTOTI,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHFFR(ISYMJ),
     *              ONE,XFAI(KOFF1),NTOTA,DFIJ(KOFF2),NTOTI,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
  100 CONTINUE
C
C-----------------------------------------------------------
C     Do the terms involving one T1AM and loop over virtual.
C-----------------------------------------------------------
C
      DO 110 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMC = ISYMI
         ISYMJ = MULD2H(ISYMA,ISYM)
C
         LSCR  = NVIR(ISYMA)*NVIR(ISYMC)
C
         KSCR1 = 1
         KEND1 = KSCR1 + LSCR
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory in CCFRETAI')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
         KOFF3 = IT1FRO(ISYMA,ISYMJ) + 1
         KOFF4 = IT1FRO(ISYMC,ISYMJ) + 1
         KOFFT = IT1AM(ISYMC,ISYMI)  + 1
         KOFFR = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMC),NRHFFR(ISYMJ),
     *              ONE,DFIA(KOFF3),NTOTA,XFIA(KOFF4),NTOTC,ZERO,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMC),NRHFFR(ISYMJ),
     *              -ONE,XFAI(KOFF3),NTOTA,DFAI(KOFF4),NTOTC,ONE,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),
     *              ONE,WORK(KSCR1),NTOTA,T1AM(KOFFT),NTOTC,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
  110 CONTINUE
C
C--------------------------------------------------------------
C     Do the terms involving one T1AM and loop over correlated.
C--------------------------------------------------------------
C
      DO 120 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMK = ISYMA
         ISYMJ = MULD2H(ISYMI,ISYM)
C
         LSCR  = NRHF(ISYMK)*NRHF(ISYMI)
C
         KSCR2 = 1
         KEND2 = KSCR2 + LSCR
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory in CCFRETAI')
         ENDIF
C
         CALL DZERO(WORK(KSCR2),LSCR)
C
         KOFF5 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFF6 = ICOFRO(ISYMI,ISYMJ) + 1
         KOFFT = IT1AM(ISYMA,ISYMK)  + 1
         KOFFR = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTA = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','T',NRHF(ISYMK),NRHF(ISYMI),NRHFFR(ISYMJ),
     *              ONE,DFJI(KOFF5),NTOTK,XFJI(KOFF6),NTOTI,ZERO,
     *              WORK(KSCR2),NTOTK)
C
         CALL DGEMM('N','T',NRHF(ISYMK),NRHF(ISYMI),NRHFFR(ISYMJ),
     *              -ONE,XFIJ(KOFF5),NTOTK,DFIJ(KOFF6),NTOTI,ONE,
     *              WORK(KSCR2),NTOTK)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *              -ONE,T1AM(KOFFT),NTOTA,WORK(KSCR2),NTOTK,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
  120 CONTINUE
C
C--------------------------------------
C     Do the terms involving two T1AMs.
C--------------------------------------
C
      DO 130 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMK = ISYMA
         ISYMC = ISYMI
         ISYMJ = MULD2H(ISYMK,ISYM)
C
         LSCR3 = NRHF(ISYMK)*NVIR(ISYMC)
         LSCR4 = NRHF(ISYMK)*NRHF(ISYMI)
C
         KSCR3 = 1
         KSCR4 = KSCR3 + LSCR3
         KEND3 = KSCR4 + LSCR4
         LWRK3 = LWORK - KEND3
C
         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
            CALL QUIT('Insufficient memory in CCFRETAI')
         ENDIF
C
         CALL DZERO(WORK(KSCR3),KEND3)
C
         KOFF7 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFF8 = IT1FRO(ISYMC,ISYMJ) + 1
         KOFT1 = IT1AM(ISYMC,ISYMI)  + 1
         KOFT2 = IT1AM(ISYMA,ISYMK)  + 1
         KOFFR = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTA = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','T',NRHF(ISYMK),NVIR(ISYMC),NRHFFR(ISYMJ),
     *              ONE,DFJI(KOFF7),NTOTK,XFIA(KOFF8),NTOTC,ZERO,
     *              WORK(KSCR3),NTOTK)
C
         CALL DGEMM('N','T',NRHF(ISYMK),NVIR(ISYMC),NRHFFR(ISYMJ),
     *              -ONE,XFIJ(KOFF7),NTOTK,DFAI(KOFF8),NTOTC,ONE,
     *              WORK(KSCR3),NTOTK)
C
         CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     *              ONE,WORK(KSCR3),NTOTK,T1AM(KOFT1),NTOTC,ZERO,
     *              WORK(KSCR4),NTOTK) 
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *              ONE,T1AM(KOFT2),NTOTA,WORK(KSCR4),NTOTK,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
  130 CONTINUE
C
      RETURN
      END
C  /* Deck ccfretab */
      SUBROUTINE CCFRETAB(ETAAB,XFIJ,XFJI,XFAI,XFIA,DFIJ,DFJI,
     *                    DFAI,DFIA,T1AM,WORK,LWORK,ISYM)
C
C     Written by Asger Halkier 14/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to eta-ab-0 from
C              intermediate looping over frozen.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAB(*), XFIJ(*), XFJI(*), XFAI(*), XFIA(*), DFIJ(*)
      DIMENSION DFJI(*), DFAI(*), DFIA(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C------------------------------
C     Do the four direct terms.
C------------------------------
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMB = ISYMA
         ISYMJ = MULD2H(ISYMA,ISYM)
C
         IF (NRHFFR(ISYMJ) .EQ. 0) GOTO 100
C
         KOFF1 = IT1FRO(ISYMA,ISYMJ) + 1
         KOFF2 = IT1FRO(ISYMB,ISYMJ) + 1
         KOFFR = NMATIJ(1) + IMATAB(ISYMA,ISYMB) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHFFR(ISYMJ),
     *              ONE,XFIA(KOFF1),NTOTA,DFIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFR),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHFFR(ISYMJ),
     *              -ONE,DFAI(KOFF1),NTOTA,XFAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFR),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHFFR(ISYMJ),
     *              -ONE,DFIA(KOFF1),NTOTA,XFIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFR),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHFFR(ISYMJ),
     *              ONE,XFAI(KOFF1),NTOTA,DFAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFR),NTOTA)
C
  100 CONTINUE
C
C--------------------------------------
C     Do the two first terms with T1AM.
C--------------------------------------
C
      DO 110 ISYMA = 1,NSYM
C
         ISYMB = ISYMA
         ISYMK = ISYMB
         ISYMJ = MULD2H(ISYMA,ISYM)
C
         LSCR  = NVIR(ISYMA)*NRHF(ISYMK)
C
         KSCR1 = 1
         KEND1 = KSCR1 + LSCR
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory in CCFRETAB')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
         KOFF3 = IT1FRO(ISYMA,ISYMJ) + 1
         KOFF4 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFFT = IT1AM(ISYMB,ISYMK)  + 1
         KOFFR = NMATIJ(1) + IMATAB(ISYMA,ISYMB) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTB = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMK),NRHFFR(ISYMJ),
     *              ONE,XFIA(KOFF3),NTOTA,DFJI(KOFF4),NTOTK,ZERO,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMK),NRHFFR(ISYMJ),
     *              -ONE,DFAI(KOFF3),NTOTA,XFIJ(KOFF4),NTOTK,ONE,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),
     *              ONE,WORK(KSCR1),NTOTA,T1AM(KOFFT),NTOTB,ONE,
     *              ETAAB(KOFFR),NTOTA)
C
  110 CONTINUE
C
C-------------------------------------
C     Do the last two terms with T1AM.
C-------------------------------------
C
      DO 120 ISYMA = 1,NSYM
C
         ISYMB = ISYMA
         ISYMK = ISYMA
         ISYMJ = MULD2H(ISYMA,ISYM)
C
         LSCR  = NVIR(ISYMB)*NRHF(ISYMK)
C
         KSCR2 = 1
         KEND2 = KSCR2 + LSCR
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory in CCFRETAB')
         ENDIF
C
         CALL DZERO(WORK(KSCR2),LSCR)
C
         KOFF5 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFF6 = IT1FRO(ISYMB,ISYMJ) + 1
         KOFFT = IT1AM(ISYMA,ISYMK)  + 1
         KOFFR = NMATIJ(1) + IMATAB(ISYMA,ISYMB) + 1
C
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTB = MAX(NVIR(ISYMB),1)
         NTOTA = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','T',NRHF(ISYMK),NVIR(ISYMB),NRHFFR(ISYMJ),
     *              ONE,DFJI(KOFF5),NTOTK,XFIA(KOFF6),NTOTB,ZERO,
     *              WORK(KSCR2),NTOTK)
C
         CALL DGEMM('N','T',NRHF(ISYMK),NVIR(ISYMB),NRHFFR(ISYMJ),
     *              -ONE,XFIJ(KOFF5),NTOTK,DFAI(KOFF6),NTOTB,ONE,
     *              WORK(KSCR2),NTOTK)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),
     *              -ONE,T1AM(KOFFT),NTOTA,WORK(KSCR2),NTOTK,ONE,
     *              ETAAB(KOFFR),NTOTA)
C
  120 CONTINUE
C
      RETURN
      END
C  /* Deck ccfretij */
      SUBROUTINE CCFRETIJ(ETAIJ,XFIJ,XFJI,XFAI,XFIA,DFIJ,DFJI,
     *                    DFAI,DFIA,T1AM,WORK,LWORK,ISYM)
C
C     Written by Asger Halkier 14/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to eta-ij-0 from
C              intermediate looping over frozen.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), XFIJ(*), XFJI(*), XFAI(*), XFIA(*), DFIJ(*)
      DIMENSION DFJI(*), DFAI(*), DFIA(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C------------------------------
C     Do the four direct terms.
C------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMK = MULD2H(ISYMI,ISYM)
C
         IF (NRHFFR(ISYMK) .EQ. 0) GOTO 100
C
         KOFF1 = ICOFRO(ISYMI,ISYMK) + 1
         KOFF2 = ICOFRO(ISYMJ,ISYMK) + 1
         KOFFR = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFJI(KOFF1),NTOTI,DFJI(KOFF2),NTOTJ,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHFFR(ISYMK),
     *              -ONE,DFIJ(KOFF1),NTOTI,XFIJ(KOFF2),NTOTJ,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHFFR(ISYMK),
     *              -ONE,DFJI(KOFF1),NTOTI,XFJI(KOFF2),NTOTJ,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFIJ(KOFF1),NTOTI,DFIJ(KOFF2),NTOTJ,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
  100 CONTINUE
C
C--------------------------------------
C     Do the two first terms with T1AM.
C--------------------------------------
C
      DO 110 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMC = ISYMJ
         ISYMK = MULD2H(ISYMI,ISYM)
C
         LSCR  = NRHF(ISYMI)*NVIR(ISYMC)
C
         KSCR1 = 1
         KEND1 = KSCR1 + LSCR
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory in CCFRETIJ')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
         KOFF3 = ICOFRO(ISYMI,ISYMK) + 1
         KOFF4 = IT1FRO(ISYMC,ISYMK) + 1
         KOFFT = IT1AM(ISYMC,ISYMJ)  + 1
         KOFFR = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTC = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NVIR(ISYMC),NRHFFR(ISYMK),
     *              ONE,DFJI(KOFF3),NTOTI,XFIA(KOFF4),NTOTC,ZERO,
     *              WORK(KSCR1),NTOTI)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NVIR(ISYMC),NRHFFR(ISYMK),
     *              -ONE,XFIJ(KOFF3),NTOTI,DFAI(KOFF4),NTOTC,ONE,
     *              WORK(KSCR1),NTOTI)
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),
     *              ONE,WORK(KSCR1),NTOTI,T1AM(KOFFT),NTOTC,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
  110 CONTINUE
C
C-------------------------------------
C     Do the last two terms with T1AM.
C-------------------------------------
C
      DO 120 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMC = ISYMI
         ISYMK = MULD2H(ISYMJ,ISYM)
C
         LSCR  = NRHF(ISYMJ)*NVIR(ISYMC)
C
         KSCR2 = 1
         KEND2 = KSCR2 + LSCR
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory in CCFRETIJ')
         ENDIF
C
         CALL DZERO(WORK(KSCR2),LSCR)
C
         KOFF5 = IT1FRO(ISYMC,ISYMK) + 1
         KOFF6 = ICOFRO(ISYMJ,ISYMK) + 1
         KOFFT = IT1AM(ISYMC,ISYMI)  + 1
         KOFFR = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)
         NTOTI = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('N','T',NVIR(ISYMC),NRHF(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFIA(KOFF5),NTOTC,DFJI(KOFF6),NTOTJ,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('N','T',NVIR(ISYMC),NRHF(ISYMJ),NRHFFR(ISYMK),
     *              -ONE,DFAI(KOFF5),NTOTC,XFIJ(KOFF6),NTOTJ,ONE,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),
     *              -ONE,T1AM(KOFFT),NTOTC,WORK(KSCR2),NTOTC,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
  120 CONTINUE
C
      RETURN
      END
C  /* Deck ccfd1ii */
      SUBROUTINE CCFD1II(DFII)
C
C     Written by Asger Halkier 14/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the simple frozen-frozen 1 electron density.
C
#include <implicit.h>
      PARAMETER(TWO = 2.0D0)
      DIMENSION DFII(*)
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CALL DZERO(DFII,NFROFR(1))
C
      DO 100 ISYM = 1,NSYM
         DO 110 I = 1,NRHFFR(ISYM)
            NII = IFROFR(ISYM,ISYM) + NRHFFR(ISYM)*(I - 1) + I
            DFII(NII) = TWO
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_etijf */
      SUBROUTINE CC_ETIJF(ETAIJ,DIJ,DAB,DAI,DIA,DIJT,DFII,DFIJ,DFJI,
     *                    DFAI,DFIA,XIJ,XAI,XIA,XIJT,XABT,XFIJ,XFJI,
     *                    XFAI,XFIA,XFII,T1AM,WORK,LWORK,ISYM,IOPT)
C
C     Written by Asger Halkier 14/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate eta-iJ-0. ISYM is the symmetry of the
C              integrals and densities and IOPT controls the one
C              and two electron contributions.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), DIJ(*), DAB(*), DAI(*), DIA(*), DIJT(*)
      DIMENSION DFII(*), DFIJ(*), DFJI(*), DFAI(*), DFIA(*), XIJ(*)
      DIMENSION XAI(*), XIA(*), XIJT(*), XABT(*), XFIJ(*), XFJI(*)
      DIMENSION XFAI(*), XFIA(*), XFII(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C-------------------------------
C     Do the direct terms first.
C-------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMC = MULD2H(ISYMI,ISYM)
         ISYMK = ISYMC
C
         IF (NRHFFR(ISYMJ) .EQ. 0 ) GOTO 100
C
         KOFFR = ICOFRO(ISYMI,ISYMJ) + 1
         KOFF1 = IT1AM(ISYMC,ISYMI)  + 1
         KOFF2 = IT1FRO(ISYMC,ISYMJ) + 1
C
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTC = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *              -ONE,DIA(KOFF1),NTOTC,XFIA(KOFF2),NTOTC,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *              -ONE,DAI(KOFF1),NTOTC,XFAI(KOFF2),NTOTC,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         KOFF3 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF4 = ICOFRO(ISYMK,ISYMJ) + 1
C
         NTOTK = MAX(NRHF(ISYMK),1)
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *              -ONE,DIJ(KOFF3),NTOTI,XFJI(KOFF4),NTOTK,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *              -ONE,DIJT(KOFF3),NTOTI,XFIJ(KOFF4),NTOTK,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         KOFF5 = ICOFRO(ISYMI,ISYMK) + 1
         KOFF6 = IFROFR(ISYMK,ISYMJ) + 1
         KOFF7 = IFROFR(ISYMJ,ISYMK) + 1
C
         NTOTK = MAX(NRHFFR(ISYMK),1)
         NTOTJ = MAX(NRHFFR(ISYMJ),1)
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFJI(KOFF5),NTOTI,DFII(KOFF6),NTOTK,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFIJ(KOFF5),NTOTI,DFII(KOFF7),NTOTJ,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
C------------------------------------------------------------
C        The terms that only has contributions from 2e- part.
C------------------------------------------------------------
C
         IF (IOPT .EQ. 2) THEN
C
            CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *                 ONE,XAI(KOFF1),NTOTC,DFAI(KOFF2),NTOTC,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *                 ONE,XIA(KOFF1),NTOTC,DFIA(KOFF2),NTOTC,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            NTOTK = MAX(NRHF(ISYMK),1)
C
            CALL DGEMM('N','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *                 ONE,XIJT(KOFF3),NTOTI,DFIJ(KOFF4),NTOTK,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('N','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *                 ONE,XIJ(KOFF3),NTOTI,DFJI(KOFF4),NTOTK,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            NTOTK = MAX(NRHFFR(ISYMK),1)
C
            CALL DGEMM('N','T',NRHF(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *                 -ONE,DFIJ(KOFF5),NTOTI,XFII(KOFF7),NTOTJ,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('N','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *                 -ONE,DFJI(KOFF5),NTOTI,XFII(KOFF6),NTOTK,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
         ENDIF
  100 CONTINUE
C
C---------------------------------
C     Do the terms involving T1AM.
C---------------------------------
C
      DO 110 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMC = ISYMI
         ISYMD = MULD2H(ISYMJ,ISYM)
         ISYMK = MULD2H(ISYMJ,ISYM)
C
         IF (NRHFFR(ISYMJ) .EQ. 0 ) GOTO 110
C
         LSCR  = NRHFFR(ISYMJ)*NVIR(ISYMC)
C
         KSCR1 = 1
         KEND1 = KSCR1 + LSCR
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory in CC_ETIJF')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
         KOFFT  = IT1AM(ISYMC,ISYMI)  + 1
         KOFFR  = ICOFRO(ISYMI,ISYMJ) + 1
C
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTD  = MAX(NVIR(ISYMD),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
         NTOKF  = MAX(NRHFFR(ISYMK),1)
C
         KOFF8  = IMATAB(ISYMC,ISYMD) + 1
         KOFF9  = IT1FRO(ISYMD,ISYMJ) + 1
         KOFF10 = IT1AM(ISYMC,ISYMK)  + 1
         KOFF11 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFF12 = IT1FRO(ISYMC,ISYMK) + 1
         KOFF13 = IFROFR(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHFFR(ISYMJ),NVIR(ISYMD),
     *              -ONE,DAB(KOFF8),NTOTC,XFIA(KOFF9),NTOTD,ZERO,
     *              WORK(KSCR1),NTOTC)
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHFFR(ISYMJ),NRHF(ISYMK),
     *              -ONE,DAI(KOFF10),NTOTC,XFJI(KOFF11),NTOTK,ONE,
     *              WORK(KSCR1),NTOTC)
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFIA(KOFF12),NTOTC,DFII(KOFF13),NTOKF,ONE,
     *              WORK(KSCR1),NTOTC)
C
C------------------------------------------------------------
C        The terms that only has contributions from 2e- part.
C------------------------------------------------------------
C
         IF (IOPT .EQ. 2) THEN
C
            CALL DGEMM('N','N',NVIR(ISYMC),NRHFFR(ISYMJ),NVIR(ISYMD),
     *                 ONE,XABT(KOFF8),NTOTC,DFAI(KOFF9),NTOTD,ONE,
     *                 WORK(KSCR1),NTOTC)
C
            CALL DGEMM('N','N',NVIR(ISYMC),NRHFFR(ISYMJ),NRHF(ISYMK),
     *                 ONE,XIA(KOFF10),NTOTC,DFIJ(KOFF11),NTOTK,ONE,
     *                 WORK(KSCR1),NTOTC)
C
            KOFF14 = IFROFR(ISYMJ,ISYMK) + 1
            NTOTJ  = MAX(NRHFFR(ISYMJ),1)
C
            CALL DGEMM('N','T',NVIR(ISYMC),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *                 -ONE,DFAI(KOFF12),NTOTC,XFII(KOFF14),NTOTJ,ONE,
     *                 WORK(KSCR1),NTOTC)
C
         ENDIF
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *              -ONE,T1AM(KOFFT),NTOTC,WORK(KSCR1),NTOTC,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
  110 CONTINUE
C
      RETURN
      END
C  /* Deck cc_etifjf */
      SUBROUTINE CC_ETIFJF(ETAIJ,DFII,DFIJ,DFJI,DFAI,DFIA,XFIJ,XFJI,
     *                     XFAI,XFIA,XFII,ISYM,IOPT)
C
C     Written by Asger Halkier 15/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate eta-IJ-0. ISYM is the symmetry of the
C              integrals and densities and IOPT controls the one
C              and two electron contributions.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), DFII(*), DFIJ(*), DFJI(*), DFAI(*), DFIA(*)
      DIMENSION XFIJ(*), XFJI(*), XFAI(*), XFIA(*), XFII(*)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMK = MULD2H(ISYMI,ISYM)
         ISYMC = MULD2H(ISYMI,ISYM)
C
         IF ((NRHFFR(ISYMI) .EQ. 0).OR.(NRHFFR(ISYMJ) .EQ. 0)) GOTO 100
C
         KOFFR = IFROFR(ISYMI,ISYMJ) + 1
         KOFF1 = IFROFR(ISYMK,ISYMI) + 1
         KOFF2 = IFROFR(ISYMK,ISYMJ) + 1
         KOFF3 = IFROFR(ISYMI,ISYMK) + 1
         KOFF4 = IFROFR(ISYMJ,ISYMK) + 1
C
         NTOTI = MAX(NRHFFR(ISYMI),1)
         NTOTK = MAX(NRHFFR(ISYMK),1)
         NTOTJ = MAX(NRHFFR(ISYMJ),1)
C
         CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFII(KOFF1),NTOTK,DFII(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('N','T',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *              -ONE,DFII(KOFF3),NTOTI,XFII(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *              -ONE,DFII(KOFF1),NTOTK,XFII(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
         CALL DGEMM('N','T',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHFFR(ISYMK),
     *              ONE,XFII(KOFF3),NTOTI,DFII(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFR),NTOTI)
C
C------------------------------------------------------------
C        The terms that only has contributions from 2e- part.
C------------------------------------------------------------
C
         IF (IOPT .EQ. 2) THEN
C
            KOFF5 = ICOFRO(ISYMK,ISYMI) + 1
            KOFF6 = ICOFRO(ISYMK,ISYMJ) + 1
C
            NTOTK = MAX(NRHF(ISYMK),1)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *                 ONE,XFIJ(KOFF5),NTOTK,DFIJ(KOFF6),NTOTK,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *                 -ONE,DFJI(KOFF5),NTOTK,XFJI(KOFF6),NTOTK,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *                 -ONE,DFIJ(KOFF5),NTOTK,XFIJ(KOFF6),NTOTK,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *                 ONE,XFJI(KOFF5),NTOTK,DFJI(KOFF6),NTOTK,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            KOFF7 = IT1FRO(ISYMC,ISYMI) + 1
            KOFF8 = IT1FRO(ISYMC,ISYMJ) + 1
C
            NTOTC = MAX(NVIR(ISYMC),1)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *                 ONE,XFAI(KOFF7),NTOTC,DFAI(KOFF8),NTOTC,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *                 -ONE,DFIA(KOFF7),NTOTC,XFIA(KOFF8),NTOTC,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *                 -ONE,DFAI(KOFF7),NTOTC,XFAI(KOFF8),NTOTC,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *                 ONE,XFIA(KOFF7),NTOTC,DFIA(KOFF8),NTOTC,ONE,
     *                 ETAIJ(KOFFR),NTOTI)
C
         ENDIF
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_zkfcb */
      SUBROUTINE CC_ZKFCB(ZKDIA,WORK,LWORK)
C
C     Written by Asger Halkier 15/8 - 1998.
C
C     To calculate the frozen-occupied block (iJ) of kappa-bar-0.
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
      DIMENSION ZKDIA(*), WORK(LWORK)
C
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KETAIJ = 1
      KFOCKD = KETAIJ + NCOFRO(1)
      KEND1  = KFOCKD + NORBTS
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in CC_ZKFCB')
      ENDIF
C
      CALL DZERO(WORK(KETAIJ),NCOFRO(1))
      CALL DZERO(WORK(KFOCKD),NORBTS)
      CALL DCOPY(NCOFRO(1),ZKDIA(1),1,WORK(KETAIJ),1)
      CALL DZERO(ZKDIA(1),NCOFRO(1))
C
C--------------------------
C     Get orbital energies.
C--------------------------
C
      CALL FOCK_ALL(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C---------------------------
C     Calculate the results.
C---------------------------
C
      DO 130 ISYMI = 1,NSYM
         ISYMJ = ISYMI
         DO 140 J = 1,NRHFFR(ISYMJ)
            DO 150 I = 1,NRHF(ISYMI)
               NIJ = ICOFRO(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
               KOI = KFOCKD + IRHFS(ISYMI) + NRHFFR(ISYMI) + I - 1
               KOJ = KFOCKD + IRHFS(ISYMJ) + J - 1
C
               DENOM    = WORK(KOJ) - WORK(KOI)
               ZKDIA(NIJ) = HALF*WORK(KETAIJ+NIJ-1)/DENOM
C
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      CALL DCOPY(NCOFRO(1),ZKDIA(1),1,ZKDIA(1+NCOFRO(1)),1)
C
      RETURN
      END
C  /* Deck cc_etaif */
      SUBROUTINE CC_ETAIF(ETAAI,DAB,DAI,DIA,DIJT,DABT,DFII,DFIJ,DFJI,
     *                    DFAI,DFIA,XIJ,XAB,XAI,XIA,XABT,XFIJ,XFJI,
     *                    XFAI,XFIA,XFII,T1AM,WORK,LWORK,ISYM,IOPT)
C
C     Written by Asger Halkier 16/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate eta-aI-0. ISYM is the symmetry of the
C              integrals and densities and IOPT controls the one
C              and two electron contributions.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAI(*), DAB(*), DAI(*), DIA(*), DIJT(*), DABT(*)
      DIMENSION DFII(*), DFIJ(*), DFJI(*), DFAI(*), DFIA(*), XIJ(*)
      DIMENSION XAB(*), XAI(*), XIA(*), XABT(*), XFIJ(*), XFJI(*)
      DIMENSION XFAI(*), XFIA(*), XFII(*), T1AM(*), WORK(LWORK)
#include <priunit.h>      
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C-------------------------------
C     Do the direct terms first.
C-------------------------------
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMC = MULD2H(ISYMA,ISYM)
         ISYMK = MULD2H(ISYMA,ISYM)
         ISYMJ = MULD2H(ISYMA,ISYM)
C
         IF (NRHFFR(ISYMI) .EQ. 0) GOTO 100
C
         KOFFR = IT1FRO(ISYMA,ISYMI) + 1
         KOFF1 = IMATAB(ISYMA,ISYMC) + 1
         KOFF2 = IT1FRO(ISYMC,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF4 = ICOFRO(ISYMK,ISYMI) + 1
         KOFF5 = IT1FRO(ISYMA,ISYMJ) + 1
         KOFF6 = IFROFR(ISYMJ,ISYMI) + 1
         KOFF7 = IFROFR(ISYMI,ISYMJ) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTJ = MAX(NRHFFR(ISYMJ),1)
         NTOTI = MAX(NRHFFR(ISYMI),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NVIR(ISYMC),
     *              -ONE,DAB(KOFF1),NTOTA,XFIA(KOFF2),NTOTC,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NVIR(ISYMC),
     *              -ONE,DABT(KOFF1),NTOTA,XFAI(KOFF2),NTOTC,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *              -ONE,DAI(KOFF3),NTOTA,XFJI(KOFF4),NTOTK,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *              -ONE,DIA(KOFF3),NTOTA,XFIJ(KOFF4),NTOTK,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHFFR(ISYMJ),
     *              ONE,XFIA(KOFF5),NTOTA,DFII(KOFF6),NTOTJ,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHFFR(ISYMI),NRHFFR(ISYMJ),
     *              ONE,XFAI(KOFF5),NTOTA,DFII(KOFF7),NTOTI,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
C------------------------------------------------------------
C        The terms that only has contributions from 2e- part.
C------------------------------------------------------------
C
         IF (IOPT .EQ. 2) THEN
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NVIR(ISYMC),
     *                 ONE,XABT(KOFF1),NTOTA,DFAI(KOFF2),NTOTC,ONE,
     *                 ETAAI(KOFFR),NTOTA)
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NVIR(ISYMC),
     *                 ONE,XAB(KOFF1),NTOTA,DFIA(KOFF2),NTOTC,ONE,
     *                 ETAAI(KOFFR),NTOTA)
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *                 ONE,XIA(KOFF3),NTOTA,DFIJ(KOFF4),NTOTK,ONE,
     *                 ETAAI(KOFFR),NTOTA)
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *                 ONE,XAI(KOFF3),NTOTA,DFJI(KOFF4),NTOTK,ONE,
     *                 ETAAI(KOFFR),NTOTA)
C
            CALL DGEMM('N','T',NVIR(ISYMA),NRHFFR(ISYMI),NRHFFR(ISYMJ),
     *                 -ONE,DFAI(KOFF5),NTOTA,XFII(KOFF7),NTOTI,ONE,
     *                 ETAAI(KOFFR),NTOTA)
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHFFR(ISYMJ),
     *                 -ONE,DFIA(KOFF5),NTOTA,XFII(KOFF6),NTOTJ,ONE,
     *                 ETAAI(KOFFR),NTOTA)
C
         ENDIF
  100 CONTINUE
C
C---------------------------------
C     Do the terms involving T1AM.
C---------------------------------
C
      DO 110 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMK = ISYMA
         ISYMC = MULD2H(ISYMI,ISYM)
         ISYML = MULD2H(ISYMI,ISYM)
         ISYMJ = MULD2H(ISYMI,ISYM)
C
         IF (NRHFFR(ISYMI) .EQ. 0) GOTO 110
C
         LSCR = NRHF(ISYMK)*NRHFFR(ISYMI)
C
         KSCR1 = 1
         KEND1 = KSCR1 + LSCR
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory in CC_ETAIF')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),LSCR)
C
         KOFFR  = IT1FRO(ISYMA,ISYMI) + 1
         KOFFT  = IT1AM(ISYMA,ISYMK)  + 1
         KOFF8  = IT1AM(ISYMC,ISYMK)  + 1
         KOFF9  = IT1FRO(ISYMC,ISYMI) + 1
         KOFF10 = IMATIJ(ISYMK,ISYML) + 1
         KOFF11 = ICOFRO(ISYML,ISYMI) + 1
         KOFF12 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFF13 = IFROFR(ISYMI,ISYMJ) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
         NTOTL  = MAX(NRHF(ISYML),1)
         NTOTI  = MAX(NRHFFR(ISYMI),1)
C
         CALL DGEMM('T','N',NRHF(ISYMK),NRHFFR(ISYMI),NVIR(ISYMC),
     *              ONE,DAI(KOFF8),NTOTC,XFAI(KOFF9),NTOTC,ZERO,
     *              WORK(KSCR1),NTOTK)
C
         CALL DGEMM('N','N',NRHF(ISYMK),NRHFFR(ISYMI),NRHF(ISYML),
     *              ONE,DIJT(KOFF10),NTOTK,XFIJ(KOFF11),NTOTL,ONE,
     *              WORK(KSCR1),NTOTK)
C
         CALL DGEMM('N','T',NRHF(ISYMK),NRHFFR(ISYMI),NRHFFR(ISYMJ),
     *              -ONE,XFIJ(KOFF12),NTOTK,DFII(KOFF13),NTOTI,ONE,
     *              WORK(KSCR1),NTOTK)
C
C------------------------------------------------------------
C        The terms that only has contributions from 2e- part.
C------------------------------------------------------------
C
         IF (IOPT .EQ. 2) THEN
C
            CALL DGEMM('T','N',NRHF(ISYMK),NRHFFR(ISYMI),NVIR(ISYMC),
     *                 -ONE,XIA(KOFF8),NTOTC,DFIA(KOFF9),NTOTC,ONE,
     *                 WORK(KSCR1),NTOTK)
C
            CALL DGEMM('N','N',NRHF(ISYMK),NRHFFR(ISYMI),NRHF(ISYML),
     *                 -ONE,XIJ(KOFF10),NTOTK,DFJI(KOFF11),NTOTL,ONE,
     *                 WORK(KSCR1),NTOTK)
C
            KOFF14 = IFROFR(ISYMJ,ISYMI) + 1
            NTOTJ  = MAX(NRHFFR(ISYMJ),1)
C
            CALL DGEMM('N','N',NRHF(ISYMK),NRHFFR(ISYMI),NRHFFR(ISYMJ),
     *                 ONE,DFJI(KOFF12),NTOTK,XFII(KOFF14),NTOTJ,ONE,
     *                 WORK(KSCR1),NTOTK)
C
         ENDIF
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *              -ONE,T1AM(KOFFT),NTOTA,WORK(KSCR1),NTOTK,ONE,
     *              ETAAI(KOFFR),NTOTA)
C
  110 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcoc1 */
      SUBROUTINE CC_FRCOC1(XINTAI,XINOFO,ISYDIS)
C
C     Written by Asger Halkier 17/8 - 1998.
C
C     To calculate the coulomb part and the first exchange part of the
C     intermediate needed for the frozen-correlated correction to eta-ai.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0, EIGHT = 8.0D0)
      DIMENSION XINTAI(*), XINOFO(*)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDIS
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHF(ISYMI) .EQ. 0) RETURN
C
C-------------------------------
C     Copy the coulomb part out.
C-------------------------------
C
      KOFF1 = IOFROO(ISYMJK,ISYMI) + 1
      KOFF2 = IOFROO(ISYMJK,ISYMI) + 1
      LEN   = NCOFRO(ISYMJK)*NRHF(ISYMI)
      CALL DAXPY(LEN,EIGHT,XINOFO(KOFF1),1,XINTAI(KOFF2),1)
C
C--------------------------------
C     Do the first exchange part.
C--------------------------------
C
      DO 100 ISYMJ = 1,NSYM
         ISYMK  = MULD2H(ISYMJ,ISYMJK)
         ISYMIK = MULD2H(ISYMI,ISYMK)
         DO 110 J = 1,NRHF(ISYMJ)
            DO 120 K = 1,NRHFFR(ISYMK)
C
               KOFF2 = IOFROO(ISYMIK,ISYMJ) + NCOFRO(ISYMIK)*(J - 1)
     *               + ICOFRO(ISYMI,ISYMK)  + NRHF(ISYMI)*(K - 1) + 1
               KOFF3 = IOFROO(ISYMJK,ISYMI) + ICOFRO(ISYMJ,ISYMK)
     *               + NRHF(ISYMJ)*(K - 1)  + J
C
               CALL DAXPY(NRHF(ISYMI),-TWO,XINOFO(KOFF2),1,
     *                    XINTAI(KOFF3),NCOFRO(ISYMJK))
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcoc1 */
      SUBROUTINE CC_FRINDU(XINTAI,XINAIF,IDEL,ISYMD,LUN1,LUN2)
C
C     Written by Asger Halkier 17/8 - 1998.
C
C     To dump intermediates needed for the frozen-correlated corrections
C     to eta-ai & eta-aI to disc.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0, EIGHT = 8.0D0)
      DIMENSION XINTAI(*), XINAIF(*)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CHARACTER NAME1*8
      CHARACTER NAME2*8
C
      NAME1 = 'CCFRO1IN'
      NAME2 = 'CCFRO2IN'
C
C------------------------------------------
C     Write integral intermediates to disc.
C------------------------------------------
C
      D     = IDEL - IBAS(ISYMD)
      NTOT  = NOFROO(ISYMD)
      IOFF  = IOFOAO(ISYMD,ISYMD) + NTOT*(D - 1) + 1
C
      IF (NTOT .GT. 0) THEN
         CALL PUTWA2(LUN1,NAME1,XINTAI,IOFF,NTOT)
      ENDIF
C
      D     = IDEL - IBAS(ISYMD)
      NTOT  = NCOFRF(ISYMD)
      IOFF  = IOFFAO(ISYMD,ISYMD) + NTOT*(D - 1) + 1
C
      IF (NTOT .GT. 0) THEN
         CALL PUTWA2(LUN2,NAME2,XINAIF,IOFF,NTOT)
      ENDIF
C
      RETURN
      END
C  /* Deck cc_etacor */
      SUBROUTINE CC_ETACOR(ETAAI,ETAAIF,ZETAJK,CMO,LUN1,LUN2,WORK,LWORK)
C
C     Written by Asger Halkier 17/8 - 1998.
C
C     To calculate the frozen-correlated corrections to eta-ai &
C     eta-aI from intermediates on disc.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0, ONE = 1.0D0)
      DIMENSION ETAAI(*), ETAAIF(*), ZETAJK(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CHARACTER NAME1*8
      CHARACTER NAME2*8
C
      NAME1 = 'CCFRO1IN'
      NAME2 = 'CCFRO2IN'
C
C-----------------------------
C     First we do the ai part.
C-----------------------------
C
      DO 100 ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 100
C
         ISYJKI = ISYMD
         ISYMA  = ISYMD
         ISYMI  = ISYMA
         ISYMJK = 1
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCR1 = 1
         KINIM = KSCR1 + NBAS(ISYMD)*NOFROO(ISYJKI)
         KEND1 = KINIM + NVIR(ISYMA)*NOFROO(ISYJKI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory in CC_ETACOR')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND1)
C
C---------------------------------
C        Read integrals from disc.
C---------------------------------
C
         NTOT = NOFROO(ISYJKI)*NBAS(ISYMD)
         IOFF = IOFOAO(ISYJKI,ISYMD) + 1
C
         IF (NTOT .GT. 0) THEN
           CALL GETWA2(LUN1,NAME1,WORK(KSCR1),IOFF,NTOT)
         END IF
C
C--------------------------------------
C        Transform AO index to virtual.
C--------------------------------------
C
         KOFF1  = ILVISI(ISYMD) + 1
C
         NTOJKI = MAX(NOFROO(ISYJKI),1)
         NTOTD  = MAX(NBAS(ISYMD),1)
C
         CALL DGEMM('N','N',NOFROO(ISYJKI),NVIR(ISYMA),NBAS(ISYMD),
     *              ONE,WORK(KSCR1),NTOJKI,CMO(KOFF1),NTOTD,ZERO,
     *              WORK(KINIM),NTOJKI)
C
C-----------------------------------
C        Calculate the contribution.
C-----------------------------------
C
         DO 110 A = 1,NVIR(ISYMA)
            DO 120 I = 1,NRHF(ISYMI)
               KOFF2 = KINIM + NOFROO(ISYJKI)*(A - 1)
     *               + IOFROO(ISYMJK,ISYMI) + NCOFRO(ISYMJK)*(I - 1)
               KOFF3 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
               FACT  = DDOT(NCOFRO(1),WORK(KOFF2),1,ZETAJK(1),1)
               ETAAI(KOFF3) = ETAAI(KOFF3) + FACT
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C----------------------------
C     Then we do the aI part.
C----------------------------
C
      DO 130 ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 130
C
         ISYJKI = ISYMD
         ISYMA  = ISYMD
         ISYMI  = ISYMA
         ISYMJK = 1
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KSCR2 = 1
         KINIM = KSCR2 + NBAS(ISYMD)*NCOFRF(ISYJKI)
         KEND2 = KINIM + NVIR(ISYMA)*NCOFRF(ISYJKI)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory in CC_ETACOR')
         ENDIF
C
         CALL DZERO(WORK(KSCR2),KEND2)
C
C---------------------------------
C        Read integrals from disc.
C---------------------------------
C
         NTOT = NCOFRF(ISYJKI)*NBAS(ISYMD)
         IOFF = IOFFAO(ISYJKI,ISYMD) + 1
C
         IF (NTOT .GT. 0) THEN
           CALL GETWA2(LUN2,NAME2,WORK(KSCR2),IOFF,NTOT)
         END IF
C
C--------------------------------------
C        Transform AO index to virtual.
C--------------------------------------
C
         KOFF4  = ILVISI(ISYMD) + 1
C
         NTOJKI = MAX(NCOFRF(ISYJKI),1)
         NTOTD  = MAX(NBAS(ISYMD),1)
C
         CALL DGEMM('N','N',NCOFRF(ISYJKI),NVIR(ISYMA),NBAS(ISYMD),
     *              ONE,WORK(KSCR2),NTOJKI,CMO(KOFF4),NTOTD,ZERO,
     *              WORK(KINIM),NTOJKI)
C
C-----------------------------------
C        Calculate the contribution.
C-----------------------------------
C
         DO 140 A = 1,NVIR(ISYMA)
            DO 150 I = 1,NRHFFR(ISYMI)
               KOFF5 = KINIM + NCOFRF(ISYJKI)*(A - 1)
     *               + ICOFRF(ISYMJK,ISYMI) + NCOFRO(ISYMJK)*(I - 1)
               KOFF6 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
               FACT  = DDOT(NCOFRO(1),WORK(KOFF5),1,ZETAJK(1),1)
               ETAAIF(KOFF6) = ETAAIF(KOFF6) + FACT
C
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcomi */
      SUBROUTINE CC_FRCOMI(XINTAI,XINAIF,DSFRO,CMO,WORK,LWORK,
     *                     ISYDEL)
C
C     Written by Asger Halkier 17/8 - 1998.
C
C     To calculate the second exchange part of the intermediates needed
C     for the frozen-correlated correction to eta-ai and eta-aI.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINTAI(*), XINAIF(*), DSFRO(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
C---------------------------------------------------------------------------
C     Loops over third index in integral batch & final wok space allocation.
C---------------------------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMJK,ISYMK)
         ISYMIJ = MULD2H(ISYMI,ISYMJ)
         ISALBE = ISYMIJ
         ISYMAL = ISYMI
         ISYMBE = ISYMJ
C
         KAOINT = 1
         KSCR1  = KAOINT + N2BST(ISALBE)
         KMOIN1 = KSCR1  + NBAS(ISYMAL)*NRHF(ISYMJ)
         KMOIN2 = KMOIN1 + NRHF(ISYMI)*NRHF(ISYMJ)
         KEND1  = KMOIN2 + NRHFFR(ISYMI)*NRHF(ISYMJ)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space for allocation '//
     &                'in CC_FRCOMI')
         ENDIF
C
         CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
         CALL DZERO(WORK(KMOIN1),NRHF(ISYMI)*NRHF(ISYMJ))
         CALL DZERO(WORK(KMOIN2),NRHFFR(ISYMI)*NRHF(ISYMJ))
C
         DO 110 K = 1,NRHFFR(ISYMK)
C
C-----------------------------------------------------------
C           Square up integral distribution (al be | K del).
C-----------------------------------------------------------
C
            KOFF1 = IDSFRO(ISALBE,ISYMK) + NNBST(ISALBE)*(K - 1) + 1
C
            CALL CCSD_SYMSQ(DSFRO(KOFF1),ISALBE,WORK(KAOINT))
C
C---------------------------------------------
C           Calculate integrals (i j | K del).
C---------------------------------------------
C
            KOFF2  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF3  = ILRHSI(ISYMBE) + NBAS(ISYMBE)*NRHFFR(ISYMJ) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF2),NTOTAL,CMO(KOFF3),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF4  = ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMI) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTI  = MAX(NRHF(ISYMI),1)
C
            CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF4),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KMOIN1),NTOTI)
C
C-----------------------------------------------
C           Add contribution to ai-intermediate.
C-----------------------------------------------
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               KOFF5 = KMOIN1 + NRHF(ISYMI)*(J - 1)
               KOFF6 = IOFROO(ISYMJK,ISYMI) + ICOFRO(ISYMJ,ISYMK)
     *               + NRHF(ISYMJ)*(K - 1)  + J
               CALL DAXPY(NRHF(ISYMI),-TWO,WORK(KOFF5),1,
     *                    XINTAI(KOFF6),NCOFRO(ISYMJK))
C
  120       CONTINUE
C
C---------------------------------------------
C           Calculate integrals (I j | K del).
C---------------------------------------------
C
            KOFF7  = ILRHSI(ISYMAL) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTI  = MAX(NRHFFR(ISYMI),1)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHF(ISYMJ),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF7),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KMOIN2),NTOTI)
C
C-----------------------------------------------
C           Add contribution to aI-intermediate.
C-----------------------------------------------
C
            DO 130 J = 1,NRHF(ISYMJ)
C
               KOFF8 = KMOIN2 + NRHFFR(ISYMI)*(J - 1)
               KOFF9 = ICOFRF(ISYMJK,ISYMI) + ICOFRO(ISYMJ,ISYMK)
     *               + NRHF(ISYMJ)*(K - 1)  + J
               CALL DAXPY(NRHFFR(ISYMI),-TWO,WORK(KOFF8),1,
     *                    XINAIF(KOFF9),NCOFRO(ISYMJK))
C
  130       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcof1 */
      SUBROUTINE CC_FRCOF1(XINAIF,DSFRO,CMO,WORK,LWORK,ISYDEL)
C
C     Written by Asger Halkier 17/8 - 1998.
C
C     To calculate the coulomb part of the intermediates needed
C     for the frozen-correlated correction to eta-aI.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, EIGHT = 8.0D0)
      DIMENSION XINAIF(*), DSFRO(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
      ISALBE = ISYMJK
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
      KAOINT = 1
      KMOINT = KAOINT + N2BST(ISALBE)
      KEND1  = KMOINT + NCOFRO(ISYMJK)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CC_FRCOF1')
      ENDIF
C
      CALL DZERO(WORK(KAOINT),N2BST(ISALBE))
C
      DO 100 I = 1,NRHFFR(ISYMI)
C
         CALL DZERO(WORK(KMOINT),NCOFRO(ISYMJK))
C
C--------------------------------------------------------
C        Square up integral distribution (al be | I del).
C--------------------------------------------------------
C
         KOFF1 = IDSFRO(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
         CALL CCSD_SYMSQ(DSFRO(KOFF1),ISALBE,WORK(KAOINT))
C
C------------------------------------------------------
C        Inner symmetry loop and final work allocation.
C------------------------------------------------------
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMJK,ISYMJ)
            ISYMAL = ISYMJ
            ISYMBE = ISYMK
C
            KSCR1 = KEND1
            KEND2 = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMK)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT('Insufficient memory for allocation '//
     &                   'in CC_FRCOF1')
            ENDIF
C
            CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NRHFFR(ISYMK))
C
C---------------------------------------------
C           Calculate integrals (j K | I del).
C---------------------------------------------
C
            KOFF2  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF3  = ILRHSI(ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMK),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF2),NTOTAL,CMO(KOFF3),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF4  = ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMJ) + 1
            KOFF5  = KMOINT + ICOFRO(ISYMJ,ISYMK)
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTJ  = MAX(NRHF(ISYMJ),1)
C
            CALL DGEMM('T','N',NRHF(ISYMJ),NRHFFR(ISYMK),NBAS(ISYMAL),
     *                 ONE,CMO(KOFF4),NTOTAL,WORK(KSCR1),NTOTAL,ZERO,
     *                 WORK(KOFF5),NTOTJ)
C
  110    CONTINUE
C
C-------------------------------------------------------
C        Calculate the contribution to the intermediate.
C-------------------------------------------------------
C
         KOFF6 = ICOFRF(ISYMJK,ISYMI) + NCOFRO(ISYMJK)*(I - 1) + 1
         CALL DAXPY(NCOFRO(ISYMJK),EIGHT,WORK(KMOINT),1,XINAIF(KOFF6),1)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_frcof2 */
      SUBROUTINE CC_FRCOF2(XINAIF,DSRHF,CMO,WORK,LWORK,ISYDEL)
C
C     Written by Asger Halkier 17/8 - 1998.
C
C     To calculate the first exchange part of the intermediates
C     needed for the frozen-correlated correction to eta-aI.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINAIF(*), DSRHF(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
C--------------------------------------------------
C     Initial symmetries and work space allocation.
C--------------------------------------------------
C
      ISYMA  = ISYDEL
      ISYMI  = ISYMA
      ISYMJK = 1
C
      IF (NRHFFR(ISYMI) .EQ. 0) RETURN
C
C---------------------------------------------------------------------
C     Loops over third index in integral batch & work space allocation.
C---------------------------------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMK  = MULD2H(ISYMJK,ISYMJ)
         ISYMIK = MULD2H(ISYMI,ISYMK)
         ISALBE = ISYMIK
         ISYMAL = ISYMI
         ISYMBE = ISYMK
C
         IF (NRHFFR(ISYMK) .EQ. 0) GOTO 100
C
         KAOINT = 1
         KSCR1  = KAOINT + N2BST(ISALBE)
         KMOINT = KSCR1  + NBAS(ISYMAL)*NRHFFR(ISYMK)
         KEND1  = KMOINT + NRHFFR(ISYMI)*NRHFFR(ISYMK)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space for allocation '//
     &                'in CC_FRCOF2')
         ENDIF
C
         CALL DZERO(WORK(KAOINT),KEND1)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
C-----------------------------------------------------------
C           Square up integral distribution (al be | j del).
C-----------------------------------------------------------
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE)*(J - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KAOINT))
C
C---------------------------------------------
C           Calculate integrals (I K | j del).
C---------------------------------------------
C
            KOFF2  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF3  = ILRHSI(ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMK),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF2),NTOTAL,CMO(KOFF3),NTOTBE,ZERO,
     *                 WORK(KSCR1),NTOTAL)
C
            KOFF4  = ILRHSI(ISYMAL) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTI  = MAX(NRHFFR(ISYMI),1)
C
            CALL DGEMM('T','N',NRHFFR(ISYMI),NRHFFR(ISYMK),
     *                 NBAS(ISYMAL),ONE,CMO(KOFF4),NTOTAL,WORK(KSCR1),
     *                 NTOTAL,ZERO,WORK(KMOINT),NTOTI)
C
C----------------------------------------------------------
C           Calculate the contribution to the intermediate.
C----------------------------------------------------------
C
            DO 120 K = 1,NRHFFR(ISYMK) 
C
               KOFF5 = KMOINT + NRHFFR(ISYMI)*(K - 1)
               KOFF6 = ICOFRF(ISYMJK,ISYMI) + ICOFRO(ISYMJ,ISYMK)
     *               + NRHF(ISYMJ)*(K - 1) + J
               CALL DAXPY(NRHFFR(ISYMI),-TWO,WORK(KOFF5),1,
     *                    XINAIF(KOFF6),NCOFRO(ISYMJK))
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
