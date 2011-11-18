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
C  /* Deck erifck */
      SUBROUTINE ERIFCK(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     *                  CCFBT,INDXBT,FCKTMP,WORK,LWORK)
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <dummy.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
C
C	I don't know if IFCTYP and ISYMDM is used properly in this 
C	routine at the current stage, K.Ruud-jan96
C
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <veclen.h>
#include <odbtch.h>
#include <symmet.h>
#include <infpar.h>
      DIMENSION FMAT(*), DMAT(*), WORK(LWORK),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*), 
     &	        FCKTMP(*)
C
      CALL QENTER('ERIFCK')
C
      IF (SLAVE) THEN
         IPRINT = IPRFCK
      ELSE
         CALL TIMER('START ',TIMSTR,TIMEND)
C
C        Initialization in ER2INI
C
         CALL ER2INI
C
         IPRINT = MAX(IPRERI,IPRFCK)
      END IF
C
      THRSH  = MAX(THRS,1.00D-15)
      NDMAT  = NDMT
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
      WRTINT = .FALSE.
      FCKINT = .TRUE.
C
C     AO batches
C     ==========
C
      CALL SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
C
C     Density matrix (not needed yet)
C     ===============================
C
C     KDAO  = 1
C     KLAST = KDAO + N2BASX*NDMAT
C     IF(KLAST.GT.LWORK) CALL STOPIT('ERIFCK','DSOTAO',KLAST,LWORK)
C     LWRK  = LWORK - KLAST + 1
C     DO IMAT = 1,NDMAT
C        JDAO = KDAO + (IMAT-1)*N2BASX
C        CALL DSOTAO(DMAT(1,IMAT),WORK(JDAO),NBAST,ISYMDM,IPRINT)
C     END DO
C
C     OD batches
C     ==========
C
C     This subroutine returns several arrays for each electron
C     starting at addresses K????1 and K????2. These are to be
C     transferred to ODCDRV.
C
      CALL ODBCHS(KODCL1,KODCL2,
     &            KODBC1,KODBC2,KRDBC1,KRDBC2,
     &            KODPP1,KODPP2,KRDPP1,KRDPP2,
     &            KFREE,LFREE,CCFBT,WORK,
     &            LWORK,IPRINT)
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWORK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C
      ICALL = 0
      CALL GETDST(ICALL,ICALL,IPRINT)
C
C     Select integrals to be calculated
C     =================================
C
      CALL PICKAO(IPRINT)
C
C     Information about distributions
C     ===============================
C
      CALL ERIDSI(INDXBT,IPRINT)
#if defined (VAR_VECTOR)
      ICHUNK = MAX(IVECLN/NDMT,1)
      CALL DZERO(FCKTMP,ICHUNK*NDMT*(NBASE + NODD)*NBASE)
#endif
C
      KLAST = KFREE
      LWRK  = LFREE
C
C     Calculate integrals
C     ===================
C
      IF (SLAVE) THEN
#if defined (VAR_VECTOR)
         CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               FCKTMP,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &	             INDXBT,WORK(KLAST),LWRK,IPRINT)
#else
         CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               FMAT,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &               INDXBT,WORK(KLAST),LWRK,IPRINT)
#endif
      ELSE
         IF (.NOT.INTSKP) THEN
#if defined (VAR_VECTOR)
            CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &                  WORK(KODBC1),WORK(KODBC2),
     &                  WORK(KRDBC1),WORK(KRDBC2),
     &                  WORK(KODPP1),WORK(KODPP2),
     &                  WORK(KRDPP1),WORK(KRDPP2),
     &                  FCKTMP,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &	                INDXBT,WORK(KLAST),LWRK,IPRINT)
#else
            CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &                  WORK(KODBC1),WORK(KODBC2),
     &                  WORK(KRDBC1),WORK(KRDBC2),
     &                  WORK(KODPP1),WORK(KODPP2),
     &                  WORK(KRDPP1),WORK(KRDPP2),
     &                  FMAT,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &                  INDXBT,WORK(KLAST),LWRK,IPRINT)
#endif
C
C           Error message in case of insufficient memory
C
            IF (.NOT.MEMOK) THEN
               WRITE (LUPRI,'(//A,3(/A,I12))')
     &            ' Not enough memory for this run of ERIFCK.',
     &            ' Available memory in ERIFCK:',LWORK,
     &            ' Required memory for ERIFCK:',LWORK + MEMADD,
     &            ' Increase memory (LWORK) by:',MEMADD
               WRITE (LUPRI,'(/A,2I5)')
     &            ' Memory requirements largest for OD classes :',
     &              MODAB,MODCD
               CALL QUIT('Insufficient memory in ERIFCK.')
            END IF
         END IF
#if defined (VAR_VECTOR)
         IOFF = 0
         DO I = 1, ICHUNK
            DO J = 1, NDMT
               DO L = 1, NBASE
                  DO K = 1, NBASE
C                    FMAT(K,L,J) = FMAT(K,L,J) + FCKTMP(IOFF + K)
                     FMAT(K+(L-1)*NBASE+(J-1)*NBASE*NBASE) =
     &                    FMAT(K+(L-1)*NBASE+(J-1)*NBASE*NBASE) +
     &                    FCKTMP(IOFF + K)
                  END DO
                  IOFF = IOFF + NBASE+NODD
               END DO
            END DO
         END DO
#endif
C
C        Symmetrize Fock matrix (Not needed!)
C        ====================================
C
C NECgh980314  Kenneth found on 98/03/14 that this is the IFCTYP=3 bug.
C NECgh980314  Now, we scale with 2.0 instead of symmetrizing.
C NECgh980314  CALL ERISFK(FMAT,NBASE,NDMT)
C NECgh980505 Instead of doing this call, we adjust the FAC in the FCKCON call.
C NECgh980505  call dscal(nbase*nbase*ndmt,2.0d0,fmat,1)
C
C        Print densities and Fock matrix
C        ===============================
C
         IF (IPRINT.GT.4) THEN
            CALL HEADER('Density and Fock matrices in ERIFCK',-1)
            KSTR = 1
            DO I = 1, NDMT
               WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
               CALL OUTPUT(DMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                     NBASE,1,LUPRI)
               WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
               CALL OUTPUT(FMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                     NBASE,1,LUPRI)
               KSTR = KSTR + NBASE*NBASE
            END DO
         END IF
C
         CALL TIMER('ERIFCK',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
      END IF
C
      CALL QEXIT('ERIFCK')
      RETURN
      END
C  /* Deck erifok */
      SUBROUTINE ERIFOK(SO,IPNTCR,IODDCC,IPNTUV,FMAT,DMAT,NDMT,
     &                  IFCTYP,CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <eridst.h>
      DIMENSION SO(*), IPNTCR(MAXBCH,4), CCFBT(*), INDXBT(*),
     &          IPNTUV(KC2MAX,0:NRDER,2), WORK(LWORK), IODDCC(NRTOP),
     &          FMAT(NBASE,NBASE), DMAT(NBASE,NBASE), IFCTYP(NDMT)
#include <cbieri.h>
#include <aobtch.h>
#include <ericom.h>
#include <eribuf.h>
#include <symmet.h>
#include <hertop.h>
C
C     Allocation for ERIFOK
C
      LBIN   = NCCT*KHKTA*KHKTB*KHKTC*KHKTD
      KBIN   = 1
      KIBIN  = KBIN   +  LBIN
      KINDEX = KIBIN  + (4*LBIN - 1)/IRAT + 1
      KLAST  = KINDEX + (4*LBIN - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIFOK',' ',KLAST,LWORK)
      CALL ERIFO1(SO,WORK(KINDEX),IPNTCR,IODDCC,IPNTUV,
     &            WORK(KBIN),WORK(KIBIN),LBIN,FMAT,DMAT,NDMT,
     &            IFCTYP,CCFBT,INDXBT,IPRINT)
C
      RETURN
      END
C  /* Deck erifo1 */
      SUBROUTINE ERIFO1(SO,INDEX,IPNTCR,IODDCC,IPNTUV,
     &                  BIN,IBIN,LBIN,FMAT,DMAT,NDMT,IFCTYP,
     &                  CCFBT,INDXBT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
C
#include <eridst.h>
      DIMENSION SO(*), INDEX(LBIN,4), CCFBT(*), INDXBT(*),
     &          IPNTCR(MAXBCH,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          BIN(LBIN), IBIN(LBIN,4),
     &          FMAT(NBASE,NBASE,NDMT), DMAT(NBASE,NBASE,NDMT),
     &          IFCTYP(NDMT)
#include <cbieri.h>
#include <ericom.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIFO1',-1)
C
C     Collect (non-zero) integrals and attach indices
C     ===============================================
C
      CALL ERINDF(SO,INDEX,IPNTCR,IODDCC,IPNTUV,
     &            BIN,IBIN,LBIN,CCFBT,INDXBT,INT,IPRINT)
C
C     Construct Fock matrix contribution
C     ==================================
C
      CALL FOKDI1(FMAT,DMAT,NDMT,IFCTYP,BIN,IBIN,LBIN,INT,
     &            ABS(IFITDM),IPRINT)
C
      RETURN
      END
C  /* Deck erindf */
      SUBROUTINE ERINDF(SO,INDEX,IPNTCR,IODDCC,IPNTUV,
     &                  BIN,IBIN,LBIN,CCFBT,INDXBT,INT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      PARAMETER (D1 = 1.0D0, DP5=0.50D0)
      INTEGER A, B, C, D, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ,
     &        AEQB, CEQD, PEQQ
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD),
     &          INDEX(NCCS,4), CCFBT(MXPRIM*MXCONT),
     &          IPNTCR(MAXBCH,4), INDXBT(MXSHEL*MXCONT,0:7),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          BIN(LBIN), IBIN(LBIN,4),
     &          IPNRST(0:7,3),
     &          IADCMP(MXAQN,MXAQN,2)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
C
#include <ibtfun.h>
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERINDF',-1)
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
      CALL PRPREP(DOREP(0,3),NHKTC,KHKTC,ISTBLC)
      CALL PRPREP(DOREP(0,4),NHKTD,KHKTD,ISTBLD)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP C  ',(DOREP(I,3),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP D  ',(DOREP(I,4),I=0,MAXREP)
      END IF
C
      NCCS1 = (NPQBCS - NPPBCS)*NCTFAB*NCTFCD
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      INT = 0
C
      DO A = 0, MAXREP
      IF (DOREP(A,1)) THEN
      DO B = 0, MAXREP
      IF (DOREP(B,2)) THEN
      DO C = 0, MAXREP
      IF (DOREP(C,3) .AND. DOREP(IBTXOR(IBTXOR(A,B),C),4)) THEN
         D = IBTXOR(IBTXOR(A,B),C)
         CD = IBTXOR(C,D)
C
         IF (DIAGAB .AND. B.GT.A) GO TO 100
         IF (DIAGCD .AND. D.GT.C) GO TO 100
C
         CTRIAB = DIAGAB .AND. A.EQ.B
         CTRICD = DIAGCD .AND. C.EQ.D
         CTRIAC = A.EQ.C
         CTRIBD = B.EQ.D
         CTRIPQ = A.EQ.C .AND. B.EQ.D
C
         R = IPNRST(B,1)
         S = IPNRST(D,2)
         T = IPNRST(CD,3)
C
         CALL ERIPNT(INDEX,A,B,C,D,IPNTCR,INDXBT,1)
C
         IF (NCCS1 .GT. 0) THEN
C
            IA   = -1
            MAXB = KHKTB
            MAXD = KHKTD
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  AEQB = CTRIAB .AND. ICMPA.EQ.ICMPB
                  DO ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
                     IF (CTRICD) MAXD = ICMPC
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        IF (IODDAB .EQ. IODDCD) THEN
                           CEQD = CTRICD .AND. ICMPC.EQ.ICMPD
                           FAC = D1
                           IF (AEQB) FAC = DP5*FAC
                           IF (CEQD) FAC = DP5*FAC
                           DO I = 1, NCCS1
                              SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD)
                              IF (ABS(SOABCD) .GT. THRSH) THEN
                                 INT = INT + 1
                                 BIN(INT) = FAC*SOABCD
                                 IBIN(INT,1) = INDEX(I,1) + IA
                                 IBIN(INT,2) = INDEX(I,2) + IB
                                 IBIN(INT,3) = INDEX(I,3) + IC
                                 IBIN(INT,4) = INDEX(I,4) + ID
                              END IF
                           END DO
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
         END IF
C
         IF (NCCS2 .GT. 0) THEN
C
            IF (C.GT.A .OR. (C.EQ.A .AND. D.GT.B)) GO TO 200
C
            MAXB = KHKTB
            MAXC = KHKTC
            MAXD = KHKTD
C
            IA = -1
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  AEQB = CTRIAB .AND. ICMPA.EQ.ICMPB
                  IF (CTRIAC) MAXC = ICMPA
                  DO ICMPC = 1, MAXC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
                     IF (CTRIPQ .AND. ICMPA.EQ.ICMPC) THEN
                        MAXD = ICMPB
                     ELSE
                        MAXD = KHKTD
                        IF (CTRICD) MAXD = ICMPC
                     END IF
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        IF (IODDAB .EQ. IODDCD) THEN
                           CEQD = CTRICD .AND. ICMPC.EQ.ICMPD
                           PEQQ = CTRIPQ .AND. ICMPA.EQ.ICMPC
     &                                   .AND. ICMPB.EQ.ICMPD 
                           FAC = D1
                           IF (AEQB) FAC = DP5*FAC
                           IF (CEQD) FAC = DP5*FAC
                           IF (PEQQ) FAC = DP5*FAC
                           DO I = NCCS1 + 1, NCCS
                              SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD)
                              IF (ABS(SOABCD) .GT. THRSH) THEN
                                 INT = INT + 1
                                 BIN(INT) = FAC*SOABCD
                                 IBIN(INT,1) = INDEX(I,1) + IA
                                 IBIN(INT,2) = INDEX(I,2) + IB
                                 IBIN(INT,3) = INDEX(I,3) + IC
                                 IBIN(INT,4) = INDEX(I,4) + ID
                              END IF
                           END DO
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
  200       CONTINUE
         END IF
  100    CONTINUE
      END IF
      END DO
      END IF
      END DO
      END IF
      END DO
      RETURN
      END
C  /* Deck fokdi1 */
      SUBROUTINE FOKDI1(FMAT,DMAT,NDMT,IFCTYP,BUF,IBUF,
     &                  LBIN,LENGTH,IFIT_DMAT,IPRINT)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C
C     This subroutine adds derivative two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
      PARAMETER (D4 = 4.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      INTEGER P, Q, R, S
#include <priunit.h>
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST,NDMT,*), DMAT(NBAST,NBAST,NDMT),
     &          BUF(LBIN), IBUF(LBIN,4), IFCTYP(NDMT)
C
      DO I = 1, NDMT
	 IX = IFCTYP(I) / 10
	 IY = MOD (IFCTYP(I),10)
C   FAC account for the different integrals in eri and twoint.
C   twoint are 4times larger, but FAC seems only to be 2(?), at least for
C   IFCTYP = 13
	 IF      (IFCTYP(I).EQ.13) THEN
C NECgh980505  We adjust the FAC, since we do not symmetrize/scal anymore.
C NECgh980505  FAC = D2
	       FAC = D4
         ELSE IF (IFCTYP(I).EQ.3)  THEN
C NECgh980505  FAC = DP5 
 	       FAC = D1
         ELSE
            FAC = D1
             WRITE(LUPRI,*) '*** WARNING!!! This value for IFCTYP is '//
     &               'probably not correctly implemented! ***'
         END IF
         IF (IFIT_DMAT.EQ.0) THEN
C...........Ordinary Fock matrix with exact density
            CALL FCKCON(FMAT,DMAT,NDMT,I,BUF,IBUF,LBIN,LBIN,LENGTH,IX,
     &                  IY,FAC)
         ELSE
C...........Construction of special Fock matrix for density fitting
C...........The multiplication factor should be one in this case.
            FAC = D1
            CALL DF_FCKCON (FMAT,DMAT,FMAT,DMAT,NDMT,I,BUF,IBUF,LBIN,
     &                      LBIN,LENGTH,IX,IY,FAC,IFIT_DMAT)
         ENDIF
      END DO
      RETURN
      END
C  /* Deck erisfk */
      SUBROUTINE ERISFK(FMAT,NBASE,NDMT)
#include <implicit.h>
#include <priunit.h>
      DIMENSION FMAT(NBASE,NBASE,NDMT)
      DO IDM = 1, NDMT
         DO I = 1, NBASE
         DO J = 1, I
            FMT = FMAT(I,J,IDM) + FMAT(J,I,IDM)
            FMAT(I,J,IDM) = FMT
            FMAT(J,I,IDM) = FMT
         END DO
         END DO
      END DO
      RETURN
      END
