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
C  /* Deck intexp */
      SUBROUTINE INTEXP(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,
     &                  NINDAB,NINDCD,NCCINT,NINTYP,WORK,LWORK,
     &                  ISYMR,ISYMS,ISYMT,ICORBA,ICORBB,ICORBC,ICORBD,
     &                  THRESH,SYMFAC,
     &                  IPRINT,NOPMAT,NODV,EXPECT,LONDON,SUSCEP,DIA2SO,
     &                  ZFS2EL,DDFOCK,SPNORB,IXPAB,IXPCD,IFCTYP)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL NOPMAT, NODV, LONDON, SUSCEP, DDFOCK, EXPECT, DIA2SO,
     &        ZFS2EL, SPNORB
      DIMENSION PSO(*), PSA(*), DMAT(*), AOINT(NCCINT,NINTYP),
     &          ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &          WORK(LWORK), NINDAB(*), NINDCD(*), FT(*), FV(*),
     &          HESSEE(*), IXPAB(*), IXPCD(*), IFCTYP(*)
#include <twocom.h>
#include <twosta.h>
C
C     Allocations
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK

      CALL MEMGET('REAL',KDERIV,NINTYP,WORK,KFREE,LFREE)
      LPAO = NORBA*NORBB*NORBC*NORBD
      CALL MEMGET('REAL',KPAO  ,LPAO  ,WORK,KFREE,LFREE)
      IF (SUSCEP) THEN
         CALL MEMGET('REAL',KPAA  ,LPAO  ,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KPAA  ,0     ,WORK,KFREE,LFREE)
      END IF
      IF (LONDON .AND. DDFOCK) THEN
         CALL MEMGET('REAL',KPINT ,3*NCCINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KQINT ,3*NCCINT,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KPINT ,0       ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KQINT ,0       ,WORK,KFREE,LFREE)
      END IF
      MWINTE = MAX(MWINTE,KFREE)
      LWTOT  = LWTOT + KFREE
      MWTOT  = MAX(MWTOT,LWTOT)
      CALL INTEX1(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,WORK(KDERIV),
     &            WORK(KPAO),WORK(KPAA),WORK(KPINT),
     &            WORK(KQINT),WORK(KFREE),LFREE,
     &            NCCINT,NINTYP,ISYMR,ISYMS,ISYMT,
     &            ICORBA,ICORBB,ICORBC,ICORBD,THRESH,
     &            SYMFAC,IPRINT,NOPMAT,NODV,NUCABQ,NUCCDQ,
     &            NINDAB,NINDCD,EXPECT,LONDON,SUSCEP,DIA2SO,ZFS2EL,
     &            DDFOCK,SPNORB,IXPAB,IXPCD,IFCTYP)
      LWTOT  = LWTOT - KFREE
      CALL MEMREL('INTEXP',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck intex1 */
      SUBROUTINE INTEX1(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,FT,FV,DERIV,PAO,
     &                  PAA,PINT,QINT,WORK,LWORK,NCCINT,NINTYP,
     &                  ISYMR,ISYMS,ISYMT,ICORBA,ICORBB,ICORBC,ICORBD,
     &                  THRESH,SYMFAC,IPRINT,NOPMAT,NODV,NUCABQ,NUCCDQ,
     &                  NINDAB,NINDCD,EXPECT,LONDON,SUSCEP,DIA2SO,
     &                  ZFS2EL,DDFOCK,SPNORB,IXPAB,IXPCD,IFCTYP)
C
C     Calculates expectation values and Fock matrices
C     of differentiated integrals
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (DP25 = 0.25 D00, DP5 = 0.5 D00,
     &           D1 = 1.0 D00, D1P5 = 1.5 D00, D2 = 2.0 D00, 
     &           ZERADD = 1.D-15,
     &           D0 = 0.0 D00, DP125 = 0.125D00, D3=3.0D0, D4=4.0D0)
      INTEGER A, B, C, D, X, Y, Z, Y2, Z2
      LOGICAL NODER, DCMPAB, DCMPCD, NOPMAT, NODV, DV, SUSCEP, DDFOCK,
     &        EXPECT, LONDON, DIA2SO, ZFS2EL, SPNORB
#include <nuclei.h>
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &          DERIV(NINTYP), PSO(*), PSA(*),
     &          PAO(NORBA,NORBB,NORBC,NORBD),
     &          PAA(NORBA,NORBB,NORBC,NORBD),
     &          AOINT(NCCINT,NINTYP),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          DMAT(NBASIS,NBASIS,NDMAT), WORK(LWORK),
     &          FT(NBASIS,NBASIS,*), FV(NBASIS,NBASIS,*), DINT(12),
     &          NCART(12), PINT(NCCINT,3), QINT(NCCINT,3),
     &          DIFAB(3), DIFCD(3), HESSEE(*),
     &          IXPAB(NUCA*NUCB,2), IXPCD(NUCC*NUCD,2),
     &          IFCTYP(NDMAT)
#include <twocom.h>
#include <dftcom.h>
#include <symmet.h>
#include <expcom.h>
#include <expopt.h>
#include <dirprt.h>
#include <infinp.h>
#include <rspprp.h>
#include <esg.h>
#include <ibtfun.h>
      XAND(I) = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IBTAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1
C
      DV = .NOT. NODV
      IF (IPRINT .GT. 9) THEN
         CALL HEADER('Output from INTEXP',-1)
         WRITE (LUPRI, '(A,6L5)') 
     &      ' NODV, NOPMAT, SUSCEP, DDFOCK, DIA2SO, ZFS2EL',
     &        NODV, NOPMAT, SUSCEP, DDFOCK, DIA2SO, ZFS2EL
         WRITE (LUPRI, '(A,3I5)') ' ISYMR/S/T  ', ISYMR,ISYMS,ISYMT
         WRITE (LUPRI, '(A,4I5)') ' NORB ', NORBA,NORBB,NORBC,NORBD
         WRITE (LUPRI, '(A,4I5)') ' ICENT ',ICENT1,ICENT2,ICENT3,ICENT4
         WRITE (LUPRI, '(A,2L5)') ' DIAGAB/CD ', DIAGAB,DIAGCD
         WRITE (LUPRI, '(A,3L5)') ' SHAEQB, SHCEQD, SHABAB ',
     &                              SHAEQB, SHCEQD, SHABAB
         WRITE (LUPRI, '(A,4I5)') ' NHKTA', NHKTA,NHKTB,NHKTC,NHKTD
         WRITE (LUPRI, '(A,4I5)') ' KHKTA', KHKTA,KHKTB,KHKTC,KHKTD
         WRITE (LUPRI, '(A,F12.6)') ' THRESH ', THRESH
         WRITE (LUPRI, '(A,F12.6)') ' SYMFAC ', SYMFAC
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals A '
         WRITE (LUPRI, '(20I5)') (ICORBA(I),I=1, NORBA)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals B '
         WRITE (LUPRI, '(20I5)') (ICORBB(I),I=1, NORBB)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals C '
         WRITE (LUPRI, '(20I5)') (ICORBC(I),I=1, NORBC)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals D '
         WRITE (LUPRI, '(20I5)') (ICORBD(I),I=1, NORBD)
      END IF
C
      IF (DDFOCK .AND. (.NOT.(LONDON .OR. DIA2SO .OR. ZFS2EL))) THEN
         IF (EXPECT) THEN
            IF (TWOCEN) THEN
               NCENTS = 1
               NCART(1) = 3*ICENT1 - 2
               NCART(2) = 3*ICENT1 - 1
               NCART(3) = 3*ICENT1
               NCART(4) = 3*ICENT2 - 2
               NCART(5) = 3*ICENT2 - 1
               NCART(6) = 3*ICENT2
            ELSE IF (THRCEN) THEN
               NCENTS = 2
               NCART(1) = 3*ICENT1 - 2
               NCART(2) = 3*ICENT1 - 1
               NCART(3) = 3*ICENT1
               NCART(4) = 3*ICENT2 - 2
               NCART(5) = 3*ICENT2 - 1
               NCART(6) = 3*ICENT2
               NCART(7) = 3*ICENT3 - 2
               NCART(8) = 3*ICENT3 - 1
               NCART(9) = 3*ICENT3
            ELSE
               NCENTS = 3
               NCART(1)  = 3*ICENT1 - 2
               NCART(2)  = 3*ICENT1 - 1
               NCART(3)  = 3*ICENT1
               NCART(4)  = 3*ICENT2 - 2
               NCART(5)  = 3*ICENT2 - 1
               NCART(6)  = 3*ICENT2
               NCART(7)  = 3*ICENT3 - 2
               NCART(8)  = 3*ICENT3 - 1
               NCART(9)  = 3*ICENT3
               NCART(10) = 3*ICENT4 - 2
               NCART(11) = 3*ICENT4 - 1
               NCART(12) = 3*ICENT4
            END IF
            NCARTS = 3*NCENTS
            NCARTZ = NCARTS + 3
         ELSE
            NCARTS = 3*NOATMS
            NCARTZ = NCARTS
            DO 5 I = 1, NCARTZ
               ICNT = (I + 2)/3
               IATM = ICNTDR(ICNT) - NUCNUM(NCNTDR(ICNT),1) + 1
               NCART(I) = 3*IATM + MOD(I-1,3) - 2
   5        CONTINUE
         END IF
      END IF
      ISYMTS   = IBTXOR(ISYMT,ISYMS)
      IF (DDFOCK .AND. LONDON) THEN
         DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
         DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
         DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
         DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
         DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
         DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0
         DO 10 N = 1, 3
            Y   = NEXT(N)
            Z   = NEXT(Y)
            ABY = DIFAB(Y)
            ABZ = DIFAB(Z)
            CDY = DIFCD(Y)
            CDZ = DIFCD(Z)
            Y2  = Y + 3
            Z2  = Z + 3
            DO 20 I = 1, NCCINT
               AOAB = ABY*AOINT(I,Z ) - ABZ*AOINT(I,Y )
               AOCD = CDY*AOINT(I,Z2) - CDZ*AOINT(I,Y2)
               PINT(I,N) = AOAB + AOCD
               QINT(I,N) = AOAB - AOCD
  20        CONTINUE
  10     CONTINUE
      END IF
C
      SFAC = SYMFAC
      IF (.NOT.SHABAB.AND..NOT.(DIA2SO.OR.SPNORB)) SFAC = D2*SFAC
      IF (.NOT.SHAEQB) SFAC = D2*SFAC
      IF (.NOT.SHCEQD) SFAC = D2*SFAC
C
      THRSH = MAX(THRESH,ZERADD)
      NODER = .TRUE.
      CALL DZERO(DERIV,NINTYP)
C
C     ***** Loop over shell components *****
C
      IAOFF = 1
      DO 100 ICOMPA = 1, KHKTA
         KHKTBB = KHKTB
         IF (DIAGAB) KHKTBB = ICOMPA
         DO 200 ICOMPB = 1, KHKTBB
            DCMPAB = DIAGAB .AND. ICOMPA .EQ. ICOMPB
            FACAB  = D1
            IF (DIAGAB .AND. ICOMPA .NE. ICOMPB) FACAB = D2*FACAB
            DO 300 ICOMPC = 1, KHKTC
               KHKTDD = KHKTD
               IF (DIAGCD) KHKTDD = ICOMPC
               DO 400 ICOMPD = 1, KHKTDD
C
C     Step 2 screening on gradient, hessian, London
C     Always do Step 2 screening /Aug.98,hjaaj
C
                  AOMAX = D0
                  DO I = 1,NINTYP
                     IAOMAX = IDAMAX(NOABCD,AOINT(IAOFF,I),1)-1
                     AOMAX = MAX(AOMAX,ABS(AOINT(IAOFF+IAOMAX,I)))
                  END DO
                  IF (AOMAX .LT. THRSH) GO TO 510
C
                  DCMPCD = DIAGCD .AND. ICOMPC .EQ. ICOMPD
                  FACCD  = D1
                  IF (DIAGCD .AND. ICOMPC .NE. ICOMPD) FACCD = D2*FACCD
                  FCABCD = FACAB*FACCD*SFAC
C
C                 Transform P-matrix block from SO basis to AO basis
C
                  IF (.NOT.NOPMAT) THEN
                     CALL PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     &                           NHKTA,NHKTB,NHKTC,NHKTD,
     &                           KHKTA,KHKTB,KHKTC,MULA,MULB,MULC,MULD,
     &                           NORBA,NORBB,NORBC,NORBD,
     &                           ISYMR,ISYMS,ISYMT)
                  END IF
                  IF (.NOT.NOPMAT .AND. SUSCEP) THEN
                     CALL PBLOCK(PSA,PAA,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     &                           NHKTA,NHKTB,NHKTC,NHKTD,
     &                           KHKTA,KHKTB,KHKTC,MULA,MULB,MULC,MULD,
     &                           NORBA,NORBB,NORBC,NORBD,
     &                           ISYMR,ISYMS,ISYMT)
                  END IF
                  IF (IPRINT .GT. 10) THEN
                     WRITE (LUPRI,'(A,I10)') ' IAOFF ', IAOFF
                     WRITE (LUPRI,'(A,4I5)') ' ICOMP ',
     *                                       ICOMPA,ICOMPB,ICOMPC,ICOMPD
                     WRITE (LUPRI,'(A,2L5)') ' DCMPAB/CD ',DCMPAB,DCMPCD
                  END IF
C
                  INT = IAOFF
C
C                 ***** Loop over contracted functions *****
C
                  DO 500 IORBAB = 1, NORBAB
                     IORBA = NINDAB(IORBAB,1)
                     IORBB = NINDAB(IORBAB,2)
                     A = ICORBA(IORBA) + ICOMPA
                     B = ICORBB(IORBB) + ICOMPB
                     DTAB  = DMAT(A,B,1)
                     IF (DV) DVAB  = DMAT(A,B,2)
                                                      FAB = D1
                     IF (TCONAB .AND. IORBA.NE.IORBB) FAB = D2
                     DO 600 IORBCD = 1, NORBCD
                        IORBC = NINDCD(IORBCD,1)
                        IORBD = NINDCD(IORBCD,2)
                        C = ICORBC(IORBC) + ICOMPC
                        D = ICORBD(IORBD) + ICOMPD
C
                        DTAC = DMAT(C,A,1)
                        DTAD = DMAT(D,A,1)
                        DTBC = DMAT(C,B,1)
                        DTBD = DMAT(D,B,1)
                        DTCD = DMAT(D,C,1)
                        IF (DV) THEN
                           DVAC = DMAT(C,A,2)
                           DVAD = DMAT(D,A,2)
                           DVBC = DMAT(C,B,2)
                           DVBD = DMAT(D,B,2)
                           DVCD = DMAT(D,C,2)
                        END IF
                                                       FCD = D1
                        IF (TCONCD.AND.IORBC.NE.IORBD) FCD = D2
C
C                       ***********************************
C                       ***** Geometrical distortions *****
C                       ***********************************
C
                        IF (.NOT.(LONDON .OR. DIA2SO .OR. ZFS2EL
     &                      .OR. SPNORB)) THEN
C
C                          Expectation values
C                          ==================
C
                           IF (EXPECT) THEN
                              IF (NOPMAT) THEN
                                 PAOVAL = D0
                              ELSE
                                 PAOVAL = PAO(IORBA,IORBB,IORBC,IORBD)
                              ENDIF
                              PFAC = DP5*FAB*FCD*FCABCD 
                              IF (DFTADD) THEN
                                 XFAC = HFXFAC
                              ELSE
                                 XFAC = 1.0D0
                              END IF
                              IF (HFXMU .NE. D0) THEN
                                 PVAL = PFAC*( - HFXFAC*DP25*
     &                                  (DTAC*DTBD + DTAD*DTBC))
                              ELSE IF (NODV) THEN
                                 PVAL = PFAC*(DTAB*DTCD
     &                                - XFAC*DP25*(DTAC*DTBD+DTAD*DTBC)
     &                                + PAOVAL)
                              ELSE IF (HSROHF .OR. 
     &                                (DFTADD .AND. .NOT. NODV)) THEN
                                 PVAL = PFAC*(DTAB*DTCD
     &                                - XFAC*DP25*(DTAC*DTBD + DVAC*DVBD
     &                                      + DTAD*DTBC + DVAD*DVBC))
                              ELSE IF (ESG) THEN 
                                 CALL ESG_PVAL(DFTADD,HFXFAC,
     &                                         PVAL,DMAT,A,B,C,D,PFAC)
                              ELSE
                                 PVAL = PFAC*(DTAB*DTCD - DVAB*DVCD
     &                                - DP25*(DTAC*DTBD - DVAC*DVBD
     &                                      + DTAD*DTBC - DVAD*DVBC)
     &                                + PAOVAL)
                              END IF
                              IF (ABS(PVAL) .GT. THRSH) THEN
                                 IF (EXPGRA) THEN
                                    NODER = .FALSE.
                                    IA = IXPAB(IORBAB,1)
                                    IB = IXPAB(IORBAB,2)
                                    IC = IXPCD(IORBCD,1)
                                    ID = IXPCD(IORBCD,2)
                                    ALPGRD(IA) = ALPGRD(IA)
     &                                         + PVAL*AOINT(INT,1)
                                    ALPGRD(IB) = ALPGRD(IB)
     &                                         + PVAL*AOINT(INT,2)
                                    ALPGRD(IC) = ALPGRD(IC)
     &                                         + PVAL*AOINT(INT,3)
                                    ALPGRD(ID) = ALPGRD(ID)
     &                                         + PVAL*AOINT(INT,4)
                                 ELSE
                                    NODER = .FALSE.
                                    DO 700 I = 1, NINTYP
                                     DERIV(I)=DERIV(I)+PVAL*AOINT(INT,I)
  700                               CONTINUE
                                 END IF
                              END IF
                           END IF
C
C                          Fock matrices
C                          =============
C
                           IF (DDFOCK) THEN
                              CALL DCOPY(NCARTS,AOINT(INT,1),NCCINT,
     &                                   DINT,1)
                              IF (EXPECT) THEN
                                 DINT(NCARTS+1)=-DSUM(NCENTS,DINT(1),3)
                                 DINT(NCARTS+2)=-DSUM(NCENTS,DINT(2),3)
                                 DINT(NCARTS+3)=-DSUM(NCENTS,DINT(3),3)
                              END IF
C
C                         950210-hjaaj: new SKLFC1 symmetrizes with
C                         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
                              DFAC = FAB*FCD*FCABCD
                              EFAC = -HFXFAC*DP25*DFAC
                              IF (HFXMU .NE. D0) DFAC = D0
C
C                             Total Fock matrix
C
                              DCD  = DFAC*DTCD
                              DAB  = DFAC*DTAB
                              DBD  = EFAC*DTBD
                              DBC  = EFAC*DTBC
                              DAD  = EFAC*DTAD
                              DAC  = EFAC*DTAC
                              DO 800 I = 1, NCARTZ
                                 N = NCART(I)
                                 FT(D,C,N) = FT(D,C,N) + DINT(I)*DAB
                                 FT(B,A,N) = FT(B,A,N) + DINT(I)*DCD
                                 FT(C,A,N) = FT(C,A,N) + DINT(I)*DBD
                                 FT(D,A,N) = FT(D,A,N) + DINT(I)*DBC
                                 FT(C,B,N) = FT(C,B,N) + DINT(I)*DAD
                                 FT(D,B,N) = FT(D,B,N) + DINT(I)*DAC
  800                         CONTINUE
C
C                             Active Fock matrix
C
                              IF (DV) THEN
                                 DCD  = DFAC*DVCD
                                 DAB  = DFAC*DVAB
                                 DBD  = EFAC*DVBD
                                 DBC  = EFAC*DVBC
                                 DAD  = EFAC*DVAD
                                 DAC  = EFAC*DVAC
                                 DO 810 I = 1, NCARTZ
                                    N = NCART(I)
                                    FV(B,A,N) = FV(B,A,N) + DINT(I)*DCD
                                    FV(D,C,N) = FV(D,C,N) + DINT(I)*DAB
                                    FV(C,A,N) = FV(C,A,N) + DINT(I)*DBD
                                    FV(D,A,N) = FV(D,A,N) + DINT(I)*DBC
                                    FV(C,B,N) = FV(C,B,N) + DINT(I)*DAD
                                    FV(D,B,N) = FV(D,B,N) + DINT(I)*DAC
  810                            CONTINUE
                              END IF
                           END IF
                           IF (IPRINT .GT. 30) THEN
                              WRITE(LUPRI,'(//A,4I5,D12.4)')
     &                              ' IORBA... PAOVAL',
     &                              IORBA,IORBB,IORBC,IORBD,PAOVAL
                              WRITE(LUPRI,'(//A,4I5)')'A,B,C,D',A,B,C,D
                              WRITE(LUPRI,'(//A,6F12.6)') 'DTAB,...',
     &                              DTAB,DTAC,DTAD,DTBC,DTBD,DTCD
                              WRITE(LUPRI,'(//A,6F12.6)') 'DVAB,...',
     &                              DVAB,DVAC,DVAD,DVBC,DVBD,DVCD
                              WRITE (LUPRI, '(A,3F12.6)')
     &                              ' Fab, Fcd, FCabab',FAB, FCD, FCABCD
                              WRITE (LUPRI, '(A,F12.6)')' PVAL ',PVAL
                              IF (DDFOCK) THEN
                                 WRITE (LUPRI,'(A,12F6.3)') ' DINT ',
     &                               (DINT(I),I=1,NCARTZ)
                              END IF
                           END IF
C
C                       **************************
C                       ***** Magnetic field *****
C                       **************************
C
                        ELSE IF (LONDON) THEN
                           FACTOR = DP125*FAB*FCD*FCABCD
                           IF (SUSCEP) THEN
                              IF (NOPMAT) THEN
                                 PAOVL1 = D0
                                 PAOVL2 = D0
                              ELSE
                                 PAOVL1 = PAO(IORBA,IORBB,IORBC,IORBD)
                                 PAOVL2 = PAA(IORBA,IORBB,IORBC,IORBD)
                              ENDIF
                              IF (HFXMU .NE. D0) THEN
                                 PVAL1  = FACTOR*(
     &                                  - HFXFAC*DP25*
     &                                    (DTAC*DTBD + DTAD*DTBC)
     &                                  + PAOVL1)
                                 PVAL2  = FACTOR*(HFXFAC*DP25*
     &                                    (DTAC*DTBD - DTAD*DTBC)
     &                                  + PAOVL2)
                              ELSE IF (NODV) THEN
                                 PVAL1  = FACTOR*(DTAB*DTCD
     &                                  - HFXFAC*DP25*
     &                                    (DTAC*DTBD + DTAD*DTBC)
     &                                  + PAOVL1)
                                 PVAL2  = FACTOR*(HFXFAC*DP25*
     &                                    (DTAC*DTBD - DTAD*DTBC)
     &                                  + PAOVL2)
                              ELSE
                                 PVAL1  = FACTOR*
     &                                    (DTAB*DTCD - DVAB*DVCD
     &                                    - HFXFAC*DP25*(DTAC*DTBD
     &                                          - DVAC*DVBD + DTAD*DTBC
     &                                          - DVAD*DVBC)
     &                                  + PAOVL1)
                                 PVAL2  = FACTOR*(HFXFAC*DP25*
     &                                    (DTAC*DTBD - DVAC*DVBD
     &                                     - DTAD*DTBC + DVAD*DVBC)
     &                                  + PAOVL2)
                              END IF
                              IF (ABS(PVAL1) .GT. THRSH) THEN
                                 NODER = .FALSE.
                                 DO 900 I = 7, 18
                                    DERIV(I)=DERIV(I)+PVAL1*AOINT(INT,I)
  900                            CONTINUE
                              END IF
                              IF (ABS(PVAL2) .GT. THRSH) THEN
                                 NODER = .FALSE.
                                 DO 910 I = 19, 27
                                    DERIV(I)=DERIV(I)+PVAL2*AOINT(INT,I)
  910                            CONTINUE
                              END IF
                           END IF
                           IF (DDFOCK) THEN
                              DFAC = DP125*FAB*FCD*FCABCD
                              EFAC = DP5*HFXFAC*DFAC
                              IF (HFXMU .NE. D0) DFAC = D0 
C
C                             Total Fock matrix
C
                              DAB  = DFAC*DTAB
                              DCD  = DFAC*DTCD
                              DAC  = EFAC*DTAC
                              DAD  = EFAC*DTAD
                              DBC  = EFAC*DTBC
                              DBD  = EFAC*DTBD
                              DO 920 N = 1, 3
                                 ABCD = PINT(INT,N)
                                 ABDC = QINT(INT,N)
                                 FT(D,C,N) = FT(D,C,N) - (ABCD-ABDC)*DAB
                                 FT(B,A,N) = FT(B,A,N) - (ABCD+ABDC)*DCD
                                 FT(C,A,N) = FT(C,A,N) + ABDC*DBD
                                 FT(D,A,N) = FT(D,A,N) + ABCD*DBC
                                 FT(C,B,N) = FT(C,B,N) - ABCD*DAD
                                 FT(D,B,N) = FT(D,B,N) - ABDC*DAC
  920                         CONTINUE
C
C                             Active Fock matrix
C
                              IF (DV) THEN
                                 DAB  = DFAC*DVAB
                                 DCD  = DFAC*DVCD
                                 DAC  = EFAC*DVAC
                                 DAD  = EFAC*DVAD
                                 DBC  = EFAC*DVBC
                                 DBD  = EFAC*DVBD
                                 DO 930 N = 1, 3
                                    ABCD = PINT(INT,N)
                                    ABDC = QINT(INT,N)
                                    FV(B,A,N) = FV(B,A,N)
     &                                        - (ABCD+ABDC)*DCD
                                    FV(D,C,N) = FV(D,C,N)
     &                                        - (ABCD-ABDC)*DAB
                                    FV(C,A,N) = FV(C,A,N) + ABDC*DBD
                                    FV(D,A,N) = FV(D,A,N) + ABCD*DBC
                                    FV(C,B,N) = FV(C,B,N) - ABCD*DAD
                                    FV(D,B,N) = FV(D,B,N) - ABDC*DAC
  930                            CONTINUE
                              END IF
                           END IF
                           IF (IPRINT .GT. 30) THEN
                              WRITE(LUPRI,'(//A,4I5,2D12.4)')
     &                             ' IORBA... PAOVAL',
     &                             IORBA,IORBB,IORBC,IORBD,PAOVL1,PAOVL2
                              WRITE(LUPRI,'(//A,4I5)')'A,B,C,D',A,B,C,D
                              WRITE(LUPRI,'(//A,6F12.6)') 'DTAB,...',
     &                              DTAB,DTAC,DTAD,DTBC,DTBD,DTCD
                              WRITE(LUPRI,'(//A,6F12.6)') 'DVAB,...',
     &                              DVAB,DVAC,DVAD,DVBC,DVBD,DVCD
                              WRITE (LUPRI, '(A,3F12.6)')
     &                              ' Fab, Fcd, FCabab',FAB, FCD, FCABCD
                              WRITE (LUPRI, '(A,2F12.6)')' PVAL ',
     &                               PVAL1,PVAL2
                              IF (DDFOCK) THEN
                                 WRITE (LUPRI,'(A,12F6.3)') ' DINT ',
     &                               (DINT(I),I=1,NCARTZ)
                              END IF
                           END IF
C
C                          Spin-orbit Fock matrices
C                          ========================
C
                        ELSE IF (SPNORB) THEN
                           FACTOR = DP25*FAB*FCD*FCABCD
                           DO N = 1, NDMAT 
                           IF (IFCTYP(N).NE.0) THEN
C
C                             Coulomb contribution
C
                              DTOTAB = DMAT(A,B,N) + DMAT(B,A,N)
                              DTOTCD = DMAT(D,C,N) - DMAT(C,D,N)
                              AOC = FACTOR*AOINT(INT,ABS(IFCTYP(N)))
C
C                                triplet case
C
                              IF (IFCTYP(N).LT.0) THEN 
                                 FT(A,B,N) = FT(A,B,N) + DTOTCD*AOC*D2
                                 FT(B,A,N) = FT(B,A,N) + DTOTCD*AOC*D2
                                 FT(C,D,N) = FT(C,D,N) + DTOTAB*AOC 
                                 FT(D,C,N) = FT(D,C,N) - DTOTAB*AOC 
C
C                                singlet case
C
                              ELSE
                                 FT(A,B,N) = FT(A,B,N) + DTOTCD*AOC
                                 FT(B,A,N) = FT(B,A,N) + DTOTCD*AOC
                                 FT(C,D,N) = FT(C,D,N) + DTOTAB*AOC*D2
                                 FT(D,C,N) = FT(D,C,N) - DTOTAB*AOC*D2
                              END IF
C
C                             exchange contribution
C
                              IF (HFXFAC .NE. 0.0D0) THEN
                                 AOX = D1P5*HFXFAC*AOC
                                 FT(A,D,N) = FT(A,D,N) - DMAT(B,C,N)*AOX
                                 FT(D,A,N) = FT(D,A,N) + DMAT(C,B,N)*AOX
                                 FT(A,C,N) = FT(A,C,N) + DMAT(B,D,N)*AOX
                                 FT(C,A,N) = FT(C,A,N) - DMAT(D,B,N)*AOX
                                 FT(B,D,N) = FT(B,D,N) - DMAT(A,C,N)*AOX
                                 FT(D,B,N) = FT(D,B,N) + DMAT(C,A,N)*AOX
                                 FT(B,C,N) = FT(B,C,N) + DMAT(A,D,N)*AOX
                                 FT(C,B,N) = FT(C,B,N) - DMAT(D,A,N)*AOX
                              END IF
                           END IF
                           END DO
                        ELSE IF (DIA2SO) THEN
C
C                          Expectation values of diamagnetic spin-orbit
C                          ============================================
C
                           IF (NOPMAT) THEN
                              PAOVAL =  3*(
     &                           2*DVAB*DVCD - DVAD*DVBC - DVAC*DVBD
     &                           )
                           ELSE
                              PAOVAL = PAO(IORBA,IORBB,IORBC,IORBD)
                           ENDIF
                           PVAL= 
     &                         D2*DTCD*DVAB + D4*DTAB*DVCD 
     &                         - D3*DP5*(
     &                              DTBC*DVAD + DTBD*DVAC
     &                            + DTAD*DVBC + DTAC*DVBD
     &                              ) + PAOVAL
                           PVAL = DP5*FAB*FCD*FCABCD*PVAL
                           IF (ABS(PVAL) .GT. THRSH) THEN
                              NODER = .FALSE.
                              DO 799 I = 1, NINTYP
                                 DERIV(I)=DERIV(I)+PVAL*AOINT(INT,I)
 799                          CONTINUE
                           END IF
                           IF (IPRINT .GT. 30) THEN
                              WRITE(LUPRI,'(//A,4I5,D12.4)')
     &                             ' IORBA... PAOVAL',
     &                             IORBA,IORBB,IORBC,IORBD,PAOVAL
                              WRITE(LUPRI,'(//A,4I5)')'A,B,C,D',A,B,C,D
C                             WRITE(LUPRI, '(//A,I5,A/3F12.6)')
C    &                           'AOINT(',INT,',*)',
C    &                           (AOINT(INT,I),I=1,NINTYP)
                              WRITE(LUPRI,'(3F12.6)')
     &                           (AOINT(INT,I),I=1,NINTYP)
                              WRITE(LUPRI,'(//A,6F12.6)') 'DTAB,...',
     &                              DTAB,DTAC,DTAD,DTBC,DTBD,DTCD
                              WRITE(LUPRI,'(//A,6F12.6)') 'DVAB,...',
     &                              DVAB,DVAC,DVAD,DVBC,DVBD,DVCD
                              WRITE (LUPRI, '(A,3F12.6)')
     &                              ' Fab, Fcd, FCabab',FAB, FCD, FCABCD
                              WRITE (LUPRI, '(A,F12.6)')' PVAL ',
     &                               PVAL
                              WRITE (LUPRI, '(A)') ' PVAL*AOINT '
                              WRITE (LUPRI,'(3F12.6)')
     &                           (PVAL*AOINT(INT,I),I=1,NINTYP)
                              WRITE (LUPRI, '(A)') ' DERIV '
                              WRITE (LUPRI,'(3F12.6)')
     &                           (DERIV(I),I=1,NINTYP)
                              IF (DDFOCK) THEN
                                 WRITE (LUPRI,'(A,12F6.3)') ' DINT ',
     &                               (DINT(I),I=1,NCARTZ)
                              END IF
                           END IF
                        ELSE
C
C                          Expectation values of spin-spin coupling
C                          ============================================
C
                           IF (NOPMAT) THEN
                              PAOVAL = D0
                           ELSE
                              PAOVAL = PAO(IORBA,IORBB,IORBC,IORBD)
                           ENDIF
                           IF (.NOT.NODV) THEN
                              PAOVAL = DP5*(DVAB*DVCD 
     &                           - DP5*(DVAD*DVBC+DVAC*DVBD))
                           END IF
                           PVAL = DP5*FAB*FCD*FCABCD*PAOVAL
                           IF (ABS(PVAL) .GT. THRSH) THEN
                              NODER = .FALSE.
                              DO 899 I = 1, NINTYP
                                 DERIV(I)=DERIV(I)+PVAL*AOINT(INT,I)
 899                          CONTINUE
                           END IF
                           IF (IPRINT .GT. 30) THEN
                              WRITE(LUPRI,'(//A,4I5,D12.4)')
     &                             ' IORBA... PAOVAL',
     &                             IORBA,IORBB,IORBC,IORBD,PAOVAL
                              WRITE(LUPRI,'(//A,4I5)')'A,B,C,D',A,B,C,D
C                             WRITE(LUPRI, '(//A,I5,A/3F12.6)')
C    &                           'AOINT(',INT,',*)',
C    &                           (AOINT(INT,I),I=1,NINTYP)
                              WRITE(LUPRI,'(3F12.6)')
     &                           (AOINT(INT,I),I=1,NINTYP)
                              WRITE(LUPRI,'(//A,6F12.6)') 'DTAB,...',
     &                              DTAB,DTAC,DTAD,DTBC,DTBD,DTCD
                              WRITE(LUPRI,'(//A,6F12.6)') 'DVAB,...',
     &                              DVAB,DVAC,DVAD,DVBC,DVBD,DVCD
                              WRITE (LUPRI, '(A,3F12.6)')
     &                              ' Fab, Fcd, FCabab',FAB, FCD, FCABCD
                              WRITE (LUPRI, '(A,F12.6)')' PVAL ',
     &                               PVAL
                              WRITE (LUPRI, '(A)') ' PVAL*AOINT '
                              WRITE (LUPRI,'(3F12.6)')
     &                           (PVAL*AOINT(INT,I),I=1,NINTYP)
                              WRITE (LUPRI, '(A)') ' DERIV '
                              WRITE (LUPRI,'(3F12.6)')
     &                           (DERIV(I),I=1,NINTYP)
                              IF (DDFOCK) THEN
                                 WRITE (LUPRI,'(A,12F6.3)') ' DINT ',
     &                               (DINT(I),I=1,NCARTZ)
                              END IF
                           END IF
                        END IF
                        INT = INT + 1
  600                CONTINUE
  500             CONTINUE
  510             IAOFF = IAOFF + NOABCD
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      IF (.NOT.NODER) THEN
         IF (.NOT.LONDON .AND. EXPECT .AND. .NOT.EXPGRA) THEN
            CALL DEROUT(DERIV,HESSEE,NINTYP,IPRINT)
         END IF
         IF (SUSCEP) THEN
            CALL LNDOUT(DERIV,ISYMR,ISYMT,ISYMTS,IPRINT)
         END IF
         IF (DIA2SO) THEN
            CALL SO2OUT(DERIV,IPRINT)
         END IF
         IF (ZFS2EL) THEN
            CALL ZFSOUT(DERIV,IPRINT)
         END IF
      END IF
      RETURN
      END
C  /* Deck derout */
      SUBROUTINE DEROUT(DERIV,HESSEE,NINTYP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      INTEGER AX, AY, AZ, BX, BY, BZ, CX, CY, CZ, DX, DY, DZ
      DIMENSION DERIV(NINTYP), HESSEE(MXCOOR,MXCOOR)
C
#include <nuclei.h>
#include <expcom.h>
#include <energy.h>
C
      IF (IPRINT .GE.10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) ICENT1, ICENT2, ICENT3, ICENT4
      END IF
C
C     ***************************
C     ***** Two-Center Case *****
C     ***************************
C
      IF (TWOCEN) THEN
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
C
C           A electron-repulsion gradient elements:
C
            GRADEE(AX) = GRADEE(AX) + DAX
            GRADEE(AY) = GRADEE(AY) + DAY
            GRADEE(AZ) = GRADEE(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            GRADEE(BX) = GRADEE(BX) - DAX
            GRADEE(BY) = GRADEE(BY) - DAY
            GRADEE(BZ) = GRADEE(BZ) - DAZ
         END IF
C
C        ***** Hessian *****
C
         IF (DERTWO) THEN
C
            AAXX = DERIV(4)
            AAXY = DERIV(5)
            AAXZ = DERIV(6)
            AAYY = DERIV(7)
            AAYZ = DERIV(8)
            AAZZ = DERIV(9)
C
C           A-A electron-repulsion Hessian elements:
C
            HESSEE(AX,AX) = HESSEE(AX,AX) + AAXX
            HESSEE(AX,AY) = HESSEE(AX,AY) + AAXY
            HESSEE(AX,AZ) = HESSEE(AX,AZ) + AAXZ
            HESSEE(AY,AY) = HESSEE(AY,AY) + AAYY
            HESSEE(AY,AZ) = HESSEE(AY,AZ) + AAYZ
            HESSEE(AZ,AZ) = HESSEE(AZ,AZ) + AAZZ
C
C           A-B electron-repulsion Hessian elements:
C
            HESSEE(AX,BX) = HESSEE(AX,BX) - AAXX
            HESSEE(AX,BY) = HESSEE(AX,BY) - AAXY
            HESSEE(AX,BZ) = HESSEE(AX,BZ) - AAXZ
            HESSEE(AY,BX) = HESSEE(AY,BX) - AAXY
            HESSEE(AY,BY) = HESSEE(AY,BY) - AAYY
            HESSEE(AY,BZ) = HESSEE(AY,BZ) - AAYZ
            HESSEE(AZ,BX) = HESSEE(AZ,BX) - AAXZ
            HESSEE(AZ,BY) = HESSEE(AZ,BY) - AAYZ
            HESSEE(AZ,BZ) = HESSEE(AZ,BZ) - AAZZ
C
C           B-B electron-repulsion Hessian elements:
C
            HESSEE(BX,BX) = HESSEE(BX,BX) + AAXX
            HESSEE(BX,BY) = HESSEE(BX,BY) + AAXY
            HESSEE(BX,BZ) = HESSEE(BX,BZ) + AAXZ
            HESSEE(BY,BY) = HESSEE(BY,BY) + AAYY
            HESSEE(BY,BZ) = HESSEE(BY,BZ) + AAYZ
            HESSEE(BZ,BZ) = HESSEE(BZ,BZ) + AAZZ
         END IF
C
C     *****************************
C     ***** Three-Center Case *****
C     *****************************
C
      ELSE IF (THRCEN) THEN
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
         CX = 3*ICENT3 - 2
         CY = 3*ICENT3 - 1
         CZ = 3*ICENT3
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
C
C           A electron-repulsion gradient elements:
C
            GRADEE(AX) = GRADEE(AX) + DAX
            GRADEE(AY) = GRADEE(AY) + DAY
            GRADEE(AZ) = GRADEE(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            GRADEE(BX) = GRADEE(BX) + DBX
            GRADEE(BY) = GRADEE(BY) + DBY
            GRADEE(BZ) = GRADEE(BZ) + DBZ
C
C           C electron-repulsion gradient elements:
C
            GRADEE(CX) = GRADEE(CX) - (DAX + DBX)
            GRADEE(CY) = GRADEE(CY) - (DAY + DBY)
            GRADEE(CZ) = GRADEE(CZ) - (DAZ + DBZ)
         END IF
C
C        ***** Hessian *****
C
         IF (DERTWO) THEN
C
            AAXX = DERIV(7)
            AAXY = DERIV(8)
            AAXZ = DERIV(9)
            AAYY = DERIV(10)
            AAYZ = DERIV(11)
            AAZZ = DERIV(12)
C
            BBXX = DERIV(13)
            BBXY = DERIV(14)
            BBXZ = DERIV(15)
            BBYY = DERIV(16)
            BBYZ = DERIV(17)
            BBZZ = DERIV(18)
C
            ABXX = DERIV(19)
            ABXY = DERIV(20)
            ABXZ = DERIV(21)
            ABYX = DERIV(22)
            ABYY = DERIV(23)
            ABYZ = DERIV(24)
            ABZX = DERIV(25)
            ABZY = DERIV(26)
            ABZZ = DERIV(27)
C
C           A-A electron-repulsion Hessian elements:
C
            HESSEE(AX,AX) = HESSEE(AX,AX) + AAXX
            HESSEE(AX,AY) = HESSEE(AX,AY) + AAXY
            HESSEE(AX,AZ) = HESSEE(AX,AZ) + AAXZ
            HESSEE(AY,AY) = HESSEE(AY,AY) + AAYY
            HESSEE(AY,AZ) = HESSEE(AY,AZ) + AAYZ
            HESSEE(AZ,AZ) = HESSEE(AZ,AZ) + AAZZ
C
C           A-B electron-repulsion Hessian elements:
C
            HESSEE(AX,BX) = HESSEE(AX,BX) + ABXX
            HESSEE(AX,BY) = HESSEE(AX,BY) + ABXY
            HESSEE(AX,BZ) = HESSEE(AX,BZ) + ABXZ
            HESSEE(AY,BX) = HESSEE(AY,BX) + ABYX
            HESSEE(AY,BY) = HESSEE(AY,BY) + ABYY
            HESSEE(AY,BZ) = HESSEE(AY,BZ) + ABYZ
            HESSEE(AZ,BX) = HESSEE(AZ,BX) + ABZX
            HESSEE(AZ,BY) = HESSEE(AZ,BY) + ABZY
            HESSEE(AZ,BZ) = HESSEE(AZ,BZ) + ABZZ
C
C           A-C electron-repulsion Hessian elements:
C
            HESSEE(AX,CX) = HESSEE(AX,CX) - (AAXX + ABXX)
            HESSEE(AX,CY) = HESSEE(AX,CY) - (AAXY + ABXY)
            HESSEE(AX,CZ) = HESSEE(AX,CZ) - (AAXZ + ABXZ)
            HESSEE(AY,CX) = HESSEE(AY,CX) - (AAXY + ABYX)
            HESSEE(AY,CY) = HESSEE(AY,CY) - (AAYY + ABYY)
            HESSEE(AY,CZ) = HESSEE(AY,CZ) - (AAYZ + ABYZ)
            HESSEE(AZ,CX) = HESSEE(AZ,CX) - (AAXZ + ABZX)
            HESSEE(AZ,CY) = HESSEE(AZ,CY) - (AAYZ + ABZY)
            HESSEE(AZ,CZ) = HESSEE(AZ,CZ) - (AAZZ + ABZZ)
C
C           B-B electron-repulsion Hessian elements:
C
            HESSEE(BX,BX) = HESSEE(BX,BX) + BBXX
            HESSEE(BX,BY) = HESSEE(BX,BY) + BBXY
            HESSEE(BX,BZ) = HESSEE(BX,BZ) + BBXZ
            HESSEE(BY,BY) = HESSEE(BY,BY) + BBYY
            HESSEE(BY,BZ) = HESSEE(BY,BZ) + BBYZ
            HESSEE(BZ,BZ) = HESSEE(BZ,BZ) + BBZZ
C
C           B-C electron-repulsion Hessian elements:
C
            HESSEE(BX,CX) = HESSEE(BX,CX) - (BBXX + ABXX)
            HESSEE(BX,CY) = HESSEE(BX,CY) - (BBXY + ABYX)
            HESSEE(BX,CZ) = HESSEE(BX,CZ) - (BBXZ + ABZX)
            HESSEE(BY,CX) = HESSEE(BY,CX) - (BBXY + ABXY)
            HESSEE(BY,CY) = HESSEE(BY,CY) - (BBYY + ABYY)
            HESSEE(BY,CZ) = HESSEE(BY,CZ) - (BBYZ + ABZY)
            HESSEE(BZ,CX) = HESSEE(BZ,CX) - (BBXZ + ABXZ)
            HESSEE(BZ,CY) = HESSEE(BZ,CY) - (BBYZ + ABYZ)
            HESSEE(BZ,CZ) = HESSEE(BZ,CZ) - (BBZZ + ABZZ)
C
C           C-C electron-repulsion Hessian elements:
C
            HESSEE(CX,CX) = HESSEE(CX,CX) + (AAXX + ABXX + ABXX + BBXX)
            HESSEE(CX,CY) = HESSEE(CX,CY) + (AAXY + ABXY + ABYX + BBXY)
            HESSEE(CX,CZ) = HESSEE(CX,CZ) + (AAXZ + ABXZ + ABZX + BBXZ)
            HESSEE(CY,CY) = HESSEE(CY,CY) + (AAYY + ABYY + ABYY + BBYY)
            HESSEE(CY,CZ) = HESSEE(CY,CZ) + (AAYZ + ABYZ + ABZY + BBYZ)
            HESSEE(CZ,CZ) = HESSEE(CZ,CZ) + (AAZZ + ABZZ + ABZZ + BBZZ)
         END IF
C
C     ****************************
C     ***** Four-Center Case *****
C     ****************************
C
      ELSE
         AX = 3*ICENT1 - 2
         AY = 3*ICENT1 - 1
         AZ = 3*ICENT1
         BX = 3*ICENT2 - 2
         BY = 3*ICENT2 - 1
         BZ = 3*ICENT2
         CX = 3*ICENT3 - 2
         CY = 3*ICENT3 - 1
         CZ = 3*ICENT3
         DX = 3*ICENT4 - 2
         DY = 3*ICENT4 - 1
         DZ = 3*ICENT4
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
            DCX = DERIV(7)
            DCY = DERIV(8)
            DCZ = DERIV(9)
C
C           A electron-repulsion gradient elements:
C
            GRADEE(AX) = GRADEE(AX) + DAX
            GRADEE(AY) = GRADEE(AY) + DAY
            GRADEE(AZ) = GRADEE(AZ) + DAZ
C
C           B electron-repulsion gradient elements:
C
            GRADEE(BX) = GRADEE(BX) + DBX
            GRADEE(BY) = GRADEE(BY) + DBY
            GRADEE(BZ) = GRADEE(BZ) + DBZ
C
C           C electron-repulsion gradient elements:
C
            GRADEE(CX) = GRADEE(CX) + DCX
            GRADEE(CY) = GRADEE(CY) + DCY
            GRADEE(CZ) = GRADEE(CZ) + DCZ
C
C           D electron-repulsion gradient elements:
C
            GRADEE(DX) = GRADEE(DX) - (DAX + DBX + DCX)
            GRADEE(DY) = GRADEE(DY) - (DAY + DBY + DCY)
            GRADEE(DZ) = GRADEE(DZ) - (DAZ + DBZ + DCZ)
         END IF
C
C        ***** Hessian *****
C
         IF (DERTWO) THEN
C
            AAXX = DERIV(10)
            AAXY = DERIV(11)
            AAXZ = DERIV(12)
            AAYY = DERIV(13)
            AAYZ = DERIV(14)
            AAZZ = DERIV(15)
C
            BBXX = DERIV(16)
            BBXY = DERIV(17)
            BBXZ = DERIV(18)
            BBYY = DERIV(19)
            BBYZ = DERIV(20)
            BBZZ = DERIV(21)
C
            CCXX = DERIV(22)
            CCXY = DERIV(23)
            CCXZ = DERIV(24)
            CCYY = DERIV(25)
            CCYZ = DERIV(26)
            CCZZ = DERIV(27)
C
            ABXX = DERIV(28)
            ABXY = DERIV(29)
            ABXZ = DERIV(30)
            ABYX = DERIV(31)
            ABYY = DERIV(32)
            ABYZ = DERIV(33)
            ABZX = DERIV(34)
            ABZY = DERIV(35)
            ABZZ = DERIV(36)
C
            ACXX = DERIV(37)
            ACXY = DERIV(38)
            ACXZ = DERIV(39)
            ACYX = DERIV(40)
            ACYY = DERIV(41)
            ACYZ = DERIV(42)
            ACZX = DERIV(43)
            ACZY = DERIV(44)
            ACZZ = DERIV(45)
C
            BCXX = DERIV(46)
            BCXY = DERIV(47)
            BCXZ = DERIV(48)
            BCYX = DERIV(49)
            BCYY = DERIV(50)
            BCYZ = DERIV(51)
            BCZX = DERIV(52)
            BCZY = DERIV(53)
            BCZZ = DERIV(54)
C
C           A-A electron-repulsion Hessian elements:
C
            HESSEE(AX,AX) = HESSEE(AX,AX) + AAXX
            HESSEE(AX,AY) = HESSEE(AX,AY) + AAXY
            HESSEE(AX,AZ) = HESSEE(AX,AZ) + AAXZ
            HESSEE(AY,AY) = HESSEE(AY,AY) + AAYY
            HESSEE(AY,AZ) = HESSEE(AY,AZ) + AAYZ
            HESSEE(AZ,AZ) = HESSEE(AZ,AZ) + AAZZ
C
C           A-B electron-repulsion Hessian elements:
C
            HESSEE(AX,BX) = HESSEE(AX,BX) + ABXX
            HESSEE(AX,BY) = HESSEE(AX,BY) + ABXY
            HESSEE(AX,BZ) = HESSEE(AX,BZ) + ABXZ
            HESSEE(AY,BX) = HESSEE(AY,BX) + ABYX
            HESSEE(AY,BY) = HESSEE(AY,BY) + ABYY
            HESSEE(AY,BZ) = HESSEE(AY,BZ) + ABYZ
            HESSEE(AZ,BX) = HESSEE(AZ,BX) + ABZX
            HESSEE(AZ,BY) = HESSEE(AZ,BY) + ABZY
            HESSEE(AZ,BZ) = HESSEE(AZ,BZ) + ABZZ
C
C           A-C electron-repulsion Hessian elements:
C
            HESSEE(AX,CX) = HESSEE(AX,CX) + ACXX
            HESSEE(AX,CY) = HESSEE(AX,CY) + ACXY
            HESSEE(AX,CZ) = HESSEE(AX,CZ) + ACXZ
            HESSEE(AY,CX) = HESSEE(AY,CX) + ACYX
            HESSEE(AY,CY) = HESSEE(AY,CY) + ACYY
            HESSEE(AY,CZ) = HESSEE(AY,CZ) + ACYZ
            HESSEE(AZ,CX) = HESSEE(AZ,CX) + ACZX
            HESSEE(AZ,CY) = HESSEE(AZ,CY) + ACZY
            HESSEE(AZ,CZ) = HESSEE(AZ,CZ) + ACZZ
C
C           A-D electron-repulsion Hessian elements:
C
            HESSEE(AX,DX) = HESSEE(AX,DX) - (AAXX + ABXX + ACXX)
            HESSEE(AX,DY) = HESSEE(AX,DY) - (AAXY + ABXY + ACXY)
            HESSEE(AX,DZ) = HESSEE(AX,DZ) - (AAXZ + ABXZ + ACXZ)
            HESSEE(AY,DX) = HESSEE(AY,DX) - (AAXY + ABYX + ACYX)
            HESSEE(AY,DY) = HESSEE(AY,DY) - (AAYY + ABYY + ACYY)
            HESSEE(AY,DZ) = HESSEE(AY,DZ) - (AAYZ + ABYZ + ACYZ)
            HESSEE(AZ,DX) = HESSEE(AZ,DX) - (AAXZ + ABZX + ACZX)
            HESSEE(AZ,DY) = HESSEE(AZ,DY) - (AAYZ + ABZY + ACZY)
            HESSEE(AZ,DZ) = HESSEE(AZ,DZ) - (AAZZ + ABZZ + ACZZ)
C
C           B-B electron-repulsion Hessian elements:
C
            HESSEE(BX,BX) = HESSEE(BX,BX) + BBXX
            HESSEE(BX,BY) = HESSEE(BX,BY) + BBXY
            HESSEE(BX,BZ) = HESSEE(BX,BZ) + BBXZ
            HESSEE(BY,BY) = HESSEE(BY,BY) + BBYY
            HESSEE(BY,BZ) = HESSEE(BY,BZ) + BBYZ
            HESSEE(BZ,BZ) = HESSEE(BZ,BZ) + BBZZ
C
C           B-C electron-repulsion Hessian elements:
C
            HESSEE(BX,CX) = HESSEE(BX,CX) + BCXX
            HESSEE(BX,CY) = HESSEE(BX,CY) + BCXY
            HESSEE(BX,CZ) = HESSEE(BX,CZ) + BCXZ
            HESSEE(BY,CX) = HESSEE(BY,CX) + BCYX
            HESSEE(BY,CY) = HESSEE(BY,CY) + BCYY
            HESSEE(BY,CZ) = HESSEE(BY,CZ) + BCYZ
            HESSEE(BZ,CX) = HESSEE(BZ,CX) + BCZX
            HESSEE(BZ,CY) = HESSEE(BZ,CY) + BCZY
            HESSEE(BZ,CZ) = HESSEE(BZ,CZ) + BCZZ
C
C           B-D electron-repulsion Hessian elements:
C
            HESSEE(BX,DX) = HESSEE(BX,DX) - (ABXX + BBXX + BCXX)
            HESSEE(BX,DY) = HESSEE(BX,DY) - (ABYX + BBXY + BCXY)
            HESSEE(BX,DZ) = HESSEE(BX,DZ) - (ABZX + BBXZ + BCXZ)
            HESSEE(BY,DX) = HESSEE(BY,DX) - (ABXY + BBXY + BCYX)
            HESSEE(BY,DY) = HESSEE(BY,DY) - (ABYY + BBYY + BCYY)
            HESSEE(BY,DZ) = HESSEE(BY,DZ) - (ABZY + BBYZ + BCYZ)
            HESSEE(BZ,DX) = HESSEE(BZ,DX) - (ABXZ + BBXZ + BCZX)
            HESSEE(BZ,DY) = HESSEE(BZ,DY) - (ABYZ + BBYZ + BCZY)
            HESSEE(BZ,DZ) = HESSEE(BZ,DZ) - (ABZZ + BBZZ + BCZZ)
C
C           C-C electron-repulsion Hessian elements:
C
            HESSEE(CX,CX) = HESSEE(CX,CX) + CCXX
            HESSEE(CX,CY) = HESSEE(CX,CY) + CCXY
            HESSEE(CX,CZ) = HESSEE(CX,CZ) + CCXZ
            HESSEE(CY,CY) = HESSEE(CY,CY) + CCYY
            HESSEE(CY,CZ) = HESSEE(CY,CZ) + CCYZ
            HESSEE(CZ,CZ) = HESSEE(CZ,CZ) + CCZZ
C
C           C-D electron-repulsion Hessian elements:
C
            HESSEE(CX,DX) = HESSEE(CX,DX) - (ACXX + BCXX + CCXX)
            HESSEE(CX,DY) = HESSEE(CX,DY) - (ACYX + BCYX + CCXY)
            HESSEE(CX,DZ) = HESSEE(CX,DZ) - (ACZX + BCZX + CCXZ)
            HESSEE(CY,DX) = HESSEE(CY,DX) - (ACXY + BCXY + CCXY)
            HESSEE(CY,DY) = HESSEE(CY,DY) - (ACYY + BCYY + CCYY)
            HESSEE(CY,DZ) = HESSEE(CY,DZ) - (ACZY + BCZY + CCYZ)
            HESSEE(CZ,DX) = HESSEE(CZ,DX) - (ACXZ + BCXZ + CCXZ)
            HESSEE(CZ,DY) = HESSEE(CZ,DY) - (ACYZ + BCYZ + CCYZ)
            HESSEE(CZ,DZ) = HESSEE(CZ,DZ) - (ACZZ + BCZZ + CCZZ)
C
C           D-D electron-repulsion Hessian elements:
C
            HESSEE(DX,DX) = HESSEE(DX,DX) + (AAXX + ABXX + ACXX
     &                + ABXX + BBXX + BCXX + ACXX + BCXX + CCXX)
            HESSEE(DX,DY) = HESSEE(DX,DY) + (AAXY + ABXY + ACXY
     &                + ABYX + BBXY + BCXY + ACYX + BCYX + CCXY)
            HESSEE(DX,DZ) = HESSEE(DX,DZ) + (AAXZ + ABXZ + ACXZ
     &                + ABZX + BBXZ + BCXZ + ACZX + BCZX + CCXZ)
            HESSEE(DY,DY) = HESSEE(DY,DY) + (AAYY + ABYY + ACYY
     &                + ABYY + BBYY + BCYY + ACYY + BCYY + CCYY)
            HESSEE(DY,DZ) = HESSEE(DY,DZ) + (AAYZ + ABYZ + ACYZ
     &                + ABZY + BBYZ + BCYZ + ACZY + BCZY + CCYZ)
            HESSEE(DZ,DZ) = HESSEE(DZ,DZ) + (AAZZ + ABZZ + ACZZ
     &                + ABZZ + BBZZ + BCZZ + ACZZ + BCZZ + CCZZ)
         END IF
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) ICENT1, ICENT2, ICENT3, ICENT4
         WRITE (LUPRI,'(A,3I5)') ' AX/Y/Z ',AX,AY,AZ
         WRITE (LUPRI,'(A,3I5)') ' BX/Y/Z ',BX,BY,BZ
         IF (.NOT.TWOCEN) THEN
         WRITE (LUPRI,'(A,3I5)') ' CX/Y/Z ',CX,CY,CZ
         IF (.NOT.THRCEN) THEN
         WRITE (LUPRI,'(A,3I5)') ' DX/Y/Z ',DX,DY,DZ
         END IF
         END IF
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, NINTYP)
         NCDEP3 = 3*NUCDEP
         IF (IPRINT .GE. 10) THEN
            IF (DERONE) THEN
               WRITE(LUPRI, 4010)
               WRITE(LUPRI, 4020) (GRADEE(I),I=1,NCDEP3)
            END IF
            IF (DERTWO) THEN
               WRITE(LUPRI, 4030)
               DO 800 I = 1,NCDEP3
                  WRITE (LUPRI, 4020)
     *                     (HESSEE(I,J) + HESSEE(J,I), J = 1,I)
  800          CONTINUE
            END IF
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine DEROUT >>>>>>>>',/)
 1110 FORMAT (2X,'TWOCEN,...:  ',2L5)
 1120 FORMAT (2X,'ICENT1/2/3/4:',4I7)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
 4010 FORMAT (//,1X,' Two-electron integral gradient ',/)
 4020 FORMAT (6F12.6)
 4030 FORMAT (//,1X,' Two-electron integral Hessian ',/)
      END
C  /* Deck so2out */
      SUBROUTINE SO2OUT(DERIV,IPRINT)
C
C     O.Vahtras and K.Ruud, May 1998
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>      
      DIMENSION DERIV(*)
#include <symmet.h>
#include <gtensor.h>
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Subroutine SO2OUT',-1)
      END IF
      IX = 1 !IPTAX(1,1)
      IY = 2 !IPTAX(2,1)
      IZ = 3 !IPTAX(3,1)
      GC2(IX,IX) = GC2(IX,IX) + DERIV(5) + DERIV(9)
      GC2(IY,IY) = GC2(IY,IY) + DERIV(1) + DERIV(9)
      GC2(IZ,IZ) = GC2(IZ,IZ) + DERIV(1) + DERIV(5)
      IF (ISYMAX(IX,1).EQ.ISYMAX(IY,1)) THEN
         GC2(IX,IY) = GC2(IX,IY) + DERIV(2)
         GC2(IY,IX) = GC2(IY,IX) + DERIV(4)
      END IF
      IF (ISYMAX(IX,1).EQ.ISYMAX(IZ,1)) THEN
         GC2(IX,IZ) = GC2(IX,IZ) + DERIV(3)
         GC2(IZ,IX) = GC2(IZ,IX) + DERIV(7)
      END IF
      IF (ISYMAX(IY,1).EQ.ISYMAX(IZ,1)) THEN
         GC2(IY,IZ) = GC2(IY,IZ) + DERIV(6)
         GC2(IZ,IY) = GC2(IZ,IY) + DERIV(8)
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, 9)
         WRITE(LUPRI,'(/A)') ' Unfinished two-electron spin-orbit '//
     &                       'contributions'
         CALL OUTPUT(GC2,1,3,1,3,3,3,1,LUPRI)
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine SO2OUT >>>>>>>>',/)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
      END
C  /* Deck lndout */
      SUBROUTINE LNDOUT(DERIV,ISYMR,ISYMT,ISYMTS,IPRINT)
C
C     tuh Nov 24 92
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      INTEGER X1, Y1, Z1, X2, Y2, Z2, X1S, X2S
      DIMENSION DERIV(*), G11(3,3), G22(3,3), G12(3,3)
      DIMENSION DIFAB(3), DIFCD(3)
#include <twocom.h>
#include <symmet.h>
#include <suscpt.h>
#include <ibtfun.h>
      XAND(I) = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IBTAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine LNDOUT',-1)
         WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
C
      DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
      DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
      DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
      DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
      DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
      DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0
C
      G11(1,1) = DERIV (7)
      G11(1,2) = DERIV (8)
      G11(1,3) = DERIV (9)
      G11(2,1) = DERIV (8)
      G11(2,2) = DERIV(10)
      G11(2,3) = DERIV(11)
      G11(3,1) = DERIV (9)
      G11(3,2) = DERIV(11)
      G11(3,3) = DERIV(12)
C
      G22(1,1) = DERIV(13)
      G22(1,2) = DERIV(14)
      G22(1,3) = DERIV(15)
      G22(2,1) = DERIV(14)
      G22(2,2) = DERIV(16)
      G22(2,3) = DERIV(17)
      G22(3,1) = DERIV(15)
      G22(3,2) = DERIV(17)
      G22(3,3) = DERIV(18)
C
      CALL DCOPY(9,DERIV(19),1,G12,1)
C
      DO 100 X1 = 1, 3
         Y1   = NEXT (X1)
         Z1   = NEXT (Y1)
         ABY1 = DIFAB(Y1)
         ABZ1 = DIFAB(Z1)
         CDY1 = DIFCD(Y1)
         CDZ1 = DIFCD(Z1)
         X1S  = IPTAX(X1,2)
         DO 200 X2 = X1, 3
         IF (ISYMAX(X1,2).EQ.ISYMAX(X2,2)) THEN
            Y2   = NEXT (X2)
            Z2   = NEXT (Y2)
            ABY2 = DIFAB(Y2)
            ABZ2 = DIFAB(Z2)
            CDY2 = DIFCD(Y2)
            CDZ2 = DIFCD(Z2)
            X2S  = IPTAX(X2,2)
            SUS2EL(X1S,X2S) = SUS2EL(X1S,X2S)
     &                      - ABY1*ABY2*G11(Z1,Z2)
     &                      + ABY1*ABZ2*G11(Z1,Y2)
     &                      + ABZ1*ABY2*G11(Y1,Z2)
     &                      - ABZ1*ABZ2*G11(Y1,Y2)
C
     &                      - CDY1*CDY2*G22(Z1,Z2)
     &                      + CDY1*CDZ2*G22(Z1,Y2)
     &                      + CDZ1*CDY2*G22(Y1,Z2)
     &                      - CDZ1*CDZ2*G22(Y1,Y2)
C
     &                      - ABY1*CDY2*G12(Z2,Z1)
     &                      + ABY1*CDZ2*G12(Y2,Z1)
     &                      + ABZ1*CDY2*G12(Z2,Y1)
     &                      - ABZ1*CDZ2*G12(Y2,Y1)
C
     &                      - ABY2*CDY1*G12(Z1,Z2)
     &                      + ABY2*CDZ1*G12(Y1,Z2)
     &                      + ABZ2*CDY1*G12(Z1,Y2)
     &                      - ABZ2*CDZ1*G12(Y1,Y2)
         END IF
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, 27)
         WRITE(LUPRI,'(/A)') ' Unfinished two-electron susceptibilities'
         CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,'(/A)') ' G11 '
            CALL OUTPUT(G11,1,3,1,3,3,3,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' G22 '
            CALL OUTPUT(G22,1,3,1,3,3,3,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' G12 '
            CALL OUTPUT(G12,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine LNDOUT >>>>>>>>',/)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
      END
C  /* Deck zfsout */
      SUBROUTINE ZFSOUT(DERIV,IPRINT)
C
C     K.Ruud, September 2000
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>      
      DIMENSION DERIV(*)
#include <symmet.h>
#include <zfs.h>
      PARAMETER (D2=2.0D0)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Subroutine ZFSOUT',-1)
      END IF
      IX = 1 !IPTAX(1,1)
      IY = 2 !IPTAX(2,1)
      IZ = 3 !IPTAX(3,1)
      ZFS(IX,IX) = ZFS(IX,IX) + DERIV(1) + DERIV(7) + D2*DERIV(13) 
      ZFS(IY,IY) = ZFS(IY,IY) + DERIV(4) + DERIV(10)+ D2*DERIV(17) 
      ZFS(IZ,IZ) = ZFS(IZ,IZ) + DERIV(6) + DERIV(12)+ D2*DERIV(21) 
      IF (ISYMAX(IZ,2).EQ.0) THEN
         ZFS(IX,IY) = ZFS(IX,IY) + DERIV(2) + DERIV(8) + DERIV(14)
     &                        + DERIV(16)
         ZFS(IY,IX) = ZFS(IY,IX) + DERIV(2) + DERIV(8) + DERIV(14)
     &                        + DERIV(16)
      END IF
      IF (ISYMAX(IY,2).EQ.0) THEN
         ZFS(IX,IZ) = ZFS(IX,IZ) + DERIV(3) + DERIV(9) + DERIV(15)
     &                        + DERIV(19)
         ZFS(IZ,IX) = ZFS(IZ,IX) + DERIV(3) + DERIV(9) + DERIV(15)
     &                        + DERIV(19)
      END IF
      IF (ISYMAX(IX,2).EQ.0) THEN
         ZFS(IY,IZ) = ZFS(IY,IZ) + DERIV(5) + DERIV(11) + DERIV(18)
     &                        + DERIV(20)
         ZFS(IZ,IY) = ZFS(IZ,IY) + DERIV(5) + DERIV(11) + DERIV(18)
     &                        + DERIV(20)
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, 9)
         WRITE(LUPRI,'(/A)') ' Unfinished two-electron zero-field '//
     &                       'splittings'
         CALL OUTPUT(ZFS,1,3,1,3,3,3,1,LUPRI)
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine SO2OUT >>>>>>>>',/)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
      END
C  /* Deck pblock */
      SUBROUTINE PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     &                  NHKTA,NHKTB,NHKTC,NHKTD,
     &                  KHKTA,KHKTB,KHKTC,MULA,MULB,MULC,MULD,
     &                  NORBA,NORBB,NORBC,NORBD,ISYMR,ISYMS,ISYMT)
#include <implicit.h>
C
C     *********************************************************
C     ***  Transform a block of the P-matrix from SO basis  ***
C     ***  to AO basis.  This is done for fixed component   ***
C     ***  indices ICOMPx and fixed symmetry operations.    ***
C     *********************************************************
C
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <symmet.h>
      DIMENSION PSO(*), PAO(*)
#include <ibtfun.h>
      IAR(I,J,K,L) = KHKTA*(KHKTB*(KHKTC*(L-1)+K-1)+J-1)+I
      IRR(I,J,K,L) = MULTA*(MULTB*(MULTC*L+K)+J)+I+1
      ISYMTS = IBTXOR(ISYMT,ISYMS)
      MULTA = MULT(MULA)
      MULTB = MULT(MULB)
      MULTC = MULT(MULC)
      MULTD = MULT(MULD)
      NCABCD = NORBA*NORBB*NORBC*NORBD
      MULTOT = IRR(MULTA-1,MULTB-1,MULTC-1,MULTD-1)
      ITYNA = ISYMAO(NHKTA,ICOMPA)
      ITYNB = ISYMAO(NHKTB,ICOMPB)
      ITYNC = ISYMAO(NHKTC,ICOMPC)
      ITYND = ISYMAO(NHKTD,ICOMPD)
      CALL DZERO(PAO,NCABCD)
C
C     Loop over irreps
C
      IRCNTA = -1
      DO 100 IREPA = 0, MAXREP
      IF(IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
         IRCNTA = IRCNTA + 1
         IRCNTB = -1
         DO 200 IREPB = 0, MAXREP
         IF(IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
            SIGNB = PT(IBTAND(ISYMR, IBTXOR(IREPB,ITYNB)))
            IRCNTB = IRCNTB + 1
            IRCNTC = -1
            DO 300 IREPC = 0, MAXREP
            IF(IBTAND(MULC,IBTXOR(IREPC,ITYNC)) .EQ. 0) THEN
               SIGNBC = SIGNB*PT(IBTAND(ISYMT, IBTXOR(IREPC,ITYNC)))
               IRPABC = IBTXOR(IREPA,IBTXOR(IREPB,IREPC))
               IRCNTC = IRCNTC + 1
               IRCNTD = -1
               DO 400 IREPD = 0,MAXREP
               IF (IBTAND(MULD,IBTXOR(IREPD,ITYND)) .EQ. 0) THEN
                  IRCNTD = IRCNTD + 1
                  IF (IBTXOR(IREPD,IRPABC) .EQ. 0) THEN
                     SIGN =SIGNBC*PT(IBTAND(ISYMTS,IBTXOR(IREPD,ITYND)))
                     IOFF = NCABCD*(MULTOT
     &                      *(IAR(ICOMPA,ICOMPB,ICOMPC,ICOMPD) - 1)
     &                      + IRR(IRCNTA,IRCNTB,IRCNTC,IRCNTD) - 1)
                     DO 500 I = 1,NCABCD
                         PAO(I) = PAO(I) + SIGN*PSO(IOFF+I)
  500                CONTINUE
                  END IF
               END IF
  400          CONTINUE
            END IF
  300       CONTINUE
         END IF
  200    CONTINUE
      END IF
  100 CONTINUE
      RETURN
      END
