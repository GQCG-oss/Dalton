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
C  /* Deck rxr */
      SUBROUTINE RXR(R2AM,V2AM,SF,TF,Q11,NOCDIM,NSPAIR)
C
C     This subroutine evaluates the matrix elements
C
C     Q(ijkl) = Sum_ab (ia|r12|jb) * (ka|(K1 + K2)r12|lb)
C
C     On input, the array R2AM contains the integrals (ia|r12|jb).
C     NOCDIM is the number of (nonfrozen) occupied orbitals and NSPAIR
C     is the number of pairs of (nonfrozen) occupied orbitals.
C
C     On output, the arrays SF and TF contain the matrices Q(ijkl) for
C     singlet and triplet coupled pairs, respectively.
C
C     V2AM (of length 2*NT2AMX) and Q11 (of length 2*NOCDIM**4) are used
C     for scratch space.
C
C     The one-particle matrix X is read from disk (file name AUXFCK).
C     This is the primary+secondary exchange matrix in the orthonormal basis.
C     It is assumed that this matrix is diagonal in the secondary space.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>  
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0,
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      PARAMETER (THRDIA = 1D-9)
      DIMENSION R2AM(*), SF(*), TF(*), ISB(8)
      REAL*8  Q11(NOCDIM,NOCDIM,NOCDIM,NOCDIM), V2AM(*),
     *        RR, VV, XBD, XAC, SFAC
      LOGICAL AVIRT, BVIRT, LDUM
      INTEGER IDUM
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      ISB(1) = 0
      DO 100 ISYM = 2, NSYM
         NBASF = NORB1(ISYM-1) + NORB2(ISYM-1)
         NNBASF = NBASF * (NBASF + 1) / 2
         ISB(ISYM) = ISB(ISYM-1) + NNBASF
  100 CONTINUE
      NBASF = NORB1(NSYM) + NORB2(NSYM)
      NNBASF = ISB(NSYM) + NBASF * (NBASF + 1) / 2
      LUMULB = -34
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ (LUMULB,'(4E30.20)') (SF(I), I = 1, NNBASF)
      CALL GPCLOSE (LUMULB,'KEEP')
      IF (.NOT. (R12NOB .OR. NORXR)) THEN
      DO 200 ISYM = 1, NSYM
         DO 210 I = NORB1(ISYM) + 2,  NVIR(ISYM)
            DO 220 J = NORB1(ISYM) + 1 , I - 1
               IJ = ISB(ISYM) + INDEX(I,J)
               IF (ABS(SF(IJ)) .GT. THRDIA) THEN
                  WRITE(LUPRI,'(/A/A,3I5,E20.10/)')
     *            '@ WARNING : Exchange matrix not diagonal', 
     *            '            Nondiagonal element is :',
     *                         ISYM,I,J,SF(IJ)
                  IF (ABS(SF(IJ)) .GT. SQRT(THRDIA))
     *               CALL QUIT('Exchange matrix not diagonal')
               END IF
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
      END IF
C  
C     Compute (ak|bl) = Sum_c X(ac)*(ck|bl) + Sum_d X(bd)*(ak|dl)
C     The sum over c runs over the auxiliary basis only
C
      DO 301 ISYMKA = 1, NSYM
         ISYMLB = ISYMKA
         DO 302 ISYMK = 1, NSYM
            ISYMA = MULD2H(ISYMK,ISYMKA)
            ISYMC = ISYMA
            ISYMKC = ISYMKA
            DO 303 ISYML = 1, NSYM
               ISYMB = MULD2H(ISYML,ISYMLB)
               ISYMD = ISYMB
               ISYMLD = ISYMLB
               DO 304 K = 1, NRHF(ISYMK)
                  DO 305 L = 1, NRHF(ISYML)
                     DO 306 A = 1, NVIR(ISYMA)
                        IF (R12CBS) THEN
			   NSTRC = 1
			   NENDC = NVIR(ISYMC)
			ELSE
                           IF (A. LE. NORB1(ISYMA)) THEN
                              NSTRC = NORB1(ISYMC) + 1
                              NENDC = NVIR(ISYMC)
                           ELSE
                              NSTRC = A                
                              NENDC = A                
                           END IF
                        END IF
                        NAK = IT1AM(ISYMA,ISYMK)
     *                      + NVIR(ISYMA)*(K-1) + A
                        DO 307 B = 1, NVIR(ISYMB)
                           IF (R12CBS) THEN
			      NSTRD = 1
			      NENDD = NVIR(ISYMD)
			   ELSE
                              IF (B. LE. NORB1(ISYMB)) THEN
                                 NSTRD = NORB1(ISYMD) + 1
                                 NENDD = NVIR(ISYMD)
                              ELSE
                                 NSTRD = B
                                 NENDD = B
                              END IF
                           END IF
                           NBL = IT1AM(ISYMB,ISYML)
     *                         + NVIR(ISYMB)*(L-1) + B
                           NAKBL = IT2AM(ISYMKA,ISYMLB)
     *                           + INDEX(NAK,NBL)
C
                           V2AM(NAKBL) = D0
C
                           DO 308 C = NSTRC, NENDC       
                              NCK = IT1AM(ISYMC,ISYMK)
     *                            + NVIR(ISYMC)*(K-1) + C 
                              NCKBL = IT2AM(ISYMKC,ISYMLB)
     *                              + INDEX(NCK,NBL) 
                              RR = R2AM(NCKBL)
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              V2AM(NAKBL) = V2AM(NAKBL) + XAC * RR
  308                      CONTINUE
                           DO 309 D = NSTRD, NENDD
                              NDL = IT1AM(ISYMD,ISYML)
     *                            + NVIR(ISYMD)*(L-1) + D
                              NAKDL = IT2AM(ISYMKA,ISYMLD)
     *                              + INDEX(NAK,NDL)
                              RR = R2AM(NAKDL)
                              XBD = SF(ISB(ISYMB)+INDEX(B,D))
                              V2AM(NAKBL) = V2AM(NAKBL) + XBD * RR
  309                      CONTINUE
C
  307                   CONTINUE
  306                CONTINUE
  305             CONTINUE
  304          CONTINUE
  303       CONTINUE
  302    CONTINUE
  301 CONTINUE
C   
      DO 999 IFLAG = 1, 2
C
      DO 110 I = 1, NOCDIM**4
         Q11(I,1,1,1) = D0
  110 CONTINUE
      DO 401 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         DO 402 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 403 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
               DO 404 I = 1, NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  DO 405 J = 1, NRHF(ISYMJ)
                     KOFFJ = IRHF(ISYMJ) + J
                     DO 406 A = 1, NVIR(ISYMA)                 
                        ISYMJA = MULD2H(ISYMJ,ISYMA)
                        DO 407 B = 1, NVIR(ISYMB)
C
                           AVIRT = A .GT. NRHFSA(ISYMA) .AND.
     *                             A .LE. NORB1(ISYMA)
                           BVIRT = B .GT. NRHFSA(ISYMB) .AND.
     *                             B .LE. NORB1(ISYMB)
                           
                           IF (R12CBS) THEN
                              SFAC = D1
			      IF (IFLAG .EQ. 1) THEN
                                 IF (A .LE. NORB1(ISYMA) .OR. 
     *                               B .LE. NORB1(ISYMB)) GOTO 407
                              ELSE IF (IFLAG .EQ. 2) THEN
                                 IF (A .LE. NRHFSA(ISYMA) .OR. 
     *                               B .LE. NRHFSA(ISYMB)) GOTO 407
                              ELSE IF (IFLAG .EQ. 3) THEN
                                 IF (A .LE. NRHFSA(ISYMA) .OR. 
     *                               B .LE. NRHFSA(ISYMB) .OR.
     *                              (AVIRT .AND. BVIRT)) GOTO 407
                              END IF
                           ELSE
                              IF (IFLAG .EQ. 2) THEN
                                 IF (AVIRT .OR. BVIRT) GOTO 407 
                              END IF
                              IF ((A .GT. NORB1(ISYMA) .AND.
     *                             B .GT. NORB1(ISYMB)) .OR.
     *                            (A .LE. NORB1(ISYMA) .AND.
     *                             B .LE. NORB1(ISYMB))) THEN
                                 SFAC = D1
                              ELSE
                                 SFAC = - D1
                              END IF
                              IF (IFLAG .EQ. 3) THEN
                                IF ((A .GT. NORB1(ISYMA) .AND.
     *                               B .GT. NORB1(ISYMB)) .OR.
     *                              (A .LE. NRHFSA(ISYMA) .AND.
     *                               B .LE. NRHFSA(ISYMB))) THEN
                                    SFAC = D1
                                ELSE
                                    SFAC = - D1
                                END IF
                                IF (AVIRT) THEN
                                   IF (.NOT. BVIRT) GOTO 407
                                ELSE IF (BVIRT) THEN
                                   IF (.NOT. AVIRT) GOTO 407
                                END IF
                              END IF
                           END IF
C           
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I-1) + A
                           NBJ = IT1AM(ISYMB,ISYMJ)
     *                         + NVIR(ISYMB)*(J-1) + B
                           NAIBJ = IT2AM(ISYMIA,ISYMJB)
     *                           + INDEX(NAI,NBJ)
                           DO 408 ISYMKA = 1, NSYM
                              ISYMLB = ISYMKA
                              ISYMK = MULD2H(ISYMA,ISYMKA)
                              ISYML = MULD2H(ISYMB,ISYMLB)
                              DO 409 K = 1, NRHF(ISYMK)
                                 KOFFK = IRHF(ISYMK) + K
                                 DO 410 L = 1, NRHF(ISYML)
                                    KOFFL = IRHF(ISYML) + L  
                                    NAK = IT1AM(ISYMA,ISYMK)
     *                                  + NVIR(ISYMA)*(K-1) + A
                                    NBL = IT1AM(ISYMB,ISYML)
     *                                  + NVIR(ISYMB)*(L-1) + B
                                    NAKBL = IT2AM(ISYMKA,ISYMLB)
     *                                    + INDEX(NAK,NBL)
                                    RR = R2AM(NAIBJ)
                                    VV = V2AM(NAKBL)
                                    Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                              SFAC * RR * VV 
  410                            CONTINUE
  409                         CONTINUE
  408                      CONTINUE
  407                   CONTINUE
  406                CONTINUE
  405             CONTINUE
  404          CONTINUE   
  403       CONTINUE
  402    CONTINUE
  401 CONTINUE   
C
      IJKL = 0
      FF = D1 / SQRT(D2)
      DO 600 K = 1, NOCDIM
         DO 601 L = 1, K
            DO 602 I = 1, NOCDIM
               DO 603 J = 1, I
                  IJKL = IJKL + 1
                  RR = Q11(I,J,K,L) + Q11(I,J,L,K)
                  VV = Q11(I,J,K,L) - Q11(I,J,L,K)
                  SF(IJKL) = RR
                  TF(IJKL) = VV
                  IF (I .EQ. J) SF(IJKL) = FF * SF(IJKL)
                  IF (K .EQ. L) SF(IJKL) = FF * SF(IJKL)
                  SF(IJKL) = SF(IJKL) * DP25 * DP5
                  TF(IJKL) = TF(IJKL) * DP75 * DP5
  603          CONTINUE
  602       CONTINUE
  601    CONTINUE
  600 CONTINUE
C
      CALL ERISFK(SF,NSPAIR,1)
      CALL ERISFK(TF,NSPAIR,1)
C
      IF (IPRINT .GT.  3) THEN
         GOTO (703,704,705), IFLAG
         CALL QUIT('INVALID IFLAG IN RXR')
  703    CALL AROUND('Singlet <RXR> matrix')
         CALL OUTPUT(SF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <RXR> matrix')
         CALL OUTPUT(TF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         GOTO 706
  704    CALL AROUND('Singlet <RXR@> matrix')
         CALL OUTPUT(SF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <RXR@> matrix')
         CALL OUTPUT(TF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         GOTO 706
  705    CALL AROUND('Singlet <RXR3> matrix')
         CALL OUTPUT(SF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <RXR3> matrix')
         CALL OUTPUT(TF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
  706    WRITE(LUPRI,*)
      END IF
C
      LUMUL = -34
      GOTO (701,702,707), IFLAG
      CALL QUIT('INVALID IFLAG IN RXR')
  701 CALL GPOPEN(LUMUL,'XMATSN','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') (SF(KL), KL = 1, NSPAIR*NSPAIR)
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'XMATTN','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') (TF(KL), KL = 1, NSPAIR*NSPAIR)
      CALL GPCLOSE(LUMUL,'KEEP')  
      GOTO 999    
  702 CALL GPOPEN(LUMUL,'XMATSP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') (SF(KL), KL = 1, NSPAIR*NSPAIR)
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN (LUMUL,'XMATTP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') (TF(KL), KL = 1, NSPAIR*NSPAIR)
      CALL GPCLOSE(LUMUL,'KEEP')  
cWK   GOTO 999
  707 CALL GPOPEN(LUMUL,'XMATSQ','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') (SF(KL), KL = 1, NSPAIR*NSPAIR)
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'XMATTQ','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') (TF(KL), KL = 1, NSPAIR*NSPAIR)
      CALL GPCLOSE(LUMUL,'KEEP')
  999 CONTINUE
      RETURN
      END 
C  /* Deck makegx */   
      SUBROUTINE MAKEGX(V2AM,GX,F2AM,WORK,LWORK,IPBAS)
C
C     This subroutine computes the X(p'k) matrix (exchange operator)
C     by summing two-electron integrals (ip'|ik) in the orthonormal
C     basis. The orbitals i and k belong to the primary basis, 
C     whereas p' belongs to the secondary basis. 
C
C     On input, V2AM contains the (ip'|ik) integrals and F2AM
C     contains the (Ip'|Ik) integrals where the orbitals I
C     belong to the frozen core. The X(p'k) exchange matrix
C     is returned as GX.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
C     Modified for natural virtual orbitals NATVIR (WK/UniKA/22-11-2005).
C
#include <implicit.h>
#include <priunit.h>    
      DIMENSION V2AM(*), GX(*), F2AM(*), WORK(*)
      integer   iv1fro(8,8),it1vm(8,8),it2vm(8,8),iv2fro(8,8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfro.h>
#include <r12int.h> 
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      IF (NATVIR) THEN
         CALL DZERO(GX,NT1AMX)
         NNBASF = 0
         DO ISYM = 1, NSYM
            NONV = NVIR(ISYM) + NRHFS(ISYM)
            NNBASF = NNBASF + NONV * (NONV + 1) / 2 
         END DO
         LUNT = -1
         IF (NNBASF .GT. LWORK) CALL QUIT('Not enough memory in MAKEGX')
         CALL GPOPEN(LUNT,'AUXFCKN','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
         READ(LUNT,'(4E30.20)') (WORK(I), I = 1, NNBASF)
         CALL GPCLOSE(LUNT,'KEEP')
         IOFF = 0
         DO ISYMP = 1, NSYM
            ISYMK =  ISYMP
            NPK = IT1AM(ISYMP,ISYMK)
            DO K =1 + NRHFFR(ISYMK), NRHFS(ISYMK)
               DO P = 1 + NRHFS(ISYMP), NVIR(ISYMP) + NRHFS(ISYMP)
                  NPK = NPK + 1
                  GX(NPK) = WORK(IOFF + P*(P-1)/2 + K)
               END DO
            END DO
            NONV = NVIR(ISYMP) + NRHFS(ISYMP)
            IOFF = IOFF + NONV * (NONV + 1) / 2
         END DO
c        LUNT = -1
c        CALL GPOPEN(LUNT,'AUXFCK','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
c        READ(LUNT,'(4E30.20)') (WORK(I), I = 1, NNBASF)
c        CALL GPCLOSE(LUNT,'KEEP')
c        IOFF = 0
c        DO ISYMP = 1, NSYM
c           ISYMK =  ISYMP
c           NPK = IT1AM(ISYMP,ISYMK)
c           DO K = NRHFFR(ISYMK)+1, NRHFSA(ISYMK)
c              DO P = 1, NVIR(ISYMP)
c                 NPK = NPK + 1
c                 GX(NPK) = WORK(IOFF + P*(P-1)/2 + K)
c              END DO
c           END DO
c           IOFF = IOFF + NVIR(ISYMP) * (NVIR(ISYMP) + 1) / 2
c        END DO
         GOTO 100
      END IF
 
      IF (R12PRP) THEN
         CALL DZERO(GX,NT1VMX)
         DO ISYMK = 1,NSYM
            NV1FRO(ISYMK) = 0
            DO ISYMJ = 1,NSYM
               ISYMI  = MULD2H(ISYMK,ISYMJ)
               NV1FRO(ISYMK) = NV1FRO(ISYMK) + NRHFFR(ISYMJ)*
     &                          (NORB1(ISYMI)-NRHFFR(ISYMI))
            ENDDO
         ENDDO
         
         DO ISYMK = 1, NSYM
            ICOUN1 = 0
            ICOUN3 = 0 
            ICOUN4 = 0 
            DO ISYMJ = 1,NSYM
               ISYMI  = MULD2h(ISYMK,ISYMJ)
               IT1VM(ISYMI,ISYMJ)  = ICOUN1
               ICOUN1 = ICOUN1 + (NORB1(ISYMJ)-NRHFFR(ISYMJ))
     &                     *NVIR(ISYMI)
               IV2FRO(ISYMI,ISYMJ) = ICOUN3
               IV1FRO(ISYMI,ISYMJ) = ICOUN4
               ICOUN3 = ICOUN3 + NT1FRO(ISYMI)*NV1FRO(ISYMJ)
               ICOUN4 = ICOUN4 + (NRHFFR(ISYMI))
     &                   *(NORB1(ISYMJ)-NRHFFR(ISYMJ))
            ENDDO
         ENDDO  
      ELSE
         CALL DZERO(GX,NT1AMX)
      ENDIF
C
C     Compute SUM_i ( i p' | i k )
C
celena
      IF (R12PRP) THEN
         DO ISYMIP = 1, NSYM
            ISYMIK = ISYMIP
            DO ISYMI = 1, NSYM
               ISYMP = MULD2H(ISYMI,ISYMIP)
               ISYMK = ISYMP
               DO I = 1, NRHF(ISYMI)
                  NPK = IT1vM(ISYMP,ISYMK)
                  DO K = NRHFFR(ISYMK)+1, NORB1(ISYMK)
                     DO P = 1, NVIR(ISYMP)
                        IF (.NOT. ONEAUX) THEN
                           NPK = NPK + 1
                           NPI = IT1AM(ISYMP,ISYMI)
     *                         + NVIR(ISYMP)*(I-1) + P
                           NKI = IT1AM(ISYMK,ISYMI)
     *                         + NVIR(ISYMK)*(I-1) + K
                           NPIKI = IT2AM(ISYMIP,ISYMIK)
     *                           + INDEX(NPI,NKI)
                           GX(NPK) = GX(NPK) + V2AM(NPIKI)
                        ENDIF
                     ENDDO   
                  ENDDO    
               ENDDO   
            ENDDO    
         ENDDO    
ccelena
      ELSE
         DO 101 ISYMIP = 1, NSYM
            ISYMIK = ISYMIP
            IKOFF = NH1AM(ISYMIK) * (NH1AM(ISYMIK) + 1) / 2
            DO 102 ISYMI = 1, NSYM
               ISYMP = MULD2H(ISYMI,ISYMIP)
               ISYMK = ISYMP  
               DO 201 I = 1, NRHFA(ISYMI)
                  NPK = IT1AM(ISYMP,ISYMK)
                  DO 202 K = NRHFFR(ISYMK)+1, NRHFS(ISYMK)
                     DO 203 P = 1, NVIR(ISYMP)
                        NPK = NPK + 1
                        IF (ONEAUX) THEN
                           IF (P .LE. NORB1(ISYMP)) THEN
                              NPI = IH1AM(ISYMP,ISYMI)
     *                            + NORB1(ISYMP)*(I-1) + P
                              NKI = IH1AM(ISYMK,ISYMI)
     *                            + NORB1(ISYMK)*(I-1) + K
                              NPIKI = IH2AM(ISYMIP,ISYMIK)
     *                              + INDEX(NPI,NKI)  
                           ELSE
                              NPI = IG1AM(ISYMP,ISYMI)
     *                           + NORB2(ISYMK)*(I-1) + P - NORB1(ISYMP)
                              NKI = IH1AM(ISYMK,ISYMI)
     *                            + NORB1(ISYMK)*(I-1) + K
                              NPIKI = IKOFF + IH2AM(ISYMIP,ISYMIK)
     *                               + NH1AM(ISYMIK) * (NPI - 1) + NKI
                           END IF
                        ELSE
                           NPI = IT1AM(ISYMP,ISYMI)
     *                         + NVIR(ISYMP)*(I-1) + P
                           NKI = IH1AM(ISYMK,ISYMI)
     *                         + NVIR(ISYMK)*(I-1) + K
                           NPIKI = IT2AM(ISYMIP,ISYMIK)
     *                           + INDEX(NPI,NKI)  
                        END IF
                        GX(NPK) = GX(NPK) + V2AM(NPIKI)
  203                CONTINUE
  202             CONTINUE
  201          CONTINUE 
  102       CONTINUE
  101    CONTINUE  
      ENDIF
C
      IF (FROIMP) THEN
C
C     Add contribution from frozen orbitals
C
         IF (R12PRP) THEN
           DO ISYMIP = 1, NSYM
              ISYMIK = ISYMIP
              DO ISYMI = 1, NSYM
                 ISYMP = MULD2H(ISYMI,ISYMIP)
                 ISYMK = ISYMP
                 DO I = 1, NRHFFR(ISYMI)
                    NPK = IT1VM(ISYMP,ISYMK)
                    DO K = 1, NORB1(ISYMK)-NRHFFR(ISYMK)
                       DO P = 1, NVIR(ISYMP)
                          NPK = NPK + 1
                          NKI = IV1FRO(ISYMI,ISYMK)
     *                        + NRHFFR(ISYMI)*(K-1) + I
                          NPI = IT1FRO(ISYMP,ISYMI)
     *                        + NVIR(ISYMP)*(I-1) + P
                          NPIKI = IV2FRO(ISYMIP,ISYMIK)
     *                          + NT1FRO(ISYMIP)*(NKI-1)
     *                          + NPI
                          GX(NPK) = GX(NPK) + F2AM(NPIKI)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
         ELSE
           DO 301 ISYMIP = 1, NSYM
              ISYMIK = ISYMIP
              DO 302 ISYMI = 1, NSYM
                 ISYMP = MULD2H(ISYMI,ISYMIP)
                 ISYMK = ISYMP 
                 DO 401 I = 1, NRHFFR(ISYMI)
                    NPK = IT1AM(ISYMP,ISYMK)  
                    DO 402 K = 1, NRHF(ISYMK)
                       DO 403 P = 1, NVIR(ISYMP)
                          NPK = NPK + 1
                          NKI = IF1FRO(ISYMI,ISYMK)
     *                        + NRHFFR(ISYMI)*(K-1) + I 
                          NPI = IT1FRO(ISYMP,ISYMI)
     *                        + NVIR(ISYMP)*(I-1) + P
                          NPIKI = IF2FRO(ISYMIP,ISYMIK)
     *                          + NT1FRO(ISYMIP)*(NKI-1)   
     *                          + NPI
                          GX(NPK) = GX(NPK) + F2AM(NPIKI)
  403                  CONTINUE
  402              CONTINUE
  401           CONTINUE 
  302        CONTINUE
  301       CONTINUE  
         ENDIF
      END IF   
C
  100 CONTINUE
C
      IF (R12OLD .OR. R12CBS) GOTO 700
C
C     Zero primary block 
C
      IF (.NOT. R12HYB .OR. (R12HYB .AND. IPBAS .EQ. 2)) THEN
         DO 500 ISYM = 1, NSYM
            IF (R12PRP) THEN 
               NPK = IT1vM(ISYM,ISYM) + 1
               DO K = 1, NORB1(ISYM)-NRHFFR(ISYM)
                  CALL DZERO(GX(NPK),NORB1(ISYM))
                  NPK = NPK + NVIR(ISYM)
               ENDDO   
            ELSE
               NPK = IT1AM(ISYM,ISYM) + 1
               DO 501 K = 1, NRHF(ISYM)
                  CALL DZERO(GX(NPK),NORB1(ISYM))
                  NPK = NPK + NVIR(ISYM)
  501          CONTINUE
            ENDIF 
  500    CONTINUE
      END IF
C
  700 CONTINUE
C
C     Zero secondary block 
C
      IF (R12HYB .AND. IPBAS .EQ. 1) THEN
         DO 600 ISYM = 1, NSYM
            IF (R12PRP) THEN
               NPK = IT1vM(ISYM,ISYM) + NORB1(ISYM)+ 1
               DO K = 1, NORB1(ISYM)-NRHFFR(ISYM)
                  CALL DZERO(GX(NPK),NORB2(ISYM))
                  NPK = NPK + NVIR(ISYM)
               ENDDO   
            ELSE 
               NPK = IT1AM(ISYM,ISYM) + NORB1(ISYM) + 1
               DO 601 K = 1, NRHF(ISYM)
                  CALL DZERO(GX(NPK),NORB2(ISYM))
                  NPK = NPK + NVIR(ISYM)
  601          CONTINUE
            ENDIF 
  600    CONTINUE
      END IF 
C
C     Print section
C
      IF (IPRINT .GT. 3) THEN
         DO 300 ISYM = 1, NSYM
            IF (NRHF(ISYM)*NVIR(ISYM) .NE. 0) 
     *      WRITE(LUPRI,'(/A,I2)') ' X-MATRIX/ISYM =',ISYM
            IF (R12PRP) THEN
               NPK = IT1vM(ISYM,ISYM) + 1
               CALL OUTPUT(GX(NPK),1,NVIR(ISYM),1,NORB1(ISYM),
     *                     NVIR(ISYM),NORB1(ISYM),1,LUPRI) 
            ELSE 
               NPK = IT1AM(ISYM,ISYM) + 1
               CALL OUTPUT(GX(NPK),1,NVIR(ISYM),1,NRHF(ISYM),
     *                     NVIR(ISYM),NRHF(ISYM),1,LUPRI) 
            ENDIF
  300    CONTINUE
      END IF
      RETURN
      END
C  /* Deck r12drv */ 
      SUBROUTINE R12DRV(V2AM,R2AM,U2AM,S2AM,T2AM,
     *                  EV,RXRS,RXRT,WRK,LWRK,
     *                  QQ2,QQ4,QQ6)
C
C     Driver routine for the MP2-R12 calculation.
C
C              V2AM = ( ip | 1/r12 | jq)
C 
C              R2AM = ( ip | r12 | jq)
C 
C              U2AM = ( ip | [T1+T2,r12] | jq)
C
C              S2AM = ( ip | [r12,K1+K2] | jq) 
C
C              T2AM = ( ip | (K1+K2) r12 | jq)
C
C              QQ2  = ( ik | (1/r12)*f(12) | jl )
C
C              QQ4  = ( ik | f(12)**2 | jl )
C
C              QQ6  = ( ik | [[f12,T1],f12] | jl )
C 
C     p and q denote orbitals of both the primary and
C     secondary basis. Integrals where both p and q belong
C     to the secondary basis are not needed. The numerical
C     values at the corresponding places in the arrays are, 
C     meaningless!             
C
C     EV is an array with orbital energies and RXRS and RXRT
C     contain the matrices computed in the subroutine RXR.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>  
#include <maxorb.h>
#include <infinp.h>
      DIMENSION R2AM(*), V2AM(*), U2AM(*), EV(*), WRK(*),
     *          RXRS(*), RXRT(*), S2AM(*), T2AM(*),
     *          QQ2(*), QQ4(*), QQ6(*),ISB(8)
      LOGICAL LGLOCAL,FRSWRT
C
      LGLOCAL = BOYORB.OR.PIPORB
C
      NOCDIM = NRHFT 
      NICDIM = NRHFTS
      NACDIM = 0
      DO ISYM = 1, NSYM
         NACDIM = NACDIM + NRHFA(ISYM)
      END DO
      NAPAIR = NACDIM * (NACDIM + 1) / 2
      NSPAIR = NOCDIM * (NOCDIM + 1) / 2
      NTPAIR = NOCDIM * (NOCDIM - 1) / 2
      NPAIR2 = NSPAIR ** 2
      NODIM4 = 2 * NOCDIM ** 4
      NIDIM4 = NICDIM ** 4
      IF (LGLOCAL) THEN
         NCT = 0
         NRT = 0
         DO ISYM = 1, NSYM
            NCT = NCT + NRHFFR(ISYM)
            NRT = NRT + NORB1(ISYM)
         END DO
         NBASM = (NRT-NCT)*(NRT-NCT+1)/2
      ELSE
         NBASM = NOCDIM*(NOCDIM+1)/2
      END IF
C
      IES = 1
      IET = IES + NSPAIR
      IFS = IET + NSPAIR
      IFT = IFS + NSPAIR
      IEVS = IFT + NSPAIR 
      ININV = IEVS + NSPAIR
      IQQ11 = ININV + NSPAIR*10
C
      IV11 = IQQ11 + NIDIM4
      IU11 = IV11 + NODIM4
      IB11 = IU11 + NODIM4
      IW11 = IB11 + NODIM4 * NSPAIR
      IQ11 = IW11 + NODIM4
      IR11 = IQ11 + NODIM4
C
      IVS11 = IR11 + NODIM4
      IUS11 = IVS11 + NPAIR2
      IBS11 = IUS11 + NPAIR2
      IWS11 = IBS11 + NPAIR2 * NSPAIR
      IQS11 = IWS11 + NPAIR2
      IRS11 = IQS11 + NPAIR2
C
      IVT11 = IRS11 + NPAIR2
      IUT11 = IVT11 + NPAIR2
      IBT11 = IUT11 + NPAIR2
      IWT11 = IBT11 + NPAIR2 * NSPAIR
      IQT11 = IWT11 + NPAIR2
      IRT11 = IQT11 + NPAIR2
C
      IFLMAT= IRT11 + NPAIR2
      ININ  = IFLMAT+ NBASM  
      III   = ININ  + NOCDIM ** 2
      IJJ   = III   + NSPAIR
      ICS11 = IJJ   + NSPAIR
      ICT11 = ICS11 + NPAIR2
      IOLWS = ICT11 + NPAIR2
      IOLWT = IOLWS + NPAIR2
      IOLYS = IOLWT + NPAIR2
      IOLYT = IOLYS + NPAIR2
      IOLZS = IOLYT + NPAIR2
      IOLZT = IOLZS + NPAIR2
      IEVL  = IOLZT + NPAIR2
C
      IENDW = IEVL  + NSPAIR
C
c    Elena  work space allocation for r12grad
      IF (R12PRP) THEN
         ISB(1) = 0
         DO ISYM = 2, NSYM
            NBASF = NORB1(ISYM-1) + NORB2(ISYM-1)
            NNBASF = NBASF * (NBASF + 1) / 2
            ISB(ISYM) = ISB(ISYM-1) + NNBASF
         ENDDO
         NBASF = NORB1(NSYM) + NORB2(NSYM)
         NNBASF = ISB(NSYM) + NBASF * (NBASF + 1) / 2

         IGRAD = IENDW + NNBASF
         LWRK1 = LWRK  - IGRAD 
         IF (IGRAD .GT. LWRK) THEN
            WRITE(LUPRI,'(A,I20/A,I20)') 
     *      ' WORK SPACE REQUIRED =  ',IGRAD,
     *      ' WORK SPACE AVAILABLE = ',LWRK
            CALL QUIT('INSUFFICIENT WORK SPACE IN R12DRV')
         ENDIF
c    Elena
C   
      ELSE
         IF (IENDW .GT. LWRK) THEN
             WRITE(LUPRI,'(A,I20/A,I20)')
     *       ' WORK SPACE REQUIRED =  ',IENDW,
     *       ' WORK SPACE AVAILABLE = ',LWRK
             CALL QUIT('INSUFFICIENT WORK SPACE IN R12DRV')
         ENDIF
      ENDIF

      IF (R12ECO) THEN
         CALL ECODRV(V2AM,R2AM,EV,WRK(IV11),WRK(IB11),WRK(IQ11),
     *               WRK(IVS11),WRK(IVT11),WRK(IWS11),WRK(IWT11),
     *               WRK(IQS11),WRK(IQT11),
     *               WRK(IBS11),WRK(IBT11),WRK(IRS11),
     *               WRK(IRT11),NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *               WRK(IES),WRK(IET),WRK(IFS),WRK(IFT),
     *               WRK(IEVS),WRK(ININV),WRK(IENDW),WRK(ICS11),
     *               WRK(ICT11))
          RETURN
      ENDIF 
      IF (.NOT. NOTONE .OR. R12OLD .OR. LGLOCAL) THEN 
         CALL AROUND('MP2-R12 ansatz := (1 - P1 - P2 + P1*P2) * R12')
         CALL MAKEVR(V2AM,R2AM,U2AM,EV,WRK(IQQ11),
     *        WRK(IV11),WRK(IU11),WRK(IB11),
     *        WRK(IW11),WRK(IQ11),WRK(IR11),
     *        WRK(IVS11),WRK(IUS11),WRK(IBS11),
     *        WRK(IWS11),WRK(IQS11),WRK(IRS11),
     *        WRK(IVT11),WRK(IUT11),WRK(IBT11),
     *        WRK(IWT11),WRK(IQT11),WRK(IRT11),
     *        NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *        WRK(IES),WRK(IET),WRK(IFS),WRK(IFT),
     *        WRK(IEVS),WRK(ININV),RXRS,RXRT,QQ2,QQ4,QQ6, 
     *        R2AM,WRK(IFLMAT),WRK(ININ),WRK(III),WRK(IJJ),
     *        WRK(ICS11),WRK(ICT11),WRK(IOLWS),WRK(IOLWT),
     *        WRK(IOLYS),WRK(IOLYT),WRK(IOLZS),WRK(IOLZT),
     *        WRK(IEVL),LGLOCAL,WRK(IGRAD),WRK(IENDW),LWRK-IENDW,
     *        NACDIM,NAPAIR)
      END IF
      IF (R12OLD .OR. LGLOCAL .OR. NOTTWO) GOTO 999
      IANSAZ = 2
      FRSWRT = .TRUE.
      CALL AROUND('MP2-R12 ansatz := (1 - O1 - O2 + O1*O2) * R12')
      CALL MAKERV(V2AM,R2AM,U2AM,S2AM,T2AM,EV,WRK(IQQ11),
     *     WRK(IV11),WRK(IU11),WRK(IB11),
     *     WRK(IW11),WRK(IQ11),WRK(IR11),
     *     WRK(IVS11),WRK(IUS11),WRK(IBS11),
     *     WRK(IWS11),WRK(IQS11),WRK(IRS11),
     *     WRK(IVT11),WRK(IUT11),WRK(IBT11),
     *     WRK(IWT11),WRK(IQT11),WRK(IRT11),
     *     NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *     WRK(IES),WRK(IET),WRK(IFS),WRK(IFT),
     *     WRK(IEVS),WRK(ININV),QQ2,QQ4,QQ6,WRK(ICS11),WRK(ICT11),
     *     IANSAZ,FRSWRT,WRK(IENDW),LWRK-IENDW,
     *     NACDIM,NAPAIR)
 999  CONTINUE
      IF (NOTTRE .OR. R12OLD .OR. LGLOCAL) RETURN
      IANSAZ = 3
      FRSWRT = .TRUE.
      CALL AROUND('MP2-R12 ansatz := '//
     *            '(1 - V1*V2) * (1 - O1 - O2 + O1*O2) * R12')
      CALL MAKERV(V2AM,R2AM,U2AM,S2AM,T2AM,EV,WRK(IQQ11),
     *     WRK(IV11),WRK(IU11),WRK(IB11),
     *     WRK(IW11),WRK(IQ11),WRK(IR11),
     *     WRK(IVS11),WRK(IUS11),WRK(IBS11),
     *     WRK(IWS11),WRK(IQS11),WRK(IRS11),
     *     WRK(IVT11),WRK(IUT11),WRK(IBT11),
     *     WRK(IWT11),WRK(IQT11),WRK(IRT11),
     *     NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *     WRK(IES),WRK(IET),WRK(IFS),WRK(IFT),
     *     WRK(IEVS),WRK(ININV),QQ2,QQ4,QQ6,WRK(ICS11),WRK(ICT11),
     *     IANSAZ,FRSWRT,WRK(IENDW),LWRK-IENDW,
     *     NACDIM,NAPAIR)
      RETURN
      END
C  /* Deck makevr */
      SUBROUTINE MAKEVR(V2AM,R2AM,U2AM,EV,QQ11,V11,U11,B11,W11,Q11,
     *                  R11,
     *                  VS11,US11,BS11,WS11,QS11,RS11,VT11,UT11,BT11,
     *                  WT11,QT11,RT11,NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *                  ES,ET,FS,FT,EVS,CNINV,RXRS,RXRT,QQ2,QQ4,QQ6,
     *                  TIJAB,FLMAT,NIND,III,JJJ,CS11,CT11,OLWS,OLWT,
     *                  OLYS,OLYT,OLZS,OLZT,EVL,LGLOCAL,GRAD,
     *                  WORK,LWORK,NACDIM,NAPAIR)
C
C     This subroutine computes the V(klmn), X(klmn), and B(klmn)
C     matrices and evaluates the MP2-R12 energies.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C     Modified by Claire C.M. Samson (University of Karlsruhe, 28 April 2003)
C                                                          for local orbitals.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0,
     *           THRES = 1.0D-10)
      DIMENSION R2AM(*), V2AM(*), U2AM(*), EV(*)
      DIMENSION TIJAB(*), FLMAT(*), WORK(*), GRAD(*)
      DIMENSION NIND(NOCDIM,NOCDIM), III(NSPAIR), JJJ(NSPAIR)
      REAL*8  R1, R2, V1, U1, U2, VR, UR, RR,
     *        V11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        U11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        B11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        W11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        Q11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        R11(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
      DIMENSION QQ2(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *          QQ4(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *          QQ6(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
      DIMENSION VS11(NSPAIR,NSPAIR)
      DIMENSION US11(NSPAIR,NSPAIR)
      DIMENSION BS11(NSPAIR,NSPAIR)
      DIMENSION WS11(NSPAIR,NSPAIR)
      DIMENSION QS11(NSPAIR,NSPAIR)
      DIMENSION RS11(NSPAIR,NSPAIR)
      DIMENSION VT11(NSPAIR,NSPAIR),WT11(NSPAIR,NSPAIR)
      DIMENSION UT11(NSPAIR,NSPAIR),BT11(NSPAIR,NSPAIR)
      DIMENSION QT11(NSPAIR,NSPAIR),RT11(NSPAIR,NSPAIR)
      DIMENSION RXRS(NSPAIR,NSPAIR),RXRT(NSPAIR,NSPAIR)
      DIMENSION ES(NSPAIR),ET(NSPAIR),FS(NSPAIR),FT(NSPAIR)
      DIMENSION EVS(NSPAIR),CNINV(NSPAIR,10)
      DIMENSION QQ11(NICDIM,NICDIM,NICDIM,NICDIM)
      DIMENSION CS11(NSPAIR,NSPAIR),CT11(NSPAIR,NSPAIR)
      DIMENSION OLWS(NSPAIR,NSPAIR),OLWT(NSPAIR,NSPAIR)
      DIMENSION OLYS(NSPAIR,NSPAIR),OLYT(NSPAIR,NSPAIR)
      DIMENSION OLZS(NSPAIR,NSPAIR),OLZT(NSPAIR,NSPAIR)
      DIMENSION EVL(NSPAIR)
      INTEGER IDUM,LWORK
      LOGICAL LDUM
      LOGICAL LGLOCAL, LOCRST
      CHARACTER*11 CFN
      CHARACTER*3 APPROX(0:7,2), IPCC
      INTEGER NKILJ(8)
#include <r12int.h>
#include <ccr12int.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      DATA APPROX
     &/"0  ","A  ","A  ","A' ","B  ","B  ","xxx","xxx", 
     & "0  ","A  ","A  ","A' ","B' ","B' ","B  ","B  "/

      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
 1001 FORMAT('CSMAT.',I2.2,'.',I2.2)
 1002 FORMAT('CTMAT.',I2.2,'.',I2.2)
Ckr   We need to have the next two unit numbers in consecutive order, and
Ckr   thus assign a very large unit number in order to avoid conflicts with 
Ckr   already open records.\kr-231007\
Ckr
      LUMULO = 90
      LUMULN = 91
      IF (R12CBS) THEN
        FACX = - D1
      ELSE
        FACX = D1
      END IF
      DELTA = R12LEV
      FF = D1 / SQRT(D2)
      LU33 = -33
      CALL GPOPEN(LU33,'AUXQ12','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ(LU33,'(4E30.20)') QQ11
      CALL GPCLOSE(LU33,'KEEP')
      DO NPIJ = 1, NSPAIR
         IJ = 0
         DO I = 1, NOCDIM
            DO J = 1, I
               IJ = IJ + 1
               IF (IJ.EQ.NPIJ) THEN
                   III(NPIJ) = I
                   JJJ(NPIJ) = J
               END IF
            END DO
         END DO
      END DO
      DO I = 1, NOCDIM
         DO J = 1, NOCDIM
            IF (I.GT.J) THEN
                NIND(I,J) = I * (I-1) /2 + J
            ELSE
                NIND(I,J) = J * (J-1) /2 + I
            END IF
         END DO
      END DO
      IJ = 0
      DO I = 1, NOCDIM
         DO J = 1, I
            IJ = IJ + 1
            EVS(IJ) = EV(I) + EV(J)
         END DO
      END DO
      IF (LGLOCAL) THEN
         NCT = 0
         NRT = 0
         DO ISYM = 1, NSYM
            NCT = NCT + NRHFFR(ISYM)
            NRT = NRT + NORB1(ISYM)
         END DO
         LU99 = -99
         CALL GPOPEN(LU99,'FLOCA','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
         NBASM = (NRT-NCT)*(NRT-NCT+1)/2
         DO I=1,NBASM 
            READ(LU99,'(D30.20)') FLMAT(I)
         END DO
         CALL GPCLOSE(LU99,'KEEP')
         DO IJ  = 1, NSPAIR
            NI  = III(IJ)
            NJ  = JJJ(IJ)
            NFI = NIND(NI,NI)
            NFJ = NIND(NJ,NJ) 
            EVL(IJ) = FLMAT(NFI) + FLMAT(NFJ)
         END DO
      ELSE
         NFMAT = NOCDIM*(NOCDIM+1)/2
         CALL DZERO(FLMAT,NFMAT)
         DO IJ  = 1, NSPAIR
            NI  = III(IJ)
            NJ  = JJJ(IJ)
            NFI = NIND(NI,NI)
            EVL(IJ) = EVS(IJ)
            IF (NI.EQ.NJ) FLMAT(NFI) = EVL(IJ) * DP5
         END DO
      END IF
      IF (IPRINT .GE. 4 .OR. LGLOCAL) THEN
         CALL AROUND('Fock matrix in occupied space') 
         WRITE(LUPRI,*)
         CALL OUTPAK(FLMAT,NOCDIM,1,LUPRI)
         CALL AROUND('Pairs of orbital energies') 
         CALL OUTPUT(EVL,1,NSPAIR,1,1,NSPAIR,1,1,LUPRI)
      END IF
C
      DO 205 I = 1, NOCDIM**4
         R11(I,1,1,1) = D0
         Q11(I,1,1,1) = D0
         V11(I,1,1,1) = D0
         U11(I,1,1,1) = D0
         B11(I,1,1,1) = D0
         W11(I,1,1,1) = D0
  205 CONTINUE
      CALL DZERO(ES,NSPAIR)
      CALL DZERO(ET,NSPAIR)
C
C     PRINT ORBITAL ENERGIES
C
      IF (IPRINT .LE. 5) GOTO 993
      DO 51 ISYM = 1, NSYM 
         WRITE(LUPRI,'(/A,I2)') ' OCCUPIED ORBITALS OF SYMMETRY',ISYM
         DO 52 I = 1, NRHF(ISYM)
            KOFFI = IRHF(ISYM) + I 
            WRITE(LUPRI,'(2I5,G20.10)') I,KOFFI,EV(KOFFI)
  52     CONTINUE
         WRITE(LUPRI,'(/A,I2)') ' PRIMARY ORBITALS OF SYMMETRY',ISYM
         DO 53 I = 1, NORB1(ISYM) 
            KOFFI = IVIR(ISYM) + I 
            WRITE(LUPRI,'(2I5,G20.10)') I,KOFFI,EV(KOFFI)
  53     CONTINUE
         WRITE(LUPRI,'(/A,I2)') ' SECONDARY ORBITALS OF SYMMETRY',ISYM
         DO 54 I = NORB1(ISYM) + 1, NVIR(ISYM)
            KOFFI = IVIR(ISYM) + I 
            WRITE(LUPRI,'(2I5,G20.10)') I,KOFFI,EV(KOFFI)
  54     CONTINUE
  51  CONTINUE
  993 CONTINUE
C
C     COMPUTE ITERATIVELY THE AMPLITUDES
C
      IF (LGLOCAL) THEN
         ICOUNT=1
         CALL DZERO(TIJAB,NH2AMX)
  900    ICONV=0
         E2 = D0
         DO 701 ISYMIA = 1, NSYM
            ISYMJB = ISYMIA
            DO 702 ISYMI = 1, NSYM
               ISYMA = MULD2H(ISYMI,ISYMIA)
               DO 703 ISYMJ = 1, NSYM
                  ISYMB = MULD2H(ISYMJ,ISYMJB)
                  DO 801 I = 1, NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     DO 802 J = 1, NRHF(ISYMJ)
                        KOFFJ = IRHF(ISYMJ) + J
                        DO 803 A = NRHFS(ISYMA) + 1, NORB1(ISYMA)
                           ISYMJA = MULD2H(ISYMJ,ISYMA)
                           KOFFA = IVIR(ISYMA) + A
                           DO 804 B = NRHFS(ISYMB) + 1, NORB1(ISYMB)
                              ISYMIB = MULD2H(ISYMI,ISYMB)
                              KOFFB = IVIR(ISYMB) + B
                              IF (ONEAUX) THEN
                              NAI = IH1AM(ISYMA,ISYMI) +
     *                              NORB1(ISYMA)*(I-1) + A
                              NBJ = IH1AM(ISYMB,ISYMJ) +
     *                              NORB1(ISYMB)*(J-1) + B
                              NAIBJ = IH2AM(ISYMIA,ISYMJB) +
     *                                INDEX(NAI,NBJ)
                              NAJ = IH1AM(ISYMA,ISYMJ) +
     *                              NORB1(ISYMA)*(J-1) + A
                              NBI = IH1AM(ISYMB,ISYMI) +
     *                              NORB1(ISYMB)*(I-1) + B
                              NAJBI = IH2AM(ISYMJA,ISYMIB) +
     *                                INDEX(NAJ,NBI)
                              ELSE
                              NAI = IT1AM(ISYMA,ISYMI) +
     *                              NVIR(ISYMA)*(I-1) + A
                              NBJ = IT1AM(ISYMB,ISYMJ) +
     *                              NVIR(ISYMB)*(J-1) + B
                              NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                                INDEX(NAI,NBJ)
                              NAJ = IT1AM(ISYMA,ISYMJ) +
     *                              NVIR(ISYMA)*(J-1) + A
                              NBI = IT1AM(ISYMB,ISYMI) +
     *                              NVIR(ISYMB)*(I-1) + B
                              NAJBI = IT2AM(ISYMJA,ISYMIB) +
     *                                INDEX(NAJ,NBI)
                              END IF
                              VAIBJ = V2AM(NAIBJ)
                              VAJBI = V2AM(NAJBI)
                              VV = VAIBJ
C
                              TKI = D0
                              TKJ = D0
                              TKA = D0
                              TKB = D0
                              TINI = TIJAB(NAIBJ)
C
C   LOOP FOR COEFFICIENTS (KJAB)
C
                              ISYMKA = ISYMJB
                              ISYMK  = MULD2H(ISYMA,ISYMKA)
                              DO 902 K = 1, NRHF(ISYMK)
                                 IF ((K.EQ.I).AND.(ISYMI.EQ.ISYMK))
     *                                                     GOTO 902
                                 IF (ONEAUX) THEN
                                   NAK = IH1AM(ISYMA,ISYMK) +
     *                                   NORB1(ISYMA)*(K-1) + A
                                   NAKBJ = IH2AM(ISYMKA,ISYMJB) +
     *                                   INDEX(NAK,NBJ)
                                 ELSE
                                   NAK = IT1AM(ISYMA,ISYMK) +
     *                                   NVIR(ISYMA)*(K-1) + A
                                   NAKBJ = IT2AM(ISYMKA,ISYMJB) +
     *                                   INDEX(NAK,NBJ)
                                 END IF
                                 NIK = INDEX(I,K)
                                 XXX = FLMAT(NIK) * TIJAB(NAKBJ)
                                 TKI = TKI + XXX
  902                         CONTINUE
C
C   LOOP FOR COEFFICIENTS (IKAB)
C
                              ISYMKB = ISYMIA
                              ISYMK  = MULD2H(ISYMB,ISYMKB)
                              DO 904 K = 1, NRHF(ISYMK)
                                 IF ((K.EQ.J).AND.(ISYMK.EQ.ISYMJ))
     *                                                     GOTO 904
                                 IF (ONEAUX) THEN
                                   NBK = IH1AM(ISYMB,ISYMK) +
     *                                   NORB1(ISYMB)*(K-1) + B
                                   NAIBK = IH2AM(ISYMKB,ISYMIA) +
     *                                   INDEX(NBK,NAI)
                                 ELSE
                                   NBK = IT1AM(ISYMB,ISYMK) +
     *                                   NVIR(ISYMB)*(K-1) + B
                                   NAIBK = IT2AM(ISYMKB,ISYMIA) +
     *                                   INDEX(NBK,NAI)
                                 END IF
                                 NJK = INDEX(J,K)
                                 XXX = FLMAT(NJK)*TIJAB(NAIBK)
                                 TKJ = TKJ + XXX
  904                         CONTINUE
C
C   LOOP FOR COEFFICIENTS (IJCB)
C
                              ISYMIC = ISYMJB
                              ISYMC  = MULD2H(ISYMI,ISYMIC)
                              DO 906 NC = NRHFS(ISYMC) + 1, NORB1(ISYMC)
                                 IF ((NC.EQ.A).AND.(ISYMA.EQ.ISYMC))
     *                                                     GOTO 906
                                 IF (ONEAUX) THEN
                                   NCI = IH1AM(ISYMC,ISYMI) +
     *                                   NORB1(ISYMC)*(I-1) + NC
                                   NCIBJ = IH2AM(ISYMIC,ISYMJB) +
     *                                   INDEX(NCI,NBJ)
                                 ELSE
                                   NCI = IT1AM(ISYMC,ISYMI) +
     *                                   NVIR(ISYMC)*(I-1) + NC
                                   NCIBJ = IT2AM(ISYMIC,ISYMJB) +
     *                                   INDEX(NCI,NBJ)
                                 END IF
                                 NCA = INDEX(NC-NCT,A-NCT)
                                 XXX = FLMAT(NCA)*TIJAB(NCIBJ)
                                 TKA = TKA + XXX
  906                         CONTINUE
C
C   LOOP FOR COEFFICIENTS (IJAC)
C
                              ISYMJC = ISYMIA
                              ISYMC  = MULD2H(ISYMJ,ISYMJC)
                              DO 908 NC = NRHFS(ISYMC) + 1, NORB1(ISYMC)
                                 IF ((NC.EQ.B).AND.(ISYMB.EQ.ISYMC))
     *                                                     GOTO 908
                                 IF (ONEAUX) THEN
                                   NCJ = IH1AM(ISYMC,ISYMJ) +
     *                                   NORB1(ISYMC)*(J-1) + NC
                                   NAICJ = IH2AM(ISYMJC,ISYMIA) +
     *                                   INDEX(NCJ,NAI)
                                 ELSE
                                   NCJ = IT1AM(ISYMC,ISYMJ) +
     *                                   NVIR(ISYMC)*(J-1) + NC
                                   NAICJ = IT2AM(ISYMJC,ISYMIA) +
     *                                   INDEX(NCJ,NAI)
                                 END IF
                                 NCB = INDEX(NC-NCT,B-NCT)
                                 XXX = FLMAT(NCB)*TIJAB(NAICJ)
                                 TKB = TKB + XXX
  908                         CONTINUE
C
                              NAA = INDEX(A-NCT,A-NCT)
                              NBB = INDEX(B-NCT,B-NCT)
                              NII = INDEX(I,I)
                              NJJ = INDEX(J,J)
                              DENOM = D1 / (FLMAT(NAA) + FLMAT(NBB) -
     *                                      FLMAT(NII) - FLMAT(NJJ))
                              TIJAB(NAIBJ) = (TKI+TKJ-TKA-TKB-VV)
     *                                       * DENOM
                              VV = (D2 * VAIBJ - VAJBI)
			      E2 = E2 + VV * TIJAB(NAIBJ)
                              IF (ABS(TIJAB(NAIBJ)-TINI).GT.THRES) 
     *                                                         ICONV = 1
  804                      CONTINUE
  803                   CONTINUE
  802                CONTINUE
  801             CONTINUE
  703          CONTINUE
  702       CONTINUE
  701    CONTINUE
	 IF (ICOUNT .EQ. 1) WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,I3,3X,F15.9)') ' ENERGY IN ITERATION',
     &   ICOUNT,E2
         IF (ICOUNT.LE.50) THEN
             IF (ICONV.EQ.1) THEN
                 ICOUNT= ICOUNT + 1
                 GOTO 900
             ELSE
                 WRITE(LUPRI,*)
                 WRITE(LUPRI,'(A,I3)') ' FINAL ITERATION OF LOCAL MP2',
     &                                                          ICOUNT
             END IF
         ELSE
            WRITE(LUPRI,'(A)') 'MORE THAN 50 ITERATIONS'
         END IF
      END IF
C
C     CONSTRUCT PRODUCTS (IA|JB)*(KA|LB)
C
      E2 = D0
      E2S = D0
      E2T = D0
      DO 101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         JBOFF = NH1AM(ISYMJB) * (NH1AM(ISYMJB) + 1) / 2
         DO 102 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 103 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
C
C              COMPUTE MP2 ENERGY
C
               DO 201 I = 1, NRHFA(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  KOFFIA = IRHFA(ISYMI) + I
                  DO 202 J = 1, NRHFA(ISYMJ)
                     KOFFJ = IRHF(ISYMJ) + J           
                     KOFFJA = IRHFA(ISYMJ) + J           
                     DO 203 A = NRHFSA(ISYMA) + 1, NORB1(ISYMA)
                        ISYMJA = MULD2H(ISYMJ,ISYMA)
                        KOFFA = IVIR(ISYMA) + A 
                        DO 204 B = NRHFSA(ISYMB) + 1, NORB1(ISYMB)
                           ISYMIB = MULD2H(ISYMI,ISYMB)
                           KOFFB = IVIR(ISYMB) + B
                           IF (ONEAUX) THEN
                              NAI = IH1AM(ISYMA,ISYMI) +
     *                              NORB1(ISYMA)*(I-1) + A
                              NBJ = IH1AM(ISYMB,ISYMJ) +
     *                              NORB1(ISYMB)*(J-1) + B
                              NAIBJ = IH2AM(ISYMIA,ISYMJB) +
     *                                INDEX(NAI,NBJ)
                              NAJ = IH1AM(ISYMA,ISYMJ) +
     *                              NORB1(ISYMA)*(J-1) + A
                              NBI = IH1AM(ISYMB,ISYMI) +
     *                              NORB1(ISYMB)*(I-1) + B
                              NAJBI = IH2AM(ISYMJA,ISYMIB) +
     *                                INDEX(NAJ,NBI)
                           ELSE
                              NAI = IT1AM(ISYMA,ISYMI) +
     *                              NVIR(ISYMA)*(I-1) + A
                              NBJ = IT1AM(ISYMB,ISYMJ) +
     *                              NVIR(ISYMB)*(J-1) + B
                              NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                                INDEX(NAI,NBJ)
                              NAJ = IT1AM(ISYMA,ISYMJ) +
     *                              NVIR(ISYMA)*(J-1) + A
                              NBI = IT1AM(ISYMB,ISYMI) +
     *                              NVIR(ISYMB)*(I-1) + B
                              NAJBI = IT2AM(ISYMJA,ISYMIB) +
     *                                INDEX(NAJ,NBI)
                           END IF
                           VAIBJ = V2AM(NAIBJ)
                           VAJBI = V2AM(NAJBI)
C
                           IF (LGLOCAL) THEN
C LOCALIZED MP2
                              VV = (D2 * VAIBJ - VAJBI)
                              DENOM = TIJAB(NAIBJ)
                              VS = (VAIBJ + VAJBI)
                              VT = (VAIBJ - VAJBI)
                           ELSE
C CANONICAL MP2
                              VV = VAIBJ * (D2 * VAIBJ - VAJBI)
                              DENOM = D1 / (EV(KOFFI) + EV(KOFFJ) -
     *                                      EV(KOFFA) - EV(KOFFB))
                              VS = (VAIBJ + VAJBI)**2
                              VT = (VAIBJ - VAJBI)**2
                           END IF
                           E2 = E2 + VV * DENOM
                           IJ = INDEX(KOFFIA,KOFFJA)
                           ES(IJ) = ES(IJ) + VS * DENOM
                           ET(IJ) = ET(IJ) + VT * DENOM
  204                   CONTINUE
  203                CONTINUE
  202             CONTINUE
  201          CONTINUE
  103       CONTINUE
  102    CONTINUE
  101 CONTINUE
      
C
      IF (LGLOCAL) THEN
C
C        WRITE AMPLITUDES AND RESTORE R12 INTEGRALS
C
         LU43 = -43
         CALL GPOPEN(LU43,'TIJAB',
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         WRITE(LU43) (R2AM(I), I = 1,NH2AMX)
         CALL GPCLOSE(LU43,'KEEP')
         CALL GPOPEN(LU43,FR12R12,
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         READ(LU43) (R2AM(I), I = 1,NH2AMX)
         CALL GPCLOSE(LU43,'KEEP')
C
      END IF
C
      DO 1101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         JBOFF = NH1AM(ISYMJB) * (NH1AM(ISYMJB) + 1) / 2
         DO 1102 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 1103 ISYMJ = 1, NSYM
              ISYMB = MULD2H(ISYMJ,ISYMJB)
C
              IF (ONEAUX) THEN
                DO 1104 ISYMKA = 1, NSYM
                  ISYMLB = ISYMKA 
                  ISYMK = MULD2H(ISYMA,ISYMKA)
                  ISYML = MULD2H(ISYMB,ISYMLB)
                  DO 1105 I = 1, NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     DO 1106 K = 1, NRHF(ISYMK)
                        KOFFK = IRHF(ISYMK) + K
                        DO 1107 A = 1, NORB1(ISYMA)
                           NAI = IH1AM(ISYMA,ISYMI) +
     *                           NORB1(ISYMA)*(I-1) + A
                           NAK = IH1AM(ISYMA,ISYMK) +
     *                           NORB1(ISYMA)*(K-1) + A
                           DO 1108 J = 1, NRHF(ISYMJ)
                              KOFFJ = IRHF(ISYMJ) + J
                              DO 1109 L = 1, NRHF(ISYML)
                                 KOFFL = IRHF(ISYML) + L
                                 DO 1110 B = 1, NORB1(ISYMB)
                                    NBJ = IH1AM(ISYMB,ISYMJ) +
     *                                    NORB1(ISYMB)*(J-1) + B
                                    NBL = IH1AM(ISYMB,ISYML) +
     *                                    NORB1(ISYMB)*(L-1) + B
                                    NAIBJ = IH2AM(ISYMIA,ISYMJB) +
     *                                      INDEX(NAI,NBJ)
                                    NAKBL = IH2AM(ISYMKA,ISYMLB) +
     *                                      INDEX(NAK,NBL)      
                                    R1 = R2AM(NAIBJ)
                                    V1 = V2AM(NAIBJ)
                                    U1 = U2AM(NAIBJ)
                                    R2 = R2AM(NAKBL)
                                    U2 = U2AM(NAKBL)
                                    VR = V1 * R2
                                    UR = U1 * R2 + U2 * R1
                                    RR = R1 * R2 
                                    V11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              V11(KOFFI,KOFFJ,KOFFK,KOFFL) + VR
                                    U11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              U11(KOFFI,KOFFJ,KOFFK,KOFFL) + UR
                                    Q11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
 1110                            CONTINUE
 1109                         CONTINUE
 1108                      CONTINUE
 1107                   CONTINUE
 1106                CONTINUE
 1105             CONTINUE
 1104           CONTINUE
                DO 2104 ISYMKA = 1, NSYM
                  ISYMLB = ISYMKA 
                  LBOFF = NH1AM(ISYMLB) * (NH1AM(ISYMLB) + 1) / 2
                  ISYMK = MULD2H(ISYMA,ISYMKA)
                  ISYML = MULD2H(ISYMB,ISYMLB)
                  DO 2105 I = 1, NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     DO 2106 K = 1, NRHF(ISYMK)
                        KOFFK = IRHF(ISYMK) + K
                        DO 2107 A = 1, NORB2(ISYMA)
                           NAI = IG1AM(ISYMA,ISYMI) +
     *                           NORB2(ISYMA)*(I-1) + A
                           NAK = IG1AM(ISYMA,ISYMK) +
     *                           NORB2(ISYMA)*(K-1) + A
                           DO 2108 J = 1, NRHF(ISYMJ)
                              KOFFJ = IRHF(ISYMJ) + J
                              DO 2109 L = 1, NRHF(ISYML)
                                 KOFFL = IRHF(ISYML) + L
                                 DO 2110 B = 1, NORB1(ISYMB)
                                    NBJ = IH1AM(ISYMB,ISYMJ) +
     *                                    NORB1(ISYMB)*(J-1) + B
                                    NBL = IH1AM(ISYMB,ISYML) +
     *                                    NORB1(ISYMB)*(L-1) + B
                                    NAIBJ = IH2AM(ISYMIA,ISYMJB) +
     *                              NH1AM(ISYMJB)*(NAI - 1) + NBJ + 
     *                              JBOFF
                                    NAKBL = IH2AM(ISYMKA,ISYMLB) +
     *                              NH1AM(ISYMLB)*(NAK - 1) + NBL + 
     *                              LBOFF
                                    R1 = R2AM(NAIBJ)
                                    V1 = V2AM(NAIBJ)
                                    U1 = U2AM(NAIBJ)
                                    R2 = R2AM(NAKBL)
                                    U2 = U2AM(NAKBL)
                                    VR = V1 * R2
                                    UR = U1 * R2 + U2 * R1
                                    RR = R1 * R2 
                                    W11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              W11(KOFFI,KOFFJ,KOFFK,KOFFL) + VR
                                    B11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              B11(KOFFI,KOFFJ,KOFFK,KOFFL) + UR
                                    R11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              R11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
                                    W11(KOFFJ,KOFFI,KOFFL,KOFFK) = 
     *                              W11(KOFFJ,KOFFI,KOFFL,KOFFK) + VR
                                    B11(KOFFJ,KOFFI,KOFFL,KOFFK) = 
     *                              B11(KOFFJ,KOFFI,KOFFL,KOFFK) + UR
                                    R11(KOFFJ,KOFFI,KOFFL,KOFFK) = 
     *                              R11(KOFFJ,KOFFI,KOFFL,KOFFK) + RR
 2110                            CONTINUE
 2109                         CONTINUE
 2108                      CONTINUE
 2107                   CONTINUE
 2106                CONTINUE
 2105             CONTINUE
 2104           CONTINUE
              ELSE
                DO 104 ISYMKA = 1, NSYM
                  ISYMLB = ISYMKA 
                  ISYMK = MULD2H(ISYMA,ISYMKA)
                  ISYML = MULD2H(ISYMB,ISYMLB)
                  DO 105 I = 1, NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     DO 106 K = 1, NRHF(ISYMK)
                        KOFFK = IRHF(ISYMK) + K
                        DO 107 A = 1, NVIR(ISYMA)
                           NAI = IT1AM(ISYMA,ISYMI) +
     *                           NVIR(ISYMA)*(I-1) + A
                           NAK = IT1AM(ISYMA,ISYMK) +
     *                           NVIR(ISYMA)*(K-1) + A
                           DO 108 J = 1, NRHF(ISYMJ)
                              KOFFJ = IRHF(ISYMJ) + J
                              DO 109 L = 1, NRHF(ISYML)
                                 KOFFL = IRHF(ISYML) + L
                                 DO 110 B = 1, NVIR(ISYMB)
                                    NBJ = IT1AM(ISYMB,ISYMJ) +
     *                                    NVIR(ISYMB)*(J-1) + B
                                    NBL = IT1AM(ISYMB,ISYML) +
     *                                    NVIR(ISYMB)*(L-1) + B
                                    NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                                      INDEX(NAI,NBJ)
                                    NAKBL = IT2AM(ISYMKA,ISYMLB) +
     *                                      INDEX(NAK,NBL)
                                    R1 = R2AM(NAIBJ)
                                    V1 = V2AM(NAIBJ)
                                    U1 = U2AM(NAIBJ)
                                    R2 = R2AM(NAKBL)
                                    U2 = U2AM(NAKBL)
                                    VR = V1 * R2
                                    UR = U1 * R2 + U2 * R1
                                    RR = R1 * R2 
C      
                                    IF (A .GT. NORB1(ISYMA)) GOTO 112
                                    IF (B .GT. NORB1(ISYMB)) GOTO 112
C 
                                    V11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              V11(KOFFI,KOFFJ,KOFFK,KOFFL) +  VR
                                    U11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              U11(KOFFI,KOFFJ,KOFFK,KOFFL) +  UR
                                    Q11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +  RR
                                    GOTO 111
  112                               CONTINUE
C       
                                    IF (A .GT. NORB1(ISYMA) .AND.
     *                                  B .GT. NORB1(ISYMB)) GOTO 111
C
                                    W11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              W11(KOFFI,KOFFJ,KOFFK,KOFFL) +  VR
                                    B11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              B11(KOFFI,KOFFJ,KOFFK,KOFFL) +  UR
                                    R11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              R11(KOFFI,KOFFJ,KOFFK,KOFFL) +  RR
C
  111                               CONTINUE
  110                            CONTINUE
  109                         CONTINUE
  108                      CONTINUE
  107                   CONTINUE
  106                CONTINUE
  105             CONTINUE
  104           CONTINUE
              ENDIF
 1103       CONTINUE
 1102    CONTINUE
 1101 CONTINUE
C
      WRITE(LUPRI,'(/A/)') ' SECOND-ORDER PAIR ENERGIES:'
      DO 230 I=1,NAPAIR
         IF (LGLOCAL) THEN
             ES(I) = ES(I) * 0.5D0
             ET(I) = ET(I) * 1.5D0
         ELSE
             ES(I) = ES(I) * DP25
             ET(I) = ET(I) * DP75
         END IF
         E2S = E2S + ES(I)
         E2T = E2T + ET(I)
        WRITE(LUPRI,'(17X,I4,2F15.9)') I,ES(I),ET(I)
  230 CONTINUE
      WRITE(LUPRI,'(/A6,3F15.9)') ' MP2 =',E2,E2S,E2T
C
      IJ = 0
      DO 300 I = 1, NOCDIM
         DO 301 J = 1, I
            IJ = IJ + 1
            KL = 0
            DO 302 K = 1, NOCDIM
               DO 303 L = 1, K
                  KL = KL + 1
C
                  IF (R12EOR) THEN
                     RR = - D2 *(QQ6(I,J,K,L) + QQ6(I,J,L,K))
                     RR = RR - (U11(I,J,K,L) + U11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (U11(I,J,K,L) + U11(I,J,L,K))
                     DX = D2
                  END IF
                  US11(KL,IJ) = RR 
                  IF (I .EQ. J) US11(KL,IJ) = FF * US11(KL,IJ)
                  IF (K .EQ. L) US11(KL,IJ) = FF * US11(KL,IJ)
                  IF (IJ .EQ. KL) US11(KL,IJ) = US11(KL,IJ) - DX
                  IF (R12EOR) THEN
                     RR = - D2 *(QQ6(I,J,K,L) - QQ6(I,J,L,K))
                     RR = RR - (U11(I,J,K,L) - U11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (U11(I,J,K,L) - U11(I,J,L,K))
                     DX = D2
                  END IF
                  UT11(KL,IJ) = RR
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            UT11(KL,IJ) = UT11(KL,IJ) - DX
C
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) + QQ2(I,J,L,K))
                     RR = RR - (V11(I,J,K,L) + V11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (V11(I,J,K,L) + V11(I,J,L,K))
                     DX = D1
                  END IF
                  VS11(KL,IJ) = RR
                  IF (I .EQ. J) VS11(KL,IJ) = FF * VS11(KL,IJ)
                  IF (K .EQ. L) VS11(KL,IJ) = FF * VS11(KL,IJ)
                  IF (IJ .EQ. KL) VS11(KL,IJ) = VS11(KL,IJ) + DX
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) - QQ2(I,J,L,K))
                     RR = RR - (V11(I,J,K,L) - V11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (V11(I,J,K,L) - V11(I,J,L,K))
                     DX = D1
                  END IF
                  VT11(KL,IJ) = RR
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            VT11(KL,IJ) = VT11(KL,IJ) + DX
C
                  IF (R12EOR) THEN
                     RR = - D2 *(QQ6(I,J,K,L) + QQ6(I,J,L,K))
                     RR = RR + (U11(I,J,K,L) + U11(I,J,L,K)) * FACX
                     RR = RR - (B11(I,J,K,L) + B11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =   (U11(I,J,K,L) + U11(I,J,L,K)) * FACX
     *                    - (B11(I,J,K,L) + B11(I,J,L,K))
                     DX = D2
                  END IF
                  BS11(KL,IJ) = RR
                  IF (I .EQ. J) BS11(KL,IJ) = FF * BS11(KL,IJ)
                  IF (K .EQ. L) BS11(KL,IJ) = FF * BS11(KL,IJ) 
                  IF (IJ .EQ. KL) BS11(KL,IJ) = BS11(KL,IJ) - DX 
                  IF (R12EOR) THEN
                     RR = - D2 *(QQ6(I,J,K,L) - QQ6(I,J,L,K))
                     RR = RR + (U11(I,J,K,L) - U11(I,J,L,K)) * FACX
                     RR = RR - (B11(I,J,K,L) - B11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =   (U11(I,J,K,L) - U11(I,J,L,K)) * FACX
     *                    - (B11(I,J,K,L) - B11(I,J,L,K))
                     DX = D2
                  END IF
                  BT11(KL,IJ) = RR
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            BT11(KL,IJ) = BT11(KL,IJ) - DX
C
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) + QQ2(I,J,L,K))
                     RR = RR + (V11(I,J,K,L) + V11(I,J,L,K)) * FACX
                     RR = RR - (W11(I,J,K,L) + W11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =    (V11(I,J,K,L) + V11(I,J,L,K)) * FACX
     *                     - (W11(I,J,K,L) + W11(I,J,L,K))
                     DX = D1
                  END IF
                  WS11(KL,IJ) = RR
                  IF (I .EQ. J) WS11(KL,IJ) = FF * WS11(KL,IJ)
                  IF (K .EQ. L) WS11(KL,IJ) = FF * WS11(KL,IJ) 
                  IF (IJ .EQ. KL) WS11(KL,IJ) = WS11(KL,IJ) + DX 
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) - QQ2(I,J,L,K))
                     RR = RR + (V11(I,J,K,L) - V11(I,J,L,K)) * FACX
                     RR = RR - (W11(I,J,K,L) - W11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =   (V11(I,J,K,L) - V11(I,J,L,K)) * FACX
     *                    - (W11(I,J,K,L) - W11(I,J,L,K))
                     DX = D1
                  END IF
                  WT11(KL,IJ) = RR
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            WT11(KL,IJ) = WT11(KL,IJ) + DX
C
                  IF (R12EOR) THEN
                     QQIJKL = QQ4(I,J,K,L)             
                     QQIJLK = QQ4(I,J,L,K)
                  ELSE
                     QQIJKL = QQ11(IQ(I),IQ(J),IQ(K),IQ(L))
                     QQIJLK = QQ11(IQ(I),IQ(J),IQ(L),IQ(K))
                  END IF
                  RR = QQIJKL + QQIJLK
                  RR = RR - (Q11(I,J,K,L) + Q11(I,J,L,K))
                  QS11(KL,IJ) = RR
                  IF (I .EQ. J) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  IF (K .EQ. L) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  RR = QQIJKL - QQIJLK
                  RR = RR - (Q11(I,J,K,L) - Q11(I,J,L,K))
                  QT11(KL,IJ) = RR
C
                  RR = QQIJKL + QQIJLK
                  RR = RR + (Q11(I,J,K,L) + Q11(I,J,L,K)) * FACX
     *                    - (R11(I,J,K,L) + R11(I,J,L,K))
ccn                  rr = + (Q11(I,J,K,L) + Q11(I,J,L,K)) * FACX
                  RS11(KL,IJ) = RR
                  IF (I .EQ. J) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  IF (K .EQ. L) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  RR = QQIJKL - QQIJLK
                  RR = RR + (Q11(I,J,K,L) - Q11(I,J,L,K)) * FACX
     *                    - (R11(I,J,K,L) - R11(I,J,L,K))
ccn                  rr = + (Q11(I,J,K,L) - Q11(I,J,L,K)) * FACX
                  RT11(KL,IJ) = RR
C 
                  US11(KL,IJ) = US11(KL,IJ) * DP5 * DP25
                  BS11(KL,IJ) = BS11(KL,IJ) * DP5 * DP25
                  VS11(KL,IJ) = VS11(KL,IJ) * DP5
                  WS11(KL,IJ) = WS11(KL,IJ) * DP5
                  QS11(KL,IJ) = QS11(KL,IJ) * DP25
                  RS11(KL,IJ) = RS11(KL,IJ) * DP25
C
                  UT11(KL,IJ) = UT11(KL,IJ) * D1P5 * DP25 
                  BT11(KL,IJ) = BT11(KL,IJ) * D1P5 * DP25 
                  VT11(KL,IJ) = VT11(KL,IJ) * D1P5
                  WT11(KL,IJ) = WT11(KL,IJ) * D1P5
                  QT11(KL,IJ) = QT11(KL,IJ) * DP75
                  RT11(KL,IJ) = RT11(KL,IJ) * DP75
c
  303          CONTINUE
  302       CONTINUE
  301    CONTINUE
  300 CONTINUE
C
      IF (R12OLD) THEN
         CALL DCOPY(NSPAIR*NSPAIR,VS11,1,WS11,1)
         CALL DCOPY(NSPAIR*NSPAIR,US11,1,BS11,1)
         CALL DCOPY(NSPAIR*NSPAIR,QS11,1,RS11,1)
         CALL DCOPY(NSPAIR*NSPAIR,VT11,1,WT11,1)
         CALL DCOPY(NSPAIR*NSPAIR,UT11,1,BT11,1)
         CALL DCOPY(NSPAIR*NSPAIR,QT11,1,RT11,1)
      ENDIF
C
      CALL VSHRNK(VS11,NSPAIR,NRHF,NRHFA,NSYM)
      CALL VSHRNK(VT11,NSPAIR,NRHF,NRHFA,NSYM)
      CALL VSHRNK(WS11,NSPAIR,NRHF,NRHFA,NSYM)
      CALL VSHRNK(WT11,NSPAIR,NRHF,NRHFA,NSYM)
C
      IF (IPRINT .LE. 3) GOTO 998
      CALL AROUND('OLD singlet <V> matrix') 
      CALL OUTPUT(VS11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD singlet <B> matrix') 
      CALL OUTPUT(US11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD singlet <S> matrix') 
      CALL OUTPUT(QS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD triplet <V> matrix') 
      CALL OUTPUT(VT11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD triplet <B> matrix') 
      CALL OUTPUT(UT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD triplet <S> matrix') 
      CALL OUTPUT(QT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      IF (R12OLD) GOTO 998
      CALL AROUND('NEW singlet <V> matrix') 
      CALL OUTPUT(WS11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW singlet <B> matrix') 
      CALL OUTPUT(BS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW singlet <S> matrix') 
      CALL OUTPUT(RS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW triplet <V> matrix') 
      CALL OUTPUT(WT11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW triplet <B> matrix') 
      CALL OUTPUT(BT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW triplet <S> matrix') 
      CALL OUTPUT(RT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
  998 CONTINUE
      CALL GPOPEN (LUMULO,'AOR12O','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      CALL GPOPEN (LUMULN,'AOR12N','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMULO,'(4E30.20)') VS11
      WRITE(LUMULO,'(4E30.20)') US11
      WRITE(LUMULO,'(4E30.20)') QS11
      WRITE(LUMULO,'(4E30.20)') VT11
      WRITE(LUMULO,'(4E30.20)') UT11
      WRITE(LUMULO,'(4E30.20)') QT11
      WRITE(LUMULN,'(4E30.20)') WS11
      WRITE(LUMULN,'(4E30.20)') BS11
      WRITE(LUMULN,'(4E30.20)') RS11
      WRITE(LUMULN,'(4E30.20)') WT11
      WRITE(LUMULN,'(4E30.20)') BT11
      WRITE(LUMULN,'(4E30.20)') RT11
C
      LU43 = -43
      LU44 = -44
      LU45 = -45
      IF (R12OLD) THEN
        LUMULA = LUMULO
        LUMULE = LUMULO  
      ELSE
        IF (R12XXL) THEN
           LUMULA = LUMULO
        ELSE
           LUMULA = LUMULN
        END IF
        LUMULE = LUMULN  
      END IF
      DO 999 LUMULX = LUMULA, LUMULE
      IF (LUMULX .EQ. LUMULO) THEN
        CALL AROUND('Original MP2-R12 method')
        IF (.NOT. R12NOB .OR. NATVIR) THEN
          CALL GPOPEN(LU43,'QMATSO','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') US11
          CALL GPCLOSE(LU43,'KEEP')
          CALL GPOPEN(LU43,'QMATTO','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') UT11
          CALL GPCLOSE(LU43,'KEEP')
        END IF
      ELSE
        CALL AROUND('Auxiliary-basis MP2-R12 method')
        IF (.NOT. R12NOB .OR. NATVIR) THEN
          CALL GPOPEN(LU43,'QMATSN','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') US11
          CALL GPCLOSE(LU43,'KEEP')
          CALL GPOPEN(LU43,'QMATTN','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') UT11
          CALL GPCLOSE(LU43,'KEEP')
          IF (.NOT. NORXR) THEN
             CALL GPOPEN(LU43,'XMATSN','UNKNOWN',' ','FORMATTED',
     &                   IDUM,LDUM)
             READ(LU43,'(4E30.20)') RXRS
             CALL GPCLOSE(LU43,'KEEP')
             CALL GPOPEN(LU43,'XMATTN','UNKNOWN',' ','FORMATTED',
     &                   IDUM,LDUM)
             READ(LU43,'(4E30.20)') RXRT
             CALL GPCLOSE(LU43,'KEEP')
          END IF
        END IF
      END IF   
      REWIND LUMULX
      READ(LUMULX,'(4E30.20)') WS11
      READ(LUMULX,'(4E30.20)') BS11
      READ(LUMULX,'(4E30.20)') RS11
      READ(LUMULX,'(4E30.20)') WT11
      READ(LUMULX,'(4E30.20)') BT11
      READ(LUMULX,'(4E30.20)') RT11
      LUM = LUMULX
      CALL GPCLOSE(LUM,'KEEP')
celena
      IF (R12PRP) THEN
          CALL GPOPEN(LU43,'CCR12_XP','UNKNOWN',' ','FORMATTED',
     &                IDUM,LDUM)
          WRITE(LU43,'(I3)') 1
          WRITE(LU43,'(4E30.20)') RS11
          WRITE(LU43,'(4E30.20)') RT11
          CALL GPCLOSE(LU43,'KEEP')
       ENDIF
celena

c     Write out V-matrix, etc. for CC2-R12 model (HF/UniKA/02-05-2003).
c     modified by C. Neiss: no not use singlet/triplet format any more,
c     make intermediates compatible with correlation factor r_12 (and 
c     not 0.5*r_12 as in the MP2-R12 code)!
c     (April 2005)
c
      IF (CCR12) THEN
        KSCR = 1
        KEND1 = KSCR + NGAMMA(1) !NGAMMA(1) in MP2-R12 >= NKILJ
c       WRITE(LUPRI,*) 'NGAMMA(1) in CC_R12: ', NGAMMA(1)
        IF (NGAMMA(1) .gt. LWORK) THEN
          CALL QUIT('Insufficient WORK space in MAKEVR')
        END IF
        CALL DSCAL(NSPAIR*NSPAIR,2.0D0,WS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,2.0D0,WT11,1)
        CALL DSCAL(NSPAIR*NSPAIR,1.0D0/3.0D0,WT11,1)
        CALL CCR12PCK(WORK(KSCR),1,WS11,WT11,NRHF,NRHFA,NKILJ)
c       WRITE(LUPRI,*) 'VIJKL written in MP2R12:'
c       CALL CC_PRPR12(WORK(KSCR),1,1,.FALSE.)
        CALL GPOPEN(LU43,FCCR12V,'UNKNOWN',' ','UNFORMATTED',
     &              IDUM,LDUM)
        WRITE(LU43) 1
        WRITE(LU43) (WORK(KSCR-1+I), I=1, NKILJ(1))
        CALL GPCLOSE(LU43,'KEEP')  
        ! undo scaling
        CALL DSCAL(NSPAIR*NSPAIR,3.0D0,WT11,1)
        CALL DSCAL(NSPAIR*NSPAIR,0.5D0,WS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,0.5D0,WT11,1)
c
c       Write out X-matrix for CC2-R12 model (HF/UniKA/02-05-2003).
c
        CALL DSCAL(NSPAIR*NSPAIR,4.0D0,RS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,4.0D0,RT11,1)
        CALL DSCAL(NSPAIR*NSPAIR,1.0D0/3.0D0,RT11,1)
        CALL CCR12PCK(WORK(KSCR),1,RS11,RT11,NRHF,NRHF,NKILJ)
        CALL GPOPEN(LU43,FCCR12X,'UNKNOWN',' ','UNFORMATTED',
     &              IDUM,LDUM)
        WRITE(LU43) 1
        WRITE(LU43) (WORK(KSCR-1+I), I=1, NKILJ(1))
        CALL GPCLOSE(LU43,'KEEP')
        ! undo scaling
        CALL DSCAL(NSPAIR*NSPAIR,3.0D0,RT11,1)
        CALL DSCAL(NSPAIR*NSPAIR,0.25D0,RS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,0.25D0,RT11,1)
c
c       Write out orbital energies for CC2-R12 model (HF/UniKA/02-05-2003).
c
        CALL GPOPEN(LU43,FCCR12E,'UNKNOWN',' ','FORMATTED',
     &           IDUM,LDUM)
        WRITE(LU43,'(4E30.20)') EVL
        CALL GPCLOSE(LU43,'KEEP')
      END IF
c        
c     Open files for B and C (HF/UniKA/25-06-2003).
c
      IF (CCR12 .AND. (LUMULA .EQ. LUMULX)) THEN
        CALL GPOPEN(LU44,FCCR12B,'UNKNOWN',' ','UNFORMATTED',
     &              IDUM,LDUM)
        CALL GPOPEN(LU45,FCCR12D,'UNKNOWN',' ','UNFORMATTED',
     &              IDUM,LDUM)
      END IF
C
      DO 999 IPDD = 0, 7
chf
        IF (R12HYB .OR. NORXR) THEN
          IPCC = APPROX(IPDD,1)
        ELSE
          IPCC = APPROX(IPDD,2)
        END IF
chf
      IF (NATVIR .AND. IPDD .EQ. 3) THEN
        DO I = 1, NSPAIR
          DO J = 1, NSPAIR
            BS11(I,J) = BS11(I,J) - US11(I,J)
            BT11(I,J) = BT11(I,J) - UT11(I,J)
          END DO
        END DO
      END IF
      IF (IPDD .EQ. 4 .AND. .NOT. NATVIR .AND. .NOT. R12NOB) THEN
        DO I = 1, NSPAIR
          DO J = 1, NSPAIR
            BS11(I,J) = BS11(I,J) - US11(I,J)
            BT11(I,J) = BT11(I,J) - UT11(I,J)
          END DO
        END DO
      END IF
      IF (IPDD .EQ. 6) THEN
        IF (.NOT. R12NOB .AND..NOT. NORXR .AND. LUMULX .EQ. LUMULN) THEN
          DO I = 1, NSPAIR
            DO J = 1, NSPAIR
              BS11(I,J) = BS11(I,J) + RXRS(I,J)
              BT11(I,J) = BT11(I,J) + RXRT(I,J)
            END DO
          END DO
        END IF
      END IF
      IF (LUMULX .EQ. LUMULO .OR. NORXR) THEN
        IF (.NOT. R12XXL .AND. 
     *      (IPDD .NE. 2 .AND. IPDD .NE. 5)) GOTO 999
      ELSE
        IF (.NOT. R12XXL .AND. 
     *      (IPDD .NE. 2 .AND. IPDD .NE. 7)) GOTO 999
      END IF
      IF (IPDD .LE. 3 .AND. R12NOA) GOTO 999
      IF (IPDD .EQ. 3 .AND. R12NOP) GOTO 999
      IF (IPDD .GE. 4 .AND. R12NOB) GOTO 999
      IF (IPDD .GE. 6 .AND. (LUMULX .EQ. LUMULO .OR. NORXR)) GOTO 999 
      IF (IPDD .EQ. 0) THEN
         XPDD = D0
         WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/0'//
     *                   ' PAIR ENERGIES WITH FIXED R12-COEFFICIENTS:'
      ELSE IF (IPDD .EQ. 1) THEN
        XPDD = D0
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/A'//
     *                   ' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 2) THEN
        XPDD = D0
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/A PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 3) THEN 
        XPDD = DP5
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/A'' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 4) THEN
        XPDD = D0    
        IF (R12HYB .OR. (LUMULX .EQ. LUMULO)) THEN
           WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B'//
     *                   ' PAIR ENERGIES:'
        ELSE
           WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B'''//
     *                   ' PAIR ENERGIES:'
        END IF
      ELSE IF (IPDD .EQ. 5) THEN
        XPDD = DP5    
        IF (R12HYB .OR. (LUMULX .EQ. LUMULO)) THEN
           WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B PAIR ENERGIES:'
        ELSE
           WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B'' PAIR ENERGIES:'
        END IF
      ELSE IF (IPDD .EQ. 6) THEN  
        XPDD = D0    
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B'//
     *                   ' PAIR ENERGIES:' 
      ELSE IF (IPDD .EQ. 7) THEN
        XPDD = DP5    
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B PAIR ENERGIES:'
      END IF
      F2S = D0
      F2T = D0
      CSALL = -D1    
      CTALL = -D1
c
c     Write out B-matrix for CC2-R12 model (HF/UniKA/02-05-2003).
c
      IF (CCR12) THEN 
        IF (LUMULX .EQ. LUMULO) THEN
          WRITE(LU44) 0,IPDD,IPCC
        ELSE
          WRITE(LU44) 1,IPDD,IPCC
        END IF
        CALL DSCAL(NSPAIR*NSPAIR,4.0D0,BS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,4.0D0,BT11,1)
        CALL DSCAL(NSPAIR*NSPAIR,1.0D0/3.0D0,BT11,1)
        CALL CCR12PCK(WORK(KSCR),1,BS11,BT11,NRHF,NRHF,NKILJ)
        WRITE(LU44) (WORK(KSCR-1+I), I=1, NKILJ(1))
        ! undo scaling
        CALL DSCAL(NSPAIR*NSPAIR,3.0D0,BT11,1)
        CALL DSCAL(NSPAIR*NSPAIR,0.25D0,BS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,0.25D0,BT11,1)
      END IF
C
C CCMS
C   
      IF (LGLOCAL .AND. (IPDD.EQ.3 .OR. IPDD.EQ.5 .OR. IPDD.EQ.7) ) THEN
      CALL R12FXL(OLWS,OLWT,FLMAT,RS11,RT11,NIND,III,JJJ,NOCDIM,NSPAIR)
      IF (IPRINT .GE. 5) THEN
         CALL AROUND('Singlet <W> matrix') 
         CALL OUTPUT(OLWS,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <W> matrix') 
         CALL OUTPUT(OLWT,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      END IF
      IF (R12RST) THEN
         WRITE(CFN,1001) IPDD,LUMULX
         CALL GPOPEN(LU43,CFN,'UNKNOWN',' ','FORMATTED',IDUM,LDUM)
         READ(LU43,'(4E30.20)') CS11
         CALL GPCLOSE(LU43,'KEEP')
         WRITE(CFN,1002) IPDD,LUMULX
         CALL GPOPEN(LU43,CFN,'UNKNOWN',' ','FORMATTED',IDUM,LDUM)
         READ(LU43,'(4E30.20)') CT11
         CALL GPCLOSE(LU43,'KEEP')
         ITER = 1
      ELSE
         ITER = 0
      END IF
  800 CONTINUE
      IF (ITER .GT. 0) THEN
         CALL R12FCL(OLZS,OLZT,OLWS,OLWT,WS11,WT11,OLYS,OLYT,RS11,RT11,
     *               CS11,CT11,FLMAT,NIND,III,JJJ,NOCDIM,NSPAIR,DELTA)
         IF (IPRINT .GE. 5) THEN
            CALL AROUND('Singlet <Z> matrix') 
            CALL OUTPUT(OLZS,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
            CALL AROUND('Triplet <Z> matrix') 
            CALL OUTPUT(OLZT,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         END IF
      END IF
      F2SI = F2S
      F2TI = F2T
      F2S = D0
      F2T = D0
      IJ = 0
      DO 500 I = 1, NACDIM
         DO 501 J = 1, I
            IJ = IJ + 1
            IF (ITER .EQ. 0) THEN
              IF (R12SVD) THEN
                CALL R12MP2(WS11(1,IJ),WT11(1,IJ),WS11(1,IJ),WT11(1,IJ),
     *          RS11,RT11,EVL(IJ),NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,
     *          QS11,QT11,BS11,BT11,EVL,XPDD,QQ11,IJ,WSMIN,WTMIN,D0,
     *                                            CS11(1,IJ),CT11(1,IJ))
              ELSE IF (R12DIA) THEN
                CALL MP2R12(WS11(1,IJ),WT11(1,IJ),WS11(1,IJ),WT11(1,IJ),
     *          RS11,RT11,EVL(IJ),NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,
     *          QS11,QT11,BS11,BT11,EVL,XPDD,QQ11,IJ,WSMIN,WTMIN,D0,
     *                                            CS11(1,IJ),CT11(1,IJ))
              ELSE
                CALL QUIT('No solver selected (R12SVD or R12DIA)')
              END IF
            ELSE
              IF (R12SVD) THEN
                CALL R12MP2(OLZS(1,IJ),OLZT(1,IJ),WS11(1,IJ),WT11(1,IJ),
     *          RS11,RT11,EVL(IJ),NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,
     *          QS11,QT11,BS11,BT11,EVL,XPDD,QQ11,IJ,WSMIN,WTMIN,DELTA,
     *                                            CS11(1,IJ),CT11(1,IJ))
              ELSE IF (R12DIA) THEN
                CALL MP2R12(OLZS(1,IJ),OLZT(1,IJ),WS11(1,IJ),WT11(1,IJ),
     *          RS11,RT11,EVL(IJ),NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,
     *          QS11,QT11,BS11,BT11,EVL,XPDD,QQ11,IJ,WSMIN,WTMIN,DELTA,
     *                                            CS11(1,IJ),CT11(1,IJ))
              ELSE
                CALL QUIT('No solver selected (R12SVD or R12DIA)')
              END IF
            END IF
            F2S = F2S + FS(IJ)
            F2T = F2T + FT(IJ) 
  501    CONTINUE
  500 CONTINUE 
      IF (IPRINT .GE. 4) THEN
         CALL AROUND('Singlet <C> matrix') 
         CALL OUTPUT(CS11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <C> matrix') 
         CALL OUTPUT(CT11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
      END IF
      CSALLI = CSALL
      CTALLI = CTALL
      CSALL = DDOT(NSPAIR*NAPAIR,CS11,1,CS11,1)
      CTALL = DDOT(NSPAIR*NAPAIR,CT11,1,CT11,1)
      CSALL = SQRT(CSALL)
      CTALL = SQRT(CTALL)
C
      IF (ITER.LT.100) THEN
C         IF (ABS(F2S-F2SI).GT.1.0d-10.OR.ABS(F2T-F2TI).GT.1.0d-10) THEN
         IF (ABS(CSALLI-CSALL).GT.1.0d-8.
     *              OR.ABS(CTALLI-CTALL).GT.1.0d-8) THEN
            WRITE(LUPRI,'(I3,A,6F15.9)') ITER,'@@',
     *                                         E2S,E2S+F2S,E2T,E2T+F2T,
     *                                         CSALL, CTALL
            ITER = ITER + 1
            GO TO 800
         ELSE
            WRITE(LUPRI,'(I3,A,6F15.9/)') ITER,'@@',
     *                                         E2S,E2S+F2S,E2T,E2T+F2T,
     *                                         CSALL, CTALL
            DO IJ = 1, NAPAIR
              WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                             ET(IJ),ET(IJ)+FT(IJ),
     *                                             WSMIN,WTMIN
            END DO
            WRITE(LUPRI,'(/A,I4,A)')' CONVERGED IN:',ITER,'  ITERATIONS'
         END IF
      ELSE
         WRITE(LUPRI,'(A)') ' ITERATIONS FOR THE ENERGY EXCEEDED 100' 
      END IF
      IF (IPRINT .GE. 4) THEN
         CALL AROUND('Singlet <C> matrix') 
         CALL OUTPUT(CS11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <C> matrix') 
         CALL OUTPUT(CT11,1,NSPAIR,1,NAPAIR,NSPAIR,NSPAIR,1,LUPRI)
      END IF
      WRITE(CFN,1001) IPDD,LUMULX
      CALL GPOPEN(LU43,CFN,'UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LU43,'(4E30.20)') CS11
      CALL GPCLOSE(LU43,'KEEP')
      WRITE(CFN,1002) IPDD,LUMULX
      CALL GPOPEN(LU43,CFN,'UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LU43,'(4E30.20)') CT11
      CALL GPCLOSE(LU43,'KEEP')
C
      ELSE 
C
      CALL DZERO(FS,NSPAIR)
      CALL DZERO(FT,NSPAIR)
      CALL DZERO(CS11,NSPAIR*NSPAIR)
      CALL DZERO(CT11,NSPAIR*NSPAIR)
      F2S = D0
      F2T = D0
      IJ = 0
      JI = 0
      II = 0
      DO 511 ISYM = 1, NSYM
      DO 511 I = 1, NRHF(ISYM)
         II = II + 1
	 JJ = 0
         DO 512 JSYM = 1, NSYM
         DO 512 J = 1, NRHF(JSYM)
	    JJ = JJ + 1
            IF (JJ .GT. II) GOTO 512
            JI = JI + 1
            IF (I .GT. NRHFA(ISYM) .OR.
     *          J .GT. NRHFA(JSYM)) GOTO 512
            IJ = IJ + 1
            CNINV(IJ,1) = WS11(JI,IJ)
            CNINV(IJ,2) = BS11(JI,JI)
            CNINV(IJ,9) = -BS11(JI,JI) / RS11(JI,JI)
            IF (IPDD .EQ. 0) THEN
               CNINV(IJ,3) = 1D0
               CNINV(IJ,4) = 2D0*WS11(JI,IJ) - BS11(JI,JI)
            ELSE
               IF (CNINV(IJ,2) .LT. -VCLTHR) THEN
                  CNINV(IJ,3) = CNINV(IJ,1) / CNINV(IJ,2)
               ELSE
                  WRITE(LUPRI,'(/A,E15.8,I4/)')
     &          ' @ WARNING: NEGATIVE DIAGONAL ELEMENT OF SINGLET'//
     &          ' B MATRIX',CNINV(IJ,2),IJ
                  CNINV(IJ,3) = D0
               END IF
               CNINV(IJ,4) = CNINV(IJ,1) * CNINV(IJ,3)
            END IF
            IF (II .NE. JJ) THEN
               CNINV(IJ,5) = WT11(JI,IJ)
               CNINV(IJ,6) = BT11(JI,JI)
               CNINV(IJ,10) = -BT11(JI,JI) / RT11(JI,JI)
               IF (IPDD .EQ. 0) THEN
                   CNINV(IJ,7) = 0.5D0
                   CNINV(IJ,8) = WT11(JI,IJ) - 0.25D0*BT11(JI,JI)
               ELSE
                  IF (CNINV(IJ,6) .LT. -VCLTHR) THEN
                     CNINV(IJ,7) = CNINV(IJ,5) / CNINV(IJ,6)
                  ELSE
                     WRITE(LUPRI,'(/A,E15.8,I4/)')
     &             ' @ WARNING: NEGATIVE DIAGONAL ELEMENT OF TRIPLET'//
     &             ' B MATRIX',CNINV(IJ,6),IJ
                     CNINV(IJ,7) = D0
                  END IF
                  CNINV(IJ,8) = CNINV(IJ,5) * CNINV(IJ,7)
               END IF
            ELSE
               CNINV(IJ,5) = D0
               CNINV(IJ,6) = D0
               CNINV(IJ,7) = D0
               CNINV(IJ,8) = D0
               CNINV(IJ,10) = D0
            END IF
            FS(IJ) = CNINV(IJ,4)
            FT(IJ) = CNINV(IJ,8)
            IF (IPDD .GT. 1 .AND. IPDD .NE. 4 .AND. IPDD .NE. 6) THEN
              IF (R12SVD) THEN
                CALL R12MP2(WS11(1,IJ),WT11(1,IJ),WS11(1,IJ),WT11(1,IJ),
     *          RS11,RT11,EVS(JI),NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,
     *          QS11,QT11,BS11,BT11,EVS,XPDD,QQ11,IJ,WSMIN,WTMIN,D0,
     *                                            CS11(1,IJ),CT11(1,IJ))
              ELSE IF (R12DIA) THEN
                CALL MP2R12(WS11(1,IJ),WT11(1,IJ),WS11(1,IJ),WT11(1,IJ),
     *          RS11,RT11,EVS(JI),NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,
     *          QS11,QT11,BS11,BT11,EVS,XPDD,QQ11,IJ,WSMIN,WTMIN,D0,
     *                                            CS11(1,IJ),CT11(1,IJ))
              ELSE
                CALL QUIT('No solver selected (R12SVD or R12DIA)')
              END IF
              WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                          ET(IJ),ET(IJ)+FT(IJ),
     *                                          WSMIN,WTMIN
            ELSE
              WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                          ET(IJ),ET(IJ)+FT(IJ)
              CS11(IJ,IJ) = CNINV(IJ,3)
              CT11(IJ,IJ) = CNINV(IJ,7)
            END IF
            F2S = F2S + FS(IJ)
            F2T = F2T + FT(IJ)
  512    CONTINUE
  511 CONTINUE
C
      END IF
C
      WRITE(LUPRI,'(/4X,6F15.9)') E2S,E2S+F2S,E2T,E2T+F2T
c
c     Write out amplitudes for CC2-R12 model (HF/UniKA/25-06-2003).
c
      IF (CCR12) THEN
        IF (LUMULX .EQ. LUMULO) THEN
          WRITE(LU45) 0,IPDD,IPCC
        ELSE
          WRITE(LU45) 1,IPDD,IPCC
        END IF
        CALL DSCAL(NSPAIR*NSPAIR,0.5D0,CS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,0.5D0,CT11,1)
        CALL CCR12PCK(WORK(KSCR),1,CS11,CT11,NRHF,NRHFA,NKILJ)
        WRITE(LU45) (WORK(KSCR-1+I), I=1, NKILJ(1))
        ! undo scaling:
        CALL DSCAL(NSPAIR*NSPAIR,2.0D0,CS11,1)
        CALL DSCAL(NSPAIR*NSPAIR,2.0D0,CT11,1)
      END IF
c
      WRITE(CFN,1001) IPDD,LUMULX
      CALL GPOPEN(LU43,CFN,'UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LU43,'(4E30.20)') CS11
      CALL GPCLOSE(LU43,'KEEP')
      WRITE(CFN,1002) IPDD,LUMULX
      CALL GPOPEN(LU43,CFN,'UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LU43,'(4E30.20)') CT11
      CALL GPCLOSE(LU43,'KEEP')
c
      IF (IPDD .EQ. 0) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)')
     * ' Noninvariant MP2-R12/0   correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO I = 1, NACDIM
            DO J = 1, I
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,8E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7),CNINV(IJ,9),
     *                                 CNINV(IJ,10)
            END DO
         END DO
      ELSE IF (IPDD .EQ. 1) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)') 
     * ' Noninvariant MP2-R12/A   correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 600 I = 1, NACDIM
            DO 601 J = 1, I       
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,8E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7),CNINV(IJ,9),
     *                                 CNINV(IJ,10)
  601       CONTINUE
  600    CONTINUE
      ELSE IF (IPDD .EQ. 2) THEN
          IF (IPRINT .GE. 4) THEN
c
c    Print MP2-R12/A -Ansatz 1- amplitudes for 
c
            CALL AROUND('Singlet <C> matrix (Ansatz 1, MP2-R12/A)') 
            CALL OUTPUT(CS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
            CALL AROUND('Triplet <C> matrix (Ansatz 1, MP2-R12/A)') 
            CALL OUTPUT(CT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
          ENDIF
c
c    Write out amplitudes for MP2-R12-Gradient(Elena) 
c
         LUNIT = -1
         IF (R12PRP .AND. LUMULX .NE. LUMULO) THEN
           LUNIT = -1
           CALL GPOPEN(LUNIT,'CCR12_C','UNKNOWN',' ','FORMATTED',
     &                 IDUM,LDUM)
           WRITE(LUNIT,'(4E30.20)') CS11
           WRITE(LUNIT,'(4E30.20)') CT11
           CALL GPCLOSE(LUNIT,'KEEP')
C
           CALL CC_R12GRAD(IPDD,EV,V2AM,GRAD,LWORK)
           IF (FROIMP) THEN
              CALL CC_R12GRADFR(IPDD,EV,V2AM,GRAD,LWORK)
           END IF
         END IF
c
cElena  

         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/A   correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 3) THEN 

         LUNIT = -1
         IF (R12PRP .AND. LUMULX .NE. LUMULO) THEN
            IF (IPRINT .GE. 4) THEN
              CALL AROUND('Singlet <C> matrix (Ansatz 1, MP2-R12/A'')')
              CALL OUTPUT(CS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
              CALL AROUND('Triplet <C> matrix (Ansatz 1, MP2-R12/A'')')
              CALL OUTPUT(CT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
           END IF
           LUNIT = -1
           CALL GPOPEN(LUNIT,'CCR12_C','UNKNOWN',' ','FORMATTED',
     &                 IDUM,LDUM)
           WRITE(LUNIT,'(4E30.20)') CS11
           WRITE(LUNIT,'(4E30.20)') CT11
           CALL GPCLOSE(LUNIT,'KEEP')
c           KEND1 = 1
c           LWRK1 = LWORK - KEND1

           CALL CC_R12GRAD(IPDD,EV,V2AM,GRAD,LWORK)
           IF (FROIMP) THEN
              CALL CC_R12GRADFR(IPDD,EV,V2AM,GRAD,LWORK)
           END IF
         END IF
c
         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/A''  correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 4) THEN
         IF (R12HYB .OR. (LUMULX .EQ. LUMULO)) THEN
            WRITE(LUPRI,'(/A,F15.9//A/)') 
     *    ' Noninvariant MP2-R12/B   correlation energy =',
     *      E2S+E2T+F2S+F2T,
     *    '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *        '    V(IJ)       U(IJ)       C(IJ)   '
         ELSE
            WRITE(LUPRI,'(/A,F15.9//A/)') 
     *    ' Noninvariant MP2-R12/B''  correlation energy =',
     *      E2S+E2T+F2S+F2T,
     *    '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *        '    V(IJ)       U(IJ)       C(IJ)   '
         END IF
         IJ = 0
         DO 602 I = 1, NACDIM
            DO 603 J = 1, I       
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  603       CONTINUE
  602    CONTINUE
      ELSE IF (IPDD .EQ. 5) THEN
         IF (R12HYB .OR. (LUMULX .EQ. LUMULO)) THEN
            IF (IPRINT. GE. 4) THEN
              CALL AROUND('Singlet <C> matrix (Ansatz 1, MP2-R12/B)')
              CALL OUTPUT(CS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
              CALL AROUND('Triplet <C> matrix (Ansatz 1, MP2-R12/B)')
              CALL OUTPUT(CT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
            END IF
            LUNIT = -1
            IF (R12PRP .AND. LUMULX .NE. LUMULO) THEN
              LUNIT = -1
              CALL GPOPEN(LUNIT,'CCR12_C','UNKNOWN',' ','FORMATTED',
     &                    IDUM,LDUM)
              WRITE(LUNIT,'(4E30.20)') CS11
              WRITE(LUNIT,'(4E30.20)') CT11
              CALL GPCLOSE(LUNIT,'KEEP')
              CALL CC_R12GRAD(IPDD,EV,V2AM,GRAD,LWORK)
              IF (FROIMP) THEN
                  CALL CC_R12GRADFR(IPDD,EV,V2AM,GRAD,LWORK)
              END IF
            ENDIF
            WRITE(LUPRI,'(/A,F15.9 )')
     *    '              MP2-R12/B   correlation energy =',
     *      E2S+E2T+F2S+F2T    
         ELSE
            WRITE(LUPRI,'(/A,F15.9 )')
     *    '              MP2-R12/B''  correlation energy =',
     *      E2S+E2T+F2S+F2T    
         END IF
      ELSE IF (IPDD .EQ. 6) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)')
     *    ' Noninvariant MP2-R12/B   correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 604 I = 1, NACDIM
            DO 605 J = 1, I
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  605       CONTINUE
  604    CONTINUE
      ELSE IF (IPDD .EQ. 7) THEN
         WRITE(LUPRI,'(/A,F15.9 )')
     *    '              MP2-R12/B   correlation energy =',
     *   E2S+E2T+F2S+F2T  
      END IF
  999 CONTINUE
      IF (CCR12) THEN
        CALL GPCLOSE(LU44,'KEEP')
        CALL GPCLOSE(LU45,'KEEP')
      ENDIF
      RETURN
      END
C  /* Deck mp2r12 */
      SUBROUTINE MP2R12(VS,VT,VOS,VOT,SS,ST,EKL,NOCDIM,NSPAIR,FS,FT,
     *                  BB,UU,PP,QQ,SF,TF,EV,XF,RR,IJPAIR,
     *                  WS,WT,DELTA,CCS,CCT)
C
C     This subroutine evaluates the MP2-R12 energy by
C     solving sets of linear equations.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h> 
      PARAMETER (QNIL = 1D-12, QTHR = 1D-10, D0 = 0D0,
     *           D1 = 1D0, D2 = 2D0, DP5 = 0.5D0, D99 = 1D99)
      DIMENSION VS(NSPAIR), VT(NSPAIR), EV(NSPAIR)
      DIMENSION SS(NSPAIR,NSPAIR), ST(NSPAIR,NSPAIR)
      DIMENSION SF(NSPAIR,NSPAIR), TF(NSPAIR,NSPAIR)
      DIMENSION BB(*), PP(*), QQ(NSPAIR,NSPAIR), UU(NSPAIR,NSPAIR)
      DIMENSION RR(NSPAIR,NSPAIR)
      DIMENSION CCS(NSPAIR),CCT(NSPAIR)
      DIMENSION VOS(NSPAIR,NSPAIR),VOT(NSPAIR,NSPAIR)
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C

      N12 = NSPAIR * (NSPAIR + 1) / 2
      N2 = NSPAIR * NSPAIR
      WS = D99
      WT = D99
C
C     SINGLET PAIRS
C
      CALL DZERO(UU,N2)
      IJ = 0
      DO 100 I = 1, NSPAIR
         UU(I,I) = D1
         DO 101 J = 1, I         
            IJ = IJ + 1
            BB(IJ) = DP5 * (SS(I,J) + SS(J,I))
  101    CONTINUE
  100 CONTINUE
      CALL JACO(BB,UU,NSPAIR,NSPAIR,NSPAIR,PP,QQ)
      CALL DZERO(RR,N2)
      II = 0
      DO 102 I = 1, NSPAIR  
         II = II + I
         BBI = BB(II)
         IF (BBI .GT. QTHR) THEN
            BBI = D1 / SQRT(BBI)
         ELSE IF (ABS(BBI) .LT. QNIL) THEN
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         ELSE
            IF (IJPAIR .EQ. 1)
     &      WRITE(LUPRI,'(/A,E15.8/)')
     &       ' @ WARNING: SMALL EIGENVALUE OF SINGLET'//
     &       ' R12 OVERLAP',BBI
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         END IF
         DO 103 L = 1, NSPAIR
            BU = BBI * UU(L,I) 
            DO 104 K = 1, NSPAIR  
               RR(K,L) = RR(K,L) + UU(K,I) * BU
  104       CONTINUE
  103    CONTINUE
  102 CONTINUE
      IJ = 0
      DO 105 I = 1, NSPAIR
         DO 106 J = 1, I
            IJ = IJ + 1
            PP(IJ) = - SF(I,J) - SF(J,I) 
     *               - XF * (D2 * EKL - EV(I) - EV(J))
     *                    * (SS(I,J) + SS(J,I))
            PP(IJ) = PP(IJ) * DP5
  106    CONTINUE
         PP(IJ) = PP(IJ) + DELTA
  105 CONTINUE  
      CALL UTHU(PP,BB,RR,UU,NSPAIR,NSPAIR)
      CALL DZERO(UU,N2)
      DO 107 I = 1, NSPAIR   
         UU(I,I) = D1
  107 CONTINUE 
      CALL JACO(BB,UU,NSPAIR,NSPAIR,NSPAIR,PP,QQ)
      CALL DZERO(QQ,N2)
      II = 0
      DO 108 I = 1, NSPAIR  
         II = II + I
         BBI = BB(II)
         IF (BBI .GT. QTHR) THEN
            WS = MIN(WS,BBI)
            BBI = D1 / BBI
         ELSE IF (ABS(BBI) .LT. QNIL) THEN
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         ELSE
            WRITE(LUPRI,'(/A,E15.8,I4/)')
     &       ' @ WARNING: NEGATIVE EIGENVALUE OF SINGLET'//
     &       ' B MATRIX',BBI,IJPAIR
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         END IF
         DO 109 L = 1, NSPAIR
            BU = BBI * UU(L,I) 
            DO 110 K = 1, NSPAIR  
               QQ(K,L) = QQ(K,L) + UU(K,I) * BU
  110       CONTINUE
  109    CONTINUE
  108 CONTINUE
      CALL MPAB(RR,NSPAIR,NSPAIR,NSPAIR,NSPAIR,
     *          VS,NSPAIR,1,NSPAIR,1,
     *          PP,NSPAIR,1)
      CALL MPAB(RR,NSPAIR,NSPAIR,NSPAIR,NSPAIR,
     *          VOS,NSPAIR,1,NSPAIR,1,
     *          BB,NSPAIR,1)
      FS = D0
      DO 111 I = 1, NSPAIR
         UU(I,1) = D0
         DO 112 J = 1, NSPAIR
            UU(I,1) = UU(I,1) - QQ(I,J) * PP(J)
  112    CONTINUE
         FS = FS + BB(I) * UU(I,1)
  111 CONTINUE
      CALL MPAB(RR,NSPAIR,NSPAIR,NSPAIR,NSPAIR,
     *          UU,NSPAIR,1,NSPAIR,1,
     *          CCS,NSPAIR,1)
C
C     TRIPLET PAIRS
C
      CALL DZERO(UU,N2)
      IJ = 0
      DO 300 I = 1, NSPAIR
         UU(I,I) = D1
         DO 301 J = 1, I
            IJ = IJ + 1
            BB(IJ) = DP5 * (ST(I,J) + ST(J,I))
  301    CONTINUE
  300 CONTINUE
      CALL JACO(BB,UU,NSPAIR,NSPAIR,NSPAIR,PP,QQ)
      CALL DZERO(RR,N2)
      II = 0
      DO 302 I = 1, NSPAIR
         II = II + I
         BBI = BB(II)
         IF (BBI .GT. QTHR) THEN
            BBI = D1 / SQRT(BBI)
         ELSE IF (ABS(BBI) .LT. QNIL) THEN
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         ELSE
            IF (IJPAIR .EQ. 1)
     &      WRITE(LUPRI,'(/A,E15.8/)')
     &       ' @ WARNING: SMALL EIGENVALUE OF TRIPLET'//
     &       ' R12 OVERLAP',BBI
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         END IF
         DO 303 L = 1, NSPAIR
            BU = BBI * UU(L,I)
            DO 304 K = 1, NSPAIR
               RR(K,L) = RR(K,L) + UU(K,I) * BU
  304       CONTINUE
  303    CONTINUE
  302 CONTINUE
      IJ = 0
      DO 305 I = 1, NSPAIR
         DO 306 J = 1, I
            IJ = IJ + 1
            PP(IJ) = - TF(I,J) - TF(J,I) 
     *               - XF * (D2 * EKL - EV(I) - EV(J))
     *                    * (ST(I,J) + ST(J,I))
            PP(IJ) = PP(IJ) * DP5
  306    CONTINUE
         IF (ABS(PP(IJ)) .GT. 1D-12) PP(IJ) = PP(IJ) + DELTA
  305 CONTINUE
      CALL UTHU(PP,BB,RR,UU,NSPAIR,NSPAIR)
      CALL DZERO(UU,N2)
      DO 307 I = 1, NSPAIR
         UU(I,I) = D1
  307 CONTINUE
      CALL JACO(BB,UU,NSPAIR,NSPAIR,NSPAIR,PP,QQ)
      CALL DZERO(QQ,N2)
      II = 0
      DO 308 I = 1, NSPAIR
         II = II + I
         BBI = BB(II)
         IF (BBI .GT. QTHR) THEN
            WT = MIN(WT,BBI)
            BBI = D1 / BBI
         ELSE IF (ABS(BBI) .LT. QNIL) THEN
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         ELSE
            WRITE(LUPRI,'(/A,E15.8,I4/)')
     &       ' @ WARNING: NEGATIVE EIGENVALUE OF TRIPLET'//
     &       ' B MATRIX',BBI,IJPAIR
            BBI = D0
            CALL DZERO(UU(1,I),NSPAIR)
         END IF
         DO 309 L = 1, NSPAIR
            BU = BBI * UU(L,I)
            DO 310 K = 1, NSPAIR
               QQ(K,L) = QQ(K,L) + UU(K,I) * BU
  310       CONTINUE
  309    CONTINUE
  308 CONTINUE
      CALL MPAB(RR,NSPAIR,NSPAIR,NSPAIR,NSPAIR,
     *          VT,NSPAIR,1,NSPAIR,1,
     *          PP,NSPAIR,1)
      CALL MPAB(RR,NSPAIR,NSPAIR,NSPAIR,NSPAIR,
     *          VOT,NSPAIR,1,NSPAIR,1,
     *          BB,NSPAIR,1)
      FT = D0
      DO 311 I = 1, NSPAIR
         UU(I,1) = D0
         DO 312 J = 1, NSPAIR
            UU(I,1) = UU(I,1) - QQ(I,J) * PP(J)
  312    CONTINUE
         FT = FT + BB(I) * UU(I,1)
  311 CONTINUE 
      CALL MPAB(RR,NSPAIR,NSPAIR,NSPAIR,NSPAIR,
     *          UU,NSPAIR,1,NSPAIR,1,
     *          CCT,NSPAIR,1)
      RETURN
      END
C  /* Deck iq */
      INTEGER FUNCTION IQ(I)
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h> 
      J = 0
      K = 0
      DO 100 ISYM = 1, NSYM
         K = K + NRHFFR(ISYM)
         DO 101 L = 1, NRHF(ISYM)
            J = J + 1
            K = K + 1
            IF (J .EQ. I) THEN
               IQ = K
               GOTO 99
            END IF
  101    CONTINUE
  100 CONTINUE
   99 RETURN
      END
C  /* Deck testuu */
      SUBROUTINE TESTUU(U2AM,R2AM)         
C
C     This subroutine generates on R2AM an array (ip|[T1+T2,r12)|jq)
C     with (ia).le.(jb) from the asymmetric integrals (ip|[T1,r12)|jq)
C     on U2AM without this triangular constraint.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION U2AM(*), R2AM(*)
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <r12int.h>
      DO 101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         NAIBJ = IT2AM(ISYMIA,ISYMJB)
         KOFFU = IU2AM(ISYMIA,ISYMJB)
         NTOTAI = NT1AM(ISYMIA) 
         DO 301 IA = 1, NTOTAI
            DO 302 JB = 1, IA
               NAIBJ = NAIBJ + 1
               R2AM(NAIBJ) = U2AM(KOFFU + (IA-1)*NTOTAI + JB) +
     &                       U2AM(KOFFU + (JB-1)*NTOTAI + IA)
  302        CONTINUE
  301    CONTINUE
  101 CONTINUE
      RETURN
      END
C  /* Deck makeug */
      SUBROUTINE MAKEUG(GX,UG,QQ,QV,WORK,LWORK,IPBAS)
C
C     This subroutine computes the orbitals i' by means
C     of the matrix product |i'> = Sum_p' X(p'i) |p'>.
C
C     It then calls R12QQ to compute the integrals
C     (i'k|r12**2|jl), which factorize into products
C     of one-particle matrices.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION GX(*), UG(*), WORK(*), QQ(*),QV(*)
      CHARACTER*7 GUNAM
      INTEGER IDUM
      LOGICAL LDUM
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J    
C
celena
      NORB1T = 0
      NORB2T = 0
      NRHFST = 0
      DO ISYM = 1, NSYM
         NORB1T = NORB1T + NORB1(ISYM)
         NORB2T = NORB2T + NORB2(ISYM)
         NRHFST = NRHFST + NRHFS(ISYM)
      ENDDO   
      LUSIRG = -1
      CALL GPOPEN(LUSIRG,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,LDUM)
      REWIND LUSIRG

      CALL MOLLAB(LABEL,LUSIRG,LUPRI)
      READ (LUSIRG)
C
      READ (LUSIRG)
      READ (LUSIRG) (WORK(i), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIRG,'KEEP')


      
celena

      NOCCT = NRHFT
      LUSIRG = -30
      KEND = NLAMDS + 2*NOCCT*NBAST + 14*NOCCT*NOCCT + 1
      IF (KEND .GT. LWORK) CALL STOPIT('MAKEUG',' ',LWORK,KEND)
      CALL GPOPEN(LUSIRG,'SIRIFC','OLD',' ','UNFORMATTED',IDUM,LDUM)
      REWIND LUSIRG

      CALL MOLLAB(LABEL,LUSIRG,LUPRI)
      READ (LUSIRG)
C
      READ (LUSIRG)
      READ (LUSIRG) (WORK(I), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIRG,'KEEP')  
C
      IWORK = 1
      IGX = 1 
      IUG = 1 
      IF (R12PRP .AND. R12HYB) THEN
         DO ISYM = 1, NSYM
            IWORK = IWORK + NBAS(ISYM)*NRHFS(ISYM) 
            CALL MPAB(WORK(IWORK),NBAS(ISYM),NVIR(ISYM),
     *                            NBAS(ISYM),NVIR(ISYM),
     *                GX(IGX),NVIR(ISYM),NORB1(ISYM),
     *                        NVIR(ISYM),NORB1(ISYM),
     *                UG(IUG),NBAS(ISYM),NORB1(ISYM))
            IF (IPRINT .GT. 3) THEN
               IF (NORB1(ISYM)*NBAS(ISYM) .NE. 0) 
     *         WRITE(LUPRI,'(/A,I2)') ' U(NORB1) x X/ISYM =',ISYM
               CALL OUTPUT(UG(IUG),1,NBAS(ISYM),1,NORB1(ISYM),
     *                     NBAS(ISYM),NORB1(ISYM),1,LUPRI)  
               WRITE(LUPRI,'(/A,I2)') ' G(NORB1) x X/ISYM =',ISYM
               CALL OUTPUT(GX(IGX),1,NBAS(ISYM),1,NORB1(ISYM),
     *                     NBAS(ISYM),NORB1(ISYM),1,LUPRI)  
            ENDIF
            IWORK = IWORK + NBAS(ISYM)*NVIRS(ISYM)
            IGX = IGX + NVIR(ISYM)*(NORB1(ISYM)-NRHFFR(ISYM))
            IUG = IUG + NBAS(ISYM)*(NORB1(ISYM)-NRHFFR(ISYM))
         ENDDO   
      ELSE 
         DO 100 ISYM = 1, NSYM
            IWORK = IWORK + NBAS(ISYM)*NRHFS(ISYM)
            CALL MPAB(WORK(IWORK),NBAS(ISYM),NVIR(ISYM),
     *                            NBAS(ISYM),NVIR(ISYM),
     *                GX(IGX),NVIR(ISYM),NRHF(ISYM),
     *                        NVIR(ISYM),NRHF(ISYM),
     *                UG(IUG),NBAS(ISYM),NRHF(ISYM))
            IF (IPRINT .GT. 3) THEN
               IF (NRHF(ISYM)*NBAS(ISYM) .NE. 0) 
     *         WRITE(LUPRI,'(/A,I2)') ' U x X/ISYM =',ISYM
               CALL OUTPUT(UG(IUG),1,NBAS(ISYM),1,NRHF(ISYM),
     *                     NBAS(ISYM),NRHF(ISYM),1,LUPRI)  
             ENDIF 
            IWORK = IWORK + NBAS(ISYM)*NVIRS(ISYM)
            IGX = IGX + NVIR(ISYM)*NRHF(ISYM)
            IUG = IUG + NBAS(ISYM)*NRHF(ISYM)
  100    CONTINUE
      ENDIF 
      NTOTGU = IUG - 1
C
C     Add contribution from natural virtual orbitals (WK/UniKA/22-11-2005).
C
      IF (NATVIR) THEN
         NR12xVIR = 0
         DO ISYM = 1, NSYM
           NR12xVIR = NR12xVIR + (NORB1(ISYM)-NRHFSA(ISYM))*NRXR12(ISYM)
         END DO
         KFCKMIX = KEND
         KEND = KFCKMIX + NR12xVIR
         KFCKHF = KEND + NR12xVIR * NBAST
         KEND = KFCKHF
         IF (KEND .GT. LWORK) CALL STOPIT('MAKEUG',' ',LWORK,KEND)
         LUNIT = -1
         CALL GPOPEN(LUNIT,'R12FVIR','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
         READ (LUNIT) (WORK(KFCKMIX+I-1), I=1, NR12xVIR)
         CALL GPCLOSE(LUNIT,'KEEP')
         XALPH = 1.0D0
         IF (R12NOB) THEN
            XBETA = 0.0D0
         ELSE
            XBETA = 1.0D0
         END IF
         IUG = 1
         JUG = 1
         IWORK = 1
         DO ISYM = 1, NSYM
           IUG = IUG + NBAS(ISYM)*NRHFA(ISYM)
           IWORK = IWORK + NBAS(ISYM)*(NRHFS(ISYM)+NRHFSA(ISYM))
           NREAVR = NORB1(ISYM)-NRHFSA(ISYM)
           IF (NRXR12(ISYM).GT.0)
     *        CALL DGEMM('N','N',NBAS(ISYM),NRXR12(ISYM),NREAVR,XALPH,
     *                    WORK(IWORK),NBAS(ISYM),WORK(KFCKMIX),NREAVR,
     *                    XBETA,UG(IUG),NBAS(ISYM))
           IF (IPRINT .GT. 4) THEN
              IF (NRXR12(ISYM)*NBAS(ISYM) .NE. 0) THEN
                 WRITE(LUPRI,'(/A,I2)') ' MO MATRIX/ISYM =',ISYM
                 CALL OUTPUT(WORK(IWORK),1,NBAS(ISYM),1,NREAVR,
     *                       NBAS(ISYM),NREAVR,1,LUPRI)  
                 WRITE(LUPRI,'(/A,I2)') ' FOCK MATRIX/ISYM =',ISYM
                 CALL OUTPUT(WORK(KFCKMIX),1,NREAVR,1,NRXR12(ISYM),
     *                       NREAVR,NRXR12(ISYM),1,LUPRI)  
              ENDIF 
              WRITE(LUPRI,'(/A,I2)') ' U x (X + F)/ISYM =',ISYM
              CALL OUTPUT(UG(JUG),1,NBAS(ISYM),1,NRHF(ISYM),
     *                    NBAS(ISYM),NRHF(ISYM),1,LUPRI)  
              JUG = JUG + NBAS(ISYM)*NRHF(ISYM)
           ENDIF 
           IUG = IUG + NBAS(ISYM)*NRXR12(ISYM)
           IWORK = IWORK + NBAS(ISYM)*(NREAVR+NORB2(ISYM))
           KFCKMIX = KFCKMIX + NREAVR*NRXR12(ISYM)
         END DO
         IF (IWORK-1 .NE. NLAMDS) 
     *   CALL STOPIT('MAKEUGa',' ',IWORK-1,NLAMDS)
         IF (NTOTGU .NE. IUG-1) 
     *   CALL STOPIT('MAKEUGb',' ',NTOTGU,IUG-1)
      END IF
C 
      NTOTGU = IUG - 1
      WRITE(GUNAM,'(A6,I1)') 'GUMAT.',IPBAS
      LU43 = -43
      CALL GPOPEN(LU43,GUNAM,'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
      REWIND(LU43)
      WRITE(LU43) NTOTGU
      WRITE(LU43) (UG(I), I = 1, NTOTGU)
      CALL GPCLOSE(LU43,'KEEP') 
C
      ISMOU = 1              
      ISMO0 = 1              
      ISMO1 = ISMO0 + NLAMDS
celena
      ISMO2 = ISMO1 + NORB1T*NBAST
celena
      CALL DZERO(WORK(ISMO1),NORB1T*NBAST)
      CALL DZERO(WORK(ISMO2),NORB1T*NBAST)
      DO 505 ISYM=1,NSYM
         ISMO0 = ISMO0 + NBAS(ISYM)*NRHFFR(ISYM)
         DO 510 I = 1, NRHF(ISYM) 
           CALL DCOPY(NBAS(ISYM),WORK(ISMO0),1,WORK(ISMO1),1)
           CALL DCOPY(NBAS(ISYM),UG(ISMOU),1,WORK(ISMO2),1)
           ISMO1 = ISMO1 + NBAST
           ISMO2 = ISMO2 + NBAST
           ISMO0 = ISMO0 + NBAS(ISYM)
           ISMOU = ISMOU + NBAS(ISYM)
  510    CONTINUE
         ISMO0 = ISMO0 + NVIR(ISYM)*NBAS(ISYM) 
         IF (R12PRP) ISMOU=ISMOU + (NORB1(ISYM)-NRHFS(ISYM))*NBAS(ISYM) 
         ISMO1 = ISMO1 + NBAS(ISYM)
         ISMO2 = ISMO2 + NBAS(ISYM)
  505 CONTINUE 
celena
      IF (R12PRP) THEN
         ISMO5 = 1 
         ISMO3 = 1 + NLAMDS + NOCCT*NBAST
         ISMO4 = 1 
         ISMO6 = ISMO3 + NORB1T*NBAST
         CALL DZERO(WORK(ISMO3),(NORB1T-NRHFST)*NBAST)
         DO  ISYM=1,NSYM
             ISMO4 = ISMO4 + 2*NBAS(ISYM)*NRHFS(ISYM)
             ISMO5 = ISMO5 + NBAS(ISYM)*NRHF(ISYM)
             DO I = 1, NORB1(ISYM)-NRHFS(ISYM)
               CALL DCOPY(NBAS(ISYM),WORK(ISMO4),1,WORK(ISMO3),1)
               CALL DCOPY(NBAS(ISYM),UG(ISMO5),1,WORK(ISMO6),1)
               ISMO3 = ISMO3 + NBAST
               ISMO4 = ISMO4 + NBAS(ISYM)
               ISMO6 = ISMO6 + NBAST
               ISMO5 = ISMO5 + NBAS(ISYM)
             ENDDO 
             ISMO3 = ISMO3 + NBAS(ISYM)
             ISMO6 = ISMO6 + NBAS(ISYM)
             ISMO4 = ISMO4 + NORB2(ISYM)*NBAS(ISYM)
          ENDDO    
      ENDIF
celena

      ISMO0 = 1              
      ISMA1 = ISMO0 + 2*NOCCT*NBAST 
      ISMO1 = ISMO0 + NLAMDS
      ISMO2 = ISMO1 + NORB1T*NBAST
      ISMS0 = ISMO2 + NORB1T*NBAST
      ISMX1 = ISMS0 + NOCCT*NOCCT*2
      ISMY1 = ISMX1 + NOCCT*NOCCT*2
      ISMZ1 = ISMY1 + NOCCT*NOCCT*2
      ISMX2 = ISMZ1 + NOCCT*NOCCT*2
      ISMY2 = ISMX2 + NOCCT*NOCCT*2
      ISMZ2 = ISMY2 + NOCCT*NOCCT*2
      IEND  = ISMZ2 + NOCCT*NOCCT*2
celena
      IF (R12PRP .AND. R12HYB) THEN
         ISMRO0  = 1 
         ISMRA1  = ISMRO0  + 2*NOCCT*NBAST
         ISMRO1  = ISMRO0  + NLAMDS
         ISMRO2  = ISMRO1  + NORB1T*NBAST
         ISMRS0  = ISMRO2  + NORB1T*NBAST
         ISMRX1  = ISMRS0  + 2*NORB1T*NORB1T
         ISMRY1  = ISMRX1  + 2*NORB1T*NORB1T
         ISMRZ1  = ISMRY1  + 2*NORB1T*NORB1T
         ISMRX2  = ISMRZ1  + 2*NORB1T*NORB1T
         ISMRY2  = ISMRX2  + 2*NORB1T*NORB1T
         ISMRZ2  = ISMRY2  + 2*NORB1T*NORB1T
         IEND    = ISMRZ2  + 2*NORB1T*NORB1T
      ENDIF
celena
      LQQ   = LWORK - IEND
      IF (IPRINT .GT. 3) THEN
         CALL AROUND('Active orbitals')
         CALL OUTPUT(WORK(ISMO1),1,NBAST,1,NOCCT,NBAST,NOCCT,1,LUPRI)
         CALL AROUND('Primed active orbitals')
         CALL OUTPUT(WORK(ISMO2),1,NBAST,1,NOCCT,NBAST,NOCCT,1,LUPRI)
      END IF
      IF (LQQ .LT. 0) CALL STOPIT('MAKEUG','R12QQ',LWORK,IEND)
      CALL R12QQ(QQ,WORK(ISMS0),WORK(ISMX1),WORK(ISMY1),WORK(ISMZ1), 
     *           WORK(ISMX2),WORK(ISMY2),WORK(ISMZ2), 
     *           WORK(ISMO2),WORK(ISMO1),
     *           WORK(IEND),LQQ,NBAST,NOCCT)
        LUNIT = -1
celena
      IF (R12PRP .AND. R12HYB) THEN
         NVIR0T = NORB1T-NRHFST

         CALL DZERO(WORK(ISMRS0),14*NORB1T*NORB1T) 
         CALL DZERO(QV,NVIR0T*NRHFT**3)
         CALL R12QV(.FALSE.,.FALSE.,QV,WORK(ISMRS0),WORK(ISMRX1),
     *              WORK(ISMRY1),
     *              WORK(ISMRZ1),
     *              WORK(ISMRX2),WORK(ISMRY2),
     *              WORK(ISMRZ2), 
     *              WORK(ISMRO2),WORK(ISMRO1),
     *              WORK(IEND),LQQ,NBAST,NORB1T,NOCCT,NVIR0T,.FALSE.)
         LUNIT = -1
         CALL GPOPEN(LUNIT,'AUXQSA12','UNKNOWN','SEQUENTIAL',
     &                      'FORMATTED',IDUMMY,LDUMMY)
         WRITE(LUNIT,'(4E30.20)') (QV(J+1), J = 0, (NORB1T-NRHFST)
     &                                     *NOCCT**3 - 1)
         CALL GPCLOSE(LUNIT,'KEEP')

         CALL DZERO(WORK(ISMRS0),14*NORB1T*NORB1T) 
         CALL DZERO(QV,NVIR0T*NRHFT**3)
         CALL R12QV(.FALSE.,.FALSE.,QV,WORK(ISMRS0),WORK(ISMRX1),
     *              WORK(ISMRY1),
     *              WORK(ISMRZ1),
     *              WORK(ISMRX2),WORK(ISMRY2),
     *              WORK(ISMRZ2), 
     *              WORK(ISMRO2),WORK(ISMRO1),
     *              WORK(IEND),LQQ,NBAST,NORB1T,NOCCT,NVIR0T,.TRUE.)
         LUNIT = -1
         CALL GPOPEN(LUNIT,'AUXQSAJ12','UNKNOWN','SEQUENTIAL',
     &                      'FORMATTED',IDUMMY,LDUMMY)
         WRITE(LUNIT,'(4E30.20)') (QV(J+1), J = 0, (NVIR0T)
     &                                     *NRHFT**3 - 1)
         CALL GPCLOSE(LUNIT,'KEEP')

         CALL DZERO(WORK(ISMRS0),14*NORB1T*NORB1T) 
         CALL DZERO(QV,NVIR0T*NRHFT**3)
         CALL R12QV(.FALSE.,.TRUE.,QV,WORK(ISMRS0),WORK(ISMRX1),
     *              WORK(ISMRY1),
     *              WORK(ISMRZ1),
     *              WORK(ISMRX2),WORK(ISMRY2),
     *              WORK(ISMRZ2), 
     *              WORK(ISMRO2),WORK(ISMRO1),
     *              WORK(IEND),LQQ,NBAST,NORB1T,NOCCT,NVIR0T,.FALSE.)
         LUNIT = -1
         CALL GPOPEN(LUNIT,'AUXQSAI12','UNKNOWN','SEQUENTIAL',
     &                      'FORMATTED',IDUMMY,LDUMMY)
         WRITE(LUNIT,'(4E30.20)') (QV(J+1), J = 0, (NVIR0T)
     &                                     *NRHFT**3 - 1)
         CALL GPCLOSE(LUNIT,'KEEP')

      ENDIF  
celena
      RETURN
      END
C  /* Deck r12qq */  
      SUBROUTINE R12QQ(QQ,S0,X1,Y1,Z1,X2,Y2,Z2,UU,ZZ,
     *                 WORK,LWORK,NBAST,NOCCT)
C
C     This subroutine computes the integrals (i'k|r12**2|jl).
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION QQ(*),X1(*),Y1(*),Z1(*),X2(*),Y2(*),Z2(*),
     *          UU(*),ZZ(*),WORK(LWORK),S0(*)
      LOGICAL FOUND
      KEND = 1 + NBAST*NBAST
      LWRK = LWORK - KEND
      CALL RDPROP('CM000000',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,S0,WORK(KEND),LWRK,NBAST,NOCCT)
      CALL RDPROP('CM010000',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,X1,WORK(KEND),LWRK,NBAST,NOCCT)
      CALL RDPROP('CM000100',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Y1,WORK(KEND),LWRK,NBAST,NOCCT)
      CALL RDPROP('CM000001',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Z1,WORK(KEND),LWRK,NBAST,NOCCT)
      CALL RDPROP('CM020000',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,X2,WORK(KEND),LWRK,NBAST,NOCCT)
      CALL RDPROP('CM000200',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Y2,WORK(KEND),LWRK,NBAST,NOCCT)
      CALL RDPROP('CM000002',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Z2,WORK(KEND),LWRK,NBAST,NOCCT)
      CALL MAKEQQ(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT)
      RETURN
      END
C  /* Deck uthz */  
      SUBROUTINE UTHZ(AA,UU,ZZ,BB,WORK,LWORK,NBAST,NOCCT) 
C
C              B(1) = Z(Transpose) * A * Z
C              B(2) = U(Transpose) * A * Z
C   
C     On input, AA is the upper triangle of the symmetric 
C     matrix A. B is returned through the array BB as 
C     square matrix.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION AA(*),UU(NBAST,NOCCT),ZZ(NBAST,NOCCT),
     *          WORK(NBAST,NOCCT),BB(NOCCT,NOCCT,2)
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J    
      IF (NBAST*NOCCT .GT. LWORK) 
     *      CALL STOPIT('UTHZ',' ',LWORK,NBAST*NOCCT)
      CALL DZERO(WORK,NBAST*NOCCT)
      DO 103 K = 1, NBAST 
         DO 102 J = 1, NOCCT
            ZZKJ = ZZ(K,J)
            DO 101 I = 1, NBAST
               IK = INDEX(I,K)
               WORK(I,J) = WORK(I,J) + AA(IK) * ZZKJ
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 202 I = 1, NOCCT
         DO 201 J = 1, NOCCT     
            BB(I,J,1) = DDOT(NBAST,ZZ(1,I),1,WORK(1,J),1)
            BB(I,J,2) = DDOT(NBAST,UU(1,I),1,WORK(1,J),1)
  201    CONTINUE
  202 CONTINUE
      RETURN
      END
C  /* Deck makeqq */
      SUBROUTINE MAKEQQ(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT)
C
C     This subroutine computes the integrals (i'k|r12**2|jl). 
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NOCCT,NOCCT,2),Y1(NOCCT,NOCCT,2),
     *          Z1(NOCCT,NOCCT,2),X2(NOCCT,NOCCT,2),
     *          Y2(NOCCT,NOCCT,2),Z2(NOCCT,NOCCT,2),
     *          S0(NOCCT,NOCCT,2),
     *          QQ(NOCCT,NOCCT,NOCCT,NOCCT) 
      DO 100 I = 1, NOCCT
        DO 200 J = 1, NOCCT
          DO 300 K = 1, NOCCT
            DO 400 L = 1, NOCCT
              QQ(I,J,K,L) =
     *        (X2(I,K,2) + Y2(I,K,2) + Z2(I,K,2)) * S0(J,L,1) +
     *        (X2(J,L,1) + Y2(J,L,1) + Z2(J,L,1)) * S0(I,K,2) -
     *        D2 * ( X1(I,K,2) * X1(J,L,1) +
     *               Y1(I,K,2) * Y1(J,L,1) +
     *               Z1(I,K,2) * Z1(J,L,1) )
 400        CONTINUE
 300      CONTINUE
 200    CONTINUE
 100  CONTINUE
      RETURN
      END 
C  /* Deck qcal */
      SUBROUTINE QCAL(QQ,R2AM,U2AM,WORK,LWORK,IFLAG,IANSAZ)
C
C     Driver routine for the computation of the Q term
C     of R12/B theory.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      DIMENSION R2AM(*), U2AM(*), WORK(*), QQ(*)
C
      NOCDIM = NRHFT
      NSPAIR = NOCDIM * (NOCDIM + 1) / 2 
C
      IR11  = 1                  
      IQ11  = IR11  + 2 * NOCDIM**4   
      IRS11 = IQ11  + 2 * NOCDIM**4 
      IQS11 = IRS11 + NSPAIR**2
      IRT11 = IQS11 + NSPAIR**2
      IQT11 = IRT11 + NSPAIR**2
      IEND  = IQT11 + NSPAIR**2 
      IF (IEND .GT. LWORK)
     *CALL STOPIT('QCAL',' ',LWORK,IEND) 
      CALL QCAL1(R2AM,U2AM,QQ,WORK(IR11),WORK(IQ11), 
     *           WORK(IRS11),WORK(IQS11),WORK(IRT11),
     *           WORK(IQT11),NOCDIM,NSPAIR,IFLAG,IANSAZ)
      RETURN
      END
C  /* Deck qcal1 */
      SUBROUTINE QCAL1(R2AM,U2AM,QQ11,R11,Q11,RS11,RT11,QS11,QT11,
     *                 NOCDIM,NSPAIR,IFLAG,IANSAZ)
C
C     This subroutine computes the Q term of R12/B:
C
C     (i'k|r12**2|jl) = Sum_pq (i'p'|r12|jq')(kp'|r12|lq')
C
C     Note that i' refers to the orbital |i'> = Sum_p' X(p'i) |p'>
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      DIMENSION R2AM(*), U2AM(*)
      REAL*8  Q11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        R11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        U1, R2, RR, QQIJKL, QQIJLK, UIJKL, UJIKL, UIJLK, UJILK
      DIMENSION QS11(NSPAIR,NSPAIR), RS11(NSPAIR,NSPAIR)
      DIMENSION QT11(NSPAIR,NSPAIR),RT11(NSPAIR,NSPAIR)
      DIMENSION QQ11(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
      INTEGER IDUM
      LOGICAL LDUM
#include <r12int.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      IF (R12NOB) THEN
         CALL DZERO(RS11,NSPAIR**2)
         CALL DZERO(RT11,NSPAIR**2)
         CALL DZERO(QS11,NSPAIR**2)
         CALL DZERO(QT11,NSPAIR**2)
         GOTO 998
      ENDIF
C
      FF = D1 / SQRT(D2)
      DO 90 I = 1, NOCDIM**4
         R11(I,1,1,1) = D0
         Q11(I,1,1,1) = D0
   90 CONTINUE
C
      IF (R12CBS) THEN
         FACX = - D1
      ELSE
         FACX = D1 
      END IF
C
C     CONSTRUCT PRODUCTS (IA|JB)*(KA|LB)
C
      DO 101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         JBOFF = NH1AM(ISYMJB) * (NH1AM(ISYMJB) + 1) / 2
         DO 102 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 103 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
               DO 104 ISYMKA = 1, NSYM
                  ISYMLB = ISYMKA 
                  LBOFF = NH1AM(ISYMLB) * (NH1AM(ISYMLB) + 1) / 2
                  ISYMK = MULD2H(ISYMA,ISYMKA)
                  ISYML = MULD2H(ISYMB,ISYMLB)
                  DO 105 I = 1, NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     DO 106 K = 1, NRHF(ISYMK)
                        KOFFK = IRHF(ISYMK) + K
                        DO 107 A = 1, NVIR(ISYMA)
                           IF (ONEAUX) THEN
                              IF (A .LE. NORB1(ISYMA)) THEN
                                 NAI = IH1AM(ISYMA,ISYMI) +
     *                                 NORB1(ISYMA)*(I-1) + A
                                 NAK = IH1AM(ISYMA,ISYMK) +
     *                                 NORB1(ISYMA)*(K-1) + A
                              ELSE
                                 NAI = IG1AM(ISYMA,ISYMI) +
     *                                 NORB2(ISYMA)*(I-1) + A -
     *                                 NORB1(ISYMA)
                                 NAK = IG1AM(ISYMA,ISYMK) +
     *                                 NORB2(ISYMA)*(K-1) + A -
     *                                 NORB1(ISYMA)
                              END IF
                           ELSE
                              NAI = IT1AM(ISYMA,ISYMI) +
     *                              NVIR(ISYMA)*(I-1) + A
                              NAK = IT1AM(ISYMA,ISYMK) +
     *                              NVIR(ISYMA)*(K-1) + A
                           END IF
                           DO 108 J = 1, NRHF(ISYMJ)
                              KOFFJ = IRHF(ISYMJ) + J
                              DO 109 L = 1, NRHF(ISYML)
                                 KOFFL = IRHF(ISYML) + L
                                 IF (ONEAUX) THEN
                                    NENDB = NORB1(ISYMB)
                                 ELSE
                                    NENDB = NVIR(ISYMB)
                                 END IF
                                 DO 110 B = 1, NENDB
                                    IF (ONEAUX) THEN
                                       NBJ = IH1AM(ISYMB,ISYMJ) +
     *                                       NORB1(ISYMB)*(J-1) + B
                                       NBL = IH1AM(ISYMB,ISYML) +
     *                                       NORB1(ISYMB)*(L-1) + B
                                       IF (A .LE. NORB1(ISYMA)) THEN
                                          NAIBJ = IH2AM(ISYMIA,ISYMJB)
     *                                          + INDEX(NAI,NBJ)
                                          NAKBL = IH2AM(ISYMKA,ISYMLB)
     *                                          + INDEX(NAK,NBL)
                                       ELSE
                                          NAIBJ = IH2AM(ISYMIA,ISYMJB) 
     *                                          + NH1AM(ISYMJB)*(NAI-1)
     *                                          + NBJ + JBOFF
                                          NAKBL = IH2AM(ISYMKA,ISYMLB) 
     *                                          + NH1AM(ISYMLB)*(NAK-1) 
     *                                          + NBL + LBOFF
                                       END IF
                                    ELSE
                                       NBJ = IT1AM(ISYMB,ISYMJ) +
     *                                       NVIR(ISYMB)*(J-1) + B
                                       NBL = IT1AM(ISYMB,ISYML) +
     *                                       NVIR(ISYMB)*(L-1) + B
                                       NAIBJ = IU2AM(ISYMIA,ISYMJB) 
     *                                       + NT1AM(ISYMJB)*(NAI-1)
     *                                       + NBJ
                                       NAKBL = IT2AM(ISYMKA,ISYMLB) +
     *                                         INDEX(NAK,NBL)
                                    END IF
                                    U1 = U2AM(NAIBJ)
                                    R2 = R2AM(NAKBL)
                                    RR = U1 * R2
C      
                                    GOTO (113,114), IFLAG
                                    CALL QUIT('INVALID IFLAG IN QCAL1')
C
  113                               CONTINUE
                                    IF (A .GT. NORB1(ISYMA)) GOTO 112
                                    IF (B .GT. NORB1(ISYMB)) GOTO 112
C 
                                    Q11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
C
                                    GOTO 111
  112                               CONTINUE
C       
                                    IF (A .GT. NORB1(ISYMA) .AND.
     *                                  B .GT. NORB1(ISYMB)) GOTO 111
C
                                    R11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              R11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
                                    IF (ONEAUX) THEN
                                       R11(KOFFJ,KOFFI,KOFFL,KOFFK) = 
     *                                 R11(KOFFJ,KOFFI,KOFFL,KOFFK) + RR
                                    END IF
C
                                    GOTO 111
C
  114                               CONTINUE
                                    IF (R12CBS) THEN
                                       NQUEA = 0
                                       NQUEB = 0
                                    ELSE
                                       NQUEA = NORB1(ISYMA)
                                       NQUEB = NORB1(ISYMB)
                                    END IF
                                    IF (A .LE. NRHFSA(ISYMA) .AND.
     *                                  B .LE. NRHFSA(ISYMB)) THEN   
C
C                                      Sum |ij><ij|
C
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) - RR
C
                                    END IF
                                    IF (A .LE. NRHFSA(ISYMA) .AND.
     *                                  B .GT. NQUEB) THEN
C
C                                      Sum |iq'><iq'|
C
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
C
                                    END IF
                                    IF (A .GT. NQUEA .AND.
     *                                  B .LE. NRHFSA(ISYMB)) THEN
C
C                                      Sum |p'j><p'j|
C
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
                                       IF (ONEAUX.AND.A.GT.NORB1(ISYMA))
     *                                 Q11(KOFFJ,KOFFI,KOFFL,KOFFK) =
     *                                 Q11(KOFFJ,KOFFI,KOFFL,KOFFK) + RR
C                                       
                                    END IF
                                    IF (A .GT. NRHFSA(ISYMA) .AND.
     *                                  A .LE. NORB1(ISYMA) .AND.
     *                                  B .GT. NRHFSA(ISYMB) .AND.
     *                                  B .LE. NORB1(ISYMB)) THEN
C
C                                      Sum |ab><ab|
C
cWK                                    IF (IANSAZ .EQ. 3)
cWK  *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
cWK  *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
C
                                    END IF
  111                               CONTINUE
  110                            CONTINUE
  109                         CONTINUE
  108                      CONTINUE
  107                   CONTINUE
  106                CONTINUE
  105             CONTINUE
  104          CONTINUE
  103       CONTINUE
  102    CONTINUE
  101 CONTINUE
C
      IJ = 0
      DO 300 I = 1, NOCDIM
         DO 301 J = 1, I
            IJ = IJ + 1
            KL = 0
            DO 302 K = 1, NOCDIM
               DO 303 L = 1, K
                  KL = KL + 1
C
                  UIJKL = QQ11(I,J,K,L)
                  UJILK = QQ11(J,I,L,K)
                  UIJLK = QQ11(I,J,L,K)
                  UJIKL = QQ11(J,I,K,L)
C
                  QQIJKL = UIJKL + UJILK
                  QQIJLK = UIJLK + UJIKL
                  
                  RR =   (QQIJKL + QQIJLK) 
     *                 - (Q11(I,J,K,L) + Q11(I,J,L,K))
                  IF (.NOT. ONEAUX)
     *            RR = RR - (Q11(J,I,L,K) + Q11(J,I,K,L))
                  QS11(KL,IJ) = RR
                  IF (I .EQ. J) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  IF (K .EQ. L) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  RR =   (QQIJKL - QQIJLK) 
     *                 - (Q11(I,J,K,L) - Q11(I,J,L,K))
                  IF (.NOT. ONEAUX)
     *            RR = RR - (Q11(J,I,L,K) - Q11(J,I,K,L))
                  QT11(KL,IJ) = RR
                  QS11(KL,IJ) = QS11(KL,IJ) * DP25 * DP5
                  QT11(KL,IJ) = QT11(KL,IJ) * DP75 * DP5
C
                  IF (IFLAG .EQ. 2) GOTO 303
                  RR =   (QQIJKL + QQIJLK)
     *                 + (Q11(I,J,K,L) + Q11(I,J,L,K)) * FACX
     *                 - (R11(I,J,K,L) + R11(I,J,L,K))
                  IF (.NOT. ONEAUX) RR = RR
     *                 + (Q11(J,I,L,K) + Q11(J,I,K,L)) * FACX
     *                 - (R11(J,I,L,K) + R11(J,I,K,L))
                  RS11(KL,IJ) = RR
                  IF (I .EQ. J) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  IF (K .EQ. L) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  RR =   (QQIJKL - QQIJLK)
     *                 + (Q11(I,J,K,L) - Q11(I,J,L,K)) * FACX
     *                 - (R11(I,J,K,L) - R11(I,J,L,K))
                  IF (.NOT. ONEAUX) RR = RR
     *                 + (Q11(J,I,L,K) - Q11(J,I,K,L)) * FACX
     *                 - (R11(J,I,L,K) - R11(J,I,K,L))
                  RT11(KL,IJ) = RR
                  RS11(KL,IJ) = RS11(KL,IJ) * DP25 * DP5
                  RT11(KL,IJ) = RT11(KL,IJ) * DP75 * DP5
  303          CONTINUE
  302       CONTINUE
  301    CONTINUE
  300 CONTINUE
C      
      CALL ERISFK(QS11,NSPAIR,1)
      CALL ERISFK(QT11,NSPAIR,1)
      IF (IFLAG .EQ. 1) THEN
         CALL ERISFK(RS11,NSPAIR,1)
         CALL ERISFK(RT11,NSPAIR,1)
      END IF
C
      IF (R12OLD) THEN
         CALL DCOPY(NSPAIR*NSPAIR,QS11,1,RS11,1)
         CALL DCOPY(NSPAIR*NSPAIR,QT11,1,RT11,1)
      ENDIF
C
      IF (IPRINT .LE. 3) GOTO 998
      GOTO (991,992), IFLAG          
      CALL QUIT('INVALID IFLAG IN QCAL1')
  991 CALL AROUND('OLD singlet <Q> matrix') 
      CALL OUTPUT(QS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD triplet <Q> matrix') 
      CALL OUTPUT(QT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      IF (R12OLD) GOTO 998
      CALL AROUND('NEW singlet <Q> matrix') 
      CALL OUTPUT(RS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW triplet <Q> matrix') 
      CALL OUTPUT(RT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      GOTO 998
  992 CALL AROUND('Singlet <Q@> matrix')
      CALL OUTPUT(QS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Triplet <Q@> matrix')
      CALL OUTPUT(QT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI) 
C
  998 LUMUL = -43
      GOTO (996,997), IFLAG          
      CALL QUIT('INVALID IFLAG IN QCAL1')
  996 CALL GPOPEN(LUMUL,'QMATSO','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QS11
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'QMATTO','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QT11
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'QMATSN','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') RS11
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'QMATTN','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') RT11
      CALL GPCLOSE(LUMUL,'KEEP')
      RETURN
  997 IF (IANSAZ .EQ. 2) THEN
      CALL GPOPEN(LUMUL,'QMATSP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QS11
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'QMATTP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QT11
      CALL GPCLOSE(LUMUL,'KEEP')
      ELSE
      CALL GPOPEN(LUMUL,'QMATSQ','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QS11
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'QMATTQ','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QT11
      CALL GPCLOSE(LUMUL,'KEEP')
      END IF
      RETURN
      END
C  /* Deck ccsd_r12 */
      SUBROUTINE CCSD_R12(WORK,LWORK,WXRK,LWXRK,CCR12RSP)
C
C     R12 driver.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <maxorb.h>
#include <priunit.h>
#include <ccorb.h>
#include <iratdef.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <ccsdio.h>
#include <cclr.h>
#include <ccfdgeo.h>
#include <cbirea.h>
#include <r12int.h>
#include <ccr12int.h>
      PARAMETER (D0 = 0D0)
      CHARACTER*6 TIMLAB, LABINT
      CHARACTER*10 APPROX
      DIMENSION WORK(LWORK), WXRK(LWXRK)
      INTEGER IDUM, LU43
      LOGICAL PROFIL, LDUM, TEMP_HERDIR, TEMP_DIRECT, TEMP_ONEAUX, LHTF,
     &        MKVAJKL, CCR12RSP,PROJVIR
      CALL QENTER('CCSD_R12')
      TEMP_DIRECT = DIRECT
      DIRECT = .TRUE.
      MKVAJKL=.FALSE.
      LUSIRG = -30
      LU43   = -43
      PROFIL = .FALSE.
      R12CAL = R12CAL .OR. LMULBS
      IF (R12CAL) THEN
         R12OLD = .TRUE.
         DO ISYM = 1, 8
            R12OLD = R12OLD .AND. NORB2(ISYM) .EQ. 0
         IF (NVIRFR(ISYM) .NE. 0)
     &      CALL QUIT('DO NOT FREEZE VIRTUALS IN R12 CALCULATIONS.')
         END DO
         NORXR = NORXR .OR. R12OLD .OR. R12HYB
         NORXR = NORXR .OR. R12OLD
         ONEAUX = (R12HYB .OR. R12NOB) .AND. .NOT. R12OLD 
     &                                 .AND. .NOT. R12EOR
     &                                 .AND. .NOT. CCR12RSP 
celena
         IF (R12PRP) ONEAUX = .FALSE.
celena
         R12HYB = R12HYB .OR. ONEAUX
         LABEL = 'FULLBAS '
         CALL CCSD_INIT1(WORK,LWORK)
         NORB1T = 0
         NORB2T = 0
         NNBASF = 0
         NSPAIR = NRHFT * (NRHFT + 1) / 2
         NOCXBS = 0
         NVCXBS = 0
         NBASTX = 0
         DO 230 ISYM = 1, NSYM
            NORB1T = NORB1T + NORB1(ISYM)
            NORB2T = NORB2T + NORB2(ISYM)
            NBASF = NORB1(ISYM) + NORB2(ISYM) 
            NNBASF = NNBASF + NBASF * (NBASF + 1) / 2
            NOCXBS = NOCXBS + NRHF(ISYM) * NBAS(ISYM)
            NVCXBS = NVCXBS + NORB1(ISYM) * NBAS(ISYM)
            NBASTX = NBASTX + NBAS(ISYM)
  230    CONTINUE
C 
         WRITE(LUPRI,'(1X,A,8I8)') 
     *     '  NRHF   :',(NRHF(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(1X,A,8I8)')
     *     '  NRHFS  :',(NRHFS(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(1X,A,8I8)') 
     *     '  NORB1  :',(NORB1(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(1X,A,8I8)') 
     *     '  NORB2  :',(NORB2(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(/1X,A,I12)')
     *     'Number of primary orbitals             : ',NORB1T
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of secondary orbitals           : ',NORB2T
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Total number of orbitals               : ',NORB2T + NORB1T
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Total number of basis functions        : ',NBASTX
         IF (ONEAUX) THEN
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (i|p) integrals              : ',NH1AMX
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (i|p'') integrals             : ',NT1AMX
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (i|u'') integrals             : ',NOCXBS
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (Ip''|jL) integrals           : ',NF2FRO(1)
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (ip''|jq) integrals           : ',NH2AMX
         ELSE
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (i|p'') integrals             : ',NT1AMX
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (i|u'') integrals             : ',NOCXBS
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (Ip''|1/r12|jL) integrals     : ',NF2FRO(1)
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (ip''|1/r12|jq'') integrals    : ',NT2AMX
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (ip''|r12|jq'') integrals      : ',NT2AMX
            WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (ip''|[T1,r12]|jq'') integrals : ',NU2AMX
         END IF
         WRITE(LUPRI,'(1X,A,G13.6)')
     *     'Threshold for R12 contributions            =',VCLTHR
         WRITE(LUPRI,'(1X,A,G13.6/)')
     *     'Threshold for singular value decomposition =',SVDTHR
         IF (R12EOR) THEN
          IF (SLATER) THEN
            WRITE(LUPRI,'(1X,A/1X,A)')  
     *        'Calculation with Gaussian-type geminals (GTGs)',
     *        'Expansion coefficients and exponents:'
          ELSE
            WRITE(LUPRI,'(1X,A/1X,A)')  
     *        'Calculation with Gaussian-damped linear-r12 terms',
     *        'Expansion coefficients and exponents:'
          END IF
            DO I = 1, NTOGAM
               WRITE(LUPRI,'(I5,2F25.15)') I, GAMAB(I), GAMAX(I)
            END DO
         ENDIF
         IF (R12DIA) THEN
            WRITE(LUPRI,'(1X,A/A/)')
     *        'The B matrix is transformed into an orthogonal basis',
     *        ' and inverted by diagonalization'
         ELSE IF (R12SVD) THEN
            WRITE(LUPRI,'(1X,A/)')
     *        'The B matrix is inverted by singular value decomposition'
         ELSE
            CALL QUIT('No solver selected')
         END IF
         IF (R12RST) WRITE(LUPRI,'(1X,A)')
     *     'The iterative local-MP2-R12 solver is restarted'
         IF (R12LEV .NE. D0) WRITE(LUPRI,'(1X,A,G13.6)')
     *     'Level-shift in local-MP2-R12 solver        =',R12LEV
         WRITE(LUPRI,'(1X,78A1)') ('*',I=1,78)
         IF (R12HYB) THEN
            WRITE(LUPRI,'(1X,A,4(/1X,A))')
     *       'The MP2-R12/A and MP2-R12/B energies'//
     *       ' are computed as described in:',
     *       'W. Klopper,'//
     *       ' J. Chem. Phys. 120, 10890-10895 (2004).',
     *       'To avoid the hybride scheme, the keyword ''.NO HYB''',
     *       'must be specified in the Dalton input file.',
     *       '----------------------------------------'//
     *       '--------------------------------------'
         ELSE
            WRITE(LUPRI,'(1X,A,2(/1X,A))')
     *       'The MP2-R12/A and MP2-R12/B energies'//
     *       ' are computed as described in:',
     *       'W. Klopper and C. C. M. Samson,'//
     *       ' J. Chem. Phys. 116, 6397-6410 (2002).',
     *       '----------------------------------------'//
     *       '--------------------------------------'
         ENDIF
         IF (R12NOA) THEN
            WRITE(APPROX,'(A10)') ' B.       '
         ELSE IF (R12NOB) THEN
            WRITE(APPROX,'(A10)') ' A.       '
         ELSE
            WRITE(APPROX,'(A10)') 's A and B.'
         END IF
         IF (NOTTWO .AND. NOTTRE) THEN
            WRITE(LUPRI,'(A)')
     *       ' Ansatz 1 is used and results are'//
     *       ' printed for approximation'//APPROX
         ELSE IF (NOTONE .AND. NOTTRE) THEN
            WRITE(LUPRI,'(A)')
     *       ' Ansatz 2 is used and results are'//
     *       ' printed for approximation'//APPROX
         ELSE IF (NOTONE .AND. NOTTWO) THEN
            WRITE(LUPRI,'(A)')
     *       ' Ansatz 3 is used and results are'//
     *       ' printed for approximation'//APPROX
         ELSE IF (NOTTRE) THEN
            WRITE(LUPRI,'(A)')
     *       ' Ansaetze 1 and 2 are used and results are'//
     *       ' printed for approximation'//APPROX
         ELSE IF (NOTTWO) THEN
            WRITE(LUPRI,'(A)')
     *       ' Ansaetze 1 and 3 are used and results are'//
     *       ' printed for approximation'//APPROX
         ELSE IF (NOTONE) THEN
            WRITE(LUPRI,'(A)')
     *       ' Ansaetze 2 and 3 are used and results are'//
     *       ' printed for approximation'//APPROX
         ELSE
            WRITE(LUPRI,'(A)')
     *       ' Ansaetze 1, 2, 3 are used and results are'//
     *       ' printed for approximation'//APPROX
         END IF
         IF (R12EOR) WRITE(LUPRI,'(A,3(/A))')
     *    '----------------------------------------'//
     *    '--------------------------------------',
     *    ' The Gaussian-damped R12 integrals'//
     *    ' are computed as described in:',
     *    ' C. C. M. Samson, W. Klopper and T. Helgaker,',
     *    ' Comp. Phys. Commun. 149, 1-10 (2002).'
         WRITE(LUPRI,'(1X,78A1/)') ('*',I=1,78)
C
         KFOCKD = 1
         KFOCKQ = KFOCKD + NORBTS
         KSF = KFOCKQ + 2 * NRHFT ** 4
         KTF = KSF + MAX(NSPAIR ** 2, NNBASF)
         KT1AM = KTF + MAX(NSPAIR ** 2, NNBASF)
         KGXAM = KT1AM + NT1AMX * 10
         IF (R12PRP) THEN
            KQQAM = KGXAM + NT1VMX
         ELSE
            KQQAM = KGXAM + NT1AMX
         ENDIF
         IF (R12EOR) THEN
            KQQ2M = KQQAM + NRHFT ** 4
            KQQ4M = KQQ2M + NRHFT ** 4
            KQQ6M = KQQ4M + NRHFT ** 4
            KCMOM = KQQ6M + NRHFT ** 4
            KQVAM = KCMOM + 2*NLAMDS
         ELSE
            KQQ2M = KQQAM
            KQQ4M = KQQ2M
            KQQ6M = KQQ4M
            KCMOM = KQQ6M
            KQVAM = KCMOM + NRHFT**4
         ENDIF
celena
         KUGAM = KQVAM + NBAST*NRHFT**3
  
celena 
         IF (R12PRP) THEN
            KFRIN = KUGAM + NVCXBS
         ELSE
            KFRIN = KUGAM + NOCXBS
         ENDIF
C
celena
         KFRIN1= KFRIN  + NFROVR(1)
         KS2AM = KFRIN1 + NFROVF(1)
celena
         IF (ONEAUX) THEN
           KT2AM = KS2AM + NH2AMX
           KU2AM = KT2AM
           KR2AM = KT2AM + NH2AMX
cWK        IF (R12XXL) THEN
              KV2AM = KR2AM + NH2AMX
              KW2AM = KV2AM + NH2AMX
cWK        ELSE
c             KV2AM = KR2AM
c             KW2AM = KV2AM
cWK        END IF
           KEND1 = KW2AM + NH2AMX
           LENT2AM = NH2AMX
         ELSE
           KT2AM = KS2AM + NT2AMX
           KU2AM = KT2AM
           KR2AM = KT2AM + NT2AMX
           KV2AM = MAX(KR2AM + NT2AMX, KU2AM + NU2AMX)
           KW2AM = KV2AM + NT2AMX
           KEND1 = KW2AM + NT2AMX
           LENT2AM = NT2AMX
         END IF
         WRITE(LUPRI,'(A,I15,A/)')' Core memory required by CCSD_R12 :',
     &                             KEND1,'  (double words)'
         LWRK1 = LWORK - KEND1
         IF (KEND1 .GT. LWORK) CALL STOPIT('CCSD_R12','R12 WORK SPACE',
     &        LWORK,KEND1)
         IF (R12EOR) THEN
C
            LUSIFC = -1
            CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
            REWIND LUSIFC
            CALL MOLLAB(LABEL,LUSIFC,LUPRI)
            READ (LUSIFC)
            READ (LUSIFC)
            READ (LUSIFC) (WORK(KCMOM+I-1),I=1,NLAMDS)
            CALL DCOPY(NLAMDS,WORK(KCMOM),1,WORK(KCMOM+NLAMDS),1)
            KCMO1 = KCMOM
            DO ISYM = 1, NSYM
               KCMO2 = KCMO1 + NRHFS(ISYM)*NBAS(ISYM)
               LENMO = NBAS(ISYM)*NRHFS(ISYM)
               CALL DCOPY(LENMO,WORK(KCMO1),1,WORK(KCMO2),1)
               KCMO1 = KCMO1 + 
     &         (NRHFS(ISYM)+NORB1(ISYM)+NORB2(ISYM))*NBAS(ISYM)
            ENDDO
            REWIND LUSIFC
            CALL MOLLAB(LABEL,LUSIFC,LUPRI)
            READ (LUSIFC)
            READ (LUSIFC)
            WRITE (LUSIFC) (WORK(KCMOM+I-1),I=1,NLAMDS)
            CALL GPCLOSE(LUSIFC,'KEEP')
C
            R12EIN = .TRUE.
C
C           Correlation factor: f12 = r12 * exp(-gamma*r12**2)
C
C           COMPUTE (IA|[[F12,T1],F12]|JB) INTEGRALS
C
            CALL TIMER('START ',TIMSTR,TIMEND)
            R12TRA = .TRUE.
            R12INT = .FALSE.
            R12SQR = .FALSE.
            U12INT = .FALSE.
            U21INT = .FALSE.
            V12INT = .FALSE.
            MBSMAX = 4
            INTGAC = 6
            CALL DZERO(WORK(KT1AM),NT1AMX)
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            CALL CCSD_IKJL(WORK(KQQ6M),NRHFT,WXRK(KS2AM))
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ik|[[f12,T1],f12]|jl) integrals done' 
            WRITE(TIMLAB,1013) '#6 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
C
C           COMPUTE (IA|F12**2|JB) INTEGRALS
C
            CALL TIMER('START ',TIMSTR,TIMEND)
            R12TRA = .TRUE.
            R12INT = .FALSE.
            R12SQR = .FALSE.
            U12INT = .FALSE.
            U21INT = .FALSE.
            V12INT = .FALSE.
            MBSMAX = 4
            INTGAC = 4
            CALL DZERO(WORK(KT1AM),NT1AMX)
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            CALL CCSD_IKJL(WORK(KQQ4M),NRHFT,WXRK(KS2AM))
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ik|f12**2|jl) integrals done' 
            WRITE(TIMLAB,1013) '#4 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
C
C           COMPUTE (IA|(1/R12)*F12|JB) INTEGRALS
C
            CALL TIMER('START ',TIMSTR,TIMEND)
            R12TRA = .TRUE.
            R12INT = .FALSE.
            R12SQR = .FALSE.
            U12INT = .FALSE.
            U21INT = .FALSE.
            V12INT = .FALSE.
            CCR12SM = .TRUE.
            MBSMAX = 4
            INTGAC = 2
            CALL DZERO(WORK(KT1AM),NT1AMX)
            LHTF = CCR12
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1) 
            CALL CCSD_IKJL(WORK(KQQ2M),NRHFT,WXRK(KS2AM))
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ik|f12/r12|jl) integrals done' 
            WRITE(TIMLAB,1013) '#2 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
            R12EIN = .FALSE.
            CCR12SM = .FALSE.
            LUSIFC = -1
            CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
            REWIND LUSIFC
            CALL MOLLAB(LABEL,LUSIFC,LUPRI)
            READ (LUSIFC)
            READ (LUSIFC)
            WRITE (LUSIFC) (WORK(NLAMDS+KCMOM+I-1),I=1,NLAMDS)
            CALL GPCLOSE(LUSIFC,'KEEP')
         ELSE
            CALL TIMER('START ',TIMSTR,TIMEND)  
         END IF
C
C        COMPUTE (IA|R12|JB) INTEGRALS
C     
         R12TRA = .TRUE.
         IF (R12EOR) THEN
            R12EIN = .TRUE.
            R12INT = .FALSE.
            INTGAC = 3
         ELSE
            R12INT = .TRUE.   
         END IF  
         R12SQR = .FALSE.
         U12INT = .FALSE.
         U21INT = .FALSE.
         V12INT = .FALSE.
         IF (R12OLD) THEN
            MBSMAX = 4
         ELSE IF ((R12NOB.AND.(.NOT.CCR12RSP)) .OR. 
     &               (NORXR.AND.(.NOT.CCR12RSP))) THEN
            MBSMAX = 5
         ELSE
            MBSMAX = 6
         END IF
         CALL DZERO(WORK(KT1AM),NT1AMX)
C........LHTF must be variable (WK/14-06-2004).
         LHTF = CCR12.AND.(CCR12RSP.OR.IANR12.EQ.2.OR.IANR12.EQ.3)
         CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),LHTF,
     &                  CCR12RSP,MKVAJKL,WORK(KEND1),LWRK1)
c------------------------------------------------------
chf      get r12 integrals on file
c------------------------------------------------------
         CALL GPOPEN(LU43,FR12R12,
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         WRITE(LU43) (WXRK(KS2AM+I-1), I = 1,NH2AMX)
         CALL GPCLOSE(LU43,'KEEP')  
         IF (R12EOR) THEN
            IF (ONEAUX) THEN
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|f12|jq) integrals done' 
            ELSE
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|f12|jq'') integrals done' 
            END IF
            WRITE(TIMLAB,1013) '#3 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            IF (ONEAUX) THEN
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|r12|jq) integrals done' 
            ELSE
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|r12|jq'') integrals done' 
            END IF
            WRITE(TIMLAB,1013) 'R  ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         END IF
C
         IF (CC2R12INT .OR. CCSDR12INT .OR. R12PRP) THEN
           IF (IANR12.EQ.1.OR.IANR12.EQ.3) THEN
             PROJVIR = (IANR12.EQ.3)
             CALL R12BATF(WXRK(KS2AM),FNBACK,DUMMY,.FALSE.,
     &                    PROJVIR,WORK(KEND1),LWRK1)
           END IF
           IF (IANR12.EQ.2.OR.IANR12.EQ.3) THEN
             CALL CC_R12BATF2(WORK(KEND1),LWRK1)
           END IF
           WRITE(LUPRI,'(1X,A)')
     *       'Backtransformation of (ip''|r12|jq'') integrals done'
           CALL TIMER('R12BATF',TIMSTR,TIMEND)  
CCN   
           CALL CC_R12MKSAK(WORK(KEND1),LWRK1)

           IF (CCR12RSP) THEN
             CALL CC_R12VXINT2(WXRK(KS2AM),WORK(KEND1),LWRK1)
           END IF
CCN
         END IF
celena
         IF (IANR12.EQ.1.AND.R12PRP) THEN
            MKVAJKL = .TRUE.
            FNVAJKL = 'CCR12XAJKL'
            CALL DZERO(WORK(KT1AM),NT1AMX)
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            MKVAJKL = .FALSE.
         ENDIF
celena
C
C        COMPUTE (IA|1/R12|JB) INTEGRALS
C
         R12TRA = .TRUE.
         R12EIN = .FALSE.
         R12INT = .FALSE. 
         R12SQR = .FALSE.
         U12INT = .FALSE.
         U21INT = .FALSE.
         V12INT = .TRUE.
         IF (PROFIL .OR. R12ECO) THEN
            MBSMAX = 6
         ELSE IF (R12OLD) THEN
            MBSMAX = 4
         ELSE
            MBSMAX = 5
         END IF
         CALL DZERO(WORK(KT1AM),NT1AMX)
         IF (CC2R12INT .OR. (R12PRP.AND.IANR12.EQ.1)) THEN
            MKVAJKL = .TRUE.            
            IF (R12PRP.AND. IANR12.EQ.1) FNVAJKL = 'CCR12VAJKL'
         END IF
C........LHTF must be variable (WK/14-06-2004).
         LHTF = CCR12.AND.(IANR12.EQ.2.OR.IANR12.EQ.3)
         CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),LHTF,
     &                  CCR12RSP,MKVAJKL,WORK(KEND1),LWRK1)
         MKVAJKL = .FALSE. 
         CALL GPOPEN(LU43,FR12V12,
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         WRITE(LU43) (WXRK(KS2AM+I-1), I = 1,NH2AMX)
         CALL GPCLOSE(LU43,'KEEP')  
C
         IF (R12PRP.AND. IANR12.EQ.1) THEN
            CALL R12BATF(WXRK(KS2AM),FV12BACK,DUMMY,.FALSE.,.FALSE.,
     &                   WORK(KEND1),LWRK1)
         END IF
C
         IF (R12ECO) GOTO 1094 
         IF (FROIMP) THEN
            call dzero(WORK(KFRIN),NFROVR(1)) 
            call dzero(WORK(KFRIN1),NFROVF(1)) 
            IF (IANR12.EQ.1.AND.R12PRP) THEN
                LUNIT = -1
                CALL GPOPEN(LUNIT,'INCJDA',
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
                REWIND(LUNIT)
                READ(LUNIT) (WORK(KFRIN+I-1), I = 1,NFROVR(1))
                CALL GPCLOSE(LUNIT,'DELETE')  

                LUNIT = -1
                CALL GPOPEN(LUNIT,'INCJDI',
     &                         'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
                REWIND(LUNIT)
                READ(LUNIT) (WORK(KFRIN1+I-1), I = 1,NFROVF(1))
                CALL GPCLOSE(LUNIT,'DELETE')  
            ELSE  
                CALL GPOPEN(LU43,'INCJDK',
     &                         'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
                REWIND(LU43)
                READ(LU43) (WORK(KFRIN+I-1), I = 1,NF2FRO(1))
                CALL GPCLOSE(LU43,'KEEP')  

            ENDIF 

         END IF
         CALL MAKEGX(WXRK(KS2AM),WORK(KGXAM),WORK(KFRIN),WORK(KEND1),
     &               LWRK1,2)
celena
         CALL MAKEUG(WORK(KGXAM),WORK(KUGAM),
     &               WORK(KQQAM),WORK(KQVAM),WORK(KEND1),LWRK1,2) 
         IF (R12HYB) THEN
            CALL MAKEGX(WXRK(KS2AM),WORK(KGXAM),WORK(KFRIN),WORK(KEND1),
     &                  LWRK1,1)
            CALL MAKEUG(WORK(KGXAM),WORK(KUGAM),
     &                  WORK(KQQAM),WORK(KQVAM),WORK(KEND1),LWRK1,1) 
celena
         END IF

         IF (IANR12.EQ.1.AND.R12PRP) CALL CC_R12MKSFR(WXRK(KS2AM),
     &                                   WORK(KFRIN1),WORK(KEND1),LWRK1)
 1094    CONTINUE
         CALL GPOPEN(LU43,FR12R12,
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         READ(LU43) (WXRK(KS2AM+I-1), I = 1,NH2AMX)
         CALL GPCLOSE(LU43,'KEEP')  
         IF (ONEAUX) THEN
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|1/r12|jq) integrals done' 
         ELSE
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|1/r12|jq'') integrals done' 
         END IF
         IF (R12EOR) THEN
            WRITE(TIMLAB,1013) '#1 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            WRITE(TIMLAB,1013) '1/R',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         END IF
         IF (PROFIL .OR. R12ECO) THEN
           GOTO 1096
         ELSE IF (R12OLD .OR. R12NOB .OR. NORXR) THEN
           GOTO 1095
         ELSE
            CALL RXR(WXRK(KS2AM),WXRK(KT2AM),
     *               WORK(KSF),WORK(KTF),
     *               WORK(KFOCKQ),NRHFT,NSPAIR)
            WRITE(LUPRI,'(1X,A)')
     *      'Computation of (ip''|r12|jq'')X(p''r'')'//
     *      '(kr''|r12|lq'') integrals done' 
            CALL TIMER('RXR   ',TIMSTR,TIMEND)  
         END IF
C
 1095    CONTINUE
C 
         IF (IANR12.EQ.1.AND.R12PRP) THEN
            MKVAJKL = .TRUE.            
chf         FNBACK = 'CCV12_BACK'
            FNVAJKL = 'CCR12VIJAL'
            R12TRA = .TRUE.
            V12INT = .FALSE.
            IF (R12EOR) THEN
               R12EIN = .TRUE.
               R12INT = .FALSE.
               INTGAC = 3
            ELSE
               R12INT = .TRUE.
            END IF
            CALL DZERO(WORK(KT1AM),NT1AMX)
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
 
            MKVAJKL = .FALSE.           
         END IF
C
C        COMPUTE (I'A|R12|JB) INTEGRALS
C     
         COMBSS = .NOT. R12HYB
         TEMP_ONEAUX = ONEAUX
         R12TRA = .TRUE.
         IF (R12EOR) THEN
            ONEAUX = .FALSE.
            R12EIN = .TRUE.
            R12INT = .FALSE.
         ELSE
            R12INT = .TRUE.   
         END IF  
         R12SQR = .TRUE. 
         U12INT = .FALSE.
         U21INT = .FALSE.
         V12INT = .FALSE.
         IF (R12EOR) THEN   
            R12EIN = .TRUE.
            IF (.NOT. R12NOB) THEN
               LUSIFC = -1
               CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &                     IDUMMY,.FALSE.)
               REWIND LUSIFC
               CALL MOLLAB(LABEL,LUSIFC,LUPRI)
               READ (LUSIFC)
               READ (LUSIFC)
               READ (LUSIFC) (WORK(KCMOM+I-1),I=1,NLAMDS)
               CALL DCOPY(NLAMDS,WORK(KCMOM),1,WORK(KCMOM+NLAMDS),1)
               KCMO1 = KCMOM
               DO ISYM = 1, NSYM
                  KCMO2 = KCMO1 + NRHFS(ISYM)*NBAS(ISYM)
                  LENMO = NBAS(ISYM)*NRHFS(ISYM)
                  CALL DCOPY(LENMO,WORK(KCMO1),1,WORK(KCMO2),1)
                  KCMO1 = KCMO1 + 
     &            (NRHFS(ISYM)+NORB1(ISYM)+NORB2(ISYM))*NBAS(ISYM)
               ENDDO
               REWIND LUSIFC
               CALL MOLLAB(LABEL,LUSIFC,LUPRI)
               READ (LUSIFC)
               READ (LUSIFC)
               WRITE (LUSIFC) (WORK(KCMOM+I-1),I=1,NLAMDS)
               CALL GPCLOSE(LUSIFC,'KEEP')
               INTGAC = 4       
               IF (R12OLD .OR. R12HYB) THEN
                  MBSMAX = 4
               ELSE
                  MBSMAX = 5
               END IF
               CALL DZERO(WORK(KT1AM),NT1AMX)
               LHTF = .FALSE.
               CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),
     &                        LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
               CALL CCSD_IKJL(WORK(KQQAM),NRHFT,WXRK(KU2AM))
               LUSIFC = -1
               CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &                     IDUMMY,.FALSE.)
               REWIND LUSIFC
               CALL MOLLAB(LABEL,LUSIFC,LUPRI)
               READ (LUSIFC)
               READ (LUSIFC)
               WRITE (LUSIFC) (WORK(NLAMDS+KCMOM+I-1),I=1,NLAMDS)
               CALL GPCLOSE(LUSIFC,'KEEP')
               WRITE(TIMLAB,1013) '#4 ',MBSMAX
               CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
            END IF
            INTGAC = 3       
         END IF
         IF (R12OLD) THEN
            IF (R12NOB .AND. .NOT. NATVIR) GOTO 1097
            MBSMAX = 4
         ELSE IF (R12NOB .OR. R12HYB) THEN
            MBSMAX = 5
         ELSE
            MBSMAX = 6
         END IF
         ONEAUX = TEMP_ONEAUX
celena
         IF (R12PRP) ONEAUX = .FALSE.
celena
         IF (IANR12.EQ.1.AND.R12PRP) THEN
            MBSMAX  = 5
            MKVAJKL = .TRUE.
            FNVAJKL= 'CCR12QAJKL' 
            R12INT = .TRUE.
            R12SQR = .FALSE.
            U12INT = .FALSE.
            U21INT = .FALSE.
            V12INT = .FALSE.
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            R12SQR = .TRUE.
            MKVAJKL = .TRUE.
            FNVAJKL = 'CCR12QIJAL' 
            LHTF = .FALSE.
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            MBSMAX  = 5
            CALL R12BATF(WXRK(KU2AM),FQ12BACK,DUMMY,.FALSE.,.FALSE.,
     &                   WORK(KEND1),LWRK1)
            FNVAJKL = 'CCR12UIJAL' 
            R12SQR = .FALSE.
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            CALL R12BATF(WXRK(KU2AM),FU12BACK,DUMMY,.FALSE.,.FALSE.,
     &                   WORK(KEND1),LWRK1)
            FNVAJKL = 'CCR12UAJKL' 
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            R12SQR = .TRUE.
            MKVAJKL = .FALSE.
         ELSE
            CALL DZERO(WORK(KT1AM),NT1AMX)
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
         ENDIF
         CALL QCAL(WORK(KQQAM),WXRK(KS2AM),WXRK(KU2AM),
     *             WORK(KEND1),LWRK1,1,1)
         IF (.NOT. R12OLD) THEN
            CALL QCAL(WORK(KQQAM),WXRK(KS2AM),WXRK(KU2AM),
     *                WORK(KEND1),LWRK1,2,2)
            CALL QCAL(WORK(KQQAM),WXRK(KS2AM),WXRK(KU2AM),
     *                WORK(KEND1),LWRK1,2,3)
            IF (R12HYB) THEN
               IF (R12EOR) THEN
                  IF (ONEAUX) THEN
                     WRITE(LUPRI,'(1X,A)')
     *                 'Computation of (i* p''|f12|jq) integrals done' 
                  ELSE
                     WRITE(LUPRI,'(1X,A)')
     *               'Computation of (i* p''|f12|jq'') integrals done' 
                  END IF
                  WRITE(TIMLAB,1013) '#3 ',MBSMAX
                  CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
               ELSE
                  IF (ONEAUX) THEN
                     WRITE(LUPRI,'(1X,A)')
     *                  'Computation of (i* p''|r12|jq) integrals done' 
                  ELSE
                     WRITE(LUPRI,'(1X,A)')
     *                'Computation of (i* p''|r12|jq'') integrals done' 
                  END IF
                  WRITE(TIMLAB,1013) 'R  ',MBSMAX
                  CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
               END IF
               COMBSS = .TRUE.  
               CALL DZERO(WORK(KT1AM),NT1AMX)
               LHTF = .FALSE.
               CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),
     &                       LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            END IF
            CALL COMKR(WXRK(KS2AM),WXRK(KU2AM),WORK(KV2AM),
     *                 WORK(KW2AM),WORK(KEND1),LWRK1)
         END IF
         IF (R12HYB) THEN
            NLBINT = 3
            LABINT(1:NLBINT) = 'i''p'
         ELSE
            NLBINT = 4
            LABINT(1:NLBINT) = 'i''p'''
         END IF
         IF (R12EOR) THEN
            IF (ONEAUX) THEN
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of ('//LABINT(1:NLBINT)//
     *           '|f12|jq) integrals done' 
            ELSE
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of ('//LABINT(1:NLBINT)//
     *           '|f12|jq'') integrals done' 
            END IF
            WRITE(TIMLAB,1013) '#3 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            IF (ONEAUX) THEN
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of ('//LABINT(1:NLBINT)//
     *           '|r12|jq) integrals done' 
            ELSE
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of ('//LABINT(1:NLBINT)//
     *           '|r12|jq'') integrals done' 
            END IF
            WRITE(TIMLAB,1013) 'R  ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         END IF
C
C        COMPUTE (IA|[T1,R12]|JB) INTEGRALS
C
 1097    CONTINUE
         R12TRA = .TRUE.
         IF (R12EOR) THEN
            R12EIN = .TRUE.
            INTGAC = 5
         END IF  
         R12SQR = .FALSE.
         R12INT = .FALSE.
         U12INT = .TRUE. 
         U21INT = .FALSE.
         V12INT = .FALSE.
         IF (R12OLD) THEN
            MBSMAX = 4
         ELSE
            MBSMAX = 5
         END IF
C
         TEMP_HERDIR = HERDIR
         TEMP_ONEAUX = ONEAUX
         IF (R12EOR) ONEAUX = .FALSE.
         IF (R12PRP .OR. ONEAUX) THEN
            HERDIR = .TRUE.
            U21INT = .TRUE. 
         END IF
         IF (IANR12.EQ.1.AND.R12PRP) THEN
            MKVAJKL = .TRUE.            
chf         FNBACK = 'CCR12_BACK'
            FNVAJKL = 'CCR12BAJKL'
         END IF
         CALL DZERO(WORK(KT1AM),NT1AMX)
         IF (ONEAUX) THEN
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
         ELSE
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            CALL TESTUU(WXRK(KU2AM),WXRK(KS2AM))
         END IF
         HERDIR = TEMP_HERDIR
         ONEAUX = TEMP_ONEAUX
         MKVAJKL = .FALSE.
         U21INT = .FALSE.
C
         IF (IANR12.EQ.1.AND.R12PRP) THEN
chf         FNBACK = 'CCT12_BACK'
            CALL R12BATF(WXRK(KS2AM),FT12BACK,DUMMY,.FALSE.,.FALSE.,
     &                   WORK(KEND1),LWRK1)
            R12TRA = .TRUE.
            IF (R12EOR) THEN
               R12EIN = .TRUE.
               R12INT = .FALSE.
               INTGAC = 3
            ELSE
               R12INT = .TRUE.
            END IF
            R12SQR = .FALSE.
            U12INT = .FALSE.
            U21INT = .FALSE.
            V12INT = .FALSE.
            IF (R12OLD) THEN
               MBSMAX = 4
            ELSE
               MBSMAX = 5
            END IF
            MKVAJKL = .TRUE.            
chf         FNBACK = 'CCT12_BACK'
            FNVAJKL = 'CCR12BIJAL'
            CALL DZERO(WORK(KT1AM),NT1AMX)
            LHTF = .FALSE.
            CALL CCSD_IAJB(WXRK(KT2AM),WORK(KT1AM),
     &                     LHTF,.FALSE.,MKVAJKL,WORK(KEND1),LWRK1)
            MKVAJKL = .FALSE.           
         END IF
C
         IF (R12EOR) THEN
            IF (ONEAUX) THEN
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|[T,f12]|jq) integrals done' 
            ELSE
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|[T,f12]|jq'') integrals done' 
            END IF
            WRITE(TIMLAB,1013) '#5 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            IF (ONEAUX) THEN
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|[T,r12]|jq) integrals done' 
            ELSE
               WRITE(LUPRI,'(1X,A)')
     *           'Computation of (ip''|[T,r12]|jq'') integrals done' 
            END IF
            WRITE(TIMLAB,1013) 'T,R',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         END IF
C
 1096    CONTINUE

         IF (CC2R12INT.AND.(IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
           ! Resort integrals on WORK(KS2AM) to the order needed in
           ! CC-R12, result returned on WXRK(KT2AM)
           CALL CC_R12GETKIAJB(WORK(KS2AM),WXRK(KT2AM),ONEAUX,
     &                NELEM,NH2AMX,IH2AM,NH1AM,IH1AM,NORB1,NRHF,NRHFSA) 
           
           CALL GPOPEN(LU43,FTR12,
     &                      'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
           WRITE(LU43) (WXRK(KT2AM+I-1), I = 1,NELEM)
           CALL GPCLOSE(LU43,'KEEP')  

           ! Resort integrals on WORK(KV2AM) to the order needed in
           ! CC-R12, result returned on WXRK(KT2AM)
           CALL CC_R12GETKIAJB(WORK(KV2AM),WXRK(KT2AM),ONEAUX,
     &                NELEM,NH2AMX,IH2AM,NH1AM,IH1AM,NORB1,NRHF,NRHFSA) 
           CALL GPOPEN(LU43,FKR12,'UNKNOWN','SEQUENTIAL',
     &                 'UNFORMATTED',IDUM,LDUM)
           WRITE(LU43) (WXRK(KT2AM+I-1), I = 1,NELEM)
           CALL GPCLOSE(LU43,'KEEP')  
         END IF
C
         IF (CCSDR12INT.AND.(IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
           !write integrals r_(bp')^(mn) to file
           CALL CC_R12GETRAMNP(WORK(KS2AM),'R12RMNAP',ONEAUX,NORB1,
     &                         NORB2,NRHF,NRHFSA,IH1AM,IH2AM,
     &                         WXRK(KT2AM),LENT2AM)
         END IF
C
         CALL GPOPEN(LU43,FR12V12,
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         READ(LU43) (WXRK(KT2AM+I-1), I = 1,NH2AMX)
         CALL GPCLOSE(LU43,'KEEP')  
         CALL GPOPEN(LU43,FR12R12,
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         READ(LU43) (WXRK(KR2AM+I-1), I = 1,NH2AMX)
         CALL GPCLOSE(LU43,'KEEP')  
         CALL GPOPEN(LUSIRG,'SIRIFC',
     &                      'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         REWIND LUSIRG
         CALL MOLLAB('SIR IPH ',LUSIRG,LUPRI)
         READ (LUSIRG) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,
     *              NACTEL,LSYM
         ESCF = EMCSCF
         CALL MOLLAB(LABEL,LUSIRG,LUPRI)
         READ (LUSIRG)
         READ (LUSIRG) (WORK(KFOCKD+I-1), I=1,NORBTS)
         CALL GPCLOSE(LUSIRG,'KEEP')
         IF (FROIMP .OR. FROEXP)
     *       CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
         CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C        COMPUTE MP2-R12 ENERGY
C
         IF (PROFIL) THEN
            CALL PRFDRV(WXRK(KT2AM),WXRK(KR2AM),WORK(KFOCKD),
     *                  WORK(KEND1),LWRK1)                    
            CALL QUIT('PROFILES ALL DONE')
         ELSE
            CALL R12DRV(WXRK(KT2AM),WXRK(KR2AM),WXRK(KS2AM),
     *                  WORK(KV2AM),WORK(KW2AM),
     *                  WORK(KFOCKD),WORK(KSF),WORK(KTF),
     *                  WORK(KEND1),LWRK1,
     *                  WORK(KQQ2M),WORK(KQQ4M),WORK(KQQ6M))
         END IF
         CALL TIMER('E(R12)',TIMSTR,TIMEND)  

         IF (CC2R12INT.AND.(IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
           ! Resort R12 integrals into WORK, sort them to the order
           ! needed in CC-R12, result returned on WXRK, and write them
           ! back to file
           CALL GPOPEN(LU43,FR12R12,
     &                      'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
           READ(LU43) (WORK(I), I = 1,NH2AMX)
           CALL CC_R12GETKIAJB(WORK,WXRK(KR2AM),ONEAUX,NELEM,
     &                NH2AMX,IH2AM,NH1AM,IH1AM,NORB1,NRHF,NRHFSA) 
           REWIND(LU43)
           WRITE(LU43) (WXRK(KR2AM + I - 1), I = 1,NELEM)
           CALL GPCLOSE(LU43,'KEEP')  
         END IF
      ENDIF
      R12TRA = .FALSE.
      R12INT = .FALSE.
      R12SQR = .FALSE.
      U12INT = .FALSE.
      U21INT = .FALSE.
      R12EIN = .FALSE.
      V12INT = .TRUE.
      INTGAC = 0
      MBSMAX = 4
      ONEAUX = .FALSE.
      DIRECT = TEMP_DIRECT
      CALL QEXIT('CCSD_R12')
      RETURN
 1013 FORMAT(A3,'[',I1,']')  
      END
C  /* Deck ccsd_ikjl */
      SUBROUTINE CCSD_IKJL(QQ,NOCDIM,XX)
C
C     This subroutine collects the (ik|jl) MO-integrals.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XX(*),QQ(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
#include <r12int.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J 
      CALL DZERO(QQ,NOCDIM**4)
      DO ISYMIK = 1, NSYM
         ISYMJL = ISYMIK
         DO ISYMI = 1, NSYM
            ISYMK = MULD2H(ISYMI,ISYMIK)
            DO ISYMJ = 1, NSYM
               ISYML = MULD2H(ISYMJ,ISYMJL)
               DO I = 1, NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  DO K = NRHFFR(ISYMK) + 1, NRHFS(ISYMK)
                     KOFFK = IRHF(ISYMK) + K - NRHFFR(ISYMK)
                     DO J = 1, NRHF(ISYMJ)
                        KOFFJ = IRHF(ISYMJ) + J
                        DO L = NRHFFR(ISYML) + 1, NRHFS(ISYML)
                           KOFFL = IRHF(ISYML) + L - NRHFFR(ISYML)
                           NKI = IT1AM(ISYMK,ISYMI) +
     *                           NVIR(ISYMK)*(I-1) + K
                           NLJ = IT1AM(ISYML,ISYMJ) +
     *                           NVIR(ISYML)*(J-1) + L
                           IF (ONEAUX) THEN
                              IF (R12SQR) THEN
                                 NKILJ = IU2AM(ISYMIK,ISYMJL)
     *                                 + NT1AM(ISYMJL)*(NKI - 1)
     *                                 + NLJ 
                              ELSE
                                 NKILJ = IH2AM(ISYMIK,ISYMJL)
     *                                 + INDEX(NKI,NLJ)  
                              END IF
                           ELSE
                              IF (R12SQR) THEN
                                 NKILJ = IU2AM(ISYMIK,ISYMJL)
     *                                 + NT1AM(ISYMJL)*(NKI - 1)
     *                                 + NLJ 
                              ELSE
                                 NKILJ = IT2AM(ISYMIK,ISYMJL)
     *                                 + INDEX(NKI,NLJ)  
                              END IF
                           END IF
                           QQ(KOFFI,KOFFJ,KOFFK,KOFFL)=XX(NKILJ)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF (IPRINT .GE. 5) THEN
         NP = NOCDIM**2
         CALL AROUND('Output from CCSD_IKJL')
         CALL OUTPUT(QQ,1,NP,1,NP,NP,NP,1,LUPRI)
      END IF
      RETURN
      END 
C  /* Deck vclean */
      SUBROUTINE VCLEAN(A,N,THR)  
C
C     This subroutine removes small matrix elements.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0)
      DIMENSION A(N,N)
      DO 100 I = 1, N
      DO 100 J = 1, N
         IF (ABS(A(I,J)) .LT. THR) A(I,J) = 0D0
  100 CONTINUE
      RETURN
      END
C  /* Deck r12mp2 */
      SUBROUTINE R12MP2(VS,VT,VOS,VOT,SS,ST,EKL,NOCDIM,NSPAIR,FS,FT,
     *            BB,UU,WW,QQ,SF,TF,EV,XF,VV,IJPAIR,WS,WT,DELTA,CCS,CCT)
C
C     This subroutine evaluates the MP2-R12 energy by
C     solving sets of linear equations using
C     singular value decomposition.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h> 
#include <r12int.h>
      PARAMETER (D0 = 0D0, DP5 = 5D-1, D2 = 2D0, D99 = 1D99)
      DIMENSION VS(NSPAIR), VT(NSPAIR), EV(NSPAIR)
      DIMENSION SS(NSPAIR,NSPAIR), ST(NSPAIR,NSPAIR)
      DIMENSION SF(NSPAIR,NSPAIR), TF(NSPAIR,NSPAIR)
      DIMENSION BB(NSPAIR,NSPAIR), QQ(NSPAIR,NSPAIR) 
      DIMENSION VV(NSPAIR,NSPAIR), WW(NSPAIR)
      DIMENSION UU(NSPAIR,NSPAIR)
      DIMENSION CCS(NSPAIR),CCT(NSPAIR)
      DIMENSION VOS(NSPAIR,NSPAIR),VOT(NSPAIR,NSPAIR)
      LOGICAL MATV, MATU
      MATV = .TRUE.
      MATU = .TRUE.
      IERR = 0
      WS = D99
      WT = D99   
C
C     SINGLET PAIRS
C
      IF (R12ECO) THEN
        DO J = 1, NSPAIR
         DO I = 1, NSPAIR
            BB(I,J) = SF(I,J) - EKL * SS(I,J)
         END DO
        END DO
      ELSE
        DO J = 1, NSPAIR
         DO I = 1, NSPAIR
            BB(I,J) = - SF(I,J) - SF(J,I) 
     *               - XF * (D2 * EKL - EV(I) - EV(J))
     *                    * (SS(I,J) + SS(J,I))
            BB(I,J) = BB(I,J) * DP5
         END DO
         BB(J,J) = BB(J,J) + DELTA
        END DO
      END IF
      CALL DZERO(WW,NSPAIR)
      CALL DZERO(UU,NSPAIR*NSPAIR)
      CALL DZERO(VV,NSPAIR*NSPAIR)
      CALL DZERO(QQ,NSPAIR*NSPAIR)
      CALL SVD(NSPAIR,NSPAIR,NSPAIR,BB,WW,MATU,UU,MATV,VV,IERR,QQ)
      DO I = 1, NSPAIR 
         IF (ABS(WW(I)) .LT. SVDTHR) WW(I) = D0
         IF (WW(I) .NE. D0) THEN
            WS = MIN(WS, WW(I))
            IF (WW(I) .LT. VCLTHR) WW(I) = D0
         END IF 
      END DO 
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,'(//A,2I10/)') ' SVD ERROR / IJPAIR, IERR =',
     *                                            IJPAIR, IERR
         CALL QUIT('SVD ERROR -- SINGLET')
      END IF
      CALL SVBKSB(UU,WW,VV,NSPAIR,NSPAIR,NSPAIR,NSPAIR,VS,BB,QQ)
      DO I = 1, NSPAIR
          CCS(I)= - BB(I,1)
      END DO
      FS = -DDOT(NSPAIR,VOS,1,BB,1)
C
C     TRIPLET PAIRS
C
      IF (R12ECO) THEN
        DO J = 1, NSPAIR
         DO I = 1, NSPAIR
            BB(I,J) = TF(I,J) - EKL * ST(I,J)
         END DO
        END DO
      ELSE
        DO J = 1, NSPAIR
         DO I = 1, NSPAIR
            BB(I,J) = - TF(I,J) - TF(J,I) 
     *               - XF * (D2 * EKL - EV(I) - EV(J))
     *                    * (ST(I,J) + ST(J,I))
            BB(I,J) = BB(I,J) * DP5
         END DO
         IF (ABS(BB(J,J)) .GT. 1D-12) BB(J,J) = BB(J,J) + DELTA
        END DO
      END IF
      CALL DZERO(WW,NSPAIR)
      CALL DZERO(UU,NSPAIR*NSPAIR)
      CALL DZERO(VV,NSPAIR*NSPAIR)
      CALL DZERO(QQ,NSPAIR*NSPAIR)
      CALL SVD(NSPAIR,NSPAIR,NSPAIR,BB,WW,MATU,UU,MATV,VV,IERR,QQ)
      DO I = 1, NSPAIR 
         IF (ABS(WW(I)) .LT. SVDTHR) WW(I) = D0
         IF (WW(I) .NE. D0) THEN
            WT = MIN(WT, WW(I))
            IF (WW(I) .LT. VCLTHR) WW(I) = D0
         END IF 
      END DO 
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,'(//A,2I10/)') ' SVD ERROR / IJPAIR, IERR =',
     *                                            IJPAIR, IERR
         CALL QUIT('SVD ERROR -- TRIPLET')
      END IF
      CALL SVBKSB(UU,WW,VV,NSPAIR,NSPAIR,NSPAIR,NSPAIR,VT,BB,QQ)
      DO I =1,NSPAIR
         CCT(I) = - BB(I,1)
      END DO
      FT = -DDOT(NSPAIR,VOT,1,BB,1)
      RETURN
      END
C  /* Deck r12fxl */
      SUBROUTINE R12FXL(XS,XT,FLM,SS,ST,NIND,II,JJ,NOCDIM,NSPAIR)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( DP5 = 0.5D0, D1 = 1.0D0, D1M = -1.0D0, D2 = 2.0D0 )
      PARAMETER ( SYMTHR = 1.0D-10 )
      DIMENSION FLM(*)
      DIMENSION II(NSPAIR), JJ(NSPAIR), NIND(NOCDIM,NOCDIM)
      DIMENSION WS(NSPAIR,NSPAIR), WT(NSPAIR,NSPAIR)
      DIMENSION XS(NSPAIR,NSPAIR), XT(NSPAIR,NSPAIR)
      DIMENSION SS(NSPAIR,NSPAIR), ST(NSPAIR,NSPAIR)
C
      DSQ2 = SQRT(D2)
      DSQ2I = D1 / DSQ2
      NSP2 = NSPAIR * NSPAIR
      CALL DZERO(XS,NSP2)
      CALL DZERO(XT,NSP2)
      DO NPIJ = 1, NSPAIR
         NI = II(NPIJ)
         NJ = JJ(NPIJ)
         DO NPKL = 1, NSPAIR
            NK = II(NPKL)
            NL = JJ(NPKL)
            DO NO = 1, NOCDIM
               NIO = NIND(NI,NO)
               NJO = NIND(NJ,NO)
               NKO = NIND(NK,NO)
               NLO = NIND(NL,NO)
               IF (NO.NE.NI) THEN
                  IF (NO.LT.NJ) THEN 
                     FAC = D1M
                  ELSE
                     FAC = D1
                  END IF
                  IF (NO.EQ.NJ) THEN
                     FA = DSQ2
                  ELSE
                     FA = D1
                  END IF
                  IF (NI.EQ.NJ) FA = FA * DSQ2I
                  XS(NPIJ,NPKL) = XS(NPIJ,NPKL) + FA *
     &                            SS(NJO,NPKL)*FLM(NIO) 
                  XT(NPIJ,NPKL) = XT(NPIJ,NPKL) + 
     &                                      FAC * ST(NJO,NPKL)*FLM(NIO)
               END IF
               IF (NO.NE.NJ) THEN
                  IF (NO.GT.NI) THEN
                     FAC = D1M
                  ELSE
                     FAC = D1
                  END IF 
                  IF (NO.EQ.NI) THEN
                     FA = DSQ2
                  ELSE
                     FA = D1
                  END IF
                  IF (NI.EQ.NJ) FA = FA * DSQ2I
                  XS(NPIJ,NPKL) = XS(NPIJ,NPKL) + FA *
     &                            SS(NIO,NPKL)*FLM(NJO)
                  XT(NPIJ,NPKL) = XT(NPIJ,NPKL) +
     &                                      FAC * ST(NIO,NPKL)*FLM(NJO)
               END IF
               IF (NO.NE.NK) THEN
                  IF (NO.LT.NL) THEN
                     FAC = D1M
                  ELSE
                     FAC = D1
                  END IF
                  IF (NO.EQ.NL) THEN
                     FA = DSQ2
                  ELSE
                     FA = D1
                  END IF
                  IF (NK.EQ.NL) FA = FA * DSQ2I
                  XS(NPIJ,NPKL) = XS(NPIJ,NPKL) + FA *
     &                            SS(NPIJ,NLO)*FLM(NKO)
                  XT(NPIJ,NPKL) = XT(NPIJ,NPKL) +
     &                                      FAC * ST(NPIJ,NLO)*FLM(NKO)
               END IF
               IF (NO.NE.NL) THEN
                  IF (NO.GT.NK) THEN
                     FAC = D1M
                  ELSE
                     FAC = D1
                  END IF
                  IF (NO.EQ.NK) THEN
                     FA = DSQ2
                  ELSE
                     FA = D1
                  END IF
                  IF (NK.EQ.NL) FA = FA * DSQ2I
                  XS(NPIJ,NPKL) = XS(NPIJ,NPKL) + FA *
     &                            SS(NPIJ,NKO)*FLM(NLO)
                  XT(NPIJ,NPKL) = XT(NPIJ,NPKL) + 
     &                                      FAC * ST(NPIJ,NKO)*FLM(NLO)
               END IF
            END DO
            XS(NPIJ,NPKL) = DP5 * XS(NPIJ,NPKL)
            XT(NPIJ,NPKL) = DP5 * XT(NPIJ,NPKL)
         END DO
      END DO
      DO NPIJ = 1, NSPAIR
         DO NPKL = 1, NSPAIR
            IF (ABS(XS(NPIJ,NPKL)-XS(NPKL,NPIJ)) .GT. SYMTHR) THEN
               WRITE (LUPRI,'(A,2I5,2E20.10)')
     &        'Singlet matrix <XF+FX> not symmetric',NPIJ,NPKL,
     &         XS(NPIJ,NPKL),XS(NPKL,NPIJ)
               CALL QUIT('Singlet matrix <XF+FX> not symmetric')
            END IF
            IF (ABS(XT(NPIJ,NPKL)-XT(NPKL,NPIJ)) .GT. SYMTHR) THEN
               WRITE (LUPRI,'(A,2I5,2E20.10)')
     &        'Triplet matrix <XF+FX> not symmetric',NPIJ,NPKL,
     &         XT(NPIJ,NPKL),XT(NPKL,NPIJ)
               CALL QUIT('Triplet matrix <XF+FX> not symmetric')
            END IF
         END DO
      END DO
      RETURN
      END
C  /* Deck r12fcl */
      SUBROUTINE R12FCL(ZS,ZT,WS,WT,VS,VT,YS,YT,SS,ST,CCS,CCT,FLM,NIND,
     &                                       II,JJ,NOCDIM,NSPAIR,DELTA)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D1M = -1.0D0 , D2 = 2D0 )
      DIMENSION FLM(*)
      DIMENSION II(NSPAIR), JJ(NSPAIR), NIND(NOCDIM,NOCDIM)
      DIMENSION ZS(NSPAIR,NSPAIR), ZT(NSPAIR,NSPAIR)
      DIMENSION WS(NSPAIR,NSPAIR), WT(NSPAIR,NSPAIR)
      DIMENSION VS(NSPAIR,NSPAIR), VT(NSPAIR,NSPAIR)
      DIMENSION YS(NSPAIR,NSPAIR), YT(NSPAIR,NSPAIR)
      DIMENSION SS(NSPAIR,NSPAIR), ST(NSPAIR,NSPAIR)
      DIMENSION CCS(NSPAIR,NSPAIR), CCT(NSPAIR,NSPAIR)
C
      DSQ2 = SQRT(D2)
      DSQ2I = D1 / DSQ2
      NSP2 = NSPAIR * NSPAIR
      CALL DZERO(YS,NSP2)
      CALL DZERO(YT,NSP2)
      DO NPIJ = 1, NSPAIR
         NI = II(NPIJ)
         NJ = JJ(NPIJ)
         DO NPKL = 1, NSPAIR
            DO NO = 1, NOCDIM
               NIO = NIND(NI,NO)
               NJO = NIND(NJ,NO)
               IF (NO.NE.NI) THEN
                  IF (NO.LT.NJ) THEN
                     FAC = D1M
                  ELSE
                     FAC = D1
                  END IF
                  IF (NO.EQ.NJ) THEN
                     FA = DSQ2
                  ELSE
                     FA = D1
                  END IF
                  YS(NPKL,NPIJ) = YS(NPKL,NPIJ) + FA *
     *                            CCS(NPKL,NJO)*FLM(NIO)
                  YT(NPKL,NPIJ) = YT(NPKL,NPIJ) +
     *                                      FAC * CCT(NPKL,NJO)*FLM(NIO)
               END IF
               IF (NO.NE.NJ) THEN
                  IF (NO.GT.NI) THEN
                     FAC = D1M
                  ELSE
                     FAC = D1
                  END IF
                  IF (NO.EQ.NI) THEN
                     FA = DSQ2
                  ELSE
                     FA = D1
                  END IF
                  YS(NPKL,NPIJ) = YS(NPKL,NPIJ) + FA *
     *                            CCS(NPKL,NIO)*FLM(NJO)
                  YT(NPKL,NPIJ) = YT(NPKL,NPIJ) +
     *                                     FAC *  CCT(NPKL,NIO)*FLM(NJO)
               END IF
            END DO
            IF (NI.EQ.NJ) YS(NPKL,NPIJ) = YS(NPKL,NPIJ) * DSQ2I
         END DO
      END DO
      DO NPIJ = 1, NSPAIR
         DO NPKL = 1, NSPAIR
             AAA = D0
             BBB = D0
             CCC = D0
             DDD = D0
             DO NPMN = 1,NSPAIR
                AAA = AAA + SS(NPIJ,NPMN) * YS(NPMN,NPKL)
                BBB = BBB + WS(NPIJ,NPMN) * CCS(NPMN,NPKL)
                CCC = CCC + ST(NPIJ,NPMN) * YT(NPMN,NPKL)
                DDD = DDD + WT(NPIJ,NPMN) * CCT(NPMN,NPKL)
             END DO
            ZS(NPIJ,NPKL)=VS(NPIJ,NPKL)-DELTA*CCS(NPIJ,NPKL)-AAA+BBB
            ZT(NPIJ,NPKL)=VT(NPIJ,NPKL)-DELTA*CCT(NPIJ,NPKL)-CCC+DDD
         END DO
      END DO
      RETURN
      END 
C  /* Deck svbksb */
      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X,TMP)
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
#include <implicit.h>
      PARAMETER (D0 = 0D0)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(N)
      DO 12 J=1,N
         S=D0
         IF (W(J).NE.D0) THEN
            DO 11 I=1,M
                S=S+U(I,J)*B(I)
   11       CONTINUE
            S=S/W(J)
         ENDIF
         TMP(J)=S
   12 CONTINUE
      DO 14 J=1,N
         S=0.0D0
         DO 13 JJ=1,N
             S=S+V(J,JJ)*TMP(JJ)
   13    CONTINUE
         X(J)=S
   14 CONTINUE
      RETURN
      END    
C  /* Deck ecodrv */
      SUBROUTINE ECODRV(V2AM,R2AM,EV,V11,B11,Q11,VS11,VT11,WS11,WT11,
     *                  QS11,QT11,BS11,BT11,RS11,RT11,
     *                  NICDIM,NOCDIM,NSPAIR,NTPAIR,ES,ET,FS,FT,EVS,
     *                  CNINV,SF,CS11,CT11)
C
C     This subroutine computes the V(klmn), X(klmn), and B(klmn)
C     matrices and evaluates the MP2-R12 energies.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0, DELTA=0.0d0)
      PARAMETER (THRDIA = 1D-9)
      DIMENSION V2AM(*), R2AM(*), EV(*), SF(*)
      REAL*8  R1, R2, V1, VR, RR, BB,
     *        V11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        B11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        Q11(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
      DIMENSION VS11(NSPAIR,NSPAIR), VT11(NSPAIR,NSPAIR)
      DIMENSION WS11(NSPAIR,NSPAIR), WT11(NSPAIR,NSPAIR)
      DIMENSION QS11(NSPAIR,NSPAIR), QT11(NSPAIR,NSPAIR)
      DIMENSION BS11(NSPAIR,NSPAIR), BT11(NSPAIR,NSPAIR)
      DIMENSION RS11(NSPAIR,NSPAIR), RT11(NSPAIR,NSPAIR) 
      DIMENSION ES(NSPAIR),ET(NSPAIR),FS(NSPAIR),FT(NSPAIR)
      DIMENSION EVS(NSPAIR),CNINV(NSPAIR,8)
      DIMENSION CS11(NSPAIR,NSPAIR),CT11(NSPAIR,NSPAIR)
      DIMENSION ISB(8)
      INTEGER IDUM
      LOGICAL LDUM
#include <r12int.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      FF = D1 / SQRT(D2)
      ISB(1) = 0
      DO 100 ISYM = 2, NSYM
         NBASF = NORB1(ISYM-1) + NORB2(ISYM-1)
         NNBASF = NBASF * (NBASF + 1) / 2
         ISB(ISYM) = ISB(ISYM-1) + NNBASF
  100 CONTINUE
      NBASF = NORB1(NSYM) + NORB2(NSYM)
      NNBASF = ISB(NSYM) + NBASF * (NBASF + 1) / 2
      LUMULB = -34
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ (LUMULB,'(4E30.20)') (SF(I), I = 1, NNBASF)
      CALL GPCLOSE(LUMULB,'KEEP')
      DO 200 ISYM = 1, NSYM
         DO 210 I = NORB1(ISYM) + 2,  NVIR(ISYM)
            DO 220 J = NORB1(ISYM) + 1 , I - 1
               IJ = ISB(ISYM) + INDEX(I,J)
               IF (ABS(SF(IJ)) .GT. THRDIA) THEN
                  WRITE(LUPRI,'(/A/A,3I5,E20.10/)')
     *            '@ WARNING : Fock matrix not diagonal',
     *            '            Nondiagonal element is :',
     *                         ISYM,I,J,SF(IJ)
                  IF (ABS(SF(IJ)) .GT. SQRT(THRDIA))
     *               CALL QUIT('Fock matrix not diagonal')
               END IF
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE 
C
      DO 205 I = 1, NOCDIM**4
         Q11(I,1,1,1) = D0
         V11(I,1,1,1) = D0
         B11(I,1,1,1) = D0
  205 CONTINUE
      CALL DZERO(ES,NSPAIR)
      CALL DZERO(ET,NSPAIR)
C
C     CONSTRUCT PRODUCTS (IA|JB)*(KA|LB)
C
      E2 = D0
      E2S = D0
      E2T = D0
      DO 101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         DO 102 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 103 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
C
C              COMPUTE MP2 ENERGY
C
               DO 201 I = 1, NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  DO 202 J = 1, NRHF(ISYMJ)
                     KOFFJ = IRHF(ISYMJ) + J           
                     DO 203 A = NRHFS(ISYMA) + 1, NORB1(ISYMA)
                        ISYMJA = MULD2H(ISYMJ,ISYMA)
                        KOFFA = IVIR(ISYMA) + A 
                        DO 204 B = NRHFS(ISYMB) + 1, NORB1(ISYMB)
                           ISYMIB = MULD2H(ISYMI,ISYMB)
                           KOFFB = IVIR(ISYMB) + B
                           NAI = IT1AM(ISYMA,ISYMI) +
     *                           NVIR(ISYMA)*(I-1) + A
                           NBJ = IT1AM(ISYMB,ISYMJ) +
     *                           NVIR(ISYMB)*(J-1) + B
                           NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                             INDEX(NAI,NBJ)
                           NAJ = IT1AM(ISYMA,ISYMJ) +
     *                           NVIR(ISYMA)*(J-1) + A
                           NBI = IT1AM(ISYMB,ISYMI) +
     *                           NVIR(ISYMB)*(I-1) + B
                           NAJBI = IT2AM(ISYMJA,ISYMIB) +
     *                             INDEX(NAJ,NBI)
                           VAIBJ = V2AM(NAIBJ)
                           VAJBI = V2AM(NAJBI)
                           VV = VAIBJ * (D2 * VAIBJ - VAJBI)
                           DENOM = D1 / (EV(KOFFI) + EV(KOFFJ) -
     *                                   EV(KOFFA) - EV(KOFFB))
                           E2 = E2 + VV * DENOM                     
                           IJ = INDEX(KOFFI,KOFFJ)
                           VS = (VAIBJ + VAJBI)**2
                           VT = (VAIBJ - VAJBI)**2
                           ES(IJ) = ES(IJ) + VS * DENOM
                           ET(IJ) = ET(IJ) + VT * DENOM
  204                   CONTINUE
  203                CONTINUE
  202             CONTINUE
  201          CONTINUE
C
               DO 104 ISYMKA = 1, NSYM
                  ISYMLB = ISYMKA 
                  ISYMK = MULD2H(ISYMA,ISYMKA)
                  ISYML = MULD2H(ISYMB,ISYMLB)
                  DO 105 I = 1, NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     DO 106 K = 1, NRHF(ISYMK)
                        KOFFK = IRHF(ISYMK) + K
                        DO 107 A = NORB1(ISYMA)+1, NVIR(ISYMA)
                           NAI = IT1AM(ISYMA,ISYMI) +
     *                           NVIR(ISYMA)*(I-1) + A
                           NAK = IT1AM(ISYMA,ISYMK) +
     *                           NVIR(ISYMA)*(K-1) + A
                           DO 108 J = 1, NRHF(ISYMJ)
                              KOFFJ = IRHF(ISYMJ) + J
                              DO 109 L = 1, NRHF(ISYML)
                                 KOFFL = IRHF(ISYML) + L
                                 DO 110 B = NORB1(ISYMB)+1, NVIR(ISYMB)
                                    NBJ = IT1AM(ISYMB,ISYMJ) +
     *                                    NVIR(ISYMB)*(J-1) + B
                                    NBL = IT1AM(ISYMB,ISYML) +
     *                                    NVIR(ISYMB)*(L-1) + B
                                    NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                                      INDEX(NAI,NBJ)
                                    NAKBL = IT2AM(ISYMKA,ISYMLB) +
     *                                      INDEX(NAK,NBL)
                                    V1 = V2AM(NAIBJ)
                                    R1 = R2AM(NAIBJ)
                                    R2 = R2AM(NAKBL)
                                    VR = V1 * R2
                                    RR = R1 * R2 
                                    BB = SF(ISB(ISYMA) + INDEX(A,A))
     *                                 + SF(ISB(ISYMB) + INDEX(B,B))
                                    BB = R1 * R2 * BB
C      
                                    V11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              V11(KOFFI,KOFFJ,KOFFK,KOFFL) +  VR
                                    Q11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +  RR
                                    B11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                              B11(KOFFI,KOFFJ,KOFFK,KOFFL) +  BB
C
  111                               CONTINUE
  110                            CONTINUE
  109                         CONTINUE
  108                      CONTINUE
  107                   CONTINUE
  106                CONTINUE
  105             CONTINUE
  104          CONTINUE
  103       CONTINUE
  102    CONTINUE
  101 CONTINUE
C
      WRITE(LUPRI,'(/A/)') ' SECOND-ORDER PAIR ENERGIES:'
      DO 230 I=1,NSPAIR
        ES(I) = ES(I) * DP25
        ET(I) = ET(I) * DP75
        E2S = E2S + ES(I)
        E2T = E2T + ET(I)
        WRITE(LUPRI,'(17X,I4,2F15.9)') I,ES(I),ET(I)
  230 CONTINUE
      WRITE(LUPRI,'(/A6,3F15.9)') ' MP2 =',E2,E2S,E2T
C
      IJ = 0
      DO 300 I = 1, NOCDIM
         DO 301 J = 1, I
            IJ = IJ + 1
            KL = 0
            DO 302 K = 1, NOCDIM
               DO 303 L = 1, K
                  KL = KL + 1
C
                  RR = V11(I,J,K,L) + V11(I,J,L,K)
                  VS11(KL,IJ) = RR
                  IF (I .EQ. J) VS11(KL,IJ) = FF * VS11(KL,IJ)
                  IF (K .EQ. L) VS11(KL,IJ) = FF * VS11(KL,IJ)
                  RR = V11(I,J,K,L) - V11(I,J,L,K)
                  VT11(KL,IJ) = RR
C
                  BB = B11(I,J,K,L) + B11(I,J,L,K)
                  BS11(KL,IJ) = BB
                  IF (I .EQ. J) BS11(KL,IJ) = FF * BS11(KL,IJ)
                  IF (K .EQ. L) BS11(KL,IJ) = FF * BS11(KL,IJ)
                  BB = B11(I,J,K,L) - B11(I,J,L,K)
                  BT11(KL,IJ) = BB
C
                  RR = Q11(I,J,K,L) + Q11(I,J,L,K)
                  RS11(KL,IJ) = RR
                  IF (I .EQ. J) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  IF (K .EQ. L) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  RR = Q11(I,J,K,L) - Q11(I,J,L,K)
                  RT11(KL,IJ) = RR
C
                  VT11(KL,IJ) = VT11(KL,IJ) * D3
                  BT11(KL,IJ) = BT11(KL,IJ) * D3
                  RT11(KL,IJ) = RT11(KL,IJ) * D3
C
  303          CONTINUE
  302       CONTINUE
  301    CONTINUE
  300 CONTINUE
C
      IF (IPRINT .LE. 3) GOTO 998
      CALL AROUND('Singlet <Vab> matrix') 
      CALL OUTPUT(VS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Singlet <inglet <V> matrixab> matrix') 
      CALL OUTPUT(BS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Singlet <Sab> matrix') 
      CALL OUTPUT(RS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Triplet <Vab> matrix') 
      CALL OUTPUT(VT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Triplet <Bab> matrix') 
      CALL OUTPUT(BT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Triplet <Sab> matrix') 
      CALL OUTPUT(RT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
  998 CONTINUE    
C
      CALL DZERO(FS,NSPAIR)
      CALL DZERO(FT,NSPAIR)
      F2S = D0
      F2T = D0
      IJ = 0
      WRITE(LUPRI,'(/A/)') ' SECOND-ORDER MP2-EC PAIR ENERGIES:'
      DO 500 I = 1, NOCDIM
         DO 501 J = 1, I
            IJ = IJ + 1
            IF (R12SVD) THEN
            CALL R12MP2(VS11(1,IJ),VT11(1,IJ),VS11(1,IJ),VT11(1,IJ),
     *      RS11,RT11,EV(I)+EV(J),NOCDIM,NSPAIR,FS(IJ),FT(IJ),WS11,WT11,
     *      QS11,QT11,BS11,BT11,EVS,D1,V11,IJ,WSMIN,WTMIN,DELTA,
     *                                            CS11(1,IJ),CT11(1,IJ))
            ELSE IF (R12DIA) THEN
            CALL MP2R12(VS11(1,IJ),VT11(1,IJ),VS11(1,IJ),VT11(1,IJ),
     *      RS11,RT11,EV(I)+EV(J),NOCDIM,NSPAIR,FS(IJ),FT(IJ),WS11,WT11,
     *      QS11,QT11,BS11,BT11,EVS,D1,V11,IJ,WSMIN,WTMIN,DELTA,
     *                                            CS11(1,IJ),CT11(1,IJ))
            ENDIF
            WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                       ET(IJ),ET(IJ)+FT(IJ),
     *                                       WSMIN,WTMIN 
            F2S = F2S + FS(IJ)
            F2T = F2T + FT(IJ) 
  501    CONTINUE
  500 CONTINUE 
      WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-EC      correlation energy =',
     *   E2S+E2T+F2S+F2T
      RETURN
      END
C  /* Deck prfdrv */ 
      SUBROUTINE PRFDRV(V2AM,R2AM,EV,WRK,LWRK)
C
C     Driver for computation of basis-set profiles.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>  
      DIMENSION R2AM(*), V2AM(*), EV(*), WRK(*)
      NOCDIM = NRHFT 
      NICDIM = NRHFTS
      NSPAIR = NOCDIM * (NOCDIM + 1) / 2
      NTPAIR = NOCDIM * (NOCDIM - 1) / 2
      NPAIR2 = NSPAIR ** 2
      NODIM4 = NOCDIM ** 4
      NIDIM4 = NICDIM ** 4
      IVS = 1             
      IVT = IVS + NOCDIM * NOCDIM   
      IENDW = IVT + NOCDIM * NOCDIM
      IF (IENDW .GT. LWRK) THEN
         WRITE(LUPRI,'(A,I20/A,I20)') 
     *    ' WORK SPACE REQUIRED =  ',IENDW,
     *    ' WORK SPACE AVAILABLE = ',LWRK
         CALL QUIT('INSUFFICIENT WORK SPACE IN R12DRV')
      END IF
      CALL PRFEVR(V2AM,R2AM,WRK(IVS),WRK(IVT),NICDIM,NOCDIM,NSPAIR,
     *            NTPAIR)
      RETURN
      END
C  /* Deck prfevr */
      SUBROUTINE PRFEVR(V2AM,R2AM,VS,VT,NICDIM,NOCDIM,NSPAIR,NTPAIR)
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      DIMENSION R2AM(*), V2AM(*)
      DIMENSION VS(NOCDIM,NOCDIM), VT(NOCDIM,NOCDIM)
      INTEGER IDUM
      LOGICAL LDUM
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL DZERO(VS,NOCDIM**2)
      CALL DZERO(VT,NOCDIM**2)
C
C     CONSTRUCT PRODUCTS (IA|JB)*(IA||JB)
C
      DO 101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         DO 102 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 103 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
               DO 105 I = 1, NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  DO 107 A = 1 + NORB1(ISYMA), NVIR(ISYMA)
                     KOFFA = IVIR(ISYMA) + A
                     NAI = IT1AM(ISYMA,ISYMI) +
     *                     NVIR(ISYMA)*(I-1) + A
                     DO 108 J = 1, NRHF(ISYMJ)
                        KOFFJ = IRHF(ISYMJ) + J
                        NAJ = IT1AM(ISYMA,ISYMJ) +
     *                        NVIR(ISYMA)*(J-1) + A
                        ISYMJA = MULD2H(ISYMJ,ISYMA)   
                        DO 110 B = 1 + NORB1(ISYMB), NVIR(ISYMB)
                           ISYMIB = MULD2H(ISYMI,ISYMB)   
                           KOFFB = IVIR(ISYMB) + B
                           NBJ = IT1AM(ISYMB,ISYMJ) +
     *                           NVIR(ISYMB)*(J-1) + B
                           NBI = IT1AM(ISYMB,ISYMI) +
     *                           NVIR(ISYMB)*(I-1) + B
                           NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                             INDEX(NAI,NBJ)
                           NAJBI = IT2AM(ISYMJA,ISYMIB) +
     *                             INDEX(NAJ,NBI)
                           VR = V2AM(NAIBJ) * R2AM(NAIBJ)
                           VP = V2AM(NAIBJ) * R2AM(NAJBI)
                           VS(KOFFI,KOFFJ) = VS(KOFFI,KOFFJ) + VR + VP
                           VT(KOFFI,KOFFJ) = VT(KOFFI,KOFFJ) + VR - VP
  110                   CONTINUE
  108                CONTINUE
  107             CONTINUE
  105          CONTINUE
  103       CONTINUE
  102    CONTINUE
  101 CONTINUE
C
      LU33 = -33
      LU34 = -34
      LU35 = -35
      CALL GPOPEN(LU33,'PROFILE.1','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      CALL GPOPEN(LU34,'PROFILE.2','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      CALL GPOPEN(LU35,'PROFILE.3','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      DO 201 I = 1, NOCDIM
         DO 202 J = 1, NOCDIM
            WRITE(LU33,'(2I4,2F20.14)') I,J,VS(I,J)       
            WRITE(LU34,'(2I4,2F20.14)') I,J,VT(I,J)  
            WRITE(LU35,'(2I4,2F20.14)') I,J,(VS(I,J)+VT(I,J))*DP5
  202    CONTINUE
         WRITE(LU33,*) 
         WRITE(LU34,*) 
         WRITE(LU35,*) 
  201 CONTINUE
      CALL GPCLOSE(33,'KEEP')
      CALL GPCLOSE(34,'KEEP')
      CALL GPCLOSE(35,'KEEP')
      RETURN
      END
C  /* Deck mpab */
      SUBROUTINE MPAB(A,NROWA,NCOLA,NRDIMA,NCDIMA,
     1                B,NROWB,NCOLB,NRDIMB,NCDIMB,
     2                C,NRDIMC,NCDIMC)
C-----------------------------------------------------------
C MATRIX PRODUCT A TIMES B = C.
C  Written by George D. Purvis 1983
C  Revised 6-Nov-1984 Hans Joergen Aa. Jensen
C-------------------------------------------------------------
#include <implicit.h>
      DIMENSION A(NRDIMA,NCDIMA),B(NRDIMB,NCDIMB),C(NRDIMC,NCDIMC)
      PARAMETER (ZERO=0.D00)
C
      IF (NCOLA.NE.NROWB) THEN
         WRITE (*,9000) NCOLA,NROWB
         CALL QUIT('ERROR, inconsistent matrix dimensions in MPAB')
      ENDIF
 9000 FORMAT(/' MPAB error: NCOLA .ne. NROWB, values =',2I10)
C
      IF (NROWB .EQ. 0) RETURN
      DO 40 J = 1,NCOLB
        DO 10 I = 1,NROWA
   10     C(I,J) = A(I,1)*B(1,J)
        DO 30 K = 2,NROWB
          IF (B(K,J).EQ.ZERO) GO TO 30
          BKJ = B(K,J)
          DO 20 I = 1,NROWA
   20       C(I,J) = BKJ*A(I,K) + C(I,J)
   30   CONTINUE
   40 CONTINUE
C
      RETURN
      END
C=====================================================================*
C  /* Deck vshrnk */
      SUBROUTINE VSHRNK(VS11,NSPAIR,NRHF,NRHFA,NSYM)
#include <implicit.h>
      DIMENSION VS11(NSPAIR,NSPAIR),NRHF(NSYM),NRHFA(NSYM)
      IJ = 0
      JI = 0
      II = 0
      DO ISYM = 1, NSYM
       DO I = 1, NRHF(ISYM)
        II = II + 1
        JJ = 0
        DO JSYM = 1, NSYM
        DO J = 1, NRHF(JSYM)
          JJ = JJ + 1
          IF (JJ .LE. II) THEN
            IJ = IJ + 1
            IF (I .LE. NRHFA(ISYM) .AND.
     *          J .LE. NRHFA(JSYM)) THEN
              JI = JI + 1
              IF (IJ .NE. JI)
     *        CALL DCOPY(NSPAIR,VS11(1,IJ),1,VS11(1,JI),1)
              END IF
            END IF
         END DO
        END DO
       END DO
      END DO
      IF (JI .LT. NSPAIR)
     *CALL DZERO(VS11(1,JI+1),NSPAIR*(NSPAIR-JI))
      RETURN
      END
