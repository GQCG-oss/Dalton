C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
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
      LUMULB = 34
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ (LUMULB,'(4E30.20)') (SF(I), I = 1, NNBASF)
      CALL GPCLOSE (LUMULB,'KEEP')
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
C  
C     Compute (ak|bl) = Sum_c X(ac)*(ck|bl) + Sum_d X(bd)*(ak|dl)
C     The sum over c runs over the auxilliary basis only
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
                        IF (A. LE. NORB1(ISYMA)) THEN
                           NSTRC = NORB1(ISYMC) + 1
                           NENDC = NVIR(ISYMC)
                        ELSE
                           NSTRC = A                
                           NENDC = A                
                        END IF
                        NAK = IT1AM(ISYMA,ISYMK)
     *                      + NVIR(ISYMA)*(K-1) + A
                        DO 307 B = 1, NVIR(ISYMB)
                           IF (B. LE. NORB1(ISYMB)) THEN
                              NSTRD = NORB1(ISYMD) + 1
                              NENDD = NVIR(ISYMD)
                           ELSE
                              NSTRD = B
                              NENDD = B
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
                           IF (IFLAG .EQ. 2) THEN
                              AVIRT = A .GT. NRHFS(ISYMA) .AND.
     *                                A .LE. NORB1(ISYMA)
                              BVIRT = B .GT. NRHFS(ISYMB) .AND.
     *                                B .LE. NORB1(ISYMB)
                              IF (AVIRT .OR. BVIRT) GOTO 407 
                           END IF
                           IF ((A .GT. NORB1(ISYMA) .AND.
     *                          B .GT. NORB1(ISYMB)) .OR.
     *                         (A .LE. NORB1(ISYMA) .AND.
     *                          B .LE. NORB1(ISYMB))) THEN
                              SFAC = D1
                           ELSE
                              SFAC = - D1
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
         GOTO (703,704), IFLAG
         CALL QUIT('INVALID IFLAG IN RXR')
  703    CALL AROUND('Singlet <RXR> matrix')
         CALL OUTPUT(SF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <RXR> matrix')
         CALL OUTPUT(TF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         GOTO 705
  704    CALL AROUND('Singlet <RXR@> matrix')
         CALL OUTPUT(SF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <RXR@> matrix')
         CALL OUTPUT(TF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
  705    WRITE(LUPRI,*)
      END IF
C
      LUMUL = 34
      GOTO (701,702), IFLAG
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
  999 CONTINUE
      RETURN
      END 
C  /* Deck makegx */   
      SUBROUTINE MAKEGX(V2AM,GX,F2AM,IPBAS)
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
#include <implicit.h>
#include <priunit.h>    
      DIMENSION V2AM(*), GX(*), F2AM(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfro.h>
#include <r12int.h> 
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL DZERO(GX,NT1AMX)
C
C     Compute SUM_i ( i p' | i k )
C
      DO 101 ISYMIP = 1, NSYM
         ISYMIK = ISYMIP
         DO 102 ISYMI = 1, NSYM
            ISYMP = MULD2H(ISYMI,ISYMIP)
            ISYMK = ISYMP  
            DO 201 I = 1, NRHF(ISYMI)
               NPK = IT1AM(ISYMP,ISYMK)
               DO 202 K = NRHFFR(ISYMK)+1, NRHFS(ISYMK)
                  DO 203 P = 1, NVIR(ISYMP)
                     NPK = NPK + 1
                     NPI = IT1AM(ISYMP,ISYMI)
     *                   + NVIR(ISYMP)*(I-1) + P
                     NKI = IT1AM(ISYMK,ISYMI)
     *                   + NVIR(ISYMK)*(I-1) + K
                     NPIKI = IT2AM(ISYMIP,ISYMIK)
     *                     + INDEX(NPI,NKI)  
                     GX(NPK) = GX(NPK) + V2AM(NPIKI)
  203             CONTINUE
  202          CONTINUE
  201       CONTINUE 
  102    CONTINUE
  101 CONTINUE  
C
      IF (FROIMP) THEN
C
C     Add contribution from frozen orbitals
C
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
                        NKI = IT1AM(ISYMI,ISYMK)
     *                      + NVIR(ISYMI)*(K-1) + I 
                        NPI = IT1FRO(ISYMP,ISYMI)
     *                      + NVIR(ISYMP)*(I-1) + P
                        NPIKI = IT2FRO(ISYMIP,ISYMIK)
     *                        + NT1FRO(ISYMIP)*(NKI-1)   
     *                        + NPI
                         GX(NPK) = GX(NPK) + F2AM(NPIKI)
  403             CONTINUE
  402          CONTINUE
  401       CONTINUE 
  302    CONTINUE
  301 CONTINUE  
      END IF   
C
      IF (R12OLD) GOTO 700
C
C     Zero primary block 
C
      IF (.NOT. R12HYB .OR. (R12HYB .AND. IPBAS .EQ. 2)) THEN
         DO 500 ISYM = 1, NSYM
            NPK = IT1AM(ISYM,ISYM) + 1
            DO 501 K = 1, NRHF(ISYM)
               CALL DZERO(GX(NPK),NORB1(ISYM))
               NPK = NPK + NVIR(ISYM)
  501       CONTINUE
  500    CONTINUE
      END IF
C
C     Zero secondary block 
C
      IF (R12HYB .AND. IPBAS .EQ. 1) THEN
         DO 600 ISYM = 1, NSYM
            NPK = IT1AM(ISYM,ISYM) + NORB1(ISYM) + 1
            DO 601 K = 1, NRHF(ISYM)
               CALL DZERO(GX(NPK),NORB2(ISYM))
               NPK = NPK + NVIR(ISYM)
  601       CONTINUE
  600    CONTINUE
      END IF 
C
C     Print section
C
  700 CONTINUE
C
      IF (IPRINT .GT. 3) THEN
         DO 300 ISYM = 1, NSYM
            IF (NRHF(ISYM)*NVIR(ISYM) .NE. 0) 
     *      WRITE(LUPRI,'(/A,I2)') ' X-MATRIX/ISYM =',ISYM
            NPK = IT1AM(ISYM,ISYM) + 1
            CALL OUTPUT(GX(NPK),1,NVIR(ISYM),1,NRHF(ISYM),
     *                  NVIR(ISYM),NRHF(ISYM),1,LUPRI) 
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
#include <r12int.h>  
      DIMENSION R2AM(*), V2AM(*), U2AM(*), EV(*), WRK(*),
     *          RXRS(*), RXRT(*), S2AM(*), T2AM(*),
     *          QQ2(*), QQ4(*), QQ6(*)
C
      NOCDIM = NRHFT 
      NICDIM = NRHFTS
      NSPAIR = NOCDIM * (NOCDIM + 1) / 2
      NTPAIR = NOCDIM * (NOCDIM - 1) / 2
      NPAIR2 = NSPAIR ** 2
      NODIM4 = 2 * NOCDIM ** 4
      NIDIM4 = NICDIM ** 4
C
      IES = 1
      IET = IES + NSPAIR
      IFS = IET + NSPAIR
      IFT = IFS + NSPAIR
      IEVS = IFT + NSPAIR 
      ININV = IEVS + NSPAIR
      IQQ11 = ININV + NSPAIR*8
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
      IENDW = IRT11 + NPAIR2
C
      IF (IENDW .GT. LWRK) THEN
         WRITE(LUPRI,'(A,I20/A,I20)') 
     *    ' WORK SPACE REQUIRED =  ',IENDW,
     *    ' WORK SPACE AVAILABLE = ',LWRK
         CALL QUIT('INSUFFICIENT WORK SPACE IN R12DRV')
      ENDIF
      IF (R12ECO) THEN
         CALL ECODRV(V2AM,R2AM,EV,WRK(IV11),WRK(IB11),WRK(IQ11),
     *               WRK(IVS11),WRK(IVT11),WRK(IWS11),WRK(IWT11),
     *               WRK(IQS11),WRK(IQT11),
     *               WRK(IBS11),WRK(IBT11),WRK(IRS11),
     *               WRK(IRT11),NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *               WRK(IES),WRK(IET),WRK(IFS),WRK(IFT),
     *               WRK(IEVS),WRK(ININV),WRK(IENDW))
          RETURN
      ENDIF 
      CALL AROUND('MP2 ansatz := (1 - P1 - P2 + P1*P2) * R12')
      CALL MAKEVR(V2AM,R2AM,U2AM,EV,WRK(IQQ11),
     *     WRK(IV11),WRK(IU11),WRK(IB11),
     *     WRK(IW11),WRK(IQ11),WRK(IR11),
     *     WRK(IVS11),WRK(IUS11),WRK(IBS11),
     *     WRK(IWS11),WRK(IQS11),WRK(IRS11),
     *     WRK(IVT11),WRK(IUT11),WRK(IBT11),
     *     WRK(IWT11),WRK(IQT11),WRK(IRT11),
     *     NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *     WRK(IES),WRK(IET),WRK(IFS),WRK(IFT),
     *     WRK(IEVS),WRK(ININV),RXRS,RXRT,QQ2,QQ4,QQ6) 
      IF (R12OLD .OR. .NOT. R12XXL) RETURN
      CALL AROUND('MP2 ansatz := (1 - O1 - O2 + O1*O2) * R12')
      CALL MAKERV(V2AM,R2AM,U2AM,S2AM,T2AM,EV,WRK(IQQ11),
     *     WRK(IV11),WRK(IU11),WRK(IB11),
     *     WRK(IW11),WRK(IQ11),WRK(IR11),
     *     WRK(IVS11),WRK(IUS11),WRK(IBS11),
     *     WRK(IWS11),WRK(IQS11),WRK(IRS11),
     *     WRK(IVT11),WRK(IUT11),WRK(IBT11),
     *     WRK(IWT11),WRK(IQT11),WRK(IRT11),
     *     NICDIM,NOCDIM,NSPAIR,NTPAIR,
     *     WRK(IES),WRK(IET),WRK(IFS),WRK(IFT),
     *     WRK(IEVS),WRK(ININV),QQ2,QQ4,QQ6)
      RETURN
      END
C  /* Deck makevr */
      SUBROUTINE MAKEVR(V2AM,R2AM,U2AM,EV,QQ11,V11,U11,B11,W11,Q11,R11,
     *                  VS11,US11,BS11,WS11,
     *                  QS11,RS11,VT11,UT11,BT11,WT11,QT11,RT11,
     *                  NICDIM,NOCDIM,NSPAIR,NTPAIR,ES,ET,FS,FT,EVS,
     *                  CNINV,RXRS,RXRT,QQ2,QQ4,QQ6)
C
C     This subroutine computes the V(klmn), X(klmn), and B(klmn)
C     matrices and evaluates the MP2-R12 energies.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      DIMENSION R2AM(*), V2AM(*), U2AM(*), EV(*)
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
      DIMENSION EVS(NSPAIR),CNINV(NSPAIR,8)
      DIMENSION QQ11(NICDIM,NICDIM,NICDIM,NICDIM)
      INTEGER IDUM
      LOGICAL LDUM
#include <r12int.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      FF = D1 / SQRT(D2)
      LU33 = 33
      CALL GPOPEN(LU33,'AUXQ12','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ(LU33,'(4E30.20)') QQ11
      CALL GPCLOSE(LU33,'KEEP')
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
                     RR = RR + (U11(I,J,K,L) + U11(I,J,L,K))
                     RR = RR - (B11(I,J,K,L) + B11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =   (U11(I,J,K,L) + U11(I,J,L,K))
     *                    - (B11(I,J,K,L) + B11(I,J,L,K))
                     DX = D2
                  END IF
                  BS11(KL,IJ) = RR
                  IF (I .EQ. J) BS11(KL,IJ) = FF * BS11(KL,IJ)
                  IF (K .EQ. L) BS11(KL,IJ) = FF * BS11(KL,IJ) 
                  IF (IJ .EQ. KL) BS11(KL,IJ) = BS11(KL,IJ) - DX 
                  IF (R12EOR) THEN
                     RR = - D2 *(QQ6(I,J,K,L) - QQ6(I,J,L,K))
                     RR = RR + (U11(I,J,K,L) - U11(I,J,L,K))
                     RR = RR - (B11(I,J,K,L) - B11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =   (U11(I,J,K,L) - U11(I,J,L,K))
     *                    - (B11(I,J,K,L) - B11(I,J,L,K))
                     DX = D2
                  END IF
                  BT11(KL,IJ) = RR
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            BT11(KL,IJ) = BT11(KL,IJ) - DX
C
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) + QQ2(I,J,L,K))
                     RR = RR + (V11(I,J,K,L) + V11(I,J,L,K))
                     RR = RR - (W11(I,J,K,L) + W11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =   (V11(I,J,K,L) + V11(I,J,L,K))
     *                    - (W11(I,J,K,L) + W11(I,J,L,K))
                     DX = D1
                  END IF
                  WS11(KL,IJ) = RR
                  IF (I .EQ. J) WS11(KL,IJ) = FF * WS11(KL,IJ)
                  IF (K .EQ. L) WS11(KL,IJ) = FF * WS11(KL,IJ) 
                  IF (IJ .EQ. KL) WS11(KL,IJ) = WS11(KL,IJ) + DX 
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) - QQ2(I,J,L,K))
                     RR = RR + (V11(I,J,K,L) - V11(I,J,L,K))
                     RR = RR - (W11(I,J,K,L) - W11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR =   (V11(I,J,K,L) - V11(I,J,L,K))
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
                  RR = RR + (Q11(I,J,K,L) + Q11(I,J,L,K))
     *                    - (R11(I,J,K,L) + R11(I,J,L,K))
                  RS11(KL,IJ) = RR
                  IF (I .EQ. J) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  IF (K .EQ. L) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  RR = QQIJKL - QQIJLK
                  RR = RR + (Q11(I,J,K,L) - Q11(I,J,L,K))
     *                    - (R11(I,J,K,L) - R11(I,J,L,K))
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
      IF (IPRINT .LE. 3) GOTO 998
      CALL AROUND('OLD singlet <V> matrix') 
      CALL OUTPUT(VS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD singlet <B> matrix') 
      CALL OUTPUT(US11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD singlet <S> matrix') 
      CALL OUTPUT(QS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD triplet <V> matrix') 
      CALL OUTPUT(VT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD triplet <B> matrix') 
      CALL OUTPUT(UT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('OLD triplet <S> matrix') 
      CALL OUTPUT(QT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      IF (R12OLD) GOTO 998
      CALL AROUND('NEW singlet <V> matrix') 
      CALL OUTPUT(WS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW singlet <B> matrix') 
      CALL OUTPUT(BS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW singlet <S> matrix') 
      CALL OUTPUT(RS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW triplet <V> matrix') 
      CALL OUTPUT(WT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW triplet <B> matrix') 
      CALL OUTPUT(BT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('NEW triplet <S> matrix') 
      CALL OUTPUT(RT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
  998 LUMULO = 34
      LUMULN = LUMULO + 1
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
      IJ = 0
      DO 490 I = 1, NOCDIM
         DO 491 J = 1, I
            IJ = IJ + 1
            EVS(IJ) = EV(I) + EV(J)
  491    CONTINUE
  490 CONTINUE
C
      LU43 = 43
      IF (R12OLD) THEN
        LUMULE = LUMULO  
      ELSE
        LUMULE = LUMULN  
      END IF
      DO 999 LUMULX = LUMULO, LUMULE
      IF (LUMULX .EQ. LUMULO) THEN
        CALL AROUND('Original MP2-R12 method')
        IF (.NOT. R12NOB) THEN
          CALL GPOPEN(LU43,'QMATSO','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') US11
          CALL GPCLOSE(LU43,'KEEP')
          CALL GPOPEN(LU43,'QMATTO','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') UT11
          CALL GPCLOSE(LU43,'KEEP')
        END IF
      ELSE
        CALL AROUND('Improved MP2-R12 method')
        IF (.NOT. R12NOB) THEN
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
      DO 999 IPDD = 1, 7
      IF (IPDD .EQ. 4) THEN
        DO I = 1, NSPAIR
           DO J = 1, NSPAIR
              BS11(I,J) = BS11(I,J) - US11(I,J)
              BT11(I,J) = BT11(I,J) - UT11(I,J)
           END DO
        END DO
      END IF
      IF (IPDD .EQ. 6) THEN
        DO I = 1, NSPAIR
           DO J = 1, NSPAIR
             BS11(I,J) = BS11(I,J) + RXRS(I,J)
             BT11(I,J) = BT11(I,J) + RXRT(I,J)
           END DO
        END DO
      END IF
      IF (LUMULX .EQ. LUMULO .OR. NORXR) THEN
        IF (.NOT. R12XXL .AND. 
     *      (IPDD .NE. 2 .AND. IPDD .NE. 3 .AND. IPDD .NE. 5)) GOTO 999
      ELSE
        IF (.NOT. R12XXL .AND. 
     *      (IPDD .NE. 2 .AND. IPDD .NE. 3 .AND. IPDD .NE. 7)) GOTO 999
      END IF
      IF (IPDD .LE. 2 .AND. R12NOA) GOTO 999
      IF (IPDD .EQ. 3 .AND. R12NOP) GOTO 999
      IF (IPDD .GE. 4 .AND. R12NOB) GOTO 999
      IF (IPDD .GE. 6 .AND. (LUMULX .EQ. LUMULO .OR. NORXR)) GOTO 999 
      IF (IPDD .EQ. 1) THEN
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
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B'//
     *                   ' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 5) THEN
        XPDD = DP5    
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 6) THEN  
        XPDD = D0    
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B'''//
     *                   ' PAIR ENERGIES:' 
      ELSE IF (IPDD .EQ. 7) THEN
        XPDD = DP5    
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B'' PAIR ENERGIES:'
      END IF
      CALL DZERO(FS,NSPAIR)
      CALL DZERO(FT,NSPAIR)
      F2S = D0
      F2T = D0
      IJ = 0
      DO 500 I = 1, NOCDIM
         DO 501 J = 1, I
            IJ = IJ + 1
            CNINV(IJ,1) = WS11(IJ,IJ)
            CNINV(IJ,2) = BS11(IJ,IJ)
            IF (ABS(CNINV(IJ,2)) .GT. VCLTHR) THEN
               CNINV(IJ,3) = CNINV(IJ,1) / CNINV(IJ,2)
            ELSE
               CNINV(IJ,3) = D0
            END IF
            CNINV(IJ,4) = CNINV(IJ,1) * CNINV(IJ,3)
            IF (I .NE. J) THEN
               CNINV(IJ,5) = WT11(IJ,IJ)
               CNINV(IJ,6) = BT11(IJ,IJ)
               IF (ABS(CNINV(IJ,6)) .GT. VCLTHR) THEN
                  CNINV(IJ,7) = CNINV(IJ,5) / CNINV(IJ,6)
               ELSE
                  CNINV(IJ,7) = D0
               END IF
               CNINV(IJ,8) = CNINV(IJ,5) * CNINV(IJ,7) 
            ELSE
               CNINV(IJ,5) = D0
               CNINV(IJ,6) = D0
               CNINV(IJ,7) = D0
               CNINV(IJ,8) = D0
            END IF 
            FS(IJ) = CNINV(IJ,4) 
            FT(IJ) = CNINV(IJ,8) 
            IF (IPDD .NE. 1 .AND. IPDD .NE. 4 .AND. IPDD .NE. 6) THEN
            CALL R12MP2(WS11(1,IJ),WT11(1,IJ),RS11,RT11,EVS(IJ),
     *      NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,QS11,QT11,
     *      BS11,BT11,EVS,XPDD,QQ11,IJ,WSMIN,WTMIN)
            WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                       ET(IJ),ET(IJ)+FT(IJ),
     *                                       WSMIN,WTMIN           
            ELSE
            WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                       ET(IJ),ET(IJ)+FT(IJ)
            END IF
            F2S = F2S + FS(IJ)
            F2T = F2T + FT(IJ) 
  501    CONTINUE
  500 CONTINUE 
      WRITE(LUPRI,'(/4X,6F15.9)') E2S,E2S+F2S,E2T,E2T+F2T
      IF (IPDD .EQ. 1) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)') 
     * ' Noninvariant MP2-R12/A   correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 600 I = 1, NOCDIM
            DO 601 J = 1, I       
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  601       CONTINUE
  600    CONTINUE
      ELSE IF (IPDD .EQ. 2) THEN
         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/A   correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 3) THEN 
         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/A''  correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 4) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)') 
     * ' Noninvariant MP2-R12/B   correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 602 I = 1, NOCDIM
            DO 603 J = 1, I       
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  603       CONTINUE
  602    CONTINUE
      ELSE IF (IPDD .EQ. 5) THEN
         WRITE(LUPRI,'(/A,F15.9 )')
     * '              MP2-R12/B   correlation energy =',
     *   E2S+E2T+F2S+F2T    
      ELSE IF (IPDD .EQ. 6) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)')
     * ' Noninvariant MP2-R12/B''  correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 604 I = 1, NOCDIM
            DO 605 J = 1, I
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  605       CONTINUE
  604    CONTINUE
      ELSE IF (IPDD .EQ. 7) THEN
         WRITE(LUPRI,'(/A,F15.9 )')
     * '              MP2-R12/B''  correlation energy =',
     *   E2S+E2T+F2S+F2T  
      END IF
  999 CONTINUE
      RETURN
      END
C  /* Deck mp2r12 */
      SUBROUTINE MP2R12(VS,VT,SS,ST,EKL,NOCDIM,NSPAIR,FS,FT,
     *                  BB,UU,PP,QQ,SF,TF,EV,XF,RR,IJPAIR,
     *                  WS,WT)
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
     &       ' @ WARNING: SMALL EIGENVALUE OF SINGLET'//
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
      FS = D0
      DO 111 I = 1, NSPAIR
         DO 112 J = 1, NSPAIR
            FS = FS - PP(I) * QQ(I,J) * PP(J)
  112    CONTINUE
  111 CONTINUE
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
     &       ' @ WARNING: SMALL EIGENVALUE OF TRIPLET'//
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
      FT = D0
      DO 311 I = 1, NSPAIR
         DO 312 J = 1, NSPAIR
            FT = FT - PP(I) * QQ(I,J) * PP(J)
  312    CONTINUE
  311 CONTINUE 
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
      SUBROUTINE MAKEUG(GX,UG,QQ,WORK,LWORK,IPBAS)
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
      DIMENSION GX(*), UG(*), WORK(*), QQ(*)
      CHARACTER*7 GUNAM
      INTEGER IDUM
      LOGICAL LDUM
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J    
C
      NOCCT = NRHFT
      LUSIRG = 30
      IF ((NLAMDS + 2*NOCCT*NBAST + 14*NOCCT*NOCCT) .GT. LWORK)
     *CALL QUIT('Not enough space in MAKEUG')
      CALL GPOPEN(LUSIRG,'SIRIFC','OLD',' ','UNFORMATTED',IDUM,LDUM)
      REWIND LUSIRG
C
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
      DO 100 ISYM = 1, NSYM
         IWORK = IWORK + NBAS(ISYM)*NRHFS(ISYM)
         CALL MPAB(WORK(IWORK),NBAS(ISYM),NVIR(ISYM),
     *                         NBAS(ISYM),NVIR(ISYM),
     *             GX(IGX),NVIR(ISYM),NRHF(ISYM),
     *                     NVIR(ISYM),NRHF(ISYM),
     *             UG(IUG),NBAS(ISYM),NRHF(ISYM))
         IF (IPRINT .GT. 3) THEN
            IF (NRHF(ISYM)*NBAS(ISYM) .NE. 0) 
     *      WRITE(LUPRI,'(/A,I2)') ' U x X/ISYM =',ISYM
            CALL OUTPUT(UG(IUG),1,NBAS(ISYM),1,NRHF(ISYM),
     *                  NBAS(ISYM),NRHF(ISYM),1,LUPRI)  
         END IF
         IWORK = IWORK + NBAS(ISYM)*NVIRS(ISYM)
         IGX = IGX + NVIR(ISYM)*NRHF(ISYM)
         IUG = IUG + NBAS(ISYM)*NRHF(ISYM)
  100 CONTINUE
      NTOTGU = IUG - 1
      WRITE(GUNAM,'(A6,I1)') 'GUMAT.',IPBAS
      LU43 = 43
      CALL GPOPEN(LU43,GUNAM,'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
      REWIND(LU43)
      WRITE(LU43) NTOTGU
      WRITE(LU43) (UG(I), I = 1, NTOTGU)
      CALL GPCLOSE(LU43,'KEEP') 
C
      ISMOU = 1              
      ISMO0 = 1              
      ISMO1 = ISMO0 + NLAMDS
      ISMO2 = ISMO1 + NOCCT*NBAST
      CALL DZERO(WORK(ISMO1),NOCCT*NBAST)
      CALL DZERO(WORK(ISMO2),NOCCT*NBAST)
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
         ISMO1 = ISMO1 + NBAS(ISYM)
         ISMO2 = ISMO2 + NBAS(ISYM)
  505 CONTINUE 
      ISMO0 = 1              
      ISMO1 = ISMO0 + NLAMDS
      ISMO2 = ISMO1 + NOCCT*NBAST
      ISMS0 = ISMO2 + NOCCT*NBAST
      ISMX1 = ISMS0 + NOCCT*NOCCT*2
      ISMY1 = ISMX1 + NOCCT*NOCCT*2
      ISMZ1 = ISMY1 + NOCCT*NOCCT*2
      ISMX2 = ISMZ1 + NOCCT*NOCCT*2
      ISMY2 = ISMX2 + NOCCT*NOCCT*2
      ISMZ2 = ISMY2 + NOCCT*NOCCT*2
      IEND  = ISMZ2 + NOCCT*NOCCT*2
      LQQ   = LWORK - IEND
      IF (IPRINT .GT. 3) THEN
         CALL AROUND('Active orbitals')
         CALL OUTPUT(WORK(ISMO1),1,NBAST,1,NOCCT,NBAST,NOCCT,1,LUPRI)
         CALL AROUND('Primed active orbitals')
         CALL OUTPUT(WORK(ISMO2),1,NBAST,1,NOCCT,NBAST,NOCCT,1,LUPRI)
      END IF
      IF (LQQ .LT. 0) CALL QUIT('No space for R12QQ')
      CALL R12QQ(QQ,WORK(ISMS0),WORK(ISMX1),WORK(ISMY1),WORK(ISMZ1), 
     *           WORK(ISMX2),WORK(ISMY2),WORK(ISMZ2), 
     *           WORK(ISMO2),WORK(ISMO1),WORK(IEND),LQQ,NBAST,NOCCT)
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
     *      CALL QUIT(' Not enough space in UTHZ')
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
      SUBROUTINE QCAL(QQ,R2AM,U2AM,WORK,LWORK,IFLAG)
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
     *CALL QUIT('Core too small in QCAL') 
      CALL QCAL1(R2AM,U2AM,QQ,WORK(IR11),WORK(IQ11), 
     *           WORK(IRS11),WORK(IQS11),WORK(IRT11),
     *           WORK(IQT11),NOCDIM,NSPAIR,IFLAG)
      RETURN
      END
C  /* Deck qcal1 */
      SUBROUTINE QCAL1(R2AM,U2AM,QQ11,R11,Q11,RS11,RT11,QS11,QT11,
     *                 NOCDIM,NSPAIR,IFLAG)
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
C     CONSTRUCT PRODUCTS (IA|JB)*(KA|LB)
C
      DO 101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         DO 102 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 103 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
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
                                    NAIBJ = IU2AM(ISYMIA,ISYMJB) 
     *                                    + NT1AM(ISYMJB)*(NAI - 1) 
     *                                    + NBJ
                                    NAKBL = IT2AM(ISYMKA,ISYMLB) +
     *                                      INDEX(NAK,NBL)
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
C
                                    GOTO 111
C
  114                               CONTINUE
                                    IF (A .LE. NRHFS(ISYMA) .AND.
     *                                  B .LE. NRHFS(ISYMB)) THEN   
C
C                                      Sum |ij><ij|
C
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) - RR
C
                                    ELSE IF (A .LE. NRHFS(ISYMA) .AND.
     *                                       B .GT. NORB1(ISYMB)) THEN
C
C                                      Sum |iq'><iq'|
C
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
C
                                    ELSE IF (A .GT. NORB1(ISYMA) .AND.
     *                                       B .LE. NRHFS(ISYMB)) THEN
C
C                                      Sum |p'j><p'j|
C
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
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
     *                 - (Q11(J,I,L,K) + Q11(J,I,K,L))
                  QS11(KL,IJ) = RR
                  IF (I .EQ. J) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  IF (K .EQ. L) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  RR =   (QQIJKL - QQIJLK) 
     *                 - (Q11(I,J,K,L) - Q11(I,J,L,K))
     *                 - (Q11(J,I,L,K) - Q11(J,I,K,L))
                  QT11(KL,IJ) = RR
                  QS11(KL,IJ) = QS11(KL,IJ) * DP25 * DP5
                  QT11(KL,IJ) = QT11(KL,IJ) * DP75 * DP5
C
                  IF (IFLAG .EQ. 2) GOTO 303
                  RR =   (QQIJKL + QQIJLK)
     *                 + (Q11(I,J,K,L) + Q11(I,J,L,K))
     *                 + (Q11(J,I,L,K) + Q11(J,I,K,L))
     *                 - (R11(I,J,K,L) + R11(I,J,L,K))
     *                 - (R11(J,I,L,K) + R11(J,I,K,L))
                  RS11(KL,IJ) = RR
                  IF (I .EQ. J) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  IF (K .EQ. L) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  RR =   (QQIJKL - QQIJLK)
     *                 + (Q11(I,J,K,L) - Q11(I,J,L,K))
     *                 + (Q11(J,I,L,K) - Q11(J,I,K,L))
     *                 - (R11(I,J,K,L) - R11(I,J,L,K))
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
  998 LUMUL = 43
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
  997 CALL GPOPEN(LUMUL,'QMATSP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QS11
      CALL GPCLOSE(LUMUL,'KEEP')
      CALL GPOPEN(LUMUL,'QMATTP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      WRITE(LUMUL,'(4E30.20)') QT11
      CALL GPCLOSE(LUMUL,'KEEP')
      RETURN
      END
C  /* Deck ccsd_r12 */
      SUBROUTINE CCSD_R12(WORK,LWORK,WXRK,LWXRK)
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
#include <r12int.h>
      CHARACTER*6 TIMLAB
      DIMENSION WORK(LWORK), WXRK(LWXRK)
      INTEGER IDUM, LU43
      LOGICAL PROFIL, LDUM
      LUSIRG = 30
      LU43 = 43
      PROFIL = .FALSE.
      R12CAL = R12CAL .OR. LMULBS
      IF (R12CAL) THEN
         LABEL = 'FULLBAS '
         CALL CCSD_INIT1(WORK,LWORK)
         NORB1T = 0
         NORB2T = 0
         NNBASF = 0
         NSPAIR = NRHFT * (NRHFT + 1) / 2
         R12OLD = .TRUE.
         NOCXBS = 0
         DO 230 ISYM = 1, NSYM
            R12OLD = R12OLD .AND. NORB2(ISYM) .EQ. 0
            NORB1T = NORB1T + NORB1(ISYM)
            NORB2T = NORB2T + NORB2(ISYM)
            NBASF = NORB1(ISYM) + NORB2(ISYM) 
            NNBASF = NNBASF + NBASF * (NBASF + 1) / 2
            NOCXBS = NOCXBS + NRHF(ISYM) * NBAS(ISYM)
  230    CONTINUE
         NORXR = NORXR .OR. R12OLD .OR. R12HYB
C
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
     *     'Number of (i|p'') integrals             : ',NT1AMX
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (i|u'') integrals             : ',NOCXBS
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (Ip''|1/r12|jq'') integrals    : ',NT2FRO(1)
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (ip''|1/r12|jq'') integrals    : ',NT2AMX
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (ip''|r12|jq'') integrals      : ',NT2AMX
         WRITE(LUPRI,'(1X,A,I12)')
     *     'Number of (ip''|[T1,r12]|jq'') integrals : ',NU2AMX
         WRITE(LUPRI,'(1X,A,G13.6)')
     *     'Threshold for R12 contributions            =',VCLTHR
         IF (R12EOR) THEN
            WRITE(LUPRI,'(1X,A,G13.6)')
     *        'Threshold for singular value decomposition =',SVDTHR
            WRITE(LUPRI,'(1X,A,G13.6,A/)')  
     *        'Calculation with Gaussian-damped linear r12 '//
     *        '(gamma =', GAMAC,')'
         ELSE
            WRITE(LUPRI,'(1X,A,G13.6/)')
     *        'Threshold for singular value decomposition =',SVDTHR
         ENDIF
         WRITE(LUPRI,'(1X,78A1)') ('*',I=1,78)
         WRITE(LUPRI,'(1X,A,3(/1X,A))')
     *    'The MP2-R12/A and MP2-R12/B energies'//
     *    ' are computed as described in:',
     *    'W. Klopper and C. C. M. Samson,'//
     *    ' J. Chem. Phys. 116, 6397-6410 (2002).',
     *    '----------------------------------------'//
     *    '--------------------------------------',
     *    'Ansatz 1 is used and results are'//
     *    ' printed for the approximations A, A'' and B.'
         IF (R12HYB) WRITE(LUPRI,'(1X,A,5(/1X,A))')
     *    'However, in the present implementation,'//
     *    ' two-electron four-center integrals',
     *    'that involve two or more auxiliary basis'//
     *    ' functions are not computed, neither',
     *    'in approximations A and A'', nor in'//
     *    ' approximation B. To reproduce the results',
     *    'of the aforementioned article not only in'//
     *    ' approximation A or A'' but also',
     *    'in approximation B, the keyword ''.NO HYB'' must'//
     *    ' be specified in the Dalton',
     *    'input file.'
         IF (R12EOR) WRITE(LUPRI,'(1X,A,3(/1X,A))')
     *    '----------------------------------------'//
     *    '--------------------------------------',
     *    'The Gaussian-damped R12 integrals'//
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
         KGXAM = KT1AM + NT1AMX
         KQQAM = KGXAM + NT1AMX 
         IF (R12EOR) THEN
            KQQ2M = KQQAM + NRHFT ** 4 
            KQQ4M = KQQ2M + NRHFT ** 4 
            KQQ6M = KQQ4M + NRHFT ** 4 
         ELSE
            KQQ2M = KQQAM 
            KQQ4M = KQQAM 
            KQQ6M = KQQAM
         ENDIF
         KUGAM = KQQ6M + NRHFT ** 4 
         KFRIN = KUGAM + NOCXBS 
         KS2AM = KFRIN + NT2FRO(1)
         KT2AM = KS2AM + NT2AMX
         KU2AM = KT2AM 
         KR2AM = KT2AM + NT2AMX  
         KV2AM = MAX(KR2AM + NT2AMX, KU2AM + NU2AMX)
         KW2AM = KV2AM + NT2AMX 
         KEND1 = KW2AM + NT2AMX
         LWRK1 = LWORK - KEND1
         IF (KEND1 .GT. LWORK) CALL QUIT('R12 WORK SPACE')
         IF (R12EOR) THEN
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
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
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
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
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
            MBSMAX = 4
            INTGAC = 2
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),WORK(KEND1),LWRK1) 
            CALL CCSD_IKJL(WORK(KQQ2M),NRHFT,WXRK(KS2AM))
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ik|f12/r12|jl) integrals done' 
            WRITE(TIMLAB,1013) '#2 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
            R12EIN = .FALSE.
         ELSE
            CALL TIMER('START ',TIMSTR,TIMEND)  
         END IF
C
C        COMPUTE (IA|1/R12|JB) INTEGRALS
C
         R12TRA = .TRUE.
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
         CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
         CALL GPOPEN(LU43,'R12V12',
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         WRITE(LU43) (WXRK(KS2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LU43,'KEEP')  
         IF (R12ECO) GOTO 1094 
         IF (FROIMP) THEN
            CALL GPOPEN(LU43,'INCJDK',
     &                       'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
            REWIND(LU43)
            READ(LU43) (WORK(KFRIN+I-1), I = 1,NT2FRO(1))
            CALL GPCLOSE(LU43,'KEEP')  
         END IF
         CALL MAKEGX(WXRK(KS2AM),WORK(KGXAM),WORK(KFRIN),2)
         CALL MAKEUG(WORK(KGXAM),WORK(KUGAM),
     &               WORK(KQQAM),WORK(KEND1),LWRK1,2) 
         IF (R12HYB) THEN
            CALL MAKEGX(WXRK(KS2AM),WORK(KGXAM),WORK(KFRIN),1)
            CALL MAKEUG(WORK(KGXAM),WORK(KUGAM),
     &                  WORK(KQQAM),WORK(KEND1),LWRK1,1) 
         END IF
 1094    CONTINUE
         IF (R12EOR) THEN
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|1/r12|jq'') integrals done' 
            WRITE(TIMLAB,1013) '#1 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|1/r12|jq'') integrals done' 
            WRITE(TIMLAB,1013) '1/R',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
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
         ELSE IF (R12NOB .OR. NORXR) THEN
            MBSMAX = 5
         ELSE
            MBSMAX = 6
         END IF
         CALL DZERO(WORK(KT1AM),NT1AMX)
         CALL CCSD_IAJB(WXRK(KS2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
         CALL GPOPEN(LU43,'R12R12',
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         WRITE(LU43) (WXRK(KS2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LU43,'KEEP')  
         IF (R12EOR) THEN
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|f12|jq'') integrals done' 
            WRITE(TIMLAB,1013) '#3 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|r12|jq'') integrals done' 
            WRITE(TIMLAB,1013) 'R  ',MBSMAX
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
C        COMPUTE (I'A|R12|JB) INTEGRALS
C     
 1095    CONTINUE
         R12TRA = .TRUE.
         IF (R12EOR) THEN
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
            IF (.NOT. R12NOB) THEN
               INTGAC = 4       
               IF (R12OLD .OR. R12HYB) THEN
                  MBSMAX = 4
               ELSE
                  MBSMAX = 5
               END IF
               CALL DZERO(WORK(KT1AM),NT1AMX)
               CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
               CALL CCSD_IKJL(WORK(KQQAM),NRHFT,WXRK(KU2AM))
               WRITE(TIMLAB,1013) '#4 ',MBSMAX
               CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
            END IF
            INTGAC = 3       
         END IF
         IF (R12OLD) THEN
            IF (R12NOB) GOTO 1097
            MBSMAX = 4
         ELSE IF (R12NOB .OR. R12HYB) THEN
            MBSMAX = 5
         ELSE
            MBSMAX = 6
         END IF
         CALL DZERO(WORK(KT1AM),NT1AMX)
         CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
         CALL QCAL(WORK(KQQAM),WXRK(KS2AM),WXRK(KU2AM),
     *             WORK(KEND1),LWRK1,1)
         IF (.NOT. R12OLD) THEN
            CALL QCAL(WORK(KQQAM),WXRK(KS2AM),WXRK(KU2AM),
     *                WORK(KEND1),LWRK1,2)
            IF (R12HYB) THEN
               IF (R12EOR) THEN
                  WRITE(TIMLAB,1013) '#3 ',MBSMAX
                  CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
               ELSE
                  WRITE(TIMLAB,1013) 'R  ',MBSMAX
                  CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
               END IF
               COMBSS = .TRUE.  
               CALL DZERO(WORK(KT1AM),NT1AMX)
               CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
            END IF
            CALL COMKR(WXRK(KS2AM),WXRK(KU2AM),WORK(KV2AM),
     *                 WORK(KW2AM),WORK(KEND1),LWRK1)
         END IF
         IF (R12EOR) THEN
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (i''p''|f12|jq'') integrals done' 
            WRITE(TIMLAB,1013) '#3 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (i''p''|r12|jq'') integrals done' 
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
         CALL DZERO(WORK(KT1AM),NT1AMX)
         CALL CCSD_IAJB(WXRK(KU2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
         CALL TESTUU(WXRK(KU2AM),WXRK(KS2AM))
         IF (R12EOR) THEN
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|[T,f12]|jq'') integrals done' 
            WRITE(TIMLAB,1013) '#5 ',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         ELSE
            WRITE(LUPRI,'(1X,A)')
     *        'Computation of (ip''|[T,r12]|jq'') integrals done' 
            WRITE(TIMLAB,1013) 'T,R',MBSMAX
            CALL TIMER(TIMLAB,TIMSTR,TIMEND)  
         END IF
C
 1096    CONTINUE
         CALL GPOPEN(LU43,'R12V12',
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         READ(LU43) (WXRK(KT2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LU43,'DELETE')  
         CALL GPOPEN(LU43,'R12R12',
     &                    'UNKNOWN',' ','UNFORMATTED',IDUM,LDUM)
         READ(LU43) (WXRK(KR2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LU43,'DELETE')  
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
                           IF (R12SQR) THEN
                              NKILJ = IU2AM(ISYMIK,ISYMJL)
     *                              + NT1AM(ISYMJL)*(NKI - 1)
     *                              + NLJ 
                           ELSE
                              NKILJ = IT2AM(ISYMIK,ISYMJL)
     *                              + INDEX(NKI,NLJ)  
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
      SUBROUTINE R12MP2(VS,VT,SS,ST,EKL,NOCDIM,NSPAIR,FS,FT,
     *                  BB,UU,WW,QQ,SF,TF,EV,XF,VV,IJPAIR,WS,WT)
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
        END DO
      END IF
      CALL DZERO(WW,NSPAIR)
      CALL DZERO(UU,NSPAIR*NSPAIR)
      CALL DZERO(VV,NSPAIR*NSPAIR)
      CALL DZERO(QQ,NSPAIR*NSPAIR)
      CALL SVD(NSPAIR,NSPAIR,NSPAIR,BB,WW,MATU,UU,MATV,VV,IERR,QQ)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,'(//A,2I10/)') ' SVD ERROR / IJPAIR, IERR =',
     *                                            IJPAIR, IERR
         CALL QUIT('SVD ERROR -- SINGLET')
      END IF
      DO I = 1, NSPAIR 
         IF (ABS(WW(I)) .LT. SVDTHR) WW(I) = D0
         IF (WW(I) .NE. D0) THEN
            WS = MIN(WS, WW(I))
            IF (WW(I) .LT. VCLTHR) WW(I) = D0
         END IF 
      END DO 
      CALL SVBKSB(UU,WW,VV,NSPAIR,NSPAIR,NSPAIR,NSPAIR,VS,BB,QQ)
      FS = -DDOT(NSPAIR,VS,1,BB,1)
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
        END DO
      END IF
      CALL DZERO(WW,NSPAIR)
      CALL DZERO(UU,NSPAIR*NSPAIR)
      CALL DZERO(VV,NSPAIR*NSPAIR)
      CALL DZERO(QQ,NSPAIR*NSPAIR)
      CALL SVD(NSPAIR,NSPAIR,NSPAIR,BB,WW,MATU,UU,MATV,VV,IERR,QQ)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,'(//A,2I10/)') ' SVD ERROR / IJPAIR, IERR =',
     *                                            IJPAIR, IERR
         CALL QUIT('SVD ERROR -- TRIPLET')
      END IF
      DO I = 1, NSPAIR 
         IF (ABS(WW(I)) .LT. SVDTHR) WW(I) = D0
         IF (WW(I) .NE. D0) THEN
            WT = MIN(WT, WW(I))
            IF (WW(I) .LT. VCLTHR) WW(I) = D0
         END IF 
      END DO 
      CALL SVBKSB(UU,WW,VV,NSPAIR,NSPAIR,NSPAIR,NSPAIR,VT,BB,QQ)
      FT = -DDOT(NSPAIR,VT,1,BB,1)
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
     *                  CNINV,SF)
C
C     This subroutine computes the V(klmn), X(klmn), and B(klmn)
C     matrices and evaluates the MP2-R12 energies.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
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
      LUMULB = 34
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
      CALL AROUND('Singlet <Bab> matrix') 
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
            CALL R12MP2(VS11(1,IJ),VT11(1,IJ),RS11,RT11,EV(I)+EV(J),
     *      NOCDIM,NSPAIR,FS(IJ),FT(IJ),WS11,WT11,QS11,QT11,
     *      BS11,BT11,EVS,D1,V11,IJ,WSMIN,WTMIN)
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
      LU33 = 33
      LU34 = 34
      LU35 = 35
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
