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
C  /* Deck getsos */
      SUBROUTINE GETSOS(GSO,NCNT,CORPX,CORPY,CORPZ,
     &                  WORK,LWORK,NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
C
C     T. Helgaker feb 01
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
      LOGICAL DOLND, DOGGA
      DIMENSION GSO(NBAST*NTYPSO), WORK(LWORK), NCNT(NBAST)
C
#include <dftinf.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (MAXREP.EQ.0) THEN
         CALL DFTAOS(GSO(KSO0),GSO(KSO1),GSO(KSO2),GSO(KSOB),GSO(KSOB1),
     &               NCNT,CORPX,CORPY,CORPZ,WORK,LWORK,NBAST,DOLND,
     &               DOGGA,DFTHRI,IPRINT)
      ELSE
         KGAO = 1
         KLST = KGAO + NTYPSO*NBAST
         IF (KLST.GT.LWORK) CALL STOPIT('DFTAOS','LWORK',KLST,LWORK)
         CALL DFTAOS(WORK(KSO0),WORK(KSO1),WORK(KSO2),WORK(KSOB),
     &               WORK(KSOB1),NCNT,CORPX,CORPY,CORPZ,WORK,LWORK,
     &               NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
         CALL DFTSOS(WORK(KSO0),GSO,NBAST,NTYPSO,NCNT,IPRINT)
      END IF
C
      RETURN
      END      
C  /* Deck getaos */
      SUBROUTINE DFTAOS(GAO,GAO1,GAO2,GAB1,GAB2,NCNT,CORPX,CORPY,CORPZ,
     &                  WORK,LWORK,NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
C
C     T. Helgaker sep 99
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <aovec.h>
C
      PARAMETER (D0 = 0.0D0, DTHRS = 20.0D0)
C   
      LOGICAL DOLND, DOGGA
      DIMENSION GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6),
     &          GAB1(NBAST,3), GAB2(NBAST,3,3),
     &          WORK(LWORK), NCNT(NBAST)
C
#include <onecom.h>
#include <lmns.h>
#include <nuclei.h>
#include <shells.h>
#include <symmet.h>
#include <primit.h>
#include <sphtrm.h>
#include <orgcom.h>
#include <dftinf.h>
C
#include <ibtfun.h>
C
      CALL DZERO(GAO,NTYPSO*NBAST)
C
      IADR = 1
      DO ISHELA = 1,KMAX
         NHKTA  = NHKT(ISHELA)
         KHKTA  = KHKT(ISHELA)
         KCKTA  = KCKT(ISHELA)
         SPHRA  = SPHR(ISHELA)
         NUMCFA = NUMCF(ISHELA)
         JSTA   = JSTRT(ISHELA)
         NUCA   = NUCO(ISHELA)
         MULA   = ISTBAO(ISHELA)
         DO ISYMOP = 0, MAXOPR
         IF (IBTAND(ISYMOP,MULA).EQ.0) THEN
C
            PAX = CORPX-PT(IBTAND(ISYMAX(1,1),ISYMOP))*CENT(ISHELA,1,1)
            PAY = CORPY-PT(IBTAND(ISYMAX(2,1),ISYMOP))*CENT(ISHELA,2,1)
            PAZ = CORPZ-PT(IBTAND(ISYMAX(3,1),ISYMOP))*CENT(ISHELA,3,1)
            PA2 = PAX**2 + PAY**2 + PAZ**2
C
C           This test has been commented out since in somce cases the 
C           last exponent may not be the smallest
C
C           FAC = PRIEXP(JSTA + NUCA)*PA2
C           IF (FAC.LT.DTHRS) THEN
               IF (NDER.GT.1 .OR. NHKTA.GT.2 .OR. DOLND) THEN
                  CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
               END IF
               IF (NDER.EQ.0) THEN
                  CALL GETGAO(GAO(IADR),CSP(ISPADR(NHKTA)),
     &                        PAX,PAY,PAZ,PA2,DFTHRI)
               ELSE IF (NDER.GT.0) THEN
                  CALL GETGA1(GAO(IADR),GAO1(IADR,1),GAO1(IADR,2),
     &                        GAO1(IADR,3),CSP(ISPADR(NHKTA)),
     &                        PAX,PAY,PAZ,PA2,DFTHRI)
                  IF (NDER.GT.1) THEN
                     CALL GETGA2(GAO2(IADR,1),GAO2(IADR,2),GAO2(IADR,3),
     &                           GAO2(IADR,4),GAO2(IADR,5),GAO2(IADR,6),
     &                           CSP(ISPADR(NHKTA)),PAX,PAY,PAZ,PA2,
     &                           DFTHRI)
                  END IF
               END IF
               IF (DOLND) THEN
                  CALL GETGB1(GAO(IADR),GAB1(IADR,1),GAB1(IADR,2),
     &                        GAB1(IADR,3),PAX,PAY,PAZ)
                  IF (DOGGA) THEN
                     DO I = 1, 3
                        CALL GETGB1(GAO1(IADR,I),GAB2(IADR,I,1),
     &                              GAB2(IADR,I,2),GAB2(IADR,I,3),
     &                              PAX,PAY,PAZ)
                     END DO
                  END IF
               END IF
C           END IF
C
C           IF (NDER.GT.0 .OR. DOLND) THEN
C              NCNTA = NCENT(ISHELA)
C              DO ICMP = 0, KHKTA - 1
C                 NCNT(IADR + ICMP) = NCNTA
C              END DO
C           END IF
            IADR = IADR + KHKTA
         END IF
         END DO
      END DO
C
C     Print section
C
      IF (IPRINT .GE. 50) THEN
          CALL HEADER('Output from DFTAOS',-1)
          WRITE (LUPRI,'(A,3F12.6)') ' CORP ', CORPX,CORPY,CORPZ 
          IF (IPRINT .GT. 200) THEN
             CALL HEADER('Undifferentiated integrals at this point:',-1)
             CALL OUTPUT(GAO,1,NBAST,1,1,NBAST,1,1,LUPRI)
             IF (NDER.GT.0) THEN
               CALL HEADER('1st derivative integrals at this point:',-1)
               CALL OUTPUT(GAO1,1,NBAST,1,3,NBAST,3,1,LUPRI)
             END IF
             IF (NDER.GT.1) THEN
               CALL HEADER('2nd derivative integrals at this point:',-1)
               CALL OUTPUT(GAO2,1,NBAST,1,6,NBAST,6,1,LUPRI)
             END IF
          END IF
          WRITE (LUPRI,'(1X,A)') ' '
      END IF
      RETURN
      END
C  /* Deck getgao */
      SUBROUTINE GETGAO(GAO,CSP,PAX,PAY,PAZ,PA2,DFTHRI)
C
C     T. Helgaker Sep 99
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
      DIMENSION GAO(KCKTA)
      DIMENSION CSP(KHKTA,KCKTA)
C
#include <lmns.h>
#include <onecom.h>
#include <primit.h>
C
C     loop over primitives
C
      GA = PRICCF(JSTA + 1,NUMCFA)*DEXP(-PRIEXP(JSTA + 1)*PA2)
      DO I = JSTA + 2, JSTA + NUCA 
         GA = GA + PRICCF(I,NUMCFA)*DEXP(-PRIEXP(I)*PA2)
      END DO
C
C     contracted orbitals
C
      IF (DABS(GA).GT.DFTHRI) THEN
         IF (NHKTA .EQ. 1) THEN
            GAO(1) = GA
         ELSE IF (NHKTA .EQ. 2) THEN
            GAO(1) = PAX*GA
            GAO(2) = PAY*GA
            GAO(3) = PAZ*GA
         ELSE IF (NHKTA .EQ. 3) THEN
            IF (SPHRA) THEN
               GAX  = PAX*GA
               GAY  = PAY*GA
               GAZ  = PAZ*GA
               GAXX = PAX*GAX
               GAYY = PAY*GAY
               GAO(1) = CSP(1,2)*PAY*GAX
               GAO(2) = CSP(2,5)*PAY*GAZ
               GAO(3) = CSP(3,1)*GAXX + CSP(3,4)*GAYY + CSP(3,6)*PAZ*GAZ
               GAO(4) = CSP(4,3)*PAX*GAZ
               GAO(5) = CSP(5,1)*GAXX + CSP(5,4)*GAYY
            ELSE
               GAX  = PAX*GA
               GAY  = PAY*GA
               GAZ  = PAZ*GA
               GAO(1) = PAX*GAX
               GAO(2) = PAY*GAX
               GAO(3) = PAZ*GAX
               GAO(4) = PAY*GAY
               GAO(5) = PAZ*GAY 
               GAO(6) = PAZ*GAZ
            END IF
         ELSE
            IF (SPHRA) THEN
               DO I = 1, KHKTA
                  GAO(I) = D0
               END DO
               DO J = 1, KCKTA
                  CINT = (PAX**LVALUA(J))*(PAY**MVALUA(J))
     &                                   *(PAZ**NVALUA(J))*GA
                  IF (DABS(CINT).GT.DFTHRI) THEN
                     DO I = 1, KHKTA
                        SPHFAC = CSP(I,J)
                        IF (ABS(SPHFAC).GT.D0) THEN
                           GAO(I) = GAO(I) + SPHFAC*CINT
                        END IF
                     END DO
                  END IF 
               END DO
            ELSE
               DO I = 1, KCKTA
                  GAO(I) = (PAX**LVALUA(I))*(PAY**MVALUA(I))*
     &                     (PAZ**NVALUA(I))*GA
               END DO
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck getga1 */
      SUBROUTINE GETGA1(GAO,GAOX,GAOY,GAOZ,CSP,PAX,PAY,PAZ,PA2,DFTHRI)
C
C     T. Helgaker Sep 99
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
#include <maxaqn.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0) 
C
      DIMENSION GAO(KCKTA), GAOX(KCKTA), GAOY(KCKTA), GAOZ(KCKTA),
     &          CAO(MXAQN), CAOX(MXAQN), CAOY(MXAQN), CAOZ(MXAQN)
      DIMENSION CSP(KHKTA,KCKTA)
C
#include <lmns.h>
#include <onecom.h>
#include <primit.h>
C
#include <sdpre.h>
C
C     loop over primitives
C
      GA = D0
      GU = D0
      DO I = JSTA + 1, JSTA + NUCA
         ALPHA = PRIEXP(I)
         FAC   = PRICCF(I,NUMCFA)*DEXP(-ALPHA*PA2)
         GA    = GA + FAC
         GU    = GU + ALPHA*FAC
      END DO 
      GU = -D2*GU
C
      IF (DABS(GA).GT.DFTHRI) THEN 
C
C        s orbitals
C
         IF (NHKTA .EQ. 1) THEN
            GAO(1)  = GA
            GAOX(1) = PAX*GU
            GAOY(1) = PAY*GU
            GAOZ(1) = PAZ*GU
C
C        p orbitals
C
         ELSE IF (NHKTA .EQ. 2) THEN
            TGX = PAX*GU
            TGY = PAY*GU
            TGZ = PAZ*GU
C
            GAO(1)  = PAX*GA
            GAO(2)  = PAY*GA
            GAO(3)  = PAZ*GA
            GAOX(1) = PAX*TGX + GA
            GAOX(2) = PAY*TGX
            GAOX(3) = PAZ*TGX
            GAOY(1) = PAX*TGY
            GAOY(2) = PAY*TGY + GA
            GAOY(3) = PAZ*TGY
            GAOZ(1) = PAX*TGZ
            GAOZ(2) = PAY*TGZ
            GAOZ(3) = PAZ*TGZ + GA
C
C        d and higher orbitals
C
         ELSE
            FX = PAX*GU
            FY = PAY*GU
            FZ = PAZ*GU
            DO ICOMPA = 1,KCKTA
               L = LVALUA(ICOMPA)
               M = MVALUA(ICOMPA)
               N = NVALUA(ICOMPA)
               P0  = (PAX**L)*(PAY**M)*(PAZ**N)
               GAX = FX*P0
               GAY = FY*P0
               GAZ = FZ*P0
               IF(L.GT.0) GAX = GAX
     &                    + SDPRE(L)*(PAX**(L-1))*(PAY**M)*(PAZ**N)*GA
               IF(M.GT.0) GAY = GAY
     &                    + SDPRE(M)*(PAX**L)*(PAY**(M-1))*(PAZ**N)*GA
               IF(N.GT.0) GAZ = GAZ
     &                    + SDPRE(N)*(PAX**L)*(PAY**M)*(PAZ**(N-1))*GA
               IF (SPHRA) THEN
                  CAO (ICOMPA) = GA*P0
                  CAOX(ICOMPA) = GAX
                  CAOY(ICOMPA) = GAY
                  CAOZ(ICOMPA) = GAZ
               ELSE
                  GAO (ICOMPA) = GA*P0
                  GAOX(ICOMPA) = GAX
                  GAOY(ICOMPA) = GAY
                  GAOZ(ICOMPA) = GAZ
               END IF
            END DO
            IF(SPHRA) THEN
                IF (NHKTA.EQ.3) THEN
                   GAO(1) = CSP(1,2)*CAO(2)
                   GAO(2) = CSP(2,5)*CAO(5)
                   GAO(3) = CSP(3,1)*CAO(1) + CSP(3,4)*CAO(4) 
     &                                      + CSP(3,6)*CAO(6)
                   GAO(4) = CSP(4,3)*CAO(3)
                   GAO(5) = CSP(5,1)*CAO(1) + CSP(5,4)*CAO(4)
                   GAOX(1) = CSP(1,2)*CAOX(2)
                   GAOX(2) = CSP(2,5)*CAOX(5)
                   GAOX(3) = CSP(3,1)*CAOX(1) + CSP(3,4)*CAOX(4) 
     &                                        + CSP(3,6)*CAOX(6)
                   GAOX(4) = CSP(4,3)*CAOX(3)
                   GAOX(5) = CSP(5,1)*CAOX(1) + CSP(5,4)*CAOX(4)
                   GAOY(1) = CSP(1,2)*CAOY(2)
                   GAOY(2) = CSP(2,5)*CAOY(5)
                   GAOY(3) = CSP(3,1)*CAOY(1) + CSP(3,4)*CAOY(4) 
     &                                        + CSP(3,6)*CAOY(6)
                   GAOY(4) = CSP(4,3)*CAOY(3)
                   GAOY(5) = CSP(5,1)*CAOY(1) + CSP(5,4)*CAOY(4)
                   GAOZ(1) = CSP(1,2)*CAOZ(2)
                   GAOZ(2) = CSP(2,5)*CAOZ(5)
                   GAOZ(3) = CSP(3,1)*CAOZ(1) + CSP(3,4)*CAOZ(4) 
     &                                        + CSP(3,6)*CAOZ(6)
                   GAOZ(4) = CSP(4,3)*CAOZ(3)
                   GAOZ(5) = CSP(5,1)*CAOZ(1) + CSP(5,4)*CAOZ(4)
                ELSE
                   DO I = 1, KHKTA
                      SPH0 = D0 
                      SPHX = D0 
                      SPHY = D0 
                      SPHZ = D0 
                      DO J = 1, KCKTA
                         SPHFAC = CSP(I,J)
                         IF (ABS(SPHFAC).GT.D0) THEN
                            SPH0 = SPH0 + SPHFAC*CAO (J)
                            SPHX = SPHX + SPHFAC*CAOX(J)
                            SPHY = SPHY + SPHFAC*CAOY(J)
                            SPHZ = SPHZ + SPHFAC*CAOZ(J)
                         END IF
                      END DO
                      GAO (I) = SPH0
                      GAOX(I) = SPHX
                      GAOY(I) = SPHY
                      GAOZ(I) = SPHZ
                   END DO
                END IF
            END IF
         END IF
      END IF
      END
C  /* Deck getga2 */
      SUBROUTINE GETGA2(GAOXX,GAOXY,GAOXZ,GAOYY,GAOYZ,GAOZZ,CSP,
     &                  PAX,PAY,PAZ,PA2,DFTHRI)
C
C     T. Helgaker Sep 99
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0) 
C
      DIMENSION GAOXX(KCKTA), GAOXY(KCKTA), GAOXZ(KCKTA),
     &          GAOYY(KCKTA), GAOYZ(KCKTA), GAOZZ(KCKTA),
     &          CAOXX(MXAQN), CAOXY(MXAQN), CAOXZ(MXAQN),
     &          CAOYY(MXAQN), CAOYZ(MXAQN), CAOZZ(MXAQN)
      DIMENSION CSP(KHKTA,KCKTA)
C
#include <lmns.h>
#include <onecom.h>
#include <primit.h>
C
#include <sdpre.h>
C
      IF (SPHRA) THEN
         DO I=1,KCKTA
            CAOXX(I) = D0
            CAOXY(I) = D0
            CAOXZ(I) = D0          
            CAOYY(I) = D0
            CAOYZ(I) = D0
            CAOZZ(I) = D0
         END DO
      END IF
      DO IPRIMA = JSTA + 1, JSTA + NUCA
         ALPHA = PRIEXP(IPRIMA)
         TALPH = -D2*ALPHA
         TAPAX = TALPH*PAX
         TAPAY = TALPH*PAY
         TAPAZ = TALPH*PAZ
         GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-ALPHA*PA2)
         IF (DABS(GA).GT.DFTHRI) THEN
            DO ICOMPA = 1, KCKTA
               L = LVALUA(ICOMPA)
               M = MVALUA(ICOMPA)
               N = NVALUA(ICOMPA)
C
               PXD = D0
               PYD = D0
               PZD = D0
               IF (L.GT.1) PXD = SDPRE(L*(L-1))*(PAX**(L-2))
               IF (M.GT.1) PYD = SDPRE(M*(M-1))*(PAY**(M-2))
               IF (N.GT.1) PZD = SDPRE(N*(N-1))*(PAZ**(N-2))
               PXM = D0
               PYM = D0
               PZM = D0
               IF (L.GT.0) PXM = SDPRE(L)*(PAX**(L-1))
               IF (M.GT.0) PYM = SDPRE(M)*(PAY**(M-1))
               IF (N.GT.0) PZM = SDPRE(N)*(PAZ**(N-1))
               PX0 = PAX**L
               PY0 = PAY**M
               PZ0 = PAZ**N
               PXP = TAPAX*PX0
               PYP = TAPAY*PY0
               PZP = TAPAZ*PZ0
               P000 = PX0*PY0*PZ0 
C
               GAXX = (TAPAX**2 + TALPH*SDPRE(2*L+1))*P000 + PXD*PY0*PZ0
               GAYY = (TAPAY**2 + TALPH*SDPRE(2*M+1))*P000 + PX0*PYD*PZ0
               GAZZ = (TAPAZ**2 + TALPH*SDPRE(2*N+1))*P000 + PX0*PY0*PZD
               GAXY = TAPAX*TAPAY*P000 + (PXP*PYM+PXM*PYP+PXM*PYM)*PZ0
               GAXZ = TAPAX*TAPAZ*P000 + (PXP*PZM+PXM*PZP+PXM*PZM)*PY0
               GAYZ = TAPAY*TAPAZ*P000 + (PYP*PZM+PYM*PZP+PYM*PZM)*PX0
C
               IF (SPHRA) THEN
                  CAOXX(ICOMPA) = CAOXX(ICOMPA) + GAXX*GA 
                  CAOXY(ICOMPA) = CAOXY(ICOMPA) + GAXY*GA 
                  CAOXZ(ICOMPA) = CAOXZ(ICOMPA) + GAXZ*GA
                  CAOYY(ICOMPA) = CAOYY(ICOMPA) + GAYY*GA
                  CAOYZ(ICOMPA) = CAOYZ(ICOMPA) + GAYZ*GA
                  CAOZZ(ICOMPA) = CAOZZ(ICOMPA) + GAZZ*GA
               ELSE
                  GAOXX(ICOMPA) = GAOXX(ICOMPA) + GAXX*GA 
                  GAOXY(ICOMPA) = GAOXY(ICOMPA) + GAXY*GA 
                  GAOXZ(ICOMPA) = GAOXZ(ICOMPA) + GAXZ*GA
                  GAOYY(ICOMPA) = GAOYY(ICOMPA) + GAYY*GA
                  GAOYZ(ICOMPA) = GAOYZ(ICOMPA) + GAYZ*GA
                  GAOZZ(ICOMPA) = GAOZZ(ICOMPA) + GAZZ*GA
               END IF
            END DO
         END IF
      END DO
      IF (SPHRA) THEN
         DO I = 1, KHKTA
            SPHXX = D0 
            SPHXY = D0 
            SPHXZ = D0 
            SPHYY = D0 
            SPHYZ = D0 
            SPHZZ = D0 
            DO J = 1, KCKTA
               SPHFAC = CSP(I,J)
               IF (ABS(SPHFAC).GT.D0) THEN
                  SPHXX = SPHXX + SPHFAC*CAOXX(J)
                  SPHXY = SPHXY + SPHFAC*CAOXY(J)
                  SPHXZ = SPHXZ + SPHFAC*CAOXZ(J)
                  SPHYY = SPHYY + SPHFAC*CAOYY(J)
                  SPHYZ = SPHYZ + SPHFAC*CAOYZ(J)
                  SPHZZ = SPHZZ + SPHFAC*CAOZZ(J)
               END IF
            END DO
            GAOXX(I) = SPHXX
            GAOXY(I) = SPHXY
            GAOXZ(I) = SPHXZ
            GAOYY(I) = SPHYY
            GAOYZ(I) = SPHYZ
            GAOZZ(I) = SPHZZ
         END DO
      END IF
      RETURN
      END
C  /* Deck getgb1 */
      SUBROUTINE GETGB1(GAO,GAOX,GAOY,GAOZ,PAX,PAY,PAZ)
C
C     T. Helgaker Oct 99
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
#include <maxaqn.h>
C
      PARAMETER (DP5 = 0.5D0)
      DIMENSION GAOX(KHKTA), GAOY(KHKTA), GAOZ(KHKTA), GAO(KHKTA)
C
#include <onecom.h>
#include <primit.h>
#include <orgcom.h>
C
      FCX = DP5*PAX
      FCY = DP5*PAY
      FCZ = DP5*PAZ
      DO ICOMPA = 1,KHKTA
         GA = GAO(ICOMPA)
         GAOX(ICOMPA) = FCX*GA
         GAOY(ICOMPA) = FCY*GA
         GAOZ(ICOMPA) = FCZ*GA
      END DO
      RETURN
      END
C /* Deck dftsph */
      SUBROUTINE DFTSPH(CARINT,SPHINT,CSP,KHKTA,KCKTA)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION CARINT(KCKTA), SPHINT(KHKTA)
      DIMENSION CSP(KHKTA,KCKTA)
C
      DO I = 1, KHKTA
         SPH = D0 
         DO J = 1, KCKTA
            SPHFAC = CSP(I,J)
            IF (ABS(SPHFAC).GT.D0) SPH = SPH + SPHFAC*CARINT(J)
         END DO
         SPHINT(I) = SPH 
      END DO
      RETURN
      END
C  /* Deck dftsos */
      SUBROUTINE DFTSOS(GAO,GSO,NBAST,NVEC,NCNT,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION GSO(NBAST,NVEC), GAO(NBAST,NVEC), NCNT(NBAST)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DFTSOS',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(GSO,NBAST*NVEC)
      DO IREPA = 0, MAXREP
         NORBA = NAOS(IREPA+1)
         DO I = ISTRA,ISTRA + NORBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA = NHKT(IA)
            KHKTA = KHKT(IA)
            MULA  = ISTBAO(IA)
            INDA  = KSTRT(IA) + NA - KHKTA
            DO ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               CALL DAXPY(NVEC,FAC,GAO(INDA,1),NBAST,GSO(I,1),NBAST)
            END IF
            END DO
            NCNT(I) = NCENT(IA)
         END DO
         ISTRA = ISTRA + NORBA
      END DO
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('AOs in DFTSOS',-1)
         CALL OUTPUT(GAO,1,NBAST,1,NVEC,NBAST,NVEC,1,LUPRI)
         CALL HEADER('SOs in DFTSOS',-1)
         CALL OUTPUT(GSO,1,NBAST,1,NVEC,NBAST,NVEC,1,LUPRI)
      END IF
      RETURN
      END
