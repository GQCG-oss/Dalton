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
C  /* Deck getsos */
      SUBROUTINE GETSOS(GSO,NCNT,COR,WORK,LWORK,NBAST,
     &                  DOLND,DOGGA,DFTHRI,IPRINT)
C
C     T. Helgaker feb 01
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
C
      LOGICAL DOLND, DOGGA
      DIMENSION GSO(NBAST*NTYPSO), WORK(LWORK), NCNT(NBAST),COR(3)
C
#include <dftinf.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (MAXREP.EQ.0) THEN
         CALL DFTAOS(GSO(KSO0),GSO(KSO1),GSO(KSO2),GSO(KSOB),GSO(KSOB1),
     &               NCNT,COR(1),COR(2),COR(3),WORK,LWORK,NBAST,DOLND,
     &               DOGGA,DFTHRI,IPRINT)
      ELSE
         KGAO = 1
         KLST = KGAO + NTYPSO*NBAST
         IF (KLST.GT.LWORK) CALL STOPIT('DFTAOS','LWORK',KLST,LWORK)
         CALL DFTAOS(WORK(KSO0),WORK(KSO1),WORK(KSO2),WORK(KSOB),
     &               WORK(KSOB1),NCNT,COR(1),COR(2),COR(3),WORK,LWORK,
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
#include <maxorb.h>
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
            PA2 = PAX*PAX + PAY*PAY + PAZ*PAZ
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
#include <maxorb.h>
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
#include <maxorb.h>
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
#include <maxorb.h>
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
#include <maxorb.h>
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
#include <maxorb.h>
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
               DO J=1, NVEC
                  GSO(I,J) = GSO(I,J) + FAC*GAO(INDA,J)
               END DO
            END IF
            END DO
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
c ===================================================================
c BLOCKED VERSION OF ORBTIAL EVALUATION ROUTINES
c Written by Pawel Salek, closely based on the above.
c ===================================================================
c RETURNS: GSO: evaluated orbitals for a batch of grid points.
c     GSO(:,:,1) contains orbital values.
c     GSO(:,:,2:4) contains first geom. derivatives. - if requested.
c     GSO(:,:,5:10) contains second derivatives - if requested.
c     After requested geometric derivatives, LND related derivatives
c     are placed.
C  /* Deck blgetsos */
      SUBROUTINE BLGETSOS(NVCLEN,GSO,COOR,NBLCNT,IBLCKS,WORK,LWORK,
     &                    NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
C
C     T. Helgaker feb 01, P. Salek 03
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
C
      LOGICAL DOLND, DOGGA
      DIMENSION GSO(NVCLEN*NBAST*NTYPSO), WORK(LWORK)
      DIMENSION COOR(3,NVCLEN)
      DIMENSION IBLCKS(2,NBLCNT)
C
#include <symmet.h>
#include <dftinf.h>
C
#include <ibtfun.h>
C
      IF (MAXREP.EQ.0) THEN
         CALL BLDFTAOS(NVCLEN,GSO,COOR,NBLCNT,IBLCKS,
     &                 WORK,LWORK,NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
      ELSE
         KGAO = 1
         KLST = KGAO + NVCLEN*NTYPSO*NBAST
         IF (KLST.GT.LWORK) CALL STOPIT('BLDFTAO','LWORK',KLST,LWORK)
         LWRK = LWORK - KLST +1
         CALL BLDFTAOS(NVCLEN,WORK(KGAO),COOR,NBLCNT,
     &                 IBLCKS,WORK(KLST),LWRK,NBAST,DOLND,DOGGA,DFTHRI,
     &                 IPRINT)
         CALL BLDFTSOS(NVCLEN,WORK(KGAO),GSO,NBAST,NTYPSO,IPRINT)
      END IF
C
      RETURN
      END      
C  /* Deck getaos */
      SUBROUTINE BLDFTAOS(NVCLEN,GAO,COOR,NBLCNT,IBLCKS,WORK,LWORK,
     &                    NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
C
C     T. Helgaker sep 99, P. Salek 03
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <aovec.h>
C
      PARAMETER (D0 = 0.0D0, DTHRS = 20.0D0)
C   
      LOGICAL DOLND, DOGGA
      DIMENSION GAO(NVCLEN,NBAST,NTYPSO),WORK(LWORK)
      DIMENSION COOR(3,NVCLEN)
      DIMENSION IBLCKS(2,NBLCNT)
C     PA2 contains distance from the basis function center to respective
c     grid point.
      DIMENSION PA(3,NVCLEN), PA2(NVCLEN)
c
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
      IF(NDER.GT.1) THEN
c        fixme - this won't scale linearly 
         CALL DZERO(GAO(1,1,5),6*NBAST*NVCLEN)
      END IF

      DO IBL = 1, NBLCNT
      DO ISHELA = IBLCKS(1,IBL),IBLCKS(2,IBL)
         IADR =   KSTRT(ISHELA)+1
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
            CENX = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CENT(ISHELA,1,1)
            CENY = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CENT(ISHELA,2,1)
            CENZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CENT(ISHELA,3,1)
            DO I=1,NVCLEN
               PA(1,i) = COOR(1,i)-CENX
               PA(2,i) = COOR(2,i)-CENY
               PA(3,i) = COOR(3,i)-CENZ
               PA2(i) = PA(1,i)**2 + PA(2,i)**2 + PA(3,i)**2
            END DO
            IF (NDER.GT.1 .OR. NHKTA.GT.2 .OR. DOLND) THEN
               CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
            END IF
            IF (NDER.EQ.0) THEN
               CALL BLGETGAO(NVCLEN,GAO(1,IADR,1),CSP(ISPADR(NHKTA)),
     &                        PA,PA2,DFTHRI)
            ELSE IF (NDER.GT.0) THEN
               CALL BLGETGA1(NVCLEN,NBAST,GAO(1,IADR,1),
     &              CSP(ISPADR(NHKTA)),PA,PA2,DFTHRI)
               IF (NDER.GT.1) THEN
                  CALL BLGETGA2(NVCLEN,
     &                 GAO(1,IADR,5),GAO(1,IADR,6),GAO(1,IADR,7),
     &                 GAO(1,IADR,8),GAO(1,IADR,9),GAO(1,IADR,10),
     &                 CSP(ISPADR(NHKTA)),PA,PA2,DFTHRI)
               END IF
            END IF
            IF (DOLND) THEN
               CALL BLGETGB1(NVCLEN,NBAST,GAO(1,IADR,1),PA,
     &              GAO(1,IADR,NSOB),GAO(1,IADR,NSOB+1),
     &              GAO(1,IADR,NSOB+2))
               IF (DOGGA) THEN
                  DO I = 1, 3
                     CALL BLGETGB1(NVCLEN,NBAST,GAO(1,IADR,1+I),PA,
     &                             GAO(1,IADR,NSOB+2+I),
     &                             GAO(1,IADR,NSOB+5+I),
     &                             GAO(1,IADR,NSOB+8+I))
                  END DO
               END IF
            END IF
            IADR = IADR + KHKTA
         END IF
         END DO
      END DO
      END DO
C     if (nder.gt.1) then
C        write (2,*) ' xx integrals from BLDFTAOS '
C        call output(gao(1,1,5),1,nvclen,1,nbast,nvclen,nbast,1,2)
C        write (2,*) ' xy integrals from BLDFTAOS '
C        call output(gao(1,1,6),1,nvclen,1,nbast,nvclen,nbast,1,2)
C        write (2,*) ' xz integrals from BLDFTAOS '
C        call output(gao(1,1,7),1,nvclen,1,nbast,nvclen,nbast,1,2)
C        write (2,*) ' yy integrals from BLDFTAOS '
C        call output(gao(1,1,8),1,nvclen,1,nbast,nvclen,nbast,1,2)
C        write (2,*) ' yz integrals from BLDFTAOS '
C        call output(gao(1,1,9),1,nvclen,1,nbast,nvclen,nbast,1,2)
C        write (2,*) ' zz integrals from BLDFTAOS '
C        call output(gao(1,1,10),1,nvclen,1,nbast,nvclen,nbast,1,2)
C     end if
      RETURN
      END
C  /* Deck blgetgao */
      SUBROUTINE BLGETGAO(NVCLEN,GAO,CSP,PA,PA2,DFTHRI)
C
C     T. Helgaker Sep 99, P. Salek 03
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <dftinf.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
      DIMENSION GAO(NVCLEN,KCKTA)
      DIMENSION PA(3,NVCLEN), PA2(NVCLEN)
      DIMENSION CSP(KHKTA,KCKTA)
c
      DIMENSION GA(NVCLEN),CINT(NVCLEN)
C
#include <lmns.h>
#include <onecom.h>
#include <primit.h>
C
      CALL DZERO(GA,NVCLEN)
C
C     loop over primitives
C
      DO I = JSTA + 1, JSTA + NUCA 
         DO K = 1, NVCLEN
            GA(K) = GA(K) + PRICCF(I,NUMCFA)*DEXP(-PRIEXP(I)*PA2(K))
         END DO
      END DO
C     
C     contracted orbitals
C
      IF (NHKTA .EQ. 1) THEN
         DO K = 1, NVCLEN
            GAO(K,1) = GA(K)
         END DO
      ELSE IF (NHKTA .EQ. 2) THEN
         DO K = 1, NVCLEN
            GAO(K,1) = PA(1,K)*GA(K)
            GAO(K,2) = PA(2,K)*GA(K)
            GAO(K,3) = PA(3,K)*GA(K)
         END DO
      ELSE IF (NHKTA .EQ. 3) THEN
         IF (SPHRA) THEN
            DO K = 1, NVCLEN
               GAX  = PA(1,K)*GA(K)
               GAY  = PA(2,K)*GA(K)
               GAZ  = PA(3,K)*GA(K)
               GAXX = PA(1,K)*GAX
               GAYY = PA(2,K)*GAY
               GAO(K,1) = CSP(1,2)*PA(2,K)*GAX
               GAO(K,2) = CSP(2,5)*PA(2,K)*GAZ
               GAO(K,3) = CSP(3,1)*GAXX + CSP(3,4)*GAYY
     &                  + CSP(3,6)*PA(3,K)*GAZ
               GAO(K,4) = CSP(4,3)*PA(1,K)*GAZ
               GAO(K,5) = CSP(5,1)*GAXX + CSP(5,4)*GAYY
            END DO
         ELSE
            DO K = 1, NVCLEN
               GAX  = PA(1,K)*GA(K)
               GAY  = PA(2,K)*GA(K)
               GAZ  = PA(3,K)*GA(K)
               GAO(K,1) = PA(1,K)*GAX
               GAO(K,2) = PA(2,K)*GAX
               GAO(K,3) = PA(3,K)*GAX
               GAO(K,4) = PA(2,K)*GAY
               GAO(K,5) = PA(3,K)*GAY 
               GAO(K,6) = PA(3,K)*GAZ
            END DO
         END IF
      ELSE
         IF (SPHRA) THEN
            DO I = 1, KHKTA
               DO K = 1, NVCLEN
                  GAO(K,I) = D0
               END DO
            END DO
            DO J = 1, KCKTA
               DO K = 1, NVCLEN
                  CINT(K) = (PA(1,K)**LVALUA(J))*(PA(2,K)**MVALUA(J))
     &                 *(PA(3,K)**NVALUA(J))*GA(K)
               END DO
c              do a dgemm here?
               DO I = 1, KHKTA
                  DO K = 1, NVCLEN
                     GAO(K,I) = GAO(K,I) + CSP(I,J)*CINT(K)
                  END DO
               END DO
            END DO
         ELSE
            DO I = 1, KCKTA
               DO K = 1, NVCLEN
                  GAO(K,I) = (PA(1,K)**LVALUA(I))*(PA(2,K)**MVALUA(I))
     &                 *(PA(3,K)**NVALUA(I))*GA(K)
               END DO
            END DO
         END IF
      END IF
      RETURN
      END
C  /* Deck blgetga1 */
      SUBROUTINE BLGETGA1(NVCLEN,NBAST,GAO,CSP,PA,PA2,DFTHRI)
C
C     T. Helgaker Sep 99, P. Salek 03
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <maxorb.h>
#include <maxaqn.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0) 
C
c     GAO(1..NVCLEN, 1..KCKTA,1..4) is set
      DIMENSION GAO(NVCLEN,NBAST,4)
      PARAMETER (I0=1,IX=2,IY=3,IZ=4)
      DIMENSION PA(3,NVCLEN), PA2(NVCLEN)
c      
      DIMENSION CAO(NVCLEN,MXAQN), CAOX(NVCLEN,MXAQN),
     &         CAOY(NVCLEN,MXAQN), CAOZ(NVCLEN,MXAQN)
      DIMENSION CSP(KHKTA,KCKTA)
      DIMENSION GA(NVCLEN), GU(NVCLEN),
     &          GAX(NVCLEN),GAY(NVCLEN), GAZ(NVCLEN), P0(NVCLEN),
     &          FX(NVCLEN), FY(NVCLEN), FZ(NVCLEN),
     &          SPH0(NVCLEN), SPHX(NVCLEN), SPHY(NVCLEN), SPHZ(NVCLEN)
C
#include <lmns.h>
#include <onecom.h>
#include <primit.h>
C
#include <sdpre.h>
C
C     loop over primitives
C
      DO K = 1, NVCLEN
         GA(K) = D0
         GU(K) = D0
      END DO
      DO I = JSTA + 1, JSTA + NUCA
         DO K = 1, NVCLEN
         FAC   = PRICCF(I,NUMCFA)*DEXP(-PRIEXP(I)*PA2(K))
         GA(K) = GA(K) + FAC
         GU(K) = GU(K) - D2*PRIEXP(I)*FAC
         END DO
      END DO 
C
C        s orbitals
C
      IF (NHKTA .EQ. 1) THEN
         DO K = 1, NVCLEN
         GAO(K,1,I0) = GA(K)
         GAO(K,1,IX) = PA(1,K)*GU(K)
         GAO(K,1,IY) = PA(2,K)*GU(K)
         GAO(K,1,IZ) = PA(3,K)*GU(K)
         END DO
C        
C        p orbitals
C        
      ELSE IF (NHKTA .EQ. 2) THEN
         DO K = 1, NVCLEN
         TGX = PA(1,K)*GU(K)
         TGY = PA(2,K)*GU(K)
         TGZ = PA(3,K)*GU(K)
         GAO(K,1,I0)  = PA(1,K)*GA(K)
         GAO(K,2,I0)  = PA(2,K)*GA(K)
         GAO(K,3,I0)  = PA(3,K)*GA(K)
         GAO(K,1,IX) = PA(1,K)*TGX + GA(K)
         GAO(K,2,IX) = PA(2,K)*TGX
         GAO(K,3,IX) = PA(3,K)*TGX
         GAO(K,1,IY) = PA(1,K)*TGY
         GAO(K,2,IY) = PA(2,K)*TGY + GA(K)
         GAO(K,3,IY) = PA(3,K)*TGY
         GAO(K,1,IZ) = PA(1,K)*TGZ
         GAO(K,2,IZ) = PA(2,K)*TGZ
         GAO(K,3,IZ) = PA(3,K)*TGZ + GA(K)
         END DO
C
C        d and higher orbitals
C
      ELSE
         DO K = 1, NVCLEN
            FX(K) = PA(1,K)*GU(K)
            FY(K) = PA(2,K)*GU(K)
            FZ(K) = PA(3,K)*GU(K)
         END DO
         DO ICOMPA = 1,KCKTA
            L = LVALUA(ICOMPA)
            M = MVALUA(ICOMPA)
            N = NVALUA(ICOMPA)
            DO K = 1, NVCLEN
               P0(K)  = (PA(1,K)**L)*(PA(2,K)**M)*(PA(3,K)**N)
               GAX(K) = FX(K)*P0(K)
               GAY(K) = FY(K)*P0(K)
               GAZ(K) = FZ(K)*P0(K)
               IF(L.GT.0) GAX(K) = GAX(K)
     &              + SDPRE(L)*(PA(1,K)**(L-1))*(PA(2,K)**M)*
     &              (PA(3,K)**N)*GA(K)
               IF(M.GT.0) GAY(K) = GAY(K)
     &              + SDPRE(M)*(PA(1,K)**L)*(PA(2,K)**(M-1))*
     &              (PA(3,K)**N)*GA(K)
               IF(N.GT.0) GAZ(K) = GAZ(K)
     &              + SDPRE(N)*(PA(1,K)**L)*(PA(2,K)**M)*
     &              (PA(3,K)**(N-1))*GA(K)
            END DO
            IF (SPHRA) THEN
               DO K = 1, NVCLEN
                  CAO (K, ICOMPA) = GA(K)*P0(K)
                  CAOX(K, ICOMPA) = GAX(K)
                  CAOY(K, ICOMPA) = GAY(K)
                  CAOZ(K, ICOMPA) = GAZ(K)
               END DO
            ELSE
               DO K = 1, NVCLEN
                  GAO(K,ICOMPA,I0) = GA(K)*P0(K)
                  GAO(K,ICOMPA,IX) = GAX(K)
                  GAO(K,ICOMPA,IY) = GAY(K)
                  GAO(K,ICOMPA,IZ) = GAZ(K)
               END DO
            END IF
         END DO
         IF(SPHRA) THEN
            IF (NHKTA.EQ.3) THEN
               DO K = 1, NVCLEN
                  GAO(K,1,I0) = CSP(1,2)*CAO(K,2)
                  GAO(K,2,I0) = CSP(2,5)*CAO(K,5)
                  GAO(K,3,I0) = CSP(3,1)*CAO(K,1) + CSP(3,4)*CAO(K,4) 
     &                        + CSP(3,6)*CAO(K,6)
                  GAO(K,4,I0) = CSP(4,3)*CAO(K,3)
                  GAO(K,5,I0) = CSP(5,1)*CAO(K,1) + CSP(5,4)*CAO(K,4)
                  GAO(K,1,IX) = CSP(1,2)*CAOX(K,2)
                  GAO(K,2,IX) = CSP(2,5)*CAOX(K,5)
                  GAO(K,3,IX) = CSP(3,1)*CAOX(K,1) + CSP(3,4)*CAOX(K,4) 
     &                        + CSP(3,6)*CAOX(K,6)
                  GAO(K,4,IX) = CSP(4,3)*CAOX(K,3)
                  GAO(K,5,IX) = CSP(5,1)*CAOX(K,1) + CSP(5,4)*CAOX(K,4)
                  GAO(K,1,IY) = CSP(1,2)*CAOY(K,2)
                  GAO(K,2,IY) = CSP(2,5)*CAOY(K,5)
                  GAO(K,3,IY) = CSP(3,1)*CAOY(K,1) + CSP(3,4)*CAOY(K,4)
     &                        + CSP(3,6)*CAOY(K,6)
                  GAO(K,4,IY) = CSP(4,3)*CAOY(K,3)
                  GAO(K,5,IY) = CSP(5,1)*CAOY(K,1) + CSP(5,4)*CAOY(K,4)
                  GAO(K,1,IZ) = CSP(1,2)*CAOZ(K,2)
                  GAO(K,2,IZ) = CSP(2,5)*CAOZ(K,5)
                  GAO(K,3,IZ) = CSP(3,1)*CAOZ(K,1) + CSP(3,4)*CAOZ(K,4) 
     &                        + CSP(3,6)*CAOZ(K,6)
                  GAO(K,4,IZ) = CSP(4,3)*CAOZ(K,3)
                  GAO(K,5,IZ) = CSP(5,1)*CAOZ(K,1) + CSP(5,4)*CAOZ(K,4)
               END DO
            ELSE
c              FIXME: take timings without use of the temporary arrays
               DO I = 1, KHKTA
                  DO K = 1, NVCLEN
                     SPH0(K) = D0 
                     SPHX(K) = D0 
                     SPHY(K) = D0 
                     SPHZ(K) = D0 
                  END DO
                  DO J = 1, KCKTA
                     SPHFAC = CSP(I,J)
                     IF (ABS(SPHFAC).GT.D0) THEN
                        DO K = 1, NVCLEN
                           SPH0(K) = SPH0(K) + SPHFAC*CAO (K,J)
                           SPHX(K) = SPHX(K) + SPHFAC*CAOX(K,J)
                           SPHY(K) = SPHY(K) + SPHFAC*CAOY(K,J)
                           SPHZ(K) = SPHZ(K) + SPHFAC*CAOZ(K,J)
                        END DO
                     END IF
                  END DO
                  DO K = 1, NVCLEN
                     GAO(K,I,I0) = SPH0(K)
                     GAO(K,I,IX) = SPHX(K)
                     GAO(K,I,IY) = SPHY(K)
                     GAO(K,I,IZ) = SPHZ(K)
                  END DO
               END DO
            END IF
         END IF
      END IF
      END
C  /* Deck blgetga2 */
      SUBROUTINE BLGETGA2(NVCLEN,
     &                    GAOXX,GAOXY,GAOXZ,GAOYY,GAOYZ,GAOZZ,CSP,
     &                    PA,PA2,DFTHRI)
C
C     T. Helgaker Sep 99
C     P. Salek    Oct 03 - vector version.
C
c     routine sets GAO__(1:NVCLEN, 1:KCKTA)
c
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0) 
C
      DIMENSION GAOXX(NVCLEN,KCKTA), GAOXY(NVCLEN,KCKTA),
     &          GAOXZ(NVCLEN,KCKTA), GAOYY(NVCLEN,KCKTA),
     &          GAOYZ(NVCLEN,KCKTA), GAOZZ(NVCLEN,KCKTA)
      DIMENSION CSP(KHKTA,KCKTA)
      DIMENSION PA(3,NVCLEN), PA2(NVCLEN)
      DIMENSION CAOXX(MXAQN), CAOXY(MXAQN),
     &          CAOXZ(MXAQN), CAOYY(MXAQN),
     &          CAOYZ(MXAQN), CAOZZ(MXAQN)
C
#include <lmns.h>
#include <onecom.h>
#include <primit.h>
C
#include <sdpre.h>
C
      DO IV = 1, NVCLEN
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
            TAPAX = TALPH*PA(1,IV)
            TAPAY = TALPH*PA(2,IV)
            TAPAZ = TALPH*PA(3,IV)
            GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-ALPHA*PA2(IV))
            IF (DABS(GA).GT.DFTHRI) THEN
               DO ICOMPA = 1, KCKTA
                  L = LVALUA(ICOMPA)
                  M = MVALUA(ICOMPA)
                  N = NVALUA(ICOMPA)
C                 
                  PXD = D0
                  PYD = D0
                  PZD = D0
                  IF (L.GT.1) PXD = SDPRE(L*(L-1))*(PA(1,IV)**(L-2))
                  IF (M.GT.1) PYD = SDPRE(M*(M-1))*(PA(2,IV)**(M-2))
                  IF (N.GT.1) PZD = SDPRE(N*(N-1))*(PA(3,IV)**(N-2))
                  PXM = D0
                  PYM = D0
                  PZM = D0
                  IF (L.GT.0) PXM = SDPRE(L)*(PA(1,IV)**(L-1))
                  IF (M.GT.0) PYM = SDPRE(M)*(PA(2,IV)**(M-1))
                  IF (N.GT.0) PZM = SDPRE(N)*(PA(3,IV)**(N-1))
                  PX0 = PA(1,IV)**L
                  PY0 = PA(2,IV)**M
                  PZ0 = PA(3,IV)**N
                  PXP = TAPAX*PX0
                  PYP = TAPAY*PY0
                  PZP = TAPAZ*PZ0
                  P000 = PX0*PY0*PZ0 
C                 
                  GAXX = (TAPAX**2 + TALPH*SDPRE(2*L+1))*P000
     &                 + PXD*PY0*PZ0
                  GAYY = (TAPAY**2 + TALPH*SDPRE(2*M+1))*P000
     &                 + PX0*PYD*PZ0
                  GAZZ = (TAPAZ**2 + TALPH*SDPRE(2*N+1))*P000
     &                 + PX0*PY0*PZD
                  GAXY = TAPAX*TAPAY*P000
     &                 + (PXP*PYM+PXM*PYP+PXM*PYM)*PZ0
                  GAXZ = TAPAX*TAPAZ*P000
     &                 + (PXP*PZM+PXM*PZP+PXM*PZM)*PY0
                  GAYZ = TAPAY*TAPAZ*P000
     &                 + (PYP*PZM+PYM*PZP+PYM*PZM)*PX0
C
                  IF (SPHRA) THEN
                     CAOXX(ICOMPA) = CAOXX(ICOMPA) + GAXX*GA 
                     CAOXY(ICOMPA) = CAOXY(ICOMPA) + GAXY*GA 
                     CAOXZ(ICOMPA) = CAOXZ(ICOMPA) + GAXZ*GA
                     CAOYY(ICOMPA) = CAOYY(ICOMPA) + GAYY*GA
                     CAOYZ(ICOMPA) = CAOYZ(ICOMPA) + GAYZ*GA
                     CAOZZ(ICOMPA) = CAOZZ(ICOMPA) + GAZZ*GA
                  ELSE
                     GAOXX(IV,ICOMPA) = GAOXX(IV,ICOMPA) + GAXX*GA 
                     GAOXY(IV,ICOMPA) = GAOXY(IV,ICOMPA) + GAXY*GA 
                     GAOXZ(IV,ICOMPA) = GAOXZ(IV,ICOMPA) + GAXZ*GA
                     GAOYY(IV,ICOMPA) = GAOYY(IV,ICOMPA) + GAYY*GA
                     GAOYZ(IV,ICOMPA) = GAOYZ(IV,ICOMPA) + GAYZ*GA
                     GAOZZ(IV,ICOMPA) = GAOZZ(IV,ICOMPA) + GAZZ*GA
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
               GAOXX(IV,I) = SPHXX
               GAOXY(IV,I) = SPHXY
               GAOXZ(IV,I) = SPHXZ
               GAOYY(IV,I) = SPHYY
               GAOYZ(IV,I) = SPHYZ
               GAOZZ(IV,I) = SPHZZ
            END DO
         END IF
      END DO
      RETURN
      END
C  /* Deck blgetgb1 */
      SUBROUTINE BLGETGB1(NVCLEN,NBAST,GAO,PA,GABX,GABY,GABZ)
C
C     T. Helgaker Oct 99
c     vectorized by Pawel Salek, Apr 03
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <dftinf.h>
C
      PARAMETER (DP5 = 0.5D0)
      DIMENSION GAO (NVCLEN,NBAST), GABX(NVCLEN,NBAST),
     &          GABY(NVCLEN,NBAST), GABZ(NVCLEN,NBAST),
     &          PA(3,NVCLEN)
      
C
#include <onecom.h>
#include <primit.h>
#include <orgcom.h>
C
      DO ICOMPA = 1,KHKTA
         DO K = 1, NVCLEN
            GA = DP5*GAO(K,ICOMPA)
            GABX(K,ICOMPA) = GA*PA(1,K)
            GABY(K,ICOMPA) = GA*PA(2,K)
            GABZ(K,ICOMPA) = GA*PA(3,K)
         END DO
      END DO
      RETURN
      END
C  /* Deck bldftsos */
      SUBROUTINE BLDFTSOS(NVCLEN,GAO,GSO,NBAST,NVEC,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      DIMENSION GAO(NVCLEN,NBAST,NVEC), GSO(NVCLEN,NBAST,NVEC)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DFTSOS',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(GSO,NVCLEN*NBAST*NVEC)
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
               DO J=1, NVEC
c                 
c                 ifc cannot vectorize it because "subscript too complex"
c                 
                  DO K=1,NVCLEN
                     GSO(K,I,J) = GSO(K,I,J) + FAC*GAO(K,INDA,J)
                  END DO
               END DO
c               print *,I,'th GSO has contribution from GAO=',INDA
            END IF
            END DO
         END DO
         ISTRA = ISTRA + NORBA
      END DO
      RETURN
      END
C
      SUBROUTINE SETUPSOS(GEODRV,DOLND,ITYPE,IOFF)
CC     sets DFTINF common block up so it is ready for subsequent
C     execution of GETSOS.
#include <implicit.h>
c     Pawel Salek, feb 2001
c     GEODRV - 0 for just orbital values, 1 for first order orbital
c     derivatives, 2 for laplacian.
c     ITYPE (out) - number of all orbitals sets/derivatives computed.
c     iOFF  (out) - array of offsets.
c     IOFF(1) contains offset of "london" derivatives.
      INTEGER GEODRV, ITYPE, IOFF(1)
c     do compute derivatives wrt magnetic field ("london" derivatives).
      LOGICAL DOLND
c     the computed orbitals are ordered as follows:
c     (orbital values=O, [dO/dx,dO/dy,dO/dz, [d0/dxx, ..]],
c      [dO/dBx, dO/dBy, dO/dBz])
#include <dftinf.h>
#include <inforb.h>
C
      NDER = GEODRV
      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2 
         IF (NDER.GT.0) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8 
         ELSE
C           assuming NSOB1 not referred.
            NSOB1 = 0
         END IF
      ELSE
c         assuming NSOB not referred.
c         assuming NSOB1 not referred?
         NSOB  = 0
         NSOB1 = 0
      END IF
c      NOSB1 defined in common block but never used? or is it NSOB1?
      KSO0  = (NSO0-1)*NBAST + 1
      KSO1  = (NSO1-1)*NBAST + 1
      KSO2  = (NSO2-1)*NBAST + 1
      KSOB  = (NSOB-1)*NBAST + 1
      KSOB1 = (NSOB1-1)*NBAST + 1
      ITYPE = NTYPSO
      IOFF(1) = NSOB
      END
C
      SUBROUTINE zeroorbs(TMP,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
#include <implicit.h>
#include <inforb.h>
      DIMENSION TMP(NVCLEN,NBAST),NBLOCKS(NSYM),IBLOCKS(2,LDAIB,NSYM)
      DO ISYM = 1, NSYM
         DO IBL=1,NBLOCKS(ISYM)
            DO IDX= IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
               DO J=1,NVCLEN
                  TMP(J,IDX) = 0D0
               END DO
            END DO
         END DO
      END DO
      END
