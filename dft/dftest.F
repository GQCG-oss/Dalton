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
C /* Deck dftest */
      SUBROUTINE DFTTST(WORK,LWORK,NBUF,LUQUAD,IPRINT)
C
C     T. Helgaker sep 99
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C 
      DIMENSION WORK(LWORK)
C
#include <inforb.h>
#include <nuclei.h>
C
      KX    = 1
      KY    = KX   + NBUF
      KZ    = KY   + NBUF
      KW    = KZ   + NBUF
      KGAO  = KW   + NBUF
      KGA1  = KW   + NBUF
      KDEN  = KGA1 + 3*NBAST
      KCNT  = KDEN + NBAST*NBAST
      KLST  = KCNT + NBAST
      LWRK  = LWORK - KLST + 1
      CALL DFTES1(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NBUF,
     &            WORK(KDEN),WORK(KGAO),WORK(KGA1),WORK(KCNT),
     &            WORK(KLST),LWRK,IPRINT,LUQUAD)
      RETURN
      END
C /* Deck dftest */
      SUBROUTINE DFTES1(CORX,CORY,CORZ,WEIGHT,NBUF,DMAT,GAO,GAO1,NCNT,
     &                  WORK,LWORK,IPRINT,LUQUAD)
C
C     T. Helgaker sep 99
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
C
#include <inforb.h>
#include <dftinf.h>
C
      DIMENSION CORX(NBUF), CORY(NBUF), CORZ(NBUF), WEIGHT(NBUF),
     &          DMAT(NBAST,NBAST), GAO(NBAST), GAO1(NBAST,3),
     &          NCNT(NBAST), WORK(LWORK)
C
      LOGICAL NODV, NODC, DOLND
C
      DOLND = .FALSE.
C
C     Calculate density matrix
C     
      NODC = .FALSE.
      NODV = .FALSE. 
      CALL GETDEN(DMAT,WORK,LWORK,NODC,NODV,IPRINT)
      ELCTRN = D0
C
      NDER = 0
C
C     Integrate
C
      REWIND LUQUAD
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         CALL REAQUA(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
         DO 300 IPNT = 1, NPOINT
            CALL DFTAOS(GAO,GAO1,NCNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &                  WORK,LWORK,NBAST,DOLND,IPRINT)
            DO 400 I = 1, NBAST
            DO 400 J = 1, NBAST
               ELCTRN = ELCTRN + WEIGHT(IPNT)*DMAT(I,J)*GAO(I)*GAO(J)
  400       CONTINUE 
  300    CONTINUE
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      ELSE
         GO TO 200
      END IF
  200 CONTINUE
C
      WRITE (LUPRI,'(2X,A,F20.12)') ' Number of electrons ',ELCTRN
C
      RETURN 
      END
C  /* Deck dftnum */
      SUBROUTINE DFTNUM(DODRC,DOVWN,DOLYP,DOBCK,RHO,RHO13,RHOGRD)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, 
     &           DP3 = D1/D3, ADD = 1.0D-06)
      LOGICAL DODRC,DOVWN,DOLYP,DOBCK
C
      ADD2 = D2*ADD
C
      R0 = RHO
      RP = RHO + ADD
      RM = RHO - ADD
      Z0 = RHOGRD
      ZP = RHOGRD + ADD 
      ZM = RHOGRD - ADD
      R013 = R0**DP3
      RP13 = RP**DP3
      RM13 = RM**DP3
C
      IF (DODRC) THEN
         CALL VDRC(G10,R013)
         CALL V1DRC(F10,F20,R0,R013)
C
         CALL EDRC(EP,RP,RP13)
         CALL VDRC(V1P,RP13)
         CALL EDRC(EM,RM,RM13)
         CALL VDRC(V1M,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' DR010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC20 test',F20,H20,F20-H20
      END IF
      IF (DOVWN) THEN
         CALL VVWN(G10,R0,R013)
         CALL V1VWN(F10,F20,R0,R013)
C
         CALL EVWN(EP,RP,RP13)
         CALL VVWN(V1P,RP,RP13)
         CALL EVWN(EM,RM,RM13)
         CALL VVWN(V1M,RM,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' VW010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN20 test',F20,H20,F20-H20
      END IF
C
      IF (DOBCK) THEN
         CALL GBCK(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1BCK(F10,F01,F20,F11,F02,R0,Z0)
C
         CALL EBCK(EP,RP,RP13,Z0)
         CALL GBCK(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL EBCK(EM,RM,RM13,Z0)
         CALL GBCK(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
C
         CALL EBCK(EP,R0,R013,ZP)
         CALL GBCK(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL EBCK(EM,R0,R013,ZM)
         CALL GBCK(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' BC010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BC001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BCK01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' BCK11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' BCK02 test',F02,H02,F02-H02
      END IF
      IF (DOLYP) THEN
         CALL GLYP(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1LYP(F10,F01,F20,F11,F02,R0,Z0)
C
         CALL ELYP(EP,RP,RP13,Z0)
         CALL GLYP(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL ELYP(EM,RM,RM13,Z0)
         CALL GLYP(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
C
         CALL ELYP(EP,R0,R013,ZP)
         CALL GLYP(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL ELYP(EM,R0,R013,ZM)
         CALL GLYP(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' LY010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LY001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LYP01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' LYP11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' LYP02 test',F02,H02,F02-H02
      END IF
      RETURN
      END
