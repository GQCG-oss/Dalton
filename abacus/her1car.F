C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
C  /* Deck oneprm */
      SUBROUTINE ONEPRM(STDER0,STDER1,STDER2,ADER,SINT0,DINT1,
     &                  RLMINT,RLMTAB,FCM,
     &                  WORK,LWORK,IPRINT,PROPTY,MAXDIF,NATOMC,TOLOG,
     &                  TOLS,SECDER,DIFDIP,FACINT,COORC,
     &                  GNUEXP,NCENTC,NCLONE,EXPDER)
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
C
      LOGICAL PROPTY, SECDER, DIFDIP, NCLONE
      DIMENSION WORK(LWORK)
      DIMENSION STDER0(KCKTAB,2), STDER1(KCKTAB,3,2),
     &          STDER2(KCKTAB,6,2),
     &          SINT0(KCKTAB), DINT1(KCKTAB,3,3),
     &          ADER(*), FACINT(*), COORC(3,*), NCENTC(*),
     &          RLMINT(*), RLMTAB(*), FCM(*), GNUEXP(*), EXPDER(*)
#include <onecom.h>
#include <expopt.h>
C
      JMAXD = 2
      IF (PROPTY .AND. .NOT.ONECEN .AND. .NOT.EXPGRA) JMAXD = 4
C
      JMAXA = NHKTA - 1
      JMAXB = NHKTB - 1
      IF (EXPGRA) THEN
         JMAXA = JMAXA + 2
         JMAXB = JMAXB + 2
      END IF
      JMAXT = JMAXA + JMAXB + JMAXD
      JMAXM = 0
C
      KAHGTF = 1
      KODC   = KAHGTF + (NATOMC + 1)*NAHGTF
      KLAST  = KODC   + 3*(JMAXA+1)*(JMAXB+1)*(JMAXT+1)*(JMAXD+1)
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ONEPRM',' ',KLAST,LWORK)
      CALL DZERO(WORK(KAHGTF),(NATOMC + 1)*NAHGTF)
      CALL ONEPR1(STDER0,STDER1,STDER2,ADER,SINT0,DINT1,
     &            RLMINT,RLMTAB,FCM,
     &            WORK(KAHGTF),WORK(KODC),JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &            IPRINT,PROPTY,MAXDIF,NATOMC,TOLOG,TOLS,SECDER,DIFDIP,
     &            FACINT,COORC,GNUEXP,NCENTC,NHKTAP,NCLONE,
     &            EXPDER,WORK(KLAST),LWRK)
      RETURN
      END
C  /* Deck onepr1 */
      SUBROUTINE ONEPR1(STDER0,STDER1,STDER2,ADER,SINT0,DINT1,
     &                  RLMINT,RLMTAB,FCM,
     &                  AHGTF,ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,IPRINT,
     &                  PROPTY,MAXDIF,NATOMC,TOLOG,TOLS,SECDER,DIFDIP,
     &                  FACINT,COORC,GNUEXP,NCENTC,NHKTAP,NCLONE,
     &                  EXPDER,WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <pi.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D3 = 3.0D0, THIRD = D1/D3,
     &           GNUTHR = 1.D-4)
C
      LOGICAL PROPTY, KINODC, SECDER, DIFDIP, NCLONE, CLASIC(MXCENT)
      DIMENSION WORK(LWORK), STDER0(KCKTAB,2), STDER1(KCKTAB,3,2),
     &          STDER2(KCKTAB,6,2), RLMINT(*), RLMTAB(*), FCM(*),
     &          SINT0(KCKTAB), DINT1(KCKTAB,3,3),
     &          ADER(*), AHGTF(*), NCENTC(*), FACINT(*), COORC(3,*),
     &          GNUEXP(*), EXPDER(KCKTAB,2,2,NUCA,NUCB)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
C
#include <cbisol.h>
#include <onecom.h>
#include <ader.h>
#include <clsfmm.h>
#include <primit.h>
#include <ecpinf.h>
#include <expopt.h>
C
      IF (IPRINT .GT. 4) CALL TITLER('Output from ONEPR1','*',103)
C
      DIFABX = CORAX - CORBX
      DIFABY = CORAY - CORBY
      DIFABZ = CORAZ - CORBZ
      DISTAB = DIFABX*DIFABX + DIFABY*DIFABY + DIFABZ*DIFABZ
C
      NATOT  = 0
      IF (PROPTY .AND..NOT.EXPGRA) THEN
         IF (ONECEN) THEN
            IA0000 =  1
            IA000X =  2
            IA000Y =  3
            IA000Z =  4
            IA00XX =  5
            IA00XY =  6
            IA00XZ =  7
            IA00YY =  8
            IA00YZ =  9
            IA00ZZ = 10
            NATOT  = 10*NATOMC*KCKTAB
            CALL DZERO(ADER,NATOT)
         ELSE
            CALL DZERO(STDER0,2*KCKTAB)
            CALL DZERO(STDER1,6*KCKTAB)
            IF (SECDER) CALL DZERO(STDER2,12*KCKTAB)
C
            IA0000 =  1
            IA0X00 =  2
            IA0Y00 =  3
            IA0Z00 =  4
            IAXX00 =  5
            IAXY00 =  6
            IAXZ00 =  7
            IAYY00 =  8
            IAYZ00 =  9
            IAZZ00 = 10
            IA000X = 11
            IA000Y = 12
            IA000Z = 13
            IA00XX = 14
            IA00XY = 15
            IA00XZ = 16
            IA00YY = 17
            IA00YZ = 18
            IA00ZZ = 19
            IA0X0X = 20
            IA0X0Y = 21
            IA0X0Z = 22
            IA0Y0X = 23
            IA0Y0Y = 24
            IA0Y0Z = 25
            IA0Z0X = 26
            IA0Z0Y = 27
            IA0Z0Z = 28
            NATOT  = 28*NATOMC*KCKTAB
            CALL DZERO(ADER,NATOT)
         END IF
      ELSE
         CALL DZERO(STDER0,2*KCKTAB)
         CALL DZERO(ADER,KCKTAB)
      END IF
      IF (PROPTY) THEN
         CALL DZERO(SINT0,  KCKTAB)
         CALL DZERO(DINT1,9*KCKTAB)
      END IF
      IF (EXPGRA) THEN
         CALL DZERO(EXPDER,4*KCKTAB*NUCA*NUCB)
      END IF
      IF (SOLVNT) THEN
         IF (MAXDIF .EQ. 0) THEN
            CALL DZERO(RLMINT,LMNTOT*KCKTAB)
         ELSE IF (MAXDIF .GE. 1) THEN
            CALL DZERO(RLMINT,7*LMNTOT*KCKTAB)
            IF (MAXDIF .GE. 2) CALL DZERO(RLMTAB,21*KCKTAB)
         END IF
      END IF
C
C     Non-classical contributions
C
      IF (NCLONE) THEN
         FAC = ERFCIV(THRCLS)
         DO IATOMC = 1, NATOMC
            CLASIC(IATOMC) = .TRUE.
         END DO
         DO IPRIMA = 1,NUCA
            JPRIMA = JSTA + IPRIMA
            CONTA = PRICCF(JPRIMA,NUMCFA)
            EXPA = PRIEXP(JPRIMA)
            DO IPRIMB = 1,NUCB
               JPRIMB = JSTB + IPRIMB
               CONTB = PRICCF(JPRIMB,NUMCFB)
               EXPB = PRIEXP(JPRIMB)
               EXPP = EXPA + EXPB
               EXPPI = D1/EXPP
               EXPAPI = EXPA*EXPPI
               EXPBPI = EXPB*EXPPI
               CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
               CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
               CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ
               EP     = FAC*SQRT(EXPPI)
               DO IATOMC = 1,NATOMC
                  RPQ = SQRT((COORC(1,IATOMC) - CORPX)**2 +
     &                       (COORC(2,IATOMC) - CORPY)**2 +
     &                       (COORC(3,IATOMC) - CORPZ)**2)
                  IF (RPQ .LE. EP) CLASIC(IATOMC) = .FALSE.
               END DO
            END DO
         END DO
      END IF 
C
C     ****************************************
C     ***** Loop over primitive orbitals *****
C     ****************************************
C
      DO 100 IPRIMA = 1,NUCA
         JPRIMA = JSTA + IPRIMA
         CONTA = PRICCF(JPRIMA,NUMCFA)
         EXPA = PRIEXP(JPRIMA)
         DO 200 IPRIMB = 1,NUCB
            JPRIMB = JSTB + IPRIMB
            CONTB = PRICCF(JPRIMB,NUMCFB)
            EXPB = PRIEXP(JPRIMB)
            EXPP = EXPA + EXPB
            EXPPI = D1/EXPP
C
C           Calculate and test square root factor
C
            EXPABQ = EXPA*EXPB*DISTAB*EXPPI
         IF (EXPABQ.GT.TOLOG) GO TO 200
            SAAB = CONTA*CONTB*EXP(-EXPABQ)
            ASAAB = ABS(SAAB)
         IF (ASAAB.LT.TOLS) GO TO 200
            SAAB13 = SIGN(ASAAB**THIRD,SAAB)
C
C           Calculate coordinates of product Gaussian P
C
            EXPAPI = EXPA*EXPPI
            EXPBPI = EXPB*EXPPI
            CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
            CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
            CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ
C
C           *********************************************
C           ***** Overlap Distribution Coefficients *****
C           *********************************************
C
C           Expansion coefficients for undifferentiated orbitals
C
            KINODC = .NOT.PROPTY
            IDUMMY = 0
            IF (.NOT.EXPGRA) THEN
               CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,PROPTY,
     &                  KINODC,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,.TRUE.,.FALSE.,
     &                  ORIGIN,IDUMMY)
            ELSE
              CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,.FALSE.,
     &                    .TRUE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                    WORK,LWORK,CORPX,CORPY,CORPZ,.TRUE.,.FALSE.,
     &                    ORIGIN,IDUMMY)
            END IF
C
C           **********************************************
C           ***** Calculation of Hermitian integrals *****
C           **********************************************
C
C           Overlap integral
C
            SHGTF = SQRT(PI*EXPPI)
C
C           Nuclear attraction integrals
C
            IADR = 1
            DO 300 IATOMC = 1,NATOMC
               ICENTC = NCENTC(IATOMC)
               IF(.NOT.(PROPTY.AND.ONECEN.AND.ICENTA.EQ.ICENTC).OR.
     &            EXPGRA) THEN
                  FACTOR = FACINT(IATOMC)
                  IF (NCLONE .AND. CLASIC(IATOMC)) FACTOR = D0
                  DIFCPX = COORC(1,IATOMC) - CORPX
                  DIFCPY = COORC(2,IATOMC) - CORPY
                  DIFCPZ = COORC(3,IATOMC) - CORPZ
                  IF(ABS(GNUEXP(IATOMC)).GT.GNUTHR) THEN
                    EXPFAC = GNUEXP(IATOMC)
                    EXPFAC = EXPFAC/(EXPFAC + EXPP)
                    EXPPGN = EXPP*EXPFAC
                    FACTOR = FACTOR*(EXPFAC**1.5D0)
                    CALL HERNAI(AHGTF,JMAX,EXPPGN,DIFCPX,DIFCPY,
     *                DIFCPZ,FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,
     *                IPRINT)
                  ELSE
                    CALL HERNAI(AHGTF,JMAX,EXPP,DIFCPX,DIFCPY,
     *                DIFCPZ,FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,
     *                IPRINT)
                  ENDIF
               END IF
               IADR = IADR + NAHGTF
  300       CONTINUE
C
C           **********************************************
C           ***** Calculation of Cartesian integrals *****
C           **********************************************
C
            IF (EXPGRA) THEN
               IF (ECP) CALL QUIT(
     &            'ECP not implemented yet in ONEPR1 for "EXPGRA" !')
               CALL CINTAB(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,EXPA,EXPB,
     &                     EXPDER(1,1,1,IPRIMA,IPRIMB),SHGTF,AHGTF,
     &                     NATOMC)
            ELSE IF (PROPTY) THEN
               IF (ECP) CALL QUIT(
     &            'ECP not implemented yet in ONEPR1 for "PROPTY" !')
               IF (ONECEN) THEN
                  CALL CINT1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,ADER,
     &                       AHGTF,NCENTC,NATOMC,SECDER)
               ELSE
                  CALL CINT2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,STDER0,
     &                       STDER1,STDER2,ADER,SHGTF,AHGTF,NATOMC,
     &                       SECDER)
               END IF
               CALL DIPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SINT0,
     &                     DINT1,SHGTF,CORPX,CORPY,CORPZ,DIFDIP)
            ELSE
               CALL CINT0(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,STDER0,
     &                    SHGTF,ADER,AHGTF,NATOMC,IPRINT)
               IF (ECP) THEN
                  CALL ECPINT(CORAX,CORAY,CORAZ,CORBX,CORBY,CORBZ,
     &                 EXPA,EXPB,JMAXA,JMAXB,CONTA,CONTB,ADER,
     &                 WORK,LWORK,IPRINT)
               END IF
            END IF
            IF (SOLVNT) THEN
               CALL SOLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                     RLMINT,RLMTAB,CORPX,CORPY,CORPZ,MAXDIF,
     &                     EXPPI,FCM,WORK,LWORK,IPRINT)
            END IF
  200    CONTINUE
  100 CONTINUE
C
C     Print
C
      IF (IPRINT .GE. 5) THEN
         IF (.NOT. ONECEN) THEN
            WRITE (LUPRI,'(/A/)') ' Undifferentiated overlap integrals '
            WRITE (LUPRI,'(1P,6D12.4)') (STDER0(I,1),I=1,KCKTAB)
         END IF
      END IF
      RETURN
      END
C  /* Deck hernai */
      SUBROUTINE HERNAI(AHGTF,JMAX,PVAL,CPX,CPY,CPZ,FACINT,IADR,
     &                  ISTEPU,ISTEPV,NRTUV,IPRINT)
#include <implicit.h>
C
C     This subroutine calculates the R integrals as defined by
C     McMurchie and Davidson in J. Comp. Phys. 26 (1978) 218.
C     The recursion formulas (4.6) - (4.8) are used.
C
C     Number of elements generated: JMAXP*(JMAXP+1)*(JMAXP+2)/6
C
C     The array AHGTF(J) must be dimensioned as
C
C     J = (NUMNUC + 1)*((MAX(JMAX) + 1)**3)
C
C     where NUCNUM Is the total number of nuclei.
C
C     TUH spring 1984
C     Modified TUH 880613 to avoid negative addressing in recursion.
C
#include <priunit.h>
#include <pi.h>
#include <maxaqn.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.D0, D2 = 2.D0, D2PI = D2*PI)
      DIMENSION AHGTF(*)
#include <gamcom.h>
#include <ibtfun.h>
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(/A/)') ' <<<<< Output from Hernai >>>>> '
         WRITE (LUPRI,'(A,I10)')   ' JMAX   ', JMAX
         WRITE (LUPRI,'(A,I10)')   ' IADR   ', IADR
         WRITE (LUPRI,'(A,I10)')   ' ISTEPU ', ISTEPU
         WRITE (LUPRI,'(A,I10)')   ' ISTEPV ', ISTEPV
         WRITE (LUPRI,'(A,I10)')   ' NRTUV  ', NRTUV
         WRITE (LUPRI,'(A,F12.6)') ' PVAL   ', PVAL
         WRITE (LUPRI,'(A,F12.6)') ' FACINT ', FACINT
         WRITE (LUPRI,'(A,3F12.6)') ' CP ', CPX, CPY, CPZ
      END IF
C
C     *************************************
C     ***** Incomplete Gamma Function *****
C     *************************************
C
      WVAL = PVAL*(CPX*CPX + CPY*CPY + CPZ*CPZ)
      JMAX0 = JMAX
      CALL GAMFUN
      IF (IPRINT .GT. 20) THEN
         CALL HEADER ('FJW after GAMFUN',-1)
         WRITE (LUPRI,'(6F12.6)') (FJW(I), I = 1, JMAX)
      END IF
C
C     **********************************
C     ***** Special Case: JMAX = 0 *****
C     **********************************
C
      IF (JMAX .EQ. 0) THEN
         AHGTF(IADR) = FACINT*D2PI*FJW(0)/PVAL
         RETURN
      END IF
C
C     **********************************
C     ***** General Case: JMAX > 0 *****
C     **********************************
C
      IF (IBTAND(JMAX,1) .EQ. 0) THEN
         ISTRTJ = IADR
         ISTEPJ = NRTUV
      ELSE
         ISTRTJ = IADR + NRTUV
         ISTEPJ = - NRTUV
      END IF
      D2PVAL = PVAL + PVAL
      FACTOR = FACINT*D2PI/PVAL
      DO 100 JVAL = 0, JMAX
         FJW(JVAL) = FACTOR*FJW(JVAL)
         FACTOR    = - D2PVAL*FACTOR
  100 CONTINUE
      IF (IPRINT .GT. 20) THEN
         CALL HEADER ('FJW after multiplication by FACTOR',-1)
         WRITE (LUPRI,'(6F12.6)') (FJW(I), I = 1, JMAX)
      END IF
C
C     ***** JVAL = 0 *****
C
      AHGTF(ISTRTJ) = FJW(JMAX)
C
C     ***** JVAL = 1 *****
C
      ISTRTJ                 =   ISTRTJ + ISTEPJ
      FJWMAX                 =   FJW(JMAX)
      AHGTF(ISTRTJ)          =   FJW(JMAX - 1)
      AHGTF(ISTRTJ +      1) = - CPX*FJWMAX
      AHGTF(ISTRTJ + ISTEPU) = - CPY*FJWMAX
      AHGTF(ISTRTJ + ISTEPV) = - CPZ*FJWMAX
C
      IF (JMAX .GT. 1) THEN
C
C        ***** JVAL > 1 *****
C
         ISTPTU = 1 - ISTEPU
         ISTEPJ =   - ISTEPJ
         DO 200 JVAL = 2, JMAX
            ISTRTJ = ISTRTJ + ISTEPJ
            ISTRTV = ISTRTJ
            ISTRTU = ISTRTV
            IPREV1 = ISTRTU - ISTEPJ - 1
            IPREV2 = IPREV1 - 1
C
C           RJ(0,0,0)
C
            AHGTF(ISTRTU) = FJW(JMAX - JVAL)
C
C           RJ(1,0,0)
C
            AHGTF(ISTRTU + 1) = - CPX*AHGTF(IPREV1 + 1)
C
C           RJ(t,0,0) for  t > 1
C
            TMIN1 = D1
            DO 300 IT = 2, JVAL
               AHGTF(ISTRTU + IT) = - CPX*AHGTF(IPREV1 + IT)
     *                            + TMIN1*AHGTF(IPREV2 + IT)
               TMIN1 = TMIN1 + D1
  300       CONTINUE
C
C           RJ(t,1,0)
C
            ISTRTU = ISTRTU + ISTEPU
            IPREV1 = IPREV1 + 1
            DO 400 IT = 0, JVAL - 1
               AHGTF(ISTRTU + IT) = - CPY*AHGTF(IPREV1 + IT)
  400       CONTINUE
C
C           RJ(t,u,0) for  u > 1
C
            UMIN1 = D1
            DO 500 IU = 2,JVAL
               ISTRTU = ISTRTU + ISTEPU
               IPREV1 = IPREV1 + ISTEPU
               IPREV2 = IPREV1 - ISTEPU
               DO 510 IT = 0, JVAL - IU
                  AHGTF(ISTRTU + IT) = - CPY*AHGTF(IPREV1 + IT)
     *                               + UMIN1*AHGTF(IPREV2 + IT)
  510          CONTINUE
               UMIN1 = UMIN1 + D1
  500       CONTINUE
C
C           RJ(t,u,1)
C
            ISTRTV = ISTRTV + ISTEPV
            ISTRTU = ISTRTV
            IPREV1 = ISTRTU - ISTEPJ - ISTEPV
            IUMAX = JVAL - 1
            DO 600 IU = 0, IUMAX
               DO 610 IT = 0, IUMAX - IU
                  AHGTF(ISTRTU + IT) = - CPZ*AHGTF(IPREV1 + IT)
  610          CONTINUE
               ISTRTU = ISTRTU + ISTEPU
               IPREV1 = IPREV1 + ISTEPU
  600       CONTINUE
C
C           RJ(t,u,v) for v > 1
C
            VMIN1 = D1
            DO 700 IV = 2,JVAL
               ISTRTV = ISTRTV + ISTEPV
               ISTRTU = ISTRTV
               IPREV1 = ISTRTU - ISTEPJ - ISTEPV
               IPREV2 = IPREV1 - ISTEPV
               IUMAX = JVAL - IV
               DO 710 IU = 0, IUMAX
                  DO 720 IT = 0, IUMAX - IU
                     AHGTF(ISTRTU + IT) = - CPZ*AHGTF(IPREV1 + IT)
     *                                  + VMIN1*AHGTF(IPREV2 + IT)
  720             CONTINUE
                  ISTRTU = ISTRTU + ISTEPU
                  IPREV1 = IPREV1 + ISTEPU
                  IPREV2 = IPREV2 + ISTEPU
  710          CONTINUE
               VMIN1 = VMIN1 + D1
  700       CONTINUE
            ISTEPJ = - ISTEPJ
  200    CONTINUE
      END IF
      IF (IPRINT .GT. 20) THEN
         CALL HEADER ('Hermite integrals in HERNAI',-1)
         WRITE (LUPRI,'(6F12.6)') (AHGTF(IADR + I - 1), I = 1, NRTUV)
      END IF
      RETURN
      END
C  /* Deck cint0 */
      SUBROUTINE CINT0(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,STDER0,SHGTF,
     &                 ADER,AHGTF,NATOMC,IPRINT)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0 D00, DP5 = 0.5 D00)
      DIMENSION STDER0(KCKTAB,2), AHGTF(*), ADER(*)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <onecom.h>
#include <lmns.h>
      SFAC = SHGTF**3
      TFAC = -DP5*SFAC
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C    **********************************************************
C    ***** CALCULATE OVERLAP AND KINETIC ENERGY INTEGRALS *****
C    **********************************************************
C
         X0 = ODC(LVALA,LVALB,0,0,0,1)
         Y0 = ODC(MVALA,MVALB,0,0,0,2)
         Z0 = ODC(NVALA,NVALB,0,0,0,3)
         X2 = ODC(LVALA,LVALB,0,2,0,1)
         Y2 = ODC(MVALA,MVALB,0,2,0,2)
         Z2 = ODC(NVALA,NVALB,0,2,0,3)
         INT = INT + 1
         STDER0(INT,1) = STDER0(INT,1)+SFAC*X0*Y0*Z0
         STDER0(INT,2) = STDER0(INT,2)+TFAC*(X2*Y0*Z0+X0*Y2*Z0+X0*Y0*Z2)
C
C     **************************************************
C     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C     **************************************************
C
         IADRAV = 1
         AINT = D0
         DO 200 IV = 0, NVALA + NVALB
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MVALA + MVALB
               EE = ODC(MVALA,MVALB,IU,0,0,2)*EV
               DO 400 IT = 0, LVALA + LVALB
                  EEE = ODC(LVALA,LVALB,IT,0,0,1)*EE
                  IADR00 = IADRAU + IT
                  IADD = - NAHGTF
                  DO 500 IATOM = 1,NATOMC
                     IADD = IADD + NAHGTF
                     AINT = AINT + EEE*AHGTF(IADR00 + IADD)
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         ADER(INT) = ADER(INT) + AINT
C
C        PRINT SECTION
C
         IF (IPRINT .GE. 10) THEN
            CALL HEADER ('Output from CINT0',-1)
            WRITE (LUPRI,'(//,A,2I5)') ' ICOMPA, ICOMPB ', ICOMPA,ICOMPB
            WRITE (LUPRI,'(/,A,2F12.6)') ' SFAC, TFAC ', SFAC, TFAC
            WRITE (LUPRI,'(/,A,3F12.6)') ' SINT, TINT, AINT ',
     *        SFAC*X0*Y0*Z0, TFAC*(X2*Y0*Z0 + X0*Y2*Z0 + X0*Y0*Z2), AINT
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck cint1 */
      SUBROUTINE CINT1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,ADER,AHGTF,
     &                 NCENTC,NATOMC,SECDER)
C
C     TUH
C
#include <implicit.h>
#include <maxaqn.h>
#include <mxcent.h>
      LOGICAL SECDER
      DIMENSION ADER(KCKTAB,NATOMC,*), AHGTF(*), NCENTC(*)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <onecom.h>
#include <ader.h>
#include <lmns.h>
C
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         INT = INT + 1
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         MAXT  = LVALA + LVALB
         MAXU  = MVALA + MVALB
         MAXV  = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EE = ODC(MVALA,MVALB,IU,0,0,2)*EV
               DO 400 IT = 0,MAXT
                  EEE = ODC(LVALA,LVALB,IT,0,0,1)*EE
                  IADR00 = IADRAU + IT
                  IADR0T = IADR00 + 1
                  IADR0U = IADR00 + ISTEPU
                  IADR0V = IADR00 + ISTEPV
                  IF (SECDER) THEN
                     IADRTT = IADR0T + 1
                     IADRTU = IADR0T + ISTEPU
                     IADRTV = IADR0T + ISTEPV
                     IADRUU = IADR0U + ISTEPU
                     IADRUV = IADR0U + ISTEPV
                     IADRVV = IADR0V + ISTEPV
                  END IF
                  IADD = - NAHGTF
C
C                 ***** LOOP OVER NUCLEI *****
C
                  DO 500 I = 1, NATOMC
                     IADD = IADD + NAHGTF
                     ICENTC = NCENTC(I)
                     IF (ICENTC .NE. ICENTA) THEN
                        AH00 = AHGTF(IADR00 + IADD)
                        AH0T = AHGTF(IADR0T + IADD)
                        AH0U = AHGTF(IADR0U + IADD)
                        AH0V = AHGTF(IADR0V + IADD)
C
C                       Undifferentiated integral:
C
                        ADER(INT,I,IA0000) = ADER(INT,I,IA0000)+EEE*AH00
C
C                       C differentiated integrals:
C
                        ADER(INT,I,IA000X) = ADER(INT,I,IA000X)-EEE*AH0T
                        ADER(INT,I,IA000Y) = ADER(INT,I,IA000Y)-EEE*AH0U
                        ADER(INT,I,IA000Z) = ADER(INT,I,IA000Z)-EEE*AH0V
C
C                       Second derivatives:
C
                        IF (SECDER) THEN
                          AHTT = AHGTF(IADRTT + IADD)
                          AHTU = AHGTF(IADRTU + IADD)
                          AHTV = AHGTF(IADRTV + IADD)
                          AHUU = AHGTF(IADRUU + IADD)
                          AHUV = AHGTF(IADRUV + IADD)
                          AHVV = AHGTF(IADRVV + IADD)
C
C                          C-C differentiated integrals:
C
                          ADER(INT,I,IA00XX)=ADER(INT,I,IA00XX)+EEE*AHTT
                          ADER(INT,I,IA00XY)=ADER(INT,I,IA00XY)+EEE*AHTU
                          ADER(INT,I,IA00XZ)=ADER(INT,I,IA00XZ)+EEE*AHTV
                          ADER(INT,I,IA00YY)=ADER(INT,I,IA00YY)+EEE*AHUU
                          ADER(INT,I,IA00YZ)=ADER(INT,I,IA00YZ)+EEE*AHUV
                          ADER(INT,I,IA00ZZ)=ADER(INT,I,IA00ZZ)+EEE*AHVV
                        END IF
                     END IF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck cint2 */
      SUBROUTINE CINT2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,STDER0,STDER1,
     &                 STDER2,ADER,SHGTF,AHGTF,NATOMC,SECDER)
C
C     TUH
C
#include <implicit.h>
#include <maxaqn.h>
#include <mxcent.h>
      PARAMETER (DP5 = 0.5 D00)
      LOGICAL SECDER
      DIMENSION STDER0(KCKTAB,2),
     &          STDER1(KCKTAB,3,2), STDER2(KCKTAB,6,2),
     &          ADER(KCKTAB,NATOMC,*), AHGTF(*)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <onecom.h>
#include <ader.h>
#include <lmns.h>
      MAXADD = 1
      IF (SECDER) MAXADD = 2
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C    **********************************************************
C    ***** CALCULATE OVERLAP AND KINETIC ENERGY INTEGRALS *****
C    **********************************************************
C
         DERX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         DERY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         DERZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DERX1 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DERY1 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DERZ1 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DERX2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DERY2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DERZ2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         DERX3 = SHGTF*ODC(LVALA,LVALB,0,3,0,1)
         DERY3 = SHGTF*ODC(MVALA,MVALB,0,3,0,2)
         DERZ3 = SHGTF*ODC(NVALA,NVALB,0,3,0,3)
         IF (SECDER) THEN
            DERX4 = SHGTF*ODC(LVALA,LVALB,0,4,0,1)
            DERY4 = SHGTF*ODC(MVALA,MVALB,0,4,0,2)
            DERZ4 = SHGTF*ODC(NVALA,NVALB,0,4,0,3)
         END IF
         INT = INT + 1
         STDER0(INT,1)   = STDER0(INT,1)   + DERX0*DERY0*DERZ0
         STDER1(INT,1,1) = STDER1(INT,1,1) + DERX1*DERY0*DERZ0
         STDER1(INT,2,1) = STDER1(INT,2,1) + DERX0*DERY1*DERZ0
         STDER1(INT,3,1) = STDER1(INT,3,1) + DERX0*DERY0*DERZ1
         STDER0(INT,2)   = STDER0(INT,2)   -(DERX2*DERY0*DERZ0
     &                                     + DERX0*DERY2*DERZ0
     &                                     + DERX0*DERY0*DERZ2)*DP5
         STDER1(INT,1,2) = STDER1(INT,1,2) -(DERX3*DERY0*DERZ0
     &                                     + DERX1*DERY2*DERZ0
     &                                     + DERX1*DERY0*DERZ2)*DP5
         STDER1(INT,2,2) = STDER1(INT,2,2) -(DERX2*DERY1*DERZ0
     &                                     + DERX0*DERY3*DERZ0
     &                                     + DERX0*DERY1*DERZ2)*DP5
         STDER1(INT,3,2) = STDER1(INT,3,2) -(DERX2*DERY0*DERZ1
     &                                     + DERX0*DERY2*DERZ1
     &                                     + DERX0*DERY0*DERZ3)*DP5
         IF (SECDER) THEN
            STDER2(INT,1,1) = STDER2(INT,1,1) + DERX2*DERY0*DERZ0
            STDER2(INT,2,1) = STDER2(INT,2,1) + DERX1*DERY1*DERZ0
            STDER2(INT,3,1) = STDER2(INT,3,1) + DERX1*DERY0*DERZ1
            STDER2(INT,4,1) = STDER2(INT,4,1) + DERX0*DERY2*DERZ0
            STDER2(INT,5,1) = STDER2(INT,5,1) + DERX0*DERY1*DERZ1
            STDER2(INT,6,1) = STDER2(INT,6,1) + DERX0*DERY0*DERZ2
            STDER2(INT,1,2) = STDER2(INT,1,2) -(DERX4*DERY0*DERZ0
     &                                        + DERX2*DERY2*DERZ0
     &                                        + DERX2*DERY0*DERZ2)*DP5
            STDER2(INT,2,2) = STDER2(INT,2,2) -(DERX3*DERY1*DERZ0
     &                                        + DERX1*DERY3*DERZ0
     &                                        + DERX1*DERY1*DERZ2)*DP5
            STDER2(INT,3,2) = STDER2(INT,3,2) -(DERX3*DERY0*DERZ1
     &                                        + DERX1*DERY2*DERZ1
     &                                        + DERX1*DERY0*DERZ3)*DP5
            STDER2(INT,4,2) = STDER2(INT,4,2) -(DERX2*DERY2*DERZ0
     &                                        + DERX0*DERY4*DERZ0
     &                                        + DERX0*DERY2*DERZ2)*DP5
            STDER2(INT,5,2) = STDER2(INT,5,2) -(DERX2*DERY1*DERZ1
     &                                        + DERX0*DERY3*DERZ1
     &                                        + DERX0*DERY1*DERZ3)*DP5
            STDER2(INT,6,2) = STDER2(INT,6,2) -(DERX2*DERY0*DERZ2
     &                                        + DERX0*DERY2*DERZ2
     &                                        + DERX0*DERY0*DERZ4)*DP5
         END IF
C
C     **************************************************
C     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C     **************************************************
C
         MAXT = LVALA + LVALB + MAXADD
         MAXU = MVALA + MVALB + MAXADD
         MAXV = NVALA + NVALB + MAXADD
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            GV = ODC(NVALA,NVALB,IV,2,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               GU = ODC(MVALA,MVALB,IU,2,0,2)
               EE = EU*EV
               FE = FU*EV
               GE = GU*EV
               EF = EU*FV
               FF = FU*FV
               EG = EU*GV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  EEE = ET*EE
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  IADR00 = IADRAU + IT
                  IADR0T = IADR00 + 1
                  IADR0U = IADR00 + ISTEPU
                  IADR0V = IADR00 + ISTEPV
                  IF (SECDER) THEN
                     GT = ODC(LVALA,LVALB,IT,2,0,1)
                     FFE = FT*FE
                     FEF = FT*EF
                     EFF = ET*FF
                     GEE = GT*EE
                     EGE = ET*GE
                     EEG = ET*EG
                     IADRTT = IADR0T + 1
                     IADRTU = IADR0T + ISTEPU
                     IADRTV = IADR0T + ISTEPV
                     IADRUU = IADR0U + ISTEPU
                     IADRUV = IADR0U + ISTEPV
                     IADRVV = IADR0V + ISTEPV
                  END IF
                  IADD = - NAHGTF
C
C                 ***** Loop over nuclei *****
C
                  DO 500 I = 1, NATOMC
C
C                    Pick up HGTF integrals
C
                     IADD = IADD + NAHGTF
                     AH00 = AHGTF(IADR00 + IADD)
                     AH0T = AHGTF(IADR0T + IADD)
                     AH0U = AHGTF(IADR0U + IADD)
                     AH0V = AHGTF(IADR0V + IADD)
C
C                    Multiply by expansion coefficients
C                    and add to appropriate CGTF integral
C
C                    Undifferentiated integral:
C
                     ADER(INT,I,IA0000) = ADER(INT,I,IA0000) + EEE*AH00
C
C                    A differentiated integrals:
C
                     ADER(INT,I,IA0X00) = ADER(INT,I,IA0X00) + FEE*AH00
                     ADER(INT,I,IA0Y00) = ADER(INT,I,IA0Y00) + EFE*AH00
                     ADER(INT,I,IA0Z00) = ADER(INT,I,IA0Z00) + EEF*AH00
C
C                    C differentiated integrals:
C
                     ADER(INT,I,IA000X) = ADER(INT,I,IA000X) - EEE*AH0T
                     ADER(INT,I,IA000Y) = ADER(INT,I,IA000Y) - EEE*AH0U
                     ADER(INT,I,IA000Z) = ADER(INT,I,IA000Z) - EEE*AH0V
C
C                    Second derivatives
C
                     IF (SECDER) THEN
                        AHTT = AHGTF(IADRTT + IADD)
                        AHTU = AHGTF(IADRTU + IADD)
                        AHTV = AHGTF(IADRTV + IADD)
                        AHUU = AHGTF(IADRUU + IADD)
                        AHUV = AHGTF(IADRUV + IADD)
                        AHVV = AHGTF(IADRVV + IADD)
C
C                       A-A differentiated integrals:
C
                        ADER(INT,I,IAXX00) = ADER(INT,I,IAXX00)+GEE*AH00
                        ADER(INT,I,IAXY00) = ADER(INT,I,IAXY00)+FFE*AH00
                        ADER(INT,I,IAXZ00) = ADER(INT,I,IAXZ00)+FEF*AH00
                        ADER(INT,I,IAYY00) = ADER(INT,I,IAYY00)+EGE*AH00
                        ADER(INT,I,IAYZ00) = ADER(INT,I,IAYZ00)+EFF*AH00
                        ADER(INT,I,IAZZ00) = ADER(INT,I,IAZZ00)+EEG*AH00
C
C                       A-C differentiated integrals:
C
                        ADER(INT,I,IA0X0X) = ADER(INT,I,IA0X0X)-FEE*AH0T
                        ADER(INT,I,IA0X0Y) = ADER(INT,I,IA0X0Y)-FEE*AH0U
                        ADER(INT,I,IA0X0Z) = ADER(INT,I,IA0X0Z)-FEE*AH0V
                        ADER(INT,I,IA0Y0X) = ADER(INT,I,IA0Y0X)-EFE*AH0T
                        ADER(INT,I,IA0Y0Y) = ADER(INT,I,IA0Y0Y)-EFE*AH0U
                        ADER(INT,I,IA0Y0Z) = ADER(INT,I,IA0Y0Z)-EFE*AH0V
                        ADER(INT,I,IA0Z0X) = ADER(INT,I,IA0Z0X)-EEF*AH0T
                        ADER(INT,I,IA0Z0Y) = ADER(INT,I,IA0Z0Y)-EEF*AH0U
                        ADER(INT,I,IA0Z0Z) = ADER(INT,I,IA0Z0Z)-EEF*AH0V
C
C                       C-C differentiated integrals:
C
                        ADER(INT,I,IA00XX) = ADER(INT,I,IA00XX)+EEE*AHTT
                        ADER(INT,I,IA00XY) = ADER(INT,I,IA00XY)+EEE*AHTU
                        ADER(INT,I,IA00XZ) = ADER(INT,I,IA00XZ)+EEE*AHTV
                        ADER(INT,I,IA00YY) = ADER(INT,I,IA00YY)+EEE*AHUU
                        ADER(INT,I,IA00YZ) = ADER(INT,I,IA00YZ)+EEE*AHUV
                        ADER(INT,I,IA00ZZ) = ADER(INT,I,IA00ZZ)+EEE*AHVV
                     END IF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck dipint */
       SUBROUTINE DIPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SINT0,
     &                   DINT1,SHGTF,CORPX,CORPY,CORPZ,DIFDIP)
C
C     tuh 1985
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
      LOGICAL DIFDIP
      DIMENSION SINT0(KCKTAB), DINT1(KCKTAB,3,3)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <onecom.h>
#include <lmns.h>
      INT1 = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT1 = INT1 + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,1,0,0,1) + CORPX*SX0
         DY0 = SHGTF*ODC(MVALA,MVALB,1,0,0,2) + CORPY*SY0
         DZ0 = SHGTF*ODC(NVALA,NVALB,1,0,0,3) + CORPZ*SZ0
C
         SINT0(INT1)   = SINT0(INT1)   + SX0*SY0*SZ0
         IF (DIFDIP) THEN
            IF (.NOT. ONECEN) THEN
               SX1 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
               SY1 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
               SZ1 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
               DX1 = SHGTF*ODC(LVALA,LVALB,1,1,0,1) + CORPX*SX1
               DY1 = SHGTF*ODC(MVALA,MVALB,1,1,0,2) + CORPY*SY1
               DZ1 = SHGTF*ODC(NVALA,NVALB,1,1,0,3) + CORPZ*SZ1
               DINT1(INT1,1,1) = DINT1(INT1,1,1) + DX1*SY0*SZ0
               DINT1(INT1,2,1) = DINT1(INT1,2,1) + SX1*DY0*SZ0
               DINT1(INT1,3,1) = DINT1(INT1,3,1) + SX1*SY0*DZ0
               DINT1(INT1,1,2) = DINT1(INT1,1,2) + DX0*SY1*SZ0
               DINT1(INT1,2,2) = DINT1(INT1,2,2) + SX0*DY1*SZ0
               DINT1(INT1,3,2) = DINT1(INT1,3,2) + SX0*SY1*DZ0
               DINT1(INT1,1,3) = DINT1(INT1,1,3) + DX0*SY0*SZ1
               DINT1(INT1,2,3) = DINT1(INT1,2,3) + SX0*DY0*SZ1
               DINT1(INT1,3,3) = DINT1(INT1,3,3) + SX0*SY0*DZ1
            END IF
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck cintab */
      SUBROUTINE CINTAB(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,EXPA,EXPB,
     &                  EXPDER,SHGTF,AHGTF,NATOMC)
C
C     TUH
C
#include <implicit.h>
#include <maxaqn.h>
#include <mxcent.h>
      PARAMETER (DP5 = 0.5 D00, D1 = 1.0D00, DP75 = 0.75 D00)
      DIMENSION EXPDER(KCKTAB,2,2), AHGTF(*)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <onecom.h>
#include <ader.h>
#include <lmns.h>
C
      AINV = D1/EXPA
      BINV = D1/EXPB
C
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
C    ***************************************
C    ***** CALCULATE OVERLAP INTEGRALS *****
C    ***************************************
C
         S0X = SHGTF*ODC(LVALA,  LVALB  ,0,0,0,1)
         S0Y = SHGTF*ODC(MVALA,  MVALB  ,0,0,0,2)
         S0Z = SHGTF*ODC(NVALA,  NVALB  ,0,0,0,3)
         SAX = SHGTF*ODC(LVALA+2,LVALB  ,0,0,0,1)
         SAY = SHGTF*ODC(MVALA+2,MVALB  ,0,0,0,2)
         SAZ = SHGTF*ODC(NVALA+2,NVALB  ,0,0,0,3)
         SBX = SHGTF*ODC(LVALA,  LVALB+2,0,0,0,1)
         SBY = SHGTF*ODC(MVALA,  MVALB+2,0,0,0,2)
         SBZ = SHGTF*ODC(NVALA,  NVALB+2,0,0,0,3)
C
         AFAC = (DP5*(LVALA + MVALA + NVALA) + DP75)*AINV
         BFAC = (DP5*(LVALB + MVALB + NVALB) + DP75)*BINV
         EXPDER(INT,1,1) = EXPDER(INT,1,1) + AFAC*S0X*S0Y*S0Z 
     &                   - SAX*S0Y*S0Z - S0X*SAY*S0Z - S0X*S0Y*SAZ 
         EXPDER(INT,2,1) = EXPDER(INT,2,1) + BFAC*S0X*S0Y*S0Z
     &                   - SBX*S0Y*S0Z - S0X*SBY*S0Z - S0X*S0Y*SBZ 
C
C     ***************************************
C     ***** CALCULATE KINETIC INTEGRALS *****
C     ***************************************
C
         THGTF = - DP5*SHGTF
         T0X = THGTF*ODC(LVALA,  LVALB  ,0,2,0,1)
         T0Y = THGTF*ODC(MVALA,  MVALB  ,0,2,0,2)
         T0Z = THGTF*ODC(NVALA,  NVALB  ,0,2,0,3)
         TAX = THGTF*ODC(LVALA+2,LVALB  ,0,2,0,1)
         TAY = THGTF*ODC(MVALA+2,MVALB  ,0,2,0,2)
         TAZ = THGTF*ODC(NVALA+2,NVALB  ,0,2,0,3)
         TBX = THGTF*ODC(LVALA,  LVALB+2,0,2,0,1)
         TBY = THGTF*ODC(MVALA,  MVALB+2,0,2,0,2)
         TBZ = THGTF*ODC(NVALA,  NVALB+2,0,2,0,3)
C
         AINT = (T0X*S0Y*S0Z + S0X*T0Y*S0Z + S0X*S0Y*T0Z)*AFAC
     &        - (TAX*S0Y*S0Z + SAX*T0Y*S0Z + SAX*S0Y*T0Z)
     &        - (T0X*SAY*S0Z + S0X*TAY*S0Z + S0X*SAY*T0Z)
     &        - (T0X*S0Y*SAZ + S0X*T0Y*SAZ + S0X*S0Y*TAZ)
         BINT = (T0X*S0Y*S0Z + S0X*T0Y*S0Z + S0X*S0Y*T0Z)*BFAC
     &        - (TBX*S0Y*S0Z + SBX*T0Y*S0Z + SBX*S0Y*T0Z)
     &        - (T0X*SBY*S0Z + S0X*TBY*S0Z + S0X*SBY*T0Z)
     &        - (T0X*S0Y*SBZ + S0X*T0Y*SBZ + S0X*S0Y*TBZ)
C
C     **************************************************
C     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C     **************************************************
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO IV = 0, MAXV + 2
            E0Z = ODC(NVALA,  NVALB  ,IV,0,0,3)
            EAZ = ODC(NVALA+2,NVALB  ,IV,0,0,3)
            EBZ = ODC(NVALA,  NVALB+2,IV,0,0,3)
            IADRAU = IADRAV
            DO IU = 0, MAXU + 2
               E0Y = ODC(MVALA,  MVALB  ,IU,0,0,2)
               EAY = ODC(MVALA+2,MVALB  ,IU,0,0,2)
               EBY = ODC(MVALA,  MVALB+2,IU,0,0,2)
               DO IT = 0, MAXT + 2
               IF (IT+IU+IV .LE. MAXT+MAXU+MAXV+2) THEN
                  E0X = ODC(LVALA,  LVALB  ,IT,0,0,1)
                  EAX = ODC(LVALA+2,LVALB  ,IT,0,0,1)
                  EBX = ODC(LVALA,  LVALB+2,IT,0,0,1)
                  E0XY = E0X*E0Y
                  E0XZ = E0X*E0Z
                  E0YZ = E0Y*E0Z
                  IF      (IT.GT.MAXT) THEN
                     AE = - EAX*E0YZ 
                     BE = - EBX*E0YZ 
                  ELSE IF (IU.GT.MAXU) THEN
                     AE = - EAY*E0XZ 
                     BE = - EBY*E0XZ 
                  ELSE IF (IV.GT.MAXV) THEN
                     AE = - EAZ*E0XY
                     BE = - EBZ*E0XY
                  ELSE
                     AE = AFAC*E0XY*E0Z - EAX*E0YZ - EAY*E0XZ - EAZ*E0XY
                     BE = BFAC*E0XY*E0Z - EBX*E0YZ - EBY*E0XZ - EBZ*E0XY
                  END IF
                  IADR00 = IADRAU + IT
                  IADD = - NAHGTF
                  DO IATOM = 1,NATOMC
                     IADD = IADD + NAHGTF
                     AINT = AINT + AE*AHGTF(IADR00 + IADD)
                     BINT = BINT + BE*AHGTF(IADR00 + IADD)
                  END DO
               END IF
               END DO
               IADRAU = IADRAU + ISTEPU
            END DO
            IADRAV = IADRAV + ISTEPV
         END DO
         EXPDER(INT,1,2) = EXPDER(INT,1,2) + AINT
         EXPDER(INT,2,2) = EXPDER(INT,2,2) + BINT
  100 CONTINUE
      RETURN
      END
