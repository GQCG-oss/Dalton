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
C /* Deck dftmag */
      SUBROUTINE DFTMAG(EXCMAT,COORX,COORY,COORZ,GAO,GAO1,GAB1,
     &                  GAB2,VXC,VXB,RH,DOGGA,FROMVX)
C
C     T. Helgaker sep 99/ oct 00 / feb 01
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C     differentiated with respect to magnetic field
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
      PARAMETER (D2 = 2.0D0)
C
#include <inforb.h>
#include <nuclei.h>
#include <dfterg.h>
#include <energy.h>
#include <dftcom.h>
#include <dftinf.h>
#include <orgcom.h>
#include <symmet.h>
C
      LOGICAL FROMVX, DOGGA
      DIMENSION GAO(NBAST), GAO1(NBAST,3), 
     &          GAB1(NBAST,3), GAB2(NBAST,3,3), 
     &          EXCMAT(NBAST,NBAST,3), RH(3)
C
      POX = COORX - ORIGIN(1)
      POY = COORY - ORIGIN(2)
      POZ = COORZ - ORIGIN(3)
C
      DO K = 1, 3
         KSYM = ISYMAX(K,2) + 1
         IF (K.EQ.1) THEN
            KY = 2
            KZ = 3
            PY = POY
            PZ = POZ
            RY = RH(2)
            RZ = RH(3)
         ELSE IF (K.EQ.2) THEN
            KY = 3
            KZ = 1
            PY = POZ
            PZ = POX
            RY = RH(3)
            RZ = RH(1)
         ELSE
            KY = 1
            KZ = 2
            PY = POX
            PZ = POY
            RY = RH(1)
            RZ = RH(2)
         END IF 
         IF (.NOT.(DOGGA .AND. .NOT.FROMVX)) THEN
            G10Z = D2*VXC*PY
            G10Y = D2*VXC*PZ
            DO JSYM = 1, NSYM
               ISYM = MULD2H(JSYM,KSYM) 
               JSTR = IBAS(JSYM) + 1
               JEND = IBAS(JSYM) + NBAS(JSYM)
               ISTR = IBAS(ISYM) + 1
               IEND = IBAS(ISYM) + NBAS(ISYM)
               DO J = JSTR, JEND
                  FYZ = G10Y*GAB1(J,KY) - G10Z*GAB1(J,KZ)
                  DO I = ISTR, IEND
                     EXCMAT(I,J,K) = EXCMAT(I,J,K) + FYZ*GAO(I)
                  END DO
               END DO
            END DO
         ELSE 
            G10Z = D2*VXC*PY + VXB*RY
            G10Y = D2*VXC*PZ + VXB*RZ
            G21Z = VXB*PY*RH(1)
            G21Y = VXB*PZ*RH(1)
            G22Z = VXB*PY*RH(2)
            G22Y = VXB*PZ*RH(2)
            G23Z = VXB*PY*RH(3)
            G23Y = VXB*PZ*RH(3)
            IF (NSYM.EQ.1) THEN
               DO J = 1, NBAST
                  FYZ = G10Y*GAB1(J,  KY) - G10Z*GAB1(J,  KZ)
     &                + G21Y*GAB2(J,1,KY) - G21Z*GAB2(J,1,KZ)
     &                + G22Y*GAB2(J,2,KY) - G22Z*GAB2(J,2,KZ)
     &                + G23Y*GAB2(J,3,KY) - G23Z*GAB2(J,3,KZ)
                  GD = VXB*(RH(1)*GAO1(J,1)+RH(2)*GAO1(J,2)
     $                 +RH(3)*GAO1(J,3))
                  FY = PY*GD
                  FZ = PZ*GD
                  DO I = 1, NBAST 
                     EXCMAT(I,J,K) = EXCMAT(I,J,K) + FYZ*GAO(I)
     &                             + FY*GAB1(I,KZ) - FZ*GAB1(I,KY)
                  END DO
               END DO
            ELSE
               DO JSYM = 1, NSYM
                  ISYM = MULD2H(JSYM,KSYM) 
                  JSTR = IBAS(JSYM) + 1
                  JEND = IBAS(JSYM) + NBAS(JSYM)
                  ISTR = IBAS(ISYM) + 1
                  IEND = IBAS(ISYM) + NBAS(ISYM)
                  DO J = JSTR, JEND
                     FYZ = G10Y*GAB1(J,  KY) - G10Z*GAB1(J,  KZ)
     &                   + G21Y*GAB2(J,1,KY) - G21Z*GAB2(J,1,KZ)
     &                   + G22Y*GAB2(J,2,KY) - G22Z*GAB2(J,2,KZ)
     &                   + G23Y*GAB2(J,3,KY) - G23Z*GAB2(J,3,KZ)
                     GD = VXB*(RH(1)*GAO1(J,1) + RH(2)*GAO1(J,2) 
     &                                       + RH(3)*GAO1(J,3))
                     FY = PY*GD
                     FZ = PZ*GD
                     DO I = ISTR, IEND 
                        EXCMAT(I,J,K) = EXCMAT(I,J,K) + FYZ*GAO(I)
     &                                + FY*GAB1(I,KZ) - FZ*GAB1(I,KY)
                     END DO
                  END DO
               END DO
            END IF
         END IF
      END DO
      RETURN
      END
