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
C /* Deck dfthed */
      SUBROUTINE DFTHED(JATOM,FCDX,FCDY,FCDZ,FDFTDX,FDFTDY,FDFTDZ,
     &                WORK,LWORK,IPRINT)
C
C     add DFT component of Fock derivative matrix
C
C     O B Lutn?s, D Wilson  - Jan 04
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <priunit.h>
#include <dummy.h>
#include <nuclei.h>
#include <inforb.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0)
      LOGICAL FOUND
      INTEGER P, Q, R, S, U, V, RS
      CHARACTER*4  KEY
      DIMENSION FCDX(N2BASX), FCDY(N2BASX), FCDZ(N2BASX), WORK(LWORK),
     &          FDFTDX(N2BASX),FDFTDY(N2BASX),FDFTDZ(N2BASX)
#include <symmet.h>
#include <symind.h>
#include <dftcom.h>
#include "dftmolhes.h"
C
      EXTERNAL DFTKSD
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      DOGGA = DFT_ISGGA()
C setting IATOM in commonblock dftmolhes.h
      IATOM = JATOM
      IF (IPRINT .GE. 5) CALL TITLER('Output from DFTHED','*',103)
C
      KTMAT  = 1
      KQMAT  = KTMAT  + 3*N2BASX
      KDMAT  = KQMAT  +   N2BASX
      KV1MAT = KDMAT  +   N2BASX
      KFCD   = KV1MAT + 3*N2BASX
      KDFTFD = KFCD   + 3*N2BASX
      KSOMO  = KDFTFD + 3*N2BASX
      KCMO   = KSOMO  +   N2BASX
      KLAST  = KCMO   +   NCMOT
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTHED',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
C Construct V1MAT = Sa * DMAT in AO basis from Sa in MO basis
C
      CALL GETV1MAT(IATOM,.TRUE.,WORK(KCMO),WORK(KDMAT),WORK(KQMAT),
     &              WORK(KV1MAT),WORK(KTMAT),WORK(KLAST),LWRK,IPRINT)
C
      CALL DZERO(WORK(KFCD),3*N2BASX)
      CALL DZERO(WORK(KDFTFD),3*N2BASX)
C
      KACTORB = 0
      KNOORB = 0
      IVMATS = 0
      CALL DFTINT(WORK(KDMAT),1,1,.FALSE.,.FALSE.,WORK(KLAST+1),LWRK,
     &            DFTKSD,WORK(KDMAT),ELE)
      KACTORB = KACTORB + KNOORB
C
      IF (IPRINT.GT.10) THEN
         CALL HEADER('DFT FCD Fock derivative matrix (x direction)',-1)
         CALL OUTPUT(WORK(KFCD),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('DFT FCD Fock derivative matrix (y direction)',-1)
         CALL OUTPUT(WORK(KFCD+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &        1,LUPRI)
         CALL HEADER('DFT FCD Fock derivative matrix (z direction)',-1)
         CALL OUTPUT(WORK(KFCD+2*N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &        1,LUPRI)
         CALL HEADER('DFT FD matrix SO basis (x direction)',-1)
         CALL OUTPUT(WORK(KDFTFD),1,NBAST,1,NBAST,NBAST,NBAST,
     &        1,LUPRI)
         CALL HEADER('DFT FD matrix SO basis (y direction)',-1)
         CALL OUTPUT(WORK(KDFTFD+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &        1,LUPRI)
         CALL HEADER('DFT FD matrix SO basis (z direction)',-1)
         CALL OUTPUT(WORK(KDFTFD+2*N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &        1,LUPRI)
      END IF
C
C  Transforming FDFTD to MO basis
C
      DO I = 1,3
         CALL DZERO(WORK(KSOMO),N2BASX)
         CALL DCOPY(N2BASX,WORK(KDFTFD+(I-1)*N2BASX),1,WORK(KSOMO),1)
         CALL ONETRA(WORK(KCMO),WORK(KDFTFD+(I-1)*N2BASX),WORK(KSOMO),
     &        WORK(KLAST),LWRK,IATOM,I,IPRINT)
      END DO
      IF (IPRINT.GT.10) THEN
         CALL HEADER('DFT FD MO basis (x direction)',-1)
         CALL OUTPUT(WORK(KDFTFD),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &        LUPRI)  
         CALL HEADER('DFT FD MO basis (y direction)',-1)
         CALL OUTPUT(WORK(KDFTFD+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &        LUPRI)  
         CALL HEADER('DFT FD MO basis (z direction)',-1)
         CALL OUTPUT(WORK(KDFTFD+2*N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &        1,LUPRI)
      END IF
C     
      KXC    = KFCD   - 1
      KYC    = KFCD   - 1 +   N2BASX
      KZC    = KFCD   - 1 + 2*N2BASX
      KDFTFX = KDFTFD - 1 
      KDFTFY = KDFTFD - 1 +   N2BASX
      KDFTFZ = KDFTFD - 1 + 2*N2BASX
      DO I = 1, N2BASX
         FCDX(I)   = FCDX(I)   +  WORK(KXC + I)
         FCDY(I)   = FCDY(I)   +  WORK(KYC + I)
         FCDZ(I)   = FCDZ(I)   +  WORK(KZC + I)
         FDFTDX(I) = FDFTDX(I) +  WORK(KDFTFX + I)
         FDFTDY(I) = FDFTDY(I) +  WORK(KDFTFY + I)
         FDFTDZ(I) = FDFTDZ(I) +  WORK(KDFTFZ + I)
      END DO
C     
      RETURN
      END  
C    
C /* Deck dftksd */
      SUBROUTINE DFTKSD(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                  COORD,WGHT,DATMAT)
C
C     This subroutine calculates the part of the geometrical derivatives
C     of the Kohn-Sham matrix, needed for molecular hessians. 
C        
C     O. B. Lutn?s and D. Wilson Jan 2004
C    
C     DATMAT(KDFCD) for static contributions to molecular hessians,
C                        that is, contribution to the derivative of the
C                        generalized fock matrix. Reorthonormalization
C                        contributions are calculated elsewhere
C                        (in ABARHS-ABATR1). 
C        
C     DATMAT(KDFTFD) remaining term for total derivative fock matrix
C                        for calculation of diffferentiated electronic
C                        gradients (needed for response calculation).
C    
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <inforb.h>
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &         RHOA(NBLEN), GRADA(3,NBLEN), NBLCNT(8),
     &         NBLOCKS(2,LDAIB,8), DATMAT(10*N2BASX)
C
      KDMAT  = 1
      KV1MAT = KDMAT   +   N2BASX
      KDFCD  = KV1MAT  + 3*N2BASX
      KDFTFD = KDFCD   + 3*N2BASX
      KLAST  = KDFTFD  + 3*N2BASX
      CALL DFTKSD1(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,COORD,
     &             WGHT, DATMAT(KDMAT), DATMAT(KV1MAT),
     &             DATMAT(KDFCD), DATMAT(KDFTFD))
C
      RETURN
      END
C
C /* Deck DFTKSD1 */
      SUBROUTINE DFTKSD1(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                   COORD,WGHT,DMAT,V1MAT,DFTFCD,DFTFDX)
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <symmet.h>
#include <shells.h>
#include <nuclei.h>
#include <inforb.h>
#include "dftmolhes.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0,
     &           D2 = 2.0D0, DP5 = 0.5D0)
      LOGICAL ACTIVE(NBAST)
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN), WGHT(NBLEN),
     &      RHOA(NBLEN), GRADA(3,NBLEN), NBLCNT(8), NBLOCKS(2,LDAIB,8)
      DIMENSION DMAT(NBAST,NBAST), V1MAT(NBAST,NBAST,3),
     &          GRADAT(NBLEN,3)
      DIMENSION DFTFCD(NBAST,NBAST,3), DFTFDX(NBAST,NBAST,3)
      DIMENSION RHODER(NBLEN,3*NUCDEP), ZETDER(NBLEN,3*NUCDEP),
     &          GDRHOD(NBLEN,3*NUCDEP,3)
      DIMENSION TMP1(NBLEN,3*NUCDEP), TMP2(NBLEN,3*NUCDEP)
      DIMENSION GD1(NBLEN), GD2(NBLEN), GD3(NBLEN), GD4(NBLEN),
     &          TMP3(NBLEN), TMP4(NBLEN), TMP5(NBLEN)
      DIMENSION VXCR(NBLEN), VXCRR(NBLEN), DENS(5), VX(5),
     &          VXCZ(NBLEN), VXCRZ(NBLEN), VXCZZ(NBLEN)
      DIMENSION IACTCOR(NBAST), IXY2(9)
      DATA IXY2 / 5, 6, 7, 6, 8, 9, 7, 9, 10/
#include <ibtfun.h>
C
C Set active orbitals
C
      DO I = 1, NBAST
         ACTIVE(I) = .FALSE.
      END DO
      DO ISYM = 1, NSYM
         DO IBLA = 1, NBLCNT(ISYM)
            DO I = NBLOCKS(1,IBLA,ISYM), NBLOCKS(2,IBLA,ISYM)
               ACTIVE(I) = .TRUE.
            END DO
         END DO
      END DO
C
C Setting up array of active orbitals centered on IATOM
C
      IA = 0   
      KACTCOR = 0
      ISCORN = IPTCNT(3*(IATOM - 1) + 1,0,1)
      DO ISHELA = 1, KMAX  
         ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + 1,0,1)
         DO ICOMPA = 1 ,KHKT(ISHELA)
            IA = IA + 1
            IF (ACTIVE(IA).AND.ISCORN.EQ.ISCORA) THEN
               KACTCOR = KACTCOR + 1
               IACTCOR(KACTCOR) = IA
            END IF
         END DO
      END DO
      IF (KACTCOR.EQ.0) THEN
         KNOORB = KNOORB + 1
      ELSE
         KACTORB = KACTORB + 1
      END IF
C
C If no active orbitals on IATOM, skip routine
C
      IF (KACTCOR.EQ.0) GOTO 100
C
C Construct transpose of GRADA
C
      DO I = 1, 3
      DO J = 1, NBLEN
         GRADAT(J,I) = GRADA(I,J)
      END DO
      END DO     
C     
C DFT potential 
C    
      CALL DZERO(DENS,5)
      CALL DZERO(VX,5)
      IF (DOGGA) THEN
         DO I = 1, NBLEN
            GRDNRM = SQRT(GRADAT(I,1)**2+GRADAT(I,2)**2+GRADAT(I,3)**2)
            DENS(1) = DP5*RHOA(I)
            DENS(2) = DP5*RHOA(I)
            DENS(3) = DP5*GRDNRM 
            DENS(4) = DP5*GRDNRM 
            DENS(5) = DENS(3)*DENS(4)
            CALL DFTPOT1(VX,WGHT(I),DENS,.FALSE.)
            GDNRMI   = 1.0D0/GRDNRM 
            VXCR(I)  = VX(1)
            VXCZ(I)  = VX(2)*GDNRMI
            VXCRR(I) = VX(3)
            VXCRZ(I) = VX(4)*GDNRMI
            VXCZZ(I) = VX(5)*DP5*(GDNRMI**2)
         END DO
      ELSE
         DO I = 1, NBLEN
           DENS(1) = DP5*RHOA(I)
           DENS(2) = DP5*RHOA(I)
           CALL DFTPOT1(VX,WGHT(I),DENS(1),.FALSE.)
           VXCR(I)  = VX(1)
           VXCRR(I) = VX(3)
         END DO
      END IF
C 
C Zero arrays
C
      DO IX = 1, 3
         ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
         CALL DZERO(RHODER(1,ISCORN),NBLEN)
         IF (DOGGA)THEN
            CALL DZERO(GDRHOD(1,ISCORN,1),NBLEN)
            CALL DZERO(GDRHOD(1,ISCORN,2),NBLEN)
            CALL DZERO(GDRHOD(1,ISCORN,3),NBLEN)
         END IF
      END DO
C
C construct rho and zeta derivatives
C
      IF (.NOT. DOGGA) THEN
C lda
         DO ISYM = 1, NSYM
         DO LAC = 1, KACTCOR
            IA = IACTCOR(LAC)
            CALL DZERO(GD1,NBLEN)
            DO IBLB = 1, NBLCNT(ISYM)
               ISTRB = NBLOCKS(1,IBLB,ISYM)
               IENDB = NBLOCKS(2,IBLB,ISYM)
               DO IB = ISTRB, IENDB
                  CALL DAXPY(NBLEN,DMAT(IB,IA),GAO(1,IB,1),1,GD1(1),1)
               END DO
            END DO
            DO IX = 1, 3
               ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
               DO I = 1,NBLEN
                  RHODER(I,ISCORN) = RHODER(I,ISCORN)
     &                    - GAO(I,IA,IX+1)*GD1(I)
               END DO
            END DO
         END DO
         END DO
      ELSE
C  gga
         DO ISYM = 1, NSYM
            DO LAC = 1, KACTCOR
               IA = IACTCOR(LAC)
               CALL DZERO(GD1,NBLEN)
               CALL DZERO(GD2,NBLEN)
               CALL DZERO(GD3,NBLEN)
               CALL DZERO(GD4,NBLEN)
               DO IBLB = 1, NBLCNT(ISYM)
                  ISTRB = NBLOCKS(1,IBLB,ISYM)
                  IENDB = NBLOCKS(2,IBLB,ISYM)
                  DO IB = ISTRB, IENDB
                     DBA = DMAT(IB,IA)
                     CALL DAXPY(NBLEN,DBA,GAO(1,IB,1),1,GD1(1),1)
                     CALL DAXPY(NBLEN,DBA,GAO(1,IB,2),1,GD2(1),1)
                     CALL DAXPY(NBLEN,DBA,GAO(1,IB,3),1,GD3(1),1)
                     CALL DAXPY(NBLEN,DBA,GAO(1,IB,4),1,GD4(1),1)
                  END DO
               END DO
               DO IX = 1, 3
                  IP1 = IX + 1
                  IXA1 = IXY2(3*IX - 2)
                  IXA2 = IXY2(3*IX - 1)
                  IXA3 = IXY2(3*IX)
                  ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
                  DO I = 1, NBLEN
                     RHODER(I,ISCORN) = RHODER(I,ISCORN)
     &                                - GD1(I)*GAO(I,IA,IP1)
                     GDRHOD(I,ISCORN,1) = GDRHOD(I,ISCORN,1)
     &                                  - GD1(I)*GAO(I,IA,IXA1)
     &                                  - GD2(I)*GAO(I,IA,IP1)
                     GDRHOD(I,ISCORN,2) = GDRHOD(I,ISCORN,2)
     &                                  - GD1(I)*GAO(I,IA,IXA2)
     &                                  - GD3(I)*GAO(I,IA,IP1)
                     GDRHOD(I,ISCORN,3) = GDRHOD(I,ISCORN,3)
     &                                  - GD1(I)*GAO(I,IA,IXA3)
     &                                  - GD4(I)*GAO(I,IA,IP1)
                  END DO 
               END DO
            END DO
         END DO
      END IF
C        
      DO IX = 1, 3
         ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
         CALL DSCAL(NBLEN,D2,RHODER(1,ISCORN),1)
         IF (DOGGA) THEN 
            CALL DSCAL(NBLEN,D2,GDRHOD(1,ISCORN,1),1)
            CALL DSCAL(NBLEN,D2,GDRHOD(1,ISCORN,2),1)
            CALL DSCAL(NBLEN,D2,GDRHOD(1,ISCORN,3),1)
         END IF
      END DO
C           
C Contribution to derivatives from connection matrix
C              
C lda             
      IF (.NOT. DOGGA) THEN 
         DO ISYM = 1, NSYM
            DO IBLA = 1, NBLCNT(ISYM)
               ISTRA = NBLOCKS(1,IBLA,ISYM)
               IENDA = NBLOCKS(2,IBLA,ISYM)
               DO IA = ISTRA, IENDA
                  DO IXA = 1, 3
                     ISCORN = IPTCNT(3*(IATOM - 1) + IXA,0,1)
                     DO IBLB = 1, NBLCNT(ISYM)
                        ISTRB = NBLOCKS(1,IBLB,ISYM)
                        IENDB = NBLOCKS(2,IBLB,ISYM)
                        DO IB = ISTRB, IENDB
                           VBA = V1MAT(IB,IA,IXA)
                           DO I = 1, NBLEN
                              RHODER(I,ISCORN) = RHODER(I,ISCORN)
     &                          - GAO(I,IB,1)*VBA*GAO(I,IA,1)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO 
            END DO   
         END DO      
      ELSE        
C gga          
         DO ISYM = 1, NSYM
            DO IX = 1, 3
               ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
               DO IBLA = 1, NBLCNT(ISYM)
                  ISTRA = NBLOCKS(1,IBLA,ISYM)
                  IENDA = NBLOCKS(2,IBLA,ISYM)
                  DO IA = ISTRA, IENDA
                     CALL DZERO(TMP3,NBLEN)
                     CALL DZERO(TMP4,NBLEN)
                     DO IBLB = 1, NBLCNT(ISYM)
                        ISTRB = NBLOCKS(1,IBLB,ISYM)
                        IENDB = NBLOCKS(2,IBLB,ISYM)
                        DO IB = ISTRB, IENDB
                           VBA = V1MAT(IB,IA,IX)
                           VAB = V1MAT(IA,IB,IX)
                           DO I = 1, NBLEN
                              TMP3(I) = TMP3(I) - VBA*GAO(I,IB,1)
                              TMP4(I) = TMP4(I) - VAB*GAO(I,IB,1)
                           END DO
                        END DO
                     END DO
                     DO I = 1, NBLEN
                        TMP4(I) = TMP4(I) + TMP3(I)
                        RHODER(I,ISCORN) = RHODER(I,ISCORN)
     &                                   + TMP3(I)*GAO(I,IA,1)
                        GDRHOD(I,ISCORN,1) = GDRHOD(I,ISCORN,1)
     &                                     + TMP4(I)*GAO(I,IA,2)
                        GDRHOD(I,ISCORN,2) = GDRHOD(I,ISCORN,2)
     &                                     + TMP4(I)*GAO(I,IA,3)
                        GDRHOD(I,ISCORN,3) = GDRHOD(I,ISCORN,3)
     &                                     + TMP4(I)*GAO(I,IA,4)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C Nuclear derivative of rho and zeta now finished
C
      IF (DOGGA) THEN
         DO I = 1, NBLEN
            ZETA = GRADAT(I,1)**2 + GRADAT(I,2)**2 + GRADAT(I,3)**2
            TMP3(I) = VXCZZ(I) - DP5*VXCZ(I) / ZETA
         END DO
         DO IX = 1, 3
            ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
            DO I = 1, NBLEN
               ZETDER(I,ISCORN) = D2*(GDRHOD(I,ISCORN,1)*GRADAT(I,1)
     &                              + GDRHOD(I,ISCORN,2)*GRADAT(I,2)
     &                              + GDRHOD(I,ISCORN,3)*GRADAT(I,3))
               TMP1(I,ISCORN) = RHODER(I,ISCORN)*VXCRR(I)
     &                        + ZETDER(I,ISCORN)*VXCRZ(I)*DP5
               TMP2(I,ISCORN) = ZETDER(I,ISCORN)*TMP3(I)
     &                        + RHODER(I,ISCORN)*VXCRZ(I)
               GDRHOD(I,ISCORN,1) = GDRHOD(I,ISCORN,1)*VXCZ(I)
               GDRHOD(I,ISCORN,2) = GDRHOD(I,ISCORN,2)*VXCZ(I)
               GDRHOD(I,ISCORN,3) = GDRHOD(I,ISCORN,3)*VXCZ(I)
            END DO
         END DO
      ELSE
         DO IX = 1, 3
            ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
            DO I = 1,NBLEN
               TMP1(I,ISCORN) = RHODER(I,ISCORN)*VXCRR(I)
            END DO
         END DO
      END IF
C                 
      IF (DOGGA) THEN
         DO ISYM = 1, NSYM
         DO IX = 1, 3
            IP1 = IX + 1
            ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
            DO I = 1, NBLEN
               TMP3(I) = TMP2(I,ISCORN)*GRADAT(I,1) + GDRHOD(I,ISCORN,1)
               TMP4(I) = TMP2(I,ISCORN)*GRADAT(I,2) + GDRHOD(I,ISCORN,2)
               TMP5(I) = TMP2(I,ISCORN)*GRADAT(I,3) + GDRHOD(I,ISCORN,3)
            END DO            
            DO IBLA = 1, NBLCNT(ISYM) 
               ISTRA = NBLOCKS(1,IBLA,ISYM)
               IENDA = NBLOCKS(2,IBLA,ISYM)
               DO IA = ISTRA, IENDA
                  FAC = D0  
                  DO I = 1, NBLEN 
                     GD1(I) = TMP1(I,ISCORN)*GAO(I,IA,1)
     &                        + GAO(I,IA,2)*TMP3(I)
     &                        + GAO(I,IA,3)*TMP4(I)
     &                        + GAO(I,IA,4)*TMP5(I)
                     GD2(I) = GAO(I,IA,1)*TMP3(I)
                     GD3(I) = GAO(I,IA,1)*TMP4(I)
                     GD4(I) = GAO(I,IA,1)*TMP5(I) 
                     FAC = FAC + GD1(I)*GAO(I,IA,1)
     &                         + GD2(I)*GAO(I,IA,2)
     &                         + GD3(I)*GAO(I,IA,3)
     &                         + GD4(I)*GAO(I,IA,4)
                  END DO
                  DFTFDX(IA,IA,IX) = DFTFDX(IA,IA,IX) + FAC
                  DO IBLB = 1, NBLCNT(ISYM)
                     ISTRB = NBLOCKS(1,IBLB,ISYM)
                     IENDB = NBLOCKS(2,IBLB,ISYM)
                     DO IB = ISTRB, IENDB
                        IF (IB.GT.IA) THEN
                           FAC = D0
                           DO I = 1, NBLEN
                              FAC = FAC + GD1(I)*GAO(I,IB,1)
     &                                  + GD2(I)*GAO(I,IB,2)
     &                                  + GD3(I)*GAO(I,IB,3)
     &                                  + GD4(I)*GAO(I,IB,4)
                           END DO
                           DFTFDX(IA,IB,IX) = DFTFDX(IA,IB,IX) + FAC
                        END IF      
                     END DO         
                  END DO
               END DO         
            END DO
            DO IB = 1, NBAST  
               DO IA = IB+1, NBAST
                  DFTFDX(IA,IB,IX) = DFTFDX(IB,IA,IX)
               END DO
            END DO
         END DO
C     
         IF (KACTCOR.GT.0) THEN
C           
         DO I = 1,NBLEN
            TMP3(I) = VXCZ(I)*GRADAT(I,1)
            TMP4(I) = VXCZ(I)*GRADAT(I,2)
            TMP5(I) = VXCZ(I)*GRADAT(I,3)
         END DO
         DO IX = 1,3
            IP1 = IX + 1
            IXA1 = IXY2(3*IX - 2)
            IXA2 = IXY2(3*IX - 1)
            IXA3 = IXY2(3*IX)
C           ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
            DO LAC = 1,KACTCOR
               IA = IACTCOR(LAC)
               FAC = D0
               DO I = 1, NBLEN
                  GD1(I)  = VXCR(I)*GAO(I,IA,IP1)
     &                               + TMP3(I)*GAO(I,IA,IXA1)
     &                               + TMP4(I)*GAO(I,IA,IXA2)
     &                               + TMP5(I)*GAO(I,IA,IXA3)
                  GD2(I) = TMP3(I)*GAO(I,IA,IP1)
                  GD3(I) = TMP4(I)*GAO(I,IA,IP1)
                  GD4(I) = TMP5(I)*GAO(I,IA,IP1)
                  FAC = FAC - D2*(GD1(I)*GAO(I,IA,1)
     &                          + GD2(I)*GAO(I,IA,2)
     &                          + GD3(I)*GAO(I,IA,3)
     &                          + GD4(I)*GAO(I,IA,4))
               END DO
               DFTFCD(IA,IA,IX) = DFTFCD(IA,IA,IX) + FAC
               DO IBLB = 1, NBLCNT(ISYM)
                  ISTRB = NBLOCKS(1,IBLB,ISYM)
                  IENDB = NBLOCKS(2,IBLB,ISYM)
                  DO IB = ISTRB, IENDB
                     IF (IB.LT.IA) THEN
                        FAC = D0
                        DO I = 1, NBLEN
                           FAC = FAC - GD1(I) *GAO(I,IB,1)
     &                               - GD2(I)*GAO(I,IB,2)
     &                               - GD3(I)*GAO(I,IB,3)
     &                               - GD4(I)*GAO(I,IB,4)
                        END DO
                        DFTFCD(IB,IA,IX) = DFTFCD(IB,IA,IX) + FAC
                     ELSE IF (IB.GT.IA) THEN
                        FAC = D0
                        DO I = 1, NBLEN
                           FAC = FAC - GD1(I)*GAO(I,IB,1)
     &                               - GD2(I)*GAO(I,IB,2)
     &                               - GD3(I)*GAO(I,IB,3)
     &                               - GD4(I)*GAO(I,IB,4)
                        END DO
                        DFTFCD(IA,IB,IX) = DFTFCD(IA,IB,IX) + FAC
                     END IF
                  END DO
               END DO
            END DO
            DO IB = 1,NBAST
               DO IA = IB+1, NOCCT
                  DFTFCD(IA,IB,IX) = DFTFCD(IB,IA,IX)
               END DO
            END DO
         END DO
         END IF
         END DO
      ELSE
C lda
         DO ISYM = 1, NSYM
         DO IX = 1, 3
            IP1 = IX + 1
            ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
            IA = 0
            DO ISHELA = 1, KMAX
               ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + IX,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IA = IA + 1
                  IF (ACTIVE(IA)) THEN
                     IB = 0
                     DO ISHELB = 1, KMAX
                     ISCORB = IPTCNT(3*(NCENT(ISHELB) - 1) + IX,0,1)
                     DO ICOMPB = 1, KHKT(ISHELB)
                        IB = IB + 1  
                        IF (ACTIVE(IB) .AND. IB.GE.IA) THEN  
                           FAC = D0  
                           DO I = 1, NBLEN
                              FAC = FAC + TMP1(I,ISCORN)*
     &                                          GAO(I,IA,1)*GAO(I,IB,1)
                           END DO
                           DFTFDX(IA,IB,IX) = DFTFDX(IA,IB,IX) + FAC
                           IF (ISCORB.EQ.ISCORN) THEN
                              FAC = D0 
                              DO I = 1, NBLEN
                                 FAC = FAC - VXCR(I)*
     &                                       GAO(I,IA,1)*GAO(I,IB,IP1)
                              END DO
                              DFTFCD(IA,IB,IX) = DFTFCD(IA,IB,IX) + FAC
                           END IF
                           IF (ISCORA.EQ.ISCORN) THEN
                              FAC = D0
                              DO I = 1, NBLEN
                                 FAC = FAC - VXCR(I)*
     &                                       GAO(I,IA,IP1)*GAO(I,IB,1)
                              END DO 
                              DFTFCD(IA,IB,IX) = DFTFCD(IA,IB,IX) + FAC
                           END IF
                        END IF
                     END DO
                  END DO
               END IF   
            END DO         
         END DO                      
         DO IB = 1, NBAST            
            DO IA = IB+1, NBAST      
               DFTFDX(IA,IB,IX) = DFTFDX(IB,IA,IX)
            END DO      
            DO IA = IB+1, NOCCT
               DFTFCD(IA,IB,IX) = DFTFCD(IB,IA,IX)
            END DO 
         END DO 
      END DO
      END DO   
      END IF      
C              
 100  CONTINUE
      RETURN 
      END 

