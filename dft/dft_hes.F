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
C /* Deck dftmolhes */
      SUBROUTINE DFTMOLHES(LURD,WORK,LWORK,IPRINT)
C
C     This subroutine calculates the static DFT contributions
C     to the molecular hessian
C
C     O.B. Lutnaes, D.J. Wilson  Jan 2004 tuh Sep 2004
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
#include <nuclei.h>
#include <inforb.h>
      CHARACTER*4 KEY
      LOGICAL FOUND 
#include <dftcom.h>
#include "dftmolhes.h"
C
      LOGICAL  DFT_ISGGA           ! , DOKAPPA
      EXTERNAL DFT_ISGGA, DFTSTATHES
C     
      CALL QENTER('DFTMOLHES')
      DOGGA = DFT_ISGGA()
      IF (IPRINT .GE. 5) CALL TITLER('Output from DFTMOLHES','*',103)
C     
      KCMO   = 1 
      KDMAT  = KCMO   + NCMOT
      KV1MAT = KDMAT  + N2BASX
      KDFTHS = KV1MAT + 3*N2BASX*NUCDEP
      KLAST  = KDFTHS + MXCOOR*MXCOOR
C     IF (DOKAPPA) THEN 
C        KTRMAT = KDFTHS + MXCOOR*MXCOOR
C        KTKAP  = KTRMAT + 3*NUCDEP*N2BASX 
C        KDKAP  = KTKAP  + N2ORBT
C        KLAST  = KDKAP  + 3*NVARPT
C     END IF 
      LWRK   = LWORK  -  KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTMOLHES',' ',KLAST,LWORK)
C 
C     Calculate V1MAT = Sa * Dmat in AO basis
C     
      CALL GETV1MAT(0,WORK(KCMO),WORK(KDMAT),WORK(KV1MAT),
     &              WORK(KLAST),LWRK,IPRINT)
C
C Calculate kappa derivative matrix for debugging of rho derivative
C NOTE: not currently checked to work for gga
C     IF (DOKAPPA) THEN
C         WRITE(LUPRI,'(/,1X,A)') 'WARNING!! Kappa derivative matrix to', 
C    &        ' be calculated and added to HESDFT'
C         CALL GETKAPMAT(LURD,WORK(KCMO),WORK(KTRMAT),WORK(KDKAP),
C    &                   WORK(KTKAP),WORK(KLAST),LWRK,IPRINT)
C     END IF
C
      CALL DZERO(WORK(KDFTHS),MXCOOR*MXCOOR)
      CALL KICK_SLAVES_HESSTAT(NBAST,WORK(KDMAT),WORK(KV1MAT),NUCDEP,
     &                         IPRINT)
      CALL DFTINT(WORK(KDMAT),1,2,.FALSE.,WORK(KLAST),LWRK,
     &            DFTSTATHES,WORK(KDMAT),ELE)
      CALL DFT_HESSTAT_COLLECT(WORK(KDFTHS),WORK(KLAST),LWRK)
      CALL DFTHSSYM(WORK(KDFTHS),MXCOOR)
      CALL ADDHES(WORK(KDFTHS))
C
      IF (IPRINT.GT.1) THEN
        KCSTRA = KLAST
        KSCTRA = KCSTRA +  MXCOOR*MXCOOR
        KLAST  = KSCTRA +  MXCOOR*MXCOOR
        CALL HEADER('DFT static contribution to molecular Hessian',-1)
        CALL PRIHES(WORK(KDFTHS),'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
      END IF
      CALL QEXIT('DFTMOLHES')
      RETURN
      END
C
C /* Deck getv1mat */
      SUBROUTINE GETV1MAT(IATOM,CMO,DMAT,V1MAT,WORK,LWRK,IPRINT)
C
C     Constructs V1MAT = derivative of orbital connection matrix
C                        multiplied by density matrix.
C                        Output matrices are in AO-basis.
C
C     O. B. Lutnaes and D. Wilson Jan 2004
C     tuh Sep 2004 (revision)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION CMO(NCMOT),DMAT(NBAST,NBAST),QMAT(NBAST,NBAST),
     &          V1MAT(N2BASX,3,NUCDEP),TMAT(N2BASX,3),
     &          WORK(LWRK)
      PARAMETER( D0 = 0.0D0 , D1 = 1.0D0 , D2 = 2.0D0)
      LOGICAL FOUND
      INTEGER R,S,RS
      CHARACTER*4 KEY
#include <oneadr.h>
#include <inftap.h>
#include <abainf.h>
C
C     Construct DMAT and QMAT from MOs
C
      CALL DZERO(CMO,NCMOT)
      CALL RD_SIRIFC('CMO',FOUND,CMO,WORK,LWRK)
      IF(.NOT.FOUND)CALL QUIT('GETV1MAT error: CMO not found on SIRIFC')
      IF (NSYM.GT.1)CALL QUIT(
     &     'DFT analytical Hessians do not work with symmetry')
C    
      DO ISYM = 1, NSYM
         NISHI  = NISH(ISYM)
         NORBI  = NORB(ISYM)
         NBASI  = NBAS(ISYM)
         ICEND  = ICMO(ISYM)
         DO R = 1, NBASI
         DO S = 1, R
            DTRS = D0
            ICENDI = ICEND
            DO I = 1, NISHI 
               DTRS = DTRS + CMO(ICENDI+R)*CMO(ICENDI+S)
               ICENDI = ICENDI + NBASI
            END DO
            DMAT(R,S) = D2*DTRS
            DMAT(S,R) = D2*DTRS
            DO I = NISHI+1, NORBT
               DTRS = DTRS + CMO(ICENDI+R)*CMO(ICENDI+S)
               ICENDI = ICENDI + NBASI
            END DO
            QMAT(R,S) = DTRS 
            QMAT(S,R) = DTRS 
         END DO          
         END DO          
      END DO
C     
C     CONSTRUCT V1MAT = Q * Sa * Dmat
C
      IF (LWRK.LT.N2BASX) CALL STOPIT('GETV1MAT',' ',LWRK,N2BASX)
      IF (NODIFC) THEN
         KEY = 'OMAT'
         FAC = D1
         IF(IATOM.EQ.0) CALL GPOPEN(LUDA1,'ABACUS.DA1','OLD',
     &                              'DIRECT',' ',LABUFI,OLDDX)
      ELSE
         KEY = 'DMAT'
         FAC = D2
      END IF
      IF (IATOM.NE.0) THEN
         JSTRT = IATOM
         JSTOP = IATOM
      ELSE
         JSTRT = 1
         JSTOP = NUCIND
      END IF
      KATOM = 0
      DO JATOM = JSTRT, JSTOP 
         KATOM = KATOM + 1
         CALL DZERO(TMAT,3*N2BASX)
         CALL ONEDRL(KEY,TMAT(1,1),TMAT(1,2),TMAT(1,3),JATOM,
     &               .TRUE.,.TRUE.,.TRUE.,WORK,LWRK,0,0)
         DO I = 1, 3
            CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,TMAT(1,I),NBAST,
     &                 DMAT,NBAST,D0,WORK,NBAST)
            CALL DGEMM('N','N',NBAST,NBAST,NBAST,FAC,QMAT,NBAST,
     &                 WORK,NBAST,D0,V1MAT(1,I,KATOM),NBAST)
            IF (IPRINT.GT.20) THEN
               CALL HEADER('Derivative connection matrix, V1MAT',-1)
               WRITE (LUPRI,'(2X,A,2I5)') ' atom and coor:',JATOM,I
               CALL OUTPUT(V1MAT(1,I,KATOM),1,NBAST,1,NBAST,
     &                     NBAST,NBAST,1,LUPRI)
            END IF
         END DO 
      END DO
      RETURN        
      END      
C
C /* Deck dftstathes */
      SUBROUTINE DFTSTATHES(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                      COORD,WGHT,HESDAT)
C
C     Exchange-correlation contribution to molecular hessian - static
C     This subroutine splits up "cbdata" needed for the calculation
C
C     O.B.Lutnaes, D.J.Wilson, T. Helgaker  Jan 04
C
C     HESDAT includes:
C        dft molecular hessian contribution
C        DMAT (density matrix)
C        V1MAT (overlap derivative contribution matrix)
C        DKAPMT (response solution matrix - derivative kappa matrix).
C        [DKAPMT is only used for debugging]
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN), WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8), NBLOCKS(2,LDAIB,8),
     &          HESDAT(MXCOOR*MXCOOR + 3*N2BASX*NUCDEP + N2BASX)
      KDMAT  = 1
      KV1MAT = KDMAT  + N2BASX
      KDFTSH = KV1MAT + N2BASX*NUCDEP*3
      CALL DFTSTATHE1(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &            COORD,WGHT,HESDAT(KDMAT),HESDAT(KV1MAT),
     &            HESDAT(KDFTSH))
C
      RETURN
      END
C
C /* Deck dftstathe1 */
      SUBROUTINE DFTSTATHE1(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                      COORD,WGHT,DMAT,V1MAT,HESDFT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <taymol.h>
#include <nuclei.h>
#include <inforb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
      LOGICAL   ACTIVE(NBAST)
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN), WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN), NBLCNT(8), 
     &          NBLOCKS(2,LDAIB,8)
      DIMENSION GRADAT(NBLEN,3) 
      DIMENSION HESDFT(MXCOOR,MXCOOR) 
      DIMENSION DMAT(NBAST,NBAST),V1MAT(NBAST,NBAST,3,NUCDEP)
      DIMENSION RHODER(NBLEN,3*NUCDEP), ZETDER(NBLEN,3*NUCDEP),
     &          GDRHOD(NBLEN,3*NUCDEP,3)
      DIMENSION GD1(NBLEN,NBAST), GD2(NBLEN,NBAST), GD3(NBLEN,NBAST),
     &          GD4(NBLEN,NBAST)     
      DIMENSION GAGDT(NBLEN,NBAST,6), GAGDT2(NBLEN,NBAST),
     &          GAGDT3(NBLEN,NBAST) 
      DIMENSION TMP1(NBLEN,3*NUCDEP), TMP2(NBLEN,3*NUCDEP), 
     &          VXCT1(NBLEN), VXCT2(NBLEN)
      DIMENSION KACTCOR(NUCDEP), IACTCOR(NBAST,NUCDEP)
      DIMENSION IX2(6), IY2(6), IXY2(9), IXY3(18)
      DATA IX2 /1, 1, 1, 2, 2, 3/
      DATA IY2 /1, 2, 3, 2, 3, 3/
      DATA IXY2 /5, 6, 7, 6, 8, 9, 7, 9, 10/
      DATA IXY3 /11,12,13,12,14,15,13,15,16,14,17,18,15,18,19,16,19,20/
#include <dftcom.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
#include "dftmolhes.h"
      DIMENSION VXCR(NBLEN), VXCRR(NBLEN), VX(5), DENS(5),
     &          VXCZ(NBLEN), VXCRZ(NBLEN), VXCZZ(NBLEN)
C
C     Set Active orbitals 
C
      CALL DFTHESACT(0,NACTIVE,KACTCOR,IACTCOR,ACTIVE,NBLCNT,NBLOCKS,
     &               LDAIB)
      IF (NACTIVE.EQ.0) RETURN
C
      DO I = 1, NBLEN
         GRADAT(I,1) = GRADA(1,I)
         GRADAT(I,2) = GRADA(2,I)
         GRADAT(I,3) = GRADA(3,I)
      END DO
C
      CALL DFTHESRHO(KACTCOR,IACTCOR,RHODER,ZETDER,GDRHOD,
     &               GD1,GD2,GD3,GD4,
     &               NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,DMAT,V1MAT,GRADAT)
C
C     potential
C
      CALL DFTHESVX(VXCR,VXCRR,VXCZ,VXCRZ,VXCZZ,GRADAT,RHOA,WGHT,NBLEN)
C
      IF (DOGGA) THEN 
         DO I = 1, 3*NUCDEP 
         DO J = I, 3*NUCDEP
            TRM = D0 
            DO K = 1, NBLEN
               TRM = TRM + VXCRR(K)*RHODER(K,I)*RHODER(K,J)
     &               + DP5*VXCZZ(K)*ZETDER(K,I)*ZETDER(K,J)
     &               + DP5*VXCRZ(K)*(RHODER(K,I)*ZETDER(K,J)
     &                             + RHODER(K,J)*ZETDER(K,I))
     &                   + VXCZ(K)*(GDRHOD(K,I,1)*GDRHOD(K,J,1)
     &                            + GDRHOD(K,I,2)*GDRHOD(K,J,2)
     &                            + GDRHOD(K,I,3)*GDRHOD(K,J,3))
            END DO        
            HESDFT(I,J) = HESDFT(I,J) + TRM
         END DO        
         END DO           
      ELSE
         DO I = 1, 3*NUCDEP 
         DO J = I, 3*NUCDEP
            TRM = D0 
            DO K = 1, NBLEN
               TRM = TRM + VXCRR(K)*RHODER(K,I)*RHODER(K,J)
            END DO        
            HESDFT(I,J) = HESDFT(I,J) + TRM 
         END DO        
         END DO           
      END IF
C                         
C     dF/dp . dp2/(dxa.dxb)
C     dF/dz . dz2/(dxa.dxb) 
C     =====================
C
C     GGA
C
      IF (DOGGA) THEN 
         DO ISYM = 1, NSYM
C
            DO IXA = 1, 3       
               IXA1 = IXY2(3*IXA - 2)
               IXA2 = IXY2(3*IXA - 1)
               IXA3 = IXY2(3*IXA)
               DO IA = 1, NBAST
               IF (ACTIVE(IA)) THEN
                  DO I  = 1, NBLEN
                     GAGDT(I,IA,IXA) = GAO(I,IA,IXA1)*GRADAT(I,1)
     &                               + GAO(I,IA,IXA2)*GRADAT(I,2)
     &                               + GAO(I,IA,IXA3)*GRADAT(I,3)
                  END DO
               END IF
               END DO
            END DO
C
            DO IXA = 1, 3       
               IA = 0               
               DO ISHELA = 1, KMAX 
               DO ICOMPA = 1, KHKT(ISHELA)
                  ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + IXA,0,1)
                  IA = IA + 1
                  IF (ACTIVE(IA)) THEN
                     CALL DZERO(TMP1,NBLEN*3*NUCDEP)
                     CALL DZERO(TMP2,NBLEN*3*NUCDEP)
                     DO IXB = 1, 3
                        IB = 0
                        DO ISHELB = 1, KMAX 
                        DO ICOMPB = 1, KHKT(ISHELB)
                           ISCORB=IPTCNT(3*(NCENT(ISHELB)-1)+IXB,0,1)
                           IB = IB + 1
                           IF (ACTIVE(IB).AND.ISCORB.GE.ISCORA) THEN
                              DBA = DMAT(IB,IA) 
                              CALL DAXPY(NBLEN,DBA,GAO(1,IB,IXB+1),1,
     &                                              TMP1(1,ISCORB),1)
                              CALL DAXPY(NBLEN,DBA,GAGDT(1,IB,IXB),1,
     &                                              TMP2(1,ISCORB),1)
                           END IF
                        END DO
                        END DO
                     END DO
                     DO J = 1, NBLEN
                        VXCT1(J) = VXCZ(J)*GAGDT(J,IA,IXA)
     &                           + VXCR(J)*GAO(J,IA,IXA+1)
                        VXCT2(J) = VXCZ(J)*GAO(J,IA,IXA+1) 
                     END DO
                     DO I = ISCORA, 3*NUCDEP
                        TRM = D0 
                        DO J = 1, NBLEN 
                           TRM = TRM + VXCT1(J)*TMP1(J,I) 
     &                               + VXCT2(J)*TMP2(J,I)
                        END DO
                        HESDFT(ISCORA,I) = HESDFT(ISCORA,I) + D2*TRM
                     END DO
                  END IF
               END DO
            END DO
            END DO
         END DO
C
C     LDA
C
      ELSE
         DO ISYM = 1, NSYM 
            DO IXA = 1, 3
               IA = 0
               DO ISHELA = 1, KMAX
               DO ICOMPA = 1, KHKT(ISHELA)
                  ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + IXA,0,1)
                  IA = IA + 1
                  IF (ACTIVE(IA)) THEN
                     CALL DZERO(TMP1,NBLEN*3*NUCDEP)
                     DO IXB = 1, 3
                        IB = 0
                        DO ISHELB = 1, KMAX
                        DO ICOMPB = 1, KHKT(ISHELB)
                           ISCORB = IPTCNT(3*(NCENT(ISHELB)-1)+IXB,0,1)
                           IB = IB + 1 
                           IF (ACTIVE(IB).AND.ISCORB.GE.ISCORA) THEN
                              DBA = DMAT(IB,IA) 
                              CALL DAXPY(NBLEN,DBA,GAO(1,IB,IXB+1),1,
     &                                              TMP1(1,ISCORB),1)
                           END IF
                        END DO
                        END DO
                     END DO
                     DO J = 1, NBLEN
                        VXCT1(J) = VXCR(J)*GAO(J,IA,IXA+1)
                     END DO
                     DO I = ISCORA, 3*NUCDEP
                        HESDFT(ISCORA,I) = HESDFT(ISCORA,I) 
     &                         + D2*DDOT(NBLEN,VXCT1,1,TMP1(1,I),1)
                     END DO
                  END IF
               END DO
               END DO
            END DO
         END DO
      END IF
C
C Add GGA: 3rd derivative GAO contribution  
C     LDA: 2nd derivative GAO contribution
C 
C     GGA
C
      IF (DOGGA) THEN 
         DO ISYM = 1, NSYM
            DO IA = 1, NBAST
            IF (ACTIVE(IA)) THEN
               DO I  = 1, NBLEN
                  GAGDT2(I,IA) = VXCZ(I)*(GD2(I,IA)*GRADAT(I,1)
     &                                  + GD3(I,IA)*GRADAT(I,2)
     &                                  + GD4(I,IA)*GRADAT(I,3))
     &                                  + VXCR(I)*GD1(I,IA)
                  GAGDT3(I,IA) = VXCZ(I)*GD1(I,IA)
               END DO
            END IF
            END DO 
            DO IXY = 1, 6
               NXY0 = IXY + 4
               NXY1 = IXY3(3*IXY - 2)
               NXY2 = IXY3(3*IXY - 1)
               NXY3 = IXY3(3*IXY)
               IA = 0
               DO ISHELA = 1, KMAX
                  TRM1 = D0
                  DO ICOMPA = 1, KHKT(ISHELA)
                     IA = IA + 1
                     IF (ACTIVE(IA)) THEN
                        DO I  = 1, NBLEN
                           TRM1 = TRM1 + GAO(I,IA,NXY0)*GAGDT2(I,IA)
     &                                 + GAGDT3(I,IA)*
     &                                   (GAO(I,IA,NXY1)*GRADAT(I,1)
     &                                  + GAO(I,IA,NXY2)*GRADAT(I,2)
     &                                  + GAO(I,IA,NXY3)*GRADAT(I,3))
                        END DO
                     END IF
                  END DO
                  ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + IX2(IXY),0,1)
                  ISCORB = IPTCNT(3*(NCENT(ISHELA) - 1) + IY2(IXY),0,1)
                  HESDFT(ISCORA,ISCORB) = HESDFT(ISCORA,ISCORB)+D2*TRM1
               END DO
            END DO
         END DO
C
C     LDA
C
      ELSE
         DO ISYM = 1, NSYM
            DO IXY = 1, 6
               IP4 = IXY + 4
               IA = 0
               DO ISHELA = 1, KMAX
                  TRM1 = D0
                  DO ICOMPA = 1, KHKT(ISHELA)
                     IA = IA + 1
                     IF (ACTIVE(IA)) THEN
                        DO I = 1, NBLEN
                           TRM1 = TRM1 + VXCR(I)*GD1(I,IA)*GAO(I,IA,IP4)
                        END DO
                     END IF
                  END DO
                  ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + IX2(IXY),0,1)
                  ISCORB = IPTCNT(3*(NCENT(ISHELA) - 1) + IY2(IXY),0,1)
                  HESDFT(ISCORA,ISCORB) = HESDFT(ISCORA,ISCORB)+D2*TRM1
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
C
C /* Deck DFTHSSYM */
      SUBROUTINE DFTHSSYM(HESDFT,NDIM)
#include <implicit.h>
#include <priunit.h>
      DIMENSION HESDFT(NDIM,NDIM)
C
C     Symmetrize Hessian
C
      DO J = 1, NDIM
         DO I = J+1, NDIM 
            HESDFT(I,J) = HESDFT(J,I)
         END DO
      END DO
C
      RETURN
      END
C
C /* Deck dfthed */
      SUBROUTINE DFTHED(JATOM,FCD,FDFTD,WORK,LWORK,IPRINT)
C
C     add DFT component of Fock derivative matrix
C
C     O B Lutnaes, D Wilson  - Jan 04
C     tuh Sep 04
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <priunit.h>
#include <dummy.h>
#include <nuclei.h>
#include <inforb.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION FCD(N2BASX,3), FDFTD(N2BASX,3), WORK(LWORK)
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
C
      IATOM = JATOM
      IF (IPRINT .GE. 5) CALL TITLER('Output from DFTHED','*',103)
C
      KCMO   = 1
      KDMAT  = KCMO   +   NCMOT
      KV1MAT = KDMAT  +   N2BASX
      KFCD   = KV1MAT + 3*N2BASX
      KDFTFD = KFCD   + 3*N2BASX
      KLAST  = KDFTFD + 3*N2BASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTHED',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
C     Construct V1MAT = Sa * DMAT in AO basis from Sa in MO basis
C
      CALL GETV1MAT(IATOM,WORK(KCMO),WORK(KDMAT),WORK(KV1MAT),
     &              WORK(KLAST),LWRK,IPRINT)
C
      CALL DZERO(WORK(KFCD),  3*N2BASX)
      CALL DZERO(WORK(KDFTFD),3*N2BASX)
      CALL KICK_SLAVES_HED(IATOM,NBAST,WORK(KDMAT),WORK(KV1MAT),IPRINT)
      CALL DFTINT(WORK(KDMAT),1,1,.FALSE.,WORK(KLAST),LWRK,
     &            DFTKSD,WORK(KDMAT),ELE)
      CALL DFT_HED_COLLECT(WORK(KFCD),WORK(KDFTFD),N2BASX,
     &                     WORK(KLAST),LWRK)
      CALL DFTFDSYM(WORK(KFCD),WORK(KDFTFD),NBAST)
      CALL DAXPY(3*N2BASX,D1,WORK(KFCD),1,FCD,1)
      DO I = 1, 3
         CALL ONETRA(WORK(KCMO),FDFTD(1,I),WORK(KDFTFD+(I-1)*N2BASX),
     &               WORK(KLAST),LWRK,IATOM,I,IPRINT)
      END DO
C
      IF (IPRINT.GT.10) THEN
         CALL HEADER('DFT FCD Fock derivative matrix (x direction)',-1)
         CALL OUTPUT(FCD(1,1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('DFT FCD Fock derivative matrix (y direction)',-1)
         CALL OUTPUT(FCD(1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('DFT FCD Fock derivative matrix (z direction)',-1)
         CALL OUTPUT(FCD(1,3),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('DFT FD MO basis (x direction)',-1)
         CALL OUTPUT(FDFTD(1,1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)  
         CALL HEADER('DFT FD MO basis (y direction)',-1)
         CALL OUTPUT(FDFTD(1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)  
         CALL HEADER('DFT FD MO basis (z direction)',-1)
         CALL OUTPUT(FDFTD(1,3),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)  
      END IF
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
C     O. B. Lutnaes and D. Wilson Jan 2004
C     tuh Sep 2004
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
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, DP5 = 0.5D0)
      LOGICAL ACTIVE(NBAST)
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN), WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN), NBLCNT(8), 
     &          NBLOCKS(2,LDAIB,8)
      DIMENSION DMAT(NBAST,NBAST), V1MAT(NBAST,NBAST,3), GRADAT(NBLEN,3)
      DIMENSION DFTFCD(NBAST,NBAST,3), DFTFDX(NBAST,NBAST,3)
      DIMENSION RHODER(NBLEN,3*NUCDEP), ZETDER(NBLEN,3*NUCDEP),
     &          GDRHOD(NBLEN,3*NUCDEP,3)
      DIMENSION TMP1(NBLEN), TMP2(NBLEN), TMP3(NBLEN), TMP4(NBLEN)
      DIMENSION GMP1(NBLEN), GMP2(NBLEN), GMP3(NBLEN), GMP4(NBLEN)
      DIMENSION GD1(NBLEN,NBAST), GD2(NBLEN,NBAST), GD3(NBLEN,NBAST),
     &          GD4(NBLEN,NBAST)
      DIMENSION VXCR(NBLEN), VXCRR(NBLEN), 
     &          VXCZ(NBLEN), VXCRZ(NBLEN), VXCZZ(NBLEN)
      DIMENSION KACTCOR(NUCDEP), IACTCOR(NBAST,NUCDEP), IXY2(9)
      DATA IXY2 / 5, 6, 7, 6, 8, 9, 7, 9, 10/
#include <ibtfun.h>
C
C     Set Active orbitals
C
      CALL DFTHESACT(IATOM,NACTIVE,KACTCOR,IACTCOR,ACTIVE,NBLCNT,
     &               NBLOCKS,LDAIB)
      IF (NACTIVE.EQ.0) RETURN
C
      DO I = 1, 3
      DO J = 1, NBLEN
         GRADAT(J,I) = GRADA(I,J)
      END DO
      END DO     
C
      CALL DFTHESRHO(KACTCOR,IACTCOR,RHODER,ZETDER,GDRHOD,
     &               GD1,GD2,GD3,GD4,NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &               DMAT,V1MAT(1,1,1-3*(IATOM-1)),GRADAT)
C
C     potential 
C    
#ifdef MOD_XCFUN
      CALL xcfun_HESVX(VXCR,VXCRR,VXCZ,VXCRZ,VXCZZ,GRADAT,RHOA,WGHT,
     &     NBLEN,DOGGA)
#else
      CALL DFTHESVX(VXCR,VXCRR,VXCZ,VXCRZ,VXCZZ,GRADAT,RHOA,WGHT,NBLEN)
#endif
C
C     GGA
C
      IF (DOGGA) THEN
         DO ISYM = 1, NSYM
         DO IX = 1, 3
            ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
            DO I = 1, NBLEN
               TMP1(I) = RHODER(I,ISCORN)*VXCRR(I)
     &                        + ZETDER(I,ISCORN)*VXCRZ(I)*DP5
               FAC = ZETDER(I,ISCORN)*VXCZZ(I)+RHODER(I,ISCORN)*VXCRZ(I)
               TMP2(I) = FAC*GRADAT(I,1) + GDRHOD(I,ISCORN,1)*VXCZ(I)
               TMP3(I) = FAC*GRADAT(I,2) + GDRHOD(I,ISCORN,2)*VXCZ(I)
               TMP4(I) = FAC*GRADAT(I,3) + GDRHOD(I,ISCORN,3)*VXCZ(I)
            END DO            
            DO IBLA = 1, NBLCNT(ISYM) 
            DO IA = NBLOCKS(1,IBLA,ISYM), NBLOCKS(2,IBLA,ISYM)
               DO I = 1, NBLEN 
                  GMP1(I) = GAO(I,IA,1)*TMP1(I)
     &                    + GAO(I,IA,2)*TMP2(I)
     &                    + GAO(I,IA,3)*TMP3(I)
     &                    + GAO(I,IA,4)*TMP4(I)
                  GMP2(I) = GAO(I,IA,1)*TMP2(I)
                  GMP3(I) = GAO(I,IA,1)*TMP3(I)
                  GMP4(I) = GAO(I,IA,1)*TMP4(I) 
               END DO
               DO IBLB = 1, NBLCNT(ISYM)
               DO IB = NBLOCKS(1,IBLB,ISYM), NBLOCKS(2,IBLB,ISYM)
                  IF (IB.GE.IA) THEN
                     TRM = D0
                     DO I = 1, NBLEN
                        TRM = TRM + GMP1(I)*GAO(I,IB,1)
     &                            + GMP2(I)*GAO(I,IB,2)
     &                            + GMP3(I)*GAO(I,IB,3)
     &                            + GMP4(I)*GAO(I,IB,4)
                     END DO
                     DFTFDX(IA,IB,IX) = DFTFDX(IA,IB,IX) + TRM
                  END IF      
               END DO         
               END DO
            END DO         
            END DO
         END DO
C     
         DO I = 1,NBLEN
            TMP2(I) = VXCZ(I)*GRADAT(I,1)
            TMP3(I) = VXCZ(I)*GRADAT(I,2)
            TMP4(I) = VXCZ(I)*GRADAT(I,3)
         END DO
         DO IX = 1,3
            IAX0 = IX + 1
            IXA1 = IXY2(3*IX - 2)
            IXA2 = IXY2(3*IX - 1)
            IXA3 = IXY2(3*IX)
            DO LAC = 1, KACTCOR(IATOM)
               IA = IACTCOR(LAC,IATOM)
               DO I = 1, NBLEN
                  GMP1(I) = VXCR(I)*GAO(I,IA,IAX0)
     &                    + TMP2(I)*GAO(I,IA,IXA1)
     &                    + TMP3(I)*GAO(I,IA,IXA2)
     &                    + TMP4(I)*GAO(I,IA,IXA3)
                  GMP2(I) = TMP2(I)*GAO(I,IA,IAX0)
                  GMP3(I) = TMP3(I)*GAO(I,IA,IAX0)
                  GMP4(I) = TMP4(I)*GAO(I,IA,IAX0)
               END DO
               DO IBLB = 1, NBLCNT(ISYM)
               DO IB = NBLOCKS(1,IBLB,ISYM), NBLOCKS(2,IBLB,ISYM)
                  TRM = D0
                  DO I = 1, NBLEN
                     TRM = TRM - GMP1(I)*GAO(I,IB,1)
     &                         - GMP2(I)*GAO(I,IB,2)
     &                         - GMP3(I)*GAO(I,IB,3)
     &                         - GMP4(I)*GAO(I,IB,4)
                  END DO
                  DFTFCD(IA,IB,IX) = DFTFCD(IA,IB,IX) + D2*TRM
               END DO
               END DO
            END DO
         END DO
         END DO
C
C     LDA
C
      ELSE
         DO ISYM = 1, NSYM
         DO IX = 1, 3
            ISCORN = IPTCNT(3*(IATOM - 1) + IX,0,1)
            DO I = 1, NBLEN
               TMP1(I) = RHODER(I,ISCORN)*VXCRR(I)
            END DO            
            DO IBLA = 1, NBLCNT(ISYM) 
            DO IA = NBLOCKS(1,IBLA,ISYM), NBLOCKS(2,IBLA,ISYM)
               DO I = 1, NBLEN 
                  GMP1(I) = GAO(I,IA,1)*TMP1(I)
               END DO
               DO IBLB = 1, NBLCNT(ISYM)
               DO IB = NBLOCKS(1,IBLB,ISYM), NBLOCKS(2,IBLB,ISYM)
                  IF (IB.GE.IA) THEN
                     DFTFDX(IA,IB,IX) = DFTFDX(IA,IB,IX) + 
     &                    DDOT(NBLEN,GMP1,1,GAO(1,IB,1),1)
                  END IF      
               END DO         
               END DO
            END DO         
            END DO
C
            DO LAC = 1, KACTCOR(IATOM)
               IA = IACTCOR(LAC,IATOM)
               DO I = 1, NBLEN
                  GMP1(I) = VXCR(I)*GAO(I,IA,IX + 1)
               END DO
               DO IBLB = 1, NBLCNT(ISYM)
               DO IB = NBLOCKS(1,IBLB,ISYM), NBLOCKS(2,IBLB,ISYM)
                  DFTFCD(IA,IB,IX) = DFTFCD(IA,IB,IX) 
     &                 - D2*DDOT(NBLEN,GMP1,1,GAO(1,IB,1),1)
               END DO
               END DO
            END DO
         END DO
         END DO
      END IF      
      RETURN 
      END 
C
C /* Deck DFTFDSYM */
      SUBROUTINE DFTFDSYM(DFTFCD,DFTFDX,NBAST)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.5D0)
      DIMENSION DFTFCD(NBAST,NBAST,3), DFTFDX(NBAST,NBAST,3)
C
      DO IX = 1, 3         
         DO IB = 1, NBAST            
         DO IA = 1, IB - 1 
            AVE = DP5*(DFTFCD(IA,IB,IX) + DFTFCD(IB,IA,IX))
            DFTFCD(IA,IB,IX) = AVE 
            DFTFCD(IB,IA,IX) = AVE 
         END DO
         END DO
         DO IB = 1, NBAST            
         DO IA = IB + 1, NBAST      
            DFTFDX(IA,IB,IX) = DFTFDX(IB,IA,IX)
         END DO      
         END DO 
      END DO
      RETURN
      END    
C
C /* Deck dfthesact */
      SUBROUTINE DFTHESACT(KATOM,NACTIV,KACTCOR,IACTCOR,ACTIVE,
     &                     NBLCNT,NBLOCKS,LDAIB)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
      LOGICAL   ACTIVE(NBAST)
      DIMENSION NBLCNT(8), NBLOCKS(2,LDAIB,8)
      DIMENSION KACTCOR(NUCDEP), IACTCOR(NBAST,NUCDEP) 
#include <dftcom.h>
#include <symmet.h>
#include <shells.h>
#include "dftmolhes.h"
C
C     Set active orbitals 
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
C     Set up array of active orbitals centered on atoms
C
      NACTIV = 0
      CALL IZERO(KACTCOR,NUCDEP) 
      DO LATOM = 1, NUCDEP
      IF (KATOM.EQ.0 .OR. KATOM.EQ.LATOM) THEN
         IA = 0
         ISCORN = IPTCNT(3*(LATOM - 1) + 1,0,1)
         DO ISHELA = 1, KMAX
            ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + 1,0,1)
            DO ICOMPA = 1 ,KHKT(ISHELA)
               IA = IA + 1
               IF (ACTIVE(IA).AND.ISCORN.EQ.ISCORA) THEN
                  NACTIV = NACTIV + 1
                  KACTCOR(LATOM) = KACTCOR(LATOM) + 1
                  IACTCOR(KACTCOR(LATOM),LATOM) = IA
               END IF
            END DO
         END DO
      END IF
      END DO
      RETURN
      END
C
C /* Deck dfthesrho */
      SUBROUTINE DFTHESRHO(KACTCOR,IACTCOR,RHODER,ZETDER,GDRHOD,
     &                     GD1,GD2,GD3,GD4,NBLEN,NBLCNT,NBLOCKS,LDAIB,
     &                     GAO,DMAT,V1MAT,GRADAT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
      PARAMETER (D2 = 2.0D0)
      DIMENSION GAO(NBLEN,NBAST,*), NBLCNT(8), NBLOCKS(2,LDAIB,8)
      DIMENSION GRADAT(NBLEN,3)
      DIMENSION DMAT(NBAST,NBAST),V1MAT(NBAST,NBAST,3,NUCDEP)
      DIMENSION RHODER(NBLEN,3*NUCDEP), ZETDER(NBLEN,3*NUCDEP),
     &          GDRHOD(NBLEN,3*NUCDEP,3)
      DIMENSION GD1(NBLEN,NBAST), GD2(NBLEN,NBAST), GD3(NBLEN,NBAST),
     &          GD4(NBLEN,NBAST)     
      DIMENSION TMP1(NBLEN,3*NUCDEP), TMP2(NBLEN,3*NUCDEP) 
      DIMENSION IXY2(9), IDONEB(3*NUCDEP) 
      DIMENSION KACTCOR(NUCDEP), IACTCOR(NBAST,NUCDEP) 
      DATA IXY2 /5, 6, 7, 6, 8, 9, 7, 9, 10/
#include <dftcom.h>
#include <symmet.h>
#include <shells.h>
#include "dftmolhes.h"
C
C     GD = DMAT*GAO  
C     =============
C
C     GGA
C
      IF (DOGGA) THEN 
         DO ISYM = 1, NSYM
            CALL DZERO(GD1,NBLEN*NBAST)
            CALL DZERO(GD2,NBLEN*NBAST)
            CALL DZERO(GD3,NBLEN*NBAST)
            CALL DZERO(GD4,NBLEN*NBAST)
            DO IBLA = 1, NBLCNT(ISYM)
               ISTRA = NBLOCKS(1,IBLA,ISYM)
               IENDA = NBLOCKS(2,IBLA,ISYM)
               DO IA = ISTRA, IENDA
                  DO IBLB = 1, NBLCNT(ISYM)
                     ISTRB = NBLOCKS(1,IBLB,ISYM)
                     IENDB = NBLOCKS(2,IBLB,ISYM)
                     DO IB = ISTRB, IENDB
                        DBA = DMAT(IB,IA)
                        CALL DAXPY(NBLEN,DBA,GAO(1,IB,1),1,GD1(1,IA),1)
                        CALL DAXPY(NBLEN,DBA,GAO(1,IB,2),1,GD2(1,IA),1)
                        CALL DAXPY(NBLEN,DBA,GAO(1,IB,3),1,GD3(1,IA),1)
                        CALL DAXPY(NBLEN,DBA,GAO(1,IB,4),1,GD4(1,IA),1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C     LDA
C
      ELSE 
         DO ISYM = 1, NSYM
            CALL DZERO(GD1,NBLEN*NBAST)
            DO IBLA = 1, NBLCNT(ISYM)
               ISTRA = NBLOCKS(1,IBLA,ISYM)
               IENDA = NBLOCKS(2,IBLA,ISYM)
               DO IA = ISTRA, IENDA
                  DO IBLB = 1, NBLCNT(ISYM)
                     ISTRB = NBLOCKS(1,IBLB,ISYM)
                     IENDB = NBLOCKS(2,IBLB,ISYM)
                     DO IB = ISTRB, IENDB
                        CALL DAXPY(NBLEN,DMAT(IB,IA),GAO(1,IB,1),1,
     &                                               GD1(1,IA),1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C     Construct rho and zeta derivatives
C     Also construct grad*(rho derivative)
C     ====================================
C
      CALL DZERO(RHODER,3*NUCDEP*NBLEN)
C
C     GGA
C
      IF (DOGGA) THEN 
         CALL DZERO(GDRHOD,9*NBLEN*NUCDEP)
         DO ISYM = 1, NSYM
            DO LATOM = 1, NUCDEP
            IF (KACTCOR(LATOM).GT.0) THEN 
               DO IXA = 1, 3
                  ISCORN = IPTCNT(3*(LATOM - 1) + IXA,0,1)
                  IXA1 = IXY2(3*IXA - 2)
                  IXA2 = IXY2(3*IXA - 1)
                  IXA3 = IXY2(3*IXA)
                  DO J = 1, KACTCOR(LATOM)
                     IA = IACTCOR(J,LATOM)
                     DO I = 1, NBLEN 
                       GA = GAO(I,IA,IXA+1)
                       RHODER(I,ISCORN)   = RHODER(I,ISCORN) 
     &                                    - GD1(I,IA)*GA
                       GDRHOD(I,ISCORN,1) = GDRHOD(I,ISCORN,1)
     &                                    - GD1(I,IA)*GAO(I,IA,IXA1) 
     &                                    - GD2(I,IA)*GA
                       GDRHOD(I,ISCORN,2) = GDRHOD(I,ISCORN,2)
     &                                    - GD1(I,IA)*GAO(I,IA,IXA2) 
     &                                    - GD3(I,IA)*GA
                       GDRHOD(I,ISCORN,3) = GDRHOD(I,ISCORN,3)
     &                                    - GD1(I,IA)*GAO(I,IA,IXA3) 
     &                                    - GD4(I,IA)*GA
                     END DO
                  END DO
               END DO 
            END IF
            END DO
         END DO
C
C     LDA
C
      ELSE
         DO ISYM = 1, NSYM
            DO LATOM = 1, NUCDEP
            IF (KACTCOR(LATOM).GT.0) THEN 
               DO IXA = 1, 3
                  ISCORN = IPTCNT(3*(LATOM - 1) + IXA,0,1)
                  DO J = 1, KACTCOR(LATOM)
                     IA = IACTCOR(J,LATOM)
                     DO I = 1,NBLEN
                        RHODER(I,ISCORN) = RHODER(I,ISCORN)
     &                                   - GD1(I,IA)*GAO(I,IA,IXA+1)
                     END DO
                  END DO
               END DO
            END IF
            END DO
         END DO
      END IF
C
C     Add orbital connection matrix contribution 
C     ==========================================
C
      DO LATOM = 1, NUCDEP
      DO IX = 1, 3
         ISCORN = IPTCNT(3*(LATOM - 1) + IX,0,1) 
         CALL DSCAL(NBLEN,D2,RHODER(1,ISCORN),1)
         IF (DOGGA) THEN
            CALL DSCAL(NBLEN,D2,GDRHOD(1,ISCORN,1),1)
            CALL DSCAL(NBLEN,D2,GDRHOD(1,ISCORN,2),1)
            CALL DSCAL(NBLEN,D2,GDRHOD(1,ISCORN,3),1)
         END IF
      END DO          
      END DO          
C
C     GGA
C
      IF (DOGGA) THEN
         DO ISYM = 1, NSYM
         DO LATOM = 1, NUCDEP
            IF (KACTCOR(LATOM).GT.0) THEN 
               DO IX = 1, 3
                  ISCORN = IPTCNT(3*(LATOM - 1) + IX,0,1)
                  DO IBLA = 1, NBLCNT(ISYM)
                     ISTRA = NBLOCKS(1,IBLA,ISYM)
                     IENDA = NBLOCKS(2,IBLA,ISYM)
                     DO IA = ISTRA, IENDA
                        CALL IZERO(IDONEB,3*NUCDEP)
                        DO IBLB = 1, NBLCNT(ISYM)
                           ISTRB = NBLOCKS(1,IBLB,ISYM)
                           IENDB = NBLOCKS(2,IBLB,ISYM)
                           DO IB = ISTRB, IENDB
                              VBA = - V1MAT(IB,IA,IX,LATOM)
                              VAB = - V1MAT(IA,IB,IX,LATOM)
                              IF (IDONEB(ISCORN).EQ.0) THEN
                                 IDONEB(ISCORN) = 1
                                 DO I = 1, NBLEN
                                    TMP1(I,ISCORN) = VBA*GAO(I,IB,1)
                                    TMP2(I,ISCORN) = VAB*GAO(I,IB,1)
                                 END DO
                              ELSE
                                 CALL DAXPY(NBLEN,VBA,GAO(1,IB,1),1,
     &                                                TMP1(1,ISCORN),1)
                                 CALL DAXPY(NBLEN,VAB,GAO(1,IB,1),1,
     &                                                TMP2(1,ISCORN),1)
                              END IF
                           END DO
                        END DO
                        DO I = 1, NBLEN 
                           TMP2(I,ISCORN)=TMP2(I,ISCORN)+TMP1(I,ISCORN)
                           RHODER(I,ISCORN)   = RHODER(I,ISCORN)
     &                                 + TMP1(I,ISCORN)*GAO(I,IA,1)
                           GDRHOD(I,ISCORN,1) = GDRHOD(I,ISCORN,1)
     &                                 + TMP2(I,ISCORN)*GAO(I,IA,2)
                           GDRHOD(I,ISCORN,2) = GDRHOD(I,ISCORN,2)
     &                                 + TMP2(I,ISCORN)*GAO(I,IA,3)
                           GDRHOD(I,ISCORN,3) = GDRHOD(I,ISCORN,3)
     &                                 + TMP2(I,ISCORN)*GAO(I,IA,4)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END DO
         END DO
         DO I = 1, 3*NUCDEP
            DO K = 1, NBLEN
               ZETDER(K,I) = D2*(GDRHOD(K,I,1)*GRADAT(K,1)
     &                         + GDRHOD(K,I,2)*GRADAT(K,2)
     &                         + GDRHOD(K,I,3)*GRADAT(K,3))
            END DO
         END DO
C
C     LDA 
C
      ELSE
         DO ISYM = 1, NSYM
         DO LATOM = 1 , NUCDEP
            IF (KACTCOR(LATOM).GT.0) THEN  
               DO IX = 1, 3
                  ISCORN = IPTCNT(3*(LATOM - 1) + IX,0,1)
                  DO IBLA = 1, NBLCNT(ISYM)
                     ISTRA = NBLOCKS(1,IBLA,ISYM)
                     IENDA = NBLOCKS(2,IBLA,ISYM)
                     DO IA = ISTRA, IENDA
                        CALL IZERO(IDONEB,3*NUCDEP)
                        DO IBLB = 1, NBLCNT(ISYM)
                           ISTRB = NBLOCKS(1,IBLB,ISYM)
                           IENDB = NBLOCKS(2,IBLB,ISYM)
                           DO IB = ISTRB, IENDB
                              VBA = V1MAT(IB,IA,IX,LATOM)
                              IF (IDONEB(ISCORN).EQ.0) THEN 
                                 IDONEB(ISCORN) = 1
                                 DO I = 1, NBLEN 
                                    TMP1(I,ISCORN) = VBA*GAO(I,IB,1)
                                 END DO 
                              ELSE
                                 CALL DAXPY(NBLEN,VBA,GAO(1,IB,1),1,
     &                                                TMP1(1,ISCORN),1)
                              END IF
                           END DO
                        END DO
                        DO I = 1, NBLEN
                           RHODER(I,ISCORN) = RHODER(I,ISCORN)
     &                                      - TMP1(I,ISCORN)*GAO(I,IA,1)
                        END DO
                     END DO
                  END DO
               END DO
            END IF 
         END DO
         END DO
      END IF
      RETURN
      END
C
C /* Deck dfthesvx */
      SUBROUTINE DFTHESVX(VXCR,VXCRR,VXCZ,VXCRZ,VXCZZ,GRADAT,RHOA,WGHT,
     &                    NBLEN)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <taymol.h>
#include <nuclei.h>
#include <inforb.h>
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      DIMENSION WGHT(NBLEN), RHOA(NBLEN), GRADAT(NBLEN,3) 
      DIMENSION VXCR(NBLEN), VXCRR(NBLEN), 
     &          VXCZ(NBLEN), VXCRZ(NBLEN), VXCZZ(NBLEN)
      DIMENSION VX(9), DENS(5) 
#include <dftcom.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
#include "dftmolhes.h"
      CALL DZERO(DENS,5)
      CALL DZERO(VX,9)
      IF (DOGGA) THEN
         DO I = 1, NBLEN
            GRDNRM =  SQRT(GRADAT(I,1)**2+GRADAT(I,2)**2+GRADAT(I,3)**2)
            IF(GRDNRM.LE.1D-40) GDRNRM = 1D-40
            DENS(1) = DP5*RHOA(I)
            DENS(2) = DP5*RHOA(I)
            DENS(3) = DP5*GRDNRM
            DENS(4) = DP5*GRDNRM
            DENS(5) = DENS(3)*DENS(4)
            CALL DFTPOT1(VX,WGHT(I),DENS,.FALSE.)
            GDNRMI   = D1/GRDNRM
            VXCR(I)  = VX(1)
            VXCZ(I)  = VX(2)*GDNRMI + VX(9)
            VXCRR(I) = VX(3)
            VXCRZ(I) = VX(4)*GDNRMI + VX(6)
            VXCZZ(I) = DP5*(VX(5)*(GDNRMI**2)-VX(2)*(GDNRMI**3))
     &                 + DP5*VX(8) + GDNRMI* VX(7)
         END DO
      ELSE
         DO I = 1, NBLEN
            DENS(1) = DP5*RHOA(I)
            DENS(2) = DP5*RHOA(I)
            CALL DFTPOT1(VX,WGHT(I),DENS,.FALSE.)
            VXCR(I)  = VX(1)
            VXCRR(I) = VX(3)
         END DO
      END IF
      RETURN
      END
C    
C /* Deck getkapmat */
      SUBROUTINE GETKAPMAT(LURD,CMO,TRMAT,DKAPPA,TOTKAP,
     &                       WORK,LWRK,IPRINT)
C           
C     Construct KAPMAT (derivative of orbital variational param matrix) 
C     Output matrices are in AO-basis.
C     Only needed for debugging purposes.
C     
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
      DIMENSION TRMAT(NBAST,NBAST,NUCDEP*3), TOTKAP(NORBT,NORBT),
     &          CMO(NCMOT), DKAPPA(NVARPT), WORK(LWRK)
#include <inftap.h> 
#include <abainf.h> 
#include <gdvec.h>
#include <infvar.h> 
#include <iratdef.h>
#include <inflin.h>
C
      IF (IPRINT .GE. 5) CALL TITLER('Output from GETKAP','*',103)
C    
C  Read derivative kappa matrix and transform to AO-basis
C
      CALL DZERO(TRMAT,N2BASX*3*NUCDEP)
      DO IOP = 1, NGDVEC(1,1)
         CALL DZERO(TOTKAP,N2ORBT)
         CALL DZERO(DKAPPA,NVARPT)
         IREC = 2*IGDREC(IOP,LSYMPT,1)-1
         CALL READDX(LURD,IREC,IRAT*NVARPT,DKAPPA)
C
C  Construct matrix kappa^a(NORBT,NORBT)
C           
         DO I = 1, NWOPT
            K = JWOP(1,I)
            L = JWOP(2,I)
            TOTKAP(K,L) = - DKAPPA(I)
            TOTKAP(L,K) = + DKAPPA(I)
         END DO
         IF (IPRINT.GE.5) THEN
            CALL HEADER('Total kappa matrix (SO basis)',-1)
            CALL OUTPUT(TOTKAP(1,1),1,NORBT,1,NORBT,NORBT,
     &        NORBT,1,LUPRI)
         END IF
C
C  Transforming kappa derivative matrix to AO-basis
C
         CALL TR1DEN(CMO,TOTKAP,DUMMY,TRMAT(1,1,IOP),DUMMY,
     &               WORK,LWRK)
         IF (IPRINT.GE.5) THEN
            CALL HEADER('Total kappa matrix (AO basis)',-1)
            CALL OUTPUT(TRMAT(1,1,IOP),1,NBAST,1,NBAST,NBAST,
     &                  NBAST,1,LUPRI)
         END IF
      END DO
C
      RETURN
      END
C
      SUBROUTINE KICK_SLAVES_HESSTAT(NBAST,DMAT,V1MAT,NUCDEP,IPRINT)
#include "implicit.h"
      DIMENSION DMAT(NBAST,NBAST), V1MAT(3*NBAST*NBAST*NUCDEP)
C
#if defined (VAR_MPI)
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
C
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = 8
         CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL KSMSYNC(DMAT)
         ISTART = 1
         DO I = 1, NUCDEP
            CALL MPI_BCAST(V1MAT(ISTART),3*NBAST*NBAST,
     &                     MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                     IERR)
            ISTART = ISTART + 3*NBAST*NBAST
         END DO
      END IF
      RETURN
#endif
      END
C
#if defined (VAR_MPI)
      SUBROUTINE DFT_HESSTAT_SLAVE(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "dftmolhes.h"
#include "nuclei.h"
C
      DIMENSION WORK(LWORK)
#include "dftcom.h"
#include "mpif.h"
      LOGICAL  DFT_ISGGA
      EXTERNAL DFT_ISGGA, DFTSTATHES
C     
      KCMO   = 1 
      KDMAT  = KCMO   + NCMOT
      KV1MAT = KDMAT  + N2BASX
      KDFTHS = KV1MAT + 3*N2BASX*NUCDEP
      KLAST  = KDFTHS + MXCOOR*MXCOOR
      LWRK   = LWORK  -  KLAST + 1
      IF (KLAST .GT. LWORK) CALL 
     &     STOPIT('DFT_HESSTAT_SLAVE',' ',KLAST,LWORK)
C
      CALL DFTINTBCAST
      CALL KSMSYNC(WORK(KDMAT))
      ISTART = KV1MAT
      DO I = 1, NUCDEP
         CALL MPI_BCAST(WORK(ISTART),3*N2BASX,
     &                  Mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
         ISTART = ISTART + 3*N2BASX
      END DO
      CALL DZERO(WORK(KDFTHS),MXCOOR*MXCOOR)
      CALL DFTINT(WORK(KDMAT),1,2,.FALSE.,WORK(KLAST),LWRK,
     &            DFTSTATHES,WORK(KDMAT),ELE)
      CALL DFT_HESSTAT_COLLECT(WORK(KDFTHS),WORK(KLAST),LWRK)
      RETURN
      END
#endif
C
      SUBROUTINE DFT_HESSTAT_COLLECT(DFTHES,WORK,LWORK)
#include "implicit.h"
#include "mxcent.h"
      DIMENSION DFTHES(MXCOOR,MXCOOR), WORK(LWORK)
C
#if defined (VAR_MPI)
#include "mpif.h"
C
      CALL DCOPY(MXCOOR*MXCOOR,DFTHES,1,WORK,1)
      CALL MPI_REDUCE(WORK,DFTHES,MXCOOR*MXCOOR,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
#endif
      RETURN
      END
C
      SUBROUTINE KICK_SLAVES_HED(IATOM,NBAST,DMAT,V1MAT,IPRINT)
#include "implicit.h"
      DIMENSION DMAT(NBAST,NBAST), V1MAT(3*NBAST*NBAST)
C
#if defined (VAR_MPI)
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
C
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = 9
         CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IATOM,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL KSMSYNC(DMAT)
         CALL MPI_BCAST(V1MAT,3*NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &                  MASTER,MPI_COMM_WORLD,IERR)
      END IF
      RETURN
#endif
      END
#if defined (VAR_MPI)
C /* Deck dfthed */
      SUBROUTINE DFTHED_SLAVE(WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <priunit.h>
#include <dummy.h>
#include <nuclei.h>
#include <inforb.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION FCD(N2BASX,3), FDFTD(N2BASX,3), WORK(LWORK)
#include <symmet.h>
#include <symind.h>
#include <dftcom.h>
#include "infpar.h"
#include "dftmolhes.h"
#include <mpif.h>
C
      EXTERNAL DFTKSD
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      KDMAT  = 1
      KV1MAT = KDMAT  +   N2BASX
      KFCD   = KV1MAT + 3*N2BASX
      KDFTFD = KFCD   + 3*N2BASX
      KLAST  = KDFTFD + 3*N2BASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTHED_SLAVES',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      DOGGA = DFT_ISGGA()
C
C     Construct V1MAT = Sa * DMAT in AO basis from Sa in MO basis
C
      CALL MPI_BCAST(IATOM,1,my_MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL DFTINTBCAST
      CALL KSMSYNC(WORK(KDMAT))
      CALL MPI_BCAST(WORK(KV1MAT),3*N2BASX,MPI_DOUBLE_PRECISION,
     &               0,MPI_COMM_WORLD,IERR)
C
      CALL DZERO(WORK(KFCD),  3*N2BASX)
      CALL DZERO(WORK(KDFTFD),3*N2BASX)
      CALL DFTINT(WORK(KDMAT),1,1,.FALSE.,WORK(KLAST),LWRK,
     &            DFTKSD,WORK(KDMAT),ELE)
      CALL DFT_HED_COLLECT(WORK(KFCD),WORK(KDFTFD),N2BASX,
     &                     WORK(KLAST),LWRK)
      RETURN
      END
#endif
C
      SUBROUTINE DFT_HED_COLLECT(FCD,DFTFD,N2BASX,WORK,LWORK)
#include "implicit.h"
      DIMENSION FCD(3*N2BASX), DFTFD(3*N2BASX), WORK(LWORK)
C
#if defined (VAR_MPI)
#include "mxcent.h"
#include "maxorb.h"
#include "mpif.h"
C
      CALL DCOPY(3*N2BASX,FCD,1,WORK,1)
      CALL MPI_REDUCE(WORK,FCD,3*N2BASX,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL DCOPY(3*N2BASX,DFTFD,1,WORK,1)
      CALL MPI_REDUCE(WORK,DFTFD,3*N2BASX,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      RETURN
#endif
      END
