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
C /* Deck dftmolhes */
      SUBROUTINE DFTMOLHES(LURD,WORK,LWORK,IPRINT)
C
C     This subroutine calculates the static DFT contributions
C     to the molecular hessian
C
C     O.B. Lutnaes, D.J. Wilson  Jan 2004
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
C#include <ccom.h>
C#include <symmet.h>
C#include <symind.h>
C#include <dftcom.h>
C#include <cbires.h>
C#include <dummy.h>
C#include <abainf.h>
#include "dftmolhes.h"
C#include <inflin.h>
C
      EXTERNAL DFTSTATHES
      LOGICAL DFT_ISGGA, DOKAPPA
      EXTERNAL DFT_ISGGA
C     
      DOGGA = DFT_ISGGA()
      IF (IPRINT .GE. 5) CALL TITLER('Output from DFTMOLHES','*',103)
C     
      KCMO   = 1 
      KQMAT  = KCMO   +  NCMOT
      KTMAT  = KQMAT  +  N2BASX
      KDMAT  = KTMAT  +  3*N2BASX
      KV1MAT = KDMAT  +  N2BASX
      KDFTHS = KV1MAT +  3*N2BASX*NUCDEP
      KLAST  = KDFTHS +  MXCOOR*MXCOOR
C     IF (DOKAPPA) THEN 
C        KTRMAT = KDFTHS +  MXCOOR*MXCOOR
C        KTKAP  = KTRMAT +  3*NUCDEP*N2BASX 
C        KDKAP  = KTKAP  +  N2ORBT
C        KLAST  = KDKAP  +  3*NVARPT
C     END IF 
      LWRK   = LWORK  -  KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTMOLHES',' ',KLAST,LWORK)
C 
C calculate V1MAT = Sa * Dmat in AO basis
C     
      CALL GETV1MAT(IDUMMY,.FALSE.,WORK(KCMO),WORK(KDMAT),WORK(KQMAT),
     &              WORK(KV1MAT),WORK(KTMAT),WORK(KLAST),LWRK,IPRINT)
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
      CALL DFTINT(WORK(KDMAT),1,3,.FALSE.,WORK(KLAST),LWRK,
     &            DFTSTATHES,WORK(KDMAT),ELE)
      CALL ADDHES(WORK(KDFTHS))
C
      IF (IPRINT.GT.1) THEN
        KCSTRA = KLAST
        KSCTRA = KCSTRA +  MXCOOR*MXCOOR
        KLAST  = KSCTRA +  MXCOOR*MXCOOR
        CALL HEADER('DFT static contribution to molecular Hessian',-1)
        CALL PRIHES(WORK(KDFTHS),'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
      END IF
      RETURN
      END
C
C /* Deck getv1mat */
      SUBROUTINE GETV1MAT(IATOM,DOONE,CMO,DMAT,QMAT,V1MAT,TMAT,WORK,
     &                LWRK,IPRINT)
C
C     Constructs V1MAT = derivative of orbital connection matrix
C                        multiplied by density matrix.
C                        Output matrices are in AO-basis.
C
C     O. B. Lutn?s and D. Wilson Jan 2004
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION CMO(NCMOT),DMAT(NBAST,NBAST),QMAT(NBAST,NBAST),
     &          V1MAT(NBAST,NBAST,3,NUCDEP),TMAT(NBAST,NBAST,3),
     &          WORK(LWRK)
      PARAMETER( D0 = 0.0D0 , D1 = 1.0D0 , D2 = 2.0D0)
      LOGICAL FOUND, DOONE
      INTEGER R,S,RS
      CHARACTER*4 KEY
#include <oneadr.h>
#include <inftap.h>
#include <abainf.h>
C
C  Construct DMAT and QMAT from MOs
C
      CALL DZERO(CMO,NCMOT)
      CALL DZERO(DMAT,N2BASX)
      CALL DZERO(QMAT,N2BASX)
      NORBI = 0
      NISHI = 0
      NBASI = 0
C  
      CALL RDSIFC('CMO',FOUND,CMO,WORK,LWRK)
      IF(.NOT.FOUND) CALL QUIT('GETV1MAT error: CMO not
     &                                   found on SIRIFC')
C    
      DO ISYM = 1, NSYM
         NORBI  = NORB(ISYM)
         NISHI  = NISH(ISYM)
         NBASI  = NBAS(ISYM)
         RS     = 0
         ICEND  = 0
         ICENDI = 0
         DO R = 1, NBASI
         DO S = 1, NBASI
            RS = RS + 1
            DTRS = D0
            ICENDI = ICEND
            DO I = 1, NISHI 
               DTRS = DTRS + CMO(ICENDI+R)*CMO(ICENDI+S)
               ICENDI = ICENDI + NBASI
            END DO
            DMAT(R,S) = D2*DTRS
            DTUV = D0
            ICENDI = ICEND
            DO I = 1, NORBT
               DTUV = DTUV + CMO(ICENDI+R)*CMO(ICENDI+S)
               ICENDI = ICENDI + NBASI
            END DO
            QMAT(R,S) = DTUV 
         END DO          
         END DO          
      END DO
C     
C*****************************
C  CONSTRUCT V1MAT = Q * Sa * Dmat
C*****************************
C     DOONE = TRUE: Constructs V1MAT for only one atom (IATOM)
C             Needed for subroutine DFTHED
C     DOONE = FALSE: Constructs V1MAT for all atoms
C             Needed for subroutine DFTMOLHES
C     KEY: 'OMAT' with .NODIFC, otherwise 'DMAT'
C
      IF (DOONE) THEN
         IF (IPRINT.GT.50) WRITE(LUPRI,'(/,A,I5)') ' V1MAT for atom ',
     &                                              IATOM
         CALL DZERO(V1MAT(1,1,1,1),3*N2BASX)
         CALL DZERO(TMAT,3*N2BASX)
         IF (NODIFC) THEN
            KEY = 'OMAT'
c           CALL GPOPEN(LUDA1,'ABACUS.DA1','OLD','DIRECT',' ',
c    &           LABUFI,OLDDX)
            CALL ONEDRL(KEY,TMAT(1,1,1),TMAT(1,1,2),TMAT(1,1,3),IATOM,
     &           .TRUE.,.TRUE.,.TRUE.,WORK,LWRK,0,0)
            IF (LWRK .LT. N2BASX) CALL STOPIT('DFTMHES',' ',LWRK,N2BASX)
            DO I = 1, 3
               CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,TMAT(1,1,I),
     &              NBAST,DMAT,NBAST,D0,WORK,NBAST)
               CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,QMAT,NBAST,
     &              WORK,NBAST,D0,V1MAT(1,1,I,1),NBAST)
            END DO
         ELSE
            KEY = 'DMAT'
            CALL ONEDRL(KEY,TMAT(1,1,1),TMAT(1,1,2),TMAT(1,1,3),IATOM,
     &           .TRUE.,.TRUE.,.TRUE.,WORK,LWRK,0,0)
            IF (LWRK .LT. N2BASX) CALL STOPIT('DFTMHES',' ',LWRK,N2BASX)
            DO I = 1, 3
               CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,TMAT(1,1,I),
     &              NBAST,DMAT,NBAST,D0,WORK,NBAST)
               CALL DGEMM('N','N',NBAST,NBAST,NBAST,D2,QMAT,NBAST,
     &              WORK,NBAST,D0,V1MAT(1,1,I,1),NBAST)
            END DO
         END IF
         IF (IPRINT.GT.20) THEN
            DO I = 1, 3
               CALL HEADER('Derivative connection matrix, V1MAT',-1)
               CALL OUTPUT(V1MAT(1,1,I,1),1,NBAST,1,NBAST,NBAST,
     &                     NBAST,1,LUPRI)
            END DO
         END IF
      ELSE
         CALL DZERO(V1MAT,3*N2BASX*NUCDEP)
         IF (NODIFC) THEN
            KEY = 'OMAT'
            CALL GPOPEN(LUDA1,'ABACUS.DA1','OLD','DIRECT',' ',
     &           LABUFI,OLDDX)
            IF (LWRK .LT. N2BASX) CALL STOPIT('GETV1MAT',' ',LWRK,
     &           N2BASX)
            DO JATOM = 1, NUCIND
               IF(IPRINT.GT.100)WRITE(LUPRI,'(/,A,I5)') 'V1MAT for atom'
     &              ,JATOM
               CALL DZERO(TMAT,3*N2BASX)
               CALL ONEDRL(KEY,TMAT(1,1,1),TMAT(1,1,2),TMAT(1,1,3),JATOM
     &              ,.TRUE.,.TRUE.,.TRUE.,WORK,LWRK,0,0)
               DO I = 1, 3
                  CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,TMAT(1,1,I),
     &                 NBAST,DMAT,NBAST,D0,WORK,NBAST)
                  CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,QMAT,NBAST,
     &                 WORK,NBAST,D0,V1MAT(1,1,I,JATOM),NBAST)
               END DO
            END DO
         ELSE
            KEY = 'DMAT'
            DO JATOM = 1, NUCDEP
               CALL DZERO(TMAT,3*N2BASX)            
               CALL ONEDRL(KEY,TMAT(1,1,1),TMAT(1,1,2),TMAT(1,1,3),JATOM
     &              ,.TRUE.,.TRUE.,.TRUE.,WORK,LWRK,0,0)
               DO I = 1, 3
                  CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,TMAT(1,1,I),
     &                 NBAST,DMAT,NBAST,D0,WORK,NBAST)
                  CALL DGEMM('N','N',NBAST,NBAST,NBAST,D2,QMAT,NBAST,
     &                 WORK,NBAST,D0,V1MAT(1,1,I,JATOM),NBAST)
               END DO 
            END DO
         END IF  
         IF (IPRINT.GT.20) THEN
            DO J = 1, NUCDEP
            DO I = 1, 3 
               CALL HEADER('Derivative connection matrix, V1MAT',-1)
               CALL OUTPUT(V1MAT(1,1,I,J),1,NBAST,1,NBAST,
     &                     NBAST,NBAST,1,LUPRI)
            END DO
            END DO
         END IF
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
C /* Deck dftstathes */
      SUBROUTINE DFTSTATHES(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                       COORD,WGHT,HESDAT)
C
C     Exchange-correlation contribution to molecular hessian - static
C     This subroutine splits up "cbdata" needed for the calculation
C
C     O.B.Lutn?s, D.J.Wilson, T. Helgaker  Jan 04
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
     &                       COORD,WGHT,DMAT,V1MAT,HESDFT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <taymol.h>
#include <nuclei.h>
#include <inforb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D2 = 2.0D0, D4 = 4.0D0)
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
      DIMENSION GAGDT(NBLEN,NBAST,6), GAGDT2(NBLEN,NBAST) 
      DIMENSION TMP1(NBLEN,3*NUCDEP), TMP2(NBLEN,3*NUCDEP), 
     &          VXCT1(NBLEN), VXCT2(NBLEN)
      DIMENSION IX2(6), IY2(6), IXY2(9), IXY3(18), IDONEB(3*NUCDEP) 
      DIMENSION KACTCOR(NUCDEP), IACTCOR(NBAST,NUCDEP) 
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
#include <ibtfun.h>
C
C Set Active orbitals 
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
C Setting up array of active orbitals centered on atoms
C
      CALL IZERO(KACTCOR,NUCDEP) 
      DO LATOM = 1, NUCDEP
         IA = 0
         ISCORN = IPTCNT(3*(LATOM - 1) + 1,0,1)
         DO ISHELA = 1, KMAX
            ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + 1,0,1)
            DO ICOMPA = 1 ,KHKT(ISHELA)
               IA = IA + 1
               IF (ACTIVE(IA).AND.ISCORN.EQ.ISCORA) THEN
                  KACTCOR(LATOM) = KACTCOR(LATOM) + 1
                  IACTCOR(KACTCOR(LATOM),LATOM) = IA
               END IF
            END DO
         END DO
      END DO
C
C  Construct transpose of GRADA 
C
      DO I = 1, NBLEN
         GRADAT(I,1) = GRADA(1,I)
         GRADAT(I,2) = GRADA(2,I)
         GRADAT(I,3) = GRADA(3,I)
      END DO
C
C Kappa only used for debugging 
C     DOKAPPA = .FALSE.
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
            CALL DFTPOT1(VX,WGHT(I),DENS(1),.FALSE.)
            VXCR(I)  = D2*VX(1)
            VXCZ(I)  = VX(2)/GRDNRM
            VXCRR(I) = VX(3)
            VXCRZ(I) = DP5*VX(4)/GRDNRM
            VXCZZ(I) = DP5*DP5*(VX(5)/GRDNRM**2 - VX(2)/GRDNRM**3)
         END DO
      ELSE
         DO I = 1, NBLEN
            DENS(1) = DP5*RHOA(I)
            DENS(2) = DP5*RHOA(I)
            CALL DFTPOT1(VX,WGHT(I),DENS(1),.FALSE.)
            VXCRR(I) = VX(3)
            VXCR(I) = D2*VX(1)
         END DO
      END IF
C
C Setup temp array GD = DMAT*GAO  
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
     &                                                  GD1(1,IA),1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C Construct rho and zeta derivatives
C Also construct grad*(rho derivative)
C
      CALL DZERO(RHODER,3*NUCDEP*NBLEN)
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
     &                                       - GD1(I,IA)*GA
                          GDRHOD(I,ISCORN,1) = GDRHOD(I,ISCORN,1)
     &                                       - GD1(I,IA)*GAO(I,IA,IXA1) 
     &                                       - GD2(I,IA)*GA
                          GDRHOD(I,ISCORN,2) = GDRHOD(I,ISCORN,2)
     &                                       - GD1(I,IA)*GAO(I,IA,IXA2) 
     &                                       - GD3(I,IA)*GA
                          GDRHOD(I,ISCORN,3) = GDRHOD(I,ISCORN,3)
     &                                       - GD1(I,IA)*GAO(I,IA,IXA3) 
     &                                       - GD4(I,IA)*GA
                        END DO
                     END DO
                  END DO 
               END IF
            END DO
         END DO
C lda
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
     &                                      - GD1(I,IA)*GAO(I,IA,IXA+1)
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
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
C Add orbital connection matrix contribution 
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
      ELSE
C lda
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
C 
C Nuclear derivative of rho, gradrho, and zeta now finished
C           
      IF (DOGGA) THEN 
         DO I = 1, 3*NUCDEP 
            DO J = I, 3*NUCDEP
               IF (J.EQ.I) THEN 
                  FAC = D0
                  DO K = 1, NBLEN
                     FAC = FAC + VXCRR(K)*RHODER(K,J)**2
     &                         + VXCZZ(K)*ZETDER(K,J)**2
     &                         + VXCRZ(K)*(D2*RHODER(K,J)*ZETDER(K,J))
     &                         + VXCZ(K)*(GDRHOD(K,J,1)**2
     &                                  + GDRHOD(K,J,2)**2
     &                                  + GDRHOD(K,J,3)**2)
                  END DO
               ELSE
                  FAC = D0 
                  DO K = 1, NBLEN
                     FAC = FAC + VXCRR(K)*RHODER(K,I)*RHODER(K,J)
     &                         + VXCZZ(K)*ZETDER(K,I)*ZETDER(K,J)
     &                         + VXCRZ(K)*(RHODER(K,I)*ZETDER(K,J)
     &                                   + RHODER(K,J)*ZETDER(K,I))
     &                         + VXCZ(K)*(GDRHOD(K,I,1)*GDRHOD(K,J,1)
     &                                  + GDRHOD(K,I,2)*GDRHOD(K,J,2)
     &                                  + GDRHOD(K,I,3)*GDRHOD(K,J,3))
                  END DO        
               END IF
               HESDFT(I,J) = HESDFT(I,J) + FAC
            END DO        
         END DO           
         CALL DSCAL(NBLEN,D2,VXCZ(1),1)
      ELSE
         DO I = 1, 3*NUCDEP 
            DO J = I, 3*NUCDEP
               FAC = D0 
               DO K = 1, NBLEN
                  FAC = FAC + VXCRR(K)*RHODER(K,I)*RHODER(K,J)
               END DO        
                  HESDFT(I,J) = HESDFT(I,J) + FAC 
            END DO        
         END DO           
      END IF
C                         
C calculate dF/dp . dp2/(dxa.dxb)
C           dF/dz . dz2/(dxa.dxb) where z = zeta
C Static part doesn't include terms with Sa or Sab.
C
      IF (DOGGA) THEN 
         DO ISYM = 1, NSYM
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
         END DO
      END IF
C
      IF (DOGGA) THEN 
      DO ISYM = 1, NSYM
         DO IXA = 1, 3       
         IA = 0               
         DO ISHELA = 1, KMAX 
            ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + IXA,0,1)
            DO ICOMPA = 1, KHKT(ISHELA)
               IA = IA + 1
               IF (ACTIVE(IA)) THEN
                  CALL DZERO(TMP1,NBLEN*3*NUCDEP)
                  CALL DZERO(TMP2,NBLEN*3*NUCDEP)
                  DO IXB = 1, 3
                     IB = 0
                     DO ISHELB = 1, KMAX 
                     ISCORB = IPTCNT(3*(NCENT(ISHELB) - 1) + IXB,0,1)
                     DO ICOMPB = 1, KHKT(ISHELB)
                        IB = IB + 1
                        IF (ACTIVE(IB).AND.ISCORB.GE.ISCORA) THEN
                           DBA = DMAT(IB,IA) 
                           CALL DAXPY(NBLEN,DBA,GAO(1,IB,IXB+1),1,
     &                                           TMP1(1,ISCORB),1)
                           CALL DAXPY(NBLEN,DBA,GAGDT(1,IB,IXB),1,
     &                                           TMP2(1,ISCORB),1)
                        END IF
                     END DO
                     END DO
                  END DO
                  DO J = 1, NBLEN
                     VXCT1(J) = VXCZ(J)*GAGDT(J,IA,IXA)
     &                        + VXCR(J)*GAO(J,IA,IXA+1)
                     VXCT2(J) = VXCZ(J)*GAO(J,IA,IXA+1) 
                  END DO
                  DO I = ISCORA, 3*NUCDEP
                     FCT = D0 
                     DO J = 1, NBLEN 
                        FCT = FCT + VXCT1(J)*TMP1(J,I) 
     &                            + VXCT2(J)*TMP2(J,I)
                     END DO
                     HESDFT(ISCORA,I) = HESDFT(ISCORA,I) + FCT
                  END DO
               END IF
            END DO
         END DO
      END DO
      END DO
C LDA
      ELSE
         DO ISYM = 1, NSYM 
         DO IXA = 1, 3
            IA = 0
            DO ISHELA = 1, KMAX
               ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + IXA,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IA = IA + 1
                  IF (ACTIVE(IA)) THEN
                     CALL DZERO(TMP1,NBLEN*3*NUCDEP)
                     DO IXB = 1, 3
                        IB = 0
                        DO ISHELB = 1, KMAX
                        ISCORB = IPTCNT(3*(NCENT(ISHELB)-1) + IXB,0,1)
                        DO ICOMPB = 1, KHKT(ISHELB)
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
                        FCT = D0 
                        DO J = 1, NBLEN
                           FCT = FCT + VXCT1(J)*TMP1(J,I)
                        END DO       
                        HESDFT(ISCORA,I) = HESDFT(ISCORA,I) + FCT
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
      IF (DOGGA) THEN 
      DO ISYM = 1, NSYM
         DO IXY = 1, 6
            NXY1 = IXY3(3*IXY - 2)
            NXY2 = IXY3(3*IXY - 1)
            NXY3 = IXY3(3*IXY)
            DO IA = 1, NBAST
               IF (ACTIVE(IA)) THEN
                  DO I  = 1, NBLEN
                     GAGDT(I,IA,IXY)  = GAO(I,IA,NXY1)*GRADAT(I,1)
     &                                + GAO(I,IA,NXY2)*GRADAT(I,2)
     &                                + GAO(I,IA,NXY3)*GRADAT(I,3)
                     IF (IXY.EQ.1) THEN 
                        GAGDT2(I,IA) = GD2(I,IA)*GRADAT(I,1)
     &                               + GD3(I,IA)*GRADAT(I,2)
     &                               + GD4(I,IA)*GRADAT(I,3)
                     END IF
                  END DO
               END IF
            END DO 
         END DO
      END DO
      END IF
C           
      DO ISYM = 1, NSYM
      DO IXY = 1, 6
         NXY1 = IXY3(3*IXY - 2)
         NXY2 = IXY3(3*IXY - 1)
         NXY3 = IXY3(3*IXY)
         JX2  = IX2(IXY) 
         JY2  = IY2(IXY) 
         IP4 = IXY + 4
         IA = 0
         DO ISHELA = 1, KMAX
C gga
            IF (DOGGA) THEN
               TRM1 = D0
               DO ICOMPA = 1, KHKT(ISHELA)
                  IA = IA + 1
                  IF (ACTIVE(IA)) THEN
                     DO I = 1, NBLEN 
                        TRM1 = TRM1 + GAO(I,IA,IP4)*
     &                       (VXCR(I)*GD1(I,IA) + VXCZ(I)*GAGDT2(I,IA))
     &                              + VXCZ(I)*GD1(I,IA)*GAGDT(I,IA,IXY)
                     END DO
                  END IF
               END DO
               ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + JX2,0,1)
               ISCORB = IPTCNT(3*(NCENT(ISHELA) - 1) + JY2,0,1)
               HESDFT(ISCORA,ISCORB) = HESDFT(ISCORA,ISCORB) + TRM1 
C lda
            ELSE
               TRM1 = D0
               DO ICOMPA = 1, KHKT(ISHELA)
                  IA = IA + 1
                  IF (ACTIVE(IA)) THEN
                     DO I = 1, NBLEN
                        TRM1 = TRM1 + VXCR(I)*GD1(I,IA)*GAO(I,IA,IP4)
                     END DO
                  END IF
               END DO
               ISCORA = IPTCNT(3*(NCENT(ISHELA) - 1) + JX2,0,1)
               ISCORB = IPTCNT(3*(NCENT(ISHELA) - 1) + JY2,0,1)
               HESDFT(ISCORA,ISCORB) = HESDFT(ISCORA,ISCORB) + TRM1
            END IF
         END DO
      END DO
      END DO
C
      DO J = 1,3*NUCDEP
         DO I = J+1, 3*NUCDEP
            HESDFT(I,J) = HESDFT(J,I)
         END DO
      END DO
C
      RETURN
      END







