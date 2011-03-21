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
      SUBROUTINE DFTMOLGRAD(WORK,LWORK,IPRINT)
C
C     T. Helgaker and P. Salek, October 2003
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <inforb.h>
#include <priunit.h>
#include "dftmolgrad.h"
      DIMENSION WORK(LWORK)
#include <energy.h>
      EXTERNAL DFTMOLGRD
      EXTERNAL DFTMOLGRDAB 
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      DOGGA = DFT_ISGGA()
C      
      CALL DZERO(GRADFT,MXCOOR)
C
      IF (NASHT .GT. 0) THEN 
         NDMAT = 2
      ELSE 
         NDMAT = 1
      END IF  
      KDMAT = 1
      KLAST = KDMAT + NDMAT*NBAST*NBAST
      LWRK  = LWORK - KLAST +1
      IF(KLAST.GT.LWORK) CALL QUIT('NOMEM IN DFTGRAD')
      IF (NASHT .GT. 0) THEN
         CALL DFTDNSAB(WORK(KDMAT),WORK(KDMAT+NBAST*NBAST), 
     &                 WORK(KLAST),LWRK,0) 
      ELSE    
         CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      END IF   
      CALL KICK_SLAVES_GRAD(NBAST,NDMAT,WORK(KDMAT),IPRINT)
      IF (NASHT .GT. 0) THEN 
         CALL DFTINT(WORK(KDMAT),NDMAT,1,.FALSE.,WORK(KLAST),LWRK,
     &               DFTMOLGRDAB,WORK(KDMAT),ELE)
      ELSE    
         CALL DFTINT(WORK(KDMAT),NDMAT,1,.FALSE.,WORK(KLAST),LWRK,
     &               DFTMOLGRD,WORK(KDMAT),ELE)
      END IF  
      CALL GRADSLAVE_COLLECT(GRADFT,WORK(KLAST),LWRK)
      RETURN
      END
C
      SUBROUTINE DFTMOLGRD(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                     COORD,WGHT,DMAT)
C
C     T. Helgaker oct 2003
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
#include "dftmolgrad.h"
#include <inforb.h>
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          DMAT(NBAST,NBAST)
C
#include <dftinf.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
c
      LOGICAL ACTIVE(0:NBAST)
      DIMENSION VXC(NBLEN), VXB(NBLEN), VX(5)
      DIMENSION TMP(NBLEN), KVALS(3,3)
      DATA ((KVALS(I,J), I = 1,3), J = 1,3) /1, 2, 3, 2, 4, 5, 3, 5, 6/
C
#include <ibtfun.h>
C 
      DO I = 0, NBAST
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
C     Exchange-correlation contribution to molecular gradient 
C
      IF (DOGGA) THEN
         DO I = 1, NBLEN
            GRDNRM = SQRT(GRADA(1,I)**2 + GRADA(2,I)**2 + GRADA(3,I)**2)
            CALL DFTPTF0(RHOA(I),GRDNRM,WGHT(I),VX)
            VXC(I) = D2*VX(1)
            VXB(I) = D2*VX(2)/GRDNRM
         END DO
      ELSE
         DO I = 1, NBLEN
            CALL DFTPTF0(RHOA(I),D0,WGHT(I),VX)
            VXC(I) = D2*VX(1)
         END DO
      END IF
C
      DO IX = 1, 3
         DO ISYM = 1, NSYM
            IORBA = 0
            DO ISHELA = 1, KMAX
               ISCOOR = IPTCNT(3*(NCENT(ISHELA) - 1) + IX,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IORBA = IORBA + 1
                  IA = IPTSYM(IORBA,ISYM-1)
                  KA = IPTSYM(IORBA,IBTXOR(ISYM-1,ISYMAX(IX,1)))
                  IF (ACTIVE(IA) .AND. KA.GT.0) THEN
                     IF (.NOT.DOGGA) THEN
                        DO I = 1, NBLEN
                           TMP(I) = D0
                        END DO
                        DO IBLB = 1, NBLCNT(ISYM)
                        DO IB=NBLOCKS(1,IBLB,ISYM),NBLOCKS(2,IBLB,ISYM) 
                           DO I = 1, NBLEN
                              TMP(I) = TMP(I) + GAO(I,IB,1)*DMAT(IB,IA)
                           END DO
                        END DO
                        END DO
                        FRC = D0
                        DO I = 1, NBLEN
                           FRC = FRC + VXC(I)*TMP(I)*GAO(I,KA,IX+1)
                        END DO
                     ELSE
                        FRC = D0
                        K1 = KVALS(1,IX) + 4
                        K2 = KVALS(2,IX) + 4
                        K3 = KVALS(3,IX) + 4
                        DO I = 1, NBLEN
                           GA  = GAO(I,KA,IX+1)
                           GAX = GRADA(1,I)*GA
                           GAY = GRADA(2,I)*GA
                           GAZ = GRADA(3,I)*GA
                           GA2 = GRADA(1,I)*GAO(I,KA,K1) 
     &                         + GRADA(2,I)*GAO(I,KA,K2)
     &                         + GRADA(3,I)*GAO(I,KA,K3)
                           GD = D0
                           GF = D0
                           DO IBLB = 1, NBLCNT(ISYM)
                              ISTRB = NBLOCKS(1,IBLB,ISYM)
                              IENDB = NBLOCKS(2,IBLB,ISYM) 
                              DO IB = ISTRB, IENDB
                                 GD = GD + DMAT(IB,IA)*GAO(I,IB,1)
                                 GF = GF + DMAT(IB,IA)*(GAO(I,IB,1)*GA2
     &                                                + GAO(I,IB,2)*GAX
     &                                                + GAO(I,IB,3)*GAY
     &                                                + GAO(I,IB,4)*GAZ)
                              END DO
                           END DO
                           FRC = FRC + VXC(I)*GD*GA + VXB(I)*GF
                        END DO
                     END IF
                     GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC 
                   END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
C
      SUBROUTINE DFTMOLGRDAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,RHOB, 
     &                       GRADA,GRADB,COORD,WGHT,DMAT)
C                           
C     Open-shell adaptation of DFTMOLGRD 
C                                        
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
#include "dftmolgrad.h"
#include <inforb.h>
C                                                                                                                                                                                                                
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          RHOB(NBLEN), GRADB(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          DMAT(NBAST,NBAST,2)
C                                                                                                                                                                                                                
#include <dftinf.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
c                                                                                                                                                                                                                
      LOGICAL ACTIVE(0:NBAST)
      DIMENSION VXA(NBLEN), VXB(NBLEN), VXGA(NBLEN), VXGB(NBLEN), 
     &          VXZ(NBLEN),VX(5)  
      DIMENSION TMPA(NBLEN), TMPB(NBLEN), KVALS(3,3)
      DATA ((KVALS(I,J), I = 1,3), J = 1,3) /1, 2, 3, 2, 4, 5, 3, 5, 6/
C                                                                          
#include <ibtfun.h>
C
       DO I = 0, NBAST
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
C     Exchange-correlation contribution to molecular gradient
C                           
      IF (DOGGA) THEN
         DO I = 1, NBLEN
            GRDA  = SQRT(GRADA(1,I)**2+GRADA(2,I)**2+GRADA(3,I)**2)
            GRDB  = SQRT(GRADB(1,I)**2+GRADB(2,I)**2+GRADB(3,I)**2)
            GRDAB   = GRADA(1,I)*GRADB(1,I)+GRADA(2,I)*GRADB(2,I) 
     &              + GRADA(3,I)*GRADB(3,I)
            CALL VXCFAB(RHOA(I),RHOB(I),GRDA,GRDB,GRDAB,WGHT(I),VX)  
            VXA(I)  = D2*VX(1)
            VXB(I)  = D2*VX(2)
            VXGA(I) = D2*VX(3)/GRDA
            VXGB(I) = D2*VX(4)/GRDB
            VXZ(I)  = D2*VX(5)
         END DO
      ELSE
         DO I = 1, NBLEN
             CALL VXCFAB(RHOA(I),RHOB(I),D0,D0,D0,WGHT(I),VX)
             VXA(I) = D2*VX(1) 
             VXB(I) = D2*VX(2)     
         END DO
      END IF
C 
      DO IX = 1, 3
         DO ISYM = 1, NSYM
            IORBA = 0
            DO ISHELA = 1, KMAX
               ISCOOR = IPTCNT(3*(NCENT(ISHELA) - 1) + IX,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IORBA = IORBA + 1
                  IA = IPTSYM(IORBA,ISYM-1)
                  KA = IPTSYM(IORBA,IBTXOR(ISYM-1,ISYMAX(IX,1)))
                  IF (ACTIVE(IA) .AND. KA.GT.0) THEN
                     IF (.NOT.DOGGA) THEN
                        DO I = 1, NBLEN
                           TMPA(I) = D0
                           TMPB(I) = D0
                        END DO
                        DO IBLB = 1, NBLCNT(ISYM)
                        DO IB=NBLOCKS(1,IBLB,ISYM),NBLOCKS(2,IBLB,ISYM)
                           DO I = 1, NBLEN
                              TMPA(I)=TMPA(I)+GAO(I,IB,1)*DMAT(IB,IA,1)
                              TMPB(I)=TMPB(I)+GAO(I,IB,1)*DMAT(IB,IA,2)
                           END DO
                        END DO
                        END DO
                        FRC = D0
                        DO I = 1, NBLEN
                           FRC = FRC + VXA(I)*TMPA(I)*GAO(I,KA,IX+1)
     &                         + VXB(I)*TMPB(I)*GAO(I,KA,IX+1)      
                        END DO
                     ELSE
                        FRC = D0
                        K1 = KVALS(1,IX) + 4
                        K2 = KVALS(2,IX) + 4
                        K3 = KVALS(3,IX) + 4
                        DO I = 1, NBLEN
                           GA  = GAO(I,KA,IX+1)
                           GAX = GRADA(1,I)*GA
                           GAY = GRADA(2,I)*GA
                           GAZ = GRADA(3,I)*GA
                           GA2 = GRADA(1,I)*GAO(I,KA,K1)
     &                         + GRADA(2,I)*GAO(I,KA,K2)
     &                         + GRADA(3,I)*GAO(I,KA,K3)
                           GBX = GRADB(1,I)*GA
                           GBY = GRADB(2,I)*GA
                           GBZ = GRADB(3,I)*GA
                           GB2 = GRADB(1,I)*GAO(I,KA,K1)
     &                         + GRADB(2,I)*GAO(I,KA,K2)
     &                         + GRADB(3,I)*GAO(I,KA,K3) 
                           GDA = D0
                           GDB = D0
                           GFA = D0
                           GFB = D0
                           GGA = D0
                           GGB = D0
                           DO IBLB = 1, NBLCNT(ISYM)
                              ISTRB = NBLOCKS(1,IBLB,ISYM)
                              IENDB = NBLOCKS(2,IBLB,ISYM)
                              DO IB = ISTRB, IENDB
                                 GDA = GDA + DMAT(IB,IA,1)*GAO(I,IB,1)
                                 GDB = GDB + DMAT(IB,IA,2)*GAO(I,IB,1) 
                                 GFA = GFA + DMAT(IB,IA,1)*(GAO(I,IB,1)
     &                               * GA2 + GAO(I,IB,2)*GAX
     &                               + GAO(I,IB,3)*GAY 
     &                               + GAO(I,IB,4)*GAZ)
                                 GFB = GFB + DMAT(IB,IA,2)*(GAO(I,IB,1)
     &                               * GB2 + GAO(I,IB,2)*GBX
     &                               + GAO(I,IB,3)*GBY
     &                               + GAO(I,IB,4)*GBZ)
                                 GGA = GGA + DMAT(IB,IA,1)*(GAO(I,IB,1)
     &                               * GB2 + GAO(I,IB,2)*GBX
     &                               + GAO(I,IB,3)*GBY
     &                               + GAO(I,IB,4)*GBZ) 
                                 GGB = GGB + DMAT(IB,IA,2)*(GAO(I,IB,1)
     &                               * GA2 + GAO(I,IB,2)*GAX
     &                               + GAO(I,IB,3)*GAY
     &                               + GAO(I,IB,4)*GAZ)
                              END DO
                           END DO
                           FRC = FRC + GA*(VXA(I)*GDA + VXB(I)*GDB)
     &                         + VXGA(I)*GFA + VXGB(I)*GFB
     &                         + VXZ(I)*(GGA+GGB) 
                        END DO
                     END IF
                     GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC
                   END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END  
C
      SUBROUTINE KICK_SLAVES_GRAD(NBAST,NDMAT,DMAT,IPRINT)
#if defined (VAR_MPI)
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "infpar.h"
#include "mpif.h"
C
      DIMENSION DMAT(NBAST,NBAST,NDMAT)
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = 6
         CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL MPI_BCAST(NDMAT,1,my_MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DMAT,NDMAT*NBAST*NBAST,MPI_DOUBLE_PRECISION,0,
     &                  MPI_COMM_WORLD,IERR)        
      END IF
      RETURN
#endif
      END
#if defined (VAR_MPI)
      SUBROUTINE DFT_GRAD_SLAVE(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "dftmolgrad.h"
#include "energy.h"
#include "mpif.h"
      DIMENSION WORK(LWORK)
      EXTERNAL DFTMOLGRD 
      EXTERNAL DFTMOLGRDAB
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      DOGGA = DFT_ISGGA()
      KDMAT  = 1
      KFREE  = KDMAT + 2*N2BASX
      IF (KFREE .GT. LWORK)CALL STOPIT('DFT_GRAD_SLAVE',' ',KFREE,LWORK)
      LFREE = LWORK - KFREE + 1
      CALL DFTINTBCAST
      CALL MPI_BCAST(NDMAT,1,my_MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WORK(KDMAT),NDMAT*NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &               0,MPI_COMM_WORLD,IERR)       
      CALL DZERO(GRADFT,MXCOOR)
      IF (NDMAT .EQ. 2 ) THEN 
         CALL DFTINT(WORK(KDMAT),NDMAT,1,.FALSE.,WORK(KFREE),LFREE,
     &               DFTMOLGRDAB,WORK(KDMAT),ELE)
      ELSE 
        CALL DFTINT(WORK(KDMAT),NDMAT,1,.FALSE.,WORK(KFREE),LFREE,
     &               DFTMOLGRD,WORK(KDMAT),ELE)
      END IF 
      CALL GRADSLAVE_COLLECT(GRADFT,WORK(KFREE),LFREE)
      RETURN
      END
#endif
      SUBROUTINE GRADSLAVE_COLLECT(GRADMOL,WORK,LWORK)
#if defined (VAR_MPI)
#include "implicit.h"
#include "mxcent.h"
#include "mpif.h"
      DIMENSION GRADMOL(MXCOOR), WORK(LWORK)
      CALL DCOPY(MXCOOR,GRADMOL,1,WORK(1),1)
      CALL MPI_Reduce(WORK,GRADMOL,MXCOOR,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      RETURN
#endif
      END


