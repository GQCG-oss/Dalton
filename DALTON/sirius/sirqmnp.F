!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C  /* Deck qmnpmm_fock */
      SUBROUTINE QMNPMM_FOCK(DCAO,DVAO,RCPMAT,RMMMAT,FMAT,EQMNP,WORK,
     &                       LWORK)
C                    
C Purpose:
C     Computes QM/NP/MM contribution to Fock matrix. 
C
C Input:
C   DCAO   - Inactive density matrix.  
C   DVAO   - Active density matrix. 
C   RCPMAT - Inverted Relay matrix or initial quess for MQ vector (NP)
C   RMMMAT - Inverter Relay matrix or initial quess for M vector (MM)  
C   WORK   - Temporary memory array 
C   LWORK  - Size of temporary memory array  
C Output:
C   FMAT   - QM/NP/MM contribution to Fock matrix
C
C Last updated: 22/03/2013 by Z. Rinkevicius.  
C
#include "implicit.h"
#include "dummy.h"
#include "qmnpmm.h"
#include "priunit.h"
#include "inforb.h"
C
      DIMENSION DCAO(*), DVAO(*), RCPMAT(*), RMMMAT(*), FMAT(*), 
     &          WORK(LWORK)
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
      KFREE = 1 
      LFREE = LWORK 
C
      CALL  GETDIM_RELMAT(IDIM,.FALSE.)          
C
      EQMNP = D0
C
      IF (.NOT.MQITER) THEN
         IF (.NOT.(DOMMSUB.AND.DOMMPOL)) THEN 
            CALL MEMGET('REAL',KMQVEC,IDIM,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KFVVEC,IDIM,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KFAO,NNBASX,WORK,KFREE,LFREE)   
C           Determine electric field/potential vector  
            CALL GET_FVVEC(DCAO,DVAO,WORK(KFVVEC),IDIM,WORK(KFREE),
     &                     LFREE)
C           Determine induced momemnts/charges 
            CALL DGEMV('N',IDIM,IDIM,D1,RCPMAT,IDIM,WORK(KFVVEC),1,D0,
     &                 WORK(KMQVEC),1)
            IF (IPRTLVL.GE.15) THEN 
               WRITE(LUPRI,'(/,2X,A)') '*** Computed MQ vector ***'
               CALL OUTPUT(WORK(KMQVEC),1,IDIM,1,1,IDIM,1,1,LUPRI)
            END IF    
C           Compute induced dipoles & charges contribution to 
C           Fock/Kohn-Sham matrix 
            CALL GET_INDMQ_FOCK(DCAO,DVAO,WORK(KMQVEC),IDIM,WORK(KFAO),
     &                          WORK(KFREE),LFREE)
C           Add energy contributions 
            IF (DONPPOL) EQMNP = EQMNP+EESOLMNP+ENSOLMNP 
            IF (DONPCAP) EQMNP = EQMNP+EESOLQNP+ENSOLQNP
C           Computer permanent charges in MM region contribution to 
C           Fock/Kohn-Sham matrix and add energy contributions 
            IF (DOMMSUB) THEN
               CALL GET_PERMQ_FOCK(DCAO,DVAO,WORK(KFAO),WORK(KFREE),
     &                             LFREE)
               EQMNP = EQMNP+EESOLQMM+ENSOLQMM 
            END IF
         ELSE  
c            CALL MEMGET('REAL',KMQVEC,IDIM,WORK,KFREE,LFREE)
c            CALL MEMGET('REAL',KFVVEC,IDIM,WORK,KFREE,LFREE)
cC           Determine electric field/potential vector from QM  
c            CALL GET_FVVEC(DCAO,DVAO,WORK(KFVVEC),IDIM,WORK(KFREE),
c     &                     LFREE)
C           Determine guess for momemnts/charges 
c            CALL DGEMV('N',IDIM,IDIM,D1,RCPMAT,IDIM,WORK(KFVVEC),1,D0,
c     &                 WORK(KMQVEC),1)
c            IF (IPRTLVL.GE.15) THEN 
c               WRITE(LUPRI,'(/,2X,A)') '*** Guess for MQ(NP) vector ***'
c               CALL OUTPUT(WORK(KMQVEC),1,IDIM,1,1,IDIM,1,1,LUPRI)
c            END IF  
c            stop  
C           Compute permanent charges contribution to Fock/Kohn matrix 
         END IF              
C        Pack matrix and add it to Fock matrix 
         CALL PKSYM1(WORK(KFAO),FMAT,NBAS,NSYM,+1)     
C        Compute induced dipoles contribution to Fock/Kohn-Sham matrix 
         CALL MEMREL('QMNPMM_FOCK',WORK,1,1,KFREE,LFREE)                
      ELSE   
 
      END IF   
C
      RETURN
      END
C  /* Deck get_fvvec */
      SUBROUTINE GET_FVVEC(DCAO,DVAO,FVVEC,IDIM,WORK,LWORK)
C                    
C Purpose:
C     Computes electric field/potential vector. 
C
C Input:
C   DCAO   - Inactive density matrix.  
C   DVAO   - Active density matrix. 
C   WORK   - Temporary memory array 
C   LWORK  - Size of temporary memory array.  
C Output:
C   FVVEC  - Electric field/potential vector at NP/MM centers. 
C   IDIM   - Size of electric field/potential vector
C
C Last updated: 22/03/2013 by Z. Rinkevicius.  
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      DIMENSION DCAO(*), DVAO(*), FVVEC(IDIM), WORK(LWORK)
C
      KFREE = 1 
      LFREE = LWORK 
C
      CALL DZERO(FVVEC,IDIM)
C     determine QM region contributions to FV vector 
      CALL GET_QMNUCFV(FVVEC,IDIM)    
      CALL GET_QMELEFV(FVVEC,IDIM,DCAO,DVAO,WORK(KFREE),LFREE)
C       
      CALL GET_QLAGRAN(FVVEC,IDIM)
C     Print final FV vector  
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN 
         WRITE(LUPRI,'(/,2X,A)') '*** Computed FV vector ***'
         CALL OUTPUT(FVVEC,1,IDIM,1,1,IDIM,1,1,LUPRI)
      END IF  
C     Add MM region contribution to FV vector 
      IF (DOMMSUB.AND.(.NOT.DOMMPOL)) THEN
         CALL GET_MMFV(FVVEC,IDIM)
c         Print final FV vector with MM contribution  
          IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN 
            WRITE(LUPRI,'(/,2X,A)') '*** Computed FV+MM vector ***'
           CALL OUTPUT(FVVEC,1,IDIM,1,1,IDIM,1,1,LUPRI)
          END IF    
      END IF                 
C
      RETURN
      END
C  /* Deck get_qmnucfv */
      SUBROUTINE GET_QMNUCFV(FVVEC,IDIM)
C                    
C Purpose:
C     Computes contribution to electric field/potential vector from 
C     nuclei in QM region. 
C
C Output:
C   FVVEC  - Electric field/potential vector at NP/MM centers. 
C   IDIM   - Size of electric field/potential vector
C
C Last updated: 22/03/2013 by Z. Rinkevicius.  
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "mxcent.h"
#include "nuclei.h"
C
      DIMENSION FVVEC(IDIM)
C
      DIMENSION RIJ(3)  
C     Electric field due to nuclei in QM region    
      IF (DONPPOL) THEN 
         DO I=1,TNPATM 
            IOFF = 3*(I-1) 
            DO J=1,NUCIND 
              RIJ(1) = NPCORD(1,I)-CORD(1,J) 
              RIJ(2) = NPCORD(2,I)-CORD(2,J) 
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD 
              FACT = CHARGE(J)/RAD3 
              FVVEC(IOFF+1) = FVVEC(IOFF+1)+FACT*RIJ(1) 
              FVVEC(IOFF+2) = FVVEC(IOFF+2)+FACT*RIJ(2) 
              FVVEC(IOFF+3) = FVVEC(IOFF+3)+FACT*RIJ(3)   
            END DO
         END DO          
      END IF    
C     Determine potential if capacitancy is used 
      IF (DONPCAP) THEN
         ISTART = 0 
         IF (DONPPOL) ISTART = 3*TNPATM 
C        Fix me MM region shift 
         DO I=1,TNPATM 
            IOFF = ISTART+I 
            DO J=1,NUCIND 
              RIJ(1) = NPCORD(1,I)-CORD(1,J) 
              RIJ(2) = NPCORD(2,I)-CORD(2,J) 
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              FVVEC(IOFF) = FVVEC(IOFF)+CHARGE(J)/RAD 
            END DO
C           Fix me: add MM region polarizable points 
         END DO   
      END IF       
C
      RETURN
      END
C  /* Deck get_qmelefv */
      SUBROUTINE GET_QMELEFV(FVVEC,IDIM,DCAO,DVAO,WORK,LWORK)
C                    
C Purpose:
C     Computes contributio to electric field/potential vector 
C     for electrons in QM regions. 
C
C Input:
C   DCAO   - Inactive density matrix.  
C   DVAO   - Active density matrix. 
C   WORK   - Temporary memory array 
C   LWORK  - Size of temporary memory array.  
C Output:
C   FVVEC  - Electric field/potential vector at NP/MM centers. 
C   IDIM   - Size of electric field/potential vector
C
C Last updated: 22/03/2013 by Z. Rinkevicius.  
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"

C
      DIMENSION DCAO(*), DVAO(*), FVVEC(IDIM), WORK(LWORK)      
C
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3)
C
      KFREE = 1 
      LFREE = LWORK 
C     Save origin coordinates 
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)   
C     Electric field due to electrons in QM region
      IF (DONPPOL) THEN  
         CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
C        Set integrals evaluation flags 
         KPATOM = 0
         NOCOMP = 3
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE. 
C 
         DO I=1,TNPATM
             IOFF = (I-1)*3
             DIPORG(1) = NPCORD(1,I) 
             DIPORG(2) = NPCORD(2,I) 
             DIPORG(3) = NPCORD(3,I) 
             CALL DZERO(WORK(KINTAO),3*NNBASX) 
             CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             IF (NISHT.GT.0) THEN 
                RELFLD = DDOT(NNBASX,DCAO,1,WORK(KINTAO),1) 
                FVVEC(IOFF+1) = FVVEC(IOFF+1)+RELFLD 
                RELFLD = DDOT(NNBASX,DCAO,1,WORK(KINTAO+NNBASX),1)  
                FVVEC(IOFF+2) = FVVEC(IOFF+2)+RELFLD  
                RELFLD = DDOT(NNBASX,DCAO,1,WORK(KINTAO+2*NNBASX),1)  
                FVVEC(IOFF+3) = FVVEC(IOFF+3)+RELFLD  
             END IF  
             IF (NASHT.GT.0) THEN 
                RELFLD = DDOT(NNBASX,DVAO,1,WORK(KINTAO),1) 
                FVVEC(IOFF+1) = FVVEC(IOFF+1)+RELFLD 
                RELFLD = DDOT(NNBASX,DVAO,1,WORK(KINTAO+NNBASX),1)  
                FVVEC(IOFF+2) = FVVEC(IOFF+2)+RELFLD  
                RELFLD = DDOT(NNBASX,DVAO,1,WORK(KINTAO+2*NNBASX),1)  
                FVVEC(IOFF+3) = FVVEC(IOFF+3)+RELFLD  
             END IF 
         END DO    
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_QMELEFV',WORK,1,1,KFREE,LFREE)  
      END IF 
C     Potential due to electrons in QM region
      IF (DONPCAP) THEN 
         ISTART = 0 
         IF (DONPPOL) ISTART = 3*TNPATM 
C        Fix me MM region shift          
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
C        Set integrals evaluation flags 
         KPATOM = 0
         NOCOMP = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE. 
C         
         DO I=1,TNPATM
             IOFF = ISTART+I 
             DIPORG(1) = NPCORD(1,I) 
             DIPORG(2) = NPCORD(2,I) 
             DIPORG(3) = NPCORD(3,I) 
             CALL DZERO(WORK(KINTAO),NNBASX) 
             CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             IF (NISHT.GT.0) THEN 
                RVPOT = DDOT(NNBASX,DCAO,1,WORK(KINTAO),1) 
                FVVEC(IOFF) = FVVEC(IOFF)+RVPOT 
             END IF  
             IF (NASHT.GT.0) THEN 
                RVPOT = DDOT(NNBASX,DVAO,1,WORK(KINTAO),1) 
                FVVEC(IOFF) = FVVEC(IOFF)+RVPOT   
             END IF 
         END DO    
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_QMELEFV',WORK,1,1,KFREE,LFREE)  
      END IF
C     Restore origin coordinates 
      DIPORG(1) = RSAVORG(1) 
      DIPORG(2) = RSAVORG(2) 
      DIPORG(3) = RSAVORG(3)   
C
      RETURN
      END
C  /* Deck get_qlag */
      SUBROUTINE GET_QLAGRAN(FVVEC,IDIM)
C                    
C Purpose:
C     Determines charge contrain fro electric/field potential vector. 
C
C Output:
C   FVVEC  - Electric field/potential vector at NP/MM centers. 
C   IDIM   - Size of electric field/potential vector
C
C Last updated: 22/03/2013 by Z. Rinkevicius.  
C
#include "implicit.h"
#include "qmnpmm.h"
C
      DIMENSION FVVEC(IDIM)
C
      IF (DONPCAP) THEN
         DO I=1,TNPBLK 
            FVVEC(IDIM) = FVVEC(IDIM)+NPCHRG(I)
         END DO 
      END IF       
C
      RETURN
      END
C  /* Deck get_mmfv */
      SUBROUTINE GET_MMFV(FVVEC,IDIM)
C                    
C Purpose:
C     Computes contribution to electric field/potential vector from 
C     permanent point charges in MM region. 
C
C Output:
C   FVVEC  - Electric field/potential vector at NP/MM centers. 
C   IDIM   - Size of electric field/potential vector
C
C Last updated: 22/03/2013 by Z. Rinkevicius.  
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      DIMENSION FVVEC(IDIM)
C
      DIMENSION RIJ(3)  
C     Electric field due to point charges in MM region    
      IF (DONPPOL) THEN 
         DO I=1,TNPATM 
            IOFF = 3*(I-1) 
            DO J=1,TMMATM 
              RIJ(1) = NPCORD(1,I)-MMCORD(1,J) 
              RIJ(2) = NPCORD(2,I)-MMCORD(2,J) 
              RIJ(3) = NPCORD(3,I)-MMCORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD 
              FACT = MMFM0(MMFTYP(J))/RAD3 
              FVVEC(IOFF+1) = FVVEC(IOFF+1)+FACT*RIJ(1) 
              FVVEC(IOFF+2) = FVVEC(IOFF+2)+FACT*RIJ(2) 
              FVVEC(IOFF+3) = FVVEC(IOFF+3)+FACT*RIJ(3)   
            END DO
         END DO          
      END IF    
C     Determine potential if capacitancy is used 
      IF (DONPCAP) THEN
         ISTART = 0 
         IF (DONPPOL) ISTART = 3*TNPATM 
C        Fix me MM region shift 
         DO I=1,TNPATM 
            IOFF = ISTART+I 
            DO J=1,TMMATM 
              RIJ(1) = NPCORD(1,I)-MMCORD(1,J) 
              RIJ(2) = NPCORD(2,I)-MMCORD(2,J) 
              RIJ(3) = NPCORD(3,I)-MMCORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              FVVEC(IOFF) = FVVEC(IOFF)+MMFM0(MMFTYP(J))/RAD 
            END DO
C           Fix me: add MM region polarizable points 
         END DO   
      END IF       
C
      RETURN
      END
C  /* Deck get_fvvec */  
      SUBROUTINE GET_INDMQ_FOCK(DCAO,DVAO,MQVEC,IDIM,FCAO,WORK,LWORK)
C 
C Purpose:
C     Computes electric field/potential vector.
C
C Input:
C   DCAO   - Inactive density matrix. 
C   DVAO   - Active density matrix.     
C   MQVEC  - Induced dipoles & charges vector. 
C   IDIM   - Size of induced dipoles & charges vector.    
C   WORK   - Temporary memory array                                                                        
C   LWORK  - Size of temporary memory array.                                                               
C Output:                                                                                                  
C   FCAO   - Computed Fock matrix (AO basis).                                            
C                                                                                                          
C Last updated: 22/03/2013 by Z. Rinkevicius.                                                             
C                                                                                                         
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "pi.h"
C
      DOUBLE PRECISION MQVEC
C
      DIMENSION DCAO(*), DVAO(*), MQVEC(IDIM), FCAO(*), WORK(LWORK)   
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, DMP5 = -0.5D0 )
      PARAMETER (DP5 = 0.5D0, D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (D13 = 1.0D0/3.0D0)
C
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3), RIJ(3)
C
      KFREE = 1
      LFREE = LWORK
C
      CALL DZERO(FCAO,NNBASX)          
C
C     Save origin coordinates 
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3) 
C
      ENSOLQNP = D0 
      EESOLQNP = D0 
      ENSOLMNP = D0 
      EESOLMNP = D0
C     Induced dipole moment in NP region interaction with QM region 
      IF (DONPPOL.AND.NOVDAMP) THEN 
C        Electronic interaction part 
         CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
C        Set integrals evaluation flags 
         KPATOM = 0
         NOCOMP = 3
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.          
         DO I=1,TNPATM
             IOFF = 3*(I-1) 
             DIPORG(1) = NPCORD(1,I) 
             DIPORG(2) = NPCORD(2,I) 
             DIPORG(3) = NPCORD(3,I)
             CALL DZERO(WORK(KINTAO),3*NNBASX) 
             CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             CALL DSCAL(NNBASX,MQVEC(IOFF+1),WORK(KINTAO),1)
             CALL DSCAL(NNBASX,MQVEC(IOFF+2),WORK(KINTAO+NNBASX),1)
             CALL DSCAL(NNBASX,MQVEC(IOFF+3),WORK(KINTAO+2*NNBASX),1)
             CALL DAXPY(NNBASX,DM1,WORK(KINTAO),1,FCAO,1)
             CALL DAXPY(NNBASX,DM1,WORK(KINTAO+NNBASX),1,FCAO,1)
             CALL DAXPY(NNBASX,DM1,WORK(KINTAO+2*NNBASX),1,FCAO,1)
             IF (NISHT.GT.0) THEN 
               ELEX = DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)
               ELEY = DDOT(NNBASX,DCAO,1,WORK(KINTAO+NNBASX),1)
               ELEZ = DDOT(NNBASX,DCAO,1,WORK(KINTAO+2*NNBASX),1)
               EESOLMNP = EESOLMNP+ELEX+ELEY+ELEZ
             END IF
             IF (NASHT.GT.0) THEN 
               ELEX = DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)
               ELEY = DDOT(NNBASX,DVAO,1,WORK(KINTAO+NNBASX),1)
               ELEZ = DDOT(NNBASX,DVAO,1,WORK(KINTAO+2*NNBASX),1)
               EESOLMNP = EESOLMNP+ELEX+ELEY+ELEZ
             END IF     
         END DO     
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_INDMQ_FOCK',WORK,1,1,KFREE,LFREE)
         EESOLMNP = DMP5*EESOLMNP
C        Nuclear interaction part   
         DO I=1,TNPATM
             IOFF = 3*(I-1) 
             ENUCFX = D0 
             ENUCFY = D0
             ENUCFZ = D0
             DO J=1,NUCIND 
              RIJ(1) = NPCORD(1,I)-CORD(1,J) 
              RIJ(2) = NPCORD(2,I)-CORD(2,J) 
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD 
              FACT = CHARGE(J)/RAD3 
              ENUCFX = ENUCFX + FACT*RIJ(1) 
              ENUCFY = ENUCFY + FACT*RIJ(2)
              ENUCFZ = ENUCFZ + FACT*RIJ(3)
            END DO
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+1)*ENUCFX 
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+2)*ENUCFY
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+3)*ENUCFZ   
         END DO 
         ENSOLMNP = DMP5*ENSOLMNP                      
      END IF 
C     Induced charges in NP region interaction with QM region 
      IF (DONPCAP.AND.NOVDAMP) THEN 
         ISTART = 0 
         IF (DONPPOL) ISTART = 3*TNPATM 
C        Fix me MM region shift   
C        Electronic interaction part 
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
C        Set integrals evaluation flags 
         KPATOM = 0
         NOCOMP = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE. 
C 
         DO I=1,TNPATM
             IOFF = ISTART+I 
             DIPORG(1) = NPCORD(1,I) 
             DIPORG(2) = NPCORD(2,I) 
             DIPORG(3) = NPCORD(3,I) 
             CALL DZERO(WORK(KINTAO),NNBASX) 
             CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)

             WRITE(LUPRI,'(/,2X,A)') '*** jaimes AO matrix ***'
             CALL OUTPAK(WORK(KINTAO),NBAST,1,LUPRI) 
             CALL DSCAL(NNBASX,MQVEC(IOFF),WORK(KINTAO),1)
             CALL DAXPY(NNBASX,D1,WORK(KINTAO),1,FCAO,1)
             IF (NISHT.GT.0) THEN
               EESOLQNP = EESOLQNP+DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)    
             END IF 
             IF (NASHT.GT.0) THEN
               EESOLQNP = EESOLQNP+DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)    
             END IF 
         END DO    
         RUNQM3 = .FALSE. 
         CALL MEMREL('GET_INDMQ_FOCK',WORK,1,1,KFREE,LFREE)
C   
         EESOLQNP = DP5*EESOLQNP        
C        Nuclear interaction part   
         DO I=1,TNPATM
            IOFF = ISTART+I  
            DO J=1,NUCIND 
              RIJ(1) = NPCORD(1,I)-CORD(1,J) 
              RIJ(2) = NPCORD(2,I)-CORD(2,J) 
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              ENSOLQNP = ENSOLQNP+MQVEC(IOFF)*CHARGE(J)/RAD 
            END DO   
         END DO 
         ENSOLQNP = DP5*ENSOLQNP
      END IF 
C        
C     Restore origin coordinates 
      DIPORG(1) = RSAVORG(1) 
      DIPORG(2) = RSAVORG(2) 
      DIPORG(3) = RSAVORG(3)   
C 
C     DAMPED INTERACTION CASE 
C
C     Induced dipole moment in NP region interaction (damped) with QM region 
      IF (DONPPOL.AND.(.NOT.NOVDAMP)) THEN
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
         CALL COMP_DAMPVMAT(WORK(KINTAO), MQVEC, WORK(KFREE), LFREE)
         IF (NISHT.GT.0) THEN
            EESOLQNP = EESOLQNP+DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)    
         END IF 
         IF (NASHT.GT.0) THEN
            EESOLQNP = EESOLQNP+DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)    
         END IF  
         CALL DAXPY(NNBASX,D1,WORK(KINTAO),1,FCAO,1) 
         CALL MEMREL('GET_INDMQ_FOCK',WORK,1,1,KFREE,LFREE)
         EESOLQNP = DP5*EESOLQNP 
C        Nuclear interaction part: induced dipoles and induced charges    
         RDIM = DSQRT(D2)/SQRTPI
C        Nuclear interaction part   
         DO I=1,TNPATM
             IOFF = 3*(I-1) 
             ENUCFX = D0 
             ENUCFY = D0
             ENUCFZ = D0
             RIPOL = NPFPOL(NPFTYP(I))/D3
             RDVAL = (RDIM*RIPOL)**D13 
             DO J=1,NUCIND 
              RIJ(1) = NPCORD(1,I)-CORD(1,J) 
              RIJ(2) = NPCORD(2,I)-CORD(2,J) 
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD 
              FACT = CHARGE(J)/RAD3
              RFACT = D2*RAD*DEXP(-(RAD*RAD)/(RDVAL*RDVAL))
              RFACT = -RFACT/(SQRTPI*RDVAL) 
              RFACT = RFACT + DERF(RAD/RDVAL)
              ENUCFX = ENUCFX + FACT*RIJ(1)*RFACT 
              ENUCFY = ENUCFY + FACT*RIJ(2)*RFACT
              ENUCFZ = ENUCFZ + FACT*RIJ(3)*RFACT
            END DO
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+1)*ENUCFX 
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+2)*ENUCFY
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+3)*ENUCFZ   
         END DO 
         ENSOLMNP = DMP5*ENSOLMNP
C  
         ISTART = 0 
         IF (DONPPOL) ISTART = 3*TNPATM
         DO I=1,TNPATM
            IOFF = ISTART+I
            RQVAL = RDIM*NPFCAP(NPFTYP(I))      
            DO J=1,NUCIND 
              RIJ(1) = NPCORD(1,I)-CORD(1,J) 
              RIJ(2) = NPCORD(2,I)-CORD(2,J) 
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RFACT = DERF(RAD/RQVAL)
              ENSOLQNP = ENSOLQNP+MQVEC(IOFF)*CHARGE(J)*RFACT/RAD 
            END DO   
         END DO 
         ENSOLQNP = DP5*ENSOLQNP
      END IF   
      RETURN 
      END 
C  /* Deck get_permq_fock */  
      SUBROUTINE GET_PERMQ_FOCK(DCAO,DVAO,FCAO,WORK,LWORK)
C 
C Purpose:
C     Computes permanent charges in MM region contribution to Fock matrix. 
C
C Input:
C   DCAO   - Inactive density matrix. 
C   DVAO   - Active density matrix.     
C   WORK   - Temporary memory array                                                                        
C   LWORK  - Size of temporary memory array.                                                               
C Output:                                                                                                  
C   FCAO   - Computed Fock matrix (AO basis).                                            
C                                                                                                          
C Last updated: 22/03/2013 by Z. Rinkevicius.                                                             
C                                                                                                         
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
C
      DIMENSION DCAO(*), DVAO(*), FCAO(*), WORK(LWORK)   
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, DMP5 = -0.5D0 )
      PARAMETER (DP5 = 0.5D0)
C
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3), RIJ(3)
C
      KFREE = 1
      LFREE = LWORK
C     Save origin coordinates 
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3) 
C
      ENSOLQMM = D0 
      EESOLQMM = D0 
C     Induced charges in NP region interaction with QM region 
C     Fix me MM region shift   
C     Electronic interaction part 
      CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
C     Set integrals evaluation flags 
      KPATOM = 0
      NOCOMP = 1
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      RUNQM3 = .TRUE. 
C 
      DO I=1,TMMATM
         DIPORG(1) = MMCORD(1,I) 
         DIPORG(2) = MMCORD(2,I) 
         DIPORG(3) = MMCORD(3,I) 
         CALL DZERO(WORK(KINTAO),NNBASX) 
         CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),
     &               LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     &               TRIMAT,DUMMY,EXP1VL,DUMMY,0)
          CALL DSCAL(NNBASX,MMFM0(MMFTYP(I)),WORK(KINTAO),1) 
          CALL DAXPY(NNBASX,D1,WORK(KINTAO),1,FCAO,1)
          IF (NISHT.GT.0) THEN
              EESOLQMM = EESOLQMM+DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)    
          END IF 
          IF (NASHT.GT.0) THEN
              EESOLQMM = EESOLQMM+DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)    
          END IF  
      END DO    
      RUNQM3 = .FALSE. 
      CALL MEMREL('GET_PERMQ_FOCK',WORK,1,1,KFREE,LFREE)
C     Nuclear interaction part   
      DO I=1,TMMATM
         DO J=1,NUCIND 
            RIJ(1) = MMCORD(1,I)-CORD(1,J) 
            RIJ(2) = MMCORD(2,I)-CORD(2,J) 
            RIJ(3) = MMCORD(3,I)-CORD(3,J)
            RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
            ENSOLQMM = ENSOLQMM+MMFM0(MMFTYP(I))*CHARGE(J)/RAD 
         END DO   
      END DO 
C     Restore origin coordinates 
      DIPORG(1) = RSAVORG(1) 
      DIPORG(2) = RSAVORG(2) 
      DIPORG(3) = RSAVORG(3)   
C 
      RETURN 
      END 
C     /* Deck scf_qmnpmm_out */  
      SUBROUTINE SCF_QMNPMM_OUT() 
C Purpose:
C     Prints QM/NP/MM contribution to SCF energy and it's decomposition into 
C     components.
C
C Last updated: 22/03/2013 by Z. Rinkevicius.        
C     
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "infopt.h"
C
      IF (IPRTLVL.GE.15) THEN
        WRITE(LUPRI,'(/5X,A/)') 
     &        'QM/NP/MM calculation converged     :'   
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'NP charges nuclear contribution    :', ENSOLQNP
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'NP charges electronic contribution :', EESOLQNP
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'NP dipoles nuclear contribution    :', ENSOLMNP
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'NP dipoles electronic contribution :', EESOLMNP
        IF (DOMMSUB) THEN 
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'MM charges nuclear contribution    :', ENSOLQMM
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'MM charges electronic contribution :', EESOLQMM
        END IF   
        WRITE(LUPRI,'(5X,A,F20.12/)') 
     &        'Total QM/NP/MM energy              :', ESOLT       
      ELSE 
        WRITE(LUPRI,'(/5X,A/)') 
     &        'QM/NP/MM calculation converged :'   
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'NP charges contribution        :', 
     &         ENSOLQNP+EESOLQNP  
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'NP dipoles contribution        :', 
     &        ENSOLMNP+EESOLMNP
        IF (DOMMSUB) THEN 
        WRITE(LUPRI,'(5X,A,F20.12)') 
     &        'MM charges contribution        :', 
     &        ENSOLQMM+EESOLQMM   
        END IF 
        WRITE(LUPRI,'(5X,A,F20.12/)') 
     &        'Total QM/NP/MM energy          :', ESOLT  
      END IF
C
      RETURN 
      END    



