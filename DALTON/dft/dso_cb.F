!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
      SUBROUTINE DSOCB(DSOMAT,NBLLEN,COOR,RHO,WGHT,
     &     RVEC,R3I)
c     Pawel Salek, 20030717
c     computes Dipole Spin-Orbit contribution to spin-spin couplings.
c     WARNING: this does not work for MAXOPR!=0, i.e. when symmetry
C     is present.
c
c     It computes ((r_K . r_L)I -r_K r_L^T)/(r_K^3 r_L^3)
c
c     The elements are only computed for atoms where K>L because the
c     matrix is symmetric. Also, abamag multiplies by ALPHAC in some
c     other place.  The correctness of the code has been verified by
c     comparing against the old DSO code.
c
c     Contribution selection deserves attention since we compute
c     only upper triangle of the matrix.
c
c     Therefore:
c     1. upper triangle of the DSOMAT need to be set.
c     This is done in numdso_finish() f77 routine.
c
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "spnout.h"
      DIMENSION DSOMAT(MXCOOR,MXCOOR),RHO(NBLLEN),WGHT(NBLLEN)
      DIMENSION COOR(3,NBLLEN)
c     Temporary arrays:
      DIMENSION RVEC(NBLLEN,3,*), R3I(NBLLEN,*)
c
      DIMENSION RKL(NBLLEN), RW(NBLLEN)
c
c    fill in rvec and r3
      DO IK = 1, NUCIND
         DO K = 1, NBLLEN
            RVEC(K,1,IK) = COOR(1,K)-CORD(1,IK)
            RVEC(K,2,IK) = COOR(2,K)-CORD(2,IK)
            RVEC(K,3,IK) = COOR(3,K)-CORD(3,IK)
            R3I(K,IK) = (RVEC(K,1,IK)*RVEC(K,1,IK)+
     &           RVEC(K,2,IK)*RVEC(K,2,IK)+
     &           RVEC(K,3,IK)*RVEC(K,3,IK))**(-1.5D0)
         END DO
      END DO
c     compute contributions.
      DO IK = 1, NUCIND
         DO IL = 1, IK-1
         IF(.NOT.NCSPNI(IK).AND..NOT.NCSPNI(IL)) GO TO 10
            DO K = 1, NBLLEN
               RKL(K) = RVEC(K,1,IK)*RVEC(K,1,IL)+
     &                  RVEC(K,2,IK)*RVEC(K,2,IL)+
     &                  RVEC(K,3,IK)*RVEC(K,3,IL)
               RW(K) = RHO(K)*WGHT(K)*R3I(K,IK)*R3I(K,IL)
            END DO
            DO IX=1,3
               IDK = 3*(IK-1)+IX
               DO IY=1,3
                  IDL = 3*(IL-1)+IY
c                 VERIFY ME: IY should go with IL but it does not
c                 if one wants to generate same numbers as the 
c                 analytical version.
                  IF(IX.NE.IY) THEN
                     DO K = 1, NBLLEN
                        DSOMAT(IDL,IDK) = DSOMAT(IDL,IDK)
     &                       -RW(K)*RVEC(K,IY,IK)*RVEC(K,IX,IL)
                     END DO
                  ELSE
                     DO K = 1, NBLLEN
                        DSOMAT(IDL,IDK) = DSOMAT(IDL,IDK)
     &                  + RW(K)*(RKL(K)-RVEC(K,IX,IK)*RVEC(K,IY,IL))
                     END DO
                  END IF
               END DO
            END DO
 10      continue
         END DO
      END DO
      END
C  /* Deck numdso_finish */
      SUBROUTINE NUMDSO_FINISH(SPNDSO)
c     see above for description.
#include "implicit.h"
#include "mxcent.h"
      DIMENSION SPNDSO(MXCOOR,MXCOOR)
#include "nuclei.h"
      DO K = 1, NUCIND*3
         DO L = 1, K-1
            SPNDSO(K,L) = SPNDSO(L  ,K)
         END DO
      END DO
c      write(2,*)'SPNDSO:'
c      CALL OUTPUT(SPNDSO,1,NUCIND*3,1,NUCIND*3,MXCOOR,MXCOOR,1,2)
      END
      SUBROUTINE getdsosz(nDSODIM)
#include "implicit.h"
#include "mxcent.h"
      nDSODIM = MXCOOR
      END
#ifdef VAR_MPI
      SUBROUTINE dsosyncslaves(dmat,nind,work,lwork)
#include "implicit.h"
#include "mpif.h"
#include "codata.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "spnout.h"
#include "inforb.h"
#include "infpar.h"
      DIMENSION DMAT(NBAST,NBAST)
      IF(MYNUM.EQ.0) CALL dftdns(dmat, work, lwork, 0)
      CALL MPI_Bcast(dmat, N2BASX, MPI_DOUBLE_PRECISION,
     &     0, MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nind, 1, my_MPI_INTEGER, 0, MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nucind, 1, my_MPI_INTEGER, 0, MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(NCSPNI,MXCENT, MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      END
#endif
