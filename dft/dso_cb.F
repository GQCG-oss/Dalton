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
#include <implicit.h>
#include <codata.h>
#include <mxcent.h>
#include <nuclei.h>
#include <spnout.h>
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
#include <implicit.h>
#include <mxcent.h>
      DIMENSION SPNDSO(MXCOOR,MXCOOR)
#include <nuclei.h>
      DO K = 1, NUCIND*3
         DO L = 1, K-1
            SPNDSO(K,L) = SPNDSO(L  ,K)
         END DO
      END DO
      END
      SUBROUTINE getdsosz(nDSODIM)
#include <implicit.h>
#include <mxcent.h>
      nDSODIM = MXCOOR
      END
