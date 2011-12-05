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
      SUBROUTINE DIISWEIGHTS(METHOD,ENERG,H1AO,SAOUNP,FAOS,DAOS,INDEVC,
     &                       NEVC,GVEC,BMAT,CVEC,MXERRV,KFREE,LFREE,WRK)
c
c     Computes weights that old density matrices are summed up with
c     to create new trial density matrix.
c     GVEC and BMAT store information between steps.
c     Exact content of these date structures depends on the method used
c     to optimize (EDIIS, total energy based, h1-energy based, 
c     DIIS, with lagrange multiplier, without, etc).
c
c     The methods are assumed to optimize a quadratic problem:
c        E = G*c + c'*B*c
c     where vector G and matrix B are specific to the method.
c
c     No matter which method is implemented, it can be divided into three
c     phases: a). update GVec. b). update B matrix. c). find coefficients c.
c
c     The convention is that the weights vector C includes an extra entry
c     corresponding to the quality of the solution. It can be lambda in
c     the case of DIIS, or interpolated energy (possibly the electronic
c     part, depending on the formulation of the problem) in the case of
c     EDIIS.
c
c     nitpicking: function pointers would be ideal here...
#include <implicit.h>
#include <maxorb.h>
#include <infinp.h>
#include <priunit.h>
#include <inforb.h>
      CHARACTER*10 METHOD
      DIMENSION H1AO(NNBAST),SAOUNP(NBAST,NBAST)
      DIMENSION FAOS(NNBAST,MXERRV,2),DAOS(NNBAST,MXERRV,2)
      INTEGER   INDEVC, NEVC
      DIMENSION GVEC(NEVC),BMAT(NEVC,NEVC),CVEC(NEVC+1)

      DIMENSION GVECTMP(NEVC),WRK(*)
      LOGICAL   FEASIBLE
      CALL QENTER('DIISWEIGHTS')
c     call outpkb(DAOS(1,INDEVC,2),NBAS,NSYM,1,0)
c     setup the default that acts only when NEVC == 1
      CVEC(1) = 1D0
      IF (METHOD(1:4).EQ.'DIIS') THEN
c        GVEC is not accessed for DIIS
         IF(NEVC.GT.1) THEN
            CALL diis_update_bmat(H1AO,SAOUNP,FAOS,DAOS,BMAT,NEVC,
     &           MXERRV,KFREE,LFREE,WRK)
            CALL diis_solve_for_coefs(NEVC,GVEC,BMAT,CVEC)
         ENDIF
      ELSE IF (METHOD(1:6).EQ.'C2DIIS') THEN
c        GVEC is not accessed for C2DIIS
         IF(NEVC.GT.1) THEN
            CALL diis_update_bmat(H1AO,SAOUNP,FAOS,DAOS,BMAT,NEVC,
     &           MXERRV,KFREE,LFREE,WRK)
            CALL c2diis_solve_for_coefs(NEVC,BMAT,CVEC,WRK,KFREE,LFREE)
         ENDIF

      ELSE IF (METHOD.EQ.'EDIIS') THEN
         CALL ediis_update_gvec(nevc,energ,H1AO,DAOS,GVEC,INDEVC)
         IF(NEVC.GT.1) THEN
            CALL DCOPY(NEVC,GVEC,1,GVECTMP,1)
            CALL ediis_update_bmat(NEVC,H1AO,SAOUNP,FAOS,DAOS,BMAT,
     &                             GVECTMP,INDEVC)
            CALL ediis_solve_for_coefs(NEVC,GVECTMP,BMAT,CVEC,INDEVC,
     &                                 FEASIBLE)
            IF(.NOT.FEASIBLE) THEN
               CALL diis_update_bmat(H1AO,SAOUNP,FAOS,DAOS,BMAT,NEVC,
     &        MXERRV,KFREE,LFREE,WRK)
               CALL diis_solve_for_coefs(NEVC,GVEC,BMAT,CVEC)
            END IF
         END IF
      ELSE
         WRITE(LUPRI,*) 'ERROR: unknown DIIS optimization method: ',
     &        METHOD
         CALL QUIT('Unknown DIIS optimization method.')
      END IF
      CALL QEXIT('DIISWEIGHTS')
      RETURN
      END
c     -------------------------------------------------------------------
c     DIIS-specific chunk of code
c     -------------------------------------------------------------------
      SUBROUTINE diis_update_bmat(H1AO,SAOUNP,FAOS,DAOS,BMAT,NEVC,
     &     MXERRV,KFRSAV,LFRSAV,WRK)
#include <implicit.h>
#include <inforb.h>
      DIMENSION H1AO(*),SAOUNP(*)
      DIMENSION FAOS(NNBAST,MXERRV,2),DAOS(NNBAST,MXERRV,2)
      DIMENSION BMAT(NEVC,NEVC),WRK(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)

c      call output(d,1,nbast,1,nbast,nbast,nbast,1,0) 
c      call outpkb(FAOS(1,INDEVC,1),NBAS,NSYM,1,0)
      KFREE=KFRSAV
      LFREE=LFRSAV

      CALL MEMGET('REAL',KTMP1 ,NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOMMi ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOMMj ,N2BASX,WRK,KFREE,LFREE)
      IF (NASHT.GT.0) THEN
         CALL MEMGET('REAL',KTMP2,N2BASX,WRK,KFREE,LFREE)
      ENDIF
      DO I = 1,NEVC
C     S*DV*FD-FD*DV*S
         CALL DCOPY(NNBAST,FAOS(1,I,1),1,WRK(KTMP1),1)
         CALL DAXPY(NNBAST,D1,H1AO,1,WRK(KTMP1),1)
         CALL SDFCOMM(WRK(KTMP1),DAOS(1,I,1),
     &        SAOUNP,WRK(KCOMMi),KFREE,LFREE,WRK)
C     S*DC*FV-FV*DC*S
         IF (NASHT.GT.0) THEN
            CALL DAXPY(NNBAST,DM1,FAOS(1,I,2),1,WRK(KTMP1),1)
            CALL SDFCOMM(WRK(KTMP1),DAOS(1,I,2),
     &           SAOUNP,WRK(KTMP2),KFREE,LFREE,WRK)
            CALL DAXPY(N2BASX,D1,WRK(KTMP2),1,WRK(KCOMMi),1)
         ENDIF
         DO J = 1,NEVC
         CALL DCOPY(NNBAST,FAOS(1,J,1),1,WRK(KTMP1),1)
         CALL DAXPY(NNBAST,D1,H1AO,1,WRK(KTMP1),1)
         CALL SDFCOMM(WRK(KTMP1),DAOS(1,J,1),
     &        SAOUNP,WRK(KCOMMj),KFREE,LFREE,WRK)
         IF (NASHT.GT.0) THEN
            CALL DAXPY(NNBAST,DM1,FAOS(1,J,2),1,WRK(KTMP1),1)
            CALL SDFCOMM(WRK(KTMP1),DAOS(1,J,2),
     &           SAOUNP,WRK(KTMP2),KFREE,LFREE,WRK)
            CALL DAXPY(N2BASX,D1,WRK(KTMP2),1,WRK(KCOMMj),1)
         ENDIF
         BMAT(I,J)=2*DDOT(N2BASX,WRK(KCOMMi),1,WRK(KCOMMj),1)
         ENDDO
      ENDDO
      CALL MEMREL('DIIS.FDCOMM',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
c     conditioning of the bmatrix to get better lagrangian multipliers.
      BDIV = SQRT(BMAT(1,1)*BMAT(NEVC,NEVC))
      DO I = 1,NEVC
        DO J = 1,NEVC
            BMAT(I,J) = BMAT(I,J)/BDIV
         ENDDO
      ENDDO
      END
c     -------------------------------------------------------------------
      SUBROUTINE diis_solve_for_coefs(NEVC,BMAT,CVEC)
c     find optimal weights for DIIS using precomputed B matrix
#include <implicit.h>
      DIMENSION BMAT(NEVC,NEVC),CVEC(NEVC+1)
      DIMENSION CONSTR(NEVC+1,NEVC+1)
      INTEGER PIV(NEVC+1)
#include <priunit.h>
      NPAR = NEVC+1
      DO I = 1,NEVC
         CONSTR(NPAR,I) = 1.D0
         CONSTR(I,NPAR) = 1.D0
         DO J = 1,NEVC
            CONSTR(J,I) = BMAT(J,I)
            CONSTR(I,J) = BMAT(J,I)
         ENDDO
      ENDDO
      CONSTR(NPAR,NPAR) = 0D0
      CALL DZERO(CVEC,NEVC)
      CVEC(NPAR) = 1.D0
      CALL DGESV(NPAR,1,CONSTR,NPAR,PIV,CVEC,NPAR,INFO)
      IF (INFO.NE.0) THEN
         WRITE(LUPRI,*)'ERROR: DIIS matrix cannot be inverted'
         CALL QUIT('ERROR: DIIS matrix cannot be inverted')
      END IF
      RETURN
      END

c     
c     -------------------------------------------------------------------
c     C2DIIS-specific chunk of code
c     -------------------------------------------------------------------
      SUBROUTINE c2diis_solve_for_coefs(NEVC,BMAT,CVEC,
     &       WRK,KFRSAV,LFRSAV)
c     find optimal weights for C2DIIS using precomputed B matrix
#include <implicit.h>
      DIMENSION BMAT(NEVC,NEVC),CVEC(NEVC+1),WRK(*)
      PARAMETER (THREVC = 0.01D0, D1=1.0D0, DP5=0.5D0)
      IROW(I) = I*(I-1)/2

      KFREE=KFRSAV
      LFREE=LFRSAV

      LBMAT=NEVC*(NEVC+1)*DP5
      CALL MEMGET('REAL',KBTMP,LBMAT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVEC,NEVC*NEVC,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KWJ  ,NEVC     ,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIWJ ,NEVC     ,WRK,KFREE,LFREE)

      
      CALL DUNIT(WRK(KEVEC),NEVC)
      CALL DGETSP(NEVC,BMAT,WRK(KBTMP))
      CALL JACO(WRK(KBTMP),WRK(KEVEC),NEVC,NEVC,NEVC,WRK(KWJ),WRK(KIWJ))

      DO I = 1,NEVC
         II = KBTMP-1+IROW(I+1)
         WRK(KBTMP-1+I) = WRK(II)
      ENDDO
      CALL ORDER (WRK(KEVEC),WRK(KBTMP),NEVC,NEVC)
      IOK = 0
      DO I = 1,NEVC
         EVCSUM = DSUM(NEVC,WRK(KEVEC+(I-1)*NEVC),1)
         IF (ABS(EVCSUM) .GT. THREVC) THEN
            IOK = I
            XLMBDA = WRK(KBTMP-1+I)
            CALL DCOPY(NEVC,WRK(KEVEC+(I-1)*NEVC),1,CVEC,1)
            CALL DSCAL(NEVC,(D1/EVCSUM),CVEC,1)
            GO TO 100
         ENDIF
      ENDDO
 100  CONTINUE
      CALL MEMREL('C2DIIS.coefs',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      END

c     
c     -------------------------------------------------------------------
c     EDIIS-specific chunk of code
c     -------------------------------------------------------------------
      SUBROUTINE ediis_update_gvec(NEVC,ENERGTOT,H1AO,DAOS,GVEC,INDEVC)
#include <implicit.h>
#include <inforb.h>
#include <priunit.h>
       DIMENSION H1AO(*),DAOS(NNBAST,NEVC,2),GVEC(*)
#ifdef EDIIS_WITH_DIFFERENCES_FORMULAS
      GVEC(INDEVC) = ENERGTOT
#else /* use other, possible more numerically stable formalism */
         GVEC(INDEVC) = DDOT(NNBAST,DAOS(1,INDEVC,1),1,H1AO,1)
#endif
      END
c     -------------------------------------------------------------------
      SUBROUTINE ediis_update_bmat(NEVC,H1AO,SAOUNP,FAOS,DAOS,BMAT,
     &                             GVEC,INDEVC)
#include <implicit.h>
#include <inforb.h>
#include <priunit.h>
      INTEGER NEVC,I,J
      DIMENSION H1AO(NNBAST),SAOUNP(NBAST,NBAST)
      DIMENSION FAOS(NNBAST,NEVC,2),DAOS(NNBAST,NEVC,2)
      DIMENSION GVEC(NEVC), BMAT(NEVC,NEVC)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)

      CALL DZERO(BMAT,NEVC*NEVC)
#ifdef EDIIS_WITH_DIFFERENCES_FORMULAS
      DO I = 1,NEVC-1
         DO J = I+1,NEVC
            CALL DCOPY(NNBAST,FAOS(1,J,1),1,FTMP,1)
            CALL DCOPY(NNBAST,DAOS(1,J,1),1,DTMP,1)
            CALL DAXPY(NNBAST,DM1,FAOS(1,I,1),1,FTMP,1)
            CALL DAXPY(NNBAST,DM1,DAOS(1,I,1),1,DTMP,1)
            BMAT(I,J) = DDOT(NNBAST,DTMP,1,FTMP,1)
            BMAT(J,I)=BMAT(I,J)
         ENDDO
      ENDDO

#else /* do EDIIS with formulas that do not substract similar matrices */
      ierr = 0
      DO I = 1,NEVC
         DO J = I,NEVC
            BMAT(I,J) = -DDOT(NNBAST,DAOS(1,I,1),1,FAOS(1,J,1),1)
#ifdef NO_BMAT_SYM_TEST
            BMAT(J,I) = BMAT(I,J)
#else
            BMAT(J,I) = -DDOT(NNBAST,DAOS(1,J,1),1,FAOS(1,I,1),1)
            IF (ABS(BMAT(J,I) - BMAT(I,J)) .gt. 1.0D-10) ierr = ierr + 1
#endif
         ENDDO
      ENDDO
      DO J = 1,NEVC
         IF (J .NE. INDEVC) THEN
            DO I = 1,NEVC
               IF (I .NE. INDEVC) THEN
                  BMAT(I,J) = BMAT(I,J)-BMAT(INDEVC,J)-BMAT(I,INDEVC)
     &                 + BMAT(INDEVC,INDEVC)
               END IF
            ENDDO
            GVEC(J) = GVEC(J) - GVEC(INDEVC)
     &           - BMAT(INDEVC,J) + BMAT(INDEVC,INDEVC)
         END IF
      ENDDO

      GVEC(INDEVC) = D0
      DO J = 1,NEVC
         BMAT(J,INDEVC) = D0
         BMAT(INDEVC,J) = D0
      ENDDO
      BMAT(INDEVC,INDEVC) = 99999.D0
#endif /* EDIIS_WITH_DIFFERENCES_FORMULAS */
      END 
c     -------------------------------------------------------------------
      SUBROUTINE ediis_solve_for_coefs(NEVC,GVEC,BMAT,WEIGHTVEC,
     &                                 INDEVC,FEASIBLE)
#include <implicit.h>
#include <inforb.h>      
#include <priunit.h>

      INTEGER NEVC
      DIMENSION GVEC(*),BMAT(NEVC,NEVC),WEIGHTVEC(*)
      INTEGER INDEVC
      LOGICAL FEASIBLE

C     local, temporary arrays
      DIMENSION CVEC(NEVC)
C     The B-matrix with the constraint sum(weights)=1.
      DIMENSION CONSTR(NEVC,NEVC)
      INTEGER PIV(NEVC)

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)
C
      IF (NEVC .EQ. 1) THEN
         CALL QUIT('EDIIS called with NEVC .eq. 1')
      END IF
c#define EDIIS_SCAN_WHOLE_SPACE 1
#ifdef EDIIS_SCAN_WHOLE_SPACE 
         CALL ediis_scan_whole_space(
     &           NEVC,GVEC,BMAT,WEIGHTVEC,INDEVC,FEASIBLE)
         GO TO 9999
#endif

      N2EVC = NEVC*NEVC
      CALL DCOPY(N2EVC,BMAT,1,CONSTR,1)
c        WRITE(LUPRI,*) ' CONSTR matrix:'
c        CALL OUTPUT(CONSTR,1,NEVC,1,NEVC,NEVC,NEVC,1,LUPRI)
c
C     Put the weights where none is forced zero last in CVEC.

      CALL DCOPY(NEVC,GVEC,1,CVEC,1)
c        WRITE(LUPRI,*) 'DIISenergies. ', NEVC, 'error vectors'
c        WRITE(LUPRI,'(8F10.5)') (CVEC(I),i=1,NEVC)      
      CALL DGESV(NEVC,1,CONSTR,NEVC,PIV,CVEC,NEVC,INFO)
      FEASIBLE = INFO.EQ.0

c     Check out quality of the initial EDIIS solution.
C     (reuse CONSTR as work array)
      CALL DGEMV('N',NEVC,NEVC,1.D0,BMAT,NEVC,cvec,1, 0.D0,CONSTR,1)
      EK1 = DDOT(NEVC,CVEC,1,GVEC,1)
      EK2 = - 0.5D0*DDOT(NEVC,CVEC,1,CONSTR,1)
      EK  = EK1 + EK2
      SUMC = DSUM(NEVC,CVEC,1)
      CVEC(INDEVC) = D1 - SUMC
      DO I = 1,NEVC
         IF (CVEC(I).LT.D0) THEN
            FEASIBLE = .FALSE.
         ENDIF
      ENDDO
      WRITE(LUPRI,*) 'Computed EDIIS weights:'
      WRITE(LUPRI,'(8F10.5)') (CVEC(I),i=1,NEVC)      
      WRITE(LUPRI,*) '(full) interpolated Delta energy: ',EK, FEASIBLE
c
c     look for other, constrained solution when using ediis.
c
      IF (.NOT.FEASIBLE) THEN
c        ELSE No solution lowering the energy found
      ELSE
         CALL DCOPY(NEVC,CVEC,1,WEIGHTVEC,1)
      ENDIF
C
 9999 CONTINUE
      write(LUPRI,*) 'Final "ediis" weights, feasible=', FEASIBLE
      write(LUPRI,*) (weightvec(i),i=1,nevc)
      write(0,*) 'Final "ediis" weights, feasible=', FEASIBLE
      write(0,*) (weightvec(i),i=1,nevc)
      RETURN
      END
      SUBROUTINE ediis_scan_whole_space(
     &           NEVC,GVEC,BMAT,WEIGHTVEC,INDEVC,FEASIBLE)
#include <implicit.h>
#include <inforb.h>      
#include <priunit.h>

      INTEGER NEVC
      DIMENSION GVEC(*),BMAT(NEVC,NEVC),WEIGHTVEC(*)
      INTEGER INDEVC
      LOGICAL FEASIBLE

C     All the 2**NEVC-1 possibilities of a weightvector including energy.
      DIMENSION CVEC(NEVC+1,2**NEVC-1)
      DIMENSION CVECTMP(NEVC+1)
C     The B-matrix with the constraint sum(weights)=1.
      DIMENSION CONSTR(NEVC,NEVC)
      INTEGER I,PIV(NEVC),ACTSIZE
      LOGICAL ZEROWEIGHT(NEVC)
      INTEGER ZERCNT(2**NEVC-1)
      CHARACTER*2 CTYPE

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)
C     size (nr of possibilties) in active set
      ACTSIZE=2**NEVC-1
      N2EVC = NEVC*NEVC
      CALL DCOPY(N2EVC,BMAT,1,CONSTR,1)

c      WRITE(LUPRI,*) ' CONSTR matrix:'
c      CALL OUTPUT(CONSTR,1,NEVC,1,NEVC,NEVC,NEVC,1,LUPRI)
c
C     Put the weights where none is forced zero last in CVEC.

      CALL DCOPY(NEVC,GVEC,1,CVEC(1,ACTSIZE),1)
c      WRITE(LUPRI,*) 'DIISenergies. ', NEVC, 'error vectors'
c      WRITE(LUPRI,'(8F10.5)') (CVEC(I,ACTSIZE),i=1,NEVC)      
      CALL DGESV(NEVC,1,CONSTR,NEVC,PIV,CVEC(1,ACTSIZE),NEVC,INFO)
      FEASIBLE = INFO.EQ.0

c     Check out quality of the initial EDIIS solution.
      CALL DGEMV('N',NEVC,NEVC,1.D0,BMAT,NEVC,cvec(1,ACTSIZE),
     &     1, 0.D0,CVECTMP,1)
      EK1 = DDOT(NEVC,CVEC(1,ACTSIZE),1,GVEC,1)
      EK2 = - 0.5D0*DDOT(NEVC,CVEC(1,ACTSIZE),1,CVECTMP,1)
      EK = EK1 + EK2
      CVEC(NEVC+1,ACTSIZE) = EK
      SUMC = DSUM(NEVC,CVEC(1,ACTSIZE),1)
      CVEC(INDEVC,ACTSIZE) = D1 - SUMC
      DO I = 1,NEVC
         IF (CVEC(I,ACTSIZE).LT.D0) THEN
            FEASIBLE = .FALSE.
         ENDIF
      ENDDO
      WRITE(LUPRI,*) 'Computed EDIIS weights:'
      WRITE(LUPRI,'(8F10.5)') (CVEC(I,ACTSIZE),i=1,NEVC)      
      WRITE(LUPRI,*) '(full) interpolated Delta energy: ',EK, FEASIBLE
c
c     look for other, constrained solution when using ediis.
c
      IF (.NOT.FEASIBLE) THEN
         CALL DZERO(CONSTR,N2EVC)
         DO K = 1,ACTSIZE-1
            NROFZEROS = 0
            DO I = 1,NEVC
               IF (BTEST(K,I-1)) THEN
                  ZEROWEIGHT(I) = .TRUE.
                  NROFZEROS = NROFZEROS + 1
               ELSE
                  ZEROWEIGHT(I) = .FALSE.
               ENDIF
            ENDDO
            ZERCNT(K) = NROFZEROS
            IOFF = 0
            DO I = 1,NEVC-NROFZEROS
               IOFF = IOFF + 1
               DO WHILE (BTEST(K,IOFF-1).AND.IOFF.LE.NEVC)
                  IOFF = IOFF + 1
               ENDDO
               JOFF = 0
               DO J = 1,NEVC-1-NROFZEROS
                  JOFF = JOFF + 1
                  DO WHILE (BTEST(K,JOFF-1).AND.JOFF.LE.NEVC)
                     JOFF = JOFF + 1
                  ENDDO
                  CONSTR(I,J)=BMAT(IOFF,JOFF)
                  CONSTR(J,I)=CONSTR(I,J)
               ENDDO
               CVEC(I,K)=GVEC(IOFF)
            ENDDO
            NPAR = NEVC-NROFZEROS
c            WRITE(LUPRI,*) 'CONSTR with NROFZEROS=',NROFZEROS
c            CALL OUTPUT(CONSTR,1,NPAR,1,NPAR,NEVC,NEVC,1,LUPRI)
            CALL DGESV(NPAR,1,CONSTR,NEVC,PIV,CVEC(1,ACTSIZE),
     &                 NEVC,INFO)
            IF(INFO.NE.0) THEN
c     Failure can be perfectly normal due to the way we construct
c     the B matrix. We just have to skip this combination.
               FEASIBLE = .FALSE.
            ELSE
               FEASIBLE = .TRUE.
            ENDIF

            IOFF=0
            DO I = 1,NEVC
               IF (ZEROWEIGHT(I)) THEN
                  CVECTMP(I) = D0
                  IOFF = IOFF + 1
               ELSE
                  IF (CVEC(I-IOFF,K).LT.D0) THEN
                     FEASIBLE = .FALSE.
                  ENDIF
                  CVECTMP(I) = CVEC(I-IOFF,K)
               ENDIF
            ENDDO
            CALL DCOPY(NEVC,CVECTMP,1,CVEC(1,K),1)
            CALL DGEMV('N',NEVC,NEVC,1D0,BMAT,NEVC,cvec(1,k),
     &                 1, 0D0,CVECTMP,1)
            EK = DDOT(NEVC,CVEC(1,k),1,GVEC,1)
     &           - 0.5*DDOT(NEVC,CVEC(1,k),1,CVECTMP,1)
            SUMC = DSUM(NEVC,CVEC(1,K),1)
            CVEC(INDEVC,K) = D1 - SUMC
            FEASIBLE = FEASIBLE.AND.CVEC(INDEVC,K).GT.0D0
            IF (FEASIBLE) THEN
               CTYPE = '  '
               cvec(nevc+1,k) = EK
            ELSE
               CTYPE = 'UN'
               CVEC(NEVC+1,k) = 1.D0  !9999.D0
            ENDIF
            WRITE(LUPRI,*) CTYPE,'feasible solution. Interpolated E=',
     &           ek,' weights:'
            WRITE(LUPRI,'(5F10.4)') (cvec(i,k),i=1,nevc)
         ENDDO
         WLAMBDA = 10.0D0
         J= 99 999 999
         DO I = 1,ACTSIZE-1
c            IF (ZERCNT(I).LT.2.AND.CVEC(NEVC+1,I).LT.WLAMBD) THEN
            IF (CVEC(NEVC+1,I).LT.WLAMBD) THEN
C     p 8259 (if most recent fockmatrix' weight is zero should be handled
C     different than just skipping it
               IF (CVEC(INDEVC,I).NE.D0) THEN
                  WLAMBD = CVEC(NEVC+1,I)
                  J=I
               ENDIF
            ENDIF
         ENDDO
         IF(J.LE.ACTSIZE) THEN
            CALL DCOPY(NEVC,CVEC(1,J),1,WEIGHTVEC,1)
         ELSE
            FEASIBLE = .FALSE.
         ENDIF
c        ELSE No solution lowering the energy found
      ELSE
         CALL DCOPY(NEVC,CVEC(1,ACTSIZE),1,WEIGHTVEC,1)
      ENDIF
      write(LUPRI,*) 'Final weights, feasible=', FEASIBLE
      write(LUPRI,*) (weightvec(i),i=1,nevc)
      print *, 'Final weights, feasible=', FEASIBLE
      print *, (weightvec(i),i=1,nevc)
      END
C --- end of opt-solvers.F ---
