C
C...   Copyright (c) 2015 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program,
C...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C...   For further information, including how to get a licence, see:
C...      http://daltonprogram.org
C
C
C ============================================================================
C     newabs.F: Main author Joanna Kauczor
C     This implementation is published in:
C
C        Reference C:
C        J. Kauczor, P. Jorgensen, and P. Norman, 
C        "On the efficiency of algorithms for solving Hartree-Fock and Kohn-Sham
C        response equations",
C        J. Chem. Theory Comput. 7 (2011) 1610.
C ============================================================================
C
      SUBROUTINE ABS_CTL(LABEL,KZVAR,GD,NGD,XSOL,FREQS,NFREQS,LU,
     &                   MJWOP,CMO,UDV,FC,FCAC,FV,PV,XINDX,
     &                   RESLRF,WRK,LWRK)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
#include "abscrs.h"
C
C PURPOSE:
C
C     Solve the complex LR equation 
C
C     ( E[2] - {w+iW}*S[2] )* (NR + iNI) = B[1]
C
C     using the new rsp solver with symmetrized trial vectors
C
      CHARACTER LABEL*8
      DIMENSION XSOL(2*KZVAR,2,NFREQS),GD(4*KZVAR,NGD)
      DIMENSION MJWOP(*),XINDX(*)
      DIMENSION CMO(*),UDV(*),FC(*),FCAC(*),FV(*),PV(*)
      DIMENSION FREQS(NFREQS)
      DIMENSION RESLRF(2,NFREQS,3,3,4),WRK(LWRK)
      DIMENSION KLRED(2),KNVEC(2)
      LOGICAL CONV(NFREQS),CONVERGED,EX1,EX2,EX3
C
C
      NEQ2S=NFREQS
      KZYVAR=2*KZVAR
      KGD    = 1
      KXTRL  = KGD   + 4*NGD*KZVAR
      KXRED  = KXTRL + 4*NEQ2S*KZVAR
      KRESD  = KXRED + 4*NEQ2S*ABS_MAXRM
      KREDE  = KRESD + 4*NEQ2S*KZVAR
      KREDS  = KREDE  + 2*ABS_MAXRM*ABS_MAXRM
      KREDGD = KREDS  + ABS_MAXRM*ABS_MAXRM
      KFREE  = KREDGD + 4*ABS_MAXRM*NGD
      LFREE  = LWRK - KFREE

      IF (LFREE.LT.0) STOP 'ABS_CTL: out of memory'
c
C      IF (.NOT. ABS_IDRI_XY) THEN
C      CALL ABS_CHKONFILE(LU,CONVERGED,LABEL,ABS_GRADSYM,
C     &     ABS_DAMP,NEQ2S,FREQS,ABS_THCLR,WRK(KFREE))
C      IF (CONVERGED) THEN
C         WRITE(LUABSPRI,'(/A,I3)')
C     &        ' ABSCTL: converged response vectors found on file'//
C     &        ' RSPVEC with unit=',LUABSVECS
C         CALL GET_LRF(LU,LABEL,NEQ2S,FREQS,FC,FV,CMO,UDV,PV,
C     &               XINDX,RESLRF,KZVAR,WRK(KFREE),LFREE)
C         GOTO 900
C      END IF
C      END IF
c
      WRK(KGD:KGD+4*NGD*KZVAR-1) = 0.0d0
      WRK(KXRED:KXRED+4*NEQ2S*ABS_MAXRM-1) = 0.0d0
      WRK(KRESD:KRESD+4*NEQ2S*KZVAR-1) = 0.0d0
c
      DO I=1,NEQ2S
        CONV(I)=.FALSE.
      ENDDO
C
      IF ((MAX(ABS_KLRED(1),ABS_KLRED(2)).NE.0)
     &    .AND. .NOT. ABS_IDRI_XY.AND. .NOT. ABS_NGD
     &    .AND. .NOT. ABS_GDCOMPLEX) THEN
C     NHL 04/16: Needs to be checked for ABS_NGD=TRUE and
C               accomodated for complex rhs.
c        write(luabspri,*) 'ABS_KLRED >0, restart'
       DO J=1,NGD
        DO I=1,KZVAR
          WRK(KGD-1+(J-1)*4*KZVAR+I)=0.5d0*(GD(I,J)+GD(KZVAR+I,J))
          WRK(KGD+KZVAR-1+(J-1)*4*KZVAR+I)=0.5d0*(GD(I,J)-GD(KZVAR+I,J))
C         IF (ABS_GDCOMPLEX) THEN
C           take care of the imaginary part  
C         END IF
        END DO
       END DO
        CALL GPINQ('ABS_E1RED','EXIST',EX1)
        CALL GPINQ('ABS_E2RED','EXIST',EX2)
        CALL GPINQ('ABS_SRED','EXIST',EX3)
        IF (EX1.AND.EX2.AND.EX3) THEN
          write(luabspri,*)
     &          'Use ABS_READIN_RED ',ABS_KLRED(1),ABS_KLRED(2)
          CALL ABS_READIN_RED(KZVAR,WRK(KREDE),WRK(KREDS),WRK(KGD),
     &                        WRK(KREDGD),KNVEC,WRK(KFREE),LFREE)
        ELSE
          write(luabspri,*)
     &    'Use ABS_REDSPACE_REBUILD',ABS_KLRED(1),ABS_KLRED(2),
     &                               EX1,EX2,EX3
          CALL ABS_REDSPACE_REBUILD(KZVAR,WRK(KGD),WRK(KREDE),
     &                            WRK(KREDS),WRK(KREDGD),KNVEC,
     &                            WRK(KFREE),LFREE)
        ENDIF
        KLRED(1)=ABS_KLRED(1)
        KLRED(2)=ABS_KLRED(2)
        NNLT=KNVEC(1)+KNVEC(2)
        IF (NNLT .GT.0) THEN
          CALL ABS_GETRIAL(KZVAR,GD,WRK(KGD),NGD,FREQS,NEQ2S,WRK(KXTRL),
     &               WRK(KRESD),
     &               KLRED,KNVEC,NNLT,CONV,WRK(KFREE),LFREE)
        ENDIF
      ELSE
c       CALL DZERO(WRK(KRESD),4*NEQ2S*KZVAR)
c       CALL DZERO(WRK(KREDE),3*ABS_MAXRM*ABS_MAXRM)
c       CALL DZERO(WRK(KREDGD),2*ABS_MAXRM)
        WRK(KREDE:KREDE+3*ABS_MAXRM*ABS_MAXRM-1) = 0.0d0
        WRK(KREDGD:KREDGD+4*ABS_MAXRM*NGD-1) = 0.0d0
c
        KNVEC(1)=2*NEQ2S
        KNVEC(2)=2*NEQ2S
        NNLT=4*NEQ2S
        KLRED(1)=0
        KLRED(2)=0
        IF ((.NOT. ABS_REBD) .AND. (.NOT. ABS_IDRI_XY)) THEN
          REWIND(LUSB)
          REWIND(LUAB)
          REWIND(LUSS)
          REWIND(LUAS)
        ENDIF
        
        CALL ABS_GETRIAL(KZVAR,GD,WRK(KGD),NGD,FREQS,NEQ2S,WRK(KXTRL),
     &               WRK(KRESD),KLRED,KNVEC,NNLT,CONV,WRK(KFREE),LFREE)
c
      ENDIF
c
c       CALL DZERO(WRK(KRESD),4*NEQ2S*KZVAR)
c      WRK(KRESD:KRESD+4*NEQ2S*KZVAR-1) = 0.0d0
      WRK(KRESD:KRESD+4*NEQ2S*KZVAR-1) = 0.0d0
      
      ITER = 1
 100  CONTINUE
C     
      IF (IPRABSLRS.GE.2) THEN
       WRITE(LUABSPRI,'(/A,/A,I3,A,/A)')
     &        ' =======================================',
     &        ' ----  Macro iteration number',ITER,'  ------',
     &        ' ======================================='
      END IF
C
C     Perform linear transformation and build reduce spaces
C
      CALL ABS_RED(ITER,KZVAR,WRK(KGD),WRK(KXTRL),WRK(KRESD),
     &             WRK(KREDGD),WRK(KREDE),WRK(KREDS),NEQ2S,NNLT,KLRED,
     &             KNVEC,NGD,MJWOP,CMO,UDV,FC,FCAC,FV,PV,XINDX,
     &             WRK(KFREE),LFREE)
C
C     Solve RPA equation in a reduced space
C
      KZLRED=2*KLRED(1)+2*KLRED(2)
      IF (ABS_IMFREQ.AND. .NOT. ABS_NGD .AND. .NOT. ABS_GDCOMPLEX) THEN
        CALL ABS_IMGSOLV(FREQS,NEQ2S,WRK(KREDE),WRK(KREDS),WRK(KREDGD),
     &                 WRK(KXRED),KLRED,KZLRED,CONV,WRK(KFREE),LFREE)
      ELSE
        CALL ABS_SOLV(FREQS,NEQ2S,WRK(KREDE),WRK(KREDS),WRK(KREDGD),
     &                WRK(KXRED),NGD,KLRED,KZLRED,CONV,WRK(KFREE),LFREE)
      ENDIF
C
C     Compute the residual in this iteration
C     
      WRK(KRESD:KRESD+4*NEQ2S*KZVAR-1) = 0.0d0
c       call DZERO(WRK(KRESD),4*NEQ2S*KZVAR)
      CALL ABS_RSD(WRK(KXRED),LABEL,KZVAR,NEQ2S,FREQS,KLRED,
     &             KNVEC,NNLT,NGD,WRK(KGD),
     &             WRK(KRESD),ITER,CONVERGED,
     &             CONV,WRK(KXTRL),XSOL,LU,WRK(KFREE),LFREE)
C   
C     Calculate polarizability
C
      DO J=1,NEQ2S
        IF (LABEL(2:8).EQ.'DIPLEN') THEN
          ITYPE = 1
        ELSE IF (LABEL(2:8).EQ.'ANGMOM') THEN
          ITYPE = 2
        ELSE IF (LABEL(3:7).EQ.'THETA') THEN
          ITYPE = 3
        ELSE IF (LABEL(2:8).EQ.'DIPVEL') THEN
          ITYPE = 4
        ELSE
!         WRITE(LUABSPRI,'(2A)') ' Warning : Unknown operator! ',LABEL
!         - operator might be PSO (for .NSCD)
          GO TO 778
        END IF
C
        CALL DIPLAB(LABEL,INDEX)
        RESLRF(1,J,INDEX,INDEX,ITYPE) = 
     &         DDOT(KZYVAR,GD,1,XSOL(1,1,J),1)
        RESLRF(2,J,INDEX,INDEX,ITYPE) = 
     &         DDOT(KZYVAR,GD,1,XSOL(1,2,J),1)
C     
      ENDDO
 778  CONTINUE
c
C     Check status
C     
      IF (CONVERGED) GOTO 900
       
      IF (ITER.GE.ABS_MAXITER) GOTO 910
C     
C     Iteration completed
C     
      ITER = ITER + 1
      GOTO 100
C     
C     Final messages from ABSCTL
C     
 900  CONTINUE
      IF (IPRABSLRS.GE.0) THEN
         WRITE(LUABSPRI,'(/,1X,A,I4,A,I4,I4)')
     &        '*** ABSCTL: REQUESTED',NFREQS,
     &        ' SOLUTION VECTORS CONVERGED IN RED SPACES OF',KLRED(1),
     &        KLRED(2)
      END IF
      GOTO 777
 910  CONTINUE
      IF (IPRABSLRS.GE.0) THEN
         WRITE(LUABSPRI,'(/,1X,A,I4,A,I4,I4)')
     &        '*** ABSCTL: REQUESTED',NFREQS,
     &        ' SOLUTION VECTORS NOT CONVERGED IN RED SPACES OF',
     &        KLRED(1),KLRED(2)
         WRITE(LUABSPRI,'(/,1X,A)')
     &        '--- MAXIMUM NUMBER OF ITERATIONS REACHED ---'
      END IF
      GOTO 777
C
 777  CONTINUE
C
      IF (ABS_REBD) THEN 
        ABS_KLRED(1)=KLRED(1)
        ABS_KLRED(2)=KLRED(2)
C
        CALL GPINQ('ABS_E1RED','EXIST',EX1)
        CALL GPINQ('ABS_E2RED','EXIST',EX2)
        CALL GPINQ('ABS_SRED','EXIST',EX3)
c
        IF (EX1.AND.EX2.AND.EX3) THEN
          REWIND(LUE1RED)
          REWIND(LUE2RED)
          REWIND(LUSRED)
          DO I=1,KLRED(1)
            CALL WRITE_VEC(LUE1RED,KLRED(1),WRK(KREDE+(I-1)*ABS_MAXRM))
          ENDDO
          DO I=1,KLRED(2)
            CALL WRITE_VEC(LUE2RED,KLRED(2),WRK(KREDE+(I-1)*ABS_MAXRM+
     &                                           ABS_MAXRM*ABS_MAXRM))
            CALL WRITE_VEC(LUSRED,KLRED(1),WRK(KREDS+(I-1)*ABS_MAXRM))
          ENDDO
        ENDIF
       ELSE
         ABS_KLRED(1)=0
         ABS_KLRED(2)=0
       ENDIF
C
      END
      SUBROUTINE ABS_SOLV(FREQS,NFREQS,REDE,REDS,REDGD,RDX,NGD,KLRED,
     &                    KZLRED,CONV,WRK,LWRK)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
C
C PURPOSE:
C     Solve the equation in a reduced space 4nx4n
C
      INTEGER NFREQS,KZLRED, NGD
      DIMENSION REDE(ABS_MAXRM,ABS_MAXRM,2),REDS(ABS_MAXRM,ABS_MAXRM)
      DIMENSION REDGD(ABS_MAXRM,4,NGD)
      DIMENSION RDX(ABS_MAXRM,4,NFREQS),FREQS(NFREQS),WRK(LWRK)
      INTEGER KLRED(2),KIPIV(KZLRED)
      LOGICAL CONV(NFREQS)
C     local
      INTEGER I,J,K,INFO,NOS,NOA
c      DOUBLE PRECISION FREQI
C
      INFO=0
      KREDZGD = 1
      KREDZ   = KREDZGD + KZLRED
      KFREE   = KREDZ   + KZLRED*KZLRED
      LFREE   = LWRK    - KFREE
      IF (LFREE.LT.0) STOP 'ABS_SOLV: out of memory'
C      
      NOS=KLRED(1)
      NOA=KLRED(2)
C
      DO K=1,NFREQS
        IF (ABS_NGD) THEN
           IGD=K
        ELSE 
           IGD=1
        END IF
        IF (.NOT.CONV(K)) THEN
           WRK(KREDZ:KREDZ+KZLRED*KZLRED-1)=0.0d0
c
c          dgesv overwrites REDGD, therefore reduced RHS must be constructed 
C          for each frequency
c
           WRK(KREDZGD:KREDZGD+KZLRED-1)=0.0d0
           CALL DCOPY(NOS,REDGD(1,1,IGD),1,WRK(KREDZGD),1)
           CALL DCOPY(NOA,REDGD(1,2,IGD),1,WRK(KREDZGD+NOS),1)
           IF (ABS_GDCOMPLEX) THEN
              CALL DCOPY(NOA,REDGD(1,3,IGD),1,WRK(KREDZGD+NOS+NOA),1)
              CALL DCOPY(NOS,REDGD(1,4,IGD),1,WRK(KREDZGD+2*NOA+NOS),1)
           END IF
           DO I=1,KZLRED
             KIPIV(I)=0
           ENDDO
           FREQI=FREQS(K)
           DO I=1,NOS
C               E1
             DO J=1,NOS
                WRK(KREDZ-1+(I-1)*KZLRED+J)=REDE(I,J,1)
             ENDDO
C              -w*S1     
             DO J=1,NOA
                WRK(KREDZ-1+NOS*KZLRED+(J-1)*KZLRED+I)=
     &          -FREQI*REDS(I,J)
             ENDDO
C               g*S2
             DO J=1,NOA
                WRK(KREDZ-1+(NOS+NOA)*KZLRED+(J-1)*KZLRED+I)=
     &          ABS_DAMP*REDS(I,J)
             ENDDO
           ENDDO
           DO I=1,NOA
C               E2
             DO J=1,NOA
                WRK(KREDZ-1+NOS*(KZLRED+1)+(I-1)*KZLRED+J)=REDE(I,J,2)
             ENDDO
C               -w*S1T
             DO J=1,NOS
                WRK(KREDZ-1+NOS+(J-1)*KZLRED+I)=-FREQI*REDS(J,I)
             ENDDO
C               g*S3
             DO J=1,NOS
                WRK(KREDZ-1+(NOS+2*NOA)*KZLRED+NOS+(J-1)*KZLRED+I)
     &          =ABS_DAMP*REDS(J,I)
             ENDDO
           ENDDO
           DO I=1,NOA
C              -E3
             DO J=1,NOA
                WRK(KREDZ-1+(NOS+NOA)*(KZLRED+1)+(I-1)*KZLRED+J)=
     &          -REDE(I,J,2)
             ENDDO
C               g*S2T
             DO J=1,NOS
                WRK(KREDZ-1+(NOS+NOA)+(J-1)*KZLRED+I)=
     &          ABS_DAMP*REDS(J,I)
             ENDDO
C               w*S4
             DO J=1,NOS
                WRK(KREDZ-1+(NOS+2*NOA)*KZLRED+(NOS+NOA)
     &          +(J-1)*KZLRED+I)=FREQI*REDS(J,I)
             ENDDO
           ENDDO
           DO I=1,NOS
C              -E4
             DO J=1,NOS
                WRK(KREDZ-1+(NOS+2*NOA)*(KZLRED+1)+(I-1)*KZLRED+J)=
     &          -REDE(I,J,1)
             ENDDO
C               g*S3T
             DO J=1,NOA
                WRK(KREDZ-1+NOS*(KZLRED)+(NOS+2*NOA)+
     &          (J-1)*KZLRED+I)=ABS_DAMP*REDS(I,J)
             ENDDO
C               w*S4T
             DO J=1,NOA
                WRK(KREDZ-1+(NOS+NOA)*KZLRED+(NOS+2*NOA)+
     &          (J-1)*KZLRED+I)=FREQI*REDS(I,J)
             ENDDO
           ENDDO
C
           IF (IPRABSLRS.GT.10) THEN
               WRITE(LUABSPRI,'(/,5X,A,/,5X,2(A,D11.4),/,5X,A)') 
     &             ' Reduced ( E[2] - {w+ig}*S[2] )  matrix',
     &             ' with w =', FREQI,' and g =', ABS_DAMP,
     &             '========================================' 
               CALL OUTPUT(WRK(KREDZ),1,KZLRED,1,KZLRED,KZLRED,KZLRED,
     &              1,LUABSPRI)
               write(luabspri,'(2(/,5X,A))')'Reduced B[1]',
     &            '=========================================='
               CALL OUTPUT(WRK(KREDZGD),1,KZLRED,1,1,KZLRED,1,1,
     &              LUABSPRI)
          END IF
C
C          Solve set of linear equations Ax = b:
C
          CALL DGESV(KZLRED, 1, WRK(KREDZ), KZLRED,KIPIV, 
     &                       WRK(KREDZGD),KZLRED, INFO)
           IF (INFO /= 0) THEN
               WRITE(LUABSPRI,'(/A, i4)') 
     &            'Problem in DGESV, IERR = ', INFO
               CALL QUIT(' Problem in DGESV')
           ENDIF
C
C          Solution vector is found in RHS.
C
           CALL DCOPY(NOS,WRK(KREDZGD),1,RDX(1,1,K),1)
           CALL DCOPY(NOA,WRK(KREDZGD+NOS),1,RDX(1,2,K),1)
           CALL DCOPY(NOA,WRK(KREDZGD+NOS+NOA),1,RDX(1,3,K),1)
           CALL DCOPY(NOS,WRK(KREDZGD+NOS+2*NOA),1,RDX(1,4,K),1)
C
           IF (IPRABSLRS.GE.5) THEN
              WRITE(LUABSPRI,'(/,5X,A,/,5X,2(A,D11.4),/,5X,A)') 
     &           ' Reduced ( E[2] - {w+ig}*S[2] )-1 * B[1] vector',
     &           ' with w =', FREQI,' and g =', ABS_DAMP,
     &           '================================================' 
              CALL OUTPUT(RDX(1,1,K),1,NOS,1,1,NOS,1,1,LUABSPRI)
              CALL OUTPUT(RDX(1,2,K),1,NOA,1,1,NOA,1,1,LUABSPRI)
              CALL OUTPUT(RDX(1,3,K),1,NOA,1,1,NOA,1,1,LUABSPRI)
              CALL OUTPUT(RDX(1,4,K),1,NOS,1,1,NOS,1,1,LUABSPRI)
           END IF
        END IF
      END DO
C
      RETURN
      END
      SUBROUTINE ABS_IMGSOLV(FREQS,NFREQS,REDE,REDS,REDGD,RDX,KLRED,
     &                    KZLRED,CONV,WRK,LWRK)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
C
C PURPOSE:
C     Solve the equation for purely imaginary freq in a reduced space 4nx4n
C
      INTEGER NFREQS,KZLRED
      DIMENSION REDE(ABS_MAXRM,ABS_MAXRM,2),REDS(ABS_MAXRM,ABS_MAXRM)
      DIMENSION REDGD(ABS_MAXRM,2)
      DIMENSION RDX(ABS_MAXRM,4,NFREQS),FREQS(NFREQS),WRK(LWRK)
      INTEGER KLRED(2),KIPIV(KZLRED)
      LOGICAL CONV(NFREQS)
C     local
      INTEGER I,J,K,NOS,NOA,INFO
c      DOUBLE PRECISION FREQI 
C     
      !initialize INFO to prevent problems with DGESV on 64bit
      INFO=0
      KREDZGD = 1
      KREDZ   = KREDZGD + KZLRED
      KFREE   = KREDZ   + KZLRED*KZLRED
      LFREE   = LWRK    - KFREE
      IF (LFREE.LT.0) STOP 'ABS_IMAGSOLV: out of memory'
C     
      NOS=KLRED(1)
      NOA=KLRED(2)
C
      DO K=1,NFREQS
        IF (.NOT.CONV(K)) THEN
           WRK(KREDZ:KREDZ+KZLRED*KZLRED-1)=0.0d0
c
c          dgesv overwrites REDGD, therefore reduced RHS must be constructed 
C          for each frequency
c
           WRK(KREDZGD:KREDZGD+KZLRED-1)=0.0d0
           CALL DCOPY(NOS,REDGD(1,1),1,WRK(KREDZGD),1)
           CALL DCOPY(NOA,REDGD(1,2),1,WRK(KREDZGD+NOS),1)
           DO I=1,KZLRED
             KIPIV(I)=0
           ENDDO
           FREQI=FREQS(K)
           DO I=1,NOS
C               E1
             DO J=1,NOS
                WRK(KREDZ-1+(I-1)*KZLRED+J)=REDE(I,J,1)
             ENDDO
C               g*S2
             DO J=1,NOA
                WRK(KREDZ-1+(NOS+NOA)*KZLRED+(J-1)*KZLRED+I)=
     &          FREQI*REDS(I,J)
             ENDDO
           ENDDO
           DO I=1,NOA
C               E2
             DO J=1,NOA
                WRK(KREDZ-1+NOS*(KZLRED+1)+(I-1)*KZLRED+J)=REDE(I,J,2)
             ENDDO
C               g*S3
             DO J=1,NOS
                WRK(KREDZ-1+(NOS+2*NOA)*KZLRED+NOS+(J-1)*KZLRED+I)
     &          =FREQI*REDS(J,I)
             ENDDO
           ENDDO
           DO I=1,NOA
C              -E3
             DO J=1,NOA
                WRK(KREDZ-1+(NOS+NOA)*(KZLRED+1)+(I-1)*KZLRED+J)=
     &          -REDE(I,J,2)
             ENDDO
C               g*S2T
             DO J=1,NOS
                WRK(KREDZ-1+(NOS+NOA)+(J-1)*KZLRED+I)=
     &          FREQI*REDS(J,I)
             ENDDO
           ENDDO
           DO I=1,NOS
C              -E4
             DO J=1,NOS
                WRK(KREDZ-1+(NOS+2*NOA)*(KZLRED+1)+(I-1)*KZLRED+J)=
     &          -REDE(I,J,1)
             ENDDO
C               g*S3T
             DO J=1,NOA
                WRK(KREDZ-1+NOS*(KZLRED)+(NOS+2*NOA)+
     &          (J-1)*KZLRED+I)=FREQI*REDS(I,J)
             ENDDO
           ENDDO
C
           IF (IPRABSLRS.GT.10) THEN
               WRITE(LUABSPRI,'(/,5X,A,/,5X,2(A,D11.4),/,5X,A)') 
     &             ' Reduced ( E[2] - {w+ig}*S[2] )  matrix',
     &             ' with w =', FREQI,' and g =', 0.d0,
     &             '========================================' 
               CALL OUTPUT(WRK(KREDZ),1,KZLRED,1,KZLRED,KZLRED,KZLRED,
     &              1,LUABSPRI)
               write(luabspri,'(2(/,5X,A))')'Reduced B[1]',
     &            '=========================================='
               CALL OUTPUT(WRK(KREDZGD),1,KZLRED,1,1,KZLRED,1,1,
     &              LUABSPRI)
          END IF
C
C          Solve set of linear equations Ax = b:
C
          CALL DGESV(KZLRED, 1, WRK(KREDZ), KZLRED,KIPIV, 
     &                       WRK(KREDZGD),KZLRED, INFO)
           IF (INFO /= 0) THEN
               WRITE(LUABSPRI,'(/A, i4)') 
     &            'Problem in DGESV, IERR = ', INFO
               CALL QUIT(' Problem in DGESV')
           ENDIF
C
C          Solution vector is found in RHS.
C
           CALL DCOPY(NOS,WRK(KREDZGD),1,RDX(1,1,K),1)
           CALL DCOPY(NOA,WRK(KREDZGD+NOS),1,RDX(1,2,K),1)
           CALL DCOPY(NOA,WRK(KREDZGD+NOS+NOA),1,RDX(1,3,K),1)
           CALL DCOPY(NOS,WRK(KREDZGD+NOS+2*NOA),1,RDX(1,4,K),1)
C
           IF (IPRABSLRS.GE.5) THEN
              WRITE(LUABSPRI,'(/,5X,A,/,5X,2(A,D11.4),/,5X,A)') 
     &           ' Reduced ( E[2] - {w+ig}*S[2] )-1 * B[1] vector',
     &           ' with w =', FREQI,' and g =', 0d0,
     &           '================================================' 
              CALL OUTPUT(RDX(1,1,K),1,NOS,1,1,NOS,1,1,LUABSPRI)
              CALL OUTPUT(RDX(1,2,K),1,NOA,1,1,NOA,1,1,LUABSPRI)
              CALL OUTPUT(RDX(1,3,K),1,NOA,1,1,NOA,1,1,LUABSPRI)
              CALL OUTPUT(RDX(1,4,K),1,NOS,1,1,NOS,1,1,LUABSPRI)
           END IF
        END IF
      END DO
C
      RETURN
      END
      SUBROUTINE ABS_RSD(RDX,LABEL,KZVAR,NFREQS,FREQS,
     &           KLRED,KNVEC,NNLT,NGD,
     &           GD,RES,ITER,CONVERGED,CONV,XVEC,XSOL,LU,WRK,LWRK)
C
C PURPOSE:
C     Compute the complex residual vector to the linear response 
C     equation.
C
C     R = B[1] - ( E[2] - {w+iW}*S[2] )*(NR + iNI)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
C
      CHARACTER*8 LABEL, BLANK
      PARAMETER (BLANK='        ')
      INTEGER NFREQS,KLRED(2),KNVEC(2),NNLT,ITER,KZVAR,NGD
      LOGICAL CONV(NFREQS)
      PARAMETER (D0 = 0.0D0, D1=1.0D0)
      DIMENSION RDX(ABS_MAXRM,4,NFREQS),FREQS(NFREQS),GD(KZVAR,4,NGD),
     &     XVEC(KZVAR*NFREQS*4),XSOL(2*KZVAR,2,NFREQS),
     &     RES(KZVAR,NFREQS,4),WRK(LWRK)
      LOGICAL CONVERGED,ORTHO_QUIT1,ORTHO_QUIT2,SOLVER_QUIT
C     local
      INTEGER KZYVAR,K,J,NOS,NOA,NCONV,NSR,NSI,NAR,NAI,N1LT,N2LT
c      DOUBLE PRECISION DAMP,FREQI,FR1,FR2,FR3,FR4,DNORM_RR,DNORM_RI,
c     &                 DNORM_RT,DNORM_GR,DNORM_GI,DNORM_GT
C
      KZYVAR=2*KZVAR
C
      KBVEC = 1
      KTMP  = KBVEC + KZVAR
      KSVEC = KTMP  + KZYVAR
      KFREE = KSVEC + KZVAR
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) STOP 'ABS_RSD: out of memory'
C      
      CONVERGED = .TRUE.
      SOLVER_QUIT=.FALSE.
      DAMP=ABS_DAMP
C
C     Initialize residual vector
C  
      DO K=1,NFREQS
        IF (ABS_NGD) THEN
           IGD=K
        ELSE
           IGD=1
        END IF
        IF (.NOT.CONV(K)) THEN
           CALL DCOPY(KZVAR,GD(1,1,IGD),1,RES(1,K,1),1)
           CALL DCOPY(KZVAR,GD(1,2,IGD),1,RES(1,K,2),1)
           CALL DSCAL(KZVAR,-1.0d0,RES(1,K,1),1)
           CALL DSCAL(KZVAR,-1.0d0,RES(1,K,2),1)
           IF (ABS_GDCOMPLEX) THEN
              CALL DCOPY(KZVAR,GD(1,3,IGD),1,RES(1,K,3),1)
              CALL DCOPY(KZVAR,GD(1,4,IGD),1,RES(1,K,4),1)
              CALL DSCAL(KZVAR,-1.0d0,RES(1,K,3),1)
              CALL DSCAL(KZVAR,-1.0d0,RES(1,K,4),1)
           END IF
        ENDIF
      ENDDO
C
C     Read trial and sigma vectors and add contributions to residual vector
C
      REWIND(LUSB)
      REWIND(LUAB)
      REWIND(LUSS) 
      REWIND(LUAS) 
c
      DO J=1,NFREQS 
        IF (.NOT.CONV(J)) THEN
c          XSOL(1,1,J:J+2*KZYVAR-1)=0.0d0
          CALL DZERO(XSOL(1,1,J),2*KZYVAR)
        ENDIF
      ENDDO
c
      NOS=KLRED(1)
      NOA=KLRED(2)
c
C     Sigma vectors used that equal E[2]*N
C
      DO J = 1, NOS
         CALL READ_VEC(LUSS,KZVAR,WRK(KBVEC))
c         CALL READT(LUSS,KZVAR,WRK(KBVEC))
         DO K=1,NFREQS
            IF (.NOT.CONV(K)) THEN
              CALL DAXPY(KZVAR,RDX(J,1,K),WRK(KBVEC),1,RES(1,K,1),1)
              CALL DAXPY(KZVAR,-RDX(J,4,K),WRK(KBVEC),1,RES(1,K,4),1)
            ENDIF
         ENDDO
      ENDDO
      DO J = 1, NOA
         CALL READ_VEC(LUAS,KZVAR,WRK(KBVEC))
c          CALL READT(LUAS,KZVAR,WRK(KBVEC))
         DO K=1,NFREQS
            IF (.NOT.CONV(K)) THEN
                CALL DAXPY(KZVAR,RDX(J,2,K),WRK(KBVEC),1,RES(1,K,2),1)
                CALL DAXPY(KZVAR,-RDX(J,3,K),WRK(KBVEC),1,RES(1,K,3),1)
            ENDIF
         ENDDO
      ENDDO
C     
C     Trial vectors used to perform S[2]*N
C
      DO J = 1, NOS
         CALL READ_VEC(LUSB,KZVAR,WRK(KBVEC))
c          CALL READT(LUSB,KZVAR,WRK(KBVEC))
         CALL GETSVEC(KZVAR,1,WRK(KBVEC),WRK(KSVEC))
         DO K=1,NFREQS
            IF (.NOT.CONV(K)) THEN
              IF (ABS_IMFREQ) THEN
                 DAMP=FREQS(K)
                 FREQI=0.0d0
              ELSE
                 FREQI=FREQS(K)
              ENDIF
              CALL DAXPY(KZVAR,RDX(J,1,K),WRK(KBVEC),1,XSOL(1,1,K),1)
              CALL DAXPY(KZVAR,RDX(J,1,K),WRK(KBVEC),1,
     &             XSOL(KZVAR+1,1,K),1)
              FR1=-FREQI*RDX(J,1,K)
              FR2=DAMP*RDX(J,1,K)
              CALL DAXPY(KZVAR,FR1,WRK(KSVEC),1,RES(1,K,2),1)
              CALL DAXPY(KZVAR,FR2,WRK(KSVEC),1,RES(1,K,3),1)
c
              CALL DAXPY(KZVAR,RDX(J,4,K),WRK(KBVEC),1,XSOL(1,2,K),1)
              CALL DAXPY(KZVAR,RDX(J,4,K),WRK(KBVEC),1,
     &             XSOL(KZVAR+1,2,K),1)
              FR3=FREQI*RDX(J,4,K)
              FR4=DAMP*RDX(J,4,K)
              CALL DAXPY(KZVAR,FR3,WRK(KSVEC),1,RES(1,K,3),1)
              CALL DAXPY(KZVAR,FR4,WRK(KSVEC),1,RES(1,K,2),1)
            ENDIF
         ENDDO
      ENDDO
C     
      DO J = 1, NOA
         CALL READ_VEC(LUAB,KZVAR,WRK(KBVEC))
c          CALL READT(LUAB,KZVAR,WRK(KBVEC))
         CALL GETSVEC(KZVAR,1,WRK(KBVEC),WRK(KSVEC))
         DO K=1,NFREQS
            IF (.NOT.CONV(K)) THEN
              IF (ABS_IMFREQ) THEN
                 DAMP=FREQS(K)
                 FREQI=0.0d0
              ELSE
                 FREQI=FREQS(K)
              ENDIF
              CALL DAXPY(KZVAR,RDX(J,2,K),WRK(KBVEC),1,XSOL(1,1,K),1)
              CALL DAXPY(KZVAR,-RDX(J,2,K),WRK(KBVEC),1,
     &             XSOL(KZVAR+1,1,K),1)
              FR1=-FREQI*RDX(J,2,K)
              FR2=DAMP*RDX(J,2,K)
              CALL DAXPY(KZVAR,FR1,WRK(KSVEC),1,RES(1,K,1),1)
              CALL DAXPY(KZVAR,FR2,WRK(KSVEC),1,RES(1,K,4),1)
c
              CALL DAXPY(KZVAR,RDX(J,3,K),WRK(KBVEC),1,XSOL(1,2,K),1)
              CALL DAXPY(KZVAR,-RDX(J,3,K),WRK(KBVEC),1,
     &             XSOL(KZVAR+1,2,K),1)
              FR3=FREQI*RDX(J,3,K)
              FR4=DAMP*RDX(J,3,K)
              CALL DAXPY(KZVAR,FR3,WRK(KSVEC),1,RES(1,K,4),1)
              CALL DAXPY(KZVAR,FR4,WRK(KSVEC),1,RES(1,K,1),1)
            ENDIF
         ENDDO
      ENDDO
C     
      IF (IPRABSLRS.GE.200) THEN
         DO K=1,NFREQS
            IF (.NOT.CONV(K)) THEN
              WRITE(LUABSPRI,'(2(/,5X,A))') ' REAL SYMM. RES. COMP.',
     &           '========================'
              CALL OUTPUT(RES(1,K,1),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
c
              WRITE(LUABSPRI,'(2(/,5X,A))')' REAL ANTISYMM. RES. COMP.',
     &          '========================'
              CALL OUTPUT(RES(1,K,2),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
c
              WRITE(LUABSPRI,'(2(/,5X,A))') ' IMG ANTISYMM. RES. COMP.',
     &          '========================'
              CALL OUTPUT(RES(1,K,3),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
c
              WRITE(LUABSPRI,'(2(/,5X,A))') ' IMG SYMM. RES. COMP.',
     &          '========================'
              CALL OUTPUT(RES(1,K,4),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
            ENDIF
         ENDDO
      END IF
C
      NCONV=0
      DO K=1,NFREQS
         IF (ABS_NGD) THEN
            IGD=K
         ELSE
            IGD=1
         ENDIF
         IF (.NOT.CONV(K)) THEN
            DNORM_X=SQRT(DNRM2(KZYVAR,XSOL(1,1,K),1)+
     &                   DNRM2(KZYVAR,XSOL(1,2,K),1))
            CALL DCOPY(KZVAR,RES(1,K,1),1,WRK(KTMP),1)
            CALL DCOPY(KZVAR,RES(1,K,1),1,WRK(KTMP+KZVAR),1)
            CALL DAXPY(KZVAR,1.0d0,RES(1,K,2),1,WRK(KTMP),1)
            CALL DAXPY(KZVAR,-1.0d0,RES(1,K,2),1,WRK(KTMP+KZVAR),1)
            DNORM_RR = DNRM2(KZYVAR,WRK(KTMP),1)
c            CALL DZERO(WRK(KTMP),KZYVAR)
            WRK(KTMP:KTMP+KZYVAR-1) = 0.0d0
            CALL DCOPY(KZVAR,RES(1,K,4),1,WRK(KTMP),1)
            CALL DCOPY(KZVAR,RES(1,K,4),1,WRK(KTMP+KZVAR),1)
            CALL DAXPY(KZVAR,1.0d0,RES(1,K,3),1,WRK(KTMP),1)
            CALL DAXPY(KZVAR,-1.0d0,RES(1,K,3),1,WRK(KTMP+KZVAR),1)
            DNORM_RI = DNRM2(KZYVAR,WRK(KTMP),1)
            DNORM_RT = SQRT(DNORM_RR**2 + DNORM_RI**2)
            DNORM_GR = DNRM2(KZVAR,GD(1,1,IGD),1)
            DNORM_GI = DNRM2(KZVAR,GD(1,2,IGD),1)
            DNORM_GT = SQRT(2.0d0*(DNORM_GR**2 + DNORM_GI**2))
            IF (ABS_GDCOMPLEX) THEN
             DNORM_GRR = DNRM2(KZVAR,GD(1,3,IGD),1)
             DNORM_GII = DNRM2(KZVAR,GD(1,4,IGD),1)
             DNORM_GT = DNORM_GT+SQRT(2.0d0*(DNORM_GRR**2+DNORM_GII**2))
            END IF
             
C     
            DNORM_REL=DNORM_GT
            IF (ABS_XREL) THEN
              DNORM_REL=DNORM_X
            ENDIF
            RESID=DNORM_RT/DNORM_REL
C    
            IF (DNORM_RT/DNORM_REL.LE.ABS_THCLR) THEN
               NCONV=NCONV+1
               CONV(K)=.TRUE.
C               CALL WRITE_XVEC(LU,2*KZYVAR,XSOL(1,1,K),LABEL,
C     &              FREQS(K),RESID)
               CALL WRITE_XVEC2(LU,2*KZYVAR,XSOL(1,1,K),LABEL,
     &              BLANK,FREQS(K),0.0d0,RESID)
            ELSEIF (ITER.EQ.ABS_MAXITER) THEN
               CONVERGED=.FALSE.
               SOLVER_QUIT=.TRUE.
               IF (RESID .LE. 10.0d0*ABS_THCLR) THEN
C                  CALL WRITE_XVEC(LU,2*KZYVAR,XSOL(1,1,K),LABEL,
C     &                 FREQS(K),RESID)
                  CALL WRITE_XVEC2(LU,2*KZYVAR,XSOL(1,1,K),LABEL,
     &                 BLANK,FREQS(K),0.0d0,RESID)
                  WRITE(LUABSPRI,*) 'WARNING! For',FREQS(K)
                  WRITE(LUABSPRI,*) 'CPP equations not converged'
                 WRITE(LUABSPRI,*) 'Threshold in iteration',ITER,RESID
                 WRITE(LUABSPRI,*) 'Requested threshold', ABS_THCLR
               ELSE
                 WRITE(LUABSPRI,*) 'Maximum number of iter. exceeded'
                 WRITE(LUABSPRI,*) 'and CPP equations not converged!'
                 WRITE(LUABSPRI,*) 'Threshold in iteration',ITER,RESID
                 WRITE(LUABSPRI,*) 'Requested threshold', ABS_THCLR
          WRITE(LUABSPRI,*) 'Increase the iteration number or threshold'
                 CALL QUIT('CPP equations not converged!')
               ENDIF
              
            ELSE
               CONVERGED=.FALSE.
            ENDIF
C
C     Print norm of residual vector
C
            IF (IPRABSLRS.GE.2) THEN
               WRITE(LUABSPRI,*)'-------------------------------------'
               write(luabspri,*) 'dimension of reduced spaces:',NOS,NOA
                write(luabspri,*)
     &            'Frequency:',FREQS(K),
     &            'Damping:',ABS_DAMP,
     &            'Residual:',DNORM_RT/DNORM_GT
            END IF
            WRITE(LUABSPRI,'(A,2(I4),2X,A,2X,F7.4,2X,A,2(D15.7))')
     &             'ITER, IVEC:',ITER,K, LABEL,FREQS(K),
     &             ' RESID, RESID/X:',DNORM_RT/DNORM_GT,DNORM_RT/DNORM_X
c            write(luabspri,*)'residual ',LABEL,ITER,FREQS(K),
c     &        DNORM_RT/DNORM_GT
         ELSE
            NCONV=NCONV+1
            WRITE(LUABSPRI,'(A,2(I3),2X,A,2X,F7.4,A,2(L1))')
     *           'ITER, IVEC:',ITER,K,LABEL,FREQS(K),
     *           ' RESID, RESID/X:',CONV(K),CONV(K)
c            write(luabspri,*)'residual ',LABEL,ITER,FREQS(K),
c     &        CONV(K)
c           write(luabspri,*) 'Freq skipped in ABS_RSD',FREQS(K),CONV(K)
         END IF
       END DO
       IF (.NOT. CONVERGED .AND. (.NOT. SOLVER_QUIT)) THEN
          CALL ABS_PRECON(KZVAR,RES,FREQS,NFREQS,XVEC,CONV,
     &                    WRK(KFREE),LFREE)
            KNVEC(1)=2*NFREQS
            KNVEC(2)=2*NFREQS
c          RES(1:KZVAR,1:NFREQS,1:4)=0.0d0
          CALL DZERO(RES,4*NFREQS*KZVAR)
          ORTHO_QUIT1=.TRUE.
          ORTHO_QUIT2=.TRUE.
          NSR=0
          NSI=0
          NAR=0
          NAI=0
          CALL ABS_ORTHNORM(KZVAR,KLRED(1),KNVEC(1),
     &         XVEC(1),XVEC(3*NFREQS*KZVAR+1),LUSB,RES,NSR,NSI,
     &         CONV,ORTHO_QUIT1,WRK(KFREE),LFREE)
          N1LT=KNVEC(1)
          CALL ABS_ORTHNORM(KZVAR,KLRED(2),KNVEC(2),
     &         XVEC(NFREQS*KZVAR+1),XVEC(2*NFREQS*KZVAR+1),
     &         LUAB,RES(1,1,3),NAR,NAI,CONV,ORTHO_QUIT2,WRK(KFREE),
     &         LFREE)
          N2LT=KNVEC(2)
          CALL DCOPY(N1LT*KZVAR,RES,1,XVEC,1)
          CALL DCOPY(N2LT*KZVAR,RES(1,1,3),1,XVEC(1+N1LT*KZVAR),1)
          IF (ORTHO_QUIT1 .AND. ORTHO_QUIT2)
c     &        STOP 'No new trial vectors!'
     &        CALL QUIT('ERROR IN ABS_ORTHO, no new trial vectors!')
c
          NNLT=KNVEC(1)+KNVEC(2)
c       write(luabspri,*) KZRED,NSR,NSI,NAR,NAI,NCONV
       ENDIF
C
      RETURN
      END
      SUBROUTINE ABS_PRECON(KZVAR,RR,FREQS,NFREQS,XVEC,CONV,WRK,LWRK)
C
C PURPOSE:
C      perform preconditioning 
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
C
      INTEGER NFREQS,KZVAR
      LOGICAL CONV(NFREQS)
      DIMENSION RR(KZVAR,NFREQS,4),XVEC(KZVAR,NFREQS,4),
     &          FREQS(NFREQS),WRK(LWRK)
c     local
      INTEGER K,I,KFREE,LFREE,KDIAE,KS2
C
      KDIAE = 1
      KS2   = KDIAE + KZVAR
      KFREE = KS2   + KZVAR
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) STOP 'ABS_PRECON: out of memory'
c
      DAMP=ABS_DAMP
      WRK(KDIAE:KDIAE+KZVAR-1) = 0.0d0
      WRK(KS2:KS2+KZVAR-1) = 0.0d0
      CALL RSPEDG(WRK(KDIAE))
      CALL ABS_GETS2(WRK(KS2),KZVAR)
      IF (IPRABSLRS.GE.150) THEN
         WRITE(LUABSPRI,*)' DIAGONAL HESSIAN MATRIX WVAL=0.0'
         CALL OUTPUT(WRK(KDIAE),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
         WRITE(LUABSPRI,*)' DIAGONAL S2 MATRIX WVAL=0.0'
         CALL OUTPUT(WRK(KS2),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
      ENDIF
         
      DO K=1,NFREQS
         DO I=1,KZVAR
            IF (ABS(RR(I,K,1)) .LT. 1.0d-12) RR(I,K,1)=0.0d0 
            IF (ABS(RR(I,K,2)) .LT. 1.0d-12) RR(I,K,2)=0.0d0 
            IF (ABS(RR(I,K,3)) .LT. 1.0d-12) RR(I,K,3)=0.0d0 
            IF (ABS(RR(I,K,4)) .LT. 1.0d-12) RR(I,K,4)=0.0d0
         ENDDO 
      ENDDO
      DO K=1,NFREQS
         IF (.NOT.CONV(K)) THEN
            IF (ABS_IMFREQ) THEN
               FREQI=0.0d0
               DAMP=FREQS(K)
            ELSE
               FREQI=FREQS(K)
            ENDIF
            DO I=1,KZVAR
               D2=WRK(KS2-1+I)*WRK(KS2-1+I)
               PR1=WRK(KDIAE-1+I)*WRK(KDIAE-1+I)-
     &              D2*(FREQI*FREQI-DAMP*DAMP)
               DP=PR1*PR1+4.0d0*D2*D2*FREQI*FREQI*DAMP*DAMP
               PRA=WRK(KDIAE-1+I)*PR1
               PRB=WRK(KS2-1+I)*FREQI*(WRK(KDIAE-1+I)*WRK(KDIAE-1+I)-
     &              D2*(FREQI*FREQI+DAMP*DAMP))
               PRC=WRK(KS2-1+I)*DAMP*(WRK(KDIAE-1+I)*WRK(KDIAE-1+I)+
     &              D2*(FREQI*FREQI+DAMP*DAMP))
               PRD=2.0d0*D2*FREQI*DAMP*WRK(KDIAE-1+I)
C     
               XVEC(I,K,1)=(PRA*RR(I,K,1)+PRB*RR(I,K,2)+
     &              PRC*RR(I,K,3)+PRD*RR(I,K,4))/DP
C     
               XVEC(I,K,2)=(PRB*RR(I,K,1)+PRA*RR(I,K,2)+
     &              PRD*RR(I,K,3)+PRC*RR(I,K,4))/DP
C     
               XVEC(I,K,3)=(PRC*RR(I,K,1)+PRD*RR(I,K,2)-
     &              PRA*RR(I,K,3)-PRB*RR(I,K,4))/DP
C     
               XVEC(I,K,4)=(PRD*RR(I,K,1)+PRC*RR(I,K,2)-
     &              PRB*RR(I,K,3)-PRA*RR(I,K,4))/DP
C     
            END DO
            IF (IPRABSLRS.GE.100) THEN
               WRITE(LUABSPRI,*)'UPPER HALFS OF TRIAL VECTORS
     & AFTER PRECONDITIONING'
               WRITE(LUABSPRI,'(2(/,5X,A))') ' REAL SYMMETRIC VECTOR',
     &              '========================'
               CALL OUTPUT(XVEC(1,K,1),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
c     
               WRITE(LUABSPRI,'(2(/,5X,A))')
     &              ' REAL ANTISYMMETRIC VECTOR',
     &              '========================'
               CALL OUTPUT(XVEC(1,K,2),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
c     
               WRITE(LUABSPRI,'(2(/,5X,A))')' IMG ANTISYMMETRIC VECTOR',
     &              '========================'
               CALL OUTPUT(XVEC(1,K,3),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
c     
               WRITE(LUABSPRI,'(2(/,5X,A))') ' IMG SYMMETRIC VECTOR',
     &              '========================'
               CALL OUTPUT(XVEC(1,K,4),1,KZVAR,1,1,KZVAR,1,1,LUABSPRI)
            END IF
         ENDIF
      END DO
C
      RETURN
      END 
      SUBROUTINE ABS_GETS2(S2M,KZVAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER KZVAR
      DIMENSION S2M(KZVAR)
C     local
      INTEGER NKD,ISYM1,ISYM2,NVRT2,NEL
C
#include "abslrs.h"
C
C loop over symmetries
C
      NKD=0
      DO ISYM1 = 1,ABS_NSYM
         ISYM2 = ABS_MULD2H(ISYM1,ABS_GRADSYM)
         NVRT2 = ABS_NORB(ISYM2)-ABS_NISH(ISYM2)-ABS_NASH(ISYM2)
c         write(luabspri,*)'isym1,isym2,na1,na2,ni1,ni2,nvrt2',isym1,isym2,
c     &                  nash1,nash2,nish1,nish(isym2),norb(isym1),
c     &                  norb(isym2),nvrt2
          NEL=ABS_NISH(ISYM1)*ABS_NASH(ISYM2)
          IF (NEL.GT.0) THEN
            DO I=1,NEL
              S2M(NKD+I)=1.0d0
            ENDDO
            NKD=NKD+NEL
          ENDIF
C
          NEL=ABS_NISH(ISYM1)*NVRT2
          IF (NEL.GT.0) THEN
            DO I=1,NEL
              S2M(NKD+I)=2.0d0
            ENDDO
            NKD=NKD+NEL
          ENDIF
C
          NEL=ABS_NASH(ISYM1)*NVRT2
          IF (NEL.GT.0) THEN
            DO I=1,NEL
              S2M(NKD+I)=1.0d0
            ENDDO
            NKD=NKD+NEL
          ENDIF
       ENDDO
C
      RETURN
      END
      SUBROUTINE GETSVEC(KZVAR,NOSIM,VECB,SVEC)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOSIM,KZVAR
      DIMENSION VECB(KZVAR,NOSIM),SVEC(KZVAR,NOSIM)
C     local
      INTEGER NKD,K,ISYM1,ISYM2,NVRT2,NEL
C
#include "abslrs.h"
C
C loop over symmetries
C
      CALL DCOPY(KZVAR*NOSIM,VECB,1,SVEC,1)
      DO K=1,NOSIM
        NKD=0
        DO ISYM1 = 1,ABS_NSYM
          ISYM2 = ABS_MULD2H(ISYM1,ABS_GRADSYM)
          NVRT2 = ABS_NORB(ISYM2)-ABS_NISH(ISYM2)-ABS_NASH(ISYM2)
          NEL=ABS_NISH(ISYM1)*ABS_NASH(ISYM2)
          NKD=NKD+NEL
          NEL=ABS_NISH(ISYM1)*NVRT2
            IF (NEL .GT.0) THEN
              CALL DSCAL(NEL,2.0d0,SVEC(NKD+1,K),1)
            ENDIF
          NKD=NKD+NEL
          NEL=ABS_NASH(ISYM1)*NVRT2
          NKD=NKD+NEL
        ENDDO
      ENDDO
C 
      RETURN
      END
      SUBROUTINE ABS_RED(ITER,KZVAR,GD,XVEC,XTMP,REDGD,REDE,REDS,NFREQS,
     &                   NNLT,KLRED,KNVEC,NGD,MJWOP,CMO,UDV,FC,FCAC,FV,
     &                   PV,
     &                   XINDX,WRK,LWRK)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
C
C PURPOSE:
C     extend reduce subspace
C
      INTEGER NFREQS,NNLT,KZVAR,NGD
      INTEGER KLRED(2),KNVEC(2)
      DIMENSION GD(4*KZVAR,NGD),XVEC(KZVAR,NNLT)
      DIMENSION REDE(ABS_MAXRM,ABS_MAXRM,2),REDS(ABS_MAXRM,ABS_MAXRM),
     &          REDGD(ABS_MAXRM,4,NGD),E2XF(KZVAR,NNLT),
     &          XTMP(4*NFREQS*KZVAR)
      DIMENSION MJWOP(*),XINDX(*)
      DIMENSION CMO(*),UDV(*),FC(*),FCAC(*),FV(*),PV(*)
      DIMENSION WRK(LWRK)
      PARAMETER (D1 = -1.0D0)
C     local
      INTEGER NNS,NNA,NOS,NOA,NOLT,NNMAX,I,J,K
C
       NNS=KNVEC(1)
       NNA=KNVEC(2)
       NOS=KLRED(1)
       NOA=KLRED(2)
       NOLT=NOS+NOA
       NNMAX=MAX(NNS,NNA)
C
      KBVEC = 1
      KE2XF = KBVEC + KZVAR
      KSVEC = KE2XF + NNLT*KZVAR
      KFREE = KSVEC + 4*NFREQS*KZVAR
      LFREE = LWRK - KFREE
      IF (LFREE.LT.0) STOP 'ABS_RED: out of memory'
c
      IF (((KLRED(1)+KNVEC(1)).GT.ABS_MAXRM).OR.
     &     ((KLRED(2)+KNVEC(2)).GT.ABS_MAXRM)) THEN 
          CALL QUIT('Error in ABS_RED, Increase ABS_MAXRM!')
      ENDIF
      IF (NNLT.GT.0) THEN
c
C       perform linear transformations
C
Cpanor
c         IF (ABS_DALINP) THEN
         CALL ABS_LINE2_INTER(KZVAR,WRK(KE2XF),XVEC,XTMP,WRK(KSVEC),
     &        ABS_GRADSYM,NNLT,NNMAX,KNVEC,CMO,FC,FCAC,FV,UDV,XINDX,
     &         MJWOP,WRK(KFREE),LFREE)
C         ELSE
C         get linear transformations A(x) -> WRK(KE2XF)
C         ENDIF
C
        WRK(KSVEC:KSVEC+4*NFREQS*KZVAR-1)=0.0d0
        CALL GETSVEC(KZVAR,NNLT,XVEC,WRK(KSVEC))
        IF (IPRABSLRS.GE.200) THEN
           WRITE(LUABSPRI,'(//A)') 
     &          ' RHO VECTORS (IN ORDER REAL-SYM'//
     &          ' REAL-ANTYSYM., IMG-ANTYSYM., IMG-SYM.'
           CALL OUTPUT(WRK(KSVEC),1,KZVAR,1,NNLT,KZVAR,NNLT,1,LUABSPRI)
        ENDIF
C
C       Construct the reduced gradient with trial vectors
C
       DO J=1,NGD 
        DO I=1,NNS
          REDGD(NOS+I,1,J) = 2.0d0*DDOT(KZVAR,XVEC(1,I),1,
     &                     GD(1,J),1)
        END DO
        DO I=1,NNA
          REDGD(NOA+I,2,J) = 2.0d0*DDOT(KZVAR,XVEC(1,NNS+I),1,
     &                     GD(1+KZVAR,J),1)
        END DO 
        IF (ABS_GDCOMPLEX) THEN
         DO I=1,NNA
          REDGD(NOA+I,3,J) = 2.0d0*DDOT(KZVAR,XVEC(1,NNS+I),1,
     &                     GD(1+2*KZVAR,J),1)
         END DO
         DO I=1,NNS
          REDGD(NOS+I,4,J) = 2.0d0*DDOT(KZVAR,XVEC(1,I),1,
     &                     GD(1+3*KZVAR,J),1)
         END DO 
        END IF
       END DO
C
C        extend reduced spaces
c
        IF (NNS.GT.0) THEN
          DO I=1,NNS
            DO J=I,NNS
              REDE(NOS+I,NOS+J,1) = 2.0d0*DDOT(KZVAR,
     &            XVEC(1,I),1,WRK(KE2XF+(J-1)*KZVAR),1)
               IF (I.NE.J) REDE(NOS+J,NOS+I,1)=REDE(NOS+I,NOS+J,1)
            ENDDO
          ENDDO
        END IF
        IF (NNA.GT.0) THEN
          DO I=1,NNA
            DO J=I,NNA
              REDE(NOA+I,NOA+J,2) = 2.0d0*DDOT(KZVAR,
     &           XVEC(1,NNS+I),1,WRK(KE2XF+(NNS+J-1)*KZVAR),1)
               IF (I.NE.J) REDE(NOA+J,NOA+I,2)=REDE(NOA+I,NOA+J,2)
            ENDDO
          ENDDO
        END IF
c
        IF ((NNS.GT.0).AND.(NNA.GT.0)) THEN
          DO I=1,NNS
            DO J=1,NNA
              REDS(NOS+I,NOA+J) = 2.0D0*DDOT(KZVAR,
     &             XVEC(1,I),1,WRK(KSVEC+(NNS+J-1)*KZVAR),1)
            END DO
          ENDDO
        END IF
c
        IF (NOLT .GT. 0) THEN
          REWIND(LUSB)
          REWIND(LUAB)
          DO I=1,NOS
            CALL READ_VEC(LUSB,KZVAR,WRK(KBVEC))
c             CALL READT(LUSB,KZVAR,WRK(KBVEC))
            DO J=1,NNS
              REDE(I,NOS+J,1)= 2.0d0*DDOT(KZVAR,WRK(KBVEC),1,
     &           WRK(KE2XF+(J-1)*KZVAR),1)
              REDE(NOS+J,I,1)=REDE(I,NOS+J,1)
            ENDDO
            DO J=1,NNA
              REDS(I,NOA+J)=2.0d0*DDOT(KZVAR,WRK(KBVEC),1,
     &            WRK(KSVEC+(NNS+J-1)*KZVAR),1)
            ENDDO
            IF ((ITER.EQ.1) .AND. (NOLT.GT.1)) THEN 
              DO J=1,NGD
               REDGD(I,1,J) = 2.0d0*DDOT(KZVAR,WRK(KBVEC),1,GD(1,J),1)
              END DO
            ENDIF
          END DO
          DO I=1,NOA
            CALL READ_VEC(LUAB,KZVAR,WRK(KBVEC))
c             CALL READT(LUAB,KZVAR,WRK(KBVEC))
            DO J=1,NNA
              REDE(I,NOA+J,2)=2.0d0*DDOT(KZVAR,WRK(KBVEC),1,
     &           WRK(KE2XF+(NNS+J-1)*KZVAR),1)
              REDE(NOA+J,I,2)=REDE(I,NOA+J,2)
            ENDDO
            DO J=1,NNS
              REDS(NOS+J,I)=2.0d0*DDOT(KZVAR,WRK(KBVEC),1,
     &             WRK(KSVEC+(J-1)*KZVAR),1)
            ENDDO
            IF ((ITER.EQ.1) .AND. (NOLT.GT.1)) THEN 
              DO J=1,NGD
               REDGD(I,2,J) = 2.0d0*DDOT(KZVAR,WRK(KBVEC),
     &                             1,GD(1+KZVAR,J),1)
              END DO
            ENDIF
          END DO
        END IF
C
        DO I=1,NNS
c          CALL WRITT(LUSB,KZVAR,XVEC(1,I))
c          CALL WRITT(LUSS,KZVAR,WRK(KE2XF+(I-1)*KZVAR))
          CALL WRITE_VEC(LUSB,KZVAR,XVEC(1,I))
          CALL WRITE_VEC(LUSS,KZVAR,WRK(KE2XF+(I-1)*KZVAR))
        ENDDO
        DO I=1,NNA
c          CALL WRITT(LUAB,KZVAR,XVEC(1,NNS+I))
c          CALL WRITT(LUAS,KZVAR,WRK(KE2XF+(NNS+I-1)*KZVAR))
          CALL WRITE_VEC(LUAB,KZVAR,XVEC(1,NNS+I))
          CALL WRITE_VEC(LUAS,KZVAR,WRK(KE2XF+(NNS+I-1)*KZVAR))
        ENDDO
C
        KLRED(1)=KLRED(1)+KNVEC(1)
        KLRED(2)=KLRED(2)+KNVEC(2)
      ELSE
        REWIND(LUSB) 
        REWIND(LUAB)
        DO I=1,NOS
          CALL READ_VEC(LUSB,KZVAR,WRK(KBVEC))
c           CALL READT(LUSB,KZVAR,WRK(KBVEC))
          IF ((ITER.EQ.1) .AND. (NOLT.GT.1)) THEN 
            DO J=1,NGD
             REDGD(I,1,J) = 2.0d0*DDOT(KZVAR,WRK(KBVEC),1,GD(1,J),1)
            ENDDO
          ELSE
            CALL QUIT('ERROR IN ABS_RED, no trial vectors!')
          ENDIF
        END DO
        DO I=1,NOA
          CALL READ_VEC(LUAB,KZVAR,WRK(KBVEC))
c           CALL READT(LUAB,KZVAR,WRK(KBVEC))
          IF ((ITER.EQ.1) .AND. (NOLT.GT.1)) THEN 
            DO J=1,NGD
             REDGD(I,2,J) = 2.0d0*DDOT(KZVAR,WRK(KBVEC),
     &                   1,GD(1+KZVAR,J),1)
            END DO
          ELSE
            CALL QUIT('ERROR IN ABS_RED, no trial vectors!')
          ENDIF
        END DO
      END IF
c
      RETURN
      END
      SUBROUTINE ABS_GETRIAL(KZVAR,GD,GDSYM,NGD,FREQS,NFREQS,XVEC,TEMP,
     &                       KLRED,KNVEC,NNLT,CONV,WRK,LWRK)
C
C PURPOSE
C    1.  Split RHS to the symmetric and antisymmetric components
C    2.  Precondition RHS to obtain trial vectors in order 
C        /RS RA IA IS/ [where R(eal)S(ymmetric), w-freqi]
C    3.  Othonormalize trial vectors in order
C        /RS IS RA IA/
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "abslrs.h"
      INTEGER NFREQS,KLRED(2),KNVEC(2),NNLT,KZVAR,NGD
      LOGICAL CONV(NFREQS)
      DIMENSION GD(4*KZVAR,NGD),GDSYM(KZVAR,4,NGD),FREQS(NFREQS)
      DIMENSION XVEC(KZVAR,NFREQS,4),TEMP(4*NFREQS*KZVAR),WRK(LWRK)
      LOGICAL ORTHO_QUIT
C     local
      INTEGER KFREE,LFREE,I,J,NSR,NSI,NAR,NAI,N1LT
c     
      KFREE=1
      LFREE=LWRK
C
      DO J=1,NGD
       DO I=1,KZVAR
           GDSYM(I,1,J)=0.5d0*(GD(I,J)+GD(KZVAR+I,J))
           GDSYM(I,2,J)=0.5d0*(GD(I,J)-GD(KZVAR+I,J))
           IF (ABS_GDCOMPLEX) THEN
              GDSYM(I,3,J)=-0.5d0*(GD(2*KZVAR+I,J)-GD(3*KZVAR+I,J))
              GDSYM(I,4,J)=-0.5d0*(GD(2*KZVAR+I,J)+GD(3*KZVAR+I,J))
           END IF
       END DO
      END DO
      IF (ABS_GDCOMPLEX) THEN
          CALL DCOPY(4*NGD*KZVAR,GDSYM,1,TEMP,1)
      ELSE
      DO J=1,2
        DO I=1,NFREQS
          CALL DCOPY(KZVAR,GDSYM(1,J,1),1,
     &         TEMP(1+((J-1)*NFREQS+(I-1))*KZVAR),1)
        ENDDO
      ENDDO
      END IF
c    
      CALL ABS_PRECON(KZVAR,TEMP,FREQS,NFREQS,XVEC,CONV,
     &                WRK(KFREE),LFREE)
          NSR=0
          NSI=0
          NAR=0
          NAI=0
      ORTHO_QUIT=.TRUE.
      CALL ABS_ORTHNORM(KZVAR,KLRED(1),KNVEC(1),XVEC(1,1,1),XVEC(1,1,4),
     &            LUSB,TEMP,NSR,NSI,CONV,ORTHO_QUIT,WRK(KFREE),LFREE)
      N1LT=KNVEC(1)
      CALL ABS_ORTHNORM(KZVAR,KLRED(2),KNVEC(2),XVEC(1,1,2),XVEC(1,1,3),
     &            LUAB,TEMP(1+N1LT*KZVAR),NAR,NAI,CONV,ORTHO_QUIT,
     &            WRK(KFREE),LFREE)
c
       NNLT=KNVEC(1)+KNVEC(2)
       CALL DCOPY(NNLT*KZVAR,TEMP,1,XVEC,1)
C
      RETURN
      END
      SUBROUTINE ABS_ORTHNORM(KZVAR,NBPREV,NBNEW,XRVEC,XIVEC,LU1,XORT,
     &                     NNR,NNI,CONV,ORTHO_QUIT,WRK,LWRK)
C
C PURPOSE:
C      preconditioning 
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NBPREV,NBNEW,NNR,NNI,KZVAR,LU1
      DIMENSION XRVEC(KZVAR,NBNEW/2),XIVEC(KZVAR,NBNEW/2)
      DIMENSION XORT(NBNEW*KZVAR),WRK(LWRK)
      LOGICAL CONV(NBNEW/2)
      LOGICAL ORCHK(NBNEW),ORTHO_QUIT
c      PARAMETER (THR_LD=1.0D-2*THCLR_ABSORP)
#include "abslrs.h"
C
      KTMP =1
      KFREE = KTMP + KZVAR
      LFREE = LWRK - KFREE
      IF (LFREE.LT.0) STOP 'ABS_ORTHONORM: out of memory'
C
      THR_LD=1.0D-8
      THR_LD2=0.0d0
c
      DO I=1,NBNEW/2
        ORCHK(I)=CONV(I)
        ORCHK(NBNEW/2+I)=CONV(I)
      ENDDO
c      
      DO ITURN=1,2
C
C     Orthonormalize trial vectors against the old vectors
C   
        IF (NBPREV.GT.0) THEN
c         write(luabspri,*)'nbprev',nbprev,nbnew
          REWIND(LU1)
          DO K = 1,NBPREV
            CALL READ_VEC(LU1,KZVAR,WRK(KTMP))
c             CALL READT(LU1,KZVAR,WRK(KTMP))
            DO I = 1,NBNEW/2
              IF (.NOT.ORCHK(I)) THEN
                TT = 2.0d0*DDOT(KZVAR,WRK(KTMP),1,XRVEC(1,I),1)
                CALL DAXPY(KZVAR,-TT,WRK(KTMP),1,XRVEC(1,I),1)
              ENDIF
            ENDDO
            DO I = 1,NBNEW/2
              IF (.NOT.ORCHK(NBNEW/2+I)) THEN
                TT = 2.0d0*DDOT(KZVAR,WRK(KTMP),1,XIVEC(1,I),1)
                CALL DAXPY(KZVAR,-TT,WRK(KTMP),1,XIVEC(1,I),1)
              ENDIF
            ENDDO
          ENDDO
        END IF
        DO I=1,NBNEW/2
          IF (.NOT.ORCHK(I)) THEN
            PR1=SQRT(2.0d0*DDOT(KZVAR,XRVEC(1,I),1,XRVEC(1,I),1))
            IF (PR1.LE.THR_LD) THEN
              ORCHK(I)=.TRUE.
            ELSE
              PR2=1.0d0/PR1
              CALL DSCAL(KZVAR,PR2,XRVEC(1,I),1)
            ENDIF
          ENDIF
          IF (.NOT.ORCHK(NBNEW/2+I)) THEN
            PI1=SQRT(2.0d0*DDOT(KZVAR,XIVEC(1,I),1,XIVEC(1,I),1))
            IF (PI1.LE.THR_LD) THEN
              ORCHK(NBNEW/2+I)=.TRUE.
            ELSE
              PI2=1.0d0/PI1
              CALL DSCAL(KZVAR,PI2,XIVEC(1,I),1)
            ENDIF
          ENDIF
        ENDDO
c
C Orthogonalize new vectors against each other
c   
        DO I = 1,NBNEW/2 !index for current bvector
          IF (.NOT.ORCHK(I)) THEN
            IF ((NBNEW/2.GT.1).AND.(I.GT.1)) THEN
              DO J = 1,(I-1)
                IF (.NOT.ORCHK(J)) THEN
                  T1 = 2.0d0*DDOT(KZVAR,XRVEC(1,J),1,XRVEC(1,J),1)
                  IF (T1 .GT. THR_LD2) THEN
                    T2 = 2.0d0*DDOT(KZVAR,XRVEC(1,J),1,XRVEC(1,I),1)
                    TT = T2/T1
                    CALL DAXPY(KZVAR,-TT,XRVEC(1,J),1,XRVEC(1,I),1)
                  ENDIF
                ENDIF
              END DO
            ENDIF
          ENDIF
        END DO
        DO I = 1,NBNEW/2 !index for current bvector
          IF (.NOT.ORCHK(NBNEW/2+I)) THEN
            DO J = 1,NBNEW/2
              IF (.NOT.ORCHK(J)) THEN
                T1 = 2.0d0*DDOT(KZVAR,XRVEC(1,J),1,XRVEC(1,J),1)
                IF (T1 .GT. THR_LD2) THEN
                  T2 = 2.0d0*DDOT(KZVAR,XRVEC(1,J),1,XIVEC(1,I),1)
                  TT = T2/T1
                  CALL DAXPY(KZVAR,-TT,XRVEC(1,J),1,XIVEC(1,I),1)
                ENDIF
              ENDIF
            END DO
            IF (NBNEW/2.GT.1) THEN
              DO J = 1,(I-1)
                IF (.NOT.ORCHK(NBNEW/2+J)) THEN
                  T1 = 2.0d0*DDOT(KZVAR,XIVEC(1,J),1,XIVEC(1,J),1)
                  IF (T1 .GT. THR_LD2) THEN
                    T2 = 2.0d0*DDOT(KZVAR,XIVEC(1,J),1,
     &                      XIVEC(1,I),1)
                    TT = T2/T1
                    CALL DAXPY(KZVAR,-TT,XIVEC(1,J),1,XIVEC(1,I),1)
                  ENDIF
                ENDIF
              END DO
            ENDIF
          ENDIF
        END DO
c 
        DO I =1,NBNEW/2 
          IF (.NOT.ORCHK(I)) THEN
            PR1=SQRT(2.0d0*DDOT(KZVAR,XRVEC(1,I),1,XRVEC(1,I),1))
            IF (PR1.LE.THR_LD) THEN
              IF (ITURN==1) THEN
                ORCHK(I)=.TRUE.
              ELSE
                CALL QUIT('error in R orthonormalization')
              ENDIF 
            ELSE
              PR2=1.0d0/PR1
              CALL DSCAL(KZVAR,PR2,XRVEC(1,I),1)
            ENDIF
          ENDIF
          IF (.NOT.ORCHK(NBNEW/2+I)) THEN
            PI1=SQRT(2.0d0*DDOT(KZVAR,XIVEC(1,I),1,XIVEC(1,I),1))
            IF (PI1.LE.THR_LD) THEN
              IF (ITURN==1) THEN
                ORCHK(NBNEW/2+I)=.TRUE.
              ELSE
                CALL QUIT('error in I orthonormalization')
              ENDIF 
            ELSE
              PI2=1.0d0/PI1
              CALL DSCAL(KZVAR,PI2,XIVEC(1,I),1)
            ENDIF
          ENDIF
        END DO
      ENDDO
c
      NNR=0
      DO I=1,NBNEW/2
        IF (.NOT.ORCHK(I)) THEN
          CALL DCOPY(KZVAR,XRVEC(1,I),1,XORT(1+NNR*KZVAR),1)
          NNR=NNR+1
          ORTHO_QUIT=.FALSE.
c        ELSE
c          write(luabspri,*)'vector skipped',I
        END IF  
      ENDDO
      NNI=0
      DO I=1,NBNEW/2
        IF (.NOT.ORCHK(NBNEW/2+I)) THEN
          CALL DCOPY(KZVAR,XIVEC(1,I),1,XORT(1+(NNR+NNI)*KZVAR),1)
          NNI=NNI+1
          ORTHO_QUIT=.FALSE.
c        ELSE
c          write(luabspri,*)'vector skipped',I
        END IF  
      ENDDO
c
       NBNEW=NNR+NNI
c         write(luabspri,*)'overlap matrix'
c        do i=1,Nbnew
c          do j=1,Nbnew
c           write(luabspri,*) I,J,2.0d0*ddot(kzvar,xort(1+(I-1)*KZVAR),1,
c     &                    xort(1+(J-1)*KZVAR),1)
c          enddo
c      enddo
c
      IF (IPRABSLRS.GE.200) THEN
         WRITE(LUABSPRI,'(2(/,5X,A))') ' NORMALIZED TRIAL VECTOR',
     &        '========================'
         CALL OUTPUT(XORT,1,KZVAR,1,NBNEW,KZVAR,NBNEW,1,LUABSPRI)
      END IF
C
      RETURN
      END
C
