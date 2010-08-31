C
C  /* Deck dc_r1vec */
      SUBROUTINE DC_R1VEC(R1VEC,LR1VEC,EXVAL,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995.
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate the first order part, R(1), of the eigenvector
C              in RPA(D) theory. See eq. (16) in RPA(D) paper.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ONE    = 1.0D0, TWO = 2.0D0 )
      DIMENSION R1VEC(LR1VEC),WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <inftap.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('DC_R1VEC')
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      KSCR1   = 1
      KEND    = KSCR1 + NORBTS
      LWORK   = LWORK - KEND
C
      CALL SO_MEMMAX ('DC_R1VEC',LWORK1)
      IF (LWORK .LT. 0) CALL STOPIT('DC_R1VEC',' ',KEND,LWORK)
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(I), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
      IF (FROIMP .OR. FROEXP) 
     &   CALL CCSD_DELFRO(WORK(KSCR1),WORK(KEND),LWORK)
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 ISYMBJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)
C
         IF (ISYMAI. GT. ISYMBJ) GOTO 100
C
         DO 200 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 300 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 400 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 500 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     &                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 600 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 700 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     &                         + NVIR(ISYMA)*(I - 1) + A
C
                           IF ( (ISYMAI .EQ. ISYMBJ) .AND. 
     &                          (NAI .GT. NBJ) ) GOTO 700
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           EDIFF = - (  WORK(MA) + WORK(MB)
     &                               - WORK(MI) - WORK(MJ) - EXVAL )
C
                           EDIFF = ONE / EDIFF
C
                           R1VEC(NAIBJ) = R1VEC(NAIBJ) * EDIFF
C
  700                   CONTINUE
C
  600                CONTINUE
C
  500             CONTINUE
C
  400          CONTINUE
C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('DC_R1VEC')
C
      RETURN
      END
