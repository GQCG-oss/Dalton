C  /* Deck cc_freeze_exci */
      SUBROUTINE CC_FREEZE_core(CAM,ISYMTR,
     &           MAXCORE, 
     &           NRHFCORE,IRHFCORE)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     09-7-2015 Sonia Coriani
!     Specular to CVS separation - could be replaced by CC_CORE()
C
C     Purpose: Project out specific CORE excitations 
C              from a trial vector (by zeroing 
C              specific elements)
C     Ex1: zero all ai and aibj elements where i and j
C     are CORE orbitals 
C
C Based on cc_pram()
! CAM is the vector analyzed, of symmetry ISYMTR
! Control is passed via argument list, not via common block
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      Implicit none

#include "ccsdsym.h"
      Double precision CAM(*)
      Integer MAXCORE, NRHFCORE(8),IRHFCORE(MAXCORE,8)
      integer ISYMTR,ISYMAI,ISYMI,ISYMA,ISYMJ,ISYMB,ISYMBJ
      Double precision TWO, THR1, THR2, zero
      PARAMETER (TWO = 2.0D0,zero=0.0d0)
      Logical LOCDBG, ikeep, LBOTH
      Parameter (Locdbg = .false.)
      Integer AA, II, MA, MI, JJ, BB, NBJ, NAI, MJ, MB
      Integer KAIBJ, NAIBJ, INDEX
C
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "priunit.h"
Cholesky
#include "maxorb.h"
#include "ccdeco.h"
C
      LOGICAL CCSEFF
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CCSEFF = CCS .OR. (CHOINT.AND.CC2)
C
      THR1 = 1.0D-9
      THR2 = 1.0D-9
C
C------------------------------------------
C     Loop through single excitation part.
C------------------------------------------
C
      if (locdbg) then
      WRITE(LUPRI,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LUPRI,'(1X,A)')
     *     '| symmetry|  orbital index  |   Excitation Numbers'
     *     //'             |   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
     *     //'     NAIBJ    |              |'
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
      end if
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
      DO 100 ISYMA = 1,NSYM
         ISYMI = MULD2H(ISYMAI,ISYMA)
         DO 110 I = 1,NRHF(ISYMI)
            MI = IORB(ISYMI) + I
            DO 120 A=1,NVIR(ISYMA)
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
                 do ii = 1, NRHFCORE(ISYMI)
                   IF (I==IRHFCORE(II,ISYMI)) THEN
                      CAM(NAI) = zero
                   end if
                 end do
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL FLSHFO(LUPRI)
C
C--------------------------------------------
C     Loop through double excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCSEFF .OR. CCP2 )) THEN
C
      if (locdbg) then
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
      end if
C
      DO 200 ISYMAI = 1,NSYM
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
         IF (ISYMAI.lt.ISYMBJ) GO TO 200
         DO 210 ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            DO 220 ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO 230 J = 1,NRHF(ISYMJ)
                  MJ = IORB(ISYMJ) + J
                  DO 240 B = 1,NVIR(ISYMB)
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
                     DO 250 I = 1,NRHF(ISYMI)
                        MI = IORB(ISYMI) + I
                        DO 260 A = 1,NVIR(ISYMA)
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
                           IF ((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ))
     *                          GOTO 260
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSE
                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
         
                           ENDIF
                           KAIBJ = NAIBJ + NT1AM(ISYMTR)  !same vector contains singles and doubles
                                                          !skip the singles
                             do ii = 1, nrhfcore(isymi)
                              if (i==IRHFCORE(II,ISYMI)) then
                                 CAM(KAIBJ) = zero
                              end if
                             end do
                             do jj = 1, nrhfcore(isymj)
                              if (j==IRHFCORE(JJ,ISYMJ)) then
                                 CAM(KAIBJ) = zero
                              end if
                             end do
!========================================
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      ENDIF
C
 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
     *       ' | ',12x,' | ',1x, F15.9,'  |')
 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | ',I12,' | ',1x,F15.9,'  |')

      RETURN
      END
C  /* Deck cc_cvs_interface */
      SUBROUTINE CC_cvs_INTERFACE(MSYM)
C
C PURPOSE:
C  interface for transfer of CVS info module
C  Sonia, 2015
#include "implicit.h"
#include "priunit.h"
#include "ccexcicvs.h"
#include "ccxscvs.h"
C
      integer MSYM

      NRHFCORE = NXCORE
      IRHFCORE = IXCORE
      LCVSEXCI = LXSCVS
      LRMCORE  = LXRMCORE

      RETURN
      END

