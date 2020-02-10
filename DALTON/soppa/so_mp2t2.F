C  /* Deck so_mp2t2 */
      SUBROUTINE SO_MP2T2(FOCKD,LFOCKD,T2AM,LT2AM,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Gets (ia|jb) integrals in T2AM and
C              calculates the T² MP2 amplitudes.
C
C     Written by Lilli Irene Ør Kristensen
C
      use so_info, only: sop_dp
      implicit none
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"

      real(sop_dp), intent(in) :: fockd(lfockd)
      real(sop_dp), intent(inout) :: t2am(lt2am), work(lwork)
      integer, intent(in) :: lfockd, lt2am, lwork

      integer :: isymai, isymbj, isymi, isyma
      integer :: koffi, koffa, ioff, ioff2
      integer :: nai, nbj, naibj
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_MP2T2')
C
      DO ISYMAI = 1,NSYM
         ISYMBJ = ISYMAI

         ! Establish F_aa - F_ii
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMAI)
            DO I = 1, NRHF(ISYMI)
               KOFFI = IRHF(ISYMI) + I
               DO A = 1, NVIR(ISYMA)
                  NAI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                  KOFFA = IVIR(ISYMA) + A
                  WORK(NAI) =  FOCKD(KOFFA) - FOCKD(KOFFI)
               END DO
            END DO
         END DO
C
         IOFF = IT2AM(ISYMAI,ISYMAI)
         DO NAI = 1, NT1AM(ISYMAI)
            IOFF2 = IOFF + (NAI*(NAI-1))/2
            DO NBJ = 1, NAI
               NAIBJ = IOFF2 + NBJ
               T2AM(NAIBJ) = -T2AM(NAIBJ)/(WORK(NAI)+WORK(NBJ))
            END DO
         END DO

      END DO
C
C Printing the amplitudes
C      IF (IPRINT .GT. 10) THEN
C         CALL AROUND('MP2 amplitudes')
C         DO 250 ISYMBJ = 1,NSYM
C            ISYMAI = ISYMBJ
C            IOFF   = IT2AM(ISYMAI,ISYMBJ) + 1
C            NTOTAI = NT1AM(ISYMAI)
C            CALL OUTPUT(T2AM(KOFF),1,NTOTAI,1,1,NTOTAI,1,1,LUPRI)
C250      CONTINUE
C      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_MP2T2')
C
      RETURN
      END
