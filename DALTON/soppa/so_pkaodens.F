      SUBROUTINE SO_PKAODENS(PKDENS,LPKDENS,DENSIJ,LDENSIJ,
     &                       DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                       CMO,LCMO,WORK,LWORK)
C
C     Rasmus Faber, 2017
C    
C     Calculate the "SOPPA" density matrix in AO-basis, 
C     
      use so_info, only : sop_dp
      implicit none
C
C     Arguments
      real(sop_dp), intent(out) :: PKDENS(LPKDENS)
      real(sop_dp), intent(in) :: DENSIJ(LDENSIJ), DENSAB(LDENSAB),
     &                            DENSAI(LDENSAI), CMO(LCMO)
      real(sop_dp), intent(inout) :: WORK(LWORK)
      integer, intent(in) :: LPKDENS, LDENSIJ, LDENSAB, LDENSAI, 
     &                       LCMO, LWORK
C
C     Parameters 
      real(sop_dp), parameter :: one = 1.0D0, zero = 0.0D0, two = 2.0D0,
     &                           four = 4.0D0

C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "soppinf.h"
C
C     Local variables
      integer :: KTMP, KEND
      integer :: LTMP
      integer :: KAB, KDT, KMO, KOFFMO, KOFFAB, KOFFAI, KOFFIJ
      integer :: NRHFI, NVIRA, NBASA
      integer :: isym, IOFFOUT, AL

      IOFFOUT = 0

      CALL DZERO(PKDENS,LPKDENS)

      DO ISYM = 1, NSYM
C
         LTMP = NBAS(ISYM)*MAX(NVIR(ISYM),NRHF(ISYM))
         KTMP = 1
         KEND = KTMP + LTMP
         KOFFAI = IT1AM(ISYM,ISYM) + 1
         KOFFIJ = IIJDEN(ISYM,ISYM) + 1 
         KOFFAB = IABDEN(ISYM,ISYM) + 1
         KOFFMO = ILMVIR(ISYM) + 1

         NVIRA = MAX(NVIR(ISYM),1)
         NRHFI = MAX(NRHF(ISYM),1)
         NBASA = MAX(NBAS(ISYM),1)
C
C----------------------------------------
C        D(a,i) * U(beta,a) => DT(i,beta)
C----------------------------------------
         CALL DGEMM('T','T',NRHF(ISYM),NBAS(ISYM),NVIR(ISYM),FOUR,
     &              DENSAI(KOFFAI),NVIRA,
     &              CMO(KOFFMO),NBASA,
     &              ZERO,WORK(KTMP),NRHFI)
C
C------------------------------------------
C        D(j,i) * U(beta,j) => DT(i,beta)
C------------------------------------------
         KOFFMO = ILMRHF(ISYM) + 1
         CALL DGEMM('T','T',NRHF(ISYM),NBAS(ISYM),NRHF(ISYM),TWO,
     &              DENSIJ(KOFFIJ),NRHFI,
     &              CMO(KOFFMO),NBASA,
     &              ONE,WORK(KTMP),NRHFI)
C---------------------------------------------
C        U(beta,j) => DT(i,beta) : RHF density
C---------------------------------------------
          DO B = 1, NBAS(ISYM)
            DO I = 1, NRHF(ISYM)
               KMO = ILMRHF(ISYM) + NBAS(ISYM)*(I-1) + B
               KDT = KTMP-1 + NRHF(ISYM)*(B-1) + I
               WORK(KDT) = WORK(KDT) + TWO*CMO(KMO)
            END DO
          END DO

C         CALL DAXPY(NBAS(ISYM)*NRHF(ISYM),TWO,CMO(KOFFMO),1,
C     &              WORK(KTMP),1)
C
C----------------------------------------------------------------
C        Transform the second index of the D(i,beta) intermediate 
C----------------------------------------------------------------
C        DT(i,beta) * U(alpha,i) + DT(i,alpha) * U(beta,i) => out
C
         DO B = 1, NBAS(ISYM)
C
            DO I = 1, NRHF(ISYM)
               KOFFMO = ILMRHF(ISYM) + NBAS(ISYM)*(I-1)
               KDT = NRHF(ISYM)*(B-1) + I-1
               DO A = 1, B
                  KMO = KOFFMO + A
                  KAB = IOFFOUT + B*(B-1)/2 + A
                  PKDENS(KAB) = PKDENS(KAB) + WORK(KTMP+KDT)*CMO(KMO)
               END DO
            END DO

            DO A = 1, B-1
               KAB = IOFFOUT + B*(B-1)/2 + A
               
               DO I = 1, NRHF(ISYM)
                  KMO = ILMRHF(ISYM) + NBAS(ISYM)*(I-1) + B
                  KDT = NRHF(ISYM)*(A-1) + I-1
                  PKDENS(KAB) = PKDENS(KAB) + WORK(KTMP+KDT)*CMO(KMO)
               END DO
            END DO
         END DO
C
C----------------------------------------
C        D(b,a) * U(beta,b) => Dt(a,beta)
C----------------------------------------
C
         KOFFMO = ILMVIR(ISYM) + 1
         CALL DGEMM('T','T',NVIR(ISYM),NBAS(ISYM),NVIR(ISYM),TWO,
     &              DENSAB(KOFFAB),NVIRA,
     &              CMO(KOFFMO),NBASA,
     &              ZERO,WORK(KTMP),NVIRA)
C
C------------------------------------------
C        Dt(a,beta) * U(alpha,a) => out(alpha,beta)
C------------------------------------------
         DO B = 1, NBAS(ISYM)
C
            DO AL = 1, B-1
               KAB = IOFFOUT + B*(B-1)/2 + AL 
               DO A = 1, NVIR(ISYM)
                  KDT = NVIR(ISYM)*(B-1) + A-1
                  KMO = ILMVIR(ISYM) + NBAS(ISYM)*(A-1) + AL
                  PKDENS(KAB) = PKDENS(KAB) +TWO*WORK(KTMP+KDT)*CMO(KMO)
               END DO
            END DO
C        Handle alpha == beta
            KAB = IOFFOUT + B*(B+1)/2
            DO A = 1, NVIR(ISYM)
               KDT = NVIR(ISYM)*(B-1) + A-1
               KMO = ILMVIR(ISYM) + NBAS(ISYM)*(A-1) + B
               PKDENS(KAB) = PKDENS(KAB) + WORK(KTMP+KDT)*CMO(KMO)
            END DO
         END DO

C         WRITE (LUPRI,'(1X,A,I5)') ' Symmetry', ISYM
C         CALL HEADER('Total density matrix (SO basis)',-1)
C         CALL OUTPAK(PKDENS(IOFFOUT+1),NBASA,1,LUPRI)

         IOFFOUT = IOFFOUT + (NBAS(ISYM)*(NBAS(ISYM)+1))/2

      END DO

      RETURN

      END
