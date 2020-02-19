
      SUBROUTINE SO_ONEAVE(SPNDSO,WORK,LWORK)
      
      use so_info, only: sop_dp, sop_num_models, sop_dens_label,
     &                   so_get_active_models                        

      implicit none

#include "soppinf.h"
#include "ccsdsym.h"
#include "ccorb.h"
C  MXCOOR       
#include "mxcent.h"
C  DODSO - Should we do DSO term
#include "spnout.h"
C  IRAT - number of intergers per (DP) floats
#include "iratdef.h"
C  NUCDEP
#include "nuclei.h"
C  MPQUAD
#include "magone.h" 
!     Arguments
      real(sop_dp), intent(out) :: SPNDSO(*)
      real(sop_dp), intent(inout) :: WORK(LWORK)
      integer, intent(in) :: lwork
!     Parameters
      character(len=*), parameter :: myname = 'SO_ONEAVE'
      integer, parameter :: IPRINT = 0
!     Local variables
      integer :: LCMO, LDENSIJ, LDENSAB, LDENSAI, LAODENS
      integer :: KCMO, KDENSIJ, KDENSAB, KDENSAI, KAODENS, KDENMAT,
     &           KINTAD, KINTRP, KEXPVA, KPATOM
      integer :: LWORK3, LWORK4, NCOMP, KCOMP
      integer :: KEND1, KEND2, KEND3, KEND4
      character(len=4) :: key
      logical :: active_models(sop_num_models)
      character(len=8) :: LABINT(3*MXCOOR) !Mirrors definition in magdr1
      real(sop_dp) :: DUMMY

      CALL QENTER(myname)

      if (.not.DODSO) return
      
      LCMO = NLAMDT ! From which block? 
      LDENSIJ = NIJDEN(1)
      LDENSAB = NABDEN(1)
      LDENSAI = NAIDEN(1)
      LAODENS = NNBASX
C
C     Memory for MO density and AO-MO matrix
      KCMO = 1
      KDENSIJ = KCMO + LCMO
      KDENSAB = KDENSIJ + LDENSIJ
      KDENSAI = KDENSAB + LDENSAB
      KEND1 = KDENSAI + LDENSAI
C  
C     Memory for final AO density
      KDENMAT = 1
      KEND2 = KDENMAT + LAODENS
C
C     Memory for symmetry-adapted AO density
      KAODENS = MAX(KEND1,KEND2)
      KEND3 = KAODENS + LAODENS
      LWORK3 = LWORK - KEND3
      IF (LWORK3 .LT. 0) CALL QUIT('Insuffient memory in '//myname)
C
C------------------------
C     Get MO coefficients
C------------------------      
C
      CALL SO_GETMO(WORK(KCMO),LCMO,WORK(KEND3),LWORK3)
C
C     Find only last method!
      call so_get_active_models(active_models) 
      do I = 1, sop_num_models
         if (active_models(i)) key = sop_dens_label(i)
      end do
C
C-----------------------------------------
C     Transform density matrix to AO basis
C-----------------------------------------
C
      if (key.ne.'NONE') then
         call so_read_dens(key,work(kdensij),ldensij,
     &                     work(kdensab),ldensab,
     &                     work(kdensai),ldensai)
         call so_pkaodens(WORK(KAODENS),LAODENS,WORK(KDENSIJ),LDENSIJ,
     &                    WORK(KDENSAB),LDENSAB,WORK(KDENSAI),LDENSAI,
     &                    WORK(KCMO),LCMO,WORK(KEND3),LWORK3)

         call dsym1(work(kdenmat),work(kend3),work(kaodens),
     &              work(kend3),.FALSE.,NBAST,0)

      else
         call quit('RPA expectation-values not suported in AOSOPPA'//
     &             'module')
      end if
C
C--------------------------------------
C        Calculate DSO term
C--------------------------------------
C
      IF (DODSO) THEN
            KCOMP = 0
            KPATOM = 1
            NCOMP = (3*NUCDEP*(NUCDEP+1)/2) ! Maybe lower in some cases?
            KEXPVA = KEND2
            KINTRP = KEXPVA + NCOMP !<- Which needs to be?
            KINTAD = KINTRP + (9*MXCENT**2+1)/IRAT
            KEND4 = KINTAD +(9*MXCENT**2+1)/IRAT
            LWORK4 = LWORK - KEND4
            CALL GET1IN(DUMMY,'DSO    ',KCOMP,WORK(KEND4),LWORK4,LABINT,
     &                  WORK(KINTRP),WORK(KINTAD),MPQUAD,
     &                  .FALSE.,KPATOM,.TRUE.,WORK(KEXPVA),
     &                  .TRUE.,WORK(KDENMAT),IPRINT)

            CALL MAGAVE(SPNDSO,WORK(KEXPVA),KCOMP,WORK(KEND4),LWORK4,
     &                  IPRINT,'DSO    ',LABINT,1.0D0)
                        
      END IF

      CALL QEXIT(myname)

      END
