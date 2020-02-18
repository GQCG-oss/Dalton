C
C  /* Deck so_lrnsl */
      SUBROUTINE SO_LRNSL(MODEL,ISYM,PROPAR,DELETE_VECTORS,
     &                    FREQ,NFREQ,LABELS,
     &                    SOLV_LABELS, NLAB,
     &                    DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                    WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer: 5.12.2003
C     Rasmus Faber: Nov. 2015 --> Implemented pert. dens. approach.
C
C     PURPOSE: Calculates the frequency dependent linear response
C              properties from the perturbed density matrices and
C              appropriate property integrals with the atomic
C              integral direct SOPPA program.
C
      use so_info, only: fn_rdens, sop_stat_trh, sop_dp,
     &                   sop_use_seller, so_has_doubles
      use so_data, only: fileinf
C
C#include "implicit.h"
      implicit none
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C REP <- Symbols of the ireducible representations.
C#include "pgroup.h"
#include "iratdef.h"
C

C
      integer, intent(in) :: isym, ! symmetry of this batch of perturbations
     &       LDENSIJ, LDENSAB, LDENSAI, LWORK, ! array lengths
     &       NFREQ, NLAB
      REAL(sop_dp),INTENT(IN) :: DENSIJ(LDENSIJ), DENSAB(LDENSAB),
     &                     DENSAI(LDENSAI), FREQ(NFREQ)
      REAL(sop_dp),INTENT(INOUT) :: WORK(LWORK)
      REAL(sop_dp),intent(inout) :: PROPAR(NFREQ,NLAB,NLAB)
      CHARACTER(LEN=8), intent(in) :: LABELS(NLAB) ! Operator lables
      CHARACTER(len=5), intent(in) :: MODEL
      LOGICAL, INTENT(IN) :: DELETE_VECTORS, ! Whether to delete solution and residual vectors
     &                       SOLV_LABELS(NLAB) ! Whether response equations have been solved
C
      CHARACTER*8 LABEL1, LABEL2
      CHARACTER*8 RTNLBL(2)
      REAL(sop_dp)  ::  PROPVAL, SELLCOR
      INTEGER :: LPRP1, LWORK1, LPDENSIJ, LPDENSAB, LPDENSAI, LPDENSTOT,
     &           LWORK2, LWORK3, LSOLV, LRESI
      INTEGER :: KPRP1, KEND1, KPDENSIJ, KPDENSAB, KPDENSAI, KEND2,
     &           KEND3, KSOLV, KRESI
      INTEGER :: IDUMMY, ILAB1, ILAB2, IFREQ, LURDENS
      INTEGER :: NSTART
C
      LOGICAL   IMAGPROP, doubles, static, OLDDX, SOLVED1
      REAL(sop_dp), PARAMETER :: DP5=0.5D0
C
      CHARACTER(len=*), parameter :: myname = 'SO_LRNSL'
      real(sop_dp) :: DDOT
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER(myname)
C
      PROPAR = 0.0D0
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      KEND1 = 1
C
C----------------------------------
C     Space for MO property matrix.
C----------------------------------
C
      LPRP1 = N2BST(ISYM)
      KPRP1 = KEND1
      KEND1 = KPRP1 + LPRP1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX (myname//'.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT(myname//'.1',' ',KEND1,LWORK)
C
C-----------------------------------------
C     Open files with perturbed densities.
C-----------------------------------------
C
      LPDENSIJ = NIJDEN(ISYM)
      LPDENSAI = NAIDEN(ISYM)
      IF (MODEL.EQ.'AORPA') THEN
         LPDENSAB = 0
         LPDENSTOT = LPDENSAI
      ELSE
         LPDENSAB = NABDEN(ISYM)
         LPDENSTOT = LPDENSIJ + LPDENSAI + LPDENSAB
      ENDIF
      LURDENS = -1
      CALL GPOPEN(LURDENS, FN_RDENS, 'OLD','DIRECT','UNFORMATTED',
     &            IRAT*LPDENSTOT,OLDDX)
C
      doubles = so_has_doubles(model)
C
C=============================================
C     Loop over the first property operators.
C=============================================
C
      DO 200 ILAB1 = 1, NLAB
C
C---------------------------------------------------
C        Find label and symmetry of second operator.
C---------------------------------------------------
C
         LABEL1 = LABELS(ILAB1)
C        Check if we have solved wrt LABEL1
         SOLVED1 = SOLV_LABELS(ILAB1)
         IF (SOLVED1) THEN
            NSTART = ILAB1
         ELSE
            NSTART = 1
         END IF
C
C--------------------------------------------------
C           Get the property integrals in MO basis.
C--------------------------------------------------
C
         CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL1,ISYM,
     &                  RTNLBL,WORK(KEND1),LWORK1)
C
         IMAGPROP = RTNLBL(2).EQ.'ANTISYMM'
C
C  These allocations as such doesn't change...

         KPDENSIJ = KEND1
         KPDENSAB = KPDENSIJ + LPDENSIJ
         KPDENSAI = KPDENSAB + LPDENSAB
         KEND2    = KPDENSAI + LPDENSAI
         LWORK2   = LWORK - KEND2
C
C
         CALL SO_MEMMAX (myname//'.2',LWORK2)
         IF (LWORK2 .LT. 0)
     &            CALL STOPIT(myname//'.2',' ',KEND2,LWORK)
C
C===============================================
C           Form second order properties SNDPRP.
C===============================================
C
C
         DO 100 IFREQ = 1, NFREQ !

            STATIC = ABS(FREQ(IFREQ)).LT.sop_stat_trh
C-----------------------------------------------------------
C           If using Seller's, get solution vector of LABEL1
C-----------------------------------------------------------
            IF (SOP_USE_SELLER.AND.SOLVED1) THEN
               LRESI = NT1AM(ISYM)
               IF (DOUBLES) LRESI = LRESI + N2P2HOP(ISYM)
               IF (.NOT.STATIC) LRESI = LRESI*2
               LSOLV = LRESI
               KSOLV = KEND2
               KRESI = KSOLV + LSOLV
               KEND3 = KRESI + LRESI
               LWORK3 = LWORK - KEND3
               IF (LWORK3 .LT. 0)
     &            CALL STOPIT(myname//'.3',' ',KEND3,LWORK)
               CALL SO_GETSOLV(WORK(KSOLV),LSOLV,LABEL1,FREQ(IFREQ),1,
     &                         STATIC,DOUBLES,ISYM)
            ENDIF

            DO ILAB2 = NSTART, NLAB
               LABEL2 = LABELS(ILAB2)
               ! Only do the following, if solution vector has been
               ! found
               IF (.NOT.SOLV_LABELS(ILAB2)) CYCLE
C
C-----------------------------------------------------------
C        Get the perturbed density matrix from file
C-----------------------------------------------------------
C
               LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
               IF (MODEL.EQ.'AORPA') THEN
                  CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYM,
     &                          LABEL2,FREQ(IFREQ),LURDENS)
                  CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
               ELSE
                  CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYM,
     &                          LABEL2,FREQ(IFREQ),LURDENS)
               ENDIF
C
C---------------------------------------------------------------------
C        Calculate second order properties PROPVAL.
C---------------------------------------------------------------------
C
               CALL SO_PROPMO(ISYM,PROPVAL,
     &                        MODEL.NE.'AORPA',IMAGPROP,
     &                        WORK(KPRP1),LPRP1,
     &                        WORK(KPDENSIJ),LPDENSIJ,
     &                        WORK(KPDENSAB),LPDENSAB,
     &                        WORK(KPDENSAI),LPDENSAI)
C
               IF (SOP_USE_SELLER.AND.SOLVED1.AND.(ILAB1.NE.ILAB2)) THEN
                  CALL SO_GETSOLV(WORK(KRESI),LRESI,LABEL2,FREQ(IFREQ),
     &                            2,STATIC,DOUBLES,ISYM)
                  SELLCOR = DDOT(LRESI,WORK(KRESI),1,WORK(KSOLV),1)
                  IF (STATIC) SELLCOR = SELLCOR*2
                  PROPVAL = PROPVAL - SELLCOR
               ENDIF
C
               PROPAR(IFREQ,ILAB1,ILAB2) = PROPVAL
               PROPAR(IFREQ,ILAB2,ILAB1) = PROPVAL
C
            END DO
C
  100    CONTINUE
C
  200 CONTINUE
C
C---------------------------------------------------------
C     Delete files for this symmetry and dump fileinf list
C---------------------------------------------------------
C
      IF (DELETE_VECTORS) THEN
         CALL GPCLOSE(LURDENS,'DELETE')
         STATIC = .TRUE.
         DO IFREQ = 1, NFREQ
            STATIC = STATIC.AND.(abs(freq(ifreq)).lt.sop_stat_trh)
         END DO
         IF (SOP_USE_SELLER)CALL SO_DELVEC(STATIC,DOUBLES,ISYM)
         call fileinf%empty
      ELSE
         CALL GPCLOSE(LURDENS,'KEEP')
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_LRNSL')
      RETURN
      END

C
      SUBROUTINE SO_RSPOUT(PROPAR,ISYM,FREQ,NFREQ,LABELS,NLAB)

      USE SO_INFO, ONLY: SOP_DP, SOP_STAT_TRH
      IMPLICIT NONE
C
C Get LUPRI, output unit (probably 6)
#include "priunit.h"
C From codata we need varius conversion factors, XT*.
#include "codata.h"
C inflr.h below requires MAXLBL from rspprp.h
#include "rspprp.h"
C We need
C LBLLR <- Label of operators
C NGPLR <- Number of operators of each symmetry
#include "inflr.h"
C TRIPLET <- Triplet property flag
#include "soppinf.h"

      REAL(SOP_DP), INTENT(IN) :: PROPAR(NFREQ,NLAB,NLAB),
     &                            FREQ(NFREQ)
      INTEGER, INTENT(IN) :: ISYM, NFREQ, NLAB
      CHARACTER(LEN=8) :: LABELS(NLAB)

      INTEGER :: ILAB1, ILAB2, IFREQ

      CALL HEADER('Final output of second order properties from'//
     &            ' linear response',-1)
      IF (TRIPLET) THEN
         WRITE (LUPRI,'(/A)') '@ Spin symmetry of operators: triplet'
      ELSE
         WRITE (LUPRI,'(/A)') '@ Spin symmetry of operators: singlet'
      END IF
      WRITE (LUPRI,'(/A/A)')
     &   ' Note that minus the linear response function:'//
     &   ' - << A; B >>(omega) is printed.',
     &   ' The results are of quadratic accuracy using Sellers formula.'

      DO IFREQ = 1, NFREQ
         IF ( ABS(FREQ(IFREQ)) .LT. sop_stat_trh) THEN
            WRITE(LUPRI,'(/A/)')
     *      '@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES'
         ELSE
            WRITE(LUPRI,'(/A/,5(/A,1P,D15.7),/)')
     * '@ FREQUENCY DEPENDENT SECOND ORDER PROPERTIES WITH FREQUENCY :',
     *      '@    a.u.:',FREQ(IFREQ),
     *      '@    cm-1:',XTKAYS*FREQ(IFREQ),
     *      '@    eV  :',XTEV*FREQ(IFREQ),
     *      '@  kJ/mol:',XKJMOL*FREQ(IFREQ),
     *      '@    nm  :',XTNM/FREQ(IFREQ)
         END IF
         DO ILAB1 = 1, NLAB
            DO ILAB2 = ILAB1, NLAB
               WRITE(LUPRI,'(5A,1P,D20.12)')
     *         '@ -<< ',LABELS(ILAB1),' ; ',LABELS(ILAB2),
     &         ' >> =',PROPAR(IFREQ,ILAB1,ILAB2)
            END DO
         END DO
      END DO

      END SUBROUTINE

      SUBROUTINE SO_GETSOLV(VECT,LVECT,LABEL,FREQ,IOPT,STATIC,DOUBLES,
     &                      ISYM)
C
C     Gets a converged solution vector or corresponding residual vector
C
      use so_info, only: sop_dp

      implicit none
#include "soppinf.h"
#include "ccsdsym.h"
#include "iratdef.h"
C     Arguments
      INTEGER, INTENT(IN)       :: lvect ! length of the vector
      real(sop_dp), intent(out) :: vect(lvect) ! the vector to read
      character(len=8), intent(in) :: LABEL
      real(sop_dp), intent(in)  :: freq ! frequency of pt.
      integer, intent(in)       :: iopt,! 1 or 2 to read solution or
                                        ! residual
     &                             isym ! symmetry of pt.
      logical, intent(in)  :: static, ! only read E part if static
     &                        doubles ! include doubles excitations
C
C     Variables
      character(len=5) :: froot
      integer          :: k1e, k1d, k2e, k2d, kend
C
C
      IF (IOPT.eq.1) then
         froot = fnsv1e(1:5)
      elseif (iopt.eq.2) then
         froot = fnrv1e(1:5)
      endif
C
C     Find positions in VECT
      k1e = 1
      kend = k1e + nt1am(isym)
      call readavector(froot//'1E', isym, label,
     &                 freq, vect(k1e), nt1am(isym))
      if (.NOT.STATIC) THEN
         k1d = kend
         kend = k1d + nt1am(isym)
         call readavector(froot//'1D',isym,
     &                    label,freq,vect(k1d),nt1am(isym))
      end if
      if (doubles) then
         k2e = kend
         kend = k2e + n2p2hop(isym)
         call readavector(froot//'2E',isym,label,freq,vect(k2e),
     &                     n2p2hop(isym))
         if (.not.static) then
            k2d = kend
            kend = k2d + n2p2hop(isym)
            call readavector(froot//'2D',isym,label,freq,vect(k2d),
     &                     n2p2hop(isym))
         end if
      end if
      !if (kend .ne. lvect+1) print *, 'Odd vector length', kend, lvect
C
      contains

         subroutine readavector(fname,isym,label,freq,vect,length)
            character(len=7),intent(in)  :: fname
            character(len=8),intent(in)  ::label
            integer, intent(in) :: isym, length
            real(sop_dp), intent(in) :: freq
            real(sop_dp), intent(out) :: vect(length)

            integer :: lufile
            logical :: olddx

            lufile = -1
            call GPOPEN(lufile,fname,'OLD    ','DIRECT','UNFORMATTED',
     &                 irat*length,olddx)
            call SO_REAVE(vect,length,isym,label,freq,lufile)
            call gpclose(lufile,'KEEP')
         end subroutine

      END SUBROUTINE

      SUBROUTINE SO_DELVEC(STATIC,DOUBLES,ISYM)
C
C     Deletes the files containing solution and
C     residual vectors
C
      implicit none
#include "soppinf.h"
#include "ccsdsym.h"
#include "iratdef.h"



      integer, intent(in)  :: isym ! symmetry of pt.
      logical, intent(in)  :: static, ! only read E part if static
     &                        doubles ! include doubles excitations
C
      character(len=5) :: froot
      integer :: iopt

      do iopt = 1,2
         IF (IOPT.eq.1) then
            froot = fnsv1e(1:5)
         elseif (iopt.eq.2) then
            froot = fnrv1e(1:5)
         endif

         CALL DELAVECT(froot//'1E', nt1am(isym))
         if (.not.static) call delavect(froot//'1D', nt1am(isym))
         if (doubles) then
            CALL DELAVECT(froot//'2E', n2p2hop(isym))
            if (.not.static) CALL DELAVECT(froot//'2D', n2p2hop(isym))
         end if
      end do

      contains

         subroutine delavect(fname,length)
            integer, intent(in) :: length
            character(len=7),intent(in)  :: fname

            integer :: lufile
            logical :: olddx
            lufile = -1
            call GPOPEN(lufile,fname,'OLD    ','DIRECT','UNFORMATTED',
     &                 irat*length,olddx)
            call gpclose(lufile,'DELETE')
         end subroutine
      END SUBROUTINE
