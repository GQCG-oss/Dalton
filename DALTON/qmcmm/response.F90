module qmcmm_response

   implicit none

   public qmcmm_lr
   public qmcmm_qr

   private

contains

   subroutine qmcmm_lr(nosim,bovecs,cref,cmo,udv,dv,udvtr,   &
                            dvtr,evecs,work,lwork)
!
! Purpose:
!     Computes QM/NP/MM contribution to linear response vector.
!
! Input:
!   NOSIM  - Number of first Fock/Kohn-Sham matrices
!   BOVECs - Linear response vectors
!   CREF   - Reference state CI coeficients (Currently not used)
!   CMO    - Molecular orbitals
!   UDV    - Density matrix
!   DV     - Active density matrix
!   UDVTR  - Triplet density matrix
!   DVTR   - Active triplet density matrix
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array
! Output:
!   EVECS  - Response vectors (XY).
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

      use qmcmm, only: getdim_relmat, read_relmat

#include "inforb.h"
#include "infdim.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "priunit.h"
#include "qmnpmm.h"

      integer, intent(in) :: nosim
      integer, intent(in) :: lwork

      real(8) :: BOVECS(*), CMO(*), UDV(NASHDI,NASHDI)
      real(8) :: UDVTR(N2ASHX), DVTR(*), EVECS(KZYVAR,*)
      real(8) :: WORK(*), DV(*), CREF(*)
      integer :: i, j, idimension, ioff

      real(8), allocatable :: mqvec(:)
      real(8), allocatable :: fvvec(:)
      real(8), allocatable :: ucmo(:)
      real(8), allocatable :: bov(:)
      real(8), allocatable :: relmat(:)
      real(8), allocatable :: rxy(:)
      real(8), allocatable :: rxyt(:)

!     np/mm contribution to response matrix is zero for triplet
!     perturbations applied to singlet reference state
      if ((nasht == 0) .and. (trplet)) return

      allocate(ucmo(nbast*norbt))
      ucmo = 0.0d0
      allocate(bov(nosim*n2orbx))
      bov = 0.0d0
      allocate(rxy(nosim*n2orbx))
      rxy = 0.0d0
      if (trplet) then
          allocate(rxyt(nosim*n2orbx))
          rxyt = 0.0d0
      end if

      if (mqiter) then ! iterative method
         call quit('mqiter not implemented in rspqmnp.F90')
      else ! non iterative method
         idimension = getdim_relmat(.false.)

         ! Zero & unpack CMO and ZY vectors
         CALL UPKCMO(CMO,uCMO)
         IF (NOSIM.GT.0) THEN
            CALL RSPZYM(NOSIM,BOVECS,BOV)
            CALL DSCAL(NOSIM*N2ORBX,-1.0d0,BOV,1)
         END IF

         allocate(mqvec(nosim*idimension))
         mqvec = 0.0d0
         allocate(fvvec(nosim*idimension))
         fvvec = 0.0d0

         ! Determine electric field/potential vector for perturbed
         ! density matrices
         call get_fvvec(idimension=idimension,    &
                        nsim=nosim,   &
                        udv=udv,      &
                        cmo=ucmo,     &
                        work=work,    &
                        lwork=lwork,  &
                        fvvec1=fvvec, &
                        udvtr=udvtr,  &
                        bovecs=bov)

         ! Allocate and compute Relay matrix
         idimension = getdim_relmat(.true.)
         allocate(relmat(idimension))
         CALL READ_RELMAT(RELMAT)
         DO I=1,NOSIM
            IOFF = (I-1)*idimension
            CALL DGEMV('N',idimension,idimension,1.0d0,RELMAT,idimension,              &
                       FVVEC(IOFF + 1),1,0.0d0,MQVEC(IOFF + 1),1)
           if (iprtlvl > 14) then
              write(lupri, '(/,2x,a,i0)') &
                  '*** Computed MQ vector start 1st-order density ', i
              do j = 1, idimension
                 write(lupri, '(i8, f18.8)') j, MQVEC(IOFF + j)
              end do
              write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
           end if
         END DO
         deallocate(relmat)

         ! compute xy contributions from induced dipoles and charges
         if (trplet) then
            call get_xyvec(ucmo,  &
                           idimension,  &
                           nosim, &
                           rxyt,  &
                           work,  &
                           lwork, &
                           mqvec)
         else
            call get_xyvec(ucmo,  &
                           idimension,  &
                           nosim, &
                           rxy,   &
                           work,  &
                           lwork, &
                           mqvec)
         end if
      END IF

      ! add qm/np/mm contributions to transformed resp. vectors
      if (trplet) then
         ! fixme: i wonder whether this is right, rxy is then never computed
         call slvsor(.true.,.false.,nosim,udvtr,evecs(1,1),rxy)
         call slvsor(.true.,.true.,nosim,udv,evecs(1,1),rxyt)
      else
         call slvsor(.true.,.true.,nosim,udv,evecs(1,1),rxy)
      endif

      if (allocated(mqvec)) deallocate(mqvec)
      if (allocated(fvvec)) deallocate(fvvec)
      if (allocated(ucmo))  deallocate(ucmo)
      if (allocated(bov))   deallocate(bov)
      if (allocated(rxy))   deallocate(rxy)
      if (allocated(rxyt))  deallocate(rxyt)

   end subroutine

   subroutine qmcmm_qr(vec1,vec2,etrs,xindx,zym1,zym2,              &
                           udv,work,lwork,kzyvr,kzyv1,kzyv2,            &
                           igrsym,isymv1,isymv2,cmo,mjwop,              &
                           ispin0,ispin1,ispin2)

      use qmcmm, only: getdim_relmat, read_relmat

#include "maxorb.h"
#include "inforb.h"
#include "infdim.h"
#include "infinp.h"
#include "infvar.h"
#include "infrsp.h"
#include "infpri.h"
#include "rspprp.h"
#include "infcr.h"
#include "inftap.h"
#include "qrinf.h"
#include "mxcent.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "infpar.h"
#include "qmnpmm.h"

      integer :: kzyvr
      integer :: kzyv1
      integer :: kzyv2
      integer :: lwork,igrsym,isymv1,isymv2,ispin0,ispin1,ispin2
      real(8) :: etrs(kzyvr),xindx(*)
      real(8) :: udv(nashdi,nashdi)
      real(8) :: zym1(*),zym2(*),work(lwork),cmo(*)
      real(8) :: vec1(kzyv1),vec2(kzyv2)
      integer :: mjwop(2,maxwop,8)
      integer :: i, j, ioff
      logical   lcon, lorb, lref
      integer :: idimension
      integer :: isymt, isymv, isymst, jspin, nsim
      integer :: nzyvec, nzcvec
      integer :: ISYMDN, idum
      real(8) :: ovlap

      real(8), allocatable :: cref(:)
      real(8), allocatable :: tres(:)
      real(8), allocatable :: ucmo(:)
      real(8), allocatable :: tlma(:)
      real(8), allocatable :: tlmb(:)
      real(8), allocatable :: trmo(:)
      real(8), allocatable :: utr(:)
      real(8), allocatable :: mqvec1(:)
      real(8), allocatable :: mqvec2(:)
      real(8), allocatable :: fvvec1(:)
      real(8), allocatable :: fvvec2(:)
      real(8), allocatable :: relmat(:)

!     Allocate arrays for response
      allocate(cref(ncref))
      allocate(tres(n2orbx))
      allocate(ucmo(norbt*nbast))
      allocate(tlma(n2orbx))
      allocate(tlmb(n2orbx))
      allocate(trmo(nnorbx))
      allocate(utr(n2orbx))
!     Initialize allocated arrays
      cref = 0.0d0
      tres = 0.0d0
      ucmo = 0.0d0
      tlma = 0.0d0
      tlmb = 0.0d0
      trmo = 0.0d0
      utr = 0.0d0
!     Reset symmetry variables
      NSIM  = 1
      ISYMT = 1
!     Get the reference state
      CALL GETREF(CREF,MZCONF(1))
!     Unpack the response vectors
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
!     Unpack symmetry blocked CMO
      CALL UPKCMO(CMO,UCMO)
!     Non-iterative method
      IF (.NOT.MQITER) THEN
        idimension = getdim_relmat(.false.)
!       Allocate FV and MQ vectors for all perturbed densities
        allocate(fvvec1(nsim*idimension))
        allocate(fvvec2(nsim*idimension))

!       compute fv vectors for second order pertubed densities
        call get_fvvec(idimension=idimension,     &
                       nsim=nsim,     &
                       udv=udv,       &
                       cmo=ucmo,      &
                       work=work,     &
                       lwork=lwork,   &
                       fvvec1=fvvec1, &
                       fvvec2=fvvec2, &
                       isymt=isymt,   &
                       isymv1=isymv1, &
                       isymv2=isymv2, &
                       zym1=zym1,     &
                       zym2=zym2)

!       Allocate and compute Relay matrix
        idimension = getdim_relmat(.true.)
        allocate(relmat(idimension))
        CALL READ_RELMAT(RELMAT)
!       Determine induced induced dipoles and charges
        allocate(mqvec1(nsim*idimension))
        allocate(mqvec2(nsim*idimension))
        mqvec1 = 0.0d0
        mqvec2 = 0.0d0
        DO I=1,NSIM
           IOFF = (I-1)*idimension
           CALL DGEMV('N',idimension,idimension,1.0d0,RELMAT,idimension,              &
     &                FVVEC1(IOFF + 1),1,0.0d0,MQVEC1(IOFF + 1),1)
           CALL DGEMV('N',idimension,idimension,1.0d0,RELMAT,idimension,              &
     &                FVVEC2(IOFF + 1),1,0.0d0,MQVEC2(IOFF + 1),1)
          if (iprtlvl > 14) then
             write(lupri, '(/,2x,a,i0)') &
                 '*** Computed MQ vector start v1 1st-order density ', i
             do j = 1, idimension
                write(lupri, '(i8, f18.8)') j, MQVEC1(IOFF+j)
             end do
             write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
             write(lupri, '(/,2x,a,i0)') &
                 '*** Computed MQ vector start v2 1st-order density ', i
             do j = 1, idimension
                write(lupri, '(i8, f18.8)') j, MQVEC2(IOFF+j)
             end do
             write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
          end if
        END DO
        deallocate(mqvec1)
        deallocate(mqvec2)
        deallocate(relmat)

        ! determine mm region contribution to qm region potential from
        ! second order density
        call get_xyvec(ucmo,   &
                       idimension,   &
                       nsim,   &
                       tres,   &
                       work,   &
                       lwork,  &
                       fvvec1, &
                       fvvec2, &
                       isymt,  &
                       isymv2, &
                       zym2)

      ELSE
!       FIX ME: ITERATIVE METHOD
      END IF
!     Set up paramterers for quadratic response gradient formation
      ISYMDN = 1
      OVLAP  = 1.0d0
      JSPIN  = 0
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = NCREF
      NZCVEC = NCREF
!     Compute gradient
      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,        &
     &            CREF,NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,           &
     &            TRES,XINDX,MJWOP,WORK,lwork,            &
     &            LORB,LCON,LREF)
!

      deallocate(cref)
      deallocate(tres)
      deallocate(ucmo)
      deallocate(tlma)
      deallocate(tlmb)
      deallocate(trmo)
      deallocate(utr)

      end subroutine







   ! computes electric field/potential vector generated by first or second
   ! order perturbed density matrix
   subroutine get_fvvec(idimension,   &
                        nsim,   &
                        udv,    &
                        cmo,    &
                        work,   &
                        lwork,  &
                        fvvec1, &
                        fvvec2, &
                        udvtr,  &
                        bovecs, &
                        isymt,  &
                        isymv1, &
                        isymv2, &
                        zym1,   &
                        zym2)

      ! size of electric field/potential vector
      integer, intent(in)              :: idimension
      ! number of perturbed density matrices
      integer, intent(in)              :: nsim
      real(8), intent(in)              :: udv(nashdi, nashdi)
      real(8), intent(in)              :: cmo(*)
      real(8), intent(inout)           :: work(*)
      integer, intent(in)              :: lwork
      ! electric field/potential vector at np/mm centers
      real(8), intent(inout)           :: fvvec1(*)
      real(8), intent(inout), optional :: fvvec2(*)
      real(8), intent(in),    optional :: udvtr(n2ashx)
      real(8), intent(in),    optional :: bovecs(*)
      integer, intent(in),    optional :: isymt
      integer, intent(in),    optional :: isymv1
      integer, intent(in),    optional :: isymv2
      real(8), intent(in),    optional :: zym1(*)
      real(8), intent(in),    optional :: zym2(*)

#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "infdim.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "infrsp.h"

      logical :: tofile
      logical :: trimat
      logical :: exp1vl
      real(8) :: dipole_origin_save(3)
      integer :: i
      integer :: j
      integer :: ioff
      integer :: joff
      integer :: istart
      integer :: iblk
      integer :: nocomp
      integer :: kpatom
      integer :: isimoff
      integer :: ixyz
      real(8) :: f1val, f2val
      real(8) :: fact
      real(8) :: tac

      real(8), external :: slvtlm
      real(8), external :: slvqlm

      real(8), allocatable :: intao(:)
      real(8), allocatable :: trmo(:)
      real(8), allocatable :: utr(:)
      real(8), allocatable :: tlma(:)
      real(8), allocatable :: tlmb(:)
      real(8), allocatable :: utrx(:)
      real(8), allocatable :: urxac(:)
      integer, allocatable :: intrep(:)
      integer, allocatable :: intadr(:)
      character(8), allocatable :: labint(:)

      ! save origin coordinates
      dipole_origin_save = diporg

      call dzero(fvvec1, idimension*nsim)
      if (present(fvvec2)) then
         call dzero(fvvec2, idimension*nsim)
      end if

      allocate(intao(3*nnbasx))
      allocate(trmo(nnorbx))
      allocate(utr(n2orbx))
      if (present(fvvec2)) then
         allocate(tlma(n2orbx))
         allocate(tlmb(n2orbx))
      else
         allocate(utrx(n2orbx))
         allocate(urxac(n2ashx))
      end if
      allocate(intrep(9*mxcent))
      allocate(intadr(9*mxcent))
      allocate(labint(9*mxcent))

      kpatom = 0
      tofile = .false.
      trimat = .true.
      exp1vl = .false.
      runqm3 = .true.

      ! compute electric field if needed
      if (donppol) then

         nocomp = 3

         ! loop over perturbed densities
         do i = 1, nsim
            ioff = (i-1)*idimension
            isimoff = (i-1)*n2orbx+1

            ! loop over np centers
            do j = 1, tnpatm
               joff = ioff+(j-1)*3
               diporg = npcord(:, j)

               intao = 0.0d0
               call get1in(intao,'NEFIELD',nocomp,work,     &
                           lwork,labint,intrep,intadr,j,tofile,kpatom,    &
                           trimat,dummy,exp1vl,dummy,0)

               do ixyz = 1, 3

                  trmo = 0.0d0
                  utr = 0.0d0
                  if (.not. present(fvvec2)) then
                     utrx = 0.0d0
                     urxac = 0.0d0
                  end if

                  ! transform integrals
                  call uthu(intao((ixyz-1)*nnbasx + 1),trmo,cmo,work,nbast,norbt)
                  call dsptsi(norbt,trmo,utr)

                  if (present(fvvec2)) then
                     ! determine electric field component size
                     f1val = 0.0d0
                     f2val = 0.0d0
                     if (isymt.eq.isymv1) then
                        tlma = 0.0d0
                        call oith1(isymv1,zym1,utr,tlma,isymt)
                        call melone(tlma,1,udv,1.0d0,f1val,200,'qmnpqro')
                        fvvec1(joff+ixyz) = f1val
                     end if
                     if (isymt.eq.muld2h(isymv1,isymv2)) then
                        tlma = 0.0d0
                        tlmb = 0.0d0
                        call oith1(isymv1,zym1,utr,tlma,isymt)
                        call oith1(isymv2,zym2,tlma,tlmb,isymv2)
                        call melone(tlmb,1,udv,1.0d0,f2val,200,'qmnpqro')
                        fvvec2(joff+ixyz) = f2val
                     endif
                  else
                     call onexh1(bovecs(isimoff),utr,utrx)
                     if (nasht.gt.0) call getacq(utrx,urxac)
                     if (trplet) then
                         fvvec1(joff+ixyz) = slvtlm(udvtr,urxac,utrx,tac)
                     else
                         fvvec1(joff+ixyz) = slvqlm(udv,urxac,utrx,tac)
                     endif
                  end if
               end do
            end do
         end do
      end if

      if (donpcap) then

         istart = 0
         if (donppol) istart = 3*tnpatm

         nocomp = 1

         ! loop over perturbed first order densities
         do i = 1, nsim
            ioff = (i-1)*idimension
            isimoff = (i-1)*n2orbx+1

            ! loop over np centers
            do j = 1, tnpatm
               joff = ioff+istart+j
               diporg = npcord(:, j)

               intao = 0.0d0
               call get1in(intao,'NPETES ',nocomp,work,     &
                           lwork,labint,intrep,intadr,j,tofile,kpatom,    &
                           trimat,dummy,exp1vl,dummy,0)

               trmo = 0.0d0
               utr = 0.0d0
               if (.not. present(fvvec2)) then
                  utrx = 0.0d0
                  urxac = 0.0d0
               end if

               ! transform integrals
               call uthu(intao,trmo,cmo,work,nbast, norbt)
               call dsptsi(norbt,trmo,utr)

               if (present(fvvec2)) then
                  f1val = 0.0d0
                  f2val = 0.0d0
                  if (isymt.eq.isymv1) then
                     tlma = 0.0d0
                     call oith1(isymv1,zym1,utr,tlma,isymt)
                     call melone(tlma,1,udv,1.0d0,f1val,200,'qmnpqro')
                     fvvec1(joff) = f1val
                  end if
                  if (isymt.eq.muld2h(isymv1,isymv2)) then
                     tlma = 0.0d0
                     tlmb = 0.0d0
                     call oith1(isymv1,zym1,utr,tlma,isymt)
                     call oith1(isymv2,zym2,tlma,tlmb,isymv2)
                     call melone(tlmb,1,udv,1.0d0,f2val,200,'qmnpqro')
                     fvvec2(joff) = f2val
                  endif
               else
                  call onexh1(bovecs(isimoff),utr,utrx)
                  if (nasht.gt.0) call getacq(utrx,urxac)
                  if (trplet) then
                      fvvec1(joff) = slvtlm(udvtr,urxac,utrx,tac)
                  else
                      fvvec1(joff) = slvqlm(udv,urxac,utrx,tac)
                  endif
               end if
            end do

            ! set lagrangian for charge equilibration
            do iblk = 1, tnpblk
               fvvec1(ioff + idimension) = fvvec1(ioff + idimension) + npchrg(iblk)
            end do
            if (present(fvvec2)) then
               do iblk = 1, tnpblk
                  fvvec2(ioff + idimension) = fvvec2(ioff + idimension) + npchrg(iblk)
               end do
            end if
         end do
      end if

      ! print final fv vector
      do i = 1, nsim
         if (iprtlvl > 14 .and. .not. mqiter) then
            ioff = (i - 1)*idimension + 1
            write(lupri, '(/,2x,a,i0)') &
                '*** Computed FV vector start 1 ', i
            do j = 1, idimension
               write(lupri, '(i8, f18.8)') j, fvvec1(IOFF + j - 1)
            end do
            write(lupri, '(/,2x,a)') '*** Computed FV vector end ***'
            if (present(fvvec2)) then
               write(lupri, '(/,2x,a,i0)') &
                   '*** Computed FV vector start 2 ', i
               do j = 1, idimension
                  write(lupri, '(i8, f18.8)') j, fvvec2(IOFF + j - 1)
               end do
               write(lupri, '(/,2x,a)') '*** Computed FV vector end ***'
            end if
         end if
      end do

      if (allocated(intao))  deallocate(intao)
      if (allocated(trmo))   deallocate(trmo)
      if (allocated(utr))    deallocate(utr)
      if (allocated(tlma))   deallocate(tlma)
      if (allocated(tlmb))   deallocate(tlmb)
      if (allocated(utrx))   deallocate(utrx)
      if (allocated(urxac))  deallocate(urxac)
      if (allocated(intrep)) deallocate(intrep)
      if (allocated(intadr)) deallocate(intadr)
      if (allocated(labint)) deallocate(labint)

      runqm3 = .false.

      ! restore origin coordinates
      diporg = dipole_origin_save

   end subroutine


   ! computes contribution to xy vector from induced dipoles moments and charge
   subroutine get_xyvec(cmo,     &
                        idimension,    &
                        nsim,    &
                        fvec,    &
                        work,    &
                        lwork,   &
                        fmqvec1, &
                        fmqvec2, &
                        isymt,   &
                        isymv2,  &
                        zym2)

      ! molecular orbital coefficients
      real(8), intent(in)           :: cmo(*)
      ! size of electric field/potential vector
      integer, intent(in)           :: idimension
      ! number of perturbed density matrice
      integer, intent(in)           :: nsim
      integer, intent(in)           :: lwork
      real(8), intent(inout)        :: fvec(*)
      real(8), intent(inout)        :: work(*)
      ! induced dipole moments and charges
      real(8), intent(in)           :: fmqvec1(*)
      real(8), intent(in), optional :: fmqvec2(*)
      integer, intent(in), optional :: isymt
      integer, intent(in), optional :: isymv2
      real(8), intent(in), optional :: zym2(*)

#include "dummy.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "infdim.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "infrsp.h"

      logical :: tofile
      logical :: trimat
      logical :: exp1vl
      real(8) :: dipole_origin_save(3)
      integer :: i
      integer :: j
      integer :: ioff
      integer :: joff
      integer :: istart
      integer :: kpatom
      integer :: nocomp
      integer :: isimoff
      real(8) :: fact
      real(8) :: fact1
      real(8) :: fact2
      integer :: ixyz

      real(8), allocatable :: utr(:)
      real(8), allocatable :: trmo(:)
      real(8), allocatable :: intao(:)
      real(8), allocatable :: tlma(:)
      integer, allocatable :: intrep(:)
      integer, allocatable :: intadr(:)
      character(8), allocatable :: labint(:)

      ! save origin coordinates
      dipole_origin_save = diporg

      allocate(utr(n2orbx))
      allocate(trmo(nnorbx))
      allocate(intao(3*nnbasx))
      if (present(fmqvec2)) then
          allocate(tlma(n2orbx))
      end if
      allocate(intrep(9*mxcent))
      allocate(intadr(9*mxcent))
      allocate(labint(9*mxcent))

      kpatom = 0
      tofile = .false.
      trimat = .true.
      exp1vl = .false.
      runqm3 = .true.

      ! induced dipole moment in np region interaction with qm region
      if (donppol .and. novdamp) then
         nocomp = 3
         do i = 1, nsim
            ioff = (i - 1)*idimension
            isimoff = (i - 1)*n2orbx + 1
            do j = 1, tnpatm
               joff = ioff + (j - 1)*3
               diporg = npcord(:, j)
               intao = 0.0d0
               call get1in(intao, 'NEFIELD', nocomp, work,                   &
                           lwork, labint, intrep, intadr, j, tofile, kpatom, &
                           trimat, dummy, exp1vl, dummy, 0)

               ! loop over x, y, z
               do ixyz = 1, 3

                  ! transform integrals
                  utr = 0.0d0
                  trmo = 0.0d0
                  call uthu(intao((ixyz - 1)*nnbasx + 1), trmo, cmo, work, nbast, norbt)
                  call dsptsi(norbt, trmo, utr)

                  ! determine mm region contribution
                  fact1 = -fmqvec1(joff + ixyz)
                  if (present(fmqvec2)) then
                     tlma = 0.0d0
                     call oith1(isymv2, zym2, utr, tlma, isymt)
                     call daxpy(n2orbx, fact1, tlma, 1, fvec, 1)
                     fact2 = -0.5d0*fmqvec2(joff + ixyz)
                     call daxpy(n2orbx, fact2, tlma, 1, fvec, 1)
                  else
                     call daxpy(n2orbx, fact1, utr, 1, fvec(isimoff), 1)
                  end if
               end do
            end do
         end do
      end if

      ! induced dipole moment in np region interaction with qm region
      if (donpcap .and. novdamp) then
         istart = 0
         if (donppol) istart = 3*tnpatm
         nocomp = 1
         do i = 1, nsim
            ioff = (i - 1)*idimension
            isimoff = (i - 1)*n2orbx + 1
            do j = 1, tnpatm
               joff = ioff + istart + j
               diporg = npcord(:, j)
               intao = 0.0d0
               call get1in(intao, 'NPETES ', nocomp, work,                   &
                           lwork, labint, intrep, intadr, j, tofile, kpatom, &
                           trimat, dummy, exp1vl, dummy, 0)

               ! transform integrals
               utr = 0.0d0
               trmo = 0.0d0
               call uthu(intao, trmo, cmo, work, nbast, norbt)
               call dsptsi(norbt, trmo, utr)

               ! determine mm region contribution
               fact1 = fmqvec1(joff)
               if (present(fmqvec2)) then
                  tlma = 0.0d0
                  call oith1(isymv2, zym2, utr, tlma, isymt)
                  call daxpy(n2orbx, fact1, tlma, 1, fvec, 1)
                  fact2 = 0.5d0*fmqvec2(joff)
                  call daxpy(n2orbx, fact2, tlma, 1, fvec, 1)
               else
                  call daxpy(n2orbx, fact1, utr, 1, fvec(isimoff), 1)
               end if
            end do
         end do
      end if

      if (allocated(utr))    deallocate(utr)
      if (allocated(trmo))   deallocate(trmo)
      if (allocated(intao))  deallocate(intao)
      if (allocated(tlma))   deallocate(tlma)
      if (allocated(intrep)) deallocate(intrep)
      if (allocated(intadr)) deallocate(intadr)
      if (allocated(labint)) deallocate(labint)

      ! restore origin coordinates
      diporg = dipole_origin_save

   end subroutine

end module
