!
! DALTON-side interface routines for the Polarizable Continuum Model
! Luca Frediani, Roberto Di Remigio 2011-2013
!
! We divide the interface routines into:
!   1. cavity-related routines;
!   2. solver-related routines.
! 
! The driver routine and the solver-related routines will be provided
! through a FORTRAN90 module. The cavity-related routines are called 
! by C++. Thus we put them outside the module, to avoid confusion due
! to name mangling.
!

!
!                   CAVITY-RELATED ROUTINES
!
! We shall provide the following data to the module:
!   1. the number of nuclei in the molecule. NUCLEI variable;
!   2. the atomic number of each nucleus. CHARGES vector;
!   3. the coordinates of each nucleus. CENTERS matrix.
! PCMSolver will provide to cavity formation depending
! on the selected mode:
!   1. Mode = Implicit means that the spheres are
!      centered according to the CENTERS matrix with
!      radii taken from the internal library;
!   2. Mode = Atoms means that we specify sphere centers
!      and radii for certain atoms in the molecule. All the
!      other spheres are created as in the Implicit mode.
!

      subroutine collect_nctot(nr_nuclei)
              implicit none
              integer(4), intent(out) :: nr_nuclei

#include "mxcent.h"
#include "nuclei.h"
      
              nr_nuclei = nucdep

      end subroutine
      
      
      subroutine collect_atoms(charges, centers)
      implicit none

#include "mxcent.h"
#include "nuclei.h"
      
      real(8), intent(out) :: charges(*)
      real(8), intent(out) :: centers(3,*)
      
      integer :: i, j, k 

! Get coordinates
      call getacord(centers)
! Get charges      
      i = 0
      do j = 1, nucind
         do k = 1, nucdeg(j)
            i = i + 1
            charges(i) = charge(j)
         enddo
      enddo
      
      end subroutine

      subroutine getacord(coora)      
!*****************************************************************************
!
!    getacord : Make list atomic coordinates
!
!               Written oct.2001 by Jesper Kielberg Pedersen
!               Copied here from DIRAC by Roberto Di Remigio, February 2012
!
!*****************************************************************************
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
      real(8), intent(out) :: coora(3,*)
!
!     Make the full matrix of cartesian coordinates from CORD(NUCIND)
!
      jatom = 0
      do icent = 1, nucind
         mulcnt = istbnu(icent)
         if (mult(mulcnt) .eq. 1) then
            jatom = jatom + 1
            coora(1,jatom) = cord(1,icent)
            coora(2,jatom) = cord(2,icent)
            coora(3,jatom) = cord(3,icent)
        else
            do isymop = 0, maxopr
            if (iand(isymop,mulcnt) .eq. 0) then
                  jatom = jatom + 1
                  coora(1,jatom) = pt(iand(isymax(1,1),isymop))*cord(1,icent)
                  coora(2,jatom) = pt(iand(isymax(2,1),isymop))*cord(2,icent)
                  coora(3,jatom) = pt(iand(isymax(3,1),isymop))*cord(3,icent)
              end if
            enddo
        end if
      enddo 

     end subroutine

module pcm_interface
   
   use pcm_write
   use pcmmod_cfg

   implicit none 
  
   public energy_pcm_drv
   public oper_ao_pcm_drv
!   public pcm_energy_driver   
!   public pcm_fock_driver
!   public pcm_lintra_driver
!   public get_pcm_energy

   private

   real(8), allocatable :: charges(:)
   real(8), allocatable :: potentials(:)
   real(8), allocatable :: tess_cent(:, :)
   real(8)              :: pcm_energy
! nr_points is 32 bit to avoid mess in passing it from C++ to FORTRAN
   integer(4)           :: nr_points = -1

   contains 

      subroutine energy_pcm_drv(dcao, dvao, pol_ene, work, lfree)

#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "inftap.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "dftcom.h"

        real(8)      :: dcao(*)
        real(8)      :: dvao(*)                                                       
        real(8)      :: pol_ene
        real(8)      :: work(*)
        integer      :: lfree
! Local variables
        character(7) :: potName, chgName, chgName1, chgName2
        integer(4)   :: nts
        integer      :: kda, kdb, kpot, kcent, kfree, lwork
        real(8)      :: factor
                                                                                      
        call get_cavity_size(nts)
        kda   = 1
        kdb   = kda   + nnbasx
        kpot  = kdb   + nnbasx
        kcent = kpot  + nts
        kfree = kcent + 3 * nts
        lwork = lfree - kfree + 1
        if (lwork .lt. 0) then 
                call errwrk('energy_pcm_drv', kfree, lfree) 
        end if
                                                                                      
! 1) Get tessera data
        call get_tess_centers(work(kcent))
        
        if ((nasht.gt.0) .and. .not. dftadd) then
           call dcopy(nnbasx, dcao, 1, work(kdb), 1)
           call daxpy(nnbasx, 1.0d0, dvao, 1, work(kdb), 1)
           call pksym1(work(kda), work(kdb), nbas, nsym, -1)
        else
           call pksym1(work(kda), dcao, nbas, nsym, -1)
        end if
                                                                                      
! 2) Compute potentials
! 3) Compute charges
        potName = 'NucPot'//CHAR(0)
        chgName = 'NucChg'//CHAR(0)
        call nuc_pot_pcm(nts, work(kcent), work(kpot))
        call set_surface_function(nts, work(kpot), potName)
        call comp_chg_pcm(potName, chgName)
                                                                                      
        potName = 'ElePot'//CHAR(0)
        chgName = 'EleChg'//CHAR(0)
        call ele_pot_pcm(nts, work(kcent), work(kpot), work(kda), work(kfree), lfree)
        call set_surface_function(nts, work(kpot), potName)
        call comp_chg_pcm(potName, chgName)
        call comp_pol_ene_pcm(pol_ene, 0)
                                                                                      
        chgName  = 'TotChg'//CHAR(0)
        chgName1 = 'NucChg'//CHAR(0)
        chgName2 = 'EleChg'//CHAR(0)
                                                                                      
        call append_surf_func(chgName)
        call clear_surf_func(chgName)
        factor =  1.0
        call add_surface_function(chgName, factor, chgName1)
        call add_surface_function(chgName, factor, chgName2)

      end subroutine
      
      subroutine nuc_pot_pcm(nts, centers, potential)

#include "mxcent.h"
#include "nuclei.h"
        
        integer(4), intent(in) :: nts
        real(8), intent(in)    :: centers(3, nts)
        real(8), intent(out)   :: potential(nts)
        real(8)                :: dist
        integer                :: i, j, k
        
        do i = 1, nts                                         
           potential(i) = 0.0d0
           do j = 1, nucdep
              dist = 0.0d0
              do k = 1, 3
                 dist = dist + (centers(k,i) - cord(k, j))**2
              end do
              dist = sqrt(dist)
              potential(i) = potential(i) + charge(j) / dist
           end do
        end do
      
      end subroutine
      
!
! Calculate exp values of potentials on tesserae
! Input: symmetry packed Density matrix in AO basis
!        cavity points
! Output: expectation values of electrostatic potential on tesserae
! Caomment: the density should not be passed. I cannot assume that 
! the density is in a given format here as this routine will be called by the module
!
      subroutine ele_pot_pcm(nts, centers, potential, density, work, lwork)

        integer, intent(in)  :: nts
        real(8), intent(in)  :: centers(3, *)
        real(8), intent(out) :: potential(nts)
        real(8) :: density(*)
        real(8) :: work(*)
        integer :: lwork
      
        call j1int_pcm(potential, nts, centers, .true., density, 1, .false., 'NPETES ', 1, work, lwork)
      
     end subroutine
!
! Calculate exp values of potentials on tesserae
! Input: symmetry packed Density matrix in AO basis
!        cavity points
! Output: expectation values of electrostatic potential on tesserae
!
      subroutine oper_ao_pcm_drv(oper, charge, work, lwork)

        real(8), intent(out) :: oper(*)
        real(8)              :: work(*)
        character            :: charge(*)
        integer              :: lwork

        integer(4)           :: nts
        integer              :: kcharge, kcenters, kfree, lfree

       
        call get_cavity_size(nts)                                                      
        kcharge = 1
        kcenters = kcharge + nts
        kfree = kcenters + 3 * nts
        lfree = lwork - kfree
                                                                                  
        if(lfree .le. 0) then 
                call quit('Not enough mem in oper_ao_pcm_drv')
        end if
        
        call get_surface_function(nts, work(kcharge), charge)
        call get_tess_centers(work(kcenters))
        call j1int_pcm(work(kcharge), nts, work(kcenters), .false.,            & 
                       oper, 1, .false., 'NPETES ', 1, work(kfree), lfree)

      end subroutine
      
      subroutine j1int_pcm(expval, nts, centers, exp1vl, denmat, nosim, tofile, intlab, ksymp, work, lwork)
!
!     denmat : if (exp1vl) then input density matrix; else output matrix
!     
#include "dummy.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "orgcom.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "infpar.h"
#include "inftap.h"
     
        character(7) :: intlab
        character(8) :: labint(9*mxcent)
        integer(4)   :: nts
        integer      :: isum
        integer      :: nosim, ksymp, lwork, i, j, iaddr, iosim, iprpcm
        integer      :: iadr, iprtyp, isymd, its, jmat, kden, klast
        integer      :: kmat, kpatom, ktmp, l, lwrk, matdim, nbastold
        integer      :: ncomp, nnbasxold, ntesp, ntsirr
        real(8)      :: intrep(9*mxcent), intadr(9*mxcent), expval(nts, nosim)
        real(8)      :: denmat(*), work(lwork), centers(3, nts)
        real(8)      :: xdiporg, ydiporg, zdiporg
        logical      :: tofile, trimat, vcheck, exp1vl
!
!lf do not insert other definitions below the following include
!#include "ibtfun.h"
!rdr move contents of ibtfun.h in functions in this module
      nbastold  = nbast
      nnbasxold = nnbasx
      nbast  = isum(maxrep+1,naos,1)
      nnbasx = nbast * (nbast + 1)/2
      n2basx = nbast*nbast
      if (intlab .eq. 'PCMBSOL') then
         matdim = n2basx
         iprtyp = 13
         ncomp  = 3
         trimat = .false.
      else
         matdim = nnbasx
         iprtyp = 11
         ncomp  = 1
         trimat = .true.
      end if
!
!     we use as a quick way of transfering tessera coordinates to hermit
!     the dipole origin. need to be restored.
!
      xdiporg = diporg(1)
      ydiporg = diporg(2)
      zdiporg = diporg(3)
!
!  2) calculation of apparent charges generated by the solute's nuclei.
!
      iprpcm=0
      if (exp1vl) then
         kden = 1
         klast = kden + nnbasx
         lwrk  = lwork - klast
         isymd = ksymp - 1
         if (ksymp .eq. 1) then
            call pksym1(denmat,work(klast),nbas,nsym,1)
            call dsym1(work(kden),dummy,work(klast),dummy,.false.,nbast,iprpcm)
         else
            call dcopy(nnbasx,denmat,1,work(kden),1)
         end if
         if (nosim .gt. 1) call quit('nosim .gt. 1 and exp1vl not permitted in j1int')
      else
         kden = 1
         if (nodtot .ge. 1) then
            klast = kden + matdim*nosim
            call dcopy(matdim*nosim,denmat,1,work(kden),1)
         else
            klast = kden
         end if
         lwrk = lwork - klast
      end if
      do  its = 1, nts
         diporg(1) = centers(1,its)
         diporg(2) = centers(2,its)
         diporg(3) = centers(3,its)

         ntesp = 1
         kpatom = 0
!
!        calculates nuclear potential energy integrals (in ao basis) for
!        the given tessera
!
         l=1
         ktmp = klast
         if (.not. tofile .and. .not. exp1vl) then
            kmat = ktmp + 8
            if (iprtyp .eq. 11) then
               klast = kmat + (maxrep + 1)*matdim
            else
               klast = kmat + (maxrep + 1)*matdim*ncomp
            end if
            ncomp = nsym
         else
            kmat  = ktmp + 8
            klast = kmat
            ncomp = 0
         end if
         call get1in(work(kmat),intlab,ncomp,work(klast),lwrk,labint, &
                     intrep,intadr,l,tofile,kpatom,trimat,work(ktmp), &
                     exp1vl,work(kden),iprpcm)
         if (iprtyp .eq. 13) then
            jmat = kmat
            do iosim = 1, nosim
               call daxpy(matdim,expval(its,iosim),work(jmat),1,denmat(matdim*(iosim - 1) + 1),1)
               jmat = jmat + matdim
            end do
         else if (exp1vl) then
            do i = 1, ncomp
               expval(its+(i-1)*ntsirr,1) = -work(ktmp+i-1)
            end do
         else if (.not. tofile) then
            do iosim = 1, nosim
               iadr = kmat + (ksymp - 1)*matdim
               call daxpy(matdim,-expval(its,iosim),work(iadr),1,denmat(matdim*(iosim - 1) + 1),1)
            end do
         end if
      enddo 
      diporg(1) = xdiporg
      diporg(2) = ydiporg
      diporg(3) = zdiporg
      nbast  = nbastold
      nnbasx = nnbasxold
      
      end
      
      subroutine pot_int_tess(potint, tessera, trimat, work, lwork)

#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "orgcom.h"
#include "inforb.h"
        
        real(8), intent(out) :: potint(*)
        real(8), intent(in)  :: tessera(3)
        real(8)              :: work(*)
        logical              :: trimat, tofile, exp1vl
        character(7)         :: intlab
        character(8)         :: labint(9*mxcent)
        integer              :: intrep(9*mxcent), intadr(9*mxcent), lwork
        integer              :: il, iprint, j, kfree, kpatom, lfree, ncomp
        real(8)              :: xdiporg, ydiporg, zdiporg

        xdiporg = diporg(1)                                                      
        ydiporg = diporg(2)
        zdiporg = diporg(3)
        diporg(1) = tessera(1)
        diporg(2) = tessera(2)
        diporg(3) = tessera(3)
                                                                                 
        intlab = 'NPETES '
        ncomp = nsym
        il = 1
        tofile = .false.
        kpatom = 0
        exp1vl = .false.
        iprint = 0
                                                                                 
        kfree = 1
        lfree = lwork - kfree + 1
        if (lfree .lt. 0) call errwrk('pot_int_tess', kfree, lwork) 
                                                                                 
        call get1in(potint,intlab,ncomp,work(kfree),lfree,labint, &
                    intrep,intadr,il,tofile,kpatom,trimat,dummy,  &
                    exp1vl,dummy,iprint)
                                                                                 
        write(lupri, *) "Potentials for tessera rsp",1,(diporg(j), j=1,3),mxcent
        call outpak(potint, norbt, 1, lupri)
                                                                                 
        diporg(1) = xdiporg
        diporg(2) = ydiporg
        diporg(3) = zdiporg

      end subroutine

      function ibtand(i, j)
                
        integer :: i, j
        integer ibtand
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtand = and(i, j)
#else
        ibtand = iand(i, j)
#endif

      end function
      
      function ibtor(i, j)
                
        integer :: i, j
        integer ibtor
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtor = or(i, j)
#else
        ibtor = ior(i, j)
#endif

      end function
      
      function ibtshl(i, j)
                
        integer :: i, j
        integer ibtshl
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtshl = shiftl(i, j)
#else
        ibtshl = ishft(i, j)
#endif

      end function

      function ibtshr(i, j)
                
        integer :: i, j
        integer ibtshr
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtshr = shiftr(i, j)
#else
        ibtshr = ishft(i, j)
#endif

      end function
      
      function ibtxor(i, j)
                
        integer :: i, j
        integer ibtxor
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
        ibtxor = xor(i, j)
#else
        ibtxor = ieor(i, j)
#endif

      end function

end module
