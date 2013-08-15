C
C/* Deck energy_pcm_drv */
C
      subroutine energy_pcm_drv(dcao, dvao, esolt, work, lfree)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <thrzer.h>
#include <iratdef.h>
#include <mxcent.h>
#include <codata.h>
#include <nuclei.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <dftcom.h>

      dimension dcao(*), dvao(*), work(*)
      character*7 potName, chgName, chgName1, chgName2
      integer nts

      call get_cavity_size(nts)
      kda   = 1
      kdb   = kda   + nnbasx
      kpot  = kdb   + nnbasx
      kcent = kpot  + nts
      kfree = kcent + 3 * nts
      lwork = lfree - kfree + 1
      if (lwork .lt. 0) call errwrk('energy_pcm_drv', kfree, lfree) 

      call get_tess_centers(work(kcent))
      
      if ((nasht.gt.0) .and. .not. dftadd) then
         call dcopy(nnbasx, dcao, 1, work(kdb), 1)
         call daxpy(nnbasx, 1.0d0, dvao, 1, work(kdb), 1)
         call pksym1(work(kda), work(kdb), nbas, nsym, -1)
      else
         call pksym1(work(kda), dcao, nbas, nsym, -1)
      end if
C
C 1) Get tessera data
C 2) Compute potentials
C 3) Compute charges
C

      potName = 'NucPot'//CHAR(0)
      chgName = 'NucChg'//CHAR(0)
      call nuc_pot_pcm(nts, work(kcent), work(kpot))
      call set_surface_function(nts, work(kpot), potName)
      call comp_chg_pcm(potName, chgName)

      potName = 'ElePot'//CHAR(0)
      chgName = 'EleChg'//CHAR(0)
      call ele_pot_pcm(nts, work(kcent), work(kpot), work(kda), 
     &                 work(kfree), lfree)
      call set_surface_function(nts, work(kpot), potName)
      call comp_chg_pcm(potName, chgName)
      call comp_pol_ene_pcm(esolt, 0)

      chgName  = 'TotChg'//CHAR(0)
      chgName1 = 'NucChg'//CHAR(0)
      chgName2 = 'EleChg'//CHAR(0)

      call append_surf_func(chgName)
      call clear_surf_func(chgName)
      factor =  1.0
      call add_surface_function(chgName, factor, chgName1)
      call add_surface_function(chgName, factor, chgName2)


      return
      end
C
C/* Deck nuclear_pot */
C
      subroutine nuc_pot_pcm(nts, centers, potential)

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <maxorb.h>
#include <infpar.h>

      dimension centers(3,*), potential(*)

      do i = 1, nts
         potential(i) = 0.0d0
         do j = 1, nctot
            dist = 0.0d0
            do k = 1, 3
               dist = dist + (centers(k,i) - cord(k, j))**2
            end do
            dist = sqrt(dist)
            potential(i) = potential(i) + charge(j) / dist
         end do
      end do
      return
      end
      
C/* Deck POTEXPVAL */
C
C Calculate exp values of potentials on tesserae
C Input: symmetry packed Density matrix in AO basis
C        cavity points
C Output: expectation values of electrostatic potential on tesserae
C Caomment: the density should not be passed. I cannot assume that 
C the density is in a given format here as this routine will be called by the module
C
      Subroutine ele_pot_pcm(nts, centers, potential, density, work, 
     $                       lwork)

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
#include <inforb.h>
#include <orgcom.h>
#include <maxorb.h>
c
      dimension density(*), work(*), potential(*), centers(3,*)
      call j1int_pcm(potential, nts, centers, 
     &     .true., density, 1, .false., 'NPETES ', 1, work, lwork)
      return
      end
C/* Deck Fock_PCMModule */
C
C Calculate exp values of potentials on tesserae
C Input: symmetry packed Density matrix in AO basis
C        cavity points
C Output: expectation values of electrostatic potential on tesserae
C
      subroutine oper_ao_pcm_drv(operator, charge, work, lwork)

#include <implicit.h>
c
      double precision work(*), operator(*)
      character charge(*)
       
      call get_cavity_size(nts)
      kcharge = 1
      kcenters = kcharge + nts
      kfree = kcenters + 3 * nts
      lfree = lwork - kfree

      if(lfree .le. 0) call quit('Not enough mem in fock_pcm_drv')

      call get_surface_function(nts, work(kcharge), charge)
      call get_tess_centers(work(kcenters))
      call j1int_pcm(work(kcharge), nts, work(kcenters), .false., 
     &     operator, 1, .false., 'NPETES ', 1, work(kfree), lfree)
      end
c
c/*deck j1int_pcm*/
c
      subroutine j1int_pcm(expval, nts, centers, 
     $     exp1vl, denmat, nosim, tofile, intlab,
     &     ksymp, work, lwork)
c
c
c     denmat : if (exp1vl) then input density matrix; else output matrix
c     
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <codata.h>
#include <infpri.h>
#include <nuclei.h>
#include <inforb.h>
#include <orgcom.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include "infpar.h"
#include "inftap.h"
#include "krcom.h"
c
      character*7 intlab
      character*8 labint(9*mxcent)
      dimension intrep(9*mxcent), intadr(9*mxcent), expval(nts,nosim)
      dimension denmat(*), work(lwork), centers(3,*)
      logical tofile, trimat, vcheck, exp1vl
c
clf do not insert other definitions below the following include
#include <ibtfun.h>
c

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
c
c     we use as a quick way of transfering tessera coordinates to hermit
c     the dipole origin. need to be restored.
c
      xdiporg = diporg(1)
      ydiporg = diporg(2)
      zdiporg = diporg(3)
c
c  2) calculation of apparent charges generated by the solute's nuclei.
c
      iprpcm=0
      if (exp1vl) then
         kden = 1
         klast = kden + nnbasx
         lwrk  = lwork - klast
         isymd = ksymp - 1
         if (ksymp .eq. 1) then
            call pksym1(denmat,work(klast),nbas,nsym,1)
            call dsym1(work(kden),dummy,work(klast),dummy,.false.,
     &                 nbast,iprpcm)
         else
            call dcopy(nnbasx,denmat,1,work(kden),1)
         end if
         if (nosim .gt. 1) call quit('nosim .gt. 1 and exp1vl not '//
     &        'permitted in j1int')
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
c
c        calculates nuclear potential energy integrals (in ao basis) for
c        the given tessera
c
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
         call get1in(work(kmat),intlab,ncomp,work(klast),lwrk,labint,
     &               intrep,intadr,l,tofile,kpatom,trimat,work(ktmp),
     &               exp1vl,work(kden),iprpcm)
         if (iprtyp .eq. 13) then
            jmat = kmat
            do iosim = 1, nosim
               call daxpy(matdim,expval(its,iosim),work(jmat),1,
     &                    denmat(matdim*(iosim - 1) + 1),1)
               jmat = jmat + matdim
            end do
         else if (exp1vl) then
            do i = 1, ncomp
               expval(its+(i-1)*ntsirr,1) = -work(ktmp+i-1)
            end do
         else if (.not. tofile) then
            do iosim = 1, nosim
               iadr = kmat + (ksymp - 1)*matdim
               call daxpy(matdim,-expval(its,iosim),work(iadr),1,
     &                    denmat(matdim*(iosim - 1) + 1),1)
            end do
         end if
      enddo 
      diporg(1) = xdiporg
      diporg(2) = ydiporg
      diporg(3) = zdiporg
      nbast  = nbastold
      nnbasx = nnbasxold
      return
      end
C  /* Deck getacord */
      subroutine getacord(coora)      
C*****************************************************************************
C
C    getacord : Make list atomic coordinates
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C               Copied here from DIRAC by Roberto Di Remigio, February 2012
C
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
      dimension coora(3,*)
C
C     Make the full matrix of cartesian coordinates from CORD(NUCIND)
C
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
                  coora(1,jatom) = 
     &                 pt(iand(isymax(1,1),isymop))*cord(1,icent)
                  coora(2,jatom) = 
     &                 pt(iand(isymax(2,1),isymop))*cord(2,icent)
                  coora(3,jatom) = 
     &                 pt(iand(isymax(3,1),isymop))*cord(3,icent)
              end if
            enddo
        end if
      enddo 
      return
      end


C
C Different Cavity Modes. Mode = Implicit
C nSpheres is now NCTOT defined in nuclei.h
C
      subroutine collect_nctot(nuclei)

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      
      nuclei = nucdep

      return
      end
      
      
      subroutine collect_atoms(charges, centers, unit_flag)

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <maxorb.h>
#include <infpar.h>
      
      dimension charges(*)
      dimension centers(3,*)
      integer unit_flag
      
      call getacord(centers)

      i=0
      do j=1,nucind
         do k=1,nucdeg(j)
            i=i+1
            charges(i) = charge(j)
         enddo
      enddo
      
C If units is AU unit_flag = 0. If Angstrom unit_flag = 1.
      unit_flag = 0
      
      return
      end
      
      subroutine pot_int_tess(potint, tessera, trimat, work, lwork)

#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <priunit.h>
#include <orgcom.h>
#include <inforb.h>

      double precision potint(*), tessera(3), work(*)
      logical trimat, tofile, exp1vl
      character*7 intlab
      character*8 labint(9*mxcent)
      integer intrep(9*mxcent), intadr(9*mxcent)

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

      call get1in(potint,intlab,ncomp,work(kfree),lfree,labint,
     &     intrep,intadr,il,tofile,kpatom,trimat,dummy,
     &     exp1vl,dummy,iprint)

      write(lupri, *) "Potentials for tessera rsp",1,(diporg(j), j=1,3),
     $     mxcent
      call outpak(potint, norbt, 1, lupri)

      diporg(1) = xdiporg
      diporg(2) = ydiporg
      diporg(3) = zdiporg

      return
      end

