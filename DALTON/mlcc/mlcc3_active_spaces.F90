!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
module mlcc3_active_spaces
!
!
!  mlcc3 initialiasion of active space variables
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
   use mlcc_typedef
   use mlcc_block_import
   use mlcc_work
   use mlcc3_data
!
!   
   implicit none
!
!  Active space variables
!
   integer  :: n_spaces, cc3_spaces, include_spaces
   integer  :: n_occ_act, n_vir_act, n_occ_gen, n_vir_gen
   integer  :: n_v_a_3
   integer  :: n_occ_inact, n_vir_inact, n_gen_inact
   integer  :: n_vir_aag, n_occ_aag !active*active*general
   integer  :: n_occ_integral, n_vir_integral !number of integrals
   integer, dimension(:), pointer   :: n_occ_space, n_vir_space
!
!  CVS variables
!
   logical  :: cvs_log, rmc_log, ion_log
   integer  :: n_cvs_rmc, n_ion
   integer, dimension(:), pointer   :: cvs_rmc_index, ion_index
!
contains
!
   subroutine center_h_import
!
      implicit none
!      
      integer, parameter   :: mxbsce   = 2000, mxacat = 50, mxacbs = 5000
      integer, parameter   :: mxexbs   = 500
      integer, parameter   :: mxspa    = 5, maxspa = mxspa+1
      integer, parameter   :: mxcent   = 500
!      
      integer :: allocate_ok, i
!
      double precision thacoc, thacvi
!
!      
      integer nbcent, ibcent, nbscen, kbscen, &
     &        nacatm, lactat, lacbas, nactbs, nabsto, &
     &        nactoc, nactvi, ninaoc, ninavi, &
     &        ioract, noract, nactfr, nextbs, iextbs, &
     &        nocvec, nvivec, nacinp, lacinp, iordec, &
     &        nseld,  iseld,  nfract, &
     &        nabsoc, nabsvi, nacbsv,lacbsv, &
     &        nspace, natoac, labspa, nspaoc, nspavi, &
     &        nlfthf, nlftvi, iof2hf, iof2vi, iepshf, iepsvi, &
     &        nabso2, lacba2, nabsv2, lacbv2, iexpvi, iexpoc
!     
      logical actsel, atomic, actfre, difadd, nboexp, seldir, &
     &        opnshl, fuldec, dospread, minspr, &
     &        dialst, extern, newact, spdils, loconl, addorb, &
     &        addexp
!     
      common /center/ thacoc, thacvi, &
     &                nbcent(mxcent,8), ibcent(mxcent,8), &
     &                nbscen(mxcent), kbscen(mxbsce,mxcent), &
     &                nacatm, lactat(mxacat), lacbas(mxacbs), &
     &                nactbs(8), nabsto, nactfr, &
     &                ioract(8), noract(8), &
     &                nactoc(8), nactvi(8), ninaoc(8), ninavi(8), &
     &                nextbs(8), iextbs(mxexbs,8), &
     &                nocvec(mxcent,8), nvivec(mxcent,8), &
     &                nacinp, lacinp(mxacat), iordec(mxcent), &
     &                nseld(mxacat), iseld(mxbsce,mxacat), nfract(8), &
     &                nspace, natoac(mxspa), labspa(mxacat,mxspa), &
     &                nspaoc(8,maxspa), nspavi(8,maxspa), &
     &                nlfthf(8), nlftvi(8), iof2hf(8,maxspa), &
     &                iof2vi(8,maxspa), iepshf(8,maxspa), &
     &                iepsvi(8,maxspa), &
     &                nabsoc, nabsvi, nacbsv(8),lacbsv(mxacbs), &
     &                nabso2(mxspa), lacba2(mxacbs,mxspa), &
     &                nabsv2(mxspa), lacbv2(mxacbs,mxspa), &
     &                actsel, atomic, actfre, difadd, nboexp, seldir, &
     &                opnshl, fuldec, dospread, minspr, &
     &                dialst, extern, newact, spdils, loconl, &
     &                addorb, addexp, iexpvi, iexpoc
!
!
!     CVS, RMC, and ionisation not implemented in this version
!     Set to to false
      cvs_log = .false.
      ion_log = .false.
      rmc_log = .false.
!
      if (mlcc3_active .and. .not. mlcc3_nrg_spa) then
!         
         n_spaces = nspace + 1
!
         write(lupri,*)
         write(lupri,*) 'nspace', nspace
         write(lupri,*) 'n_spaces', n_spaces
         write(lupri,*)
!
         if(n_general .gt. n_spaces) then
!            
            call quit('n_general greater than n_spaces in mlcc3 center import')
!
         end if
!
!      
         allocate(n_occ_space(n_spaces), stat = allocate_ok)
         if(allocate_ok /= 0) call quit('Something wrong with allocation in center_h_import')
!      
         allocate(n_vir_space(n_spaces), stat = allocate_ok)
         if(allocate_ok /= 0) call quit('Something wrong with allocation in center_h_import')
!      
         do i = 1,n_spaces
            n_occ_space(i) = nspaoc(1,i)
            n_vir_space(i) = nspavi(1,i)
         enddo
!      
         n_occ_act = 0
         n_vir_act = 0
         n_occ_gen = 0
         n_vir_gen = 0
!      
         cc3_spaces = n_active       ! Spaces treated with CC3 or higher
         include_spaces = n_general   ! Spaces to include for general orbital indexes
!      
         do i = 1,cc3_spaces
            n_occ_act = n_occ_act + n_occ_space(i)
            n_vir_act = n_vir_act + n_vir_space(i)
         end do
         do i = 1,include_spaces
            n_occ_gen = n_occ_gen + n_occ_space(i)
            n_vir_gen = n_vir_gen + n_vir_space(i)
         end do
!
         n_occ_inact = n_occ - n_occ_act
         n_vir_inact = n_vir - n_vir_act
!
         n_gen_inact = n_occ - n_occ_gen
!
      else if(mlcc3_nrg_spa) then
!
         n_occ_act = n_occ_inp
         n_vir_act = n_vir_inp
!
         if(mlcc3_nrg_gen) then
            n_occ_gen = n_gen_o_inp
            n_vir_gen = n_gen_v_inp
         else
            n_occ_gen = n_occ
            n_vir_gen = n_vir
         end if
!
         n_occ_inact = n_occ - n_occ_act
         n_vir_inact = n_vir - n_vir_act
!
         n_gen_inact = n_occ - n_occ_gen
!
      else
!
         n_occ_act = n_occ
         n_occ_gen = n_occ
         n_vir_act = n_vir
         n_vir_gen = n_vir
!
         n_occ_inact = 0
         n_vir_inact = 0
!
         n_gen_inact = 0
!
      end if
!
!     Sanity check
!
      if(n_occ_inact .lt. 0) then
         call quit('n_occ_inact less than 0')
      end if
!
      if(n_vir_inact .lt. 0) then
         call quit('n_vir_inact less than 0')
      end if
!
      if(n_gen_inact .lt. 0) then
         call quit('n_gen_inact less than 0')
      end if
!
      n_occ_aag = n_occ_act*n_occ_act*n_occ_gen
      n_vir_aag = n_vir_act*n_vir_act*n_vir_gen
!
      n_v_a_3     = n_vir_act**3
!
      n_occ_integral = n_occ_aag*n_vir_act
      n_vir_integral = n_vir_aag*n_occ_act
!
      if(print_mlcc3 .ge. 3) then
!
         write(lupri,*)
         write(lupri,*) 'Output from mlcc3_active spaces'
         write(lupri,*) 'n_basis:     ', n_basis
         write(lupri,*) 'n_orbitals:  ', n_orbitals
         write(lupri,*) 'n_occ:       ', n_occ
         write(lupri,*) 'n_occ_act:   ', n_occ_act
         write(lupri,*) 'n_occ_inact: ', n_occ_inact
         write(lupri,*) 'n_gen_inact: ', n_gen_inact
         write(lupri,*) 'n_occ_gen:   ', n_occ_gen
         write(lupri,*) 'n_vir:       ', n_vir
         write(lupri,*) 'n_vir_act:   ', n_vir_act
         write(lupri,*) 'n_vir_gen:   ', n_vir_gen
         write(lupri,*)
!
      end if
!
   end subroutine center_h_import
!   
!
   subroutine cvsexci_import
!
      implicit none
!
      integer  :: allocate_ok, i
!
      integer, parameter   :: maxcore=5, maxion=1
!
      integer  nrhfcore, irhfcore, nvirion, ivirion
      logical  lcvsexci, lionizexci, lbothexci, lrmcore, lcvsptexci
!
      common /cvsepexci/ nrhfcore(8),irhfcore(maxcore,8), &
     &               lcvsexci,lrmcore, &
     &               nvirion(8), ivirion(maxion,8), &
     &               lionizexci, lbothexci, lcvsptexci
!
!     rmc or cvs
!
      if(lcvsexci .or. lrmcore) then
!
         cvs_log = lcvsexci
         rmc_log = lrmcore
!
         n_cvs_rmc = nrhfcore(1)

         allocate(cvs_rmc_index(n_cvs_rmc), stat = allocate_ok)
         if(allocate_ok /= 0) call quit('Something wrong with allocation in center_h_import')
!      
         do i = 1,n_cvs_rmc
            cvs_rmc_index(i) = irhfcore(i,1)
         end do
!      
      else
!
         cvs_log = .false.
         rmc_log = .false.
!
         n_cvs_rmc = 0
         cvs_rmc_index => null()
!
      end if
!
!     Ionisation
!
      if(lionizexci) then
!
         ion_log = lionizexci
!
         n_ion = nvirion(1)

         allocate(ion_index(n_ion), stat = allocate_ok)
         if(allocate_ok /= 0) call quit('Something wrong with allocation in center_h_import')
!      
         do i = 1,n_ion
            ion_index = ivirion(i,1)
         end do
!      
      else
!
         ion_log = .false.
!
         n_ion = 0
         ion_index => null()
!
      end if
!
!     Offset indices in case of active space
!
      if(mlcc3_active) then
         if(cvs_log .or. rmc_log) then
            do i =1,n_cvs_rmc
               cvs_rmc_index(i) = cvs_rmc_index(i) - n_occ_inact
            end do
         end if
         if(ion_log) then
            do i =1,n_cvs_rmc
               ion_index(i) = ion_index(i) - n_vir_inact
            end do
         end if
      end if
!
!     Sanity check
!
      if(cvs_log .or. rmc_log) then
         do i =1,n_cvs_rmc
            if(cvs_rmc_index(i) .lt. 1) then
!
               write(lupri,*) 'i', i
               write(lupri,*) 'index', cvs_rmc_index(i)
               call quit('CVS/RMC index less than 1')
!
            else if(cvs_rmc_index(i) .gt. n_occ_act) then
!
               write(lupri,*) 'i', i
               write(lupri,*) 'index', cvs_rmc_index(i)
               write(lupri,*) 'n_occ_act', n_occ_act
               call quit('CVS/RMC index greater than n_occ_act')
!
            end if
         end do
      end if
!
      if(ion_log) then
         do i =1,n_ion
            if(ion_index(i) .lt. 1) then
!
               write(lupri,*) 'i', i
               write(lupri,*) 'index', cvs_rmc_index(i)
               call quit('Ion index less than 1')
!
            else if(ion_index(i) .gt. n_vir_act) then
!
               write(lupri,*) 'i', i
               write(lupri,*) 'index', ion_index(i)
               write(lupri,*) 'n_vir_act', n_vir_act
               call quit('Ion index greater than n_vir_act')
!
            end if
         end do
      end if
!
!
   end subroutine cvsexci_import
!   
end module mlcc3_active_spaces
