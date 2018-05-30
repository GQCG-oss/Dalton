C
C  /* Deck cc_molden */
      subroutine cc_molden_nto(work,lwork)
C
C     R. Faber et al., 2018
C
C     based on CC_MOLDEN by
C     Rolf H. Myhre and H. Koch
C
C     Purpose: Calculate natural transition orbitals and write them in
C     molden format to file
C
      implicit none

      character(len=*), parameter :: myname = 'CC_MOLDEN_NTO'
!
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccfop.h"
#include "ccsections.h"
#include "ccexgr.h"
#include "ccexci.h"
#include "ccexcinf.h"
#include "inftap.h"

!
      integer :: kcmot, knto, kucmo, koccu, krbtr, korbv
      integer :: kvec1, kumat, kvmat, ksval, kpos
      integer :: isyma, isymi, isym
      integer :: excisym, iexci
      integer :: mindim
      integer :: nto_count(8), ncount, nfound
      integer :: i_info, idummy, ilist, ilistoff, iopt
      integer ludena
      integer :: icmo(8)
      integer :: kend1, kend2, kend3
      integer :: lwrk1, lwrk2, lwrk3
      integer, intent(in) :: lwork
!      
      character(len=20):: filename
      character(len=10):: model
      character(len=3), parameter :: list_type = 'RE '
!      
#if defined (SYS_CRAY)
      real, intent(out):: work(lwork)
      real :: one, half, ddot
#else
      double precision, intent(out):: work(lwork)
      double precision, parameter :: one = 1.d0, half = 0.5d0,
     &                               zero = 0d0,
     &                               thresh_nto = 1.0D-2
      double precision :: ddot, DNRM2
      double precision trace, dummy
#endif
!
      model ='CCSD      '
!
      nto_count = 0
!     
      ncount = 0
      do isym = 1, nsym
         icmo(isym) = ncount
         ncount = ncount + nbas(isym)*norb(isym)
      end do

!     No symmetry allowed
!     -------------------
!      if (nsym .ne. 1) then
!         stop 'Symmetry not implemented'
!      end if
!
!     Dynamic allocation
!     ------------------
!
      kcmot = 1
      knto = kcmot + nlamds
      kucmo = knto + nlamds
      krbtr = kucmo + nbast*norbt
      korbv = krbtr + nbast*nbast
      koccu = korbv + nbast
      kend1 = koccu + norbt
      lwrk1 = lwork - kend1 + 1
      write(lupri,'(6I5)') kcmot, kucmo, krbtr, korbv, koccu, kend1
!
      if (lwrk1 .lt. 0) then
         stop 'Insufficient work space in '//myname
      endif

!     Zero the nto memory
      call dzero(work(knto), nlamds)
!     
!     Close molden_MOS, we want new files
!     -----------------------------------
      call gpclose(LUMOLDEN_MOS,'KEEP')
!      
!
!     Read in MO coefficients
!     -------------------------
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUERR)
      READ (LUSIFC)
C
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMOT+I-1),I=1,nlamds)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C      
      CALL CMO_REORDER(WORK(KCMOT), WORK(KEND1), LWRK1)
C
C
C
!     Get overlap matrix and square
!     -----------------------------
!
!      call rdonel('OVERLAP ',.true.,work(kover),nbast)
!
!      do ax = 1,nbast
!         do bx = 1,ax
!
!            abu = max(ax,bx)*(max(ax,bx)-3)/2+ax+bx
!            ab  = nbast*(ax-1) + bx
!            ba  = nbast*(bx-1) + ax
!
!            work(kuver+ab-1) = work(kover+abu-1)
!            work(kuver+ba-1) = work(kover+abu-1)
!
!         enddo
!      enddo

!
!     Natural transition orbitals
!     ---------------------
!     
      do excisym = 1, nsym
         
         ! symmetry specific entries
         kvec1 = kend1
         kend2 = kvec1 + nt1am(excisym)
         lwrk2 = lwork - kend2

         ilistoff = ISYOFE(excisym)

         do iexci = 1, nccexci(excisym,1)
            iopt = 1 ! Read only singles part of vector
            ilist = ilistoff + iexci
            ! Read in actual vector
            call cc_rdrsp(list_type, ilist, excisym, iopt,
     &                     model, work(kvec1), dummy)
            
            call dzero(work(koccu), norbt)

            write(LUPRI,'(F16.6)') DNRM2(NT1AM(EXCISYM),work(kvec1),1)
            write(LUPRI,'(A,I1,I3)') 'NTO singular values for exci ',
     &            excisym, iexci
            ! Loop over symmetry blocks of the vector
            do isymi = 1, nsym
               isyma = muld2h(isymi,excisym)
               
               mindim = min(nvir(isyma),nrhf(isymi))
               ! Skip if there is no work to do
               if (mindim.lt.1) cycle

               ! Allocate memory
               kumat = kend2
               kvmat = kumat + nvir(isyma)**2
               ksval = kvmat + nrhf(isymi)**2
               kend3 = ksval + mindim 
               lwrk3 = lwork - kend3
               
               kpos = kvec1 + it1am(isyma,isymi)
               write(lupri,*) kpos, kumat, kvmat, ksval
               ! Do the SVD... Maybe we can use 'S' instead of 'A'
               CALL DGESVD(
     &            'A','A',nvir(isyma),nrhf(isymi),
     &            WORK(kpos),
     &            nvir(isyma), work(ksval),
     &            work(kumat), nvir(isyma),
     &            work(kvmat), nrhf(isymi),
     &            work(kend3), lwrk3, i_info
     &         )
               write(LUPRI,'(3I5)') nvir(isyma), nrhf(isymi), mindim
               if (i_info .ne. 0) then    
                  print *, 'ERROR in DGSVD ', i_info
               endif
               do n = 1, mindim
                  write(LUPRI,'(F16.10)') work(ksval-1+n)
               end do
               write (LUPRI,'(A)') ''
                  
               
               ! work(kumat) now has the left/virtual eigenvectors on
               ! the columns
               ! work(kvmat) now has the right/occupied eigevectors on
               ! the rows

               ! Find the number of NTOs to write
               nfound = 0
               do n = 1, mindim
                  if (work(ksval) .lt. thresh_nto) exit
                  nfound = nfound + 1
               end do
               ! Ensure that we don't write too many, i.e. more than
               ! half the number ofbitals
               nfound = min(nfound, norb(isymi)/2, norb(isyma)/2)

               ! Transform the hole NTOs to AO basis
               call dgemm('N','T', norb(isymi), nfound, nrhf(isymi),
     &                    one,
     &                    work(kcmot+ilmrhf(isymi)), norb(isymi), ! AO x occ-MO matrix
     &                    work(kvmat), nrhf(isymi), ! NTO x occ-MO matrix
     &                    zero, work(knto+icmo(isymi)), nbas(isymi))

               ! Transform the particle NTOs to AO basis
               kpos = knto + icmo(isyma)
     &              + nbas(isyma)*(norb(isyma)-nfound)
               call dgemm('N','N', norb(isyma), nfound, nvir(isyma), 
     &                    one,
     &                    work(kcmot+ilmvir(isyma)), norb(isyma),
     &                    work(kumat), nvir(isyma),
     &                    zero, work(kpos), nbas(isyma))

               ! Put the hole NTO eigenvalues in the occupation vector
               kpos = koccu + iorb(isymi)
               work(kpos:kpos+nfound-1) = - work(ksval:ksval+nfound-1)
               kpos = koccu + iorb(isyma) + (norb(isyma) - nfound)
               work(kpos:kpos+nfound-1) = work(ksval:ksval+nfound-1)

            end do
!         
!           Create filename
!           ---------------
!            
            write(filename,"(A5,I1,A,I1,A7)") 
     &         "exci_", excisym, '_', iexci, ".molden"
!            
!           Call molden_mos to write MOs in molden format
!           ---------------------------------------------
            call gpopen(LUMOLDEN_MOS,trim(filename),'NEW',' ',
     &                  'FORMATTED',idummy,.false.)
            rewind lumolden_mos
            call molden_mos(1, work(knto), work(koccu), work(krbtr),
     &                      work(kucmo),work(korbv))
            call gpclose(LUMOLDEN_MOS,'KEEP')
         end do
      end do
!
!     open original molden_mos.tmp
!     ----------------------------
      call gpopen(LUMOLDEN_MOS,'molden_MOS.tmp','UNKNOWN',' ',
     &            'FORMATTED',idummy,.false.)
!
      return
      end
