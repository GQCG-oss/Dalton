!> @file 
!> Contains the trilevel and atoms starting guess in addition to the construction of the GCbasis se PCCP 2009, 11, 5805-5813 
!> Trilevel and Atoms module
!> \author Branislav Jansik documented by \latexonly T. Kj{\ae}rgaard\endlatexonly
!> \date 2010-03-03
module trilevel_module
use precision
use typedeftype
use typedef
use molecule_type!, only: moleculeinfo, build_atomicmolecule, free_moleculeinfo
use basis_type!, only: basisinfo, basissetinfo, free_basissetinfo
use matrix_module
use Matrix_Operations
use Matrix_Operations_aux!, only: mat_density_from_orbs
use matrix_util
use scfloop_module!, only : scfloop
use lsdalton_fock_module
use BUILDAOBATCH
use lsparameters
use READMOLEFILE
use BUILDBASISSET
use LSTIMING
use memory_handling!, only: mem_alloc,mem_dealloc
use integralinterfaceMod
use II_XC_interfaceModule
use files
use diagonalization
type trilevel_atominfo
     integer :: LUPRI
     integer :: ND,NA ! number of disticnt atoms, number of atoms
     integer, pointer :: NATOM(:)!length ND : FOR A GIVEN UNIQUE ATOM THIS POINTS TO A ATOM OF SAME TYPE
     integer, pointer :: UATOMTYPE(:) !length ND : FOR A GIVEN UNIQUE ATOM THIS IS THE ATOMTYPE IN FULL 
end type trilevel_atominfo

private
public :: trilevel_atominfo_init, trilevel_gcbasis, trilevel_atominfo_free,&
     & trilevel_ATOMS_density, trilevel_full2valence, trilevel_readdens,&
     & trilevel_cmo_valence2full, trilevel_density_valence2full,&
     & trilevel_atominfo, freeVbasis

contains

!> \brief diagonalize a angular moment block of the atomic h1 matrix 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param angular moment
!> \param basis_size size of basis
!> \param bCMO MO coefficients
!> \param F full Fock matrix 
!> \param S full Overlap matrix 
!> \param nbast full dimensions 
subroutine trilevel_diag_per_ang(ang,basis_size,bCMO,F,S,nbast,nb)
implicit none
integer, intent(in)     :: ang, nbast,nb
integer, intent(in)     :: basis_size(:)
real(realk), intent(in) :: F(nbast,nbast), S(nbast,nbast)
real(realk),target      :: bCMO(nb,nb)
!
real(realk), pointer    :: bF(:,:), bS(:,:), eig(:), wrk(:)
integer                 :: i, j, k, istart, info, lwrk,itype
integer,     pointer    :: indexlist(:) 
#ifdef VAR_LSESSL
integer :: ifail(nb),iwrk(5*nb),nfound
real(realk) :: no_ref,tol
real(realk) :: Z(nb,nb)
#endif
 
 info = 0
 itype=1
 call trilevel_indexlist(ang,basis_size,indexlist)

!create subblock matrix
 call mem_alloc(bS,nb,nb)
 call mem_alloc(bF,nb,nb)

 do j=1, nb
  do i=1, nb
       bF(i,j)=F(indexlist(i),indexlist(j))
       bS(i,j)=S(indexlist(i),indexlist(j))
  enddo
 enddo

! diagonalize
!querry optimal work size

 if (nb.gt. 1) then
  lwrk=-1
  call mem_alloc(eig,nb)
  call mem_alloc(wrk,2)
#ifdef VAR_LSESSL
  no_ref=0.0E0_realk
  tol = 0.0E0_realk
  call DSYGVX( itype,'V','A','U', nb,bF,nb,bS,nb,no_ref,no_ref,no_ref,&
   &no_ref,tol,nfound,eig, Z,nb, wrk, lwrk, iwrk, ifail, info)
#else
  call dsygv(itype,'V','U',nb,bF,nb,bS,nb,eig,wrk,lwrk,info)
#endif
  
!run diagonalization
  lwrk = NINT(wrk(1))
  call mem_dealloc(wrk)
  call mem_alloc(wrk,lwrk)
#ifdef VAR_LSESSL
  call DSYGVX( itype,'V','A','U', nb,bF,nb,bS,nb,no_ref,no_ref,no_ref,&
   &no_ref,tol,nfound,eig, Z,nb, wrk, lwrk, iwrk, ifail, info)
  call dcopy(nb*nb,Z,1,bF,1)
#else
  call dsygv(itype,'V','U',nb,bF,nb,bS,nb,eig,wrk,lwrk,info)
#endif
 
  if(info/=0)then
    call lsquit("ERROR(trilevel_diag_per_ang):dsygv(x) failed to converge",-1)
  endif 

  call mem_dealloc(eig)
  call mem_dealloc(wrk)
  bCMO = bF
 else
  bCMO(1,1)=1E0_realk
 endif

!cleanup
call mem_dealloc(indexlist)
call mem_dealloc(bS)
call mem_dealloc(bF)

end subroutine trilevel_diag_per_ang

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt Contains info about SCF optimization
!> \param D density matrix
!> \param H1 one electron contribution to the fock matrix
!> \param F fock matrix
!> \param Etotal energy 
!> \param newlupri logical unit number for output
!> \param newluerr logical unit number for error output
!> \param setting setting structure containing integral info
   SUBROUTINE trilevel_get_fock(opt,D,H1,F,Etotal,newlupri,newluerr,setting)
   ! ===================================================================
   ! di_get_fock obtains total fock matrix and corresponding energy.
   ! WE have to go through the interface to dalton before the fock
   ! evaluator learns how to handle arbitrary-type arrays.
   ! ===================================================================
      use opttype
      IMPLICIT NONE
      type(optItem),intent(in)   :: opt
      TYPE(Matrix), INTENT(IN)   :: D(1) 
      TYPE(Matrix),intent(inout) :: F(1)
      type(lssetting),intent(inout) :: setting
      real(realk), INTENT(INOUT) :: Etotal
      TYPE(Matrix),intent(in)  :: H1
      integer,intent(in) :: newlupri,newluerr
      real(realk)   :: edfty(1), edfty_a, edfty_b
      integer nbast,ndmat
      logical  :: Dsym

!     Two-electron part: G(D)
!
!     Coulomb and exchange
      Dsym = .TRUE.!symmetric Density matrix
      ndmat = 1
      call II_get_Fock_mat(newlupri,newluerr,setting,D,Dsym,F,ndmat,.FALSE.)
      Etotal = trilevel_fockenergy_f(opt%cfg_unres,F(1),D(1),H1)
!     Exchange-correlation
      if (opt%calctype == opt%dftcalc) then
         nbast = D(1)%nrow
         call II_get_xc_fock_mat(newlupri,newluerr,setting,nbast,D,F,Edfty,ndmat)
         Etotal = Etotal + Edfty(1)
      ENDIF

!     Add one-electron part: F(D) = h + G(D)
      call mat_daxpy(1E0_realk,H1,F(1))

   contains      
      double precision function trilevel_fockenergy_F(unres,F,D,H1)
         !E = Tr(h + F)D + hnuc ! No factor since molecular orbitals
         implicit none
         logical, intent(in) :: unres
         TYPE(matrix), intent(in) :: F,D,H1
         double precision :: hD, FD, fac
         double precision, external :: DF_E_ROB_CORR
         integer :: ndim

         ndim=F%nrow
         fac=2E0_realk
         if(unres) fac=1E0_realk
         !Tr(hD)
         hD =       mat_dotproduct(D,H1)
         !Tr(FD)
         FD = 0.5E0_realk*mat_dotproduct(D,F)
         !E(HF) = Tr(h +F)D + hnuc
         trilevel_fockenergy_F = (hd + FD)*fac
       end function trilevel_fockenergy_F

   END SUBROUTINE trilevel_get_fock

!> \brief the grand canonical SCF loop
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt optItem containing info about scf optimization
!> \param D density matrix
!> \param CMO MO coefficients
!> \param H1 one electron contribution to the fock matrix
!> \param F fock matrix
!> \param S Overlap matrix
!> \param ai trilevel_atominfo structure containing info about the unique atoms
!> \param setting setting structure containing integral info
!> \param iatom unique atom index
!> \param newlupri logical unit number for output
!> \param newluerr logical unit number for error output
SUBROUTINE trilevel_gcscfloop(opt,D,CMO,H1,F,S,ai,setting,molecule,basis,iatom,newlupri,newluerr)
  use av_utilities
  use matrix_util
  use opttype
  USE scf_stats
  use linsca_diis 

   implicit none
   type(optItem) :: opt
   type(trilevel_atominfo),intent(in) :: ai
   integer,intent(in) :: iatom
   type(lssetting),intent(inout) :: setting
   type(moleculeinfo),intent(in) :: molecule
   type(basisinfo),intent(in)    :: basis
   TYPE(Matrix),intent(inout)   :: D(1), F(1), S, H1, CMO
   integer,intent(in) :: newlupri,newluerr
!
   TYPE(Matrix)            :: grad
   TYPE(util_HistoryStore) :: queue
   real(realk) :: E, gradnrm
   integer     :: iteration
   logical :: energy_converged,firsttime,cfg_oao_gradnrm
   integer :: itype,nbast,iAO
   type(avItem) :: av
   type(moleculeinfo),target :: atomicmolecule

   itype = ai%UATOMTYPE(iatom)
   nbast = molecule%atom(ai%NATOM(iatom))%nContOrbREG
   cfg_oao_gradnrm = opt%cfg_oao_gradnrm
   opt%cfg_oao_gradnrm = .FALSE.
   call av_set_default_config(av)
   av%lupri = opt%lupri
   av%CFG_averaging = av%CFG_AVG_DIIS
   av%trilevel_gcscf = .true.

   if (opt%calctype == opt%dftcalc) then
      av%cfg_set_type = av%CFG_THR_dft
      av%diis_history_size  = av%cfg_settings(av%CFG_SET_type)%max_history_size
   endif
   call queue_init(av,queue)

   call scf_stats_init(opt)

   call mat_init(grad,H1%nrow,H1%nrow)

!
! GCSCF iterations
!
   CALL build_atomicmolecule(molecule,atomicmolecule,ai%NATOM(iatom),opt%lupri)
   do iAO=1,4
      setting%molecule(iAO)%p => atomicmolecule
      setting%fragment(iAO)%p => atomicmolecule
   enddo

   DO iteration = 1, trilevel_max_linscf_iterations

      call trilevel_get_fock(opt, D, H1, F, E, &
           &newlupri,newluerr,setting)

      call get_AO_gradient(F(1), D(1), S, grad)

      gradnrm = sqrt(mat_sqnorm2(grad))

      call scf_stats_update(iteration,gradnrm,E,opt)

      IF(gradnrm < 7.5E-4_realk) then 
            EXIT
      ENDIF

      CALL add_to_queue(av, F(1), D(1), S, E, grad, queue) 

      call diis(av,queue,D(1),F(1))

      CALL trilevel_get_density_blocked(D(1),CMO,F(1),S,basis,itype,nbast)
   END DO
   !IF not converged then what!
   print *, "gcscf loop done"
  call queue_free(av,queue)
  CALL mat_free(grad)
  call scf_stats_shutdown
  call free_moleculeinfo(atomicmolecule)
  opt%cfg_oao_gradnrm = cfg_oao_gradnrm

END SUBROUTINE trilevel_gcscfloop

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param CMO MO coefficients 
!> \param F fock matrix
!> \param S Overlaps matrix
!> \param itype atom type of the current atom of interest
!> \param nbast number of basis functions for this atom
subroutine trilevel_get_density_blocked(D,CMO,F,S,basis,itype,nbast)
implicit none
type(Matrix)           ::  F,D,S,CMO
type(basisinfo)        :: basis
integer                ::  nAngmom, kmult,ipos, ang, nbast,nb,itype
real(realk),pointer    ::  bCMO(:,:),bD(:,:),occ(:), Fmat(:), Smat(:)
integer, pointer       ::  perm(:), iperm(:), basis_size(:)


  nAngmom = BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%nAngmom

  call trilevel_set_basis_size(BASIS,itype,basis_size) !size of basis for a given angmom

  !build list of permutation to obtain (1S,2S,3S,2Px,3Px,2Py,3Py,2Pz,3Pz,3Pz,...)
  call trilevel_blockdiagonal_permutation(nAngmom,nbast,basis_size,perm) 
  !sorting
  call trilevel_isort(perm,nbast,iperm)

  call mem_alloc(occ,nbast)
  occ = 0E0_realk;
  CMO%elms=0E0_realk;

  kmult=1; ipos =1
  do ang=0,nAngmom-1  

     nb=basis_size(ang+1)
     if (nb.ne. 0) then

     call mem_alloc(bCMO,nb,nb)
     !diagonal angular block => bCMO
     call trilevel_diag_per_ang(ang,basis_size,bCMO,F%elms,S%elms,nbast,nb)
     !set occupation numbers based on element table in ecdata
     call trilevel_set_occ(occ(ipos:nbast),ang,BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%Charge)
     !Build full MO coefficient matrix from the bCMO blocks
     call trilevel_mdiag(CMO%elms,ipos,nbast,bCMO,nb,kmult)

     call mem_dealloc(bCMO)

     endif
     kmult = kmult + 2
  enddo
  !reorder the MO coefficient matrix
  call trilevel_reorder2d(CMO%elms,nbast,iperm)
  !build Density matrix eq. 7 from article
  call trilevel_density_from_orbs(D%elms,occ,CMO%elms,nbast)
  
  call mem_dealloc(basis_size)
  call mem_dealloc(perm)
  call mem_dealloc(iperm)
  call mem_dealloc(occ)

end subroutine trilevel_get_density_blocked

!> \brief build trilevel_blockdiagonal_permutation  
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param nAngmom number of angular moments
!> \param nbast number of basis functions
!> \param basis_size size of basis
subroutine trilevel_blockdiagonal_permutation(nAngmom,nbast,basis_size,permutation)
implicit none
integer, pointer :: permutation(:)
integer          :: nAngmom, nbast
integer          :: basis_size(nAngmom)
integer          :: i,j,k,l,istart,iend,ioff

  call mem_alloc(permutation,nbast)

  k=1;l=1
  istart=1
  do i=1,nAngmom
   iend = istart  + k*basis_size(i) -1
   do ioff=0,(k -1)
   do j=(istart+ioff),iend,k
      permutation(l) = j 
      l=l+1
   enddo
   enddo
   istart= iend+1
   k=k+2
  enddo

end subroutine trilevel_blockdiagonal_permutation

!> \brief sorting routine
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param arr  vector of integers to be sorted 
!> \param n length of vector
subroutine trilevel_isort(arr,n,iperm)
implicit none
integer, pointer :: iperm(:)
integer          :: n
integer          :: arr(n),a,b
integer          :: i,j

 call mem_alloc(iperm,n)

 do i=1,n
  iperm(i)=i
 enddo

 do j=2,n
  a=arr(iperm(j)); b=iperm(j)
  do i=j-1,1,-1
    if (arr(iperm(i)).le.a) goto 10
!   arr(i+1)=arr(i)
    iperm(i+1)=iperm(i)
  enddo
  i=0
10 continue 
! arr(i+1)=a
  iperm(i+1)=b
 enddo

end subroutine trilevel_isort

!> \brief reorder routine
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param A matrix to be reordered
!> \param n dimension
!> \param perm reorder according to this vector 
subroutine trilevel_reorder2d(A,n,perm)
implicit none
integer     :: n,i
real(realk) :: A(n,n)
real(realk) :: tmp(n,n)
integer     :: perm(n)


  do i=1,n
     tmp(:,i) = A(:,perm(i))
  enddo

  do i=1,n
     A(i,:)   = tmp(perm(i),:)
  enddo


end subroutine trilevel_reorder2d

!> \brief calculates the size of the basisset
!> \author Branislav Jansik
!> \date 2010-03-03
subroutine trilevel_set_basis_size(basis,itype,basis_size)
implicit none
integer, pointer :: basis_size(:)
Type(basisinfo)  :: basis
integer          :: itype
!
integer                          :: i

call mem_alloc(basis_size,BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%nAngmom)

do i=1, BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%nAngmom
   basis_size(i)= BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%SHELL(i)%norb
enddo

end subroutine trilevel_set_basis_size

!> \brief Eq. 7 from PCCP 2009, 11, 5805-5813
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param occ occupation numbers
!> \param CMO MO coefficients
!> \param n dimension of density matrix
subroutine trilevel_density_from_orbs(D,occ,CMO,n)
implicit none
integer     :: n, i, j
real(realk) :: occ(n), CMO(n,n), D(n,n)
real(realk) :: tmp(n,n)

 do j=1, n
      tmp(:,j)=occ(j)*CMO(:,j)
 enddo

 call dgemm('n','t',n,n,n,1E0_realk,tmp,n,CMO,n,0E0_realk,D,n)


end subroutine trilevel_density_from_orbs

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D output density matrix  
!> \param istart start index
!> \param nbast dimension of D
!> \param bD 
!> \param n dimension of bD
!> \param kmult
subroutine trilevel_mdiag(D,istart,nbast,bD,n,kmult)
implicit none
integer     :: istart,nbast,n,kmult,i
real(realk) :: bD(n,n),D(nbast,nbast)

 do i=1, kmult
  D(istart:istart+n-1,istart:istart+n-1)=bD
  istart=istart+n
 enddo

end subroutine trilevel_mdiag


!> \brief set occupation from elementtabel (from ecdata)
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param occ occupation number 
!> \param ang angular moment
!> \param charge charge of the atom of interest 
subroutine trilevel_set_occ(occ,ang,charge)
use EcData
implicit none
real(realk)   :: occ(:)
integer       :: ang, charge


 select case(ang)
 case(0); occ(1:ElementTable(charge)%nocc_s)=&
         &ElementTable(charge)%occ_s(1:ElementTable(charge)%nocc_s);
 case(1); occ(1:ElementTable(charge)%nocc_p)=&
         &ElementTable(charge)%occ_p(1:ElementTable(charge)%nocc_p);
 case(2); occ(1:ElementTable(charge)%nocc_d)=&
         &ElementTable(charge)%occ_d(1:ElementTable(charge)%nocc_d);
 case(3); occ(1:ElementTable(charge)%nocc_f)=&
         &ElementTable(charge)%occ_f(1:ElementTable(charge)%nocc_f);

 case default
 end select
end subroutine trilevel_set_occ

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param 
!> \param 
!> \param 
!> \param 
!> \param 
subroutine trilevel_indexlist(ang,basis_size,indexlist)
implicit none
integer, pointer :: indexlist(:)
integer          :: ang, basis_size(:)
integer          :: nb, i,j, k, istart
 
 nb =  basis_size(ang+1)
 k=1;
 istart=1
 do i=1,ang
   istart = istart + k*basis_size(i)
   k = k+2;
 enddo

!create list of indexes of elements with same angular and lowest magnetic number (i.e. px1,px2,px3)
 call mem_alloc(indexlist,nb)

 j = 1
 do i=istart, istart -1 + (k*nb), k
   indexlist(j)=i  
   j=j+1
 enddo

end subroutine trilevel_indexlist

!> \brief change the input basis in the ls%setting to the grand canonical basis Eq. 8
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ls lsitem structure containing info about the integral eval and basis
!> \param CMO MO coefficients
!> \param ai trilevel_atominfo containing info about the unique atoms
!> \param iatom the unique atom index
subroutine trilevel_convert_ao2gcao(ls,CMO,nbast,ai,iatom,lupri)
implicit none
integer :: lupri
TYPE(lsitem) :: ls!,atomic_ls
type(trilevel_atominfo) :: ai
real(realk)         :: CMO(nbast,nbast)
real(realk),pointer :: bCMO(:,:)
integer             :: nAngmom, nbast,ang,nb,i,j,itype,icharge
integer             :: nprim, norb,iatom,nrow,ncol
integer, pointer    :: basis_size(:), indexlist(:)
real(realk), pointer:: contrCoeffs(:)
real(realk), pointer :: CCtmp(:,:)
type(shell),pointer  :: shell2
type(segment),pointer :: GCtransSegment
integer :: icont,iprim,iprimLoc,iContLoc,iseg,ielm,ip1,ic1

  itype = ai%UATOMTYPE(iatom) !type of distinct atom in full
  nAngmom = ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%nAngmom

  call trilevel_set_basis_size(ls%setting%BASIS(1)%p,itype,basis_size)

 do ang=0, nAngmom-1
   shell2 => ls%input%BASIS%BINFO(REGBASPARAM)%&
                &ATOMTYPE(itype)%SHELL(ang+1)
   nb=basis_size(ang+1)
   if (nb.eq. 0) cycle

   call trilevel_indexlist(ang,basis_size,indexlist)

   call mem_alloc(bCMO,nb,nb)

   do j=1, nb
    do i=1, nb
       bCMO(i,j)=CMO(indexlist(i),indexlist(j))
    enddo
   enddo

!   WRITE(lupri,*)'The transformation matrix itype=',itype,'ang=',ang
!   call ls_output(bCMO,1,nb,1,nb,nb,nb,1,lupri)
!   WRITE(lupri,*)'The Original block:'
!   call ls_output(shell2%segment(1)%elms,1,shell2%nprim,1,shell2%norb,&
!        &shell2%nprim,shell2%norb,1,lupri)
   IF(.NOT. ls%input%BASIS%wbasis(gctbasparam))THEN
      call lsquit('GCtrans basis not call mem_allocd in trilevel_convert_ao2gcao',lupri)
   ENDIF
   !We always save the transformation matrix in GCtrans
   !so that we transform back and forth between 
   !AO and GCAO in the integral routine
   GCtransSegment => ls%input%BASIS%binfo(gctbasparam)%&
        &ATOMTYPE(itype)%SHELL(ang+1)%segment(1)
   call dcopy(nb*nb,bCMO,1,GCtransSegment%elms,1)


   IF(ls%input%DALTON%NOGCINTEGRALTRANSFORM)THEN
      !we transform the input basis
      ls%input%basis%BINFO(REGBASPARAM)%GCbasis = .TRUE.
      IF(ls%input%basis%BINFO(REGBASPARAM)%Gcont)THEN         
         !General contracted Case or no integraltransform
         !we transform the input basis - so that 
         !we do not need to transform back and forth
         nprim     = shell2%nprim
         norb      = shell2%norb
         call mem_alloc(CCtmp,nprim,norb)
         CCtmp = reshape(shell2%segment(1)%elms(1:nprim*norb), (/ nprim, norb /))
         CCtmp = matmul(CCtmp,bCMO)   
         shell2%segment(1)%elms = reshape(CCtmp, (/ nprim*norb /))
         call mem_dealloc(CCtmp)
      ELSE
         !     We transform the segmentet input basis
         !     transform_basis_to_GCAO
         nprim     = shell2%nprim
         norb      = shell2%norb
         call mem_alloc(CCtmp,nprim,norb)
         do iCont=1,norb
            do iprim=1,nprim
               CCtmp(iprim,iCont) = 0E0_realk
            enddo
         enddo
         iPrim=0
         iCont=0      
         DO iseg = 1,shell2%nsegments
            ielm = 0
            iContLoc = iCont
            do ic1=1,shell2%segment(iseg)%ncol
               iContLoc = iContLoc+1 
               iprimLoc = iPrim            
               do ip1=1,shell2%segment(iseg)%nrow
                  iprimLoc = iprimLoc+1 
                  ielm = ielm +1
                  CCtmp(iPrimLoc,iContLoc) = shell2%segment(iseg)%elms(ielm)
               enddo
            enddo
            iPrim = iPrim + shell2%segment(iseg)%nrow
            iCont = iCont + shell2%segment(iseg)%ncol
         ENDDO
         
         CCtmp = matmul(CCtmp,bCMO)   
         call mem_dealloc(shell2%segment(1)%elms)
         call mem_alloc(shell2%segment(1)%elms,nprim*norb)
         shell2%segment(1)%elms = reshape(CCtmp, (/ nprim*norb /))
         
         iPrim=0
         DO iseg = 1,shell2%nsegments
            do ip1=1,shell2%segment(iseg)%nrow
               iprim = iprim+1 
               CCtmp(iPrim,1) = shell2%segment(iseg)%exponents(iP1)
            enddo
         enddo
         call mem_dealloc(shell2%segment(1)%exponents)
         call mem_alloc(shell2%segment(1)%exponents,nprim)
         DO iPrim=1,nPrim
            shell2%segment(1)%exponents(iPrim) = CCtmp(iPrim,1)
         ENDDO
         call mem_dealloc(shell2%segment(1)%UCCelms)
         call mem_alloc(shell2%segment(1)%UCCelms,nprim*norb)
         shell2%segment(1)%UCCelms = 0E0_realk
         DO iseg = 2,shell2%nsegments
            call mem_dealloc(shell2%segment(iseg)%elms)
            call mem_dealloc(shell2%segment(iseg)%UCCelms)
            call mem_dealloc(shell2%segment(iseg)%exponents)
         enddo
         shell2%nsegments=1      
         shell2%segment(1)%ncol = norb
         shell2%segment(1)%nrow = nprim
         call mem_dealloc(CCtmp) 
      ENDIF
   ELSE
      IF(ls%input%basis%BINFO(REGBASPARAM)%Gcont)THEN
         !General contracted Case 
         !we transform the input basis - so that 
         !we do not need to transform back and forth
         nprim     = shell2%nprim
         norb      = shell2%norb
         call mem_alloc(CCtmp,nprim,norb)
         CCtmp = reshape(shell2%segment(1)%elms(1:nprim*norb), (/ nprim, norb /))
         CCtmp = matmul(CCtmp,bCMO)   
         shell2%segment(1)%elms = reshape(CCtmp, (/ nprim*norb /))
         call mem_dealloc(CCtmp)
         !      WRITE(6,*)'The transformed block:'
         !      call ls_output(shell2%segment(1)%UCCelms,1,shell2%nprim,1,shell2%norb,&
         !           &shell2%nprim,shell2%norb,1,6)
         ls%input%basis%BINFO(REGBASPARAM)%GCbasis = .TRUE.
      ELSE
         !Segmented contracted Case
         !so we transform back and forth. which means that we do nothing
         ls%input%basis%BINFO(REGBASPARAM)%GCbasis = .FALSE.
      ENDIF
   ENDIF
   call mem_dealloc(bCMO)
   call mem_dealloc(indexlist)
 enddo 

 call mem_dealloc(basis_size)

! print*,'print BASIS  after ao2gcao'
! call print_basissetinfo(6,ls%input%BASIS%BINFO(REGBASPARAM))

end subroutine trilevel_convert_ao2gcao

subroutine trilevel_ALLOC_SYNC_VBASIS(VBASISINFO,BASISINFO,GCtrans,integralGCtrans,&
     & GCtransAlloc,GCbasis,MOLECULE,lupri)
IMPLICIT NONE
TYPE(BASISSETINFO) :: BASISINFO,VBASISINFO,GCtrans
TYPE(MOLECULEINFO) :: MOLECULE
logical            :: integralGCtrans,GCtransAlloc,GCbasis
INTEGER            :: I,J,K,L,nrow,nsize,icharge,maxcharge,ncol,lupri
real(realk),pointer:: CCtmp(:,:),Exponents(:),bCMO(:,:) 
integer :: iprim,icont,iseg,ielm,icontloc,ic1,ip1,iprimloc,nrow2,ncol2,R
  IF(integralGCtrans)THEN
     !it is not possible to transform back and forth on level 2 
     call lsquit('VBASIS ALWAYS TRANSFORMED TO GC BASIS no integraltransform',-1)
  ENDIF
  IF (BASISINFO%natomtypes.eq. 0) THEN
     VBASISINFO%natomtypes=0
     call lsquit('Vbasis error in trilevel_ALLOC_SYNC_VBASIS',lupri)
  ENDIF
  VBASISINFO%natomtypes = BASISINFO%natomtypes
  CALL MEM_ALLOC(VBASISINFO%ATOMTYPE,BASISINFO%natomtypes)
!  VBASISINFO%ATOMTYPE = BASISINFO%ATOMTYPE  
  maxcharge = 0
  VBASISINFO%nAtomtypes = BASISINFO%nAtomtypes
  DO J=1,BASISINFO%nAtomtypes
     icharge = BASISINFO%ATOMTYPE(J)%charge
     maxcharge = MAX(maxcharge,icharge)
     !NO need to call mem_alloc SHELL
     call nullifyAtomType(VBASISINFO%ATOMTYPE(J))
     VBASISINFO%ATOMTYPE(J)%nAngmom = BASISINFO%ATOMTYPE(J)%nAngmom
     VBASISINFO%ATOMTYPE(J)%ToTnorb = BASISINFO%ATOMTYPE(J)%ToTnorb
     VBASISINFO%ATOMTYPE(J)%ToTnprim = BASISINFO%ATOMTYPE(J)%ToTnprim
     VBASISINFO%ATOMTYPE(J)%Charge = BASISINFO%ATOMTYPE(J)%Charge
     VBASISINFO%ATOMTYPE(J)%NAME = BASISINFO%ATOMTYPE(J)%NAME
     DO K=1,BASISINFO%ATOMTYPE(J)%nAngmom
        !NO need to call mem_alloc segments
        !we transform the basis so we need to allocate it as if it is general contracted
        VBASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments = 1
        !this is a general contracted basis, so nsegments = 1
        L=1
        nrow=BASISINFO%ATOMTYPE(J)%SHELL(K)%nprim
        ncol=BASISINFO%ATOMTYPE(J)%SHELL(K)%norb
        nsize=nrow*ncol
        CALL MEM_ALLOC(VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms,nSIZE)
        CALL MEM_ALLOC(VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%UCCelms,nSIZE)
        CALL MEM_ALLOC(VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents,nrow)
        Exponents => VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents
        VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(1)%nrow = nrow
        VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(1)%ncol = ncol
        VBASISINFO%ATOMTYPE(J)%SHELL(K)%nprim = nrow 
        VBASISINFO%ATOMTYPE(J)%SHELL(K)%norb = ncol
        call mem_alloc(CCtmp,nrow,ncol)
        !build elms
        CCtmp = 0E0_realk
        iPrim=0
        iCont=0      
        DO iseg = 1,BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
           ielm = 0
           iContLoc = iCont
           do ic1=1,BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(iseg)%ncol
              iContLoc = iContLoc+1 
              iprimLoc = iPrim            
              do ip1=1,BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(iseg)%nrow
                 iprimLoc = iprimLoc+1 
                 ielm = ielm +1
                 CCtmp(iPrimLoc,iContLoc) = BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(iseg)%elms(ielm)
                 Exponents(iPrim+ip1) = BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(iseg)%Exponents(ip1)
              enddo
           enddo
           iPrim = iPrim + BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(iseg)%nrow
           iCont = iCont + BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(iseg)%ncol
        ENDDO
!        WRITE(6,*)'VBASIS untransfomed',nrow,ncol
!        call ls_output(CCtmp,1,nrow,1,ncol,nrow,ncol,1,6)
        IF(GCtransAlloc)THEN
           !GCbasis requested
           IF(GCbasis)THEN
              !the input basis have been transformed so no need to modify valence basis
              VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(1)%elms = reshape(CCtmp, (/ nrow*ncol /))
           ELSE
              !the input basis have not been transformed so we need to transform valence basis
              !WRITE(6,*)'VBASIS untransfomed',nrow,ncol
              !call ls_output(CCtmp,1,nrow,1,ncol,nrow,ncol,1,6)
              nrow2 = GCtrans%ATOMTYPE(J)%SHELL(K)%nprim
              ncol2 = GCtrans%ATOMTYPE(J)%SHELL(K)%norb
              call mem_alloc(bCMO,nrow2,ncol2)
              bCMO = reshape(GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%elms,(/ nrow2,ncol2 /))
              !WRITE(6,*)'transform matrix',nrow2,ncol2
              !call ls_output(bCMO,1,nrow2,1,ncol2,nrow2,ncol2,1,6)
              CCtmp = matmul(CCtmp,bCMO) 
              VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(1)%elms = reshape(CCtmp, (/ nrow*ncol /))
              !print*,'GCtransAlloc true so we transform GC basis and we keep ',VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(1)%elms
              call mem_dealloc(bCMO)           
           ENDIF
        ELSE
           !GCbasis not requested 
           VBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(1)%elms = reshape(CCtmp, (/ nrow*ncol /))
        ENDIF
!        WRITE(6,*)'VBASIS transformed'
!        call ls_output(CCtmp,1,nrow,1,ncol,nrow,ncol,1,6)
        call mem_dealloc(CCtmp)
     ENDDO
  ENDDO
  IF(BASISINFO%labelindex.EQ. 0)THEN
     VBASISINFO%labelindex = 0
     CALL MEM_ALLOC(VBASISINFO%chargeindex,BASISINFO%nchargeindex,.TRUE.)
     VBASISINFO%nchargeindex = BASISINFO%nchargeindex
     VBASISINFO%chargeindex(0:BASISINFO%nchargeindex) = &
          &BASISINFO%chargeindex(0:BASISINFO%nchargeindex)
  ELSE
     VBASISINFO%labelindex = VALBasParam
     R = BASISINFO%labelindex
     DO J=1,molecule%natoms
        molecule%atom(J)%idtype(VALBasParam) = molecule%atom(J)%idtype(R)
     ENDDO
  ENDIF
  VBASISINFO%label = 'VALENCE  '

END SUBROUTINE trilevel_ALLOC_SYNC_VBASIS

subroutine trilevel_ALLOC_SYNC_GCTRANS(GCtrans,REGULAR,MOLECULE)
implicit none
TYPE(MOLECULEINFO) :: MOLECULE
TYPE(BASISSETINFO),intent(in) :: REGULAR
TYPE(BASISSETINFO),intent(inout) :: GCtrans
INTEGER   :: I,J,K,L,nsize,icharge,maxcharge,ncol,KK,R

  IF (REGULAR%natomtypes.eq. 0) THEN
     CALL LSQUIT('Try to build GCtrans for natomtypes=0',6)
  ENDIF

  GCtrans%nChargeindex = 0
  GCtrans%natomtypes = REGULAR%natomtypes
  CALL MEM_ALLOC(GCtrans%ATOMTYPE,REGULAR%natomtypes)
  maxcharge = 0
  GCtrans%nAtomtypes = REGULAR%nAtomtypes
  DO J=1,REGULAR%nAtomtypes
     icharge = REGULAR%ATOMTYPE(J)%charge
     maxcharge = MAX(maxcharge,icharge)
     !NO need to call mem_alloc SHELL
     call nullifyAtomType(GCtrans%ATOMTYPE(J))
     GCtrans%ATOMTYPE(J)%nAngmom = REGULAR%ATOMTYPE(J)%nAngmom
     GCtrans%ATOMTYPE(J)%Charge  = REGULAR%ATOMTYPE(J)%Charge
     DO K=1,REGULAR%ATOMTYPE(J)%nAngmom
        !NO need to call mem_alloc segments
        GCtrans%ATOMTYPE(J)%SHELL(K)%nsegments = 1
        ncol=REGULAR%ATOMTYPE(J)%SHELL(K)%norb
        GCtrans%ATOMTYPE(J)%SHELL(K)%nOrb = ncol
        GCtrans%ATOMTYPE(J)%SHELL(K)%nprim = ncol
        GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%ncol = ncol
        GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%nrow = ncol
        nsize=ncol*ncol
        CALL MEM_ALLOC(GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%elms,nSIZE)
        CALL MEM_ALLOC(GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%UCCelms,nSIZE)
        do KK=1,nsize
           GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%UCCelms(KK) = 0E0_realk
        enddo
        CALL MEM_ALLOC(GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%Exponents,ncol)
        do KK=1,ncol
           GCtrans%ATOMTYPE(J)%SHELL(K)%segment(1)%Exponents(KK) = 0E0_realk
        enddo
     ENDDO
  ENDDO
  IF(REGULAR%labelindex.EQ. 0)THEN
     GCtrans%labelindex = 0
     CALL MEM_ALLOC(GCtrans%chargeindex,REGULAR%nchargeindex,.TRUE.)
     GCtrans%nchargeindex = REGULAR%nchargeindex
     GCtrans%chargeindex(0:REGULAR%nchargeindex) = &
          &REGULAR%chargeindex(0:REGULAR%nchargeindex)
  ELSE
     GCtrans%labelindex = GCTBasParam 
     R = REGULAR%labelindex
     DO J=1,molecule%natoms
        molecule%atom(J)%idtype(GCTBasParam) = molecule%atom(J)%idtype(R)
     ENDDO
  ENDIF
  GCtrans%label = 'GCTRANS  '

end subroutine trilevel_ALLOC_SYNC_GCTRANS
!> \brief loop over all atoms and construct the grand canonical basis 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt optItem containing info about scf optimization
!> \param ls lsitem structure containing integral,molecule,basis info
!> \param ai trilevel_atominfo structure containing info about the unique atoms
SUBROUTINE trilevel_gcbasis(opt,ls,ai,LUPRI,LUERR)
use READMOLEFILE
use BUILDBASISSET
use opttype
use precision
use memory_handling
use typedeftype, only: lssetting, lsitem
use molecule_typetype, only: moleculeinfo
use basis_typetype, only: basisinfo, basissetinfo
IMPLICIT NONE
type(optItem)       :: opt
INTEGER             :: I,LUPRI,LUERR,IPRINT
TYPE(lsitem),intent(inout) :: ls
type(trilevel_atominfo) :: ai
TYPE(lsitem),pointer :: atomic_ls
Type(Matrix)         :: F(1), H1, S , D(1), CMO
integer              :: nbast, len,iAO,itype,igrid
type(moleculeinfo),target :: atomicmolecule
TYPE(lssetting)           :: atomicSetting
logical :: integraltransformGC

Write(ls%lupri,*)''
Write(ls%lupri,*)'A set of atomic calculations are performed in order to construct'
Write(ls%lupri,*)'the Grand Canonical Basis set (see PCCP 11, 5805-5813 (2009))'
Write(ls%lupri,*)'as well as JCTC 5, 1027 (2009)'
Write(ls%lupri,*)'This is done in order to use the TRILEVEL starting guess and '
Write(ls%lupri,*)'perform orbital localization'
Write(ls%lupri,*)'This is Level 1 of the TRILEVEL starting guess and is performed per default.'
Write(ls%lupri,*)'The use of the Grand Canonical Basis can be deactivated using .NOGCBASIS'
Write(ls%lupri,*)'under the **GENERAL section. This is NOT recommended if you do TRILEVEL '
Write(ls%lupri,*)'or orbital localization.'
Write(ls%lupri,*)''
!> Contains the trilevel and atoms starting guess in addition to the construction of the GCbasis se PCCP 2009, 11, 5805-5813 


!atomic calc is always done in AO basis and not transformed
integraltransformGC = ls%setting%integraltransformGC
ls%setting%integraltransformGC = .FALSE.
! We call mem_alloc and initiate the ls%input%BASIS%binfo(gctbasparam) in this routine
! This is needed whenever we go from AO to GC baiss and back.
ls%input%basis%WBASIS(GCTBasParam) = .TRUE.
CALL trilevel_ALLOC_SYNC_GCTRANS(ls%input%BASIS%BINFO(GCTBasParam),&
     & ls%input%basis%BINFO(RegBasParam),ls%input%molecule)
do i=1, ai%ND
   itype = ai%UATOMTYPE(i)
   IF(ls%input%molecule%atom(ai%NATOM(i))%pointcharge)CYCLE
   call io_free(ls%setting%IO)
   call io_init(ls%setting%IO)
   !print statements
   len = len_trim(ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%NAME)
   write(lupri,*)
   write (lupri,'(1X,A,1X,A,1X,A,1X,I3)') 'Level 1 atomic calculation on', &
  & ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%NAME(1:len), 'Charge', &
  & ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%Charge
   write(lupri,*) '================================================'
   write (*,'(1X,A,1X,A,1X,A,1X,I3)') 'Level 1 atomic calculation on', &
  & ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%NAME(1:len), 'Charge', &
  & ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%Charge
 
   nbast = ls%input%molecule%atom(ai%NATOM(i))%nContOrbREG
   CALL mat_init(F(1),nbast,nbast)
   CALL mat_init(H1,nbast,nbast)
   CALL mat_init(S,nbast,nbast)
   CALL mat_init(D(1),nbast,nbast)
   CALL mat_init(CMO,nbast,nbast)
 
   !we change the setting to point to the atom of interest 
   CALL build_atomicmolecule(ls%input%molecule,atomicmolecule,ai%NATOM(i),lupri)
   CALL typedef_init_setting(atomicSetting)
   call typedef_set_default_setting(atomicSetting,ls%input)
   do iAO=1,4
      atomicSetting%molecule(iAO)%p => atomicmolecule
      atomicSetting%fragment(iAO)%p => atomicmolecule
   enddo
   !deactivate density fitting, FMM and screening - provide no speedup and just complicates things
   atomicSetting%scheme%admm_exchange = .FALSE.
   atomicSetting%scheme%densfit = .FALSE.
   atomicSetting%scheme%df_k = .FALSE.
   atomicSetting%scheme%PARI_J = .FALSE.
   atomicSetting%scheme%PARI_K = .FALSE.
   atomicSetting%scheme%MOPARI_K = .FALSE.
   atomicSetting%scheme%FMM = .FALSE.
   atomicSetting%scheme%MBIE_SCREEN = .FALSE.
   atomicSetting%scheme%CS_SCREEN = .FALSE.
   atomicSetting%scheme%PS_SCREEN = .FALSE.
   atomicSetting%scheme%LINK = .FALSE.
   atomicSetting%scheme%DFT%CS00 = .FALSE.
   atomicSetting%scheme%DFT%LB94 = .FALSE.
   atomicSetting%scheme%DFT%griddone = 0 !the DFT grid should be made for all unique atoms
   igrid = atomicSetting%scheme%DFT%igrid
   atomicSetting%scheme%DFT%GridObject(igrid)%griddone = 0

   !build atomic overlap
   CALL II_get_overlap(lupri,luerr,atomicSetting,S)
   !build atomic h1
   CALL II_get_h1(lupri,luerr,atomicSetting,H1)
    !build Density matrix eq. 7 from article
   CALL trilevel_get_density_blocked(D(1),CMO,H1,S,ls%input%basis,itype,nbast)
   !use Density as start guess for a SCF convergence for this atom
   call trilevel_gcscfloop(opt,D,CMO,H1,F,S,ai,atomicSetting,ls%input%molecule,ls%input%basis,&
      &                    i,lupri,luerr)
  !use the MO coefficients from the gcscf loop to create the gcao basis (Eq. 8)
  !WARNING: this changes the input basis in ls to the grand canonical basis
  call trilevel_convert_ao2gcao(ls,CMO%elms,nbast,ai,i,lupri)

  CALL typedef_free_setting(atomicSetting)
  call free_moleculeinfo(atomicmolecule)
  CALL mat_free(F(1))
  CALL mat_free(H1)
  CALL mat_free(S)
  CALL mat_free(D(1))
  CALL mat_free(CMO)
enddo
!reverted back
ls%setting%integraltransformGC = integraltransformGC
ls%optlevel = 3

end subroutine trilevel_gcbasis

!> \brief initialise the trilevel_atominfo containing info about the unique atoms
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ai trilevel_atominfo containing info about the unique atoms to be built
!> \param ls lsitem structure containing the entire molecule and basis
!> \param lupri output logical unit number
subroutine trilevel_atominfo_init(ai,ls,LUPRI)
IMPLICIT NONE
TYPE(trilevel_atominfo) :: ai !Atomic Info
TYPE(lsitem) :: ls
INTEGER             :: LUPRI,IPRINT,I

ai%LUPRI = LUPRI
ai%ND = ls%input%BASIS%BINFO(REGBASPARAM)%nAtomtypes
ai%NA = ls%input%MOLECULE%Natoms
CALL MEM_ALLOC(ai%NATOM,ai%ND)
CALL MEM_ALLOC(ai%UATOMTYPE,ai%ND)

IPRINT=ls%input%dalton%BASPRINT
CALL BUILD_DISTINCT_ATOMS(LUPRI,ai%ND,ai%NA,ai%NATOM,ai%UATOMTYPE,ls,IPRINT)

end subroutine trilevel_atominfo_init

!> \brief deallocation routines for freeing the trilevel_atominfo containing info about the unique atoms
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ai trilevel_atominfo containing info about the unique atoms to be built
subroutine trilevel_atominfo_free(ai)
  use daltoninfo
IMPLICIT NONE
TYPE(trilevel_atominfo) :: ai !Atomic Info
CALL MEM_DEALLOC(ai%NATOM)
CALL MEM_DEALLOC(ai%UATOMTYPE) 
end subroutine trilevel_atominfo_free

!> \brief free the valens basis - we need to set the nAngmom correctly
!> \author Thomas Kjaergaard
!> \date 2011
!> \param vbasis basisset structure contains the valens basis - the level 2 basis
!> \param lsbasis the full basis
subroutine freeVbasis(lsBASIS)
implicit none
TYPE(BASISINFO) :: lsBASIS
!
integer :: I,J,K,L,nangmom1,nAngmom2,nsegments

DO J=1,lsBASIS%BINFO(VALBASPARAM)%natomtypes
   nangmom1=lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%nAngmom
   nangmom2=lsBASIS%BINFO(REGBASPARAM)%ATOMTYPE(J)%nAngmom
   DO K=nAngmom1+1,nAngmom2
      nsegments=lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%SHELL(K)%nsegments
      IF(nsegments .NE. 0)THEN
         DO L=1,nsegments
            if (.not.ASSOCIATED(lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)) then
               print*,'memory previously released!!'
               call lsquit('Error in FREE_BASISSETINFO1 - memory previously released',-1)
            endif
            CALL MEM_DEALLOC(lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)
            if (.not.ASSOCIATED(lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%SHELL(K)%segment(L)%UCCelms)) then
               print*,'memory previously released!!'
               call lsquit('Error in FREE_BASISSETINFO1 - memory previously released',-1)
            endif
            CALL MEM_DEALLOC(lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%SHELL(K)%segment(L)%UCCelms)
            if (.not.ASSOCIATED(lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)) then
               print*,'memory previously released!!'
               call lsquit('Error in FREE_BASISSETINFO2 - memory previously released',-1)
            endif
            CALL MEM_DEALLOC(lsBASIS%BINFO(VALBASPARAM)%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)
         ENDDO
      ENDIF
   ENDDO
ENDDO
end subroutine freeVbasis

!> \brief build the valens basis as a subset of full, ls contains on input the full basis and as output the valens basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param vbasis basisset structure contains the valens basis - the level 2 basis
!> \param ls lsitem structure containing integral,molecule,basis info
!> \param ai trilevel_atominfo structure containing info about the unique atoms
subroutine trilevel_full2valence(ls,ai,lupri)
use READMOLEFILE
use BUILDBASISSET
use BUILDAOBATCH
use EcData
IMPLICIT NONE
INTEGER             :: I,LUPRI,IPRINT
TYPE(lsitem) :: ls
TYPE(trilevel_atominfo) :: ai
!TYPE(lsitem),pointer :: atomic_ls(:)
integer             :: nbast,itype,nAngmom,charge,ang,nb,nocc
integer             :: icharge,jatom,nsegments,tmpindex,seg,iprim
integer, pointer    :: basis_size(:)
integer :: len,R,J
CALL trilevel_ALLOC_SYNC_VBASIS(ls%input%basis%BINFO(VALBasParam),&
     & ls%input%basis%BINFO(RegBasParam),ls%input%basis%BINFO(GCTBasParam),&
     & ls%setting%integraltransformGC,ls%input%basis%wbasis(gctbasparam),&
     & ls%input%basis%BINFO(RegBasParam)%GCbasis,ls%input%MOLECULE,lupri)

R = ls%input%basis%BINFO(VALBasParam)%labelindex
IF(R.NE.0)THEN
   do i=1, ls%setting%nAO   
      DO J=1,ls%setting%MOLECULE(i)%p%natoms
         ls%setting%MOLECULE(i)%p%atom(J)%idtype(VALBasParam) = ls%input%MOLECULE%atom(J)%idtype(R)
      ENDDO
   enddo
ENDIF

ls%input%basis%WBASIS(VALBasParam) = .TRUE.
do i=1, ai%ND
   itype = ai%UATOMTYPE(i)
   !points to function which calculates the size of the basis
   call trilevel_set_basis_size(ls%setting%BASIS(1)%p,itype,basis_size)
   jatom = ai%NATOM(i) !an atom in the full input molecule
   nAngmom = ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%nAngmom
   charge  = ls%input%BASIS%BINFO(REGBASPARAM)%ATOMTYPE(itype)%Charge
   
   ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%ToTnorb = &
        & ElementTable(charge)%nocc_s + ElementTable(charge)%nocc_p &
        &+ElementTable(charge)%nocc_d + ElementTable(charge)%nocc_f
   if (ElementTable(charge)%nocc_s .gt. 0 ) ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%nAngmom = 1
   if (ElementTable(charge)%nocc_p .gt. 0 ) ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%nAngmom = 2
   if (ElementTable(charge)%nocc_d .gt. 0 ) ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%nAngmom = 3
   if (ElementTable(charge)%nocc_f .gt. 0 ) ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%nAngmom = 4
   
   do ang=0, nAngmom-1
      
      nb=basis_size(ang+1)
      if (nb.eq. 0) cycle
      
      select case(ang)
      case(0); nocc=ElementTable(charge)%nocc_s
      case(1); nocc=ElementTable(charge)%nocc_p/3
      case(2); nocc=ElementTable(charge)%nocc_d/5
      case(3); nocc=ElementTable(charge)%nocc_f/7
      case default; nocc=0;
      end select
      IF(ls%setting%integraltransformGC)CALL LSQUIT('the Valence basis must be transformed to GCbasis',-1)
      !the Valence basis is a GCbasis 
      ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%SHELL(ang+1)%norb = nocc
      ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%SHELL(ang+1)%segment(1)%ncol = nocc
      if (nocc.eq. 0) then
         ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%SHELL(ang+1)%nprim = 0
         ls%input%basis%BINFO(VALBASPARAM)%ATOMTYPE(itype)%SHELL(ang+1)%segment(1)%nrow = 0
      endif
   enddo
   call mem_dealloc(basis_size)
enddo
call determine_nbast(ls%input%MOLECULE,ls%input%basis%BINFO(VALBASPARAM),&
      &ls%setting%scheme%DoSpherical,ls%setting%scheme%uncont)

!print*,'trilevel_full2valence vbasis MERGE'
!call print_basissetinfo(6,ls%input%basis%BINFO(VALBASPARAM))
!print*,'trilevel_full2valence BINFO(REGBASPARAM)'
!call print_basissetinfo(6,ls%input%basis%BINFO(REGBASPARAM))

END SUBROUTINE trilevel_full2valence

!> \brief build trilevel atomic density as a diagonal matrix with occupation numbers on the diagonal
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param ls lsitem structure containing basis info
subroutine trilevel_ATOMS_density(D,ls,REGBASIS)
  use BUILDAOBATCH
implicit none
TYPE(lsitem) :: ls
TYPE(basissetinfo) :: REGbasis
Type(Matrix) :: D
real(realk),pointer  :: occ(:), tmp(:,:)
integer :: nAtoms, nAngmom, norb, ipos, charge,icharge
integer :: itype, ang, i, nbast, kmult,R

 nbast = REGBASIS%nbast
 call mem_alloc(occ,nbast)

 occ = 0E0_realk;

 nAtoms = ls%input%MOLECULE%nAtoms
 R = REGBASIS%Labelindex
 ipos = 1
 do i=1, nAtoms
    IF(ls%input%MOLECULE%ATOM(i)%phantom)CYCLE
    IF(ls%input%MOLECULE%ATOM(i)%pointcharge)CYCLE
    IF(REGBASIS%labelindex .EQ. 0)THEN
       icharge = INT(ls%input%MOLECULE%ATOM(i)%charge) 
       itype = REGBASIS%chargeindex(icharge)
    ELSE
       itype = ls%input%MOLECULE%ATOM(i)%IDtype(R)
    ENDIF
    nAngmom = REGBASIS%ATOMTYPE(itype)%nAngmom
    charge  = REGBASIS%ATOMTYPE(itype)%Charge

    kmult = 1
    do ang = 0,nAngmom-1
       norb = REGBASIS%ATOMTYPE(itype)%SHELL(ang+1)%norb
       call trilevel_set_occ(occ(ipos:nbast),ang,charge)  
       ipos = ipos + (norb*kmult)
       kmult = kmult + 2
    enddo
 enddo

 call mem_alloc(tmp,nbast,nbast)
 tmp = 0E0_realk;

 do i=1,nbast
  tmp(i,i) = occ(i)
 enddo

 call mat_set_from_full(tmp,1E0_realk,D,unres3=.true.)

 call mem_dealloc(occ)
 call mem_dealloc(tmp)
end subroutine trilevel_ATOMS_density

!> \brief convert the Density matrix from the valence basis to full basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D full density matrix
!> \param Dval valence density matrix
!> \param list
!> \param vlist
!> \param len
subroutine trilevel_density_valence2full(D,Dval,list,vlist,len,DENSELEVEL2)
implicit none
Type(Matrix) :: D,Dval
real(realk),pointer  :: tD(:,:), tDval(:,:)
integer :: i,j,k,len
integer :: l1,l2,lvi,lvj,vk
integer, pointer :: vlist(:,:), list(:,:)
logical :: DENSELEVEL2

 call mem_alloc(tD,D%ncol,D%nrow)
 call mem_alloc(tDval,Dval%ncol,Dval%nrow)
 call set_matop_timer_optlevel(2)

 IF(DENSELEVEL2)THEN
  call mat_to_full(Dval,1E0_realk,tDval,matype=mtype_dense) !Level 2 Mat 
 ELSE
  call mat_to_full(Dval,1E0_realk,tDval)                   !Level 2 Mat 
 ENDIF
 tD = 0E0_realk

 do j=1,len
 if (vlist(j,1).gt.vlist(j,2)) cycle
  vk=0
  do l1 = list(j,1),list(j,2)
  do i=1,len
    if (vlist(i,1).gt.vlist(i,2)) cycle
     lvj=vlist(j,1)+vk
     k=0
     do l2 = list(i,1),list(i,2)
      lvi =vlist(i,1)+k

      tD(l2,l1) = tDval(lvi,lvj)
      k=k+1
    !  tD(list(i,1):list(i,2),list(j,1):list(j,2)) = &
    !    & tDval(vlist(i,1):vlist(i,2),vlist(j,1):vlist(j,2))
     enddo
    enddo
     vk=vk+1

  enddo
 enddo
 call set_matop_timer_optlevel(3)
 call mat_set_from_full(tD,1E0_realk,D,'Dmat')  !Level 3 Mat

 call mem_dealloc(tD)
 call mem_dealloc(tDval)
end subroutine trilevel_density_valence2full

!> \brief convert the MO coefficients from the valence basis to full basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param CMO MO coefficients in full basis
!> \param CMO MO coefficients in valence basis
!> \param list
!> \param vlist
!> \param len
subroutine trilevel_cmo_valence2full(CMO,vCMO,list,vlist,len,DENSELEVEL2)
implicit none
Type(Matrix) :: CMO,vCMO
real(realk),pointer :: tCmo(:,:), tvCMO(:,:)
integer :: i,j,k,lv,vk,len
integer, pointer :: vlist(:,:), list(:,:)
logical :: DENSELEVEL2

 call mem_alloc(tCMO,CMO%ncol,CMO%nrow)
 call mem_alloc(tvCMO,vCMO%ncol,vCMO%nrow)

 call set_matop_timer_optlevel(2)
 IF(DENSELEVEL2)THEN
  call mat_to_full(vCMO,1E0_realk,tvCMO,matype=mtype_dense) !Level 2 Mat 
 ELSE
  call mat_to_full(vCMO,1E0_realk,tvCMO)                   !Level 2 Mat 
 ENDIF
 tCMO = 0E0_realk

 do i=1,len
 if (vlist(i,1).gt.vlist(i,2)) cycle
 vk = 0
 do j= list(i,1),list(i,2)
    lv = vlist(i,1) + vk
 do k =1,vCMO%ncol

   
    tCMO(j,k) = tvCMO(lv,k)

    !tCMO(list(i,1):list(i,2),1:vCMO%ncol) = &
   !& tvCMO(vlist(i,1):vlist(i,2),1:vCMO%ncol)
 enddo
 vk=vk+1
 enddo
 enddo


 k = vCMO%ncol + 1;
 do i=1, len-1
 if (vlist(i,1).gt.vlist(i,2)) cycle
   do j=list(i,2)+1 , list(i+1,1)-1
      tCMO(j,k) = 1E0_realk
      k = k + 1
   enddo  
 enddo

 do j=list(len,2)+1, CMO%ncol
      tCMO(j,k) = 1E0_realk
      k = k + 1
 enddo

 call set_matop_timer_optlevel(3)
 call mat_set_from_full(tCMO,1E0_realk,CMO) !Level 3 Mat 

 call mem_dealloc(tCMO)
 call mem_dealloc(tvCMO)
end subroutine trilevel_cmo_valence2full

!> \brief read density from disk
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param 
!> \param 
!> \param 
!> \param 
!> \param 
subroutine trilevel_readdens(D,lupri)
implicit none
Type(Matrix) :: D
integer      :: lu, lupri, IDUM,LDUM
logical      :: OnMaster
data lu /-1/
OnMaster = .TRUE.
call lsopen(lu,'vdens.restart','OLD','UNFORMATTED')
rewind lu
call mat_read_from_disk(lu,D,OnMaster)
call lsclose(lu,'KEEP')

WRITE(LUPRI,*)
WRITE(LUPRI,*) '*** RESTART FROM DENSITY ON DISK - READ FROM vdens.restart  ***'
WRITE(LUPRI,*)
WRITE(*,*)
WRITE(*,*) '*** RESTART FROM DENSITY ON DISK - READ FROM vdens.restart  ***'
WRITE(*,*)

end subroutine trilevel_readdens

end module trilevel_module

!> \brief main wrapper routine to obtain the trilevel_basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt optItem containing info about scf optimization
!> \param ls lsitem structure containing integral,molecule,basis info
SUBROUTINE trilevel_basis(opt,ls)
use basis_typetype
use trilevel_module
use typedeftype, only: lsitem
use io, only: io_free, io_init
use screen_mod, only: screen_free, screen_init
use ecdata
use opttype
use precision
use memory_handling
use Matrix_Operations, only: matrix_type,mtype_dense,mat_select_type
use lstiming
use typedef, only: getNbasis
use GCtransMod, only: write_GCtransformationmatrix
use lsparameters
implicit none
type(optItem), intent(inout) :: opt
type(optItem)             :: gcopt
TYPE(lsitem) :: ls
TYPE(trilevel_atominfo) :: ai
integer                 :: CFG_averaging_sav, matrix_sav,nbast
logical                 :: unres_sav

  IF(ls%input%basis%BINFO(REGBASPARAM)%GCbasis)THEN
     !the atomic calculation
     WRITE(ls%LUPRI,'(A)')' Skipping gcbasis calculation'
     WRITE(ls%LUPRI,'(A)')' This have already been done and input basis transformed.'
     !we have transformed the input basis to GCbasis
     ls%setting%integraltransformGC = .FALSE.              
  ELSE
     opt%optlevel = 1
     ls%optlevel = 1
     call set_matop_timer_optlevel(1)
     
     gcopt = opt
     
     matrix_sav  = matrix_type
     call mat_select_type(mtype_dense,gcopt%lupri)
     !  matrix_type = mtype_dense !this is not allowed due to MPI complications
     
     gcopt%cfg_unres = .false.
     
     !initialise atominfo : find # unique atoms, and map to atoms in molecule
     call trilevel_atominfo_init(ai,ls,gcopt%LUPRI)
     
     !main driver to build the grand canonical basis
     call trilevel_gcbasis(gcopt,ls,ai,gcopt%LUPRI,gcopt%LUERR)
     
     !free atominfo
     call trilevel_atominfo_free(ai)
     
     call mat_select_type(matrix_sav,gcopt%lupri,ls%input%MOLECULE%nbastREG)
     !  matrix_type = matrix_sav !this is not allowed due to MPI complications
     
     !This must be done after the mat_select_type
     IF(ls%setting%integraltransformGC)THEN
        nbast = getNbasis(AORdefault,Contractedinttype,ls%input%MOLECULE,gcopt%LUPRI)
        call write_GCtransformationmatrix(nbast,ls%setting,gcopt%lupri)
     ENDIF

     !CFG_averaging = CFG_averaging_sav
     !we have written screening matrices which should not 
     !be used in the rest of the calculations so we reset the IO  
     call io_free(ls%setting%IO)
     call io_init(ls%setting%IO)
     call set_matop_timer_optlevel(3)
  ENDIF
END SUBROUTINE trilevel_basis

!> \brief the atoms start guess routine. The gcao basis have been created by a call to trilevel_basis and now a start D is made
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param decomp Contains matrices from OAO decomposition of overlap matrix
!> \param D density matrix 
!> \param H1 one electron contribution to the fock matrix
!> \param ls lsitem structure containing all info about integrals,basis,..
SUBROUTINE atoms_start(config,D,H1,S,ls,ndmatalloc)
use decompMod
use basis_typetype
use configurationType
use trilevel_module
use typedeftype, only: lssetting, lsitem
use matrix_module
use precision
use memory_handling
use matrix_operations
use Matrix_Operations_aux!, only: mat_density_from_orbs
use matrix_util
use lsdalton_fock_module
use dal_interface
use initial_guess
implicit none
type(ConfigItem),intent(in) :: config
TYPE(lsitem),intent(inout) :: ls
Type(Matrix),target        :: H1
Type(Matrix),intent(inout) :: D(ndmatalloc),S
integer,intent(in)         :: ndmatalloc
!
Type(Matrix) :: Dval,F(1),Cmo,Dpure,Dpure2
TYPE(trilevel_atominfo) :: ai
real(realk)         :: E(1),trace
real(realk),pointer :: eival(:) 
integer      :: nbast,Nelectrons,sz,ndmat
logical      :: dalink,DiagFmat,McWeeny,purify_failed,CS00,admm_exchange
real(realk),parameter :: THRNEL=1E-3_realk
real(realk),external :: HOMO_energy
ndmat = 1
!a diagonal matrix with occupation numbers on the diagonal
call trilevel_ATOMS_density(D(1),ls,ls%input%basis%binfo(RegBasParam))

nbast = ls%input%BASIS%BINFO(REGBASPARAM)%nbast

!The matrix D is not idempotent. We have 2 options: 
!1. Do McWeeny purification
!2. Perform a Fock matrix build followed by a diagonalization
McWeeny = config%opt%MWPURIFYATOMSTART  !Option 1 
DiagFmat = .FALSE.  !Option 2 
IF(McWeeny)THEN
   write(config%lupri,'(A,I3)')'Perform McWeeny purification on the non idempotent Atoms Density',1
   Nelectrons = 2*config%decomp%nocc + config%decomp%nactive
   call mat_init(Dpure,D(1)%nrow,D(1)%ncol)
   call mat_assign(Dpure,D(1))
   call Robust_McWeeney_purify(S,Dpure,purify_failed,Nelectrons)     
   DiagFmat = purify_failed !revert to diagonalization if fail
   if (.not.purify_failed) then
      trace = mat_dotproduct(Dpure,S)
      if ( abs(2.0E0_realk*trace-Nelectrons) .le. THRNEL*Nelectrons)then
         call mat_assign(D(1),Dpure)
      else
         write(config%lupri,*)'Warning: Descrepancy between 2.0*Tr(DS) after purification and Nelectrons '
         write(config%lupri,*)'2.0*Tr(DS) = ',2.0E0_realk*trace,' Nelectrons = ', Nelectrons
         write(config%lupri,*)'We fall back to a diagonalization'
         DiagFmat = .TRUE.
      endif
   else
      write(config%lupri,*)'McWeeny purification of atomic density failed. We fall back to a diagonalization'
   endif
   call mat_free(Dpure)
!   call mat_free(Dpure2)
ELSE
   DiagFmat = .TRUE. 
ENDIF

CALL mat_init(F(1),nbast,nbast)

IF(DiagFmat)THEN
   !We cannot use DaLink in the 0'th iteration - this gives a diagonal Fock matrix
   ! => bad starting guess. If DaLink is requested, turn it off and then back on after
   ! the 0'th iteration. /Stinne, Thomas, Brano 19/11-2009
   dalink = ls%setting%scheme%DALINK
   ls%setting%scheme%DALINK = .FALSE.
   CS00 = ls%setting%scheme%DFT%CS00
   ls%setting%scheme%DFT%CS00 = .FALSE.
   !turn off ADMM exchange due missing idempotency 
   ADMM_EXCHANGE = ls%setting%scheme%ADMM_EXCHANGE
   IF (ls%setting%scheme%ADMM1) ls%setting%scheme%ADMM_EXCHANGE = .FALSE.
   ls%setting%scheme%DFT%CS00eHOMO = config%diag%eHOMO
   ls%setting%scheme%DFT%DFTELS = 1E0_realk !the density is not idempotent so it would giv e a wrong number of electrons
   ! Iteration 0 : The density matrix is not idempotent; a diagonalization gives a proper 
   ! idempotent density matrix  
   call di_get_fock_LSDALTON(D,H1,F,ndmat,E,config%decomp%lupri,config%decomp%luerr,ls)
   write(*,*) ' Iteration 0 energy:', E(1)
   write(config%decomp%lupri,*) ' Iteration 0 energy:', E(1)
   
   if (config%decomp%cfg_unres) then
      call mem_alloc(eival,2*nbast)
   else
      call mem_alloc(eival,nbast)
   endif
   
   CALL mat_init(Cmo,nbast,nbast)
   call mat_diag_f(F(1),config%decomp%S,eival,Cmo)
   !Commentet out by Thomas Kjaergaard - no testcase - code not testet ...
   !!Asymetrizing starting guess if .ASYM is in input
   !! 21.04.2010 C. Nygaard
   !!Only works if HOMO and LUMO are of different symmetry
   !if (config%decomp%cfg_unres .and. config%opt%cfg_asym) then
   !   call asymmetrize_starting_guess (Cmo, config%decomp)
   !endif
   
   call mat_density_from_orbs(Cmo,D(1),config%decomp%nocc,config%decomp%nocca,config%decomp%noccb)
   
   if (config%decomp%cfg_unres) then
      sz = 2*CMO%nrow
   else
      sz = CMO%nrow
   endif
   CALL mat_free(Cmo)

   ls%setting%scheme%DFT%CS00eHOMO = HOMO_energy(config%decomp%cfg_unres,&
        & config%decomp%nocc,config%decomp%nocca,config%decomp%noccb,eival,sz)
   
   call mem_dealloc(eival)
   
   !Turn DaLink back on, if requested:
   ls%setting%scheme%DALINK = dalink
   ls%setting%scheme%DFT%CS00 = cs00
   ls%setting%scheme%ADMM_EXCHANGE = ADMM_EXCHANGE
   
   ls%setting%scheme%DFT%DFTELS = ls%input%dalton%DFT%DFTELS      
ENDIF
CALL mat_free(F(1))

END SUBROUTINE ATOMS_START

!> \brief Main Driver routine for the trilevel start guess
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param ls structure containing info about molecule and integral evaluation
!> \param config structure containing basicly all info
SUBROUTINE trilevel_start(D,ls,config)
use decompMod
use configurationType
use dal_interface
use trilevel_module
use typedeftype, only: lssetting, lsitem
use screen_mod, only: screen_free,screen_init
use precision
use memory_handling
use files
use ecdata
use lsdalton_fock_module
use integralinterfaceMod
use typedef
use dal_interface
use ks_settings, only: ks_init_incremental_fock, ks_free_incremental_fock
use direct_dens_util
use daltoninfo
use matrix_module
use matrix_operations
use Matrix_Operations_aux!, only: mat_density_from_orbs
use matrix_util
use lsparameters
use IntegralInterfaceMOD
use diagonalization
use scfloop_module
use molecule_module
use optimlocMOD, only: optimloc
implicit none
TYPE(lsitem),target :: ls
TYPE(trilevel_atominfo) :: ai
Type(Matrix) :: Dval(1),F(1),D(1),Cmo, fCmo,Dpure,Dpure2
Type(Matrix),target :: H1,S
Type(lsint_fock_data_type) :: lsint_fock_data_sav
real(realk)         :: E(1), mx
real(realk)         :: maxelm_save, maxstep_save,thr_save,trace
real(realk),pointer :: eival(:) 
integer      :: nbast, len, nocc, Nelectrons
integer, pointer :: vlist(:,:), list(:,:)
integer :: lun, idum, ldum,iAO,sz,ndmat,matrix_sav
integer(8) :: fperm, vperm
logical :: restart_from_dens, no_rhdiis, dalink, vdens_exists,admm_exchange
logical :: OnMaster,DiagFmat,McWeeny,purify_failed,CS00,integraltransformGC
character(len=11) :: string11
type(ConfigItem) :: config
real(realk),parameter :: THRNEL=1E-3_realk
real(realk),external :: HOMO_energy
type(LowAccuracyStartType)  :: LAStype

  ndmat = 1
  OnMaster = .TRUE.
  !initialise the IO so that old screening matrices are not used 
  call io_free(ls%setting%IO)
  call io_init(ls%setting%IO)
  call screen_free()
  call screen_init()

  !build a minls used in the valens calculation
!  call trilevel_alloc_sync_ls(minls,ls)

  call set_matop_timer_optlevel(-1)
  !build a trilevel_atominfo structure containing info for each distinct atom
  call trilevel_atominfo_init(ai,ls,config%LUPRI)

  !Level 2 integral eval must be done in the GCbasis 
  integraltransformGC = ls%setting%integraltransformGC
  ls%setting%integraltransformGC = .FALSE.

  !turn off ADMM exchange ass it does not make sense 
  !Warning it turns out this is redundant as it is turned off in 
  !SolverUtilities/dalton_interface.F90 if config%opt%optlevel .EQ. 2
  ADMM_EXCHANGE = ls%setting%scheme%ADMM_EXCHANGE
  ls%setting%scheme%ADMM_EXCHANGE = .FALSE.

  !loops over all distinct atoms and build valensbasis
  call trilevel_full2valence(ls,ai,config%LUPRI)

  IF(.NOT.ASSOCIATED(ls%SETTING%BASIS(1)%p,ls%input%basis))THEN
     DO iAO=1,ls%SETTING%nAO
        IF(ls%SETTING%BASIS(iAO)%p%BINFO(VALBASPARAM)%nAtomtypes.NE.0)THEN
           call lsquit('this should not happen setting valence basis non zero',-1)
           call free_basissetinfo(ls%SETTING%BASIS(iAO)%p%BINFO(VALBASPARAM))
        ENDIF
        call copy_basissetinfo(ls%input%basis%binfo(valbasparam),ls%SETTING%BASIS(iAO)%p%BINFO(VALBASPARAM))
        call determine_nbast(ls%setting%MOLECULE(iAO)%p,ls%SETTING%BASIS(iAO)%p%BINFO(VALBASPARAM),&
             &ls%setting%scheme%DoSpherical,ls%setting%scheme%uncont)
     ENDDO
  ENDIF
  !in the valence basis is built as a GC basis
  !so we do not need to transform when doing integrals

  !set setting basis to the valens basis
  call set_default_AOs(AOVAL,AOdfAux)

  write(config%lupri,'(1X,A)') 'Level 2 molecular calculation'
  write(   * ,'(1X,A)') 'Level 2 molecular calculation'
  write(config%lupri,'(A)') '============================='
  Write(ls%lupri,'(A)')' '
  Write(ls%lupri,'(1X,A)')'We now perform the 2. Level of the TRILEVEL starting guess '
  Write(ls%lupri,'(1X,A)')'see PCCP 11, 5805-5813 (2009) as well as JCTC 5, 1027 (2009)'
  Write(ls%lupri,'(1X,A)')'The 2. level is a full molecular calculation done in a minimal basis'
  Write(ls%lupri,'(A)')' '
  write(config%lupri,'(A)') '  The 2. Level Basis'
  CALL PRINT_LEVEL2BASIS(config%lupri,&
       & ls%input%MOLECULE,ls%input%basis%BINFO(VALBASPARAM))
!  call print_basissetinfo(config%LUPRI,ls%input%basis%VALENCE)
  ls%optlevel = 2

  nbast = ls%input%basis%BINFO(VALBASPARAM)%nbast !nbast for valence basis
  call set_matop_timer_optlevel(2)

  IF(config%opt%DENSELEVEL2)THEN
     matrix_sav  = matrix_type
     call mat_select_tmp_type(mtype_dense,config%lupri,nbast)
     Write(ls%lupri,'(A)')'Chosing to do Level 2 in Dense matrix type'
  ENDIF

  call mat_init(Dval(1),nbast,nbast) !Level 2 Mat

  !select starting from atoms or from file
  INQUIRE(file='vdens.restart',EXIST=vdens_exists) 
  restart_from_dens = config%diag%cfg_restart.AND.vdens_exists
  ! Never restart from file for geometry optimization
  ! because the level 2 density is not idempotent.
  if(config%optinfo%optimize .OR. config%dynamics%do_dynamics) then
     restart_from_dens=.false.
  end if

  if (restart_from_dens) then
     !.RESTART specified in input and if a density 
     !have been written to file
     call trilevel_readdens(Dval(1),config%lupri)
  else
     !default option
     call trilevel_ATOMS_density(Dval(1),ls,ls%input%BASIS%BINFO(VALBASPARAM))
  endif

  CALL mat_init(F(1),nbast,nbast)
  CALL mat_init(H1,nbast,nbast)
  CALL mat_init(S,nbast,nbast)
  CALL mat_init(Cmo,nbast,nbast)

  !Get Level2 screening matrices
  call II_precalc_ScreenMat(config%decomp%lupri,config%decomp%luerr,ls%SETTING)
  !Get overlap and H1
  CALL II_get_overlap(config%decomp%lupri,config%decomp%luerr,ls%setting,S)
  CALL II_get_h1(config%decomp%lupri,config%decomp%luerr,ls%setting,H1)

  config%decomp%S => S

  !save L3 inputs
  lsint_fock_data_sav = lsint_fock_data
  lsint_fock_data%ls  => ls
  lsint_fock_data%H1  => H1
  lsint_fock_data%lupri = config%lupri
  lsint_fock_data%luerr = config%decomp%luerr
  !maxstep_save = cfg_max_step    
  !maxelm_save = cfg_max_element
  !thr_save = cfg_convergence_threshold !We don't need to converge hard on level 2
  config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold*config%opt%cfg_level2_convfactor

  !The matrix D is not idempotent. We have 2 options: 
  !1. Do McWeeny purification
  !2. Perform a Fock matrix build followed by a diagonalization
  McWeeny = config%opt%MWPURIFYATOMSTART  !Option 1 
  DiagFmat = .FALSE.  !option 2
  IF(McWeeny)THEN
     write(config%lupri,'(A,I3)')'Perform McWeeny purification on the non idempotent Atoms Density',1
     Nelectrons = 2*config%decomp%nocc + config%decomp%nactive
     call mat_init(Dpure,Dval(1)%nrow,Dval(1)%ncol)
     call mat_assign(Dpure,Dval(1))
     call Robust_McWeeney_purify(S,Dpure,purify_failed,Nelectrons)
     DiagFmat = purify_failed !revert to diagonalization if fail
     if (.not.purify_failed) then
        trace = mat_dotproduct(Dpure,S)
        if ( abs(2.0E0_realk*trace) .le. THRNEL*Nelectrons)then
           call mat_assign(Dval(1),Dpure)
        else
           write(config%lupri,*)'Warning: Descrepancy between 2.0*Tr(DS) after purification and Nelectrons '
           write(config%lupri,*)'2.0*Tr(DS) = ',2.0E0_realk*trace,' Nelectrons = ', Nelectrons
           write(config%lupri,*)'We fall back to a diagonalization'
           DiagFmat = .TRUE.
        endif
     else
        write(config%lupri,*)'McWeeny purification of atomic density failed. We fall back to a diagonalization'
     endif
     call mat_free(Dpure)
  ELSE
     DiagFmat = .TRUE. 
  ENDIF
  
  IF(DiagFmat)THEN
     !We cannot use DaLink in the 0'th iteration - this gives a diagonal Fock matrix
     ! => bad starting guess. If DaLink is requested, turn it off and then back on after
     ! the 0'th iteration. /Stinne, Thomas, Brano 19/11-2009
     dalink = ls%setting%scheme%DALINK
     ls%setting%scheme%DALINK = .FALSE.
     CS00 = ls%setting%scheme%DFT%CS00
     ls%setting%scheme%DFT%CS00 = .FALSE.
     
     ls%setting%scheme%DFT%CS00eHOMO = config%diag%eHOMO
     ! Iteration 0
     ls%setting%scheme%DFT%DFTELS = 1E0_realk !the density is not idempotent so it would give a wrong number of electrons
     if (.not.restart_from_dens) then
        call di_get_fock_LSDALTON(Dval,H1,F,ndmat,E,config%lupri,config%decomp%luerr,ls)
        write(*,*) ' Iteration 0 energy:', E(1)
        write(config%lupri,*) ' Iteration 0 energy:', E(1)
        call mem_alloc(eival,nbast)
        call mat_diag_f(F(1),S,eival,Cmo)
        call mat_density_from_orbs(Cmo,Dval(1),config%decomp%nocc,config%decomp%nocca,config%decomp%noccb)
        if (config%decomp%cfg_unres) then
           sz = 2*CMO%nrow
        else
           sz = CMO%nrow
        endif
        ls%setting%scheme%DFT%CS00eHOMO = HOMO_energy(config%decomp%cfg_unres,&
             & config%decomp%nocc,config%decomp%nocca,config%decomp%noccb,eival,sz)
        
        call mem_dealloc(eival)
     endif
     ls%setting%scheme%DFT%DFTELS = ls%input%dalton%DFT%DFTELS
     !Turn DaLink back on, if requested:
     ls%setting%scheme%DALINK = dalink
     ls%setting%scheme%DFT%CS00 = CS00
  ENDIF

  !initialize incremental scheme
  if (config%opt%cfg_incremental) call ks_init_incremental_fock(nbast)
     
  !Does optimization method require overlap decomposition?
  no_rhdiis =(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
       & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
       & config%decomp%cfg_check_converged_solution .or. &
       & config%decomp%cfg_rsp_nexcit > 0 ) 
  
  !Do decomposition
  if (no_rhdiis .or. config%decomp%cfg_lcv) then   
     call decomp_init(nbast,config%decomp)
     config%decomp%lcv_basis = .false.
     call decomposition(config%decomp)
  endif
  
#if 0
  !Prepare for L2 calculation in localized CMO basis.
  !Do localization of initial orbitals and do decomposition
  !to get local orthonormal CMO basis
  if (config%decomp%cfg_lcv.and.(.not.restart_from_dens).and. no_rhdiis) then
     nocc = config%decomp%nocc
     call leastchange_lcv(config%decomp,Cmo,nocc,ls)
     
     call leastchangeOrbspreadStandalone(mx,ls,Cmo,config%decomp%lupri,config%decomp%luerr)
     write(*,*) 'Orbspread standalone: ', mx
     
     config%decomp%lcv_basis = .true.
     call mat_init(config%decomp%lcv_CMO,nbast,nbast)
     config%decomp%decompMatInit_lcv_CMO = .TRUE.
     call mat_assign(config%decomp%lcv_CMO,Cmo)
     
     call save_decomposition(config%decomp)
     call decomposition(config%decomp)
  endif
#endif
  
  if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then !FIXME: put this somewhere else!
     call mat_init(config%av%Fprev,nbast,nbast)
     call mat_init(config%av%Dprev,nbast,nbast)
  endif

  ! L2 minimization
  config%opt%optlevel = 2

  IF(config%integral%LOW_ACCURACY_START)THEN
     call set_Low_accuracy_start_settings(config%lupri,ls,config,LAStype)
     call scfloop(H1,F,Dval,S,E,ls,config)
     call revert_Low_accuracy_start_settings(config%lupri,ls,config,LAStype)
  ENDIF
  call scfloop(H1,F,Dval,S,E,ls,config)
  write(   * ,'(1X,A)') 'Level 2 minimization done! Starting Level 3 minimization...'
  write(config%lupri,'(1X,A)') 'Level 2 minimization done! Starting Level 3 minimization...'
  Write(ls%lupri,'(1X,A)')'Level 3 is the final level of the TRILEVEL algorithm '
  Write(ls%lupri,'(1X,A)')'and is the standard calculation using the TRILEVEL starting guess'
  Write(ls%lupri,'(1X,A)')'obtained from Level 2.'

  write(config%lupri,*)
  config%opt%optlevel = 3
  ls%optlevel = 3
   if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then
      call mat_free(config%av%Fprev)
      call mat_free(config%av%Dprev)
   endif

  !Dump level 2 density matrix to disk
  !densL2_lun = -1
  !call lsopen(densL2_lun,'densL2.restart','unknown','UNFORMATTED')
  !call mat_write_to_disk(densL2_lun,Dval)
  !call lsclose(densL2_lun,'KEEP')

  !release data related to incremental fock
  if (config%opt%cfg_incremental) call ks_free_incremental_fock
  !restore L3 inputs
  lsint_fock_data = lsint_fock_data_sav
  ls%setting%integraltransformGC = integraltransformGC
  config%solver%set_max_step    = config%solver%cfg_max_step
  config%solver%set_max_element = config%solver%cfg_max_element
  config%solver%set_local       = .false.
  config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold
  if (config%av%cfg_averaging == config%av%cfg_avg_van_lenthe) then
     config%av%vanlentheCounter = 0
     config%diag%CFG_lshift = diag_lshift_vanlenthe
     config%av%usequeue = .false.
   endif
  !cfg_max_step = maxstep_save
  !cfg_max_element = maxelm_save
  !cfg_dd_local = .false.
  !cfg_convergence_threshold = thr_save

  ! get Valence Cmo
  call mem_alloc(eival,nbast)
  call mat_diag_f(F(1),S,eival,Cmo)
  call mem_dealloc(eival)


  CALL mat_free(H1)
  CALL mat_free(F(1))
  CALL mat_free(S)

#if 0
  !restore overlap decomposition
  if (config%decomp%cfg_lcv.and.(.not.restart_from_dens).and. no_rhdiis) &
  & call restore_decomposition(config%decomp)
#endif

  ! localization of Valence CMO, by least change algorithm
  if (config%decomp%cfg_lcv) then
     nocc = config%decomp%nocc
     ls%setting%integraltransformGC = .false.
     call leastchange_lcv(config%decomp,Cmo,nocc,ls)
     if (config%decomp%cfg_mlo .and. (.not. config%davidOrbLoc%PFM_input%TESTCASE) .and. (.not. config%davidOrbLoc%NOL2OPT)) then
        write(ls%lupri,'(a)')'  %LOC%   '
        write(ls%lupri,'(a)')'  %LOC%   ********************************************'
        write(ls%lupri,'(a)')'  %LOC%   *       LEVEL 2 ORBITAL LOCALIZATION       *'
        write(ls%lupri,'(a)')'  %LOC%   ********************************************'
        write(ls%lupri,'(a)')'  %LOC%   '
        call optimloc(Cmo,config%decomp%nocc,config%decomp%cfg_mlo_m,ls,config%davidOrbLoc)
     endif
     ls%setting%integraltransformGC = integraltransformGC
  endif

  ! shutdown of decomposition and dd
  if (no_rhdiis .or. config%decomp%cfg_lcv) then
     call decomp_shutdown(config%decomp)
     !call dd_shutdown(config%decomp%cfg_unres)
  endif

  ! create integer list for conversion to full basis
  call typedef_setlist_valence2full(list,vlist,len,ls,ls%input%BASIS%BINFO(VALBASPARAM))
  !Due to the construction of Vbasis from regular basis we need to free some space
  !that is not used in the valence basis - But we do not actually free the basis
  IF(.NOT.ASSOCIATED(ls%SETTING%BASIS(1)%p,ls%input%basis))THEN
     DO iAO=1,ls%SETTING%nAO
        call freeVbasis(ls%SETTING%BASIS(iAO)%p)
     ENDDO
  ENDIF
  call freeVbasis(ls%input%BASIS)

  ! initialize decomp%lcv_CMO
  nbast = ls%input%BASIS%BINFO(REGBASPARAM)%nbast !full basis nbast 
  !set setting basis to the full basis
  call set_default_AOs(AORegular,AOdfAux)

  !This is a tricky part of the code where we go from CMO and Dval(1) in the 
  !Valence basis and transform to full level 3 basis. So we have both 
  !Level 2 and Level 3 matrices in Memory. TK

  IF(config%opt%DENSELEVEL2)THEN
     !This restores the proper value
     call mat_select_tmp_type(matrix_sav,config%lupri,nbast)
     call string11_mtype(matrix_sav,string11)
     Write(ls%lupri,'(A,A)')'Done Level 2 so revert to chosen matrix type=',string11
  ENDIF
  call set_matop_timer_optlevel(3)

  call mat_init(config%decomp%lcv_CMO,nbast,nbast) !Level 3 Mat
  config%decomp%decompMatInit_lcv_CMO = .TRUE.

  ! convert valence CMO to full cmo  Level 2 Mat (Cmo) => Level 3 Mat (lcv_Cmo)
  call trilevel_cmo_valence2full(config%decomp%lcv_Cmo,Cmo,list,vlist,len,&
       & config%opt%DENSELEVEL2)
  config%decomp%lcv_basis = .true.
  call set_matop_timer_optlevel(2)

  IF(config%opt%DENSELEVEL2)THEN
     call mat_free(Cmo,matype=mtype_dense)  !Level 2 Mat 
  ELSE
     call mat_free(Cmo)                    !Level 2 Mat 
  ENDIF
  ! convert valence Dval to full D   Level 2 Mat (Dval) => Level 3 Mat (D)
  call trilevel_density_valence2full(D(1),Dval(1),list,vlist,len,config%opt%DENSELEVEL2)
  call set_matop_timer_optlevel(2)
  IF(config%opt%DENSELEVEL2)THEN
     CALL mat_free(Dval(1),matype=mtype_dense) !Level 2 Mat 
  ELSE
     call mat_free(Dval(1))                   !Level 2 Mat 
  ENDIF
  call set_matop_timer_optlevel(3)

  !save lcv_CMO basis
  lun = -1
  call lsopen(lun,'lcv_basis.restart','unknown','UNFORMATTED')
  call mat_write_to_disk(lun,config%decomp%lcv_Cmo,OnMaster)
  call mat_write_info_to_disk(lun,config%decomp%cfg_lcv)
  call lsclose(lun,'KEEP')

  call mem_dealloc(list)
  call mem_dealloc(vlist)
  call set_matop_timer_optlevel(3)
  !we restore the default settings, which sets the GRDONE=0 
  !so that the dft grid is calculated with the full basis 
  !on the next kohn-sham matrix build
!  call typedef_set_default_setting(ls%setting,ls%input)
  CALL typedef_setIntegralSchemeFromInput(LS%INPUT%DALTON,LS%SETTING%SCHEME)

  IF(.NOT.config%decomp%cfg_gcbasis)then
     IF(.NOT.ls%input%basis%BINFO(REGBASPARAM)%DunningsBasis)THEN
        WRITE(config%lupri,*)'Warning: You use trilevel without a GCbasis. This is not'
        WRITE(config%lupri,*)'         recommended and you may be doing something wrong'
        print*,'Warning: You use trilevel without a GCbasis. This is not'
        print*,'         recommended and you may be doing something wrong'
     ENDIF
     !The basis is in AO and we only need AO
     ls%setting%integraltransformGC = .FALSE.              
  ELSE
     IF(ls%input%DALTON%NOGCINTEGRALTRANSFORM)THEN
        !we transform the input basis to GCbasis
        ls%setting%integraltransformGC = .FALSE.              
     ELSE
        IF(ls%input%basis%BINFO(REGBASPARAM)%Gcont)THEN
           !we transform the input basis to GCbasis
           ls%setting%integraltransformGC = .FALSE.     
        ELSE
           !Segmented contracted Case
           !so we transform back and forth.
           ls%setting%integraltransformGC = .TRUE.
           IF(.NOT.ls%input%basis%wbasis(gctbasparam))THEN
              call lsquit('GCtrans not call mem_allocd in trilevel_start',-1)
           ENDIF
        ENDIF
     ENDIF
  ENDIF
  IF(ls%setting%integraltransformGC.NEQV.integraltransformGC)&
       & CALL LSQUIT('ERROR in integraltransformGC',-1)  
  !This should be obsolete 
  call determine_nbast(ls%input%MOLECULE,ls%input%BASIS%BINFO(REGBASPARAM),&
      &ls%setting%scheme%DoSpherical,ls%setting%scheme%uncont)

  call io_free(ls%setting%IO)
  call io_init(ls%setting%IO)! we change basis, can nolonger use the mat on disk
  call screen_free()
  call screen_init()! we change basis, can nolonger use the gabs in memory
  !Get Level3 screening matrices (they have been calculated once before)
  call II_precalc_ScreenMat(config%decomp%lupri,config%decomp%luerr,ls%SETTING)

  call leastchangeOrbspreadStandalone(mx,ls,config%decomp%lcv_Cmo,config%decomp%lupri,config%decomp%luerr)
  write(*,*) 'Orbspread standalone full CMO: ', mx

  call trilevel_atominfo_free(ai)
  ls%setting%scheme%ADMM_EXCHANGE = ADMM_EXCHANGE

END SUBROUTINE trilevel_start

