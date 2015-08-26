!> @file

!> Subroutines related to the construction of the CCSD gradient
!> \author Dmytro Bykov
!> \date   May 2014

module ccsd_gradient_module

  ! ================= Declarations / Specifications =================
  use fundamental
  use precision
  use dec_typedef_module
  use lsparameters
  use tensor_interface_module
  use matrix_operations

  use dec_typedef_module
  use array2_simple_operations
  use array4_simple_operations
  use integralinterfaceMod
  use ccintegrals
  use mp2_gradient_module
  use ccsd_module
#ifdef VAR_MPI
  use infpar_module
#endif
  ! ========================= Implementation ========================
contains

  !> ----------------------------------------------------------------
  !> \brief Calculate contributions to gradient  
  !>        from the single fragment
  !> ON INPUT
  !>            MyFragment      - DEC info for a fragment
  !>            m2occ, m2virt   - Lagrange multipliers
  !>            t2occ, t2virt   - amplitudes
  !> ON OUTPUT
  !>            grad            - gradint structure 
  !> \author Dmytro Bykov
  !> \date   September 2014
  !> ----------------------------------------------------------------
  subroutine single_calculate_CCSDgradient_driver&
                          &(MyFragment,t2occ,t2virt,m2occ,m2virt,m1,&
                                    &VOOO,VOVV,VOVOocc,VOVOvirt,grad)
  
    implicit none

    !INPUT declarations
    !******************
    type(decfrag), intent(inout)  :: MyFragment
    type(tensor),intent(inout)    :: m2occ
    type(tensor),intent(inout)    :: m2virt
    type(tensor),intent(inout)    :: t2occ
    type(tensor),intent(inout)    :: t2virt
    type(tensor),intent(inout)    :: m1 
    type(tensor),intent(in)       :: VOOO
    type(tensor),intent(in)       :: VOVV
    type(tensor),intent(inout)    :: VOVOocc   
    type(tensor),intent(in)       :: VOVOvirt

    !OUTPUT
    !******
    !DEBUG try to use MP2 gradient structure
    type(mp2grad),intent(inout) :: grad

    !LOCAL INTERMEDIATES
    !*******************

    !Service variables
    integer           :: i,j,a
    integer           :: nb,noEOS,nvEOS,noAOS,nvAOS
    type(array2)      :: temp1

    nb    = MyFragment%nbasis
    noEOS = MyFragment%noccEOS
    noAOS = MyFragment%noccAOS
    nvEOS = MyFragment%nvirtEOS
    nvAOS = MyFragment%nvirtAOS

    !DEBUG Init MP2 gradient structure
    call init_mp2grad(MyFragment,grad)
 
    !DEBUG
    print *, "Debug print m2_aibj: "
    do a=1,m2occ%nelms 
      if (ABS(m2occ%elm1(a))>0.00000001)then 
        write (*,'(i5,F12.8)') a, m2occ%elm1(a)
      endif
    end do 

    !DEBUG
    print *, "Debug print t2_aibj: "
    do a=1,t2occ%nelms
      if (ABS(t2occ%elm1(a))>0.00000001)then
        write (*,'(i5,F12.8)') a, t2occ%elm1(a)
      endif
    end do

    !calculate the density
    !**********************
    !the diagonal blocks dab and dij of the one-electron density 
    call single_dab_dij_formation(MyFragment,grad,t2occ,t2virt,m2occ,m2virt)

    !the off-diagonal blocks dia and dai
    call single_dai_dia_formation(MyFragment,grad,t2occ,t2virt,m1)

    !clean up if needed


  end subroutine single_calculate_CCSDgradient_driver

  !> ----------------------------------------------------------------
  !> \brief For a given fragment, calculate Dab and Dij parts of the 
  !>        unrelaxed one-electron density 
  !>
  !> ON INPUT
  !>            MyFragment      - DEC info for a fragment
  !>            m2occ, m2virt   - Lagrange multipliers
  !>            t2occ, t2virt   - amplitudes
  !> ON OUTPUT
  !>            dab,dij         - Dab and Dij parts of the unrel. dens.
  !>            grad            - grad is updated with Dab and Dij
  !> \author Dmytro Bykov
  !> \date   August 2014
  !> ----------------------------------------------------------------
  subroutine single_dab_dij_formation(MyFragment,grad,t2occ,t2virt,m2occ,m2virt)
   
    implicit none

    !INPUT declarations
    !******************
    type(decfrag), intent(inout) :: MyFragment
    type(tensor),intent(inout)    :: m2occ
    type(tensor),intent(inout)    :: m2virt
    type(tensor),intent(inout)    :: t2occ
    type(tensor),intent(inout)    :: t2virt
    !DEBUG use MP2 gradient structure for now
    type(mp2grad),intent(inout) :: grad

    !OUTPUT
    !******
    type(array2)   :: dab 
    type(array2)   :: dij

    !LOCAL INTERMEDIATES
    !*******************

    !Service variables
    type(array2)   :: Xij
    type(array4)   :: temp1, temp2
    integer        :: i,j,nb,noEOS,nvEOS,noAOS,nvAOS
    !DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    integer        :: a,b,aFull,bFull,iFull,jFull
    type(tensor)    :: m2tmp
    !DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    logical        :: something_wrong
    real(realk)    :: tcpu,twall

    !nb    = MyFragment%nbasis
    noEOS = MyFragment%noccEOS
    noAOS = MyFragment%noccAOS
    nvEOS = MyFragment%nvirtEOS
    nvAOS = MyFragment%nvirtAOS

    !start the stopwatch
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    !Initialize SD density
    dab  = array2_init([nvAOS,nvAOS])
    dij  = array2_init([noAOS,noAOS])

    !Zero matrices in dens structure
    grad%dens%X = 0e0_realk
    grad%dens%Y = 0e0_realk

    ! Sanity check 1: Density input structure
    !****************************************
    something_wrong=.false.
    if(noAOS/=grad%dens%nocc)   something_wrong=.true.
    if(nvAOS/=grad%dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'dens%nocc  ', grad%dens%nocc
       write(DECinfo%output,*) 'dens%nvirt', grad%dens%nvirt
       write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noAOS
       write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvAOS
       call lsquit('single_dab_dij_formation: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

    ! Sanity check 2: Array dimensions
    !*********************************
    if(nvAOS /= t2occ%dims(3)) something_wrong=.true.
    if(noEOS /= t2occ%dims(4)) something_wrong=.true.
    if(nvEOS /= t2virt%dims(3)) something_wrong=.true.
    if(noAOS /= t2virt%dims(4)) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'EOS: nocc, nvirt', noEOS, nvEOS
       write(DECinfo%output,*) 'AOS: nocc, nvirt', noAOS, nvAOS
       write(DECinfo%output,*) 't2occ dims     ', t2occ%dims
       write(DECinfo%output,*) 't2virt dims    ', t2virt%dims
       call lsquit('single_dab_dij_formation: &
            & Something wrong with amplitude and Theta dimensions!',DECinfo%output)
    end if

    !******************************
    !Do the calculation           *
    !******************************
    ! virt-virt block
    !******************
    !using the occupied partitioning scheme
    !Dab = Sum_ckl (t^cb_kl * z^ca_kl) = Y_ba

    !containers for amps/mults voov shape
    temp1 = array4_init([nvAOS,noEOS,noEOS,nvAOS])
    temp2 = array4_init([nvAOS,noEOS,noEOS,nvAOS])

    !sort amps(ckbl) -> cklb
    call array_reorder_4d(1.0E0_realk,m2occ%elm4,nvAOS,noEOS,nvAOS,noEOS,[1,2,4,3],0.0E0_realk,temp1%val)
    call array_reorder_4d(1.0E0_realk,t2occ%elm4,nvAOS,noEOS,nvAOS,noEOS,[1,2,4,3],0.0E0_realk,temp2%val)

    !form virt-virt dens. block 
    call array4_contract3(temp2,temp1,dab)
    
    !fill in the grad structure 
    grad%dens%Y(1:nvAOS,1:nvAOS) = dab%val(1:nvAOS,1:nvAOS)
 
    !clean up!
    call array4_free(temp1)
    call array4_free(temp2)
    call array2_free(dab)

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "DEC SINGLES: virt-virt block of unrel. one-el. dens."
       call print_matrix_real(grad%dens%Y,nvAOS,nvAOS)
       print *, "";    endif

    !time
    call LSTIMER('Y MATRIX',tcpu,twall,DECinfo%output)

    ! occ-occ block
    !****************
    !using the virtual partitioning scheme
    !Dij = 2*delta_ij - Sum_cdk (t^cd_ki * z^cd_kj) = 2*delta_ij - X_ij

    !containers for amps/mults vovo shape
    temp1 = array4_init([nvEOS,noAOS,nvEOS,noAOS])
    temp2 = array4_init([nvEOS,noAOS,nvEOS,noAOS])
    call array_reorder_4d(1.0E0_realk,m2virt%elm4,nvEOS,noAOS,nvEOS,noAOS,[1,2,3,4],0.0E0_realk,temp1%val)
    call array_reorder_4d(1.0E0_realk,t2virt%elm4,nvEOS,noAOS,nvEOS,noAOS,[1,2,3,4],0.0E0_realk,temp2%val)

    !form X intermediate
    Xij = array2_init([noAOS,noAOS])
    call array4_contract3(temp2,temp1,Xij)

    !clean up!
    call array4_free(temp1)
    call array4_free(temp2)
     
    !add to the SD density
    do i=1,noAOS
       do j=1,noAOS

          dij%val(i,j) = dij%val(i,j) - Xij%val(i,j)

       end do
    end do

    !fill in the grad structure 
    grad%dens%X(1:noAOS,1:noAOS) = dij%val(1:noAOS,1:noAOS)

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, 'Xij intermediate'
       call print_matrix(Xij)
       print *, "";       print *, "DEC SINGLES: occ-occ block of unrel. one-el. dens."
       call print_matrix_real(grad%dens%X,noAOS,noAOS)
       print *, "";    endif

    !clean up!
    call array2_free(Xij)
    call array2_free(dij)

    !time
    call LSTIMER('X MATRIX',tcpu,twall,DECinfo%output)

  end subroutine single_dab_dij_formation

  !> ----------------------------------------------------------------
  !> \brief For a given fragment, calculate Dai and Dia off-diagonal
  !>        parts of the unrelaxed one-electron density 
  !>
  !> ON INPUT
  !>            MyFragment       - DEC info for a fragment
  !>            m2occ, m2virt,m1 - Lagrange multipliers
  !>            t2occ, t2virt    - amplitudes
  !> ON OUTPUT
  !>            dai,dia          - Dai and Dia off-diagonal parts 
  !>                              of the unrel. dens.
  !>            grad             - grad is updated with Dai and Dia
  !> \author Dmytro Bykov
  !> \date   October 2014
  !> ----------------------------------------------------------------
  subroutine single_dai_dia_formation(MyFragment,grad,t2occ,t2virt,m1)
   
    implicit none

    !INPUT declarations
    !******************
    type(decfrag), intent(inout) :: MyFragment
    type(tensor),intent(inout)    :: t2occ
    type(tensor),intent(inout)    :: t2virt
    type(tensor),intent(inout)    :: m1 
    !DEBUG use MP2 gradient structure for now
    type(mp2grad),intent(inout) :: grad

    !OUTPUT
    !******
    type(array2)   :: dai 
    type(array2)   :: dia

    !LOCAL INTERMEDIATES
    !*******************

    !Service variables
    type(array2)   :: tmp1,tmp2 
    type(array4)   :: temp1, temp2
    integer        :: nb,no,nv,noEOS,nvEOS,noAOS,nvAOS
    integer        :: iEOS,aEOS,jEOS,bEOS,iAOS,jAOS,aAOS,bAOS,iMOL,aMOL 
    logical        :: something_wrong
    real(realk)    :: tcpu,twall

    nb    = MyFragment%nbasis
    !DEBUG 
    no    = MyFragment%noccAOS
    nv    = MyFragment%nvirtAOS
    !DEBUG
    noEOS = MyFragment%noccEOS
    noAOS = MyFragment%noccAOS
    nvEOS = MyFragment%nvirtEOS
    nvAOS = MyFragment%nvirtAOS

    !start the stopwatch
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    !Initialize SD density
    dai  = array2_init([nvAOS,noAOS])
    dia  = array2_init([noAOS,nvAOS])
    tmp1 = array2_init([grad%dens%nbasis,grad%dens%nbasis])
    tmp2 = array2_init([grad%dens%nbasis,grad%dens%nbasis])
    
    !clean up first
    grad%dens%rho = 0e0_realk

    ! Sanity check 1: Density input structure
    !****************************************
    something_wrong=.false.
    if(noAOS/=grad%dens%nocc)   something_wrong=.true.
    if(nvAOS/=grad%dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'dens%nocc  ', grad%dens%nocc
       write(DECinfo%output,*) 'dens%nvirt', grad%dens%nvirt
       write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noAOS
       write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvAOS
       call lsquit('single_dai_dia_formation: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

    ! Sanity check 2: Array dimensions
    !*********************************
    if(nvAOS /= t2occ%dims(3)) something_wrong=.true.
    if(noEOS /= t2occ%dims(4)) something_wrong=.true.
    if(nvEOS /= t2virt%dims(3)) something_wrong=.true.
    if(noAOS /= t2virt%dims(4)) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'EOS: nocc, nvirt', noEOS, nvEOS
       write(DECinfo%output,*) 'AOS: nocc, nvirt', noAOS, nvAOS
       write(DECinfo%output,*) 't2occ dims     ', t2occ%dims
       write(DECinfo%output,*) 't2virt dims    ', t2virt%dims
       call lsquit('single_dai_dia_formation: &
            & Something wrong with amplitude and Theta dimensions!',DECinfo%output)
    end if
    
    !******************************
    !Do the calculation           *
    !******************************
    ! virt-occ block
    !******************
    ! m1 contains the multipliers: just need to put them into the full molecular 
    ! orbital space 
    !Dai = z^a_i
    do aEOS=1,nvEOS
       do iEOS=1,noEOS

          iAOS=MyFragment%idxo(iEOS) 
          aAOS=MyFragment%idxu(aEOS) 

          iMOL=MyFragment%occEOSidx(iEOS) 
          aMOL=MyFragment%virtEOSidx(aEOS) 

          grad%dens%Phivo(aAOS,iAOS) = m1%elm2(aAOS,iAOS)

       end do
    end do

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "DEC SINGLES: virt-occ block of unrel. one-el. dens."
       print *, "";       print *,no,nv
       print *, "";       call print_matrix_real(grad%dens%Phivo,nv,no)
       print *, "";    endif

    
    !clean up if needed


    !time
    call LSTIMER('Dai MATRIX',tcpu,twall,DECinfo%output)

    ! occ-virt block
    !****************
    !using the virtual partitioning scheme
    !D^sd_ia = Sum_bj z^b_j*(2t^ab_ij - t^ba_ij)
    do iEOS=1,noEOS
       !iMOL=MyFragment%occEOSidx(iEOS)
       iAOS=MyFragment%idxo(iEOS)
       
      do aAOS=1,nvAOS
         !aMOL=MyFragment%unoccAOSidx(aAOS)
            
         do jEOS=1,noEOS
            jAOS=MyFragment%idxo(jEOS)
            
           do bAOS=1,nvAOS

              grad%dens%Phiov(iAOS,aAOS) = grad%dens%Phiov(iAOS,aAOS) + m1%elm2(bAOS,jAOS)*&
                                    &(2E0_realk*t2occ%elm4(aAOS,iEOS,bAOS,jEOS)&
                                    & - t2occ%elm4(bAOS,iEOS,aAOS,jEOS))
 
           end do
         end do
      end do
    end do 

    if( DECinfo%PL>2 )then
       print *, "";       print *, "DEC SINGLES: occ-virt block of unrel. one-el. dens. "
       call print_matrix(dia)
       print *, "";       call print_matrix_real(grad%dens%Phiov,no,nv)
       print *, "";    endif

    !clean up!

    !time
    call LSTIMER('Dia MATRIX',tcpu,twall,DECinfo%output)

  end subroutine single_dai_dia_formation

  !> ----------------------------------------------------------------
  !> \brief Calculate contributions to gradient  
  !>        from the pair fragment
  !> ON INPUT
  !>            Fragment1 & 2    - DEC info for fragments
  !>            PairFragment     - DEC info for a pair
  !>            m2occ, m2virt,m1 - Lagrange multipliers
  !>            t2occ, t2virt    - amplitudes
  !> ON OUTPUT
  !>            grad             - gradint structure 
  !> \author Dmytro Bykov
  !> \date   October 2014
  !> ----------------------------------------------------------------
  subroutine pair_calculate_CCSDgradient_driver&
                            &(Fragment1,Fragment2,PairFragment,t2occ,t2virt,m2occ,m2virt,m1,&
                                      &VOOO,VOVV,VOVOocc,VOVOvirt,grad)
  
    implicit none

    !INPUT declarations
    !******************
    type(decfrag),intent(in)  :: Fragment1
    type(decfrag),intent(in)  :: Fragment2
    type(decfrag),intent(inout)  :: PairFragment
    type(tensor),intent(inout)    :: m2occ
    type(tensor),intent(inout)    :: m2virt
    type(tensor),intent(inout)    :: t2occ
    type(tensor),intent(inout)    :: t2virt
    type(tensor),intent(inout)    :: m1 
    type(tensor),intent(in)       :: VOOO
    type(tensor),intent(in)       :: VOVV
    type(tensor),intent(inout)    :: VOVOocc   
    type(tensor),intent(in)       :: VOVOvirt

    !OUTPUT
    !******
    !DEBUG try to use MP2 gradient structure
    type(mp2grad),intent(inout) :: grad

    !LOCAL INTERMEDIATES
    !*******************

    !Service variables
    type(array2)      :: Xij
    type(array4)      :: temp1, temp2
    integer           :: i,j,nb,noEOS,nvEOS,noAOS,nvAOS

    nb    = PairFragment%nbasis
    noEOS = PairFragment%noccEOS
    noAOS = PairFragment%noccAOS
    nvEOS = PairFragment%nvirtEOS
    nvAOS = PairFragment%nvirtAOS

    !DEBUG Init MP2 gradient structure
    call init_mp2grad(PairFragment,grad)

    !DEBUG DEBUG DEBUG
    !call pair_calculate_DEBUG(Fragment1,Fragment2,PairFragment,t2occ,t2virt,&
    !   & m2occ, m2virt, VOOO,VOVV,grad%dens)

    !calculate the density
    !**********************
    !the diagonal blocks dab and dij of the one-electron density 
    call pair_dab_dij_formation(Fragment1,Fragment2,PairFragment,grad,t2occ,t2virt,m2occ,m2virt)

    !the off-diagonal blocks dia and dai
    call pair_dai_dia_formation(Fragment1,Fragment2,PairFragment,grad,t2occ,t2virt,m1)
    
    !clean up if needed

  end subroutine pair_calculate_CCSDgradient_driver

  !> ----------------------------------------------------------------
  !> \brief For a pair fragment, calculate Dab and Dij parts of the 
  !>        unrelaxed one-electron density 
  !>
  !> ON INPUT
  !>            PairFragment      - DEC info for a pair fragment
  !>            Fragment1 & 2     - DEC info for fragments
  !>            m2occ, m2virt, m1 - Lagrange multipliers
  !>            t2occ, t2virt     - amplitudes
  !> ON OUTPUT
  !>            dab,dij           - Dab and Dij parts of the unrel. dens.
  !>            grad              - grad is updated with Dab and Dij
  !> \author Dmytro Bykov
  !> \date   October 2014
  !> ----------------------------------------------------------------
  subroutine pair_dab_dij_formation(Fragment1,Fragment2,PairFragment,grad,t2occ,t2virt,m2occ,m2virt)
   
    implicit none

    !INPUT declarations
    !******************
    type(decfrag),intent(in)  :: Fragment1
    type(decfrag),intent(in)  :: Fragment2
    type(decfrag),intent(inout)  :: PairFragment
    type(tensor),intent(inout)    :: m2occ
    type(tensor),intent(inout)    :: m2virt
    type(tensor),intent(inout)    :: t2occ
    type(tensor),intent(inout)    :: t2virt

    !OUTPUT
    !******
    !DEBUG try to use MP2 gradient structure
    type(mp2grad),intent(inout) :: grad

    !LOCAL INTERMEDIATES
    !*******************
    type(array2)   :: dab 
    type(array2)   :: dij

    !Service variables
    type(array2)     :: Xij
    type(array4)     :: temp1, temp2
    integer          :: i,j,nb,noEOS,nvEOS,noAOS,nvAOS
    integer          :: a,b,c,d,k,l
    integer          :: iFull,jFull
    logical          :: something_wrong
    real(realk)      :: tcpu,twall
    logical, pointer :: dopair_occ(:,:), dopair_virt(:,:)

    nb    = PairFragment%nbasis
    noEOS = PairFragment%noccEOS
    noAOS = PairFragment%noccAOS
    nvEOS = PairFragment%nvirtEOS
    nvAOS = PairFragment%nvirtAOS

    ! Which "interaction pairs" to include for occ and virt space (avoid double counting)
    call mem_alloc(dopair_occ,noEOS,noEOS)
    call mem_alloc(dopair_virt,nvEOS,nvEOS)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)

    !start the stopwatch
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    !Initialize SD density
    dab  = array2_init([nvAOS,nvAOS])
    dij  = array2_init([noAOS,noAOS])
    Xij  = array2_init([noAOS,noAOS])

    !Zero matrices in dens structure
    grad%dens%X = 0e0_realk
    grad%dens%Y = 0e0_realk

    ! Sanity check 1: Density input structure
    !****************************************
    something_wrong=.false.
    if(noAOS/=grad%dens%nocc)   something_wrong=.true.
    if(nvAOS/=grad%dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'dens%nocc  ', grad%dens%nocc
       write(DECinfo%output,*) 'dens%nvirt', grad%dens%nvirt
       write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noAOS
       write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvAOS
       call lsquit('pair_dab_dij_formation: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

    ! Sanity check 2: Array dimensions
    !*********************************
    if(nvAOS /= t2occ%dims(3)) something_wrong=.true.
    if(noEOS /= t2occ%dims(4)) something_wrong=.true.
    if(nvEOS /= t2virt%dims(3)) something_wrong=.true.
    if(noAOS /= t2virt%dims(4)) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'EOS: nocc, nvirt', noEOS, nvEOS
       write(DECinfo%output,*) 'AOS: nocc, nvirt', noAOS, nvAOS
       write(DECinfo%output,*) 't2occ dims     ', t2occ%dims
       write(DECinfo%output,*) 't2virt dims    ', t2virt%dims
       call lsquit('pair_dab_dij_formation: &
            & Something wrong with amplitude dimensions!',DECinfo%output)
    end if
    
    !******************************
    !Do the calculation           *
    !******************************
    ! virt-virt block
    !******************
    !using the occupied partitioning scheme
    !Dab = Sum_ckl (t^cb_kl * z^ca_kl) = Y_ba

    !form virt-virt dens. block 
    do k=1,noEOS
      do l=1,noEOS

          ! Only update for "interaction orbital pairs" - see which_pairs_occ
          if(dopair_occ(k,l)) then !YDoPair

            do c=1,nvAOS
              do b=1,nvAOS
                do a=1,nvAOS
                  dab%val(b,a) = dab%val(b,a) + t2occ%elm4(c,k,b,l)*m2occ%elm4(c,k,a,l)
                end do
              end do
            end do

          end if

      end do
    end do

    !fill in the grad structure 
    grad%dens%Y(1:nvAOS,1:nvAOS) = dab%val(1:nvAOS,1:nvAOS)

    !fill in the rho
    !do i=1, PairFragment%NvirtAOS !nvAOS
    !  iFull = PairFragment%virtAOSidx(i)
    !  do j=1,PairFragment%NvirtAOS !nvAOS
    !    jFull = PairFragment%virtAOSidx(j)
    !    grad%dens%rho(iFull,jFull) = dab%val(i,j)
    !  end do
    !end do

    !clean up!
    call array2_free(dab)

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "DEC PAIRS: virt-virt block of unrel. one-el. dens."
       call print_matrix_real(grad%dens%Y,nvAOS,nvAOS)
       print *, "";    endif

    !time
    call LSTIMER('Y MATRIX',tcpu,twall,DECinfo%output)

    ! occ-occ block
    !****************
    !using the virtual partitioning scheme
    !Dij = 2*delta_ij - Sum_cdk (t^cd_ki * z^cd_kj) = 2*delta_ij - X_ij

    do c=1,nvEOS
      do d=1,nvEOS

        ! Only update for "interaction orbital pairs" - see which_pairs_virt
        if(dopair_virt(c,d)) then !XDoPair

          do k=1,noAOS
            do i=1,noAOS
              do j=1,noAOS

                Xij%val(i,j) = Xij%val(i,j) + t2virt%elm4(c,k,d,i)*m2virt%elm4(c,k,d,j)

              end do
            end do
          end do

          end if

       end do
    end do

    !Fill in the dij
    do i=1,noAOS
     do j=1,noAOS

      dij%val(i,j) = - Xij%val(i,j)

      !if(i==j)then
      ! dij%val(i,j) = 2E0_realk/numPairFragments - Xij%val(i,j)  
      !else
      ! dij%val(i,j) = - Xij%val(i,j)
      !endif

     end do
    end do
 
    !fill in the grad structure 
    grad%dens%X(1:noAOS,1:noAOS) = dij%val(1:noAOS,1:noAOS)

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "DEC PAIR: Xij intermediate"
       call print_matrix(Xij)
       print *, "";       print *, "DEC PAIR: occ-occ block of unrel. one-el. dens."
       call print_matrix_real(grad%dens%X,noAOS,noAOS)
       print *, "";    endif

    !clean up!
    call array2_free(dij)
    call array2_free(Xij)

    !time
    call LSTIMER('X MATRIX',tcpu,twall,DECinfo%output)

  end subroutine pair_dab_dij_formation

  !> ----------------------------------------------------------------
  !> \brief For a pair fragment, calculate Dai and Dia off-diagonal
  !>        parts of the unrelaxed one-electron density 
  !>
  !> ON INPUT
  !>            PairFragment      - DEC info for a fragment
  !>            m2occ, m2virt, m1 - Lagrange multipliers
  !>            t2occ, t2virt     - amplitudes
  !> ON OUTPUT
  !>            dai,dia           - Dai and Dia off-diagonal parts 
  !>                                of the unrel. dens.
  !>            grad              - grad is updated with Dai and Dia
  !> \author Dmytro Bykov
  !> \date   October 2014
  !> ----------------------------------------------------------------
  subroutine pair_dai_dia_formation(Fragment1,Fragment2,PairFragment,grad,t2occ,t2virt,m1)
   
    implicit none

    !INPUT declarations
    !******************
    type(decfrag),intent(in)  :: Fragment1
    type(decfrag),intent(in)  :: Fragment2
    type(decfrag),intent(inout)  :: PairFragment
    type(tensor),intent(inout)    :: t2occ
    type(tensor),intent(inout)    :: t2virt
    type(tensor),intent(inout)    :: m1 
    !DEBUG try to use MP2 gradient structure
    type(mp2grad),intent(inout) :: grad

    !OUTPUT
    !******
    type(array2)   :: dai 
    type(array2)   :: dia

    !LOCAL INTERMEDIATES
    !*******************

    !Service variables
    type(array2)     :: tmp1,tmp2 
    type(array4)     :: temp1, temp2
    integer          :: nb,no,nv,noEOS,nvEOS,noAOS,nvAOS
    integer          :: iEOS,aEOS,jEOS,bAOS,bEOS,iAOS,jAOS,aAOS,iMOL,aMOL 
    logical          :: something_wrong
    real(realk)      :: tcpu,twall
    logical, pointer :: dopair_occ(:,:), dopair_virt(:,:)
    logical, pointer :: dopair_virt_occ(:,:)
    logical, pointer :: dopair_virt_occ_aos(:,:)
 
    nb    = PairFragment%nbasis
    !DEBUG - works only for simulate full mol
    no    = PairFragment%noccAOS
    nv    = PairFragment%nvirtAOS
    !DEBUG
    noEOS = PairFragment%noccEOS
    noAOS = PairFragment%noccAOS
    nvEOS = PairFragment%nvirtEOS
    nvAOS = PairFragment%nvirtAOS

    ! Which "interaction pairs" to include for occ and virt space (avoid double counting)
    call mem_alloc(dopair_occ,noEOS,noEOS)
    call mem_alloc(dopair_virt,nvEOS,nvEOS)
    call mem_alloc(dopair_virt_occ,nvEOS,noEOS)
!    call mem_alloc(dopair_virt_occ_aos,nvAOS,noAOS)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)
    call which_pairs_occ_virt(Fragment1,Fragment2,PairFragment,dopair_virt_occ)
!    call which_pairs_occ_virt_aos(Fragment1,Fragment2,PairFragment,dopair_virt_occ_aos)

    !start the stopwatch
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    !Initialize SD density
    dai  = array2_init([nvAOS,noAOS])
    dia  = array2_init([noAOS,nvAOS])
    tmp1 = array2_init([grad%dens%nbasis,grad%dens%nbasis])
    tmp2 = array2_init([grad%dens%nbasis,grad%dens%nbasis])

    !clean up first
    grad%dens%rho = 0e0_realk

    ! Sanity check 1: Density input structure
    !****************************************
    something_wrong=.false.
    if(noAOS/=grad%dens%nocc)   something_wrong=.true.
    if(nvAOS/=grad%dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'dens%nocc  ', grad%dens%nocc
       write(DECinfo%output,*) 'dens%nvirt', grad%dens%nvirt
       write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noAOS
       write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvAOS
       call lsquit('pair_dai_dia_formation: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

    ! Sanity check 2: Array dimensions
    !*********************************
    if(nvAOS /= t2occ%dims(3)) something_wrong=.true.
    if(noEOS /= t2occ%dims(4)) something_wrong=.true.
    if(nvEOS /= t2virt%dims(3)) something_wrong=.true.
    if(noAOS /= t2virt%dims(4)) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'EOS: nocc, nvirt', noEOS, nvEOS
       write(DECinfo%output,*) 'AOS: nocc, nvirt', noAOS, nvAOS
       write(DECinfo%output,*) 't2occ dims     ', t2occ%dims
       write(DECinfo%output,*) 't2virt dims    ', t2virt%dims
       call lsquit('pair_dai_dia_formation: &
            & Something wrong with amplitude dimensions!',DECinfo%output)
    end if
    
    !******************************
    !Do the calculation           *
    !******************************
    ! virt-occ block
    !******************
    ! m1 contains the multipliers: just need to put them into the full molecular 
    ! orbital space 
    !Dai = z^a_i
    do aEOS=1,nvEOS
       do iEOS=1,noEOS

        ! Only update for "interaction orbital pairs":
        ! see which_pairs_occ_virt
!        print*, dopair_virt_occ(aEOS,iEOS), 'dopair_virt_occ'
        if(dopair_virt_occ(aEOS,iEOS)) then 
          
          iAOS=PairFragment%idxo(iEOS) 
          aAOS=PairFragment%idxu(aEOS) 
               
          iMOL=PairFragment%occEOSidx(iEOS) 
          aMOL=PairFragment%virtEOSidx(aEOS) 

          !dai%val(aMOL,iMOL) = m1%elm2(aAOS,iAOS)
          grad%dens%Phivo(aAOS,iAOS) = m1%elm2(aAOS,iAOS)

        end if

       end do
    end do

    !fill in the grad structure 
    !grad%dens%rho((noAOS+1):grad%dens%nbasis,1:noAOS) = dai%val(1:nvAOS,1:noAOS)
 
    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "DEC PAIR: virt-occ block of unrel. one-el. dens. "
       print *, "";       print *,no,nv
       print *, "";       call print_matrix_real(grad%dens%Phivo,nv,no)
       print *, "";    endif

    
    !clean up if needed


    !time
    call LSTIMER('Dai MATRIX',tcpu,twall,DECinfo%output)

    ! occ-virt block
    !****************
    !D^sd_ia = Sum_bj z^b_j*(2t^ab_ij - t^ba_ij)
    do iEOS=1,noEOS
       !iMOL=PairFragment%occEOSidx(iEOS)
       iAOS=PairFragment%idxo(iEOS)

      do aAOS=1,nvAOS
         !aMOL=PairFragment%unoccAOSidx(aAOS)

           do jEOS=1,noEOS
              jAOS=PairFragment%idxo(jEOS)

             do bAOS=1,nvAOS
         
                ! Only update for "interaction orbital pairs":
                ! see which_pairs_occ_virt
                if(dopair_occ(iEOS,jEOS)) then 
          
                grad%dens%Phiov(iAOS,aAOS) = grad%dens%Phiov(iAOS,aAOS)  + m1%elm2(bAOS,jAOS)*&
                                    &(2E0_realk*t2occ%elm4(aAOS,iEOS,bAOS,jEOS)&
                                    & - t2occ%elm4(bAOS,iEOS,aAOS,jEOS))

                end if
             end do
           end do
      end do
    end do
 
    if( DECinfo%PL>2 )then
       print *, "";       print *, "DEC PAIR: occ-virt block of unrel. one-el. dens. "
       print *, "";       call print_matrix_real(grad%dens%rho,no+nv,no+nv)
       print *, "";    endif

    !clean up!

    !time
    call LSTIMER('Dia MATRIX',tcpu,twall,DECinfo%output)

  end subroutine pair_dai_dia_formation

  !> ----------------------------------------------------------------
  !> \brief For a given fragment, calculate RHS of z-vector equations
  !>
  !> ON INPUT
  !>            MyFragment      - DEC info for a fragment
  !> ON OUTPUT
  !>            eta_arr         - RHS of Z-vector equations
  !> \author Dmytro Bykov
  !> \date   May 2014
  !> ----------------------------------------------------------------
  subroutine z_vec_rhs_ccsd(MyFragment,eta,m1_arr,m2_arr,t1_arr,t2_arr)

     implicit none

     !INPUT declarations
     !******************

     !> Fragment info
     type(decfrag), intent(inout) :: MyFragment
   !> Singles multipliers m1(a,i)
   type(tensor),intent(inout), optional :: m1_arr
   !> Doubles multipliers m2(a,i,b,j)
   type(tensor),intent(inout), optional :: m2_arr
   !> Singles amplitudes t1(a,i)
   type(tensor),intent(inout), optional :: t1_arr
   !> Doubles amplitudes t2(a,i,b,j)
   type(tensor),intent(inout), optional :: t2_arr
     !> Singles multipliers m1(a,i)
!     type(array) :: m1_arr
     !> Doubles multipliers m2(a,i,b,j)
!     type(array) :: m2_arr
!     real(realk),intent(inout) :: m1_arr(:,:),m2_arr(:,:,:,:)

     !OUTPUT
     !******

     !RHS
     type(tensor)   :: eta
     type(array2)  :: eta_arr, eta_arr2

     !LOCAL INTERMEDIATES
     !*******************

     !SD density container
     type(array2) :: Dsd
     type(array4) :: d2el,d2elaibj,d2elaibj2

     !Service variables
     character(TENSOR_MSG_LEN)    :: msg
     real(realk)               :: norm,nrmm2
     real(realk), pointer      :: TEST(:)
     type(array2)              :: Xij, Yba, Eai, temp
     type(array4)              :: tbar, temp1, temp2, temp3
     type(array4)              :: Mbar, M1bar, M2bar, M3bar, M4bar, M5bar
     type(array4)              :: Yae,Xim,Xjm,Ybc
     type(array4)              :: gao
     type(array4)              :: goovv,govvv,gvovv,gvvvv
     type(array4)              :: goooo,govoo,gvooo,gvvoo
     type(array4)              :: goovo,govvo,gvovo,gvvvo
     type(array4)              :: gooov,govov,gvoov,gvvov
     integer                   :: nb,no,nv
     integer                   :: v4,o4,o3v,ov3,o2v2,ov,b2,v2,o2,o2v,ov2
     integer                   :: i,j,k,l,m,a,b,c,d,t,r,s
     integer, dimension(2)     :: dims
     type(matrix)              :: h1,tmp1,tmp2
     type(array2)              :: xocc,yocc,xvirt,yvirt
     type(array2)              :: Co,Cv,Co2,Cv2
     type(array2)              :: t1_trans
     type(array2)              :: hoo,hov,hvo,hvv
     real(realk), pointer      :: w1(:)
     type(array2)              :: Co_d,Cv_d,Uocc,Uvirt,Co_dl,Co_dr,Cv_dl,Cv_dr
     real(realk), pointer      :: focc(:),fvirt(:)

     !INTERNAL PARAMETERS
     logical :: local

     nb = MyFragment%nbasis
     no = MyFragment%noccAOS
     nv = MyFragment%nvirtAOS

     b2   = nb*nb
     v2   = nv*nv
     o2   = no*no
     o2v  = o2*nv
     ov   = no*nv
     ov2  = ov*nv
     o2v2 = ov*ov
     ov3  = ov*v2
     o3v  = ov*o2
     o4   = o2*o2
     v4   = v2*v2

     dims = [nb,nb]

     !Initialize SD density
     Dsd  = array2_init(dims)

     !memory allocations
!     call mem_alloc(temp2,max(max(no,nv)*nb**3,max(max(max(max(o2v2,ov3),v4),o2*v2),o4)))
!     call mem_alloc(temp1,max(max(max(max(o2v2,ov3),v4),o2*v2),o4))
!     call mem_alloc(Yba,  max(max(max(max(o2v2,ov3),v4),o2*v2),o4))

     call mem_alloc(TEST, 8)

     !***************************************************************
     !DEBUG Printout
     !***************************************************************
!     if( DECinfo%PL>2 )then
!
!        write (msg, *)"DEBUG z1 n**2 n"
!        call print_norm(m1_arr%elm2,int(ov,kind=8),norm,.true.)
!        print *,msg,norm,sqrt(norm)
!
!        nrmm2 = 0.0E0_realk
!        do j = 1, no
!          do b = 1, nv
!            do i = 1, no
!             do a = 1,nv
!               if(a+(i-1)*nv<=b+(j-1)*nv)then
!                 nrmm2 = nrmm2 + m2_arr%elm4(a,i,b,j)**2
!               endif
!             enddo
!           enddo
!         enddo
!       enddo
!
!       write (msg,*)"z2 n**2 n pn**2 pn"
!       call print_norm(m2_arr%elm4,int(o2v2,kind=8),norm,.true.)
!       print *,msg,norm,sqrt(norm),nrmm2,sqrt(nrmm2)
!    endif
!
!    !Print just to see that it's nicely empty
!    if( DECinfo%PL>2 )then
!       print *, "";!       call print_matrix(Dsd)
!       print *, "";!    endif

     ! Sanity check: This routine is not intended for MP2
     if(MyFragment%ccmodel == MODEL_MP2) then
        call lsquit('fragment_ccsolver cannot be used for MP2!',&
        & DECinfo%output)
     end if

     local=.true.
#ifdef VAR_MPI
     if(infpar%lg_nodtot>1)local=.false.
#endif

!*******************************************************************************
!DEBUG: put all ones to the input multipliers to be sure not to get exceptions
!from the incoming buggy arrays.

!     do a=no+1,nb
!        do i=1,no
!
!          m1_arr%elm2(a-no,i) = 1
!
!        end do
!     end do
!
!     do a=no+1,nb
!        do i=1,no
!          do b=no+1,nb
!            do j=1,no
!
!              m2_arr%elm4(a-no,i,b-no,j) = 1
!
!            end do
!          end do
!        end do
!     end do
!
     !DEBUG: hard-coded multipliers
!
!        do a=1,m1_arr%nelms 
!            m1_arr%elm1(a)=0.0
!        end do
!
!        do a=1,m2_arr%nelms 
!            m2_arr%elm1(a)=0.0
!        end do
!
!        m1_arr%elm1(1)= 0.00120460
!        m1_arr%elm1(4)= 0.00081651
!        m1_arr%elm1(5)= 0.05956286
!        m1_arr%elm1(8)=-0.00181810
!
!       m2_arr%elm1(1)= -0.00616018
!       m2_arr%elm1(4)=  0.00408643
 !      m2_arr%elm1(5)= -0.00132998
 !      m2_arr%elm1(8)= -0.01766175
 !      m2_arr%elm1(10)= -0.00332510
 !      m2_arr%elm1(14)= -0.00583353
 !      m2_arr%elm1(19)= -0.00332510
 !      m2_arr%elm1(23)= -0.00583353
 !      m2_arr%elm1(25)=  0.00408643
 !      m2_arr%elm1(28)= -0.00453717
 !      m2_arr%elm1(29)=  0.00189826
 !      m2_arr%elm1(32)= -0.00730401
 !      m2_arr%elm1(33)= -0.00132998
 !      m2_arr%elm1(36)=  0.00189826
 !      m2_arr%elm1(37)= -0.04139758
 !      m2_arr%elm1(40)= -0.08471432
 !      m2_arr%elm1(42)= -0.00583353
 !      m2_arr%elm1(46)= -0.07998981
 !      m2_arr%elm1(51)= -0.00583353
 !      m2_arr%elm1(55)= -0.07998981
 !      m2_arr%elm1(57)= -0.01766175
 !      m2_arr%elm1(60)= -0.00730401
 !      m2_arr%elm1(61)= -0.08471432
 !      m2_arr%elm1(64)= -0.19593169
 !
 !       do a=1,t1_arr%nelms 
 !           t1_arr%elm1(a)=0.0
 !       end do
 !
 !       do a=1,t2_arr%nelms 
 !           t2_arr%elm1(a)=0.0
 !       end do
 !
 !t1_arr%elm1(1)=0.00079127
 !t1_arr%elm1(4)=0.00032443
 !t1_arr%elm1(5)=0.03098741
 !t1_arr%elm1(8)=0.00050793
 !
 !t2_arr%elm1( 1)=-0.00306775  375     1 -0.00306775
 !t2_arr%elm1( 3)=-0.00166009  376     4 -0.00205805
 !t2_arr%elm1( 6)=-0.00166009  377     5 -0.00069422
 !t2_arr%elm1( 7)= 0.00205805  378     8  0.00568419
 !t2_arr%elm1(10)=-0.00229763  379    10 -0.00166009
 !t2_arr%elm1(11)=-0.00069422  380    14 -0.00296272
 !t2_arr%elm1(14)=-0.00236328  381    19 -0.00166009
 !t2_arr%elm1(15)=-0.02206757  382    23 -0.00296272
 !t2_arr%elm1(17)=-0.00296272  383    25 -0.00205805
 !t2_arr%elm1(21)=-0.04078488  384    28 -0.00229763
 !t2_arr%elm1(24)=-0.00296272  385    29  0.00236328
 !t2_arr%elm1(28)=-0.04078488  386    32 -0.00370535
 !t2_arr%elm1(29)=-0.00568418  387    33 -0.00069422
 !t2_arr%elm1(32)=-0.00370535  388    36  0.00236328
 !t2_arr%elm1(33)=-0.04318715  389    37 -0.02206754
 !t2_arr%elm1(36)=-0.09988663  390    40  0.04318712
 !t2_arr%elm1(37)=-0.00306775  391    42 -0.00296272
 !t2_arr%elm1(39)=-0.00166009  392    46 -0.04078487
 !t2_arr%elm1(42)=-0.00166009  393    51 -0.00296272
 !t2_arr%elm1(43)= 0.00205805  394    55 -0.04078487
 !t2_arr%elm1(46)=-0.00229763  395    57  0.00568419
 !t2_arr%elm1(47)=-0.00069422  396    60 -0.00370535
 !t2_arr%elm1(50)= 0.00095763  397    61  0.04318712
 !t2_arr%elm1(51)=-0.02206757  398    64 -0.09988664
 !t2_arr%elm1(53)=-0.00296272
 !t2_arr%elm1(57)=-0.04078488
 !t2_arr%elm1(60)=-0.00296272
 !t2_arr%elm1(64)=-0.04078488


!        print *, "Debug print m1_ai: "
!        do a=1,m1_arr%nelms 
!          if (ABS(m1_arr%elm1(a))>0.00000001)then
!            write (*,'(i5,F12.8)') a, m1_arr%elm1(a)
!          endif
!        end do
!
!        print *, "Debug print m2_aibj: "
!        do a=1,m2_arr%nelms 
!          if (ABS(m2_arr%elm1(a))>0.00000001)then 
!            write (*,'(i5,F12.8)') a, m2_arr%elm1(a)
!          endif
!        end do
!        
!        print *, "Debug print t1_aibj: "
!        do a=1,t1_arr%nelms 
!          if (ABS(t1_arr%elm1(a))>0.00000001)then 
!            write (*,'(i5,F12.8)') a, t1_arr%elm1(a)
!          endif
!        end do
!
!        print *, "Debug print t2_aibj: "
!        do a=1,t2_arr%nelms 
!          if (ABS(t2_arr%elm1(a))>0.00000001)then 
!            write (*,'(i5,F12.8)') a, t2_arr%elm1(a)
!          endif
!        end do

!        print *, "t2_aibj like in DALTON"
!        k = 1
!        do j=1,no 
!         do b=1,nv
!          do i=1,no 
!           do a=1,nv
!            if (ABS(t2_arr%elm4(a,i,b,j))>0.00000001)then 
!             write (*,'(i5,i5,i5,i5,i5,F12.8)') k,a,i,b,j, t2_arr%elm4(a,i,b,j)
!            endif
!           k = k + 1
!          end do
!         end do
!        end do
!        end do

!*******************************************************************************

     !***************************************************************
     ! Construct the one-electron SD density
     !***************************************************************

     !1. virt-occ block
     !*****************
     !D^sd_ai = z^a_i
     do a=no+1,nb
        do i=1,no
            Dsd%val(a,i) = m1_arr%elm2(a-no,i)
        end do
     end do

     !print debug info if needed
     if( DECinfo%PL>2 )then
        print *, "";        print *, "Debug print Dsd virt-occ block: "
        print *, "";        print *,no,nv
        call print_matrix(Dsd)
        print *, "";     endif

     !2. virt-virt block
     !******************
     !D^sd_ab = Sum_ckl (t^cb_kl * z^ca_kl) = Y_ba

     !container for amps voov shape
     temp1 = array4_init([nv,no,no,nv])
     temp2 = array4_init([nv,no,no,nv])

     !sort amps(ckbl) -> bckl
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,2,4,3],0.0E0_realk,temp1%val)
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[1,2,4,3],0.0E0_realk,temp2%val)

     !form Y_ba intermediate
     Yba = array2_init([nv,nv])
     call array4_contract3(temp1,temp2,Yba)

     !clean up!
     call array4_free(temp1)
     call array4_free(temp2)

     !D^sd_ab = Yba
     do a=no+1,nb
        do b=no+1,nb
            Dsd%val(a,b) = 1.0*Yba%val(a-no,b-no)
        end do
     end do

     !print debug info if needed
     if( DECinfo%PL>2 )then
        print *, "";        print *, 'Yba intermediate'
        print *, no
        print *, nv
        print *, nb
        call print_matrix(Yba)
        print *, "";        print *, "Debug print Dsd virt-virt block: "
        call print_matrix(Dsd)
        print *, "";     endif

     !3. occ-occ block
     !****************
     !D^sd_ij = 2*delta_ij - Sum_cdk (t^cd_ki * z^cd_kj) = 2*delta_ij - X_ij

     !container for amps vovo shape
     temp1 = array4_init([nv,no,nv,no])
     temp2 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp1%val)
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp2%val)

     !form X intermediate
     Xij = array2_init([no,no])
     call array4_contract3(temp2,temp1,Xij)

     !clean up!
     call array4_free(temp1)
     call array4_free(temp2)
     
     !add to the SD density
     do i=1,no
        do j=1,no
           if(i==j)then
              Dsd%val(i,j) = Dsd%val(i,j) + 2 - Xij%val(i,j)
           else
              Dsd%val(i,j) = Dsd%val(i,j) - Xij%val(i,j)
           endif
        end do
     end do

     !print debug info if needed
     if( DECinfo%PL>2 )then
        print *, "";        print *, 'Xij intermediate'
        call print_matrix(Xij)
        print *, "";        print *, "Debug print Dsd occ-occ block: "
        call print_matrix(Dsd)
        print *, "";     endif

     !4. occ-virt block
     !*****************
     !D^sd_ia = Sum_bj z^b_j*(2t^ab_ij - t^ba_ij) = Sum_bj z^b_j * tbar^ab_ij

     !intermediate container vovo shape
     tbar = array4_init([nv,no,nv,no])

     !Construct the tbar intermediate
     do j=1,t2_arr%dims(4)
       do b=1,t2_arr%dims(3)
          do i=1,t2_arr%dims(2)
             do a=1,t2_arr%dims(1)
                tbar%val(a,i,b,j) = tbar%val(a,i,b,j) + 2.0*t2_arr%elm4(a,i,b,j)&
                                     & - t2_arr%elm4(b,i,a,j)
             end do
          end do
       end do
     end do

     !form eai intermediate
     Eai = array2_init([nv,no])

     !container for amps
     temp = array2_init([nv,no],m1_arr%elm2)
     call array4_contract_array2(tbar,temp,Eai)
     call array2_free(temp)

!     do i=1,no
!       do a=1,nv
!          do j=1,no
!             do b=1,nv
!                Eai%val(a,i) = Eai%val(a,i) + temp%val(b,j)*&
!                                     &tbar%val(a,i,b,j)
!             end do
!          end do
!       end do
!     end do
 
    !add to the SD density
     do a=no+1,nb
        do i=1,no
            Dsd%val(i,a) = Eai%val(a-no,i)
        end do
     end do

     !print debug info if needed
     if( DECinfo%PL>2 )then
        print *, "";        print *, "Tbar intermediate"
        k = 1
        do j=1,tbar%dims(4) 
         do b=1,tbar%dims(3) 
          do i=1,tbar%dims(2) 
           do a=1,tbar%dims(1) 
            if (ABS(tbar%val(a,i,b,j))>0.00000001)then 
             write (*,'(i5,i5,i5,i5,i5,F12.8)') k,a,i,b,j, tbar%val(a,i,b,j)
            endif
            k = k + 1
           end do
          end do
         end do
        end do

        print *, "";        print *, 'Eai intermediate'
        call print_matrix(Eai)

        print *, "";        print *, "Debug print Dsd occ-virt block: "
        call print_matrix(Dsd)
        print *, "";     endif

     !***************************************************************
     ! Construct the two-electron SD density
     !***************************************************************

     !1. occ-occ-occ-occ block
     !************************
     !d^sd_ijkl = 4*delta_ij * delta_kl - 2*delta_kj * delta_il +
     !          +   delta_kj * Sum_cdm z^cd_ml * t^cd_mi -
     !          - 2*delta_ij * Sum_cdm z^cd_ml * t^cd_mk -
     !          - 2*delta_kl * Sum_cdm z^cd_mj * t^cd_mi +
     !          +   delta_il * Sum_cdm z^cd_mj * t^cd_mk +
     !          +              Sum_cd  z^cd_jl * t^cd_ik =
     !
     !          = 4*delta_ij * delta_kl - 2*delta_kj * delta_il +
     !          +   delta_kj * X_il -
     !          - 2*delta_ij * X_kl -
     !          - 2*delta_kl * X_ij +
     !          +   delta_il * X_kj +
     !          + Mbar_ijkl
 
    !container for amps vvoo shape
    temp1 = array4_init([nv,nv,no,no])
    temp2 = array4_init([nv,nv,no,no])
    temp3 = array4_init([no,no,no,no])
 
    !sort amps(cjdn) -> cdjn
    call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,3,2,4],0.0E0_realk,temp1%val)
    call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[1,3,2,4],0.0E0_realk,temp2%val)
 
    !form Mbar intermediate
    ! Mbar_ijmn = Sum_cd z^cd_jn * t^cd_im
    Mbar = array4_init([no,no,no,no])
    call array4_contract2(temp1,temp2,temp3)
    call array_reorder_4d(1.0E0_realk,temp3%val,no,no,no,no,[3,1,4,2],0.0E0_realk,Mbar%val)

     !clean up!
     call array4_free(temp1)
     call array4_free(temp2)
 
     !Initialize two-el density
    d2el  = array4_init([nb,nb,nb,nb])
    d2elaibj  = array4_init([nv,no,nv,nb])
 
     !construct d^sd_ijkl
     do i=1,no
       do j=1,no
         do k=1,no
           do l=1,no
 
             d2el%val(i,j,k,l) = d2el%val(i,j,k,l) + Mbar%val(i,j,k,l)
 
             if(i==l)then
               d2el%val(i,j,k,l) = d2el%val(i,j,k,l) + Xij%val(k,j)
               if(k==j)then
                 d2el%val(i,j,k,l) = d2el%val(i,j,k,l) - 2.0
               endif
             endif
 
             if(k==l)then
               d2el%val(i,j,k,l) = d2el%val(i,j,k,l) - 2.0*Xij%val(i,j)
               if(i==j)then
                 d2el%val(i,j,k,l) = d2el%val(i,j,k,l) + 4.0
               endif
             endif
 
             if(i==j)then
               d2el%val(i,j,k,l) = d2el%val(i,j,k,l) - 2.0*Xij%val(k,l)
             endif
 
             if(k==j)then
               d2el%val(i,j,k,l) = d2el%val(i,j,k,l) + Xij%val(i,l)
             endif
 
           end do
         end do
       end do
     end do
 
     !print debug info if needed
     if( DECinfo%PL>2 )then
 
      print *, "";      print *, "Mbar intermediate "
      m = 1
      do l=1,no
        do k=1,no
          do j=1,no
            do i=1,no
 
              if (ABS(Mbar%val(i,j,k,l))>0.00000001)then 
               write (*,'(i5,i5,i5,i5,i5,F12.8)') m,i,j,k,l, Mbar%val(i,j,k,l)
              endif
              m = m + 1
 
            end do
          end do
        end do
      end do
      print *, "";      print *, "occ-occ-occ-occ block of two-el density"
      m = 1
      do l=1,no
        do k=1,no
          do j=1,no
            do i=1,no
 
              if (ABS(d2el%val(i,j,k,l))>0.00000001)then 
               write (*,'(i5,i5,i5,i5,i5,F12.8)') m,i,j,k,l, d2el%val(i,j,k,l)
              endif
              m = m + 1
 
            end do
          end do
        end do
      end do
 
     endif
 
     !2. occ-occ-occ-virt block
     !*************************
     !d^sd_ijka = 2*delta_ij * Sum_cl z^c_l * tbar^ac_kl -
     !          -   delta_jk * Sum_cl z^c_l * tbar^ac_il -
     !          -              Sum_c  z^c_j * tbar^ac_ki =
     !
     !          = 2*delta_ij * E_ak -
     !          -   delta_jk * E_ai -
     !          -              Sum_c  z^c_j * tbar^ac_ki
 
 
     !container for amps vooo shape
     temp1 = array4_init([no,no,nv,no])
 
     !temp1 = Sum_c z^c_j * tbar^ac_ki
     temp = array2_init([nv,no],m1_arr%elm2)
     call array4_contract1(tbar,temp,temp1,.true.)
 
     !clean up!
     call array2_free(temp)
 
     !construct d^sd_ijka
     do i=1,no
       do j=1,no
         do k=1,no
           do a=no+1,nb
 
             d2el%val(i,j,k,a) = d2el%val(i,j,k,a) - temp1%val(i,j,a-no,k)
 
             if(i==j)then
               d2el%val(i,j,k,a) = d2el%val(i,j,k,a) + 2.0 * Eai%val(a-no,k)
             endif
 
             if(k==j)then
               d2el%val(i,j,k,a) = d2el%val(i,j,k,a) - Eai%val(a-no,i)
             endif
 
           end do
         end do
       end do
     end do

     !clean up!
     call array4_free(temp1)

!     !print debug info if needed
!     if( DECinfo%PL>2 )then
!     print *, "";!     print *, "occ-occ-occ-virt block of two-el density"
!     m = 1
!     do a=no+1,nb
!       do k=1,no
!         do j=1,no
!           do i=1,no
!
!             if (ABS(d2el%val(i,j,k,a))>0.00000001)then 
!              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,i,j,k,a, d2el%val(i,j,k,a)
!             endif
!             m = m + 1
!
!           end do
!         end do
!       end do
!     end do
!     endif
!

     !3. occ-occ-virt-virt block
     !**************************
     !d^sd_ijab = 2*delta_ij * Sum_ckl z^ca_kl * t^cb_kl -
     !          -              Sum_ck  z^ca_kj * t^cb_ki -
     !          -              Sum_ck  z^ca_jk * t^cb_ik =
     !
     !          = 2*delta_ij * Y_ba -
     !          -              Sum_ck  z^ca_kj * t^cb_ki -
     !          -              Sum_ck  z^ca_jk * t^cb_ik
 
 
     ! use containers for amps vovo shape
     temp1 = array4_init([nv,no,nv,no])
     temp2 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp1%val)
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp2%val)
 
     !form M1bar intermediate
     ! M1bar_ajbi = Sum_ck z^ca_kj * t^cb_ki
     M1bar = array4_init([nv,no,nv,no])
     call array4_contract2(temp1,temp2,M1bar)
     
     !clean up!
     call array4_free(temp1)
     call array4_free(temp2)
 
     !sort amps(cjak) -> ckaj
     temp1 = array4_init([no,nv,no,nv])
     temp2 = array4_init([no,nv,no,nv])
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,4,2,3],0.0E0_realk,temp1%val)
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[1,4,2,3],0.0E0_realk,temp2%val)
 
     !form M2bar intermediate
     ! M2bar_jaib = Sum_ck z^ca_jk * t^cb_ik
     M2bar = array4_init([no,nv,no,nv])
     call array4_contract2(temp1,temp2,M2bar)
     
     !clean up!
     call array4_free(temp1)
     call array4_free(temp2)
 
     !construct d^sd_ijab
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do b=no+1,nb
 
             d2el%val(i,j,a,b) = d2el%val(i,j,a,b) - M1bar%val(a-no,j,b-no,i)
             d2el%val(i,j,a,b) = d2el%val(i,j,a,b) - M2bar%val(j,a-no,i,b-no)
 
             if(i==j)then
               d2el%val(i,j,a,b) = d2el%val(i,j,a,b) + 2.0*Yba%val(b-no,a-no)
             endif
 
           end do
         end do
       end do
     end do
     
     !clean up!
     call array4_free(M1bar)
     call array4_free(M2bar)
 
     !print debug info if needed
     if( DECinfo%PL>2 )then
     print *, "";     print *, "occ-occ-virt-virt block of two-el density"
     m = 1
     do b=no+1,nb
       do a=no+1,nb
         do i=1,no
           do j=1,no
 
             if (ABS(d2el%val(i,j,a,b))>0.00000001)then 
              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,i,j,a,b, d2el%val(i,j,a,b)
             endif
             m = m + 1
 
           end do
         end do
       end do
     end do
     endif
     !clean up!

 
!     !4. occ-virt-occ-virt block
!     !**************************
!     !d^sd_iajb = 2*tbar^ab_ij + Sum_ckem    t^ea_mi *    t^ec_mk * tbar^cb_kj +
!     !          +                Sum_ckem    t^ab_mk *    t^ce_ij *    t^ce_mk -
!     !          -                Sum_ckem tbar^eb_ij *    t^ca_km *    t^ce_km -
!     !          -                Sum_ckem tbar^ab_mj *    t^ce_ki *    t^ce_km -
!     !          -                Sum_ckem tbar^ab_im *    t^ce_kj *    t^ce_km -
!     !          -                Sum_ckem tbar^ac_ij *    t^eb_km *    t^ec_km -
!     !          -                Sum_ckem    t^ca_kj * tbar^eb_mi *    t^ce_km +
!     !          +                Sum_ckem    t^ac_kj *    t^eb_mi *    t^ce_km +
!     !          +                Sum_ckem    t^ac_kj *    t^eb_im *    t^ce_mk
!     !
!     !          = 2*tbar^ab_ij + Sum_ckem tbar^cb_kj * tbar^ea_mi *    t^ec_mk +
!     !                         + Sum_km t^ab_mk * Mbar_imjk - Sum_e tbar^eb_ij * Y_ae -
!     !                         - Sum_m tbar^ab_mj * X_im - Sum_m tbar^ab_im * X_jm -
!     !                         - Sum_c tbar^ac_ij * Y_bc - Sum_ckem t^ca_kj * tbar^eb_mi * t^ce_mk +
!     !                         + Sum_ckem t^ac_kj * t^eb_mi * t^ce_km +
!     !                         + Sum_ckem t^ac_kj * t^eb_im * t^ce_mk

     !Sum_ckem tbar^cb_kj * tbar^ea_mi *    z^ec_mk
     !temporary containers
     temp1 = array4_init([nv,no,nv,no])
     M1bar = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp1%val)
     call array4_contract2(tbar,temp1,M1bar)
     !Sum_ck tbar^cb_kj * M1bar_aick
     temp2 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,M1bar%val,nv,no,nv,no,[3,4,1,2],0.0E0_realk,temp2%val)
     M1bar = array4_init([nv,no,nv,no])
     call array4_contract2(tbar,temp2,M1bar)

     !Sum_km t^ab_mk * Mbar_imjk
     !temporary containers
     temp1 = array4_init([no,no,no,no])
     temp2 = array4_init([no,no,nv,nv])
     M2bar = array4_init([nv,nv,no,no])
     call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,2,1,3],0.0E0_realk,temp1%val)

     !call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,3,1,2],0.0E0_realk,temp1%val)

     !call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,1,2,3],0.0E0_realk,temp1%val)
     !call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,1,3,2],0.0E0_realk,temp1%val)
     !call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,2,1,3],0.0E0_realk,temp1%val)
     !call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,2,3,1],0.0E0_realk,temp1%val)
     !call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,3,1,2],0.0E0_realk,temp1%val)
     !call array_reorder_4d(1.0E0_realk,Mbar%val,no,no,no,no,[4,3,2,1],0.0E0_realk,temp1%val)

     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[4,2,1,3],0.0E0_realk,temp2%val)
     call array4_contract2(temp2,temp1,M2bar)

     !temp container
     Yae = array4_init([nv,no,nv,no])
     !Yae = Sum_e tbar^eb_ij * Y_ae
     call array4_contract1(tbar,Yba,Yae,.true.)

     !temp container
     Xim = array4_init([no,nv,nv,no])
     !tbar(ambj) -> tbar(mabj)
     temp1 = array4_init([no,nv,nv,no])
     call array_reorder_4d(1.0E0_realk,tbar%val,nv,no,nv,no,[2,1,3,4],0.0E0_realk,temp1%val)
     !Xim = Sum_m tbar^ab_mj * X_im
     call array4_contract1(temp1,Xij,Xim,.true.)

     !temp container
     Xjm = array4_init([no,nv,no,nv])
     !tbar(aibm) -> tbar(maib)
     temp1 = array4_init([no,nv,no,nv])
     call array_reorder_4d(1.0E0_realk,tbar%val,nv,no,nv,no,[4,1,2,3],0.0E0_realk,temp1%val)
     !Xim = Sum_m tbar^ab_mj * X_im
     call array4_contract1(temp1,Xij,Xjm,.true.)

     !temp container
     Ybc = array4_init([nv,nv,no,no])
     !tbar(aicj) -> tbar(caij)
     temp1 = array4_init([nv,nv,no,no])
     call array_reorder_4d(1.0E0_realk,tbar%val,nv,no,nv,no,[3,1,2,4],0.0E0_realk,temp1%val)
     !Xim = Sum_m tbar^ab_mj * X_im
     call array4_contract1(temp1,Yba,Ybc,.true.)

     !Sum_ckem t^ca_kj * tbar^eb_mi * t^ce_mk
     temp1 = array4_init([nv,no,nv,no])
     temp2 = array4_init([nv,no,nv,no])
     !t(cmek) -> temp2(emck)
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[3,2,1,4],0.0E0_realk,temp2%val)
     !Sum_em tbar^eb_mi * temp2^ec_mk
     call array4_contract2(tbar,temp2,temp1)
     !Sum_ck t^ca_kj * temp1_bick
     temp2 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp2%val)
     temp3 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,temp1%val,nv,no,nv,no,[3,4,1,2],0.0E0_realk,temp3%val)
     M3bar = array4_init([nv,no,nv,no])
     call array4_contract2(temp2,temp3,M3bar)

     !Sum_ckem t^ac_kj * t^eb_mi * t^ce_km
     temp1 = array4_init([nv,no,nv,no])
     temp2 = array4_init([nv,no,nv,no])
     temp3 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp3%val)
     !t(ckem) -> t(emck)
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[3,4,1,2],0.0E0_realk,temp2%val)
     !Sum_em t^eb_mi * t^ec_mk
     call array4_contract2(temp3,temp2,temp1)
     !Sum_ck t^ac_kj * temp1_bick
     temp3 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,temp1%val,nv,no,nv,no,[3,4,1,2],0.0E0_realk,temp3%val)
     temp2 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[3,2,1,4],0.0E0_realk,temp2%val)
     M4bar = array4_init([nv,no,nv,no])
     call array4_contract2(temp2,temp3,M4bar)

     !+ Sum_ckem t^ac_kj * t^eb_im * t^ce_mk
     temp1 = array4_init([nv,no,nv,no])
     temp2 = array4_init([nv,no,nv,no])
     temp3 = array4_init([nv,no,nv,no])
     !t^eb_im (eibm) -> (embi)
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[1,4,3,2],0.0E0_realk,temp3%val)
     !t^ce_mk (cmek) -> (emck)
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[3,2,1,4],0.0E0_realk,temp2%val)
     !Sum_em t^eb_im * t^ce_mk
     call array4_contract2(temp3,temp2,temp1)
     !Sum_ck t^ac_kj * temp1_bick
     temp3 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,temp1%val,nv,no,nv,no,[3,4,1,2],0.0E0_realk,temp3%val)
     temp2 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[3,2,1,4],0.0E0_realk,temp2%val)
     M5bar = array4_init([nv,no,nv,no])
     call array4_contract2(temp2,temp3,M5bar)

     !construct d^sd_iajb
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do b=no+1,nb

             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) + 2 * tbar%val(a-no,i,b-no,j)

             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) + M1bar%val(b-no,j,a-no,i)
             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) + M2bar%val(a-no,b-no,i,j)

             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) - Yae%val(a-no,i,b-no,j)

             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) - Xim%val(i,a-no,b-no,j)
             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) - Xjm%val(j,a-no,i,b-no)

             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) - Ybc%val(b-no,a-no,i,j)

             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) - M3bar%val(a-no,j,b-no,i)
             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) + M4bar%val(a-no,j,b-no,i)
             d2el%val(i,a,j,b) = d2el%val(i,a,j,b) + M5bar%val(a-no,j,b-no,i)

           end do
         end do
       end do
     end do

!     !print debug info if needed
!     if( DECinfo%PL>2 )then
!     print *, "M1bar intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, M1bar%val(a-no,i,b-no,j)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "M2bar intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, M2bar%val(b-no,a-no,i,j)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "M3bar intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, M3bar%val(b-no,j,a-no,i)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "M4bar intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, M4bar%val(b-no,j,a-no,i)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "M5bar intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, M5bar%val(b-no,j,a-no,i)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "Yae intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, Yae%val(a-no,i,b-no,j)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "Xim intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, Xim%val(i,a-no,b-no,j)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "Xjm intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, Xjm%val(i,a-no,j,b-no)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "Ybc intermediate"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, Ybc%val(a-no,b-no,i,j)
!
!           end do
!         end do
!       end do
!     end do
!
!     print *, "numbers from occ-virt-occ-virt block of two-el density"
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             print *, d2el%val(i,a,j,b)
!
!           end do
!         end do
!       end do
!     end do
!     !ONLY if needed (huge array)
!     !   call print_4_dimensional_array_on_screen([nb,nb,nb,nb],d2el%val,"Two-electron density N4")
!     endif
!     !clean up!
!

     !5. virt-occ-virt-occ block
     !**************************
     !d^sd_aibj = z^ab_ij
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do b=no+1,nb

             d2el%val(a,i,b,j) = d2el%val(a,i,b,j) + m2_arr%elm4(a-no,i,b-no,j)

           end do
         end do
       end do
     end do

     !print debug info if needed
     if( DECinfo%PL>2 )then
     print *, "";     print *, "virt-occ-virt-occ block of two-el density"
     m = 1
     do j=1,no
       do b=no+1,nb
         do i=1,no
           do a=no+1,nb

             if (ABS(d2el%val(a,i,b,j))>0.00000001)then 
              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,i,b,j, d2el%val(a,i,b,j)
             endif
             m = m + 1

           end do
         end do
       end do
     end do
     endif
     !clean up!


     !6. virt-occ-occ-occ block
     !**************************
     !d^sd_aijk = 2*delta_jk * z^a_i - delta_ij * z^a_k
     do i=1,no
       do j=1,no
         do k=1,no
           do a=no+1,nb
 
             if(i==j)then
               d2el%val(a,i,j,k) = d2el%val(a,i,j,k) - m1_arr%elm2(a-no,k)
!               d2el%val(i,j,a,k) = d2el%val(i,j,a,k) - m1_arr%elm2(a-no,k)
             endif
 
             if(j==k)then
               d2el%val(a,i,j,k) = d2el%val(a,i,j,k) + 2.0*m1_arr%elm2(a-no,i)
!               d2el%val(i,j,a,k) = d2el%val(i,j,a,k) - m1_arr%elm2(a-no,i)
             endif
 
           end do
         end do
       end do
     end do
 
     !print debug info if needed
     if( DECinfo%PL>2 )then
     print *, "";     print *, "virt-occ-occ-occ block of two-el density"
     m = 1
     do k=1,no
       do j=1,no
         do i=1,no
           do a=no+1,nb
 
             if (ABS(d2el%val(a,i,j,k))>0.00000001)then 
              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,i,j,k, d2el%val(a,i,j,k)
             endif
             m = m + 1
 
           end do
         end do
       end do
     end do
     endif
 

     !7. virt-occ-occ-virt block
     !**************************
     !d^sd_aijb = Sum_dk z^ad_ik * tbar^db_kj - delta_ij * Sum_ckl z^ca_kl * t^cb_kl =
     !          = Sum_dk z^ad_ik * tbar^db_kj - delta_ij * Yba
 
 
     !Sum_dk z^ad_ik * tbar^db_kj
     !temporary containers
     temp1 = array4_init([nv,no,nv,no])
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[3,4,1,2],0.0E0_realk,temp1%val)
     temp2 = array4_init([nv,no,nv,no])
     call array4_contract2(temp1,tbar,temp2)
 
     !clean up!
     call array4_free(temp1)
 
     !construct d^sd_aijb
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do b=no+1,nb
 
             d2el%val(a,i,j,b) = d2el%val(a,i,j,b) + temp2%val(a-no,i,b-no,j)
 
             if(i==j)then
              ! d2el%val(a,i,j,b) = d2el%val(a,i,j,b) - Yba%val(b-no,a-no)
               d2el%val(a,i,j,b) = d2el%val(a,i,j,b) - Yba%val(a-no,b-no)
             endif
 
           end do
         end do
       end do
     end do
 
     !clean up!
     call array4_free(temp2)
 
     !print debug info if needed
     if( DECinfo%PL>2 )then
     print *, "";     print *, "virt-occ-occ-virt block of two-el density"
     m = 1
     do b=no+1,nb
       do j=1,no
         do i=1,no
           do a=no+1,nb
 
             if (ABS(d2el%val(a,i,j,b))>0.00000001)then 
              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,i,j,b, d2el%val(a,i,j,b)
             endif
             m = m + 1
 
           end do
         end do
       end do
     end do
     endif
 
     !8. virt-virt-occ-virt block
     !***************************
     !d^sd_abic = Sum_m z^a_m * tbar^bc_mi
 
     !container for amps vooo shape
     temp  = array2_init([nv,no],m1_arr%elm2)
     call array2_transpose(temp)
     temp1 = array4_init([no,nv,nv,no])
     temp2 = array4_init([nv,nv,nv,no])
 
     call array_reorder_4d(1.0E0_realk,tbar%val,nv,no,nv,no,[2,1,3,4],0.0E0_realk,temp1%val)
 
     !temp1 = Sum_m z^a_m * tbar^bc_mi
     call array4_contract1(temp1,temp,temp2,.true.)
 
     !clean up!
     call array4_free(temp1)
     call array2_free(temp)
 
     !construct d^sd_abic
     do i=1,no
       do c=no+1,nb
         do a=no+1,nb
           do b=no+1,nb
 
             d2el%val(a,b,i,c) = d2el%val(a,b,i,c) + temp2%val(a-no,b-no,c-no,i)
 
           end do
         end do
       end do
     end do
 
     !clean up!
     call array4_free(temp2)
 
     !print debug info if needed
     if( DECinfo%PL>2 )then
     print *, "";     print *, "virt-virt-occ-virt block of two-el density"
     m = 1
     do c=no+1,nb
       do i=1,no
         do b=no+1,nb
           do a=no+1,nb
 
             if (ABS(d2el%val(a,b,i,c))>0.00000001)then 
              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,b,i,c, d2el%val(a,b,i,c)
             endif
             m = m + 1
 
           end do
         end do
       end do
     end do
     endif
 

     !9. virt-virt-virt-virt block
     !****************************
     !d^sd_abcd = Sum_mn z^ac_mn * t^bd_mn
 
     !temporary containers
     temp1 = array4_init([no,no,nv,nv])
     temp2 = array4_init([no,no,nv,nv])
     temp3 = array4_init([nv,nv,nv,nv])
     call array_reorder_4d(1.0E0_realk,m2_arr%elm4,nv,no,nv,no,[4,2,3,1],0.0E0_realk,temp1%val)
     call array_reorder_4d(1.0E0_realk,t2_arr%elm4,nv,no,nv,no,[4,2,3,1],0.0E0_realk,temp2%val)
     call array4_contract2(temp1,temp2,temp3)
 
     !clean up!
     call array4_free(temp1)
     call array4_free(temp2)
 
     !construct d^sd_abcd
     do a=no+1,nb
       do b=no+1,nb
         do c=no+1,nb
           do d=no+1,nb
 
               d2el%val(a,b,c,d) = d2el%val(a,b,c,d) + temp3%val(a-no,c-no,b-no,d-no)
 
           end do
         end do
       end do
     end do
     
     !clean up!
     call array4_free(temp3)
 
     !print debug info if needed
     if( DECinfo%PL>2 )then
     print *, "";     print *, "virt-virt-virt-virt block of two-el density"
     m = 1
     do d=no+1,nb
       do c=no+1,nb
         do b=no+1,nb
           do a=no+1,nb
 
             if (ABS(d2el%val(a,b,c,d))>0.00000001)then 
              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,b,c,d, d2el%val(a,b,c,d)
             endif
             m = m + 1
 
           end do
         end do
       end do
      end do
      endif
 
     !10. redundant blocks
     !****************************
     ! ovoo = ooov
     ! vvoo = oovv
     ! oovo = vooo
     ! ovvo = voov
     ! ovvv = vvov

     !construct
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do k=1,no
 
             d2el%val(i,a,j,k) = d2el%val(i,a,j,k) + d2el%val(j,k,i,a)
 
           end do
         end do
       end do
     end do
 
     !construct
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do b=no+1,nb
 
             d2el%val(a,b,i,j) = d2el%val(a,b,i,j) + d2el%val(i,j,a,b)
 
           end do
         end do
       end do
     end do
 
     !construct
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do k=1,no
 
             d2el%val(i,j,a,k) = d2el%val(i,j,a,k) + d2el%val(a,k,i,j)
 
           end do
         end do
       end do
     end do
 
     !construct
     do i=1,no
       do j=1,no
         do a=no+1,nb
           do b=no+1,nb
 
             d2el%val(i,a,b,j) = d2el%val(i,a,b,j) + d2el%val(b,j,i,a)
 
           end do
         end do
       end do
     end do
 
     !construct
     do i=1,no
       do c=no+1,nb
         do a=no+1,nb
           do b=no+1,nb
 
             d2el%val(i,a,b,c) = d2el%val(i,a,b,c) + d2el%val(b,c,i,a)
 
           end do
         end do
       end do
     end do
 
 !     !clean up!
 !     call array2_free(Xij)
 !     call array2_free(Yba)
 !     call array2_free(Eai)
 !     call array4_free(tbar)
 !     call array4_free(Mbar)
     !***************************************************************
     ! Construct the RHS from the one- and two-el density
     !***************************************************************

     !density-based form
     !************************
     ! eta^ccsd_ai = Sum_t(1 - P_pq)
     !                      Sum_tj  (D^SD_tq   + D^SD_qt) * h_pt
     !                          Sum_tj  (D^SD_tq   + D^SD_qt) * h_pt
     !                              +Sum_tbj(d^SD_tqrs + d^SD_qtrs)*g_qtrs]

     !T1-transformed integrals first
     xocc = array2_init([nb,no])
     yocc = array2_init([nb,no])
     xvirt = array2_init([nb,nv])
     yvirt = array2_init([nb,nv])

     Co = array2_init([nb,no],MyFragment%Co)
     Cv = array2_init([nb,nv],MyFragment%Cv)
     Co2 = array2_init([nb,no],MyFragment%Co)
     Cv2 = array2_init([nb,nv],MyFragment%Cv)

!     Co%val => MyFragment%Co
!     Cv%val => MyFragment%Cv
!     Co2%val => MyFragment%Co
!     Cv2%val => MyFragment%Cv

     t1_trans = array2_init([nv,no])

     !DEBUG
     print *, "Debug print t1_trans: "
     call print_matrix(t1_trans)

     t1_trans = array2_init([nv,no],m1_arr%elm2)

     !DEBUG
     print *, "Debug print t1_trans: "
     call print_matrix(t1_trans)
          print *, "Debug print Co: "
     call print_matrix(Co)
          print *, "Debug print Cv: "
     call print_matrix(Cv)
          print *, "Debug print Co2: "
     call print_matrix(Co2)
          print *, "Debug print Cv2: "
     call print_matrix(Cv2)

     !put in temp amplitudes
!     call getT1transformation(t1_trans,xocc,xvirt,yocc,yvirt,&
!                              &Co,Cv,Co2,Cv2)

     !DEBUG
     print *, "Debug print: "
          print *, "Debug print xocc: "
     call print_matrix(xocc)
          print *, "Debug print xvirt: "
     call print_matrix(xvirt)
          print *, "Debug print yocc: "
     call print_matrix(yocc)
          print *, "Debug print yvirt: "
     call print_matrix(yvirt)

     ! Get one-electron contribution
     call mat_init(h1,nb,nb)
     call mat_zero(h1)

     call mat_init(tmp1,nb,nb)
     call mat_init(tmp2,nb,nb)
     call mat_zero(tmp1)
     call mat_zero(tmp2)

    ! call II_get_h1(DECinfo%output,DECinfo%output,MyFragment%mylsitem%setting,h1)
     
    call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyFragment%MyLsItem%setting,&
         & h1%elms,nb,nb,AORdefault,AORdefault)
     !DEBUG
     write (*,*)'The one-electron integrals:'
     print *, "";     m = 1
     do a=1,nb*nb
       if (ABS(h1%elms(a))>0.00000001)then 
         write (*,'(i5,i5,F12.8)') m,a, h1%elms(a)
       endif
       m = m + 1
     end do

     hoo = array2_init([no,no])
     hov = array2_init([no,nv])
     hvo = array2_init([nv,no])
     hvv = array2_init([nv,nv])

     call mem_alloc(w1,max(max(max(max(o2v2,ov3),v4),o2*v2),o4))

     !DEBUG: get MO coeff to compare with DALTON
     call mem_alloc( focc,     no     )
     call mem_alloc( fvirt,    nv     )
     Co_d  = array2_init([nb,no])
     Cv_d  = array2_init([nb,nv])
     Uocc  = array2_init([no,no])
     Uvirt = array2_init([nv,nv])

     write (*,*)'The MO coefficients:'
     call print_matrix(Co_d)
     call print_matrix(Cv_d)

     !DEBUG: get MO coeff to compare with DALTON
     !call get_canonical_integral_transformation_matrices(no,nv,nb,&
     !                           &MyFragment%ppfock,MyFragment%qqfock,Co%val,Cv%val,&                         
     !                           &Co_d%val,Cv_d%val,Uocc%val,Uvirt%val,focc,fvirt) 

     Co_dl  = array2_init([nb,no])
     Cv_dl  = array2_init([nb,nv])
     Co_dr  = array2_init([nb,no])
     Cv_dr  = array2_init([nb,nv])

     Co_dl%val(1,1) = -0.97995428
     Co_dl%val(2,1) = -0.01007583 
     Co_dl%val(3,1) =  0.01997469 
     Co_dl%val(4,1) =  0.0
     Co_dl%val(5,1) =  0.0
     Co_dl%val(6,1) = -0.05636319 
     Co_dl%val(1,2) =  0.29102213 
     Co_dl%val(2,2) = -0.37799340 
     Co_dl%val(3,2) = -0.40019444 
     Co_dl%val(4,2) =  0.0
     Co_dl%val(5,2) =  0.0
     Co_dl%val(6,2) = -0.54463061 

     Co_dr%val(1,1) = -0.97988625
     Co_dr%val(2,1) = -0.01105724
     Co_dr%val(3,1) =  0.02006786
     Co_dr%val(4,1) =  0.0
     Co_dr%val(5,1) =  0.0
     Co_dr%val(6,1) = -0.05575207
     Co_dr%val(1,2) =  0.29649270
     Co_dr%val(2,2) = -0.40489724
     Co_dr%val(3,2) = -0.38299817
     Co_dr%val(4,2) =  0.0
     Co_dr%val(5,2) =  0.0
     Co_dr%val(6,2) = -0.53909859

     Cv_dl%val(1,1) =  0.17207035
     Cv_dl%val(2,1) = -0.84100329
     Cv_dl%val(3,1) =  0.58553533
     Cv_dl%val(4,1) =  0.0
     Cv_dl%val(5,1) =  0.0
     Cv_dl%val(6,1) =  0.17071817
     Cv_dl%val(1,2) =  0.0
     Cv_dl%val(2,2) =  0.0
     Cv_dl%val(3,2) =  0.0
     Cv_dl%val(4,2) =  0.99144175
     Cv_dl%val(5,2) =  0.13054983
     Cv_dl%val(6,2) =  0.0
     Cv_dl%val(1,3) =  0.0
     Cv_dl%val(2,3) =  0.0
     Cv_dl%val(3,3) =  0.0
     Cv_dl%val(4,3) =  0.13054983
     Cv_dl%val(5,3) = -0.99144175
     Cv_dl%val(6,3) =  0.0
     Cv_dl%val(1,4) = -0.22990007
     Cv_dl%val(2,4) = -0.94506670
     Cv_dl%val(3,4) = -1.11052009
     Cv_dl%val(4,4) =  0.0
     Cv_dl%val(5,4) =  0.0
     Cv_dl%val(6,4) =  1.50887818

     Cv_dr%val(1,1) =  0.18031296
     Cv_dr%val(2,1) = -0.85272430
     Cv_dr%val(3,1) =  0.57315014
     Cv_dr%val(4,1) =  0.0
     Cv_dr%val(5,1) =  0.0
     Cv_dr%val(6,1) =  0.15379688
     Cv_dr%val(1,2) =  0.0
     Cv_dr%val(2,2) =  0.0
     Cv_dr%val(3,2) =  0.0
     Cv_dr%val(4,2) =  0.99144175
     Cv_dr%val(5,2) =  0.13054983
     Cv_dr%val(6,2) =  0.0
     Cv_dr%val(1,3) =  0.0
     Cv_dr%val(2,3) =  0.0
     Cv_dr%val(3,3) =  0.0
     Cv_dr%val(4,3) =  0.13054983
     Cv_dr%val(5,3) = -0.99144175
     Cv_dr%val(6,3) =  0.0
     Cv_dr%val(1,4) = -0.23007018
     Cv_dr%val(2,4) = -0.94526196
     Cv_dr%val(3,4) = -1.11071688
     Cv_dr%val(4,4) =  0.0
     Cv_dr%val(5,4) =  0.0
     Cv_dr%val(6,4) =  1.50858327

     !DEBUG
     print *, "";     print *,MyFragment%ppfock(1,1)
     print *, "";     write (*,*)'The MO coefficients:'
     m = 1
     do i=1,no
     do a=1,nb
       if (ABS(Co_d%val(a,i))>0.00000001)then 
         write (*,'(i5,F12.8)') m, Co_d%val(a,i)
       endif
       m = m + 1
     end do
     end do
     call print_matrix(Co)
     call print_matrix(Co_dl)
     call print_matrix(Co_dr)
     call print_matrix(Cv)
     call print_matrix(Cv_dl)
     call print_matrix(Cv_dr)

     temp  = array2_init([no,nb])
     !DEBUG Transform integrals to MO basis
     ! -> hoo
     call dgemm('t','n',no,nb,nb,1.0E0_realk,Co_dl%val,nb,h1%elms,nb,0.0E0_realk,w1,no)
     call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,Co_dr%val,nb,0.0E0_realk,hoo%val,no)
     ! -> hov
     call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,Cv_dr%val,nb,0.0E0_realk,hov%val,no)
     ! -> hvo
     call dgemm('t','n',nv,nb,nb,1.0E0_realk,Cv_dl%val,nb,h1%elms,nb,0.0E0_realk,w1,nv)
     call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,Co_dr%val,nb,0.0E0_realk,hvo%val,nv)
     ! -> hvv
     call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,Cv_dr%val,nb,0.0E0_realk,hvv%val,nv)

!     !Transform integrals to T1-transformed basis
!     ! -> hoo
!     call dgemm('t','n',no,nb,nb,1.0E0_realk,xocc%val,nb,h1%elms,nb,0.0E0_realk,w1,no)
!     call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,yocc%val,nb,0.0E0_realk,hoo%val,no)
!     ! -> hov
!     call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,yvirt%val,nb,0.0E0_realk,hov%val,no)
!     ! -> hvo
!     call dgemm('t','n',nv,nb,nb,1.0E0_realk,xvirt%val,nb,h1%elms,nb,0.0E0_realk,w1,nv)
!     call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,yocc%val,nb,0.0E0_realk,hvo%val,nv)
!     ! -> hvv
!     call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,yvirt%val,nb,0.0E0_realk,hvv%val,nv)

     !the T1-transformation matricies
     write (*,*)'The one-electron integrals:'
     call print_matrix(hoo)
     call print_matrix(hov)
     call print_matrix(hvo)
     call print_matrix(hvv)

     !*********************************************
     !add ONE-electron contribution to the RHS
     ! of orbital-rotation multiplier equations
     !*********************************************

     !occ block
     eta_arr = array2_init([nv,no])
!     eta_arr2 = array2_init([nv,no])
     
!     write (*,*)'ETA ai initialized:'
!     call print_matrix(eta_arr)

     call eta_ai_construction(eta_arr,t1_arr,hoo,hov,hvo,hvv,Dsd,no,nb)

!     call array2_add_to(eta_arr,1.0E0_realk,eta_arr2) 
!
!     write (*,*)'ETA ai one-electron part:'
!     call print_matrix(eta_arr)
 
     !*********************************************
     !add TWO-electron contribution to the RHS
     ! of orbital-rotation multiplier equations
     !*********************************************

     ! get two-electron integrals in ao
     write(*,*) 'debug :: calculating AO integrals'
     ! Always in the DEBUG SOLVER: calculate full 4-dimensional AO integrals
     call get_full_eri(MyFragment%mylsitem,nb,gao)
     write(*,*) 'debug :: AO integrals done'

     call array4_read(gao)
!     if( DECinfo%PL>2 )then
!     print *, "";!     print *, "AO INTEGRALS"
!     m = 1
!     do d=1,nb
!       do b=1,nb
!         do i=1,nb
!           do a=1,nb
!
!             if (ABS(gao%val(a,i,b,d))>0.00000001)then 
!              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,i,b,d, gao%val(a,i,b,d)
!             endif
!             m = m + 1
!
!          end do
!        end do
!      end do
!    end do
!    endif

!    print *,* 

!     !d^sd_aibj = z^ab_ij
!     do d=1,nb
!     do i=1,no
!       do j=1,no
!         do a=no+1,nb
!           do b=no+1,nb
!
!             d2elaibj%val(a-no,i,b-no,d) = d2elaibj%val(a-no,i,b-no,d) +&
!                                           Co_dr%val(d,j)*m2_arr%elm4(a-no,i,b-no,j)
!
!             if (ABS(d2elaibj%val(a-no,i,b-no,d))>0.00000001)then 
!              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,i,b,j, Co_dr%val(d,j)
!              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,i,b,j, m2_arr%elm4(a-no,i,b-no,j)
!             endif
!
!           end do
!         end do
!       end do
!     end do
!     end do
!
!     !print debug info if needed
!     if( DECinfo%PL>2 )then
!     print *, "";!     print *, "ONLY virt-occ-virt-occ block of two-el density"
!     m = 1
!     do b=no+1,nb
!       do d=1,nb
!         do i=1,no
!           do a=no+1,nb
!
!             if (ABS(d2elaibj%val(a-no,i,b-no,d))>0.00000001)then 
!              write (*,'(i5,i5,i5,i5,i5,F12.8)') m,a,i,b,d, d2elaibj%val(a-no,i,b-no,d)
!             endif
!             m = m + 1
!
!          end do
!        end do
!      end do
!    end do
!    endif


     goooo = array4_init([no,no,no,no])
     goooo = get_gmo_simple(gao,Co_dl,Co_dr,Co_dl,Co_dr)
  
     govoo = array4_init([no,nv,no,no])
      govoo = get_gmo_simple(gao,Co_dl,Cv_dr,Co_dl,Co_dr)
  
     gvooo = array4_init([nv,no,no,no])
      gvooo = get_gmo_simple(gao,Cv_dl,Co_dr,Co_dl,Co_dr)
  
     gvvoo = array4_init([nv,nv,no,no])
      gvvoo = get_gmo_simple(gao,Cv_dl,Cv_dr,Co_dl,Co_dr)
      
     do i=1,no
       do j=1,no

         hoo = array2_init([no,no],goooo%val(:,:,i,j))
         hov = array2_init([no,nv],govoo%val(:,:,i,j))
         hvo = array2_init([nv,no],gvooo%val(:,:,i,j))
         hvv = array2_init([nv,nv],gvvoo%val(:,:,i,j))
      
         Dsd = array2_init([nb,nb], d2el%val(:,:,i,j))
      
!         if( DECinfo%PL>2 )then
!         print *, 'in the OO block'
!           print *, 'hoo'
!           call print_matrix(hoo)
!           print *, 'hov'
!           call print_matrix(hov)
!           print *, 'hvo'
!           call print_matrix(hvo)
!           print *, 'hvv'
!           call print_matrix(hvv)
!           print *, 'Dsd'
!           call print_matrix(Dsd)
!         endif
     
         call eta_ai_construction(eta_arr,t1_arr,hoo,hov,hvo,hvv,Dsd,no,nb)

       end do
     end do

     goovv = array4_init([no,no,nv,nv])
     goovv = get_gmo_simple(gao,Co_dl,Co_dr,Cv_dl,Cv_dr)
  
     govvv = array4_init([no,nv,nv,nv])
      govvv = get_gmo_simple(gao,Co_dl,Cv_dr,Cv_dl,Cv_dr)
  
     gvovv = array4_init([nv,no,nv,nv])
      gvovv = get_gmo_simple(gao,Cv_dl,Co_dr,Cv_dl,Cv_dr)
  
     gvvvv = array4_init([nv,nv,nv,nv])
      gvvvv = get_gmo_simple(gao,Cv_dl,Cv_dr,Cv_dl,Cv_dr)
      
     do a=no+1,nb
       do b=no+1,nb

         hoo = array2_init([no,no],goovv%val(:,:,a-no,b-no))
         hov = array2_init([no,nv],govvv%val(:,:,a-no,b-no))
         hvo = array2_init([nv,no],gvovv%val(:,:,a-no,b-no))
         hvv = array2_init([nv,nv],gvvvv%val(:,:,a-no,b-no))
      
         Dsd = array2_init([nb,nb], d2el%val(:,:,a,b))
     
 !        if( DECinfo%PL>2 )then
 !        print *, 'in the VV block'
 !          print *, 'hoo'
 !          call print_matrix(hoo)
 !          print *, 'hov'
 !          call print_matrix(hov)
 !          print *, 'hvo'
 !          call print_matrix(hvo)
 !          print *, 'hvv'
 !          call print_matrix(hvv)
 !          print *, 'Dsd'
 !          call print_matrix(Dsd)
 !        endif
         call eta_ai_construction(eta_arr,t1_arr,hoo,hov,hvo,hvv,Dsd,no,nb)

       end do
     end do

     goovo = array4_init([no,no,nv,no])
     goovo = get_gmo_simple(gao,Co_dl,Co_dr,Cv_dl,Co_dr)
  
     govvo = array4_init([no,nv,nv,no])
      govvo = get_gmo_simple(gao,Co_dl,Cv_dr,Cv_dl,Co_dr)
  
     gvovo = array4_init([nv,no,nv,no])
      gvovo = get_gmo_simple(gao,Cv_dl,Co_dr,Cv_dl,Co_dr)
  
     gvvvo = array4_init([nv,nv,nv,no])
      gvvvo = get_gmo_simple(gao,Cv_dl,Cv_dr,Cv_dl,Co_dr)
      
     do a=no+1,nb
       do i=1,no

         hoo = array2_init([no,no],goovo%val(:,:,a-no,i))
         hov = array2_init([no,nv],govvo%val(:,:,a-no,i))
         hvo = array2_init([nv,no],gvovo%val(:,:,a-no,i))
         hvv = array2_init([nv,nv],gvvvo%val(:,:,a-no,i))
      
         Dsd = array2_init([nb,nb], d2el%val(:,:,a,i))
     
 !        if( DECinfo%PL>2 )then
 !        print *, 'in the VO block'
 !          print *, 'hoo'
 !          call print_matrix(hoo)
 !          print *, 'hov'
 !          call print_matrix(hov)
 !          print *, 'hvo'
 !          call print_matrix(hvo)
 !          print *, 'hvv'
 !          call print_matrix(hvv)
 !          print *, 'Dsd'
 !          call print_matrix(Dsd)
 !        endif
         call eta_ai_construction(eta_arr,t1_arr,hoo,hov,hvo,hvv,Dsd,no,nb)

       end do
     end do

     gooov = array4_init([no,no,no,nv])
     gooov = get_gmo_simple(gao,Co_dl,Co_dr,Co_dl,Cv_dr)
  
     govov = array4_init([no,nv,no,nv])
      govov = get_gmo_simple(gao,Co_dl,Cv_dr,Co_dl,Cv_dr)
  
     gvoov = array4_init([nv,no,no,nv])
      gvoov = get_gmo_simple(gao,Cv_dl,Co_dr,Co_dl,Cv_dr)
  
     gvvov = array4_init([nv,nv,no,nv])
      gvvov = get_gmo_simple(gao,Cv_dl,Cv_dr,Co_dl,Cv_dr)
      
     do i=1,no
       do a=no+1,nb

         hoo = array2_init([no,no],gooov%val(:,:,i,a-no))
         hov = array2_init([no,nv],govov%val(:,:,i,a-no))
         hvo = array2_init([nv,no],gvoov%val(:,:,i,a-no))
         hvv = array2_init([nv,nv],gvvov%val(:,:,i,a-no))
      
         Dsd = array2_init([nb,nb], d2el%val(:,:,i,a))
     
 !        if( DECinfo%PL>2 )then
 !        print *, 'in the OV block'
 !          print *, 'hoo'
 !          call print_matrix(hoo)
 !          print *, 'hov'
 !          call print_matrix(hov)
 !          print *, 'hvo'
 !          call print_matrix(hvo)
 !          print *, 'hvv'
 !          call print_matrix(hvv)
 !          print *, 'Dsd'
 !          call print_matrix(Dsd)
 !        endif
         call eta_ai_construction(eta_arr,t1_arr,hoo,hov,hvo,hvv,Dsd,no,nb)

       end do
     end do

!     !
!     do a=no+1,nb
!       do i=1,no
!         do t=1,no
!           do r=1,no
!             do s=1,no
!
!               eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + govoo%val(t,a-no,r,s)*d2el%val(t,i,r,s)
!
!             end do
!           end do
!         end do
!       end do
!     end do
!
!     gvvvv = array4_init([nv,nv,nv,nv])
!     gvvvv = get_gmo_simple(gao,xvirt,yvirt,xvirt,yvirt)
!
!     !
!     do a=no+1,nb
!       do i=1,no
!         do t=no+1,nb
!           do r=no+1,nb
!             do s=no+1,nb
!
!               eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + gvvvv%val(t-no,a-no,r-no,s-no)*d2el%val(t,i,r,s)
!
!             end do
 !          end do
 !        end do
 !      end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector 5'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    goooo = array4_init([no,no,no,no])
 !    goooo = get_gmo_simple(gao,xocc,yocc,xocc,yocc)
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - goooo%val(i,t,r,s)*d2el%val(a,t,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !    end do
!
!     govvv = array4_init([no,nv,nv,nv])
!     govvv = get_gmo_simple(gao,xocc,yvirt,xvirt,yvirt)
!
!     !
!     do a=no+1,nb
!       do i=1,no
!         do t=no+1,nb
!           do r=no+1,nb
!             do s=no+1,nb
!
!               eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - govvv%val(i,t-no,r-no,s-no)*d2el%val(t,i,r,s)
!
!             end do
 !          end do
 !        end do
 !      end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - goooo%val(t,i,r,s)*d2el%val(t,a,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !    end do
 !
 !    gvovv = array4_init([nv,no,nv,nv])
 !    gvovv = get_gmo_simple(gao,xvirt,yocc,xvirt,yvirt)
!
!     !
!     do a=no+1,nb
!       do i=1,no
!         do t=no+1,nb
!           do r=no+1,nb
!             do s=no+1,nb
!
!               eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - gvovv%val(t-no,i,r-no,s-no)*d2el%val(t,a,r,s)
!
!             end do
!           end do
!         end do
!       end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    gvooo = array4_init([nv,no,no,no])
 !    gvooo = get_gmo_simple(gao,xvirt,yocc,xocc,yocc)
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + gvooo%val(a-no,t,r,s)*d2el%val(i,t,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !    end do
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !        do t=no+1,nb
 !          do r=no+1,nb
 !            do s=no+1,nb
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + gvvvv%val(a-no,t-no,r-no,s-no)*d2el%val(i,t,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_trans%val(a-no,j)*&
 !                                                          &goooo%val(t,i,r,s)*d2el%val(t,j,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !    end do
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !        do t=no+1,nb
 !          do r=no+1,nb
 !            do s=no+1,nb
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_trans%val(a-no,j)*&
 !                                                          &gvovv%val(t-no,i,r-no,s-no)*d2el%val(t,j,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_trans%val(a-no,j)*&
 !                                                          &goooo%val(j,t,r,s)*d2el%val(i,t,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !    end do
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !        do t=no+1,nb
 !          do r=no+1,nb
 !            do s=no+1,nb
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_trans%val(a-no,j)*&
 !                                                          &govvv%val(j,t-no,r-no,s-no)*d2el%val(i,t,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do b=no+1,nb
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_trans%val(b-no,i)*&
 !                                                          &govoo%val(t,b-no,r,s)*d2el%val(t,a,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !    end do
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do b=no+1,nb
 !        do t=no+1,nb
 !          do r=no+1,nb
 !            do s=no+1,nb
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_trans%val(b-no,i)*&
 !                                                          &gvvvv%val(t-no,b-no,r-no,s-no)*d2el%val(t,a,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !      do b=no+1,nb
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_trans%val(a-no,j)*&
 !                                                         &t1_trans%val(b-no,i)*&
 !                                                          &govoo%val(t,b-no,r,s)*d2el%val(t,j,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !      end do
 !    end do
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !      do b=no+1,nb
 !        do t=no+1,nb
 !          do r=no+1,nb
 !            do s=no+1,nb
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_trans%val(a-no,j)*&
 !                                                         &t1_trans%val(b-no,i)*&
 !                                                          &gvvvv%val(t-no,b-no,r-no,s-no)*d2el%val(t,j,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !      end do
 !    end do
 !
 !    !print debug info if needed
 !    if( DECinfo%PL>2 )then
 !       print *, 'RHS for z-vector'
 !       call print_matrix(eta_arr)
 !    endif
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !      do b=no+1,nb
 !        do t=1,no
 !          do r=1,no
 !            do s=1,no
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_trans%val(a-no,j)*&
 !                                                         &t1_trans%val(b-no,i)*&
 !                                                          &goooo%val(j,t,r,s)*d2el%val(b,t,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !      end do
 !    end do
 !
 !    !
 !    do a=no+1,nb
 !      do i=1,no
 !      do j=1,no
 !      do b=no+1,nb
 !        do t=no+1,nb
 !          do r=no+1,nb
 !            do s=no+1,nb
 !
 !              eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_trans%val(a-no,j)*&
 !                                                         &t1_trans%val(b-no,i)*&
 !                                                          &govvv%val(j,t-no,r-no,s-no)*d2el%val(b,t,r,s)
 !
 !            end do
 !          end do
 !        end do
 !      end do
 !      end do
 !      end do
 !    end do
 
     !print debug info if needed
     if( DECinfo%PL>2 )then
        print *, 'RHS for z-vector FINAL'
        call print_matrix(eta_arr)
     endif


  end subroutine z_vec_rhs_ccsd


  !> \construct eta_ai
  subroutine eta_ai_construction(eta_arr,t1_arr,hoo,hov,hvo,hvv,Dsd,no,nb)

    implicit none
    type(array2), intent(inout) :: eta_arr
    type(tensor), intent(in)    :: t1_arr
    type(array2), intent(in)    :: hoo
    type(array2), intent(in)    :: hov
    type(array2), intent(in)    :: hvo
    type(array2), intent(in)    :: hvv
    type(array2), intent(in)    :: Dsd
    integer                     :: i,j,a,b,t
    integer                     :: no,nb

      !print debug info if needed
!      if( DECinfo%PL>2 )then
!         print *, 'RHS for z-vector initial'
!         call print_matrix(eta_arr)
!      endif
 
      do a=no+1,nb
        do i=1,no
          do t=1,no
 !2nd+-
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - hoo%val(i,t)*Dsd%val(a,t)
 !5th+-
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - hoo%val(t,i)*Dsd%val(t,a)
 
          end do
        end do
      end do

      !occ-virt/virt-occ block
      do a=no+1,nb
        do i=1,no
          do t=no+1,nb
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - hov%val(i,t-no)*Dsd%val(a,t)
 !3rd!
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - hvo%val(t-no,i)*Dsd%val(t,a)
 
          end do
        end do
      end do
 
      !
      do a=no+1,nb
        do i=1,no
          do t=1,no
 !1st!
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + hov%val(t,a-no)*Dsd%val(t,i)
 !6th!
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + hvo%val(a-no,t)*Dsd%val(i,t)
 
          end do
        end do
      end do
      !occ-virt/virt-occ block
      do a=no+1,nb
        do i=1,no
          do t=no+1,nb
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + hvv%val(t-no,a-no)*Dsd%val(t,i)
 !4th!
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + hvv%val(a-no,t-no)*Dsd%val(i,t)
 
          end do
        end do
      end do
 
!      !print debug info if needed
!      if( DECinfo%PL>2 )then
!         print *, 'hvo'
!         call print_matrix(hvo)
!      endif
!      !print debug info if needed
!      if( DECinfo%PL>2 )then
!         print *, 'hov'
!         call print_matrix(hov)
!      endif
!      !print debug info if needed
!      if( DECinfo%PL>2 )then
!         print *, 'Dsd'
!         call print_matrix(Dsd)
!      endif
!      !print debug info if needed
!      if( DECinfo%PL>2 )then
!         print *, 'RHS for z-vector first'
!         call print_matrix(eta_arr)
!      endif
      !
      do a=no+1,nb
        do i=1,no
          do t=1,no
            do j=1,no
            
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_arr%elm2(a-no,j)*hoo%val(t,i)*Dsd%val(t,j)
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_arr%elm2(a-no,j)*hoo%val(j,t)*Dsd%val(i,t)
                
            end do
          end do
        end do
      end do
      
      !
      do a=no+1,nb
        do i=1,no
          do t=no+1,nb
            do j=1,no
            
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_arr%elm2(a-no,j)*hvo%val(t-no,i)*Dsd%val(t,j)
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_arr%elm2(a-no,j)*hov%val(j,t-no)*Dsd%val(i,t)
                
            end do
          end do
        end do
      end do
      !
      do a=no+1,nb
        do i=1,no
          do t=1,no
            do b=no+1,nb
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_arr%elm2(b-no,i)*hov%val(t,b-no)*Dsd%val(t,a)
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_arr%elm2(b-no,i)*hvo%val(a-no,t)*Dsd%val(b,t)
 
            end do
          end do
        end do
      end do
 
      !
      do a=no+1,nb
        do i=1,no
          do t=no+1,nb
            do b=no+1,nb
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_arr%elm2(b-no,i)*hvv%val(t-no,b-no)*Dsd%val(t,a)
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_arr%elm2(b-no,i)*hvv%val(a-no,t-no)*Dsd%val(b,t)
 
            end do
          end do
        end do
      end do
!      !print debug info if needed
!      if( DECinfo%PL>2 )then
!         print *, 'RHS for z-vector second'
!         call print_matrix(eta_arr)
!      endif
 
      !
      do a=no+1,nb
        do i=1,no
          do t=1,no
            do b=no+1,nb
              do j=1,no
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_arr%elm2(a-no,j)*t1_arr%elm2(b-no,i)&
                                                      &*hov%val(t,b-no)*Dsd%val(t,j)
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_arr%elm2(a-no,j)*t1_arr%elm2(b-no,i)&
                                                      &*hoo%val(j,t)*Dsd%val(b,t)
 
              end do
            end do
          end do
        end do
      end do
      !
      do a=no+1,nb
        do i=1,no
          do t=no+1,nb
            do b=no+1,nb
              do j=1,no
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) + t1_arr%elm2(a-no,j)*t1_arr%elm2(b-no,i)&
                                                      &*hvv%val(t-no,b-no)*Dsd%val(t,j)
 
                eta_arr%val(a-no,i) = eta_arr%val(a-no,i) - t1_arr%elm2(a-no,j)*t1_arr%elm2(b-no,i)&
                                                      &*hov%val(j,t-no)*Dsd%val(b,t)
 
              end do
            end do
          end do
        end do
      end do
 
!      !print debug info if needed
!      if( DECinfo%PL>2 )then
!         print *, 'RHS for z-vector third'
!         call print_matrix(eta_arr)
!      endif

    return

  end subroutine eta_ai_construction   

  !> \Print Matrix
  subroutine print_matrix(this)

    implicit none
    type(array2), intent(in) :: this
    integer :: i,j,k

    write(*, *) ''
    do i=1,this%dims(1)
      do j=1,this%dims(2)
        write(*,'(f16.10)',advance='no') this%val(i,j)
      end do
      write(*, *) ''
    end do

    return

  end subroutine print_matrix

  !> \Print Matrix
  subroutine print_matrix_real(this,dimi,dimj)

    implicit none
    real(realk), pointer :: this(:,:) 
    integer :: i,j,dimi,dimj

    write(*, *) ''
    do i=1,dimi
      do j=1,dimj
        write(*,'(f16.10)',advance='no') this(i,j)
      end do
      write(*, *) ''
    end do

    return

  end subroutine print_matrix_real



  !> \brief Print all elements of four-dimensional array to the screen
  subroutine print_4_dimensional_array_on_screen(dims,A,label)

    implicit none
    !> Dimensions of 4-dimensional array
    integer,dimension(4),intent(in) :: dims
    !> 4-dimensional array to be printed
    real(realk),intent(in) :: A(dims(1),dims(2),dims(3),dims(4))
    !> Label for array
    character(*), intent(in) :: label
    integer :: i,j,k,l

    write(*,*)
    write(*,*)
    write(*,*) '***********************************************************'
    write(*,*) '             ARRAY LABEL: ', label
    write(*,*) '***********************************************************'
    write(*,*)
    write(*,'(8X,a,8X,a,8X,a,8X,a,12X,a)') 'i','j','k','l', 'value'

    do i=1,dims(1)
       do j=1,dims(2)
          do k=1,dims(3)
             do l=1,dims(4)
                write(*,'(4i9,5X,g18.10)') i,j,k,l,A(i,j,k,l)
             end do
          end do
       end do
    end do

    write(*,*)
    write(*,*)


  end subroutine print_4_dimensional_array_on_screen

  !> \brief Calculate contributions to MP2 density from t2 amplitudes for pair fragment.
  !> See mp2dens structure for equations.
  subroutine pair_calculate_DEBUG(fragment1,fragment2,PairFragment,t2occ,t2virt,&
       & ThetaOCC, ThetaVIRT, VOOO,VOVV,dens)


    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
    !> Pair fragment
    type(decfrag),intent(inout) :: pairfragment
    !> t2 amplitudes, only for EOS orbitals using occupied partitioning, order:  (A,I,B,J)
    type(tensor),intent(in) :: t2occ
    !> t2 amplitudes, only for EOS orbitals using virtual partitioning, order:  (A,I,B,J)
    type(tensor),intent(in) :: t2virt
    !> Theta array, only for EOS orbitals using occupied partitioning, order:  (A,I,B,J)
    type(tensor),intent(inout) :: ThetaOCC
    !> Theta array, only for EOS orbitals using virtual partitioning, order:  (A,I,B,J)
    type(tensor),intent(in) :: ThetaVIRT
    !> (C I | J L) integrals stored as (C,I,J,L)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOOO
    !> (B K | A C) integrals stored as (B,K,A,C)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOVV
    !> MP2 density matrix contribution from given pair fragment
    type(mp2dens),intent(inout) :: dens
    logical, pointer :: dopair_occ(:,:), dopair_virt(:,:)
    integer :: i,j,a,b,c,k,d,l,noccEOS,nvirtEOS,noccAOS,nvirtAOS,natoms,nocctot
    logical :: something_wrong
    real(realk) :: tcpu, twall
    real(realk),pointer :: Y_tmp(:,:),X_tmp(:,:), Phivo_tmp(:,:), Phiov_tmp(:,:)

print *, "I AM HERE 1 "
    ! Init stuff
    noccEOS = t2occ%dims(2) ! number of occupied EOS orbitals
    noccAOS = t2virt%dims(2) ! number of occupied AOS orbitals
    nvirtEOS = t2virt%dims(1) ! number of virtual EOS orbitals
    nvirtAOS = t2occ%dims(1) ! number of virtual AOS orbitals
    ! number of occupied core+valence orbitals (only different from noccAOS for frozen approx)         
!    nocctot = VOOO%dims(4)
 
    ! Which "interaction pairs" to include for occ and virt space (avoid double counting)
    call mem_alloc(dopair_occ,noccEOS,noccEOS)
    call mem_alloc(dopair_virt,nvirtEOS,nvirtEOS)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)
 
    ! Just in case, zero matrices in dens structure
    dens%X = 0e0_realk
    dens%Y = 0e0_realk
    dens%Phivo = 0e0_realk
    dens%Phiov = 0e0_realk
 
 
 !   ! Sanity check
 !   something_wrong=.false.
 !   if(noccAOS/=dens%nocc) something_wrong=.true.
 !   if(nvirtAOS/=dens%nvirt) something_wrong=.true.
 !   if(nocctot/=dens%nocctot) something_wrong=.true.
 !   if(something_wrong) then
 !      write(DECinfo%output,*) 'dens%nocc', dens%nocc
 !      write(DECinfo%output,*) 'dens%nvirt', dens%nvirt
 !      write(DECinfo%output,*) 'dens%nocctot', dens%nocctot
 !      write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noccAOS
 !      write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvirtAOS
 !      call lsquit('single_calculate_mp2density: &
 !           & Something wrong with input dimensions!',DECinfo%output)
 !   end if
 !


    ! ******************************************************************************************
    !                Calculate contribution to virt-virt block of density matrix (Y)           *
    ! ******************************************************************************************

    ! Use occupied partitioning scheme
    ! Y_{ab} = sum_{cij} t_{ji}^{ca} * mult_{ji}^{cb}           [see the mp2dens structure]
    !        = 1/2* sum_{cij} t_{ji}^{ca} * Theta_{ji}^{cb}     [see construct_theta_array]
    ! -- where we only include the contributions if i and j belong to different fragments!
    call LSTIMER('START',tcpu,twall,DECinfo%output)
       print *, "I AM HERE 1 "

call mem_TurnONThread_Memory()
!OMP PARALLEL DEFAULT(shared) PRIVATE(Y_tmp,i,j,a,b,c)
call init_threadmemvar()

    call mem_alloc(Y_tmp,nvirtAOS,nvirtAOS)
    Y_tmp = 0E0_realk
 
!OMP DO SCHEDULE(dynamic,1)
       print *, "I AM HERE 2 "
 
    do i=1,noccEOS
       do j=1,noccEOS
 
          ! Only update for "interaction orbital pairs" - see which_pairs_occ
          if(dopair_occ(i,j)) then !YDoPair
 
             do b=1,nvirtAOS
                do c=1,nvirtAOS
                   do a=1,nvirtAOS
                      Y_tmp(a,b) = Y_tmp(a,b) + t2occ%elm4(c,j,a,i)*ThetaOCC%elm4(c,j,b,i)
                   end do
                end do
             end do

          end if

       end do
    end do
!OMP END DO NOWAIT
       print *, "I AM HERE 3 "

! Total Y matrix is found by summing all thread contributions
!OMP CRITICAL
dens%Y = dens%Y + Y_tmp
!OMP END CRITICAL

call mem_dealloc(Y_tmp)
call collect_thread_memory()
!OMP END PARALLEL
call mem_TurnOffThread_Memory()
       print *, "I AM HERE 4 "

! Multiply by 1/2 due to convention for Theta array [see above]
!dens%Y = 0.5E0_realk*dens%Y
    call LSTIMER('Y MATRIX',tcpu,twall,DECinfo%output)

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Pairs: MP2 one-el density MO: "
       print *, "";       call print_matrix_real(dens%Y,((noccAOS+nvirtAOS)/2),((noccAOS+nvirtAOS))/2)
       print *, "";    endif



    ! ****************************************************************************************
    !                Calculate contribution to occ-occ block of density matrix (X)           *
    ! ****************************************************************************************

    ! Use virtual partitioning scheme
    ! X_{ij} = sum_{abk} t_{ki}^{ba} * mult_{kj}^{ba}          [see the mp2dens structure]
    !        = 1/2* sum_{abk} t_{ki}^{ba} * Theta_{kj}^{ba}    [see construct_theta_array]
    ! -- where we only include the contributions if A and B belong to different fragments!

!Call mem_TurnONThread_Memory()
!!OMP PARALLEL DEFAULT(shared) PRIVATE(X_tmp,i,j,k,a,b)
!call init_threadmemvar()
!    call mem_alloc(X_tmp,noccAOS,noccAOS)
!    X_tmp = 0E0_realk
!
!!OMP DO SCHEDULE(dynamic,1)
!
!
!    do a=1,nvirtEOS
!       do b=1,nvirtEOS
!
!          ! Only update for "interaction orbital pairs" - see which_pairs_virt
!          if(dopair_virt(a,b)) then !XDoPair
!
!             do j=1,noccAOS
!                do k=1,noccAOS
!                   do i=1,noccAOS
!                      X_tmp(i,j) = X_tmp(i,j) + t2virt%elm4(b,k,a,i)*ThetaVIRT%elm4(b,k,a,j)
!                   end do
!                end do
!             end do
!
!          end if
!
!       end do
!    end do
!!OMP END DO NOWAIT
!
!! Total X matrix is found by summing all thread contributions
!!OMP CRITICAL
!dens%X = dens%X + X_tmp
!!OMP END CRITICAL
!
!call mem_dealloc(X_tmp)
!call collect_thread_memory()
!!OMP END PARALLEL
!call mem_TurnOffThread_Memory()
!
!! Multiply by 1/2 due to convention for Theta array [see above]
!    dens%X = 0.5E0_realk*dens%X
!    call LSTIMER('X MATRIX',tcpu,twall,DECinfo%output)
!
!
!
!
!    ! ****************************************************************************************
!    !                  Calculate contribution to virt-occ block of Phi matrix                *
!    ! ****************************************************************************************
!
!
!    ! Phivo(d,l) = sum_{cij} Theta(c,i,d,j) g(c,i,j,l)        [see the mp2dens structure]
!    ! -- where we only include the contributions if I and J belong to different fragments!
!
!
!call mem_TurnONThread_Memory()
!!OMP PARALLEL DEFAULT(shared) PRIVATE(Phivo_tmp,i,j,c,l,d)
!call init_threadmemvar()
 !   call mem_alloc(Phivo_tmp,nvirtAOS,nocctot)
 !   Phivo_tmp = 0E0_realk
 !
!!OMP DO SCHEDULE(dynamic,1)
 !
 !
 !   do i=1,noccEOS
 !      do j=1,noccEOS
 !
 !         if(dopair_occ(i,j)) then !PhivoDoPair
 !
 !            do d=1,nvirtAOS
 !               do l=1,nocctot
 !                  do c=1,nvirtAOS
 !                     Phivo_tmp(d,l) = Phivo_tmp(d,l) + ThetaOCC%elm4(c,i,d,j)*VOOO%elm4(c,i,j,l)
!                   end do
!                end do
!             end do
!
!          end if
!
!       end do
!    end do
!!OMP END DO NOWAIT
!
!! Total Phivo matrix is found by summing all thread contributions
!!OMP CRITICAL
!dens%Phivo = dens%Phivo + Phivo_tmp
!!OMP END CRITICAL
!
!    call mem_dealloc(Phivo_tmp)
!call collect_thread_memory()
!!OMP END PARALLEL
!call mem_TurnOffThread_Memory()
!
!    call LSTIMER('PHIVO MATRIX',tcpu,twall,DECinfo%output)
!
!
!
!
!
!    ! ****************************************************************************************
!    !                  Calculate contribution to occ-virt block of Phi matrix                *
!    ! ****************************************************************************************
!
!    ! Phiov(l,c) = sum_{abk} Theta(b,k,a,l) g(b,k,a,c)      [see the mp2dens structure])
!    ! -- where we only include the contributions if A and B belong to different fragments!
!
!
!call mem_TurnONThread_Memory()
!!OMP PARALLEL DEFAULT(shared) PRIVATE(Phiov_tmp,a,b,c,k,l)
!call init_threadmemvar()
!    call mem_alloc(Phiov_tmp,noccAOS,nvirtAOS)
!    Phiov_tmp = 0E0_realk
!
!!OMP DO SCHEDULE(dynamic,1)
!
!    do a=1,nvirtEOS
!       do b=1,nvirtEOS
!
!          if(dopair_virt(a,b)) then ! PhiovDoPair
!
!             do l=1,noccAOS
!                do c=1,nvirtAOS
!                   do k=1,noccAOS
!                      Phiov_tmp(l,c) = Phiov_tmp(l,c) + ThetaVIRT%elm4(b,k,a,l)*VOVV%elm4(b,k,a,c)
!                   end do
!                end do
!             end do
!
!          end if
!
!       end do
!    end do
!OMP END DO NOWAIT
!
! Total Phiov matrix is found by summing all thread contributions
!OMP CRITICAL
!dens%Phiov = dens%Phiov + Phiov_tmp
!!OMP END CRITICAL
!
!    call mem_dealloc(Phiov_tmp)
!call collect_thread_memory()
!!OMP END PARALLEL
!call mem_TurnOffThread_Memory()
!
!    call mem_dealloc(dopair_occ)
!    call mem_dealloc(dopair_virt)
!
!    call LSTIMER('PHIOV MATRIX',tcpu,twall,DECinfo%output)
!
  end subroutine pair_calculate_DEBUG



  !> \brief Calculate contributions to CCSD density from t2 amplitudes for single fragment.
  !> See mp2dens structure for equations.
  !> \author Dmytro Bykov
  !> \date November 2014
  subroutine single_calculate_DEBUG(MyFragment,t2occ,t2virt,ThetaOCC,ThetaVIRT,VOOO,VOVV,dens)


    implicit none

    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> t2 amplitudes t_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(tensor),intent(in) :: t2occ
    !> t2 amplitudes t_{KL}^{AB}, only for EOS orbitals using virtual partitioning, order: (A,K,B,L)
    type(tensor),intent(in) :: t2virt
    !> Theta array, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(tensor),intent(in) :: ThetaOCC
    !> Theta array, only for EOS orbitals using virtual partitioning, order:  (A,K,B,L)
    type(tensor),intent(in) :: ThetaVIRT
    !> (C I | J L) integrals stored as (C,I,J,L)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOOO
    !> (B K | A C) integrals stored as (B,K,A,C)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOVV
    !> MP2 density matrix contribution from given single fragment
    type(mp2dens),intent(inout) :: dens
    type(array2) :: Phivo, Phiov, X, Y
    integer :: i,j,noccEOS,nvirtEOS,noccAOS,nvirtAOS,nocctot
    logical :: something_wrong
    real(realk) :: tcpu,twall
    type(array4)              :: tbar, temp1, temp2, temp3

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    noccEOS = t2occ%dims(2) ! number of occupied EOS orbitals
    noccAOS = t2virt%dims(2) ! number of occupied AOS orbitals
    nvirtEOS = t2virt%dims(1) ! number of virtual EOS orbitals
    nvirtAOS = t2occ%dims(1) ! number of virtual AOS orbitals
    ! number of occupied core+valence orbitals (only different from noccAOS for frozen approx)
!    nocctot = VOOO%dims(4)   
    ! Just in case, zero matrices in dens structure
    dens%X = 0e0_realk
    dens%Y = 0e0_realk
    dens%Phivo = 0e0_realk
    dens%Phiov = 0e0_realk



    ! Sanity check 1: Density input structure is correct
    something_wrong=.false.
    if(noccAOS/=dens%nocc) something_wrong=.true.
    if(nvirtAOS/=dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'dens%nocc', dens%nocc
       write(DECinfo%output,*) 'dens%nvirt', dens%nvirt
       write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noccAOS
       write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvirtAOS
       call lsquit('single_calculate_mp2density: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

!    ! Sanity check 2: Array dimensions match for last two indices of t2 and Theta
!    if(noccEOS /= t2occ%dims(3)) something_wrong=.true.
!    if(nvirtAOS /= t2occ%dims(4)) something_wrong=.true.
!    if(nvirtEOS /= t2virt%dims(3)) something_wrong=.true.
!    if(noccAOS /= t2virt%dims(4)) something_wrong=.true.
!    if(noccEOS /= Thetaocc%dims(3)) something_wrong=.true.
!    if(nvirtAOS /= Thetaocc%dims(4)) something_wrong=.true.
!    if(nvirtEOS /= Thetavirt%dims(3)) something_wrong=.true.
!    if(noccAOS /= Thetavirt%dims(4)) something_wrong=.true.
!    if(something_wrong) then
!       write(DECinfo%output,*) 'EOS: nocc, nvirt', noccEOS, nvirtEOS
!       write(DECinfo%output,*) 'AOS: nocc, nvirt', noccAOS, nvirtAOS
!       write(DECinfo%output,*) 't2occ dims     ', t2occ%dims
!       write(DECinfo%output,*) 'Thetaocc dims  ', Thetaocc%dims
!       write(DECinfo%output,*) 't2virt dims    ', t2virt%dims
!       write(DECinfo%output,*) 'Thetavirt dims ', Thetavirt%dims
!       call lsquit('single_calculate_mp2density: &
!            & Something wrong with amplitude and Theta dimensions!',DECinfo%output)
!    end if


    ! ******************************************************************************************
    !                Calculate contribution to virt-virt block of density matrix (Y)           *
    ! ******************************************************************************************

    ! Use occupied partitioning scheme
    ! Y_{ab} = sum_{cij} t_{ji}^{ca} * mult_{ji}^{cb}           [see the mp2dens structure]
    !        = 1/2* sum_{cij} t_{ji}^{ca} * Theta_{ji}^{cb}     [see construct_theta_array]
    Y = array2_init([nvirtAOS,nvirtAOS])

    !container for amps voov shape
    temp1 = array4_init([nvirtAOS,noccEOS,noccEOS,nvirtAOS])
    temp2 = array4_init([nvirtAOS,noccEOS,noccEOS,nvirtAOS])

    !sort amps(ckbl) -> bckl
    call array_reorder_4d(1.0E0_realk,t2occ%elm4,nvirtAOS,noccEOS,nvirtAOS,noccEOS,[1,2,4,3],0.0E0_realk,temp1%val)
    call array_reorder_4d(1.0E0_realk,ThetaOCC%elm4,nvirtAOS,noccEOS,nvirtAOS,noccEOS,[1,2,4,3],0.0E0_realk,temp2%val)

    call array4_contract3(temp1,temp2,Y)
!    call array4_contract3(t2occ,ThetaOCC,Y)
!    dens%Y(1:nvirtAOS,1:nvirtAOS) = 0.5e0_realk*Y%val(1:nvirtAOS,1:nvirtAOS)

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "MP2 Y matrix SINGLE: "
       print *, "";       write(*, *) ''
        do i=1,nvirtAOS
          do j=1,nvirtAOS
            write(*,'(f16.10)',advance='no') Y%val(i,j)
          end do
         write(*, *) ''
       end do
       print *, "";    endif

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Pairs: MP2 one-el density MO: "
       print *, "";       call print_matrix_real(dens%Y,((noccAOS+nvirtAOS))/2,((noccAOS+nvirtAOS))/2)
       print *, "";    endif

    call array2_free(Y)
    call LSTIMER('Y MATRIX',tcpu,twall,DECinfo%output)


!    ! ****************************************************************************************
!    !                Calculate contribution to occ-occ block of density matrix (X)           *
!    ! ****************************************************************************************
!
!    ! Use virtual partitioning scheme
!    ! X_{ij} = sum_{abk} t_{ki}^{ba} * mult_{kj}^{ba}          [see the mp2dens structure]
!    !        = 1/2* sum_{abk} t_{ki}^{ba} * Theta_{kj}^{ba}    [see construct_theta_array]
!    X = array2_init([noccAOS,noccAOS])
!
!    call array4_contract3(t2virt,ThetaVIRT,X)
!    dens%X(1:noccAOS,1:noccAOS) = 0.5e0_realk*X%val(1:noccAOS,1:noccAOS)
!
!    call array2_free(X)
!    call LSTIMER('X MATRIX',tcpu,twall,DECinfo%output)
!
!    ! X and Y matrices collected and transformed to AO basis (unrelaxed corr. density)
!    call get_unrelaxed_corrdens_in_AO_basis(dens%nocc,dens%nvirt,dens%nbasis,MyFragment%Co,&
!         & MyFragment%Cv,dens%X,dens%Y,dens%rho)
!
!    !print debug info if needed
!    if( DECinfo%PL>2 )then
!       print *, ""; !      print *, "MP2 Y+X  SINGLE: "
 !      print *, ""; !      write(*, *) ''
 !       do i=1,dens%nbasis
 !         do j=1,dens%nbasis
 !           write(*,'(f16.10)',advance='no') dens%rho(i,j)
 !         end do
 !        write(*, *) ''
 !      end do
 !      print *, ""; !   endif
!
!
!    ! ****************************************************************************************
!    !                  Calculate contribution to virt-occ block of Phi matrix                *
!    ! ****************************************************************************************
!
!
!    ! Phivo(D,L) = sum_{CIJ} Theta(C,I,J,D) g(C,I,J,L)     [see the mp2dens structure]
!    ! Note: L is both core+valence, also for frozen core approximation.
!    Phivo = array2_init([nvirtAOS,nocctot])
!    call array4_contract3(ThetaOCC,VOOO,Phivo)
 !   dens%Phivo(1:nvirtAOS,1:nocctot) = Phivo%val(1:nvirtAOS,1:nocctot)
 !
 !   call array2_free(Phivo)
 !   call LSTIMER('PHIVO MATRIX',tcpu,twall,DECinfo%output)
 !
 !
 !
 !
 !   ! ****************************************************************************************
 !   !                  Calculate contribution to occ-virt block of Phi matrix                *
 !   ! ****************************************************************************************
!
!
!    ! Phiov(L,C) = sum_{ABK} Theta(B,K,A,L) g(B,K,A,C)      [see the mp2dens structure]
!    Phiov = array2_init([noccAOS,nvirtAOS])
!    call array4_contract3(ThetaVIRT,VOVV,Phiov)
!    dens%Phiov(1:noccAOS,1:nvirtAOS) = Phiov%val(1:noccAOS,1:nvirtAOS)
!
!    call array2_free(Phiov)
!
!    call LSTIMER('PHIOV MATRIX',tcpu,twall,DECinfo%output)
!


  end subroutine single_calculate_DEBUG


!> Calculate CCSD gradient from fragment contributions.
!> \author Dmytro Bykov
!> \date November 2014
subroutine get_CCSDgradient_main(MyMolecule,mylsitem,DHF,mp2gradient,fullgrad)


  implicit none

  !> Full molecule info
  type(fullmolecule), intent(in) :: MyMolecule
  !> LS setting info
  type(lsitem), intent(inout) :: mylsitem
  !> HF density matrix
  type(matrix),intent(in) :: DHF
  !> MP2 gradient
  real(realk), intent(inout) :: mp2gradient(3,MyMolecule%natoms)
  !> Full MP2 gradient structure
  type(fullmp2grad),intent(inout) :: fullgrad
  integer :: nvirt,nocc,nbasis,ncore,nval
  type(matrix) :: rho,DMP2,Phivo,Phiov,Phioo
  type(matrix) :: C,F,W
  type(array2) :: DMP2AO,dia,dai,dij,dab 
  type(array2) :: diaao,daiao,dijao,dabao 
  logical :: full_equation
  real(realk), dimension(3) :: HFdipole, MP2dipole
  real(realk) :: tcpu,twall, tcpu1,tcpu2, twall1,twall2, TrSrho
  real(realk),pointer :: basis(:,:)
  type(matrix) :: Phi
  type(matrix) :: Cocc,Cvirt
  integer i,j,k

  call LSTIMER('START',tcpu,twall,DECinfo%output)
  call LSTIMER('START',tcpu1,twall1,DECinfo%output)
  if(MyMolecule%mem_distributed)then
     call lsquit("ERROR(get_CCSDgradient_main) not implemented for distributed matrices in molecule",-1)
  endif

print*, "ready to calculate the dipole moment"

  ! Easy reference to molecule info
  ! *******************************
  nbasis = MyMolecule%nbasis
  nocc = MyMolecule%nocc
  nvirt = MyMolecule%nvirt
  ncore = MyMolecule%ncore
  nval = MyMolecule%nval

  ! Get MO coefficients in matrix form
  ! **********************************
  call mat_init(Cocc,nbasis,nocc)
  call mat_init(Cvirt,nbasis,nvirt)
  call mat_set_from_full(MyMolecule%Co%elm2(1:nbasis,1:nocc), 1E0_realk, Cocc)
  call mat_set_from_full(MyMolecule%Cv%elm2(1:nbasis,1:nvirt), 1E0_realk, Cvirt)

  DMP2AO = array2_init([nbasis,nbasis])

  dij    = array2_init([nocc,nocc])
  dia    = array2_init([nocc,nvirt])
  dai    = array2_init([nvirt,nocc])
  dab    = array2_init([nvirt,nvirt])

  dijao    = array2_init([nbasis,nbasis])
  diaao    = array2_init([nbasis,nbasis])
  daiao    = array2_init([nbasis,nbasis])
  dabao    = array2_init([nbasis,nbasis])

  dij%val(1:nocc,1:nocc)                 = fullgrad%rho(1:nocc,1:nocc)
  dia%val(1:nocc,(nocc+1):nvirt)         = fullgrad%rho(1:nocc,(nocc+1):nvirt)
  dai%val((nocc+1):nvirt,1:nocc)         = fullgrad%rho((nocc+1):nvirt,1:nocc) 
  dab%val((nocc+1):nvirt,(nocc+1):nvirt) = fullgrad%rho((nocc+1):nvirt,(nocc+1):nvirt) 

  print *,'dij befor diagonal' 
  call print_matrix_real(dij%val,nocc,nocc)

  do i=1,nocc
   dij%val(i,i)= 2E0_realk + dij%val(i,i)
  end do

  print *,'dij after diagonal' 
  call print_matrix_real(dij%val,nocc,nocc)

  call dec_simple_basis_transform2(nbasis,nocc,&
       & MyMolecule%Co%elm2,dij%val,dijao%val) 

  call dec_simple_basis_transform2(nbasis,nvirt,&
       & MyMolecule%Cv%elm2,dab%val,dabao%val) 
  
  write(*, *) ''
  write(*, *) 'Co basis'
  k=1
  do i=1,nocc
    do j=1,nbasis
      write(*,'(f16.10)',advance='no') MyMolecule%Co%elm2(j,i)
      k=k+1
    end do
    write(*, *) ''
  end do

  call dec_diff_basis_transform2(nbasis,nocc,nvirt,MyMolecule%Co%elm2,&
       & MyMolecule%Cv%elm2,dia%val,diaao%val) 

  call dec_diff_basis_transform2(nbasis,nvirt,nocc,MyMolecule%Cv%elm2,&
       & MyMolecule%Co%elm2,dai%val,daiao%val) 

  DMP2AO = dijao+diaao+daiao+dabao

  call mat_init(rho,nbasis,nbasis)
  k=1
  do i=1,rho%ncol
    do j=1,rho%nrow
      rho%elms(k) = DMP2AO%val(j,i)
      k=k+1
    end do
  end do
 
  write(*, *) ''
  write(*, *) 'rho in AO basis'
  k=1
  do i=1,rho%ncol
    do j=1,rho%nrow
      write(*,'(f16.10)',advance='no') rho%elms(k)
      k=k+1
    end do
    write(*, *) ''
  end do


  ! Full density matrix = Hartree-Fock + correlation contribution
  ! *****************************************************************
  call mat_init(DMP2,nbasis,nbasis)
  DMP2 = rho

  ! Also calculate and print HF and MP2 electric dipole moments
  ! ***********************************************************
  call get_HF_and_CCSD_dipole_moments(mylsitem,DHF,DMP2,HFdipole,MP2dipole)
  call print_HF_and_CCSD_dipoles(HFdipole, MP2dipole)

  ! Save MP2 density matrices to file
  ! *********************************
  !call save_mp2density_matrices_to_file(DMP2,rho)


  GradientCalc: if(DECinfo%gradient) then ! only for gradient, simple MP2 density

!     ! Get MP2 reorthonormalization matrix
!     ! ***********************************
!     ! MO coefficient matrix and Fock matrix in type matrix form
!     call mat_init(C,nbasis,nbasis)
!     call mat_init(F,nbasis,nbasis)
!     call mem_alloc(basis,nbasis,nbasis)
!     basis(1:nbasis,1:nocc) = MyMolecule%Co(1:nbasis,1:nocc)
!     basis(1:nbasis,nocc+1:nbasis) = MyMolecule%Cv(1:nbasis,1:nvirt)
!     call mat_set_from_full(basis(1:nbasis,1:nbasis), 1E0_realk, C)
!     call mem_dealloc(basis)
!     call mat_set_from_full(MyMolecule%fock(1:nbasis,1:nbasis), 1E0_realk, F)
!
!     ! Reorthonormalization matrix W
!     call get_mp2_reorthonormalization_matrix(F,D,Phi,rho,C,MyLsitem,W)
!     call mat_free(C)
!
!
!     ! Finally, we can determine the MP2 gradient using the calculated matrices
!     ! ************************************************************************
!     call calculate_MP2_gradient(D,F,rho,W,DMP2,MyLsitem,FullGrad)
!     mp2gradient = fullgrad%mp2gradient
!
!
!     ! Free stuff
!     ! **********
!     call mat_free(F)
!     call mat_free(W)

  else
     write(DECinfo%output,*) 'Only CCSD density requested, skipping CCSD gradient part'

  end if GradientCalc


!  call mat_free(DMP2)
!  call mat_free(rho)
!  call mat_free(Phi)

  call LSTIMER('GET MP2GRAD',tcpu,twall,DECinfo%output)
  call LSTIMER('START',tcpu2,twall2,DECinfo%output)

end subroutine get_CCSDgradient_main

  !> \brief Calculate electric dipole moment at the HF and CCSD levels of theory.
  !> \author Dmytro Bykov
  !> \date December 2014
  subroutine get_HF_and_CCSD_dipole_moments(mylsitem,DHF,DMP2,HFdipole,MP2dipole)


    implicit none

    !> LS setting info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: DHF
    !> CCSD density matrix 
    type(matrix),intent(in) :: DMP2
    type(matrix),target :: DipoleIntegral(3)
    real(realk), intent(inout) :: HFdipole(3)
    real(realk), intent(inout) :: MP2dipole(3)
    real(realk),dimension(3) :: NucDipole, HF_el_dipole, MP2_el_dipole
    integer :: nbasis,i
    character(len=7) :: string

    ! Calculate pure nuclear contribution to dipole (of course the same for HF and MP2)
    ! *********************************************************************************
    call II_get_nucdip(mylsitem%setting,NucDipole)

    ! Calculate electronic dipole matrices
    ! ************************************
    string(1:7) = 'DIPLEN '
    nbasis = DHF%nrow
    do i=1,3
       call mat_init(DipoleIntegral(i),nbasis,nbasis)
    end do
    call II_get_integral(DECinfo%output, DECinfo%output ,mylsitem%setting, DipoleIntegral,3,string)


    ! Electronic contribution to HF dipole
    ! ************************************
    ! Dot product of dipole integrals and HF density (factor two, due to double occupation for closed-shell)
    do i=1,3
       HF_el_dipole(i) = 2E0_realk*mat_dotproduct(DHF,DipoleIntegral(i))
    end do


    ! Electronic contribution to CCSD dipole
    ! ************************************
    ! Dot product of dipole integrals and MP2 density
    do i=1,3
       MP2_el_dipole(i) = mat_dotproduct(DMP2,DipoleIntegral(i))
    end do


    ! Total dipoles
    ! *************
    ! "Sum" of nuclear and electronic contributions. However, the
    ! dipole by definition is MINUS the derivative of the energy/Lagrangian with respect
    ! to an electric field component. Our density matrices are associated with an HF or an MP2 energy.
    ! Therefore, the electronic contribution should be multiplied by (-1),
    ! while the nuclear contribution has the correct sign when calculated using II_get_nucdip.

    ! HF total dipole
    do i=1,3
       HFdipole(i) = -HF_el_dipole(i) + NucDipole(i)
    end do

    ! MP2 total dipole
    do i=1,3
       MP2dipole(i) = -MP2_el_dipole(i) + NucDipole(i)
    end do
   
    ! Free stuff
    ! **********
    do i=1,3
       call mat_free(DipoleIntegral(i))
    end do

  end subroutine get_HF_and_CCSD_dipole_moments

  !> \brief Print HF and CCSD dipole info to LSDALTON.OUT
  !> \author Dmytro Bykov
  !> \date December 2014
  subroutine print_HF_and_CCSD_dipoles(HFdipole, MP2dipole)


    implicit none

    !> Hartree-Fock dipole moment
    real(realk), dimension(3),intent(in) :: HFdipole
    !> MP2 dipole moment
    real(realk), dimension(3),intent(in) :: MP2dipole
    real(realk) :: au_to_debye, au_to_SI, normHFdipole,normMP2dipole
    integer :: i

    normHFdipole=0E0_realk
    normMP2dipole=0E0_realk
    do i=1,3
       normHFdipole = normHFdipole + HFdipole(i)**2
       normMP2dipole = normMP2dipole + MP2dipole(i)**2
    end do
    normHFdipole = sqrt(normHFdipole)
    normMP2dipole = sqrt(normMP2dipole)

    au_to_debye=2.54175
    au_to_SI=8.47835
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '    DIPOLE MOMENTS AT THE HARTREE-FOCK AND CCSD LEVELS OF THEORY '
    write(DECinfo%output,*) '    =========================================================== '
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(6X,A)') '                 HF: Permanent dipole moment'
    write(DECinfo%output,'(6X,A)') '                 ---------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(5X,3g18.6)') normHFDipole, au_to_debye*normHFDipole, au_to_SI*normHFDipole
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(6X,A)') '                 HF: Dipole moment components'
    write(DECinfo%output,'(6X,A)') '                 ----------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'x', HFdipole(1), &
         & au_to_debye*HFdipole(1), au_to_SI*HFdipole(1)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'y', HFdipole(2), &
         & au_to_debye*HFdipole(2), au_to_SI*HFdipole(2)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'z', HFdipole(3), &
         & au_to_debye*HFdipole(3), au_to_SI*HFdipole(3)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    write(DECinfo%output,'(6X,A)') '                 CCSD: Permanent dipole moment'
    write(DECinfo%output,'(6X,A)') '                 -----------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(5X,3g18.6)') normMP2Dipole, au_to_debye*normMP2Dipole, au_to_SI*normMP2Dipole
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(6X,A)') '                 CCSD: Dipole moment components'
    write(DECinfo%output,'(6X,A)') '                 ------------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'x', MP2dipole(1), &
         & au_to_debye*MP2dipole(1), au_to_SI*MP2dipole(1)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'y', MP2dipole(2), &
         & au_to_debye*MP2dipole(2), au_to_SI*MP2dipole(2)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'z', MP2dipole(3), &
         & au_to_debye*MP2dipole(3), au_to_SI*MP2dipole(3)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine print_HF_and_CCSD_dipoles


  !> \brief Update full CCSD gradient structure with fragment contribution.
  !> \author Dmytro Bykov
  !> \date November 2014
  subroutine update_full_CCSDgradient(fraggrad,fullgrad)

    implicit none
    !> Full MP2 gradient structure
    type(FullMP2grad),intent(inout) :: fullgrad
    !> Fragment MP2 gradient contribution
    type(mp2grad),intent(inout) :: fraggrad
    integer :: i,j,ix,jx,a,b,ax,bx,iFull,jFull
    
    !DEBUG - everything in MO
    !Update the Dab matrix
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Dab before updating rho"
       call print_matrix_real(fraggrad%dens%rho,fraggrad%dens%nvirt,fraggrad%dens%nvirt)
       print *, "";    endif

    do i=1,fraggrad%dens%nvirt 
      iFull = fraggrad%dens%basis_idx(i)
      do j=1,fraggrad%dens%nvirt 
        jFull = fraggrad%dens%basis_idx(j)
        fullgrad%rho(iFull,jFull) = fullgrad%rho(iFull,jFull)+fraggrad%dens%rho(i,j)
      end do
    end do

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Dab after updating rho"
       call print_matrix_real(fullgrad%rho,fullgrad%nbasis,fullgrad%nbasis)
       print *, "";    endif

!    !Update the Dij matrix
!    if( DECinfo%PL>2 )then
!       print *, "";!       print *, "Dij before updating rho"
!       call print_matrix_real(fraggrad%dens%Y,fraggrad%dens%nvirt,fraggrad%dens%nvirt)
!       print *, "";!    endif
!    
!    do i=1,fraggrad%dens%nvirt
!      iFull = fraggrad%dens%basis_idx(i)
!      do j=1,fraggrad%dens%nvirt 
!        jFull = fraggrad%dens%basis_idx(j)
!        fullgrad%rho(iFull,jFull) = fullgrad%rho(iFull,jFull)+fraggrad%dens%Y(i,j)
!      end do
!    end do
!    
!    !print debug info if needed
!    if( DECinfo%PL>2 )then
!       print *, "";!       print *, "Dij after updating rho"
!       call print_matrix_real(fullgrad%rho,fullgrad%nbasis,fullgrad%nbasis)
!       print *, "";!    endif

    ! Update full unrelaxed corr. density matrix rho in AO basis
    ! ----------------------------------------------------------
    do j=1,fraggrad%dens%nbasis
       jx=fraggrad%dens%basis_idx(j)
       do i=1,fraggrad%dens%nbasis
          ix=fraggrad%dens%basis_idx(i)
!          fullgrad%rho(ix,jx) = fullgrad%rho(ix,jx) + fraggrad%dens%rho(i,j)
       end do
    end do


    ! Update full Phi matrix in AO basis
    ! ----------------------------------
    do j=1,fraggrad%dens%nbasis
       jx=fraggrad%dens%basis_idx(j)
       do i=1,fraggrad%dens%nbasis
          ix=fraggrad%dens%basis_idx(i)
!          fullgrad%Phi(ix,jx) = fullgrad%Phi(ix,jx) + fraggrad%PhiAO(i,j)
       end do
    end do


    ! Update correlation energy
    ! -------------------------
!    fullgrad%Ecorr = fullgrad%Ecorr + fraggrad%dens%energy


    ! Loop over elements in fragment Ltheta vector and update full Ltheta vector
    ! --------------------------------------------------------------------------
    do j=1,fraggrad%natoms
       jx=fraggrad%atoms_idx(j)
       do i=1,3  ! Cartesion components of Ltheta
          ix=i

          ! Update full matrix Ltheta element by fragment (i,j) contribution
!          fullgrad%Ltheta(ix,jx) = fullgrad%Ltheta(ix,jx) + fraggrad%Ltheta(i,j)

       end do
    end do


  end subroutine update_full_CCSDgradient

  !> \brief Update CCSD molecular density with fragment contribution.
  !> \author Dmytro Bykov
  !> \date November 2014
  subroutine update_CCSDdensity(fragment,fraggrad,fullgrad)

    implicit none
    !> Full MP2 gradient structure
    type(FullMP2grad),intent(inout) :: fullgrad
    !> Fragment MP2 gradient contribution
    type(mp2grad),intent(inout) :: fraggrad
    type(decfrag) :: fragment
    integer :: i,j,ix,jx,a,b,ax,bx,iFull,jFull
    
    !Update the Dab matrix
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Ready to update Y:"
       call print_matrix_real(fraggrad%dens%Y,fraggrad%dens%nvirt,fraggrad%dens%nvirt)
       print *, "";    endif

    do i=1,fragment%nvirtAOS !fraggrad%dens%nvirt 
      iFull = fragment%virtAOSidx(i)
      do j=1,fragment%nvirtAOS !fraggrad%dens%nvirt 
        jFull = fragment%virtAOSidx(j)
        fullgrad%rho(fullgrad%nocc+iFull,fullgrad%nocc+jFull) = &
        fullgrad%rho(fullgrad%nocc+iFull,fullgrad%nocc+jFull) + fraggrad%dens%Y(i,j)
      end do
    end do

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Here is what I've got for Y:"
       call print_matrix_real(fullgrad%rho,fullgrad%nbasis,fullgrad%nbasis)
       print *, "";    endif

    !Update the Dij matrix
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Ready to update X:"
       call print_matrix_real(fraggrad%dens%X,fraggrad%dens%nocc,fraggrad%dens%nocc)
       print *, "";    endif

    do i=1,fraggrad%dens%nocc
      iFull = fragment%OccAOSidx(i)
      do j=1,fraggrad%dens%nocc
        jFull = fragment%OccAOSidx(j)
        fullgrad%rho(iFull,jFull) = fullgrad%rho(iFull,jFull)+fraggrad%dens%X(i,j)
      end do
    end do

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Here is what I've got for X:"
       call print_matrix_real(fullgrad%rho,fullgrad%nbasis,fullgrad%nbasis)
       print *, "";    endif

    !Update the Dai and Dia matricies
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Ready to update off-diagonal terms:"
       call print_matrix_real(fraggrad%dens%Phivo,fraggrad%dens%nvirt,fraggrad%dens%nocc)
       print *, "";       call print_matrix_real(fraggrad%dens%Phiov,fraggrad%dens%nocc,fraggrad%dens%nvirt)
       print *, "";    endif

    do i=1,fraggrad%dens%nocc
      iFull = fragment%OccAOSidx(i)
      do j=1,fraggrad%dens%nvirt
        jFull = fragment%virtAOSidx(j)
        fullgrad%rho(iFull,fullgrad%nocc+jFull) = &
        fullgrad%rho(iFull,fullgrad%nocc+jFull) + &
                                    &fraggrad%dens%Phiov(i,j)
      end do
    end do
    do i=1,fraggrad%dens%nocc
      iFull = fragment%OccAOSidx(i)
      do j=1,fraggrad%dens%nvirt
        jFull = fragment%virtAOSidx(j)
        fullgrad%rho(fullgrad%nocc+jFull,iFull) = &
        fullgrad%rho(fullgrad%nocc+jFull,iFull) + &
                                    &fraggrad%dens%Phivo(j,i)
      end do
    end do

    !print debug info if needed
    if( DECinfo%PL>2 )then
       print *, "";       print *, "Here is what I've got for the full density:"
       call print_matrix_real(fullgrad%rho,fullgrad%nbasis,fullgrad%nbasis)
       print *, "";    endif

  end subroutine update_CCSDdensity



end module ccsd_gradient_module

