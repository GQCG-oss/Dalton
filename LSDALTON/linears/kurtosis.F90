! *******************************************************
! *                                                     *
! *  File which contains necessary routines for         *
! *  fourth central moment localization (JCP, 137, 2012)* 
! *  Author: Ida-Marie Hoeyvik                          *
! *                                                     *
! *******************************************************

module kurtosis
use precision
use matrix_module
use typedef
use matrix_util
use matrix_operations
use matrix_operations_aux
use IntegralInterfaceMOD

TYPE PFMitem 

! General settings
!**********************************
! power of localization
integer :: m
! Number of orbitals to localize
integer :: norb
! Number of basis functiona
integer :: nbas
! include crossterms
logical :: crossterms
! Function value
real(realk) :: kurt_val
! used to run test case
logical :: TESTCASE

! Needed for localization
!***********************************
! Matrix of diagonal Hessian elements
type(matrix),pointer :: P

! Matrices of cartesian moments
!***********************************
type(matrix) :: carmom(35) 
! In MO basis
type(matrix) :: carmomMO(24) 
! In AO basis
type(matrix) :: carmomAO(24) 

! Prefactor vectors
!***********************************
real(realk),pointer :: omega(:)

ENDTYPE PFMitem

contains

!> \brief  runs test of gradient and linear transformation elements
!> \param cmo input matrix of coefficients
  subroutine kurtosis_test(ls,cmo,nbas,norb)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: HK,kappa
    type(lsitem) :: ls
    type(matrix) :: cmo,G,Prec
    integer :: norb,nbas
    integer :: i,p,q
    real(realk),pointer :: testmat(:,:)
    real(realk) :: Gpq,Hpqpq

    call mem_alloc(testmat,norb,norb)
    p=2
    q=4
    testmat=0.0_realk
    testmat(p,q)=1.0_realk
    testmat(q,p)=-1.0_realk
    call mat_init(G,cmo%ncol,cmo%ncol)
    call mat_init(Prec,cmo%ncol,cmo%ncol)
    call mat_init(kappa,cmo%ncol,cmo%ncol)
    call mat_set_from_full(testmat,1.0_realk,kappa)
    call mat_init(HK,cmo%ncol,cmo%ncol)
    call mem_alloc(KURT%omega,cmo%ncol)
    KURT%nbas=nbas
    KURT%norb=norb
    KURT%m = 2
    KURT%crossterms = .true.

    do i=1,24
       call mat_init(KURT%carmomMO(i),norb,norb)
       call mat_init(KURT%carmomAO(i),nbas,nbas)
    end do
    call compute_carmom(KURT,ls)
    do i=1,24
       call carmomAO_to_carmomMO(KURT%carmomAO(i),cmo,KURT%carmomMO(i))
    end do
    call compute_omega(KURT%carmomMO,KURT%omega,KURT%crossterms,cmo%ncol)
    call compute_gradient(KURT,G,G%ncol)
    call FiniteDiff_gradient(KURT,cmo,Gpq)
    write(ls%lupri,'(a,f13.5)') 'TEST: gradient                       ',G%elms(nbas*(q-1)+p)
    write(ls%lupri,'(a,f13.5)') 'TEST: finite difference gradient     ',Gpq
    call compute_LinTrans(KURT,kappa,G,HK) 
    call FiniteDiff_HK(KURT,cmo,Hpqpq)
    write(ls%lupri,'(a,f13.5)') 'TEST: linear transformations         ',HK%elms(nbas*(q-1)+p)
    write(ls%lupri,'(a,f13.5)') 'TEST: finite difference  LT/diagonal ', Hpqpq
    call kurtosis_precond_matrix(KURT,Prec)
    write(ls%lupri,'(a,f13.5)') 'TEST: diagonal elements              ', Prec%elms(nbas*(q-1)+p)

    do i=1,24
       call mat_free(KURT%carmomMO(i))
       call mat_free(KURT%carmomAO(i))
    end do
    call mem_dealloc(KURT%omega)
    call mat_free(G)
    call mat_free(HK)
    call mat_free(Prec)
    call mat_free(kappa)
    call mem_dealloc(testmat)
  end subroutine kurtosis_test

!> \brief Computes cartesian matrices in AO basis
!> \brief  Do only once per calculation (incl.core/valence/virt)
  subroutine kurt_initAO(KURT,ls,nbas)
    implicit none
    type(PFMitem) :: KURT
    type(lsitem)       :: ls
    integer :: i,nbas

    KURT%nbas = nbas

    do i=1,24
       call mat_init(KURT%carmomAO(i),KURT%nbas,KURT%nbas)
    end do

    call compute_carmom(KURT,ls) 

  end subroutine kurt_initAO

!> \brief initialize cartesian moment matrices in MO basis
!> \brief Do once for each core, valence and virtual
  subroutine kurt_initMO(KURT,cmo,Memreduced)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: cmo
    logical,optional :: Memreduced
    type(lsitem) :: ls
    integer :: i
    logical :: DoMemreduced

    DoMemreduced = .FALSE.
    IF(present(Memreduced))DoMemreduced = Memreduced

    KURT%norb = cmo%ncol
    call mem_alloc(KURT%omega,KURT%norb)
    IF(DoMemreduced)THEN
       call compute_memreduced_omega(KURT%carmomAO,KURT%omega,KURT%crossterms,KURT%norb,cmo)
    ELSE
       do i=1,24
          call mat_init(KURT%carmomMO(i),cmo%ncol,cmo%ncol)
          call carmomAO_to_carmomMO(KURT%carmomAO(i),cmo,KURT%carmomMO(i))
       end do
       call compute_omega(KURT%carmomMO,KURT%omega,KURT%crossterms,KURT%norb)
    ENDIF

  end subroutine kurt_initMO

!> \brief Free cartesian matrices in AO basis. Once per calculation
  subroutine kurt_freeAO(KURT)
    implicit none
    type(PFMitem) :: kurt
    integer :: i

    do i=1,24
       call mat_free(KURT%carmomAO(i))
    end do

  end subroutine kurt_freeAO

!> \brief Free cartesian matrices in MO basis. Done for both core/valence/virtual
  subroutine kurt_freeMO(KURT,Memreduced)
    implicit none
    type(PFMitem) :: kurt
    logical,optional :: Memreduced
    integer :: i
    logical :: DoMemreduced
    
    DoMemreduced = .FALSE.
    IF(present(Memreduced))DoMemreduced = Memreduced

    call mem_dealloc(KURT%omega)
    IF(.NOT.DoMemreduced)THEN
       do i=1,24
          call mat_free(KURT%carmomMO(i))
       end do
    ENDIF
  end subroutine kurt_freeMO

!> \brief Update CM matrices with current CMO: Mnew = CMO^T M_AO CMO
  subroutine kurt_updateAO(KURT,cmo)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: temp,cmo
    integer :: i

    call mat_init(temp,KURT%norb,KURT%nbas)
    do i=1,24
       call mat_mul(cmo,KURT%carmomAO(i),'T','n',1.0_realk,0.0_realk,temp)
       call mat_mul(temp,cmo,'n','n',1.0_realk,0.0_realk,KURT%carmomMO(i))
    end do
    call mat_free(temp)

    call compute_omega(KURT%carmomMO,KURT%omega,KURT%crossterms,KURT%norb)

  end subroutine kurt_updateAO

!> \brief Update CM matrices in MO basis with exp. M_MO = expX^T M_MP expX 
  subroutine kurt_update(KURT,expX)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: temp,expX
    integer :: i

    call mat_init(temp,KURT%norb,KURT%norb)
    do i=1,24
       call mat_mul(expX,KURT%carmomMO(i),'T','n',1.0_realk,0.0_realk,temp)
       call mat_mul(temp,expX,'n','n',1.0_realk,0.0_realk,KURT%carmomMO(i))
    end do
    call mat_free(temp)

    call compute_omega(KURT%carmomMO,KURT%omega,KURT%crossterms,KURT%norb)

  end subroutine kurt_update

!> \brief Compute function value sum_i = FM_i^m
  subroutine kurt_value(KURT)
    implicit none
    type(PFMitem) :: KURT
    integer :: i
    real(realk) ::temp(KURT%norb)

    do i=1,KURT%norb
       temp(i)=KURT%omega(i)**KURT%m
    end do

    KURT%kurt_val = sum(temp)

  end subroutine kurt_value



!> \brief Compute gradient for fourth central moment 
!> \param G gradient (output)
  subroutine compute_gradient(KURT,G,norb)
    implicit none
    type(PFMitem) :: KURT
    integer :: norb
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz
    real(realk) :: factor(norb)
    type(matrix) :: G,Gt
    integer :: i,j,m
    real(realk) ::diatemp(norb),diatemp2(norb),diatemp3(norb),temp(norb)


    call mat_zero(G)

    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24

    call compute_omega(KURT%carmomMO,KURT%omega,KURT%crossterms,KURT%norb)
    m=KURT%m
    factor=KURT%omega
    if (m>1) then
       do i=1,norb
          factor(i) = factor(i)**(m-1)
       end do
    else
       factor=1.0_realk
    end if

    call mat_zero(G)
    call mat_dmul(factor,KURT%carmomMO(xxxx),'T',-dble(2*m),1.0_realk,G)
    call mat_dmul(factor,KURT%carmomMO(yyyy),'T',-dble(2*m),1.0_realk,G)
    call mat_dmul(factor,KURT%carmomMO(zzzz),'T',-dble(2*m),1.0_realk,G)

    diatemp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
    diatemp2=diatemp*factor
    call mat_dmul(diatemp2,KURT%carmomMO(xxx),'T',dble(8*m),1.0_realk,G)


    diatemp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(y))
    diatemp2=diatemp*factor
    call mat_dmul(diatemp2,KURT%carmomMO(yyy),'T',dble(8*m),1.0_realk,G)


    diatemp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(z))
    diatemp2=diatemp*factor
    call mat_dmul(diatemp2,KURT%carmomMO(zzz),'T',dble(8*m),1.0_realk,G)


    diatemp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
    diatemp = diatemp*diatemp*factor
    call mat_dmul(diatemp,KURT%carmomMO(xx),'T',-dble(12*m),1.0_realk,G)


    diatemp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(y))
    diatemp = diatemp*diatemp*factor
    call mat_dmul(diatemp,KURT%carmomMO(yy),'T',-dble(12*m),1.0_realk,G)

    diatemp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(z))
    diatemp = diatemp*diatemp*factor
    call mat_dmul(diatemp,KURT%carmomMO(zz),'T',-dble(12*m),1.0_realk,G)


    diatemp=0.0_realk;diatemp2=0.0_realk;diatemp3=0.0_realk;temp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
    call mat_extract_diagonal(diatemp2,KURT%carmomMO(xx))
    call mat_extract_diagonal(diatemp3,KURT%carmomMO(xxx))
    temp=diatemp*diatemp*diatemp
    diatemp= (8.0_realk*diatemp3-24.0_realk*diatemp2*diatemp +24.0_realk*temp)
    temp=0.0_realk
    temp=diatemp*factor
    call mat_dmul(temp,KURT%carmomMO(x),'T',dble(m),1.0_realk,G)


    diatemp=0.0_realk;diatemp2=0.0_realk;diatemp3=0.0_realk;temp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(y))
    call mat_extract_diagonal(diatemp2,KURT%carmomMO(yy))
    call mat_extract_diagonal(diatemp3,KURT%carmomMO(yyy))
    temp=diatemp*diatemp*diatemp
    diatemp= (8.0_realk*diatemp3-24.0_realk*diatemp2*diatemp +24.0_realk*temp)
    temp=0.0_realk
    temp=diatemp*factor
    call mat_dmul(temp,KURT%carmomMO(y),'T',1.0_realk*dble(m),1.0_realk,G)



    diatemp=0.0_realk;diatemp2=0.0_realk;diatemp3=0.0_realk;temp=0.0_realk
    call mat_extract_diagonal(diatemp,KURT%carmomMO(z))
    call mat_extract_diagonal(diatemp2,KURT%carmomMO(zz))
    call mat_extract_diagonal(diatemp3,KURT%carmomMO(zzz))
    temp=diatemp*diatemp*diatemp
    diatemp= (8.0_realk*diatemp3-24.0_realk*diatemp2*diatemp +24.0_realk*temp)
    temp=0.0_realk
    temp=diatemp*factor
    call mat_dmul(temp,KURT%carmomMO(z),'T',1.0_realk*dble(m),1.0_realk,G)


    if (KURT%crossterms) then
       call mat_dmul(factor,KURT%carmomMO(xxyy),'T',-dble(4*m),1.0_realk,G)

       call mat_dmul(factor,KURT%carmomMO(xxzz),'T',-dble(4*m),1.0_realk,G)

       call mat_dmul(factor,KURT%carmomMO(yyzz),'T',-dble(4*m),1.0_realk,G)


       call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
       call mat_extract_diagonal(diatemp2,KURT%carmomMO(y))
       diatemp3=diatemp*diatemp2*factor
       call mat_dmul(diatemp3,KURT%carmomMO(xy),'T',-dble(16*m),1.0_realk,G)


       call mat_extract_diagonal(diatemp3,KURT%carmomMO(z))
       call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
       diatemp2=diatemp*diatemp3*factor
       call mat_dmul(diatemp2,KURT%carmomMO(xz),'T',-dble(16*m),1.0_realk,G)


       call mat_extract_diagonal(diatemp3,KURT%carmomMO(y))
       call mat_extract_diagonal(diatemp2,KURT%carmomMO(z))
       diatemp=diatemp3*diatemp2*factor 
       call mat_dmul(diatemp,KURT%carmomMO(yz),'T',-dble(16*m),1.0_realk,G)

       call GradTerm_xj(KURT,factor,x,y,yy,xy,xyy,G,norb)

       call GradTerm_xj(KURT,factor,y,x,xx,xy,xxy,G,norb)

       call GradTerm_xj(KURT,factor,x,z,zz,xz,xzz,G,norb)

       call GradTerm_xj(KURT,factor,z,x,xx,xz,xxz,G,norb)

       call GradTerm_xj(KURT,factor,y,z,zz,yz,yzz,G,norb)

       call GradTerm_xj(KURT,factor,z,y,yy,yz,yyz,G,norb)

       call mat_extract_diagonal(diatemp,KURT%carmomMO(y))
       diatemp2=diatemp*diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(xx),'T',-dble(4*m),1.0_realk,G)

       call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
       diatemp2=diatemp*diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(yy),'T',-dble(4*m),1.0_realk,G)


       call mat_extract_diagonal(diatemp,KURT%carmomMO(z))
       diatemp2=diatemp*diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(xx),'T',-dble(4*m),1.0_realk,G)


       call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
       diatemp2=diatemp*diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(zz),'T',-dble(4*m),1.0_realk,G)


       call mat_extract_diagonal(diatemp,KURT%carmomMO(z))
       diatemp2=diatemp*diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(yy),'T',-dble(4*m),1.0_realk,G)


       call mat_extract_diagonal(diatemp,KURT%carmomMO(y))
       diatemp2=diatemp*diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(zz),'T',-dble(4*m),1.0_realk,G)



       call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
       diatemp2=diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(xyy),'T',dble(8*m),1.0_realk,G) 

       call mat_extract_diagonal(diatemp,KURT%carmomMO(y))
       diatemp2=diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(xxy),'T',dble(8*m),1.0_realk,G) 


       call mat_extract_diagonal(diatemp,KURT%carmomMO(x))
       diatemp2=diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(xzz),'T',dble(8*m),1.0_realk,G) 

       call mat_extract_diagonal(diatemp,KURT%carmomMO(z))
       diatemp2=diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(xxz),'T',dble(8*m),1.0_realk,G) 

       call mat_extract_diagonal(diatemp,KURT%carmomMO(y))
       diatemp2=diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(yzz),'T',dble(8*m),1.0_realk,G) 

       call mat_extract_diagonal(diatemp,KURT%carmomMO(z))
       diatemp2=diatemp*factor
       call mat_dmul(diatemp2,KURT%carmomMO(yyz),'T',dble(8*m),1.0_realk,G) 
    endif


    call mat_init(Gt,G%ncol,G%nrow)
    call mat_trans(G,Gt)
    call mat_daxpy(-1.0_realk,Gt,G)
    call mat_free(Gt)

  end subroutine compute_gradient

  subroutine GradTerm_xj(KURT,factor,xj,xi,xixi,xixj,xixixj,G,norb)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: G
    integer :: norb,xi,xj,xixi,xixj,xixixj
    real(realk),intent(in) :: factor(norb) 
    real(realk) :: diatemp(norb),diatemp2(norb),diatemp3(norb)

    diatemp=0.0_realk
    call mat_extract_diagonal(diatemp3,KURT%carmomMO(xixi))
    call mat_extract_diagonal(diatemp2,KURT%carmomMO(xj))
    diatemp=-8.0_realk*diatemp3*diatemp2
    call mat_extract_diagonal(diatemp2,KURT%carmomMO(xj))
    call mat_extract_diagonal(diatemp3,KURT%carmomMO(xi))
    diatemp=diatemp +24.0_realk*diatemp3*diatemp3*diatemp2
    call mat_extract_diagonal(diatemp2,KURT%carmomMO(xixixj))
    diatemp=diatemp +8.0_realk*diatemp2
    call mat_extract_diagonal(diatemp2,KURT%carmomMO(xixj))
    call mat_extract_diagonal(diatemp3,KURT%carmomMO(xi))
    diatemp=diatemp-16.0_realk*diatemp2*diatemp3 
    diatemp=diatemp*factor
    call mat_dmul(diatemp,KURT%carmomMO(xj),'n',dble(KURT%m),1.0_realk,G)


  end subroutine GradTerm_xj

!> \brief Compute carmom in AO basis
  subroutine compute_carmom(KURT,ls)
    implicit none
    type(PFMitem) :: KURT
    type(lsitem) :: ls
    type(matrix) :: cmo
    integer :: nderiv,nmat
    integer :: i,indx,mat_indx(24)

    nderiv = 4
    nmat=(nderiv+1)*(nderiv+2)*(nderiv+3)/6

    ! Map essential matrix numbers to a vector
    do i=2,14
       mat_indx(i-1)=i
    end do
    do i=16,21
       mat_indx(i-2)=i
    end do
    mat_indx(20)=24
    mat_indx(21)=26
    mat_indx(22)=31
    mat_indx(23)=33
    mat_indx(24)=35


    do i=1,24
       call II_get_single_carmom(6,6,ls%setting,KURT%carmomAO(i),mat_indx(i),nderiv,0.0_realk,0.0_realk,0.0_realk)
    end do

  end subroutine compute_carmom

!> \brief Transform from AO basis to MO basis
!> \param carmomAO  (input)
!> \param CMO (input)
!> \param carmomMO (output)
  subroutine carmomAO_to_carmomMO(carmomAO,cmo,carmomMO)
    implicit none
    type(matrix) :: carmomAO,carmomMO,cmo
    type(matrix) :: temp

    call mat_init(temp,cmo%ncol,carmomAO%ncol)
    call mat_mul(cmo,carmomAO,'T','n',1.0_realk,0.0_realk,temp)
    call mat_mul(temp,cmo,'n','n',1.0_realk,0.0_realk,carmomMO)
    call mat_free(temp)

  end subroutine carmomAO_to_carmomMO

  subroutine compute_omega(mat,omega,crossterms,n)
    implicit none
    type(matrix),intent(in) :: mat(24) 
    integer :: n   !Number of orbitals
    logical :: crossterms
    real(realk) :: temp(n),temp2(n),temp3(n)
    real(realk) :: cross_total(n),total(n),omega(n)
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz


    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24



    total=0.0_realk
    call mat_extract_diagonal(temp,mat(xxxx))
    total=temp
    call mat_extract_diagonal(temp,mat(yyyy))
    total = total +temp
    call mat_extract_diagonal(temp,mat(zzzz))
    total = total +temp

    call mat_extract_diagonal(temp,mat(xxx))
    call mat_extract_diagonal(temp2,mat(x))
    temp3=temp*temp2
    call mat_extract_diagonal(temp,mat(yyy))
    call mat_extract_diagonal(temp2,mat(y))
    temp3=temp3+temp*temp2
    call mat_extract_diagonal(temp,mat(zzz))
    call mat_extract_diagonal(temp2,mat(z))
    temp3=temp3+temp*temp2
    total = total - 4.0_realk*temp3

    temp3 = 0.0_realk
    call mat_extract_diagonal(temp,mat(xx))
    call mat_extract_diagonal(temp2,mat(x))
    temp3=temp*temp2*temp2
    call mat_extract_diagonal(temp,mat(yy))
    call mat_extract_diagonal(temp2,mat(y))
    temp3=temp3+temp*temp2*temp2
    call mat_extract_diagonal(temp,mat(zz))
    call mat_extract_diagonal(temp2,mat(z))
    temp3=temp3+temp*temp2*temp2
    total=total+6.0_realk*temp3

    call mat_extract_diagonal(temp,mat(x))
    total= total -3.0_realk*temp*temp*temp*temp

    call mat_extract_diagonal(temp,mat(y))
    total= total -3.0_realk*temp*temp*temp*temp

    call mat_extract_diagonal(temp,mat(z))
    total= total -3.0_realk*temp*temp*temp*temp


    if (crossterms) then
       call mat_extract_diagonal(temp,mat(xxyy))
       cross_total = 2.0_realk*temp
       call mat_extract_diagonal(temp,mat(xxzz))
       cross_total = cross_total+2.0_realk*temp
       call mat_extract_diagonal(temp,mat(yyzz))
       cross_total = cross_total+2.0_realk*temp

       call mat_extract_diagonal(temp,mat(xx))
       call mat_extract_diagonal(temp2,mat(y))
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       call mat_extract_diagonal(temp,mat(yy))
       call mat_extract_diagonal(temp2,mat(x))
       cross_total=cross_total+2.0_realk*temp*temp2*temp2


       call mat_extract_diagonal(temp,mat(xx))
       call mat_extract_diagonal(temp2,mat(z))
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       call mat_extract_diagonal(temp,mat(zz))
       call mat_extract_diagonal(temp2,mat(x))
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       call mat_extract_diagonal(temp,mat(yy))
       call mat_extract_diagonal(temp2,mat(z))
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       call mat_extract_diagonal(temp,mat(zz))
       call mat_extract_diagonal(temp2,mat(y))
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       call mat_extract_diagonal(temp,mat(xxy))
       call mat_extract_diagonal(temp2,mat(y))
       cross_total = cross_total -4.0_realk*temp*temp2

       call mat_extract_diagonal(temp,mat(xyy))
       call mat_extract_diagonal(temp2,mat(x))
       cross_total = cross_total -4.0_realk*temp*temp2

       call mat_extract_diagonal(temp,mat(xxz))
       call mat_extract_diagonal(temp2,mat(z))
       cross_total = cross_total -4.0_realk*temp*temp2

       call mat_extract_diagonal(temp,mat(xzz))
       call mat_extract_diagonal(temp2,mat(x))
       cross_total = cross_total -4.0_realk*temp*temp2

       call mat_extract_diagonal(temp,mat(yyz))
       call mat_extract_diagonal(temp2,mat(z))
       cross_total = cross_total -4.0_realk*temp*temp2

       call mat_extract_diagonal(temp,mat(yzz))
       call mat_extract_diagonal(temp2,mat(y))
       cross_total = cross_total -4.0_realk*temp*temp2

       call mat_extract_diagonal(temp,mat(x))
       call mat_extract_diagonal(temp2,mat(y))
       cross_total = cross_total-6.0_realk*temp*temp*temp2*temp2

       call mat_extract_diagonal(temp,mat(x))
       call mat_extract_diagonal(temp2,mat(z))
       cross_total = cross_total-6.0_realk*temp*temp*temp2*temp2


       call mat_extract_diagonal(temp,mat(y))
       call mat_extract_diagonal(temp2,mat(z))
       cross_total = cross_total-6.0_realk*temp*temp*temp2*temp2

       call mat_extract_diagonal(temp,mat(x))
       call mat_extract_diagonal(temp2,mat(y))
       call mat_extract_diagonal(temp3,mat(xy))
       cross_total =cross_total+8.0_realk*temp*temp2*temp3


       call mat_extract_diagonal(temp,mat(x))
       call mat_extract_diagonal(temp2,mat(z))
       call mat_extract_diagonal(temp3,mat(xz))
       cross_total =cross_total+8.0_realk*temp*temp2*temp3

       call mat_extract_diagonal(temp,mat(y))
       call mat_extract_diagonal(temp2,mat(z))
       call mat_extract_diagonal(temp3,mat(yz))
       cross_total =cross_total+8.0_realk*temp*temp2*temp3

       total=total+cross_total
    end if

    omega=total

  end subroutine compute_omega

  subroutine compute_memreduced_noAOinit_omega(KURT,cmo,nbas,ls)
    implicit none
    type(PFMitem) :: KURT
    type(matrix),intent(in) :: cmo
    integer,intent(in) :: nbas
    type(lsitem) :: ls
    !
    real(realk),pointer :: matdiag(:,:)
    integer :: I,nderiv,nmat,mat_indx(24)
    type(matrix) :: carmomMO,temp,matAO

    call mem_alloc(KURT%omega,KURT%norb)
    nderiv = 4
    nmat=(nderiv+1)*(nderiv+2)*(nderiv+3)/6
    ! Map essential matrix numbers to a vector
    do i=2,14
       mat_indx(i-1)=i
    end do
    do i=16,21
       mat_indx(i-2)=i
    end do
    mat_indx(20)=24
    mat_indx(21)=26
    mat_indx(22)=31
    mat_indx(23)=33
    mat_indx(24)=35

    call mem_alloc(matdiag,KURT%norb,24)
    call mat_init(matAO,nbas,nbas)
    call mat_init(carmomMO,cmo%ncol,cmo%ncol)
    call mat_init(temp,cmo%ncol,nbas)
    do I=1,24
       call II_get_single_carmom(6,6,ls%setting,matAO,mat_indx(I),nderiv,&
            & 0.0_realk,0.0_realk,0.0_realk)
       call mat_mul(cmo,matAO,'T','n',1.0_realk,0.0_realk,temp)
       call mat_mul(temp,cmo,'n','n',1.0_realk,0.0_realk,carmomMO)
       call mat_extract_diagonal(matdiag(:,I),carmomMO)
    enddo
    call mat_free(temp)
    call mat_free(carmomMO)
    call mat_free(matAO)
    call compute_omega_diagonal(matdiag,KURT%omega,KURT%crossterms,KURT%norb)
    call mem_dealloc(matdiag)

  end subroutine compute_memreduced_noAOinit_omega

  subroutine compute_memreduced_omega(matAO,omega,crossterms,n,cmo)
    implicit none
    type(matrix),intent(in) :: matAO(24) 
    real(realk),intent(inout) :: omega(n)
    type(matrix),intent(in) :: cmo
    integer :: n   !Number of orbitals
    logical :: crossterms
    real(realk),pointer :: matdiag(:,:)
    integer :: I
    type(matrix) :: carmomMO,temp
    call mem_alloc(matdiag,n,24)
    call mat_init(carmomMO,cmo%ncol,cmo%ncol)
    call mat_init(temp,cmo%ncol,matAO(1)%ncol)
    do I=1,24
       call mat_mul(cmo,matAO(I),'T','n',1.0_realk,0.0_realk,temp)
       call mat_mul(temp,cmo,'n','n',1.0_realk,0.0_realk,carmomMO)
       call mat_extract_diagonal(matdiag(:,I),carmomMO)
    enddo
    call mat_free(temp)
    call mat_free(carmomMO)
    call compute_omega_diagonal(matdiag,omega,crossterms,n)
    call mem_dealloc(matdiag)

  end subroutine compute_memreduced_omega

  subroutine compute_omega_diagonal(matdiag,omega,crossterms,n)
    implicit none
    integer,intent(in) :: n   !Number of orbitals
    real(realk),intent(in) :: matdiag(n,24) 
    logical :: crossterms
    real(realk) :: temp(n),temp2(n),temp3(n)
    real(realk) :: cross_total(n),total(n),omega(n)
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz

    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24

    total=0.0_realk
    temp=matdiag(:,xxxx)
    total=temp
    temp=matdiag(:,yyyy)
    total = total +temp
    temp=matdiag(:,zzzz)
    total = total +temp

    temp=matdiag(:,xxx)
    temp2=matdiag(:,x)
    temp3=temp*temp2
    temp = matdiag(:,yyy)
    temp2 = matdiag(:,y)
    temp3=temp3+temp*temp2
    temp = matdiag(:,zzz)
    temp2 = matdiag(:,z)
    temp3=temp3+temp*temp2
    total = total - 4.0_realk*temp3

    temp3 = 0.0_realk
    temp = matdiag(:,xx)
    temp2 = matdiag(:,x)
    temp3=temp*temp2*temp2
    temp = matdiag(:,yy)
    temp2 = matdiag(:,y)
    temp3=temp3+temp*temp2*temp2
    temp = matdiag(:,zz)
    temp2 = matdiag(:,z)
    temp3=temp3+temp*temp2*temp2
    total=total+6.0_realk*temp3

    temp = matdiag(:,x)
    total= total -3.0_realk*temp*temp*temp*temp

    temp = matdiag(:,y)
    total= total -3.0_realk*temp*temp*temp*temp

    temp = matdiag(:,z)
    total= total -3.0_realk*temp*temp*temp*temp


    if (crossterms) then
       temp = matdiag(:,xxyy)
       cross_total = 2.0_realk*temp
       temp = matdiag(:,xxzz)
       cross_total = cross_total+2.0_realk*temp
       temp = matdiag(:,yyzz)
       cross_total = cross_total+2.0_realk*temp

       temp = matdiag(:,xx)
       temp2 = matdiag(:,y)
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       temp = matdiag(:,yy)
       temp2 = matdiag(:,x)
       cross_total=cross_total+2.0_realk*temp*temp2*temp2


       temp = matdiag(:,xx)
       temp2 = matdiag(:,z)
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       temp = matdiag(:,zz)
       temp2 = matdiag(:,x)
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       temp = matdiag(:,yy)
       temp2 = matdiag(:,z)
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       temp = matdiag(:,zz)
       temp2 = matdiag(:,y)
       cross_total=cross_total+2.0_realk*temp*temp2*temp2

       temp = matdiag(:,xxy)
       temp2 = matdiag(:,y)
       cross_total = cross_total -4.0_realk*temp*temp2

       temp = matdiag(:,xyy)
       temp2 = matdiag(:,x)
       cross_total = cross_total -4.0_realk*temp*temp2

       temp = matdiag(:,xxz)
       temp2 = matdiag(:,z)
       cross_total = cross_total -4.0_realk*temp*temp2

       temp = matdiag(:,xzz)
       temp2 = matdiag(:,x)
       cross_total = cross_total -4.0_realk*temp*temp2

       temp = matdiag(:,yyz)
       temp2 = matdiag(:,z)
       cross_total = cross_total -4.0_realk*temp*temp2

       temp = matdiag(:,yzz)
       temp2 = matdiag(:,y)
       cross_total = cross_total -4.0_realk*temp*temp2

       temp = matdiag(:,x)
       temp2 = matdiag(:,y)
       cross_total = cross_total-6.0_realk*temp*temp*temp2*temp2

       temp = matdiag(:,x)
       temp2 = matdiag(:,z)
       cross_total = cross_total-6.0_realk*temp*temp*temp2*temp2


       temp = matdiag(:,y)
       temp2 = matdiag(:,z)
       cross_total = cross_total-6.0_realk*temp*temp*temp2*temp2

       temp = matdiag(:,x)
       temp2 = matdiag(:,y)
       temp3 = matdiag(:,xy)
       cross_total =cross_total+8.0_realk*temp*temp2*temp3


       temp = matdiag(:,x)
       temp2 = matdiag(:,z)
       temp3 = matdiag(:,xz)
       cross_total =cross_total+8.0_realk*temp*temp2*temp3

       temp = matdiag(:,y)
       temp2 = matdiag(:,z)
       temp3 = matdiag(:,yz)
       cross_total =cross_total+8.0_realk*temp*temp2*temp3

       total=total+cross_total
    end if

    omega=total

  end subroutine compute_omega_diagonal

!> \brief Compute linear transformations Hkappa
!> \param kappa  trial vector (input)
!> \param G gradient (input)
!> \param Hkappa linear transformations (output)
  subroutine compute_LinTrans(KURT,kappa,G,HK)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: kappa
    type(matrix) :: HK,G
    type(matrix) :: temp,tempT
    real(realk)  :: omen(HK%ncol),KF(HK%ncol),omto(HK%ncol)
    real(realk)  :: dia(HK%ncol),omega(HK%ncol)
    real(realk)  :: mat_sqnrm2
    integer :: m,i
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz

    m=KURT%m

    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24

    if (m > 1) then
       do i=1,KURT%norb
          omen(i)=KURT%omega(i)**(m-1)
       end do
    else 
       omen=1.0_realk
    end if

    KF = 0.0_realk
    call mat_zero(HK)
    call mat_init(temp,kappa%nrow,kappa%ncol)
    call mat_init(tempT,kappa%nrow,kappa%ncol)

    call KXXXX_lintrans(kappa,KURT%carmomMO(xxxx),omen,m,HK,KF)
    call KXXXX_lintrans(kappa,KURT%carmomMO(yyyy),omen,m,HK,KF)
    call KXXXX_lintrans(kappa,KURT%carmomMO(zzzz),omen,m,HK,KF)

    call KXXX_lintrans(kappa,KURT%carmomMO(xxx),KURT%carmomMO(x),omen,m,HK,KF)
    call KXXX_lintrans(kappa,KURT%carmomMO(yyy),KURT%carmomMO(y),omen,m,HK,KF)
    call KXXX_lintrans(kappa,KURT%carmomMO(zzz),KURT%carmomMO(z),omen,m,HK,KF)

    call KXX_lintrans(kappa,KURT,omen,m,HK,.true.,.false.,.false.,KF)
    call KXX_lintrans(kappa,KURT,omen,m,HK,.false.,.true.,.false.,KF)
    call KXX_lintrans(kappa,KURT,omen,m,HK,.false.,.false.,.true.,KF)

    call KX_lintrans(kappa,KURT,omen,m,HK,.true.,.false.,.false.,KF)
    call KX_lintrans(kappa,KURT,omen,m,HK,.false.,.true.,.false.,KF)
    call KX_lintrans(kappa,KURT,omen,m,HK,.false.,.false.,.true.,KF)


    if (KURT%crossterms) then
       call KXjXjXi_lintrans(kappa,KURT%carmomMO(xxy),KURT%carmomMO(y),omen,m,HK,KF)
       call KXjXjXi_lintrans(kappa,KURT%carmomMO(xyy),KURT%carmomMO(x),omen,m,HK,KF)

       call KXjXjXi_lintrans(kappa,KURT%carmomMO(xxz),KURT%carmomMO(z),omen,m,HK,KF)
       call KXjXjXi_lintrans(kappa,KURT%carmomMO(xzz),KURT%carmomMO(x),omen,m,HK,KF)

       call KXjXjXi_lintrans(kappa,KURT%carmomMO(yyz),KURT%carmomMO(z),omen,m,HK,KF)
       call KXjXjXi_lintrans(kappa,KURT%carmomMO(yzz),KURT%carmomMO(y),omen,m,HK,KF)

       call KXjXi_lintrans(kappa,KURT%carmomMO(xy),KURT%carmomMO(x),KURT%carmomMO(y),omen,m,HK,KF)
       call KXjXi_lintrans(kappa,KURT%carmomMO(xy),KURT%carmomMO(y),KURT%carmomMO(x),omen,m,HK,KF)

       call KXjXi_lintrans(kappa,KURT%carmomMO(xz),KURT%carmomMO(x),KURT%carmomMO(z),omen,m,HK,KF)
       call KXjXi_lintrans(kappa,KURT%carmomMO(xz),KURT%carmomMO(z),KURT%carmomMO(x),omen,m,HK,KF)

       call KXjXi_lintrans(kappa,KURT%carmomMO(yz),KURT%carmomMO(y),KURT%carmomMO(z),omen,m,HK,KF)
       call KXjXi_lintrans(kappa,KURT%carmomMO(yz),KURT%carmomMO(z),KURT%carmomMO(y),omen,m,HK,KF)

       call noperm_lintrans(kappa,KURT%carmomMO(xxyy),KURT%carmomMO(xy),KURT%carmomMO(x),KURT%carmomMO(y),omen,m,HK,KF)
       call noperm_lintrans(kappa,KURT%carmomMO(xxzz),KURT%carmomMO(xz),KURT%carmomMO(x),KURT%carmomMO(z),omen,m,HK,KF)
       call noperm_lintrans(kappa,KURT%carmomMO(yyzz),KURT%carmomMO(yz),KURT%carmomMO(y),KURT%carmomMO(z),omen,m,HK,KF)
    end if

    if (m>1) then
       if (m==2) then
          omto=1.0_realk
       else
          do i=1,KURT%norb
             omto(i)=KURT%omega(i)**(m-2)
          end do
       end if
       KF = dble(m-1)*omto*KF
    end if

    ! HK terms for m>1
    if (m>1) then
       call mat_zero(temp)
       call mat_daxpy(1.0_realk,KURT%carmomMO(xxxx),temp) 
       call mat_daxpy(1.0_realk,KURT%carmomMO(yyyy),temp) 
       call mat_daxpy(1.0_realk,KURT%carmomMO(zzzz),temp) 
       call mat_dmul(KF,temp,'n',-dble(2*m),1.0_realk,HK)

       call xxx_KF(KURT%carmomMO(xxx),KURT%carmomMO(x),KF,m,HK)
       call xxx_KF(KURT%carmomMO(yyy),KURT%carmomMO(y),KF,m,HK)
       call xxx_KF(KURT%carmomMO(zzz),KURT%carmomMO(z),KF,m,HK)

       call xx_KF(KURT%carmomMO(xx),KURT%carmomMO(x),KF,m,HK)
       call xx_KF(KURT%carmomMO(yy),KURT%carmomMO(y),KF,m,HK)
       call xx_KF(KURT%carmomMO(zz),KURT%carmomMO(z),KF,m,HK)

       call x_KF(KURT%carmomMO(xxx),KURT%carmomMO(xx),KURT%carmomMO(x),KF,m,HK)
       call x_KF(KURT%carmomMO(yyy),KURT%carmomMO(yy),KURT%carmomMO(y),KF,m,HK)
       call x_KF(KURT%carmomMO(zzz),KURT%carmomMO(zz),KURT%carmomMO(z),KF,m,HK)

       cross: if (KURT%crossterms) then
          call noperm_KF(KURT%carmomMO(xxyy),KURT%carmomMO(xy),KURT%carmomMO(x),KURT%carmomMO(y),m,KF,HK)
          call noperm_KF(KURT%carmomMO(xxzz),KURT%carmomMO(xz),KURT%carmomMO(x),KURT%carmomMO(z),m,KF,HK)
          call noperm_KF(KURT%carmomMO(yyzz),KURT%carmomMO(yz),KURT%carmomMO(y),KURT%carmomMO(z),m,KF,HK)

          call xj_KF(KURT%carmomMO(x),KURT%carmomMO(y),KURT%carmomMO(yy),KURT%carmomMO(xy),KURT%carmomMO(xyy),m,KF,HK)
          call xj_KF(KURT%carmomMO(y),KURT%carmomMO(x),KURT%carmomMO(xx),KURT%carmomMO(xy),KURT%carmomMO(xxy),m,KF,HK)

          call xj_KF(KURT%carmomMO(x),KURT%carmomMO(z),KURT%carmomMO(zz),KURT%carmomMO(xz),KURT%carmomMO(xzz),m,KF,HK)
          call xj_KF(KURT%carmomMO(z),KURT%carmomMO(x),KURT%carmomMO(xx),KURT%carmomMO(xz),KURT%carmomMO(xxz),m,KF,HK)

          call xj_KF(KURT%carmomMO(y),KURT%carmomMO(z),KURT%carmomMO(zz),KURT%carmomMO(yz),KURT%carmomMO(yzz),m,KF,HK)
          call xj_KF(KURT%carmomMO(z),KURT%carmomMO(y),KURT%carmomMO(yy),KURT%carmomMO(yz),KURT%carmomMO(yyz),m,KF,HK)

          call xjxj_KF(KURT%carmomMO(xx),KURT%carmomMO(y),m,KF,HK)
          call xjxj_KF(KURT%carmomMO(yy),KURT%carmomMO(x),m,KF,HK)

          call xjxj_KF(KURT%carmomMO(xx),KURT%carmomMO(z),m,KF,HK)
          call xjxj_KF(KURT%carmomMO(zz),KURT%carmomMO(x),m,KF,HK)

          call xjxj_KF(KURT%carmomMO(yy),KURT%carmomMO(z),m,KF,HK)
          call xjxj_KF(KURT%carmomMO(zz),KURT%carmomMO(y),m,KF,HK)

          call xixixj_KF(KURT%carmomMO(xxy),KURT%carmomMO(y),m,KF,HK)
          call xixixj_KF(KURT%carmomMO(xyy),KURT%carmomMO(x),m,KF,HK)

          call xixixj_KF(KURT%carmomMO(xxz),KURT%carmomMO(z),m,KF,HK)
          call xixixj_KF(KURT%carmomMO(xzz),KURT%carmomMO(x),m,KF,HK)

          call xixixj_KF(KURT%carmomMO(yyz),KURT%carmomMO(z),m,KF,HK)
          call xixixj_KF(KURT%carmomMO(yzz),KURT%carmomMO(y),m,KF,HK)
       end if cross
    end if
    call mat_mul(G,kappa,'T','T',0.5_realk,1.0_realk,HK)
    !call mat_mul(kappa,G,'n','n',-0.5_realk,0.0_realk,temp)


    call mat_trans(HK,temp)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,temp,HK)

    call mat_free(temp);call mat_free(tempT)
  end subroutine compute_LinTrans


  subroutine xixixj_KF(xixixj,xj,m,KF,HK)
    implicit none
    type(matrix) :: xixixj,xj
    type(matrix) :: HK
    real(realk)  :: KF(xj%ncol)
    real(realk)  :: dia(xj%ncol)
    integer :: m

    call mat_extract_diagonal(dia,xj)
    dia = dia*KF
    call mat_dmul(dia,xixixj,'T',dble(8*m),1.0_realk,HK)

  end subroutine xixixj_KF



  subroutine xjxj_KF(xixi,xj,m,KF,HK)
    implicit none
    type(matrix) :: xixi,xj,HK
    real(realk)  :: KF(xj%ncol)
    real(realk)  :: dia(xj%ncol)
    integer :: m

    call mat_extract_diagonal(dia,xj)
    dia = dia*dia*KF
    call mat_dmul(dia,xixi,'T',-dble(4*m),1.0_realk,HK)

  end subroutine xjxj_KF



  subroutine xj_KF(xj,xi,xixi,xixj,xixixj,m,KF,HK) 
    implicit none
    type(matrix) :: xj,xi,xixi,xixj,xixixj 
    type(matrix) :: HK
    real(realk)  :: KF(xi%ncol)
    real(realk)  :: dia(xi%ncol)
    real(realk)  :: dia2(xi%ncol)
    real(realk)  :: dia3(xi%ncol)
    integer :: m


    call mat_extract_diagonal(dia,xixi)
    call mat_extract_diagonal(dia2,xj)
    dia = -4.0_realk*dia*dia2
    call mat_extract_diagonal(dia3,xi)
    dia =dia +12.0_realk*dia3*dia3*dia2
    call mat_extract_diagonal(dia2,xixj)
    dia =dia -8.0_realk*dia3*dia2
    call mat_extract_diagonal(dia2,xixixj)
    dia = dia+4.0_realk*dia2
    dia = dia*KF
    call mat_dmul(dia,xj,'T',dble(2*m),1.0_realk,HK)


  end subroutine xj_KF



  subroutine noperm_KF(xixixjxj,xixj,xi,xj,m,KF,HK)
    implicit none
    type(matrix) ::  xixixjxj,xixj,xi,xj
    type(matrix) :: HK
    real(realk)  ::KF(xj%ncol)
    real(realk)  ::dia(xj%ncol)
    real(realk)  ::dia2(xj%ncol)
    integer :: m

    call mat_dmul(KF,xixixjxj,'T',-dble(4*m),1.0_realk,HK)

    call mat_extract_diagonal(dia,xi)
    call mat_extract_diagonal(dia2,xj)
    dia = dia*dia2*KF
    call mat_dmul(dia,xixj,'T',-dble(16*m),1.0_realk,HK)

  end subroutine noperm_KF


  subroutine KXjXi_lintrans(k,xjxi,xj,xi,omen,m,HK,KF)
    implicit none
    type(matrix) :: xjxi,xj,xi,k,HK
    type(matrix) :: temp,tempt
    real(realk)  :: omen(k%ncol)
    real(realk)  :: KF(k%ncol)
    real(realk)  :: dia(k%ncol)
    real(realk)  :: dia2(k%ncol)
    integer :: m

    call mat_init(temp,k%ncol,k%ncol)
    call mat_init(tempt,k%ncol,k%ncol)
    call mat_mul(k,xjxi,'n','n',1.0_realk,0.0_realk,temp)
    call mat_trans(temp,tempt)
    call mat_daxpy(1.0_realk,tempt,temp)
    call mat_free(tempt)

    call mat_extract_diagonal(dia2,temp)
    call mat_extract_diagonal(dia,xi)
    dia=dia*dia2*omen
    call mat_dmul(dia,xj,'T',-dble(16*m),1.0_realk,HK)


    call mat_free(temp)

  end subroutine KXjXi_lintrans



  subroutine noperm_lintrans(k,xixixjxj,xixj,xi,xj,omen,m,HK,KF)
    implicit none
    type(matrix) :: xixixjxj,xixj,k,HK
    type(matrix) :: xi,xj
    type(matrix) :: temp,tempt
    real(realk)  :: omen(k%ncol)
    real(realk)  :: KF(k%ncol)
    real(realk)  :: dia(k%ncol)
    real(realk)  :: dia2(k%ncol)
    real(realk)  :: dia3(k%ncol)
    integer :: m

    call mat_init(temp,k%nrow,k%ncol)
    call mat_init(tempt,k%nrow,k%ncol)
    call mat_mul(k,xixixjxj,'n','n',1.0_realk,0.0_realk,temp)
    if (m>1) then
       call mat_extract_diagonal(dia,temp)
       KF=KF+4.0_realk*dia
    end if
    call mat_trans(temp,tempt)
    call mat_daxpy(1.0_realk,tempt,temp)
    call mat_dmul(omen,temp,'T',-dble(4*m),1.0_realk,HK)
    call mat_zero(temp);call mat_zero(tempt)

    call mat_mul(k,xixj,'n','n',1.0_realk,0.0_realk,temp)
    if (m>1) then
       call mat_extract_diagonal(dia,temp)
       call mat_extract_diagonal(dia2,xi)
       call mat_extract_diagonal(dia3,xj)
       KF = KF +16.0_realk*dia*dia2*dia3
    end if
    call mat_trans(temp,tempt)
    call mat_daxpy(1.0_realk,tempt,temp)
    call mat_free(tempt)
    call mat_extract_diagonal(dia,xi)
    call mat_extract_diagonal(dia2,xj)
    dia = dia*dia2*omen
    call mat_dmul(dia,temp,'T',-dble(16*m),1.0_realk,HK)
    call mat_free(temp)

  end subroutine noperm_lintrans

  subroutine KXjXjXi_lintrans(k,xjxjxi,xi,omen,m,HK,KF)
    implicit none
    type(matrix) :: xjxjxi,xi,k,HK
    type(matrix) :: temp,tempt
    real(realk)  :: omen(k%ncol)
    real(realk)  :: KF(k%ncol)
    real(realk)  :: dia(k%ncol)
    real(realk)  :: dia2(k%ncol)
    integer :: m

    call mat_init(temp,k%nrow,k%ncol)
    call mat_init(tempt,k%nrow,k%ncol)
    call mat_mul(k,xjxjxi,'n','n',1.0_realk,0.0_realk,temp)
    ! ** UPDATE KF
    if (m>1) then
       call mat_extract_diagonal(dia,temp)   
       call mat_extract_diagonal(dia2,xi)
       KF = KF-8.0_realk*dia*dia2 
    end if
    ! ** END UPDATE KF
    call mat_trans(temp,tempt)
    call mat_daxpy(1.0_realk,tempt,temp)
    call mat_free(tempt)


    call mat_extract_diagonal(dia,xi)
    dia=dia*omen
    call mat_dmul(dia,temp,'T',dble(8*m),1.0_realk,HK)

    call mat_extract_diagonal(dia,temp)
    dia=dia*omen
    call mat_dmul(dia,xi,'T',dble(8*m),1.0_realk,HK)


    call mat_free(temp)
  end subroutine KXjXjXi_lintrans


  subroutine x_KF(xxx,xx,x,KF,m,HK)
    implicit none
    type(matrix),intent(in) :: xxx,xx,x
    type(matrix) :: HK
    real(realk)  :: KF(HK%ncol)
    real(realk)  :: dia(HK%ncol),dia2(HK%ncol),dia3(HK%ncol)
    integer :: m

    call mat_extract_diagonal(dia2,x)
    call mat_extract_diagonal(dia3,xx)
    dia = 24.0_realk*dia2*dia2*dia2-24.0_realk*dia3*dia2 
    call mat_extract_diagonal(dia2,xxx)
    dia =dia+8.0_realk*dia2
    dia = KF*dia
    call mat_dmul(dia,x,'T',dble(m),1.0_realk,HK)

  end subroutine x_KF


  subroutine xxx_KF(xxx,x,KF,m,HK)
    implicit none
    type(matrix),intent(in) :: xxx,x
    type(matrix) :: HK
    real(realk)  :: KF(HK%ncol),dia(HK%ncol)
    integer :: m

    call mat_extract_diagonal(dia,x)
    dia = dia*KF
    call mat_dmul(dia,xxx,'T',dble(8*m),1.0_realk,HK)

  end subroutine xxx_KF

  subroutine xx_KF(xx,x,KF,m,HK)
    implicit none
    type(matrix),intent(in) :: xx,x
    type(matrix) :: HK
    real(realk)  :: KF(HK%ncol),dia(HK%ncol)
    integer      :: m

    call mat_extract_diagonal(dia,x)
    dia = KF*dia*dia
    call mat_dmul(dia,xx,'T',-dble(12*m),1.0_realk,HK)

  end subroutine xx_KF



  subroutine KX_lintrans(k,KURT,omen,m,HK,logx,logy,logz,KF)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: K,HK
    real(realk)  :: omen(k%ncol),KF(k%ncol)
    real(realk)  :: dia(k%ncol)
    real(realk)  :: dia2(k%ncol)
    real(realk)  :: dia3(k%ncol)
    type(matrix) :: temp,tempT
    logical :: logx,logy,logz
    integer :: m,mat1,mat2,mat3,i
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz
    integer :: xi,xi1,xi2,xixi1,xixi2,xixixj1,&
         &xixixj2,xixj1,xixj2,xixi,xixixj,xixj,xj,xj1,xj2

    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24

    if (logx) then
       mat1=x;mat2=xx;mat3=xxx
    elseif (logy) then
       mat1=y;mat2=yy;mat3=yyy
    elseif (logz) then
       mat1=z;mat2=zz;mat3=zzz
    end if

    call mat_init(temp,K%nrow,K%ncol)
    call mat_init(tempT,K%nrow,K%ncol)
    call mat_mul(K,KURT%carmomMO(mat1),'n','n',1.0_realk,0.0_realk,temp)
    !*** update KF
    if (m>1) then
       call mat_extract_diagonal(dia,KURT%carmomMO(mat3))
       dia = -8.0_realk*dia
       call mat_extract_diagonal(dia2,KURT%carmomMO(mat2))
       call mat_extract_diagonal(dia3,KURT%carmomMO(mat1))
       dia = dia+24.0_realk*dia2*dia3-24.0_realk*dia3*dia3*dia3
       call mat_extract_diagonal(dia2,temp)
       KF = KF + dia*dia2 

       if (KURT%crossterms) then
          if (logx) then
             xj=x;xi1=y;xi2=z;xixi1=yy;xixi2=zz;xixixj1=xyy;xixixj2=xzz;xixj1=xy;xixj2=xz
          elseif(logy) then
             xj=y;xi1=x;xi2=z;xixi1=xx;xixi2=zz;xixixj1=xxy;xixixj2=yzz;xixj1=xy;xixj2=yz
          elseif(logz) then
             xj=z;xi1=x;xi2=y;xixi1=xx;xixi2=yy;xixixj1=xxz;xixixj2=yyz;xixj1=xz;xixj2=yz
          end if

          do i=1,2
             if (i==1) then
                xi=xi1;xixi=xixi1;xixixj=xixixj1;xixj=xixj1
             else 
                xi=xi2;xixi=xixi2;xixixj=xixixj2;xixj=xixj2
             end if

             call mat_extract_diagonal(dia,KURT%carmomMO(xixi))
             call mat_extract_diagonal(dia2,KURT%carmomMO(xj))
             dia = 2.0_realk*dia*dia2
             call mat_extract_diagonal(dia2,KURT%carmomMO(xixixj))
             dia = dia-2.0_realk*dia2
             call mat_extract_diagonal(dia3,KURT%carmomMO(xj))
             call mat_extract_diagonal(dia2,KURT%carmomMO(xi))
             dia = dia-6.0_realk*dia2*dia2*dia3 
             call mat_extract_diagonal(dia3,KURT%carmomMO(xixj))
             dia = dia+4.0_realk*dia3*dia2
             call mat_extract_diagonal(dia3,temp)
             KF = KF +4.0_realk*dia3*dia
          end do
       end if
    end if
    !*** end update KF
    call mat_trans(temp,tempT) 
    call mat_daxpy(1.0_realk,tempT,temp)

    call mat_extract_diagonal(dia2,temp)
    dia = dia2*omen
    call mat_dmul(dia,KURT%carmomMO(mat3),'T',dble(8*m),1.0_realk,HK)

    call mat_extract_diagonal(dia,KURT%carmomMO(mat2))
    call mat_extract_diagonal(dia2,temp)
    dia = dia*dia2*omen
    call mat_dmul(dia,KURT%carmomMO(mat1),'T',-dble(24*m),1.0_realk,HK)

    call mat_extract_diagonal(dia,KURT%carmomMO(mat1))
    call mat_extract_diagonal(dia2,temp)
    dia = dia*dia*dia2*omen
    call mat_dmul(dia,KURT%carmomMO(mat1),'T',dble(72*m),1.0_realk,HK)

    call mat_extract_diagonal(dia,KURT%carmomMO(mat3))
    dia=8.0_realk*dia
    call mat_extract_diagonal(dia2,KURT%carmomMO(mat1))
    call mat_extract_diagonal(dia3,KURT%carmomMO(mat2))
    dia = dia-24.0_realk*dia2*dia3
    dia= dia +24.0_realk*dia2*dia2*dia2
    dia = dia*omen

    call mat_dmul(dia,temp,'T',dble(m),1.0_realk,HK)

    call mat_extract_diagonal(dia,temp)
    call mat_extract_diagonal(dia2,KURT%carmomMO(mat1))
    dia = dia*dia2*omen
    call mat_dmul(dia,KURT%carmomMO(mat2),'T',-dble(24*m),1.0_realk,HK)



    if (KURT%crossterms) then
       if  (logx) then
          xj=x;xi1=y;xi2=z;xixi1=yy;xixi2=zz;xixixj1=xyy;xixixj2=xzz;xixj1=xy;xixj2=xz
       elseif (logy) then 
          xj=y;xi1=x;xi2=z;xixi1=xx;xixi2=zz;xixixj1=xxy;xixixj2=yzz;xixj1=xy;xixj2=yz
       elseif (logz) then
          xj=z;xi1=x;xi2=y;xixi1=xx;xixi2=yy;xixixj1=xxz;xixixj2=yyz;xixj1=xz;xixj2=yz
       endif

       do i=1,2
          if (i==1) then
             xi=xi1;xixi=xixi1;xixj=xixj1;xixixj=xixixj1
          elseif (i==2) then
             xi=xi2;xixi=xixi2;xixj=xixj2;xixixj=xixixj2 
          endif
          dia=0.0;dia2=0.0;dia3=0.0

          call mat_extract_diagonal(dia,KURT%carmomMO(xixi))
          call mat_extract_diagonal(dia2,KURT%carmomMO(xj))
          dia = -4.0_realk*dia*dia2
          call mat_extract_diagonal(dia3,KURT%carmomMO(xi))
          dia= dia+12.0_realk*dia3*dia3*dia2
          call mat_extract_diagonal(dia2,KURT%carmomMO(xixj)) 
          dia=dia-8.0_realk*dia2*dia3
          call mat_extract_diagonal(dia2,KURT%carmomMO(xixixj))
          dia = dia+4.0_realk*dia2
          dia = dia*omen
          call mat_dmul(dia,temp,'T',dble(2*m),1.0_realk,HK)

          call mat_extract_diagonal(dia2,temp)
          call mat_extract_diagonal(dia,KURT%carmomMO(xj))
          dia = dia*dia2*omen
          call mat_dmul(dia,KURT%carmomMO(xixi),'T',-dble(8*m),1.0_realk,HK)

          call mat_extract_diagonal(dia,KURT%carmomMO(xi))
          dia = dia*dia2*omen
          call mat_dmul(dia,KURT%carmomMO(xixj),'T',-dble(16*m),1.0_realk,HK)

          dia = dia2*omen
          call mat_dmul(dia,KURT%carmomMO(xixixj),'T',dble(8*m),1.0_realk,HK)

          call mat_extract_diagonal(dia,KURT%carmomMO(xixi))
          call mat_extract_diagonal(dia2,temp)
          dia = dia*dia2*omen
          call mat_dmul(dia,KURT%carmomMO(xj),'T',-dble(8*m),1.0_realk,HK)

          call mat_extract_diagonal(dia,KURT%carmomMO(xi))
          call mat_extract_diagonal(dia3,KURT%carmomMO(xj))
          call mat_extract_diagonal(dia2,temp)
          dia = dia*dia2*dia3*omen
          call mat_dmul(dia,KURT%carmomMO(xi),'T',dble(48*m),1.0_realk,HK)

          call mat_extract_diagonal(dia,KURT%carmomMO(xixj))
          dia=dia*dia2*omen
          call mat_dmul(dia,KURT%carmomMO(xi),'T',-dble(16*m),1.0_realk,HK)

          call mat_extract_diagonal(dia,KURT%carmomMO(xi))
          call mat_extract_diagonal(dia2,temp)
          dia = dia*dia*dia2*omen
          call mat_dmul(dia,KURT%carmomMO(xj),'T',dble(24*m),1.0_realk,HK)


       end do


    end if

    call mat_free(temp);call mat_free(tempt)
  end subroutine KX_lintrans

  subroutine KXX_lintrans(k,KURT,omen,m,HK,logx,logy,logz,KF)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: K,HK
    type(matrix) :: temp,tempT,temp1
    real(realk)  :: omen(K%ncol),KF(K%ncol)
    real(realk)  :: dia(K%ncol)
    real(realk)  :: dia2(K%ncol)
    logical :: logx,logy,logz
    integer :: m,mat1,mat2,mat3,i
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz
    integer :: xi,xi1,xi2,xjxj,xixi,xixi1,xixi2
    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24
    if (logx) then
       mat1=x;mat2=xx;mat3=xxx
    elseif(logy) then
       mat1=y;mat2=yy;mat3=yyy
    elseif(logz) then
       mat1=z;mat2=zz;mat3=zzz
    end if

    call mat_init(temp,K%nrow,K%ncol)
    call mat_init(tempT,K%nrow,K%ncol)
    call mat_mul(K,KURT%carmomMO(mat2),'n','n',1.0_realk,0.0_realk,temp)
    ! *** update KF
    if (m>1) then
       call mat_extract_diagonal(dia,temp)
       call mat_extract_diagonal(dia2,KURT%carmomMO(mat1))
       KF = KF+12.0_realk*dia*dia2*dia2

       if (KURT%crossterms) then
          if (logx) then
             xjxj=xx;xi1=y;xi2=z
          elseif (logy) then
             xjxj=yy;xi1=x;xi2=z
          elseif (logz) then
             xjxj=zz;xi1=x;xi2=y
          end if

          do i=1,2
             if (i==1) xi=xi1
             if (i==2) xi=xi2
             call mat_extract_diagonal(dia,temp)
             call mat_extract_diagonal(dia2,KURT%carmomMO(xi))
             KF = KF+4.0_realk*dia*dia2*dia2 
          end do
       end if
    end if
    ! *** end update KF
    call mat_trans(temp,tempT)
    call mat_daxpy(1.0_realk,tempT,temp)
    call mat_free(tempT)

    call mat_extract_diagonal(dia,KURT%carmomMO(mat1))
    dia = dia*dia
    dia = dia*omen
    call mat_dmul(dia,temp,'T',-dble(12*m),1.0_realk,HK)

    call mat_extract_diagonal(dia,temp)
    call mat_extract_diagonal(dia2,KURT%carmomMO(mat1))
    dia = dia*dia2*omen
    call mat_dmul(dia,KURT%carmomMO(mat1),'T',-dble(24*m),1.0_realk,HK)


    if (KURT%crossterms) then
       if  (logx) then
          xi1=2;xi2=3;xixi1=yy;xixi2=zz
       elseif (logy) then 
          xi1=1;xi2=3;xixi1=xx;xixi2=zz
       elseif (logz) then
          xi1=1;xi2=2;xixi1=xx;xixi2=yy
       endif

       do i=1,2
          if (i==1) then
             xi=xi1;xixi=xixi1
          elseif (i==2) then
             xi=xi2;xixi=xixi2
          end if
          call mat_extract_diagonal(dia,KURT%carmomMO(xi))
          call mat_extract_diagonal(dia2,temp)
          dia = dia*dia2*omen
          call mat_dmul(dia,KURT%carmomMO(xi),'T',-dble(8*m),1.0_realk,HK)
          call mat_extract_diagonal(dia,KURT%carmomMO(xi))
          dia = dia*dia*omen
          call mat_dmul(dia,temp,'T',-dble(4*m),1.0_realk,HK)
       enddo
    end if


    call mat_free(temp)
  end subroutine KXX_lintrans


  subroutine KXXX_lintrans(K,xxx,x,omen,m,HK,KF)
    implicit none
    type(matrix) :: xxx,x,HK,K
    real(realk)  :: omen(x%ncol),KF(K%ncol)
    real(realk)  :: dia(x%ncol)
    real(realk)  :: dia2(x%ncol)
    integer      :: m
    type(matrix) :: temp,tempT

    call mat_init(temp,K%nrow,K%ncol)

    call mat_mul(K,xxx,'n','n',1.0_realk,0.0_realk,temp) !temp=KX^3
    ! *** update KF
    if (m>1) then
       call mat_extract_diagonal(dia,temp)
       call mat_extract_diagonal(dia2,x)
       KF =KF-8.0_realk*dia*dia2
    end if
    ! *** end update KF
    call mat_init(tempT,K%nrow,K%ncol)
    call mat_trans(temp,tempT) !tempT = (KX^3)^T=-X^3K
    call mat_daxpy(1.0_realk,tempT,temp) !temp=(KX^3-X^3K)
    call mat_free(tempT)
    !compute all terms containing KX^3
    call mat_extract_diagonal(dia,x)
    dia=dia*omen
    call mat_dmul(dia,temp,'T',dble(8*m),1.0_realk,HK)
    call mat_extract_diagonal(dia,temp)
    dia = dia*omen
    call mat_dmul(dia,x,'T',dble(8*m),1.0_realk,HK)
    call mat_free(temp)

  end subroutine KXXX_lintrans


  subroutine KXXXX_lintrans(K,XXXX,omen,m,HK,KF)
    implicit none
    type(matrix),intent(in) :: xxxx,k
    type(matrix) :: HK
    real(realk)  :: KF(k%ncol),omen(k%ncol)
    real(realk)  :: dia(k%ncol)
    integer :: m
    type(matrix) :: temp,tempt

    call mat_init(temp,HK%nrow,HK%ncol)
    call mat_init(tempt,HK%nrow,HK%ncol)

    call mat_mul(k,xxxx,'n','n',1.0_realk,0.0_realk,temp) !temp=KX^4
    !**** update KF to contain KX^4 term*******
    if (m>1) then
       call mat_extract_diagonal(dia,temp)
       KF=KF+2.0_realk*dia
    end if
    !**** end update KF to contain KX^4 term***
    call mat_trans(temp,tempT) !tempT = (KX^4)^T=-X^4K
    call mat_daxpy(1.0_realk,tempT,temp) !temp=(KX^4-X^4K)
    call mat_dmul(omen,temp,'T',-dble(2*m),1.0_realk,HK)

    call mat_free(temp)
    call mat_free(tempt)

  end subroutine KXXXX_lintrans

!> \brief compute matrix containg diagonal hessian elements
!> \param P diagonal elements (output)
  subroutine kurtosis_precond_matrix(KURT,P)
    implicit none
    type(PFMitem),intent(in) :: KURT
    type(matrix) :: P
    real(realk),pointer ::  dia(:,:) 
    type(matrix):: temp
    real(realk) :: TEMPVEC(KURT%norb)
    real(realk),pointer :: TEMPVEC2(:)
    real(realk) :: omen(KURT%norb)
    type(matrix) :: KF,Pt,tmp
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz
    integer :: i,m,norb

    call mem_alloc(dia,KURT%norb,24)

    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24

    call mat_zero(P)

    do i=1,24
       !   call mat_clone(R(i),KURT%carmomMO(i))
       call mat_extract_diagonal(dia(:,i),KURT%carmomMO(i))
    end do
    m=KURT%m
    norb=KURT%norb
    if (m>2) then
       omen=KURT%omega**(m-1)
    elseif(m==2) then
       omen=KURT%omega
    elseif(m==1) then
       omen=1.0_realk
    end if

    call mat_dger(dble(2*m),omen,dia(:,xxxx),P)
    call mat_dger(dble(2*m),omen,dia(:,yyyy),P)
    call mat_dger(dble(2*m),omen,dia(:,zzzz),P)
    TEMPVEC=1.0_realk
    call mat_dger(-dble(2*m),omen*dia(:,xxxx),TEMPVEC,P)
    call mat_dger(-dble(2*m),omen*dia(:,yyyy),TEMPVEC,P)
    call mat_dger(-dble(2*m),omen*dia(:,zzzz),TEMPVEC,P)

    call mat_dhmul(omen,KURT%carmomMO(xxx),KURT%carmomMO(x),'n','n',-dble(32*m),1.0_realk,P)
    call mat_dhmul(omen,KURT%carmomMO(yyy),KURT%carmomMO(y),'n','n',-dble(32*m),1.0_realk,P)
    call mat_dhmul(omen,KURT%carmomMO(zzz),KURT%carmomMO(z),'n','n',-dble(32*m),1.0_realk,P)

    call mat_dger(-dble(8*m),omen*dia(:,xxx),dia(:,x),P)
    call mat_dger(-dble(8*m),omen*dia(:,yyy),dia(:,y),P)
    call mat_dger(-dble(8*m),omen*dia(:,zzz),dia(:,z),P)
    TEMPVEC=1.0_realk
    call mat_dger(dble(16*m),omen*dia(:,xxx)*dia(:,x),TEMPVEC,P)
    call mat_dger(dble(16*m),omen*dia(:,yyy)*dia(:,y),TEMPVEC,P)
    call mat_dger(dble(16*m),omen*dia(:,zzz)*dia(:,z),TEMPVEC,P)


    call mat_dger(-dble(8*m),omen*dia(:,x),dia(:,xxx),P)
    call mat_dger(-dble(8*m),omen*dia(:,y),dia(:,yyy),P)
    call mat_dger(-dble(8*m),omen*dia(:,z),dia(:,zzz),P)
    TEMPVEC=1.0_realk

    call mat_dger(dble(24*m),omen*dia(:,xx)*dia(:,x),dia(:,x),P)
    call mat_dger(dble(24*m),omen*dia(:,yy)*dia(:,y),dia(:,y),P)
    call mat_dger(dble(24*m),omen*dia(:,zz)*dia(:,z),dia(:,z),P)
    call mat_dger(-dble(24*m),omen*dia(:,xx)*dia(:,x)*dia(:,x),TEMPVEC,P)
    call mat_dger(-dble(24*m),omen*dia(:,yy)*dia(:,y)*dia(:,y),TEMPVEC,P)
    call mat_dger(-dble(24*m),omen*dia(:,zz)*dia(:,z)*dia(:,z),TEMPVEC,P)

    call mat_dger(dble(12*m),omen*dia(:,x)*dia(:,x),dia(:,xx),P)
    call mat_dger(dble(12*m),omen*dia(:,y)*dia(:,y),dia(:,yy),P)
    call mat_dger(dble(12*m),omen*dia(:,z)*dia(:,z),dia(:,zz),P)
    call mat_dger(-dble(12*m),omen*dia(:,x)*dia(:,x)*dia(:,xx),TEMPVEC,P)
    call mat_dger(-dble(12*m),omen*dia(:,y)*dia(:,y)*dia(:,yy),TEMPVEC,P)
    call mat_dger(-dble(12*m),omen*dia(:,z)*dia(:,z)*dia(:,zz),TEMPVEC,P)

    !call mat_clone(temp,KURT%carmomMO(x))
    call mat_dhmul(omen*dia(:,xx),KURT%carmomMO(x),KURT%carmomMO(x),'n','n',dble(48*m),1.0_realk,P)
    !call mat_clone(temp,KURT%carmomMO(y))
    call mat_dhmul(omen*dia(:,yy),KURT%carmomMO(y),KURT%carmomMO(y),'n','n',dble(48*m),1.0_realk,P)
    !call mat_clone(temp,KURT%carmomMO(z))
    call mat_dhmul(omen*dia(:,zz),KURT%carmomMO(z),KURT%carmomMO(z),'n','n',dble(48*m),1.0_realk,P)

    call mat_dhmul(omen*dia(:,x),KURT%carmomMO(x),KURT%carmomMO(xx),'n','n',dble(96*m),1.0_realk,P)
    call mat_dhmul(omen*dia(:,y),KURT%carmomMO(y),KURT%carmomMO(yy),'n','n',dble(96*m),1.0_realk,P)
    call mat_dhmul(omen*dia(:,z),KURT%carmomMO(z),KURT%carmomMO(zz),'n','n',dble(96*m),1.0_realk,P)

    call mat_dger(-dble(24*m),omen*dia(:,x)*dia(:,x)*dia(:,x),dia(:,x),P)
    call mat_dger(-dble(24*m),omen*dia(:,y)*dia(:,y)*dia(:,y),dia(:,y),P)
    call mat_dger(-dble(24*m),omen*dia(:,z)*dia(:,z)*dia(:,z),dia(:,z),P)

    call mat_dger(dble(24*m),omen*dia(:,x)*dia(:,x)*dia(:,x)*dia(:,x),TEMPVEC,P)
    call mat_dger(dble(24*m),omen*dia(:,y)*dia(:,y)*dia(:,y)*dia(:,y),TEMPVEC,P)
    call mat_dger(dble(24*m),omen*dia(:,z)*dia(:,z)*dia(:,z)*dia(:,z),TEMPVEC,P)


    !call mat_clone(temp,KURT%carmomMO(x))
    call mat_dhmul(omen*dia(:,x)*dia(:,x),KURT%carmomMO(X),KURT%carmomMO(x),'n','n',-dble(144*m),1.0_realk,P)
    !call mat_clone(temp,KURT%carmomMO(y))
    call mat_dhmul(omen*dia(:,y)*dia(:,y),KURT%carmomMO(y),KURT%carmomMO(y),'n','n',-dble(144*m),1.0_realk,P)
    !call mat_clone(temp,KURT%carmomMO(z))
    call mat_dhmul(omen*dia(:,z)*dia(:,z),KURT%carmomMO(z),KURT%carmomMO(z),'n','n',-dble(144*m),1.0_realk,P)

    if (KURT%crossterms) then
       call mat_dhmul(omen*dia(:,xx),KURT%carmomMO(y),KURT%carmomMO(y),'n','n',dble(4*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,yy),KURT%carmomMO(x),KURT%carmomMO(x),'n','n',dble(4*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,xx),KURT%carmomMO(z),KURT%carmomMO(z),'n','n',dble(4*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,zz),KURT%carmomMO(x),KURT%carmomMO(x),'n','n',dble(4*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,yy),KURT%carmomMO(z),KURT%carmomMO(z),'n','n',dble(4*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,zz),KURT%carmomMO(y),KURT%carmomMO(y),'n','n',dble(4*4*m),1.0_realk,P) 


       call mat_dhmul(omen*dia(:,x),KURT%carmomMO(x),KURT%carmomMO(yy),'n','n',dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,y),KURT%carmomMO(y),KURT%carmomMO(xx),'n','n',dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,x),KURT%carmomMO(x),KURT%carmomMO(zz),'n','n',dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z),KURT%carmomMO(z),KURT%carmomMO(xx),'n','n',dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,y),KURT%carmomMO(y),KURT%carmomMO(zz),'n','n',dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z),KURT%carmomMO(z),KURT%carmomMO(yy),'n','n',dble(2*16*m),1.0_realk,P) 

       call mat_dhmul(omen,KURT%carmomMO(x),KURT%carmomMO(xyy),'n','n',-dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen,KURT%carmomMO(y),KURT%carmomMO(xxy),'n','n',-dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen,KURT%carmomMO(x),KURT%carmomMO(xzz),'n','n',-dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen,KURT%carmomMO(z),KURT%carmomMO(xxz),'n','n',-dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen,KURT%carmomMO(y),KURT%carmomMO(yzz),'n','n',-dble(2*16*m),1.0_realk,P) 
       call mat_dhmul(omen,KURT%carmomMO(z),KURT%carmomMO(yyz),'n','n',-dble(2*16*m),1.0_realk,P) 

       call mat_dhmul(omen*dia(:,x),KURT%carmomMO(y),KURT%carmomMO(xy),'n','n',dble(2*32*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,y),KURT%carmomMO(x),KURT%carmomMO(xy),'n','n',dble(2*32*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,x),KURT%carmomMO(z),KURT%carmomMO(xz),'n','n',dble(2*32*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z),KURT%carmomMO(x),KURT%carmomMO(xz),'n','n',dble(2*32*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,y),KURT%carmomMO(z),KURT%carmomMO(yz),'n','n',dble(2*32*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z),KURT%carmomMO(y),KURT%carmomMO(yz),'n','n',dble(2*32*m),1.0_realk,P) 

       call mat_dhmul(omen*dia(:,xy),KURT%carmomMO(x),KURT%carmomMO(y),'n','n',dble(64*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,xz),KURT%carmomMO(x),KURT%carmomMO(z),'n','n',dble(64*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,yz),KURT%carmomMO(y),KURT%carmomMO(z),'n','n',dble(64*m),1.0_realk,P) 


       call mat_dhmul(omen*dia(:,x)*dia(:,y),KURT%carmomMO(x),KURT%carmomMO(y),'n','n',-dble(24*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,y)*dia(:,x),KURT%carmomMO(y),KURT%carmomMO(x),'n','n',-dble(24*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,x)*dia(:,z),KURT%carmomMO(x),KURT%carmomMO(z),'n','n',-dble(24*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z)*dia(:,x),KURT%carmomMO(z),KURT%carmomMO(x),'n','n',-dble(24*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,y)*dia(:,z),KURT%carmomMO(y),KURT%carmomMO(z),'n','n',-dble(24*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z)*dia(:,y),KURT%carmomMO(z),KURT%carmomMO(y),'n','n',-dble(24*4*m),1.0_realk,P) 


       call mat_dhmul(omen*dia(:,y)*dia(:,y),KURT%carmomMO(x),KURT%carmomMO(x),'n','n',-dble(12*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,x)*dia(:,x),KURT%carmomMO(y),KURT%carmomMO(y),'n','n',-dble(12*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,x)*dia(:,x),KURT%carmomMO(z),KURT%carmomMO(z),'n','n',-dble(12*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z)*dia(:,z),KURT%carmomMO(x),KURT%carmomMO(x),'n','n',-dble(12*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,z)*dia(:,z),KURT%carmomMO(y),KURT%carmomMO(y),'n','n',-dble(12*4*m),1.0_realk,P) 
       call mat_dhmul(omen*dia(:,y)*dia(:,y),KURT%carmomMO(z),KURT%carmomMO(z),'n','n',-dble(12*4*m),1.0_realk,P) 


       call mem_alloc(tempvec2,norb)
       TEMPVEC2= 1.0_realk
       call mat_dger(dble(4*m),omen,dia(:,xxyy),P) 
       call mat_dger(-dble(4*m),omen*dia(:,xxyy),tempvec2,P) 
       call mat_dger(dble(4*m),omen,dia(:,xxzz),P) 
       call mat_dger(-dble(4*m),omen*dia(:,xxzz),tempvec2,P) 
       call mat_dger(dble(4*m),omen,dia(:,yyzz),P) 
       call mat_dger(-dble(4*m),omen*dia(:,yyzz),tempvec2,P) 

       call mat_dger(dble(16*m),omen*dia(:,x)*dia(:,y),dia(:,xy),P)
       call mat_dger(-dble(16*m),omen*dia(:,x)*dia(:,y)*dia(:,xy),tempvec2,P)
       call mat_dger(dble(16*m),omen*dia(:,x)*dia(:,z),dia(:,xz),P)
       call mat_dger(-dble(16*m),omen*dia(:,x)*dia(:,z)*dia(:,xz),tempvec2,P)
       call mat_dger(dble(16*m),omen*dia(:,y)*dia(:,z),dia(:,yz),P)
       call mat_dger(-dble(16*m),omen*dia(:,y)*dia(:,z)*dia(:,yz),tempvec2,P)


       TEMPVEC = -8.0_realk*dia(:,xx)*dia(:,y)+8_realk*dia(:,xxy)+24.0_realk*dia(:,x)*dia(:,x)*dia(:,y)&
            &-16.0_realk*dia(:,xy)*dia(:,x) 
       call mat_dger(-dble(m),omen*TEMPVEC,dia(:,y),P)
       call mat_dger(dble(m),omen*TEMPVEC*dia(:,y),tempvec2,P)
       TEMPVEC = -8.0_realk*dia(:,yy)*dia(:,x)+8_realk*dia(:,xyy)+24.0_realk*dia(:,y)*dia(:,y)*dia(:,x)&
            &-16.0_realk*dia(:,xy)*dia(:,y) 
       call mat_dger(-dble(m),omen*TEMPVEC,dia(:,x),P)
       call mat_dger(dble(m),omen*TEMPVEC*dia(:,x),tempvec2,P)

       TEMPVEC = -8.0_realk*dia(:,xx)*dia(:,z)+8_realk*dia(:,xxz)+24.0_realk*dia(:,x)*dia(:,x)*dia(:,z)&
            &-16.0_realk*dia(:,xz)*dia(:,x) 
       call mat_dger(-dble(m),omen*TEMPVEC,dia(:,z),P)
       call mat_dger(dble(m),omen*TEMPVEC*dia(:,z),tempvec2,P)
       TEMPVEC = -8.0_realk*dia(:,zz)*dia(:,x)+8_realk*dia(:,xzz)+24.0_realk*dia(:,z)*dia(:,z)*dia(:,x)&
            &-16.0_realk*dia(:,xz)*dia(:,z) 
       call mat_dger(-dble(m),omen*TEMPVEC,dia(:,x),P)
       call mat_dger(dble(m),omen*TEMPVEC*dia(:,x),tempvec2,P)

       TEMPVEC = -8.0_realk*dia(:,yy)*dia(:,z)+8_realk*dia(:,yyz)+24.0_realk*dia(:,y)*dia(:,y)*dia(:,z)&
            &-16.0_realk*dia(:,yz)*dia(:,y) 
       call mat_dger(-dble(m),omen*TEMPVEC,dia(:,z),P)
       call mat_dger(dble(m),omen*TEMPVEC*dia(:,z),tempvec2,P)
       TEMPVEC = -8.0_realk*dia(:,zz)*dia(:,y)+8_realk*dia(:,yzz)+24.0_realk*dia(:,z)*dia(:,z)*dia(:,y)&
            &-16.0_realk*dia(:,yz)*dia(:,z) 
       call mat_dger(-dble(m),omen*TEMPVEC,dia(:,y),P)
       call mat_dger(dble(m),omen*TEMPVEC*dia(:,y),tempvec2,P)

       call mat_dger(dble(4*m),omen*dia(:,x)*dia(:,x),dia(:,yy),P)
       call mat_dger(-dble(4*m),omen*dia(:,x)*dia(:,x)*dia(:,yy),tempvec2,P)
       call mat_dger(dble(4*m),omen*dia(:,y)*dia(:,y),dia(:,xx),P)
       call mat_dger(-dble(4*m),omen*dia(:,y)*dia(:,y)*dia(:,xx),tempvec2,P)

       call mat_dger(dble(4*m),omen*dia(:,x)*dia(:,x),dia(:,zz),P)
       call mat_dger(-dble(4*m),omen*dia(:,x)*dia(:,x)*dia(:,zz),tempvec2,P)
       call mat_dger(dble(4*m),omen*dia(:,z)*dia(:,z),dia(:,xx),P)
       call mat_dger(-dble(4*m),omen*dia(:,z)*dia(:,z)*dia(:,xx),tempvec2,P)

       call mat_dger(dble(4*m),omen*dia(:,y)*dia(:,y),dia(:,zz),P)
       call mat_dger(-dble(4*m),omen*dia(:,y)*dia(:,y)*dia(:,zz),tempvec2,P)
       call mat_dger(dble(4*m),omen*dia(:,z)*dia(:,z),dia(:,yy),P)
       call mat_dger(-dble(4*m),omen*dia(:,z)*dia(:,z)*dia(:,yy),tempvec2,P)

       call mat_dger(-dble(8*m),omen*dia(:,x),dia(:,xyy),P)
       call mat_dger(dble(8*m),omen*dia(:,x)*dia(:,xyy),tempvec2,P)
       call mat_dger(-dble(8*m),omen*dia(:,y),dia(:,xxy),P)
       call mat_dger(dble(8*m),omen*dia(:,y)*dia(:,xxy),tempvec2,P)

       call mat_dger(-dble(8*m),omen*dia(:,x),dia(:,xzz),P)
       call mat_dger(dble(8*m),omen*dia(:,x)*dia(:,xzz),tempvec2,P)
       call mat_dger(-dble(8*m),omen*dia(:,z),dia(:,xxz),P)
       call mat_dger(dble(8*m),omen*dia(:,z)*dia(:,xxz),tempvec2,P)

       call mat_dger(-dble(8*m),omen*dia(:,y),dia(:,yzz),P)
       call mat_dger(dble(8*m),omen*dia(:,y)*dia(:,yzz),tempvec2,P)
       call mat_dger(-dble(8*m),omen*dia(:,z),dia(:,yyz),P)
       call mat_dger(dble(8*m),omen*dia(:,z)*dia(:,yyz),tempvec2,P)

       call mem_dealloc(tempvec2)

    end if

    if (m>1) then
       call mat_init(Pt,P%nrow,P%ncol)
       call mat_zero(Pt)
       call mat_daxpy(dble(-2*m),KURT%carmomMO(xxxx),Pt)
       call mat_daxpy(dble(-2*m),KURT%carmomMO(yyyy),Pt)
       call mat_daxpy(dble(-2*m),KURT%carmomMO(zzzz),Pt)

       call mat_dmul(dia(:,x),KURT%carmomMO(xxx),'n',dble(8*m),1.0_realk,Pt)
       call mat_dmul(dia(:,y),KURT%carmomMO(yyy),'n',dble(8*m),1.0_realk,Pt)
       call mat_dmul(dia(:,z),KURT%carmomMO(zzz),'n',dble(8*m),1.0_realk,Pt)

       call mat_dmul(dia(:,xxx),KURT%carmomMO(x),'n',dble(8*m),1.0_realk,Pt)
       call mat_dmul(dia(:,yyy),KURT%carmomMO(y),'n',dble(8*m),1.0_realk,Pt)
       call mat_dmul(dia(:,zzz),KURT%carmomMO(z),'n',dble(8*m),1.0_realk,Pt)

       call mat_dmul(dia(:,x)*dia(:,x),KURT%carmomMO(xx),'n',-dble(12*m),1.0_realk,Pt)
       call mat_dmul(dia(:,y)*dia(:,y),KURT%carmomMO(yy),'n',-dble(12*m),1.0_realk,Pt)
       call mat_dmul(dia(:,z)*dia(:,z),KURT%carmomMO(zz),'n',-dble(12*m),1.0_realk,Pt)

       call mat_dmul(dia(:,xx)*dia(:,x),KURT%carmomMO(x),'n',-dble(24*m),1.0_realk,Pt)
       call mat_dmul(dia(:,yy)*dia(:,y),KURT%carmomMO(y),'n',-dble(24*m),1.0_realk,Pt)
       call mat_dmul(dia(:,zz)*dia(:,z),KURT%carmomMO(z),'n',-dble(24*m),1.0_realk,Pt)

       call mat_dmul(dia(:,x)*dia(:,x)*dia(:,x),KURT%carmomMO(x),'n',dble(24*m),1.0_realk,Pt)
       call mat_dmul(dia(:,y)*dia(:,y)*dia(:,y),KURT%carmomMO(y),'n',dble(24*m),1.0_realk,Pt)
       call mat_dmul(dia(:,z)*dia(:,z)*dia(:,z),KURT%carmomMO(z),'n',dble(24*m),1.0_realk,Pt)

       if (KURT%crossterms) then
          call mat_daxpy(-dble(4*m),KURT%carmomMO(xxyy),Pt)  
          call mat_daxpy(-dble(4*m),KURT%carmomMO(xxzz),Pt)  
          call mat_daxpy(-dble(4*m),KURT%carmomMO(yyzz),Pt)  

          call mat_dmul(dia(:,x)*dia(:,y),KURT%carmomMO(xy),'n',-dble(16*m),1.0_realk,Pt)
          call mat_dmul(dia(:,x)*dia(:,z),KURT%carmomMO(xz),'n',-dble(16*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z)*dia(:,y),KURT%carmomMO(yz),'n',-dble(16*m),1.0_realk,Pt)

          call mat_dmul(dia(:,xx)*dia(:,y),KURT%carmomMO(y),'n',-dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,yy)*dia(:,x),KURT%carmomMO(x),'n',-dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,xx)*dia(:,z),KURT%carmomMO(z),'n',-dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,zz)*dia(:,x),KURT%carmomMO(x),'n',-dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,yy)*dia(:,z),KURT%carmomMO(z),'n',-dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,zz)*dia(:,y),KURT%carmomMO(y),'n',-dble(8*m),1.0_realk,Pt)

          call mat_dmul(dia(:,x)*dia(:,x)*dia(:,y),KURT%carmomMO(y),'n',dble(24*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y)*dia(:,y)*dia(:,x),KURT%carmomMO(x),'n',dble(24*m),1.0_realk,Pt)
          call mat_dmul(dia(:,x)*dia(:,x)*dia(:,z),KURT%carmomMO(z),'n',dble(24*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z)*dia(:,z)*dia(:,x),KURT%carmomMO(x),'n',dble(24*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y)*dia(:,y)*dia(:,z),KURT%carmomMO(z),'n',dble(24*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z)*dia(:,z)*dia(:,y),KURT%carmomMO(y),'n',dble(24*m),1.0_realk,Pt)

          call mat_dmul(dia(:,x)*dia(:,xy),KURT%carmomMO(y),'n',-dble(16*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y)*dia(:,xy),KURT%carmomMO(x),'n',-dble(16*m),1.0_realk,Pt)
          call mat_dmul(dia(:,x)*dia(:,xz),KURT%carmomMO(z),'n',-dble(16*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z)*dia(:,xz),KURT%carmomMO(x),'n',-dble(16*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y)*dia(:,yz),KURT%carmomMO(z),'n',-dble(16*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z)*dia(:,yz),KURT%carmomMO(y),'n',-dble(16*m),1.0_realk,Pt)

          call mat_dmul(dia(:,xxy),KURT%carmomMO(y),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,xyy),KURT%carmomMO(x),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,xxz),KURT%carmomMO(z),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,xzz),KURT%carmomMO(x),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,yyz),KURT%carmomMO(z),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,yzz),KURT%carmomMO(y),'n',dble(8*m),1.0_realk,Pt)

          call mat_dmul(dia(:,x)*dia(:,x),KURT%carmomMO(yy),'n',-dble(4*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y)*dia(:,y),KURT%carmomMO(xx),'n',-dble(4*m),1.0_realk,Pt)
          call mat_dmul(dia(:,x)*dia(:,x),KURT%carmomMO(zz),'n',-dble(4*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z)*dia(:,z),KURT%carmomMO(xx),'n',-dble(4*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y)*dia(:,y),KURT%carmomMO(zz),'n',-dble(4*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z)*dia(:,z),KURT%carmomMO(yy),'n',-dble(4*m),1.0_realk,Pt)

          call mat_dmul(dia(:,x),KURT%carmomMO(xyy),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y),KURT%carmomMO(xxy),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,x),KURT%carmomMO(xzz),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z),KURT%carmomMO(xxz),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,y),KURT%carmomMO(yzz),'n',dble(8*m),1.0_realk,Pt)
          call mat_dmul(dia(:,z),KURT%carmomMO(yyz),'n',dble(8*m),1.0_realk,Pt)
       end if
       call mat_init(KF,P%nrow,P%ncol)
       call dia_KF_kl2(KURT,KF,dia)
       call mat_init(tmp,Pt%ncol,Pt%nrow)
       call mat_trans(Pt,tmp)
       call mat_free(Pt)
       call mat_hmul(1.0_realk,KF,tmp,1.0_realk,P)
       call mat_free(tmp)
       call mat_free(KF)
    end if

    call mat_init(temp,P%ncol,P%nrow)
    call mat_trans(P,temp)
    call mat_daxpy(1.0_realk,temp,P)
    call mat_scal_dia(0.0_realk,P)
    call mat_free(temp)
    call mem_dealloc(dia)

  end subroutine kurtosis_precond_matrix



  subroutine kurtosis_precond(Xout,X,mu,KURT)
    implicit none
    type(PFMitem),intent(in) :: KURT
    type(matrix) :: Xout
    type(matrix) :: X,Xtemp
    real(realk),intent(in) :: mu
    integer :: i,ne
    real(realk),pointer :: tmp(:),tmpP(:),xvec(:),yvec(:)

    select case(matrix_type)
    case(mtype_dense)
       ne = Xout%nrow*Xout%ncol
       call mat_zero(Xout)
       do i=1,ne
          if (dabs(KURT%p%elms(i)- mu) > 1d-8 ) then 
             Xout%elms(i)=X%elms(i)/(KURT%P%elms(i)-mu)
          else
             Xout%elms(i)=X%elms(i)
          end if
       end do
    case(mtype_scalapack)
       call mat_copy(1.0_realk,X,Xout)
       call mat_hdiv(Xout,KURT%P,mu)
    case default
       ne = Xout%nrow*Xout%ncol
       call mem_alloc(tmp,X%nrow*X%ncol)
       call mem_alloc(tmpP,KURT%P%nrow*KURT%P%ncol)
       call mat_to_full(X,1E0_realk,tmp)
       call mat_to_full(KURT%P,1E0_realk,tmpP)

       do i=1,ne
          if (dabs(tmpP(i) - mu)> 1d-8) tmp(i) = tmp(i)/(tmpP(i) - mu)
       enddo

       call mat_set_from_full(tmp,1E0_realk,Xout)
       call mem_dealloc(tmp)
       call mem_dealloc(tmpP) 
    end select

  end subroutine kurtosis_precond

  subroutine dia_KF_kl2(KURT,KF,dia)
    implicit none
    type(PFMitem) :: KURT
    integer :: m 
    type(matrix) :: KF!, R(24)
    real(realk)  :: dia(KURT%norb,24)
    type(matrix) :: tmp,tmp1
    real(realk),pointer :: tmpvec(:)
    integer :: i,j
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz

    m=KURT%m
    call mat_zero(KF)
    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24


    !do i=1,24
    !   call mat_clone(R(i),KURT%carmomMO(i))
    !end do

    call mat_daxpy(2.0_realk,KURT%carmomMO(xxxx),KF)
    call mat_daxpy(2.0_realk,KURT%carmomMO(yyyy),KF)
    call mat_daxpy(2.0_realk,KURT%carmomMO(zzzz),KF)

    call mat_dmul(dia(:,x),KURT%carmomMO(xxx),'n',-8.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,y),KURT%carmomMO(yyy),'n',-8.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,z),KURT%carmomMO(zzz),'n',-8.0_realk,1.0_realk,KF)

    call mat_dmul(dia(:,xxx),KURT%carmomMO(x),'n',-8.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,yyy),KURT%carmomMO(y),'n',-8.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,zzz),KURT%carmomMO(z),'n',-8.0_realk,1.0_realk,KF)

    call mat_dmul(dia(:,x)*dia(:,x),KURT%carmomMO(xx),'n',12.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,y)*dia(:,y),KURT%carmomMO(yy),'n',12.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,z)*dia(:,z),KURT%carmomMO(zz),'n',12.0_realk,1.0_realk,KF)

    call mat_dmul(dia(:,x)*dia(:,xx),KURT%carmomMO(x),'n',24.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,y)*dia(:,yy),KURT%carmomMO(y),'n',24.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,z)*dia(:,zz),KURT%carmomMO(z),'n',24.0_realk,1.0_realk,KF)

    call mat_dmul(dia(:,x)*dia(:,x)*dia(:,x),KURT%carmomMO(x),'n',-24.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,y)*dia(:,y)*dia(:,y),KURT%carmomMO(y),'n',-24.0_realk,1.0_realk,KF)
    call mat_dmul(dia(:,z)*dia(:,z)*dia(:,z),KURT%carmomMO(z),'n',-24.0_realk,1.0_realk,KF)

    if (KURT%crossterms) then
       call mat_init(tmp,KURT%norb,KURT%norb)
       call mat_zero(tmp)
       call mat_daxpy(4.0_realk,KURT%carmomMO(xxyy),tmp)  
       call mat_daxpy(4.0_realk,KURT%carmomMO(xxzz),tmp)  
       call mat_daxpy(4.0_realk,KURT%carmomMO(yyzz),tmp)  

       call mat_dmul(dia(:,x)*dia(:,y),KURT%carmomMO(xy),'n',16.0_realk,1.0_realk,tmp)
       call mat_dmul(dia(:,x)*dia(:,z),KURT%carmomMO(xz),'n',16.0_realk,1.0_realk,tmp)
       call mat_dmul(dia(:,y)*dia(:,z),KURT%carmomMO(yz),'n',16.0_realk,1.0_realk,tmp)

       call mem_alloc(tmpvec,KURT%norb)
       tmpvec =  dia(:,xx)*dia(:,y)-dia(:,xxy)-3.0_realk*dia(:,y)*dia(:,x)*dia(:,x)+2.0_realk*dia(:,x)*dia(:,xy) 
       call mat_dmul(tmpvec,KURT%carmomMO(y),'n',8.0_realk,1.0_realk,tmp)
       tmpvec =  dia(:,yy)*dia(:,x)-dia(:,xyy)-3.0_realk*dia(:,x)*dia(:,y)*dia(:,y)+2.0_realk*dia(:,y)*dia(:,xy) 
       call mat_dmul(tmpvec,KURT%carmomMO(x),'n',8.0_realk,1.0_realk,tmp)

       tmpvec =  dia(:,xx)*dia(:,z)-dia(:,xxz)-3.0_realk*dia(:,z)*dia(:,x)*dia(:,x)+2.0_realk*dia(:,x)*dia(:,xz) 
       call mat_dmul(tmpvec,KURT%carmomMO(z),'n',8.0_realk,1.0_realk,tmp)
       tmpvec =  dia(:,zz)*dia(:,x)-dia(:,xzz)-3.0_realk*dia(:,x)*dia(:,z)*dia(:,z)+2.0_realk*dia(:,z)*dia(:,xz) 
       call mat_dmul(tmpvec,KURT%carmomMO(x),'n',8.0_realk,1.0_realk,tmp)

       tmpvec =  dia(:,yy)*dia(:,z)-dia(:,yyz)-3.0_realk*dia(:,z)*dia(:,y)*dia(:,y)+2.0_realk*dia(:,y)*dia(:,yz) 
       call mat_dmul(tmpvec,KURT%carmomMO(z),'n',8.0_realk,1.0_realk,tmp)
       tmpvec =  dia(:,zz)*dia(:,y)-dia(:,yzz)-3.0_realk*dia(:,y)*dia(:,z)*dia(:,z)+2.0_realk*dia(:,z)*dia(:,yz) 
       call mat_dmul(tmpvec,KURT%carmomMO(y),'n',8.0_realk,1.0_realk,tmp)
       call mem_dealloc(tmpvec)

       call mat_dmul(dia(:,y),KURT%carmomMO(xxy),'n',-8.0_realk,1.0_realk,tmp)
       call mat_dmul(dia(:,x),KURT%carmomMO(xyy),'n',-8.0_realk,1.0_realk,tmp)
       call mat_dmul(dia(:,z),KURT%carmomMO(xxz),'n',-8.0_realk,1.0_realk,tmp)
       call mat_dmul(dia(:,x),KURT%carmomMO(xzz),'n',-8.0_realk,1.0_realk,tmp)
       call mat_dmul(dia(:,z),KURT%carmomMO(yyz),'n',-8.0_realk,1.0_realk,tmp)
       call mat_dmul(dia(:,y),KURT%carmomMO(yzz),'n',-8.0_realk,1.0_realk,tmp)


       call mat_dmul(dia(:,x)*dia(:,x),KURT%carmomMO(yy),'n',4.0_realk,1.0_realk,tmp) 
       call mat_dmul(dia(:,y)*dia(:,y),KURT%carmomMO(xx),'n',4.0_realk,1.0_realk,tmp) 

       call mat_dmul(dia(:,x)*dia(:,x),KURT%carmomMO(zz),'n',4.0_realk,1.0_realk,tmp) 
       call mat_dmul(dia(:,z)*dia(:,z),KURT%carmomMO(xx),'n',4.0_realk,1.0_realk,tmp) 

       call mat_dmul(dia(:,y)*dia(:,y),KURT%carmomMO(zz),'n',4.0_realk,1.0_realk,tmp) 
       call mat_dmul(dia(:,z)*dia(:,z),KURT%carmomMO(yy),'n',4.0_realk,1.0_realk,tmp) 

       call mat_daxpy(1.0_realk,tmp,KF)
       call mat_free(tmp)
    end if

    if (m>2) then
       call mem_alloc(tmpvec,KURT%norb)
       tmpvec=1.0_realk
       call mat_init(tmp,KF%nrow,KF%ncol)
       call mat_init(tmp1,KF%nrow,KF%ncol)
       call mat_zero(tmp)
       call mat_dger(1.0_realk,KURT%omega**(m-2),tmpvec,tmp)
       call mat_hmul(1.0_realk,tmp,KF,0.0_realk,tmp1)
       call mat_copy(1.0_realk,tmp1,KF)
       call mat_free(tmp)
       call mat_free(tmp1)
       call mem_dealloc(tmpvec)
    end if

    call mat_init(tmp,KF%nrow,KF%ncol)
    call mat_trans(KF,tmp)
    call mat_copy(-dble(m-1),tmp,KF)
    call mat_free(tmp)


  end subroutine dia_KF_kl2

! Make (f^[1])_kl
  subroutine dia_KF_kl(KURT,KF)
    implicit none
    type(PFMitem) :: KURT
    integer :: m 
    real(realk),pointer :: KF(:,:) 
    real(realk),pointer :: temp(:,:)
    real(realk),pointer :: temp2(:,:)
    real(realk) :: dia(KURT%norb),dia2(KURT%norb)
    integer :: i,j
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz
    type(matrix) :: tmp

    call mem_alloc(KF,KURT%norb,KURT%norb)
    call mem_alloc(temp,KURT%norb,KURT%norb)
    call mem_alloc(temp2,KURT%norb,KURT%norb)

    m=KURT%m
    KF = 0_realk
    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24

    !***********************************************
    call mat_to_full(KURT%carmomMO(xxxx),1.0_realk,temp)
    KF = 2.0_realk*temp
    call mat_to_full(KURT%carmomMO(yyyy),1.0_realk,temp)
    KF = KF+2.0_realk*temp
    call mat_to_full(KURT%carmomMO(zzzz),1.0_realk,temp)
    KF = KF+2.0_realk*temp

    !***********************************************
    call mat_extract_diagonal(dia,KURT%carmomMO(x))
    do i=1,KURT%norb
       temp2(i,:)=dia
    end do
    call mat_to_full(KURT%carmomMO(xxx),1.0_realk,temp)
    KF = KF -8.0_realk*temp*temp2
    call mat_extract_diagonal(dia,KURT%carmomMO(y))
    do i=1,KURT%norb
       temp2(i,:)=dia
    end do
    call mat_to_full(KURT%carmomMO(yyy),1.0_realk,temp)
    KF = KF -8.0_realk*temp*temp2
    call mat_extract_diagonal(dia,KURT%carmomMO(z))
    do i=1,KURT%norb
       temp2(i,:)=dia
    end do
    call mat_to_full(KURT%carmomMO(zzz),1.0_realk,temp)
    KF = KF -8.0_realk*temp*temp2

    !***********************************************
    call mat_extract_diagonal(dia,KURT%carmomMO(xxx))
    do i=1,KURT%norb
       temp2(i,:)=dia
    end do
    call mat_to_full(KURT%carmomMO(x),1.0_realk,temp)
    KF = KF- 8.0_realk*temp2*temp
    call mat_extract_diagonal(dia,KURT%carmomMO(yyy))
    do i=1,KURT%norb
       temp2(i,:)=dia
    end do
    call mat_to_full(KURT%carmomMO(y),1.0_realk,temp)
    KF = KF- 8.0_realk*temp2*temp
    call mat_extract_diagonal(dia,KURT%carmomMO(zzz))
    do i=1,KURT%norb
       temp2(i,:)=dia
    end do
    call mat_to_full(KURT%carmomMO(z),1.0_realk,temp)
    KF = KF- 8.0_realk*temp2*temp
    !***********************************************
    call mat_extract_diagonal(dia,KURT%carmomMO(x))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)= dia(j)*dia(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(xx),1.0_realk,temp)
    KF = KF +12.0_realk*temp2*temp
    call mat_extract_diagonal(dia,KURT%carmomMO(y))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)= dia(j)*dia(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(yy),1.0_realk,temp)
    KF = KF +12.0_realk*temp2*temp
    call mat_extract_diagonal(dia,KURT%carmomMO(z))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)= dia(j)*dia(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(zz),1.0_realk,temp)
    KF = KF +12.0_realk*temp2*temp
    !***********************************************
    call mat_extract_diagonal(dia,KURT%carmomMO(x))
    call mat_extract_diagonal(dia2,KURT%carmomMO(xx))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)=dia(j)*dia2(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(x),1.0_realk,temp)
    KF = KF + 24.0_realk*temp2*temp
    call mat_extract_diagonal(dia,KURT%carmomMO(y))
    call mat_extract_diagonal(dia2,KURT%carmomMO(yy))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)=dia(j)*dia2(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(y),1.0_realk,temp)
    KF = KF + 24.0_realk*temp2*temp
    call mat_extract_diagonal(dia,KURT%carmomMO(z))
    call mat_extract_diagonal(dia2,KURT%carmomMO(zz))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)=dia(j)*dia2(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(z),1.0_realk,temp)
    KF = KF + 24.0_realk*temp2*temp
    !***********************************************
    call mat_extract_diagonal(dia,KURT%carmomMO(x))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)=dia(j)*dia(j)*dia(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(x),1.0_realk,temp)
    KF = KF - 24.0_realk*temp*temp2
    call mat_extract_diagonal(dia,KURT%carmomMO(y))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)=dia(j)*dia(j)*dia(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(y),1.0_realk,temp)
    KF = KF - 24.0_realk*temp*temp2
    call mat_extract_diagonal(dia,KURT%carmomMO(z))
    do i=1,KURT%norb
       do j=1,KURT%norb
          temp2(i,j)=dia(j)*dia(j)*dia(j)
       end do
    end do
    call mat_to_full(KURT%carmomMO(z),1.0_realk,temp)
    KF = KF - 24.0_realk*temp*temp2
    !***********************************************

    !**** CROSSTERMS *****
    CROSSTERMS: if (KURT%crossterms) then
       call dia_crossKF_xxyy(KURT%carmomMO(xxyy),KURT%carmomMO(xxzz),KURT%carmomMO(yyzz),KF)

       call dia_crossKF_xy(KURT%carmomMO(xy),KURT%carmomMO(x),KURT%carmomMO(y),KF)
       call dia_crossKF_xy(KURT%carmomMO(xz),KURT%carmomMO(x),KURT%carmomMO(z),KF)
       call dia_crossKF_xy(KURT%carmomMO(yz),KURT%carmomMO(y),KURT%carmomMO(z),KF)

       call dia_crossKF_xx(KURT%carmomMO(xx),KURT%carmomMO(y),KF) 
       call dia_crossKF_xx(KURT%carmomMO(yy),KURT%carmomMO(x),KF) 

       call dia_crossKF_xx(KURT%carmomMO(xx),KURT%carmomMO(z),KF) 
       call dia_crossKF_xx(KURT%carmomMO(zz),KURT%carmomMO(x),KF) 

       call dia_crossKF_xx(KURT%carmomMO(yy),KURT%carmomMO(z),KF) 
       call dia_crossKF_xx(KURT%carmomMO(zz),KURT%carmomMO(y),KF) 

       call dia_crossKF_x(KURT%carmomMO(x),KURT%carmomMO(y),KF)
       call dia_crossKF_x(KURT%carmomMO(y),KURT%carmomMO(x),KF)

       call dia_crossKF_x(KURT%carmomMO(x),KURT%carmomMO(z),KF)
       call dia_crossKF_x(KURT%carmomMO(z),KURT%carmomMO(x),KF)

       call dia_crossKF_x(KURT%carmomMO(y),KURT%carmomMO(z),KF)
       call dia_crossKF_x(KURT%carmomMO(z),KURT%carmomMO(y),KF)

       call dia_crossKF_xixj(KURT%carmomMO(xy),KURT%carmomMO(x),KURT%carmomMO(y),KF)
       call dia_crossKF_xixj(KURT%carmomMO(xy),KURT%carmomMO(y),KURT%carmomMO(x),KF)

       call dia_crossKF_xixj(KURT%carmomMO(xz),KURT%carmomMO(x),KURT%carmomMO(z),KF)
       call dia_crossKF_xixj(KURT%carmomMO(xz),KURT%carmomMO(z),KURT%carmomMO(x),KF)

       call dia_crossKF_xixj(KURT%carmomMO(yz),KURT%carmomMO(y),KURT%carmomMO(z),KF)
       call dia_crossKF_xixj(KURT%carmomMO(yz),KURT%carmomMO(z),KURT%carmomMO(y),KF)

       call dia_crossKF_xxy(KURT%carmomMO(xxy),KURT%carmomMO(y),KF)
       call dia_crossKF_xxy(KURT%carmomMO(xyy),KURT%carmomMO(x),KF)

       call dia_crossKF_xxy(KURT%carmomMO(xxz),KURT%carmomMO(z),KF)
       call dia_crossKF_xxy(KURT%carmomMO(xzz),KURT%carmomMO(x),KF)

       call dia_crossKF_xxy(KURT%carmomMO(yyz),KURT%carmomMO(z),KF)
       call dia_crossKF_xxy(KURT%carmomMO(yzz),KURT%carmomMO(y),KF)

    end if  CROSSTERMS



    if (m>2) then
       do i=1,KURT%norb
          do j=1,KURT%norb
             temp(i,j)=KURT%omega(j)**(m-2)
          end do
       end do
       KF = KF*temp
    end if

    KF = -dble(m-1)*KF
    call mat_init(tmp,KURT%norb,KURT%norb)
    call mat_set_from_full(KF,1.0_realk,tmp)

    call mat_free(tmp)
    call mem_dealloc(KF)
    call mem_dealloc(temp)
    call mem_dealloc(temp2)

  end subroutine dia_KF_kl

  subroutine dia_crossKF_xxy(xxy,y,KF)
    implicit none
    type(matrix) :: xxy,y
    real(realk)  :: KF(y%ncol,y%ncol)
    real(realk)  :: temp(y%ncol,y%ncol)
    real(realk)  :: dia1(y%ncol),dia2(y%ncol)
    integer :: i,j

    call mat_to_full(y,1.0_realk,temp)
    call mat_extract_diagonal(dia1,xxy)

    do i=1,y%ncol
       do j=1,y%ncol
          KF(i,j) = KF(i,j) -8.0_realk*temp(i,j)*dia1(j) 
       end do
    end do

    call mat_to_full(xxy,1.0_realk,temp)
    call mat_extract_diagonal(dia1,y)

    do i=1,y%ncol
       do j=1,y%ncol
          KF(i,j) = KF(i,j) - 8.0_realk*temp(i,j)*dia1(j) 
       end do
    end do


  end subroutine dia_crossKF_xxy



  subroutine dia_crossKF_xixj(xy,x,y,KF)
    implicit none
    type(matrix) :: xy,x,y
    real(realk),pointer  :: KF(:,:)
    real(realk),pointer  :: temp(:,:)
    real(realk)  :: dia1(x%ncol),dia2(x%ncol)
    integer :: i,j

    call mem_alloc(KF,x%ncol,x%ncol)
    call mem_alloc(temp,x%ncol,x%ncol)

    call mat_to_full(x,1.0_realk,temp)
    call mat_extract_diagonal(dia1,xy)
    call mat_extract_diagonal(dia2,y)

    do i=1,x%ncol
       do j=1,x%ncol
          KF(i,j) = KF(i,j)+16.0_realk*temp(i,j)*dia1(j)*dia2(j)
       end do
    end do

    call mem_dealloc(KF)
    call mem_dealloc(temp)
  end subroutine dia_crossKF_xixj


  subroutine dia_crossKF_x(x,y,KF)
    implicit none
    type(matrix) :: x,y
    real(realk),pointer  :: KF(:,:)
    real(realk),pointer  :: temp(:,:)
    real(realk)  :: dia1(x%ncol),dia2(x%ncol)
    integer :: i,j

    call mem_alloc(KF,x%ncol,x%ncol)
    call mem_alloc(temp,x%ncol,x%ncol)


    call mat_to_full(x,1.0_realk,temp)
    call mat_extract_diagonal(dia1,x)
    call mat_extract_diagonal(dia2,y)

    do i=1,x%ncol
       do j=1,x%ncol
          KF(i,j) = KF(i,j) -24.0_realk*temp(i,j)*dia1(j)*dia2(j)*dia2(j)
       end do
    end do

    call mem_dealloc(KF)
    call mem_dealloc(temp)
  end subroutine dia_crossKF_x



  subroutine dia_crossKF_xx(xx,y,KF)
    implicit none
    type(matrix) :: xx,y
    real(realk),pointer  :: KF(:,:)
    real(realk),pointer  :: temp(:,:)
    real(realk)  :: dia1(y%ncol),dia2(y%ncol)
    integer :: i,j

    call mem_alloc(KF,y%ncol,y%ncol)
    call mem_alloc(temp,y%ncol,y%ncol)
    call mat_to_full(y,1.0_realk,temp)
    call mat_extract_diagonal(dia1,xx)
    call mat_extract_diagonal(dia2,y)

    do i=1,y%ncol
       do j=1,y%ncol
          KF(i,j)=KF(i,j)+8.0_realk*temp(i,j)*dia1(j)*dia2(j)
       end do
    end do

    call mat_to_full(xx,1.0_realk,temp)
    call mat_extract_diagonal(dia1,y)

    do i=1,y%ncol
       do j=1,y%ncol
          KF(i,j)=KF(i,j)+4.0_realk*temp(i,j)*dia1(j)*dia1(j)
       end do
    end do

    call mem_dealloc(KF)
    call mem_dealloc(temp)

  end subroutine dia_crossKF_xx


  subroutine dia_crossKF_xy(xy,x,y,KF)
    implicit none
    type(matrix) :: xy,x,y
    real(realk),pointer  :: KF(:,:)
    real(realk),pointer  :: temp(:,:)
    real(realk)  :: dia1(x%ncol),dia2(x%ncol)
    integer :: i,j

    call mem_alloc(KF,y%ncol,y%ncol)
    call mem_alloc(temp,y%ncol,y%ncol)
    call mat_to_full(xy,1.0_realk,temp)
    call mat_extract_diagonal(dia1,x)
    call mat_extract_diagonal(dia2,y)

    do i=1,x%ncol
       do j=1,x%ncol
          KF(i,j) = KF(i,j)+16.0_realk*temp(i,j)*dia1(j)*dia2(j)
       end do
    end do
    call mem_dealloc(KF)
    call mem_dealloc(temp)

  end subroutine dia_crossKF_xy

  subroutine dia_crossKF_xxyy(xxyy,xxzz,yyzz,KF)
    implicit none
    type(matrix) :: xxyy,xxzz,yyzz
    real(realk),pointer  :: KF(:,:)
    real(realk),pointer  :: temp(:,:)

    call mem_alloc(KF,xxyy%ncol,xxyy%ncol)
    call mem_alloc(temp,xxyy%ncol,xxyy%ncol)

    call mat_to_full(xxyy,4.0_realk,temp)
    KF = KF + temp
    call mat_to_full(xxzz,4.0_realk,temp)
    KF = KF + temp
    call mat_to_full(yyzz,4.0_realk,temp)
    KF = KF + temp

    call mem_dealloc(KF)
    call mem_dealloc(temp)

  end subroutine dia_crossKF_xxyy

! ****************************************
! *  ONLY FOR FINITE DIFFERENCE TESTING  * 
! ****************************************

  subroutine FiniteDiff_gradient(KURT,cmo,Gpq)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: carmomMOp(24),carmomMOm(24)
    type(matrix) :: expXm,expXp
    type(matrix) :: x
    type(matrix) :: cmo
    real(realk)  :: eps
    real(realk)  :: omegap(cmo%ncol),omegam(cmo%ncol)
    real(realk)  :: zetap,zetam,Gpq
    integer :: p,q,i

    do i=1,24
       call mat_init(carmomMOp(i),KURT%carmomMO(i)%nrow,KURT%carmomMO(i)%ncol)
       call mat_init(carmomMOm(i),KURT%carmomMO(i)%nrow,KURT%carmomMO(i)%ncol)
    end do
    call mat_init(expXp,cmo%ncol,cmo%ncol)
    call mat_init(expXm,cmo%ncol,cmo%ncol)
    call mat_init(X,cmo%ncol,cmo%ncol)

    eps=5.0d-5
    p=2;q=4
    call mat_zero(x)
    call expX_eps(expXp,x,eps,p,q)  
    call mat_zero(x)
    call expX_eps(expXm,x,-eps,p,q)  

    !call mat_print(expXp,1,4,1,4,6)
    do i=1,24
       call carmomAO_to_carmomMO(KURT%carmomMO(i),expXp,carmomMOp(i)) 
       call carmomAO_to_carmomMO(KURT%carmomMO(i),expXm,carmomMOm(i)) 
    end do
    call compute_omega(carmomMOp,omegap,KURT%crossterms,cmo%ncol)
    call compute_omega(carmomMOm,omegam,KURT%crossterms,cmo%ncol)

    do i=1,cmo%ncol
       omegap(i)=omegap(i)**KURT%m
       omegam(i)=omegam(i)**KURT%m
    end do
    zetap = sum(omegap)
    zetam = sum(omegam)

    Gpq = (zetap-zetam )/(2.0_realk*eps)

    print*, "Gpq  = ",Gpq, "for p,q =",p,q

    do i=1,24
       call mat_free(carmomMOp(i))
       call mat_free(carmomMOm(i))
    end do
    call mat_free(expXp)
    call mat_free(expXm)
    call mat_free(X)
  end subroutine FiniteDiff_gradient


  subroutine FiniteDiff_HK(KURT,cmo,HKpq)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: carmomMOpp(24),carmomMOmm(24)
    type(matrix) :: carmomMOpm(24),carmomMOmp(24)
    type(matrix) :: expXmm,expXpp
    type(matrix) :: expXmp,expXpm
    type(matrix) :: x
    type(matrix) :: cmo
    real(realk)  :: eps
    real(realk)  :: omegapp(cmo%ncol),omegamm(cmo%ncol)
    real(realk)  :: omegapm(cmo%ncol),omegamp(cmo%ncol)
    real(realk)  :: zetapp,zetamm,Gpq
    real(realk)  :: zetapm,zetamp,HKpq
    type(matrix) :: temp
    integer :: p,q,r,s,i

    do i=1,24
       call mat_init(carmomMOpp(i),KURT%carmomMO(i)%nrow,KURT%carmomMO(i)%ncol)
       call mat_init(carmomMOmm(i),KURT%carmomMO(i)%nrow,KURT%carmomMO(i)%ncol)
       call mat_init(carmomMOpm(i),KURT%carmomMO(i)%nrow,KURT%carmomMO(i)%ncol)
       call mat_init(carmomMOmp(i),KURT%carmomMO(i)%nrow,KURT%carmomMO(i)%ncol)
    end do
    call mat_init(expXpp,cmo%ncol,cmo%ncol)
    call mat_init(expXmm,cmo%ncol,cmo%ncol)
    call mat_init(expXpm,cmo%ncol,cmo%ncol)
    call mat_init(expXmp,cmo%ncol,cmo%ncol)
    call mat_init(temp,cmo%ncol,cmo%ncol)
    call mat_init(X,cmo%ncol,cmo%ncol)

    eps=1d-4
    p=2;q=4
    r=2;s=4
    call mat_zero(x)
    call expX_eps(expXpp,x,eps,p,q)  
    call expX_eps(expXpp,x,eps,r,s)  

    call mat_zero(x)
    call expX_eps(expXpm,x,eps,p,q)  
    call expX_eps(expXpm,x,-eps,r,s)  

    call mat_zero(x)
    call expX_eps(expXmm,x,-eps,p,q)  
    call expX_eps(expXmm,x,-eps,r,s)  

    call mat_zero(x)
    call expX_eps(expXmp,x,-eps,p,q)  
    call expX_eps(expXmp,x,eps,r,s)  
    !call mat_print(expXp,1,4,1,4,6)
    do i=1,24
       call carmomAO_to_carmomMO(KURT%carmomMO(i),expXpp,carmomMOpp(i)) 
       call carmomAO_to_carmomMO(KURT%carmomMO(i),expXmm,carmomMOmm(i)) 
       call carmomAO_to_carmomMO(KURT%carmomMO(i),expXpm,carmomMOpm(i)) 
       call carmomAO_to_carmomMO(KURT%carmomMO(i),expXmp,carmomMOmp(i)) 
    end do
    call compute_omega(carmomMOpp,omegapp,KURT%crossterms,cmo%ncol)
    call compute_omega(carmomMOmm,omegamm,KURT%crossterms,cmo%ncol)
    call compute_omega(carmomMOpm,omegapm,KURT%crossterms,cmo%ncol)
    call compute_omega(carmomMOmp,omegamp,KURT%crossterms,cmo%ncol)

    do i=1,cmo%ncol
       omegapp(i)=omegapp(i)**KURT%m
       omegamm(i)=omegamm(i)**KURT%m
       omegapm(i)=omegapm(i)**KURT%m
       omegamp(i)=omegamp(i)**KURT%m
    end do


    zetapp = sum(omegapp)
    zetamm = sum(omegamm)
    zetapm = sum(omegapm)
    zetamp = sum(omegamp)

    HKpq = ((zetapp-zetapm)-(zetamp-zetamm))/(4.0_realk*eps*eps)

    print*, "HKpq  = ",HKpq, "for p,q =",p,q, " and r,s = ",r,s 

    do i=1,24
       call mat_free(carmomMOpp(i))
       call mat_free(carmomMOpm(i))
       call mat_free(carmomMOmm(i))
       call mat_free(carmomMOmp(i))
    end do
    call mat_free(expXpp)
    call mat_free(expXmm)
    call mat_free(expXpm)
    call mat_free(expXmp)
    call mat_free(temp)
    call mat_free(X)
  end subroutine FiniteDiff_HK

  subroutine expX_eps(expX,x,eps,p,q) 
    implicit none
    type(matrix) :: X,expX
    real(realk)  :: xmat(X%nrow,X%ncol)
    real(realk)  :: eps
    integer :: p,q

    call mat_to_full(X,1.0_realk,xmat)

    xmat(p,q)=xmat(p,q)+eps
    xmat(q,p)=xmat(q,p)-eps

    call mat_set_from_full(xmat,1.0_realk,x)
    call matrix_exponential(X,expX,1.0d-12)

  end subroutine expX_eps

  subroutine debug_lintra(KURT,kappa,G,HK)
    implicit none
    type(PFMitem) :: KURT
    type(matrix) :: kappa
    type(matrix) :: HK,G
    type(matrix) :: temp,tempT
    real(realk)  :: dia(HK%ncol),omega(HK%ncol)
    integer :: m,i
    integer ::x,y,z,xx,xy,xz,yy,yz,zz,xxx,xxy,xxz,&
         &xyy,xzz,yyy,yyz,yzz,zzz,xxxx,xxyy,xxzz,yyyy,&
         &yyzz,zzzz

    m=KURT%m

    ! mapping to matrix numbers
    x=1;y=2;z=3;xx=4;xy=5;xz=6;yy=7;yz=8;zz=9
    xxx=10;xxy=11;xxz=12;xyy=13;xzz=14;yyy=15
    yyz=16;yzz=17;zzz=18;xxxx=19;xxyy=20;xxzz=21
    yyyy=22;yyzz=23;zzzz=24

    call mat_zero(HK)
    call mat_init(temp,kappa%nrow,kappa%ncol)
    call KXXX_lintrans_d(kappa,KURT%carmomMO(xxx),KURT%carmomMO(x),m,HK)

    call KXX_lintrans_d(kappa,KURT%carmomMO(xx),KURT%carmomMO(x),m,HK)

    call KX_lintrans_d(kappa,KURT%carmomMO(xxx),KURT%carmomMO(xx),KURT%carmomMO(x),m,HK)

    call mat_mul(G,kappa,'T','T',0.5_realk,1.0_realk,HK)
    !call mat_zero(temp)
    call mat_trans(HK,temp)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,temp,HK)
    print*,"G"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)
    call mat_mul(kappa,G,'n','n',-0.5_realk,0.0_realk,temp)

    !call mat_zero(temp)
    !call mat_trans(HK,temp)
    !call mat_scal(-1.0_realk,HK)
    !call mat_daxpy(1.0_realk,temp,HK)


    call mat_free(temp)

  end subroutine debug_lintra

  subroutine KX_lintrans_d(K,xxx,xx,x,m,HK)
    implicit none
    type(matrix) :: K,xxx,xx,x,HK
    real(realk)  :: dia(x%ncol)
    real(realk)  :: dia2(x%ncol)
    real(realk)  :: dia3(x%ncol)
    type(matrix) :: temp,tempT
    integer :: m

    call mat_init(temp,K%nrow,K%ncol)
    call mat_init(tempT,K%nrow,K%ncol)
    call mat_zero(temp)
    call mat_zero(tempt)
    call mat_mul(K,x,'n','n',1.0_realk,0.0_realk,temp)
    call mat_trans(temp,tempT) 
    call mat_daxpy(1.0_realk,tempT,temp)

    call mat_extract_diagonal(dia2,temp)
    call mat_dmul(dia2,xxx,'T',dble(8),1.0_realk,HK)
    !!debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KX del 1"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

    call mat_extract_diagonal(dia,x)
    call mat_extract_diagonal(dia2,temp)
    dia = dia*dia*dia2
    call mat_dmul(dia,x,'T',dble(72),1.0_realk,HK)
    !debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KX del 2"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

    call mat_extract_diagonal(dia,xxx)
    call mat_dmul(dia,temp,'T',dble(8),1.0_realk,HK)
    !!debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KX del 3a"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)
    call mat_extract_diagonal(dia2,x)
    call mat_extract_diagonal(dia3,xx)
    dia = dia2*dia3
    call mat_dmul(dia,temp,'T',-dble(24),1.0_realk,HK)
    !!debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KX del 3b"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

    call mat_extract_diagonal(dia2,x)
    dia = dia2*dia2*dia2
    call mat_dmul(dia,temp,'T',dble(24),1.0_realk,HK)
    !debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KX del 4"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

    call mat_extract_diagonal(dia,temp)
    call mat_extract_diagonal(dia2,x)
    dia = dia*dia2
    call mat_dmul(dia,xx,'T',-dble(24),1.0_realk,HK)
    !debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KX del 5"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

    call mat_free(tempt)
    call mat_free(temp)
  end subroutine KX_lintrans_d

  subroutine KXX_lintrans_d(K,xx,x,m,HK)
    implicit none
    type(matrix) :: XX,X,K,HK
    type(matrix) :: temp,tempT
    real(realk)  :: dia(K%ncol)
    real(realk)  :: dia2(K%ncol)
    integer :: m


    call mat_init(temp,K%nrow,K%ncol)
    call mat_init(tempT,K%nrow,K%ncol)
    call mat_mul(K,XX,'n','n',1.0_realk,0.0_realk,temp)
    call mat_trans(temp,tempT)
    call mat_daxpy(1.0_realk,tempT,temp)

    call mat_extract_diagonal(dia,x)
    dia = dia*dia
    call mat_dmul(dia,temp,'n',-dble(12),1.0_realk,HK)
    !debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KxX del 1"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

    call mat_extract_diagonal(dia,temp)
    call mat_extract_diagonal(dia2,x)
    dia = dia*dia2
    call mat_dmul(dia,x,'n',-dble(24),1.0_realk,HK)
    !debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KxX del 2"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

  end subroutine KXX_lintrans_d


  subroutine KXXX_lintrans_d(K,xxx,x,m,HK)
    implicit none
    type(matrix) :: xxx,x,HK,K
    real(realk)  :: dia(x%ncol)
    real(realk)  :: dia2(x%ncol)
    integer      :: m
    type(matrix) :: temp,tempT

    call mat_init(temp,K%nrow,K%ncol)

    call mat_mul(K,xxx,'n','n',1.0_realk,0.0_realk,temp) !temp=KX^3
    call mat_init(tempT,K%nrow,K%ncol)
    call mat_trans(temp,tempT) !tempT = (KX^3)^T=-X^3K
    call mat_daxpy(1.0_realk,tempT,temp) !temp=(KX^3-X^3K)
    !compute all terms containing KX^3
    call mat_extract_diagonal(dia,x)
    call mat_dmul(dia,temp,'n',dble(8),1.0_realk,HK)
    !debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KXxx del 1"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)

    call mat_extract_diagonal(dia,temp)
    call mat_dmul(dia,x,'n',dble(8),1.0_realk,HK)
    !debug
    call mat_trans(HK,tempt)
    call mat_scal(-1.0_realk,HK)
    call mat_daxpy(1.0_realk,tempt,HK)
    print*,"KxxX del 2"
    call mat_print(HK,1,4,1,4,6)
    call mat_zero(HK)



    call mat_free(temp)
    call mat_free(tempt)

  end subroutine KXXX_lintrans_d



end module kurtosis
