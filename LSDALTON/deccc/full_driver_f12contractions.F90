module full_f12contractions
  use precision
contains
!> @file
!> Full molecular CC program
!> \author Thomas Kjaergaard

#ifdef MOD_UNRELEASED
  subroutine ccsdf12_Vjiij_coupling(Vjiij,Ciajb,Taibj,Viajb,Vijja,Viaji,Tai,nocc,nvirt)
  implicit none
  real(realk),intent(INOUT) :: Vjiij(nocc,nocc)
  real(realk),intent(IN)    :: Ciajb(nocc,nvirt,nocc,nvirt)
  real(realk),intent(IN)    :: Viajb(nocc,nvirt,nocc,nvirt)
  real(realk),intent(IN)    :: Taibj(nvirt,nocc,nvirt,nocc)
  real(realk),intent(IN)    :: Vijja(nocc,nocc,nvirt)
  real(realk),intent(IN)    :: Viaji(nocc,nvirt,nocc)
  real(realk),intent(IN)    :: Tai(nvirt,nocc)

  integer,intent(IN)        :: nocc,nvirt
  !
  integer :: i,j,a,b
  real(realk) :: tmp

  do j=1,nocc
   do i=1,nocc
    tmp = 0E0_realk
    do b=1,nvirt
     do a=1,nvirt
      tmp = tmp + (Ciajb(i,a,j,b)+Viajb(i,a,j,b))*Taibj(a,j,b,i)
     enddo
    enddo
    Vjiij(i,j) = Vjiij(i,j) + tmp
   enddo
  enddo
  do j=1,nocc
   do i=1,nocc
    tmp = 0E0_realk
    do a=1,nvirt
     tmp = tmp + Viaji(i,a,j)*Tai(a,j) + Vijja(i,j,a)*Tai(a,i)
    enddo
    Vjiij(i,j) = Vjiij(i,j) + tmp
   enddo
  enddo

  end subroutine ccsdf12_Vjiij_coupling

  subroutine mp2f12_Vijij_coupling(Vijij,Ciajb,Taibj,nocc,nvirt)
  implicit none
  real(realk),intent(INOUT) :: Vijij(nocc,nocc)
  real(realk),intent(IN)    :: Ciajb(nocc,nvirt,nocc,nvirt)
  real(realk),intent(IN)    :: Taibj(nvirt,nocc,nvirt,nocc)
  integer,intent(IN)        :: nocc,nvirt
  !
  integer :: i,j,a,b
  real(realk) :: tmp

  ! Setting Ciajb = 0
  !Ciajb = 0.0E0_realk

  do j=1,nocc
   do i=1,nocc
    tmp = 0E0_realk
    do b=1,nvirt
     do a=1,nvirt
      tmp = tmp + Ciajb(i,a,j,b)*Taibj(a,i,b,j)
     enddo
    enddo
    Vijij(i,j) = Vijij(i,j) + tmp
   enddo
  enddo

  end subroutine mp2f12_Vijij_coupling


  subroutine mp2f12_Vjiij_coupling(Vjiij,Ciajb,Taibj,nocc,nvirt)
  implicit none
  real(realk),intent(INOUT) :: Vjiij(nocc,nocc)
  real(realk),intent(IN)    :: Ciajb(nocc,nvirt,nocc,nvirt)
  real(realk),intent(IN)    :: Taibj(nvirt,nocc,nvirt,nocc)
  integer,intent(IN)        :: nocc,nvirt
  !
  integer :: i,j,a,b
  real(realk) :: tmp

  ! Setting Ciajb = 0
  !Ciajb = 0.0E0_realk

  do j=1,nocc
   do i=1,nocc
    tmp = 0E0_realk
    do b=1,nvirt
     do a=1,nvirt
      tmp = tmp + Ciajb(j,a,i,b)*Taibj(a,i,b,j)
     enddo
    enddo
    Vjiij(i,j) = Vjiij(i,j) + tmp
   enddo
  enddo

  end subroutine mp2f12_Vjiij_coupling

  subroutine ccsdf12_Vijij_coupling(Vijij,Ciajb,Taibj,Viajb,Viija,Viajj,Tai,nocc,nvirt)
  implicit none
  real(realk),intent(INOUT) :: Vijij(nocc,nocc)
  real(realk),intent(IN)    :: Ciajb(nocc,nvirt,nocc,nvirt)
  real(realk),intent(IN)    :: Viajb(nocc,nvirt,nocc,nvirt)
  real(realk),intent(IN)    :: Taibj(nvirt,nocc,nvirt,nocc)
  real(realk),intent(IN)    :: Viija(nocc,nocc,nvirt)
  real(realk),intent(IN)    :: Viajj(nocc,nvirt,nocc)
  real(realk),intent(IN)    :: Tai(nvirt,nocc)

  integer,intent(IN)        :: nocc,nvirt
  !
  integer :: i,j,a,b
  real(realk) :: tmp

  do j=1,nocc
   do i=1,nocc
    tmp = 0E0_realk
    do b=1,nvirt
     do a=1,nvirt
      tmp = tmp + ( Ciajb(i,a,j,b) + Viajb(i,a,j,b) )*Taibj(a,i,b,j)
     enddo
    enddo
    Vijij(i,j) = Vijij(i,j) + tmp
   enddo
  enddo
  do j=1,nocc
   do i=1,nocc
    tmp = 0E0_realk
    do a=1,nvirt
     tmp = tmp + Viija(i,j,a) * Tai(a,j) + Viajj(i,a,j) * Tai(a,i)
    enddo
    Vijij(i,j) = Vijij(i,j) + tmp
   enddo
  enddo

  end subroutine ccsdf12_Vijij_coupling

  !> \brief Vijij contribution for MP2-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine mp2f12_Vijij(Vijij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(OUT) :: Vijij(nocc,nocc)
    Real(realk),intent(in)  :: Ripjq(nocc,nbasis,nocc,nbasis)
    Real(realk),intent(in)  :: Gipjq(nocc,nbasis,nocc,nbasis)
    Real(realk),intent(in)  :: Fijkl(nocc,nocc,nocc,nocc)
    Real(realk),intent(in)  :: Rimjc(nocc,noccfull,nocc,ncabs)
    Real(realk),intent(in)  :: Gimjc(nocc,noccfull,nocc,ncabs)
    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    !Vijij = 0.0E0_realk
    do j=1,nocc
       do i=1,nocc
          Vijij(i,j) = Fijkl(i,i,j,j)
       enddo
    enddo
    do q=1,nbasis
       do j=1,nocc
          do p=1,nbasis
             do i=1,nocc
                Vijij(i,j) = Vijij(i,j) - Ripjq(i,p,j,q)*Gipjq(i,p,j,q)
             enddo
          enddo
       enddo
    enddo
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vijij(i,j) = Vijij(i,j) - Rimjc(i,p,j,q)*Gimjc(i,p,j,q)
             enddo
          enddo
       enddo
    enddo
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vijij(i,j) = Vijij(i,j) - Rimjc(j,p,i,q)*Gimjc(j,p,i,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vijij

  !> \brief Vijij term1 contribution for MP2-f12
  !> \author Yang M. Wang
  !> \date March 2013
  subroutine mp2f12_Vijij_term1(Vijij,Fijkl,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(inout) :: Vijij(nocc,nocc)
    Real(realk),intent(in)  :: Fijkl(nocc,nocc,nocc,nocc)
    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vijij = 0.0E0_realk
    do j=1,nocc
       do i=1,nocc
          Vijij(i,j) = Fijkl(i,i,j,j)
       enddo
    enddo
  end subroutine mp2f12_Vijij_term1

  subroutine mp2f12_Vijij_term2(Vijij,Ripjq,Gipjq,nocc,noccfull,nbasis,ncabs)
    
    implicit none
    Real(realk),intent(inout) :: Vijij(nocc,nocc)
    Real(realk),intent(in)  :: Ripjq(nocc,nbasis,nocc,nbasis)
    Real(realk),intent(in)  :: Gipjq(nocc,nbasis,nocc,nbasis)

    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vijij = 0.0E0_realk
    do q=1,nbasis
       do j=1,nocc
          do p=1,nbasis
             do i=1,nocc
                Vijij(i,j) = Vijij(i,j) - Ripjq(i,p,j,q)*Gipjq(i,p,j,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vijij_term2

  subroutine mp2f12_Vijij_term3(Vijij,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(inout) :: Vijij(nocc,nocc)
    Real(realk),intent(in)  :: Rimjc(nocc,noccfull,nocc,ncabs)
    Real(realk),intent(in)  :: Gimjc(nocc,noccfull,nocc,ncabs)

    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vijij = 0.0E0_realk
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vijij(i,j) = Vijij(i,j) - Rimjc(i,p,j,q)*Gimjc(i,p,j,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vijij_term3

  subroutine mp2f12_Vijij_term4(Vijij,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(inout) :: Vijij(nocc,nocc)
    Real(realk),intent(in)  :: Rimjc(nocc,noccfull,nocc,ncabs)
    Real(realk),intent(in)  :: Gimjc(nocc,noccfull,nocc,ncabs)

    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vijij = 0.0E0_realk
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vijij(i,j) = Vijij(i,j) - Rimjc(j,p,i,q)*Gimjc(j,p,i,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vijij_term4

  !> \brief Vjiij contribution for MP2-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine mp2f12_Vjiij(Vjiij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(out) :: Vjiij(nocc,nocc)
    Real(realk),intent(in)  :: Ripjq(nocc,nbasis,nocc,nbasis)
    Real(realk),intent(in)  :: Gipjq(nocc,nbasis,nocc,nbasis)
    Real(realk),intent(in)  :: Fijkl(nocc,nocc,nocc,nocc)
    Real(realk),intent(in)  :: Rimjc(nocc,noccfull,nocc,ncabs)
    Real(realk),intent(in)  :: Gimjc(nocc,noccfull,nocc,ncabs)
    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    !Vjiij = 0.0E0_realk
    do j=1,nocc
       do i=1,nocc
          Vjiij(i,j) = Fijkl(i,j,j,i)
       enddo
    enddo
    do q=1,nbasis
       do j=1,nocc
          do p=1,nbasis
             do i=1,nocc
                Vjiij(i,j) = Vjiij(i,j) - Ripjq(i,p,j,q)*Gipjq(j,p,i,q)
             enddo
          enddo
       enddo
    enddo
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vjiij(i,j) = Vjiij(i,j) - Rimjc(i,p,j,q)*Gimjc(j,p,i,q)
             enddo
          enddo
       enddo
    enddo
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vjiij(i,j) = Vjiij(i,j) - Rimjc(j,p,i,q)*Gimjc(i,p,j,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vjiij

  !> \brief Vjiij term1 contribution for MP2-f12
  !> \author Yang M. Wang
  !> \date March 2013
  subroutine mp2f12_Vjiij_term1(Vjiij,Fijkl,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(out) :: Vjiij(nocc,nocc)
    Real(realk),intent(in)  :: Fijkl(nocc,nocc,nocc,nocc)
    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vjiij = 0.0E0_realk
    do j=1,nocc
       do i=1,nocc
          Vjiij(i,j) = Fijkl(i,j,j,i)
       enddo
    enddo
  end subroutine mp2f12_Vjiij_term1

  subroutine mp2f12_Vjiij_term2(Vjiij,Ripjq,Gipjq,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(inout) :: Vjiij(nocc,nocc)
    Real(realk),intent(in)  :: Ripjq(nocc,nbasis,nocc,nbasis)
    Real(realk),intent(in)  :: Gipjq(nocc,nbasis,nocc,nbasis)

    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vjiij = 0.0E0_realk
    do q=1,nbasis
       do j=1,nocc
          do p=1,nbasis
             do i=1,nocc
                Vjiij(i,j) = Vjiij(i,j) - Ripjq(i,p,j,q)*Gipjq(j,p,i,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vjiij_term2

  subroutine mp2f12_Vjiij_term3(Vjiij,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(inout) :: Vjiij(nocc,nocc)
    Real(realk),intent(in)  :: Rimjc(nocc,noccfull,nocc,ncabs)
    Real(realk),intent(in)  :: Gimjc(nocc,noccfull,nocc,ncabs)

    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vjiij = 0E0_realk
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vjiij(i,j) = Vjiij(i,j) - Rimjc(i,p,j,q)*Gimjc(j,p,i,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vjiij_term3

  subroutine mp2f12_Vjiij_term4(Vjiij,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    implicit none
    Real(realk),intent(inout) :: Vjiij(nocc,nocc)
    Real(realk),intent(in)  :: Rimjc(nocc,noccfull,nocc,ncabs)
    Real(realk),intent(in)  :: Gimjc(nocc,noccfull,nocc,ncabs)

    Integer,intent(in)      :: nocc,noccfull,nbasis,ncabs
    !
    Integer :: i,j,p,q
    Vjiij = 0.0E0_realk
    do q=1,ncabs
       do j=1,nocc
          do p=1,noccfull
             do i=1,nocc
                Vjiij(i,j) = Vjiij(i,j) - Rimjc(j,p,i,q)*Gimjc(i,p,j,q)
             enddo
          enddo
       enddo
    enddo
  end subroutine mp2f12_Vjiij_term4

  !> \brief Xijij contribution for MP2-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine mp2f12_Xijij(Xijij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
  implicit none
  Real(realk),intent(OUT) :: Xijij(nocc,nocc)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Tijkl(nocc,nocc,nocc,nocc)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs
  !
  Integer :: i,j,p,q

  DO j=1,nocc
    DO i=1,nocc
      Xijij(i,j) = Tijkl(i,i,j,j)
    ENDDO
  ENDDO
  DO q=1,nbasis
    DO j=1,nocc
      DO p=1,nbasis
        DO i=1,nocc
          Xijij(i,j) = Xijij(i,j) - Gipjq(i,p,j,q)*Gipjq(i,p,j,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO q=1,ncabs
    DO j=1,nocc
      DO p=1,noccfull
        DO i=1,nocc
          Xijij(i,j) = Xijij(i,j) - Gimjc(i,p,j,q)*Gimjc(i,p,j,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO q=1,ncabs
    DO j=1,nocc
      DO p=1,noccfull
        DO i=1,nocc
          Xijij(i,j) = Xijij(i,j) - Gimjc(j,p,i,q)*Gimjc(j,p,i,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine mp2f12_Xijij

  !> \brief Xijij contribution for MP2-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine mp2f12_Xijijfull(Xijkl,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
  implicit none
  Real(realk),intent(OUT) :: Xijkl(nocc,nocc,nocc,nocc)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Tijkl(nocc,nocc,nocc,nocc)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs
  !
  Integer :: i,j,p,q,l,k

  DO l=1,nocc
    DO k=1,nocc
       DO j=1,nocc
          DO i=1,nocc
             Xijkl(i,j,k,l) = Tijkl(i,j,k,l)
          ENDDO
       ENDDO
    ENDDO
 ENDDO
  DO l=1,nocc
    DO k=1,nocc
  DO q=1,nbasis
    DO j=1,nocc
      DO p=1,nbasis
        DO i=1,nocc
          Xijkl(i,j,k,l) = Xijkl(i,j,k,l) - Gipjq(i,p,k,q)*Gipjq(j,p,l,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
    ENDDO
 ENDDO
  DO l=1,nocc
    DO k=1,nocc
  DO q=1,ncabs
    DO j=1,nocc
      DO p=1,noccfull
        DO i=1,nocc
          Xijkl(i,j,k,l) = Xijkl(i,j,k,l) - Gimjc(i,p,k,q)*Gimjc(j,p,l,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
    ENDDO
 ENDDO
  DO l=1,nocc
    DO k=1,nocc
  DO q=1,ncabs
    DO j=1,nocc
      DO p=1,noccfull
        DO i=1,nocc
          Xijkl(i,j,k,l) = Xijkl(i,j,k,l) - Gimjc(k,p,i,q)*Gimjc(l,p,j,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
    ENDDO
 ENDDO
end subroutine mp2f12_Xijijfull

  !> \brief Xjiij contribution for MP2-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine mp2f12_Xjiij(Xjiij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
  implicit none
  Real(realk),intent(OUT) :: Xjiij(nocc,nocc)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Tijkl(nocc,nocc,nocc,nocc)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs
  !
  Integer :: i,j,p,q

  DO j=1,nocc
    DO i=1,nocc
      Xjiij(i,j) = Tijkl(i,j,j,i)
    ENDDO
  ENDDO
  DO q=1,nbasis
    DO j=1,nocc
      DO p=1,nbasis
        DO i=1,nocc
          Xjiij(i,j) = Xjiij(i,j) - Gipjq(i,p,j,q)*Gipjq(j,p,i,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO q=1,ncabs
    DO j=1,nocc
      DO p=1,noccfull
        DO i=1,nocc
          Xjiij(i,j) = Xjiij(i,j) - Gimjc(i,p,j,q)*Gimjc(j,p,i,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO q=1,ncabs
    DO j=1,nocc
      DO p=1,noccfull
        DO i=1,nocc
          Xjiij(i,j) = Xjiij(i,j) - Gimjc(j,p,i,q)*Gimjc(i,p,j,q)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine mp2f12_Xjiij

  !> \brief Bijij contribution for MP2-f12
  !> \author T Kjaergaard
  !> \date May 2012
  subroutine mp2f12_Bijij(Bijij,Dijkl,Tirjk,Tijkr,Girjs,hJ,K,Frr,Fpp,Fmm,Frm,Fcp,&
       & Girjm,Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,&
       & nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs)
  implicit none
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs
  Real(realk),intent(OUT) :: Bijij(nocc,nocc)
  Real(realk),intent(IN)  :: Dijkl(nocc,nocc,nocc,nocc) !The double commutator [[T,g],g]
  Real(realk),intent(IN)  :: Tirjk(nocc,ncabsAO,nocc,nocc) !The Gaussian geminal operator squared g^2
  Real(realk),intent(IN)  :: Tijkr(nocc,nocc,nocc,ncabsAO) !The Gaussian geminal operator squared g^2
  Real(realk),intent(IN)  :: Girjs(nocc,ncabsAO,nocc,ncabsAO) !The Gaussian geminal operator g
  Real(realk),intent(IN)  :: Girjm(nocc,ncabsAO,nocc,noccfull)
  Real(realk),intent(IN)  :: Grimj(ncabsAO,nocc,noccfull,nocc)
  Real(realk),intent(IN)  :: Gipja(nocc,nbasis,nocc,nvirt)
  Real(realk),intent(IN)  :: Gpiaj(nbasis,nocc,nvirt,nocc)

  Real(realk),intent(IN)  :: Gicjm(nocc,ncabs,nocc,noccfull)
  Real(realk),intent(IN)  :: Gcimj(ncabs,nocc,noccfull,nocc)
  Real(realk),intent(IN)  :: Gcirj(ncabs,nocc,ncabsAO,nocc)
  Real(realk),intent(IN)  :: Gciaj(ncabs,nocc,nvirt,nocc)

  Real(realk),intent(IN)  :: hJ(nocc,ncabsAO) !(h+J)
  Real(realk),intent(IN)  :: K(ncabsAO,ncabsAO) !exchange matrix
  Real(realk),intent(IN)  :: Frr(ncabsAO,ncabsAO) !Fock matrix in RI MO basis
  Real(realk),intent(IN)  :: Fpp(nbasis,nbasis) !Fock matrix in full MO basis
  Real(realk),intent(IN)  :: Fmm(noccfull,noccfull) !Fock matrix in full MO occ basis
  Real(realk),intent(IN)  :: Frm(ncabsAO,noccfull)
  Real(realk),intent(IN)  :: Fcp(ncabs,nbasis)
  real(realk),parameter :: D05=0.5E0_realk,D2=2.0E0_realk
  !
  Integer :: i,j,p,q,r,m,a,n,b

  DO j=1,nocc
    DO i=1,nocc
      Bijij(i,j) = Dijkl(i,i,j,j)
    ENDDO
  ENDDO
  DO j=1,nocc
     DO p=1,ncabsAO
        DO i=1,nocc
           Bijij(i,j) = Bijij(i,j) + Tirjk(i,p,j,j)*hJ(i,p)
        ENDDO
     ENDDO
  ENDDO
  DO j=1,nocc
     DO p=1,ncabsAO
        DO i=1,nocc
           Bijij(i,j) = Bijij(i,j) + Tijkr(i,i,j,p)*hJ(j,p)
        ENDDO
     ENDDO
  ENDDO
!Krr
  DO r=1,ncabsAO
     DO q=1,ncabsAO
        DO j=1,nocc
           DO p=1,ncabsAO
              DO i=1,nocc
                 Bijij(i,j) = Bijij(i,j) + (Girjs(i,r,j,p)*K(q,p)*Girjs(i,r,j,q) &
                      &  + Girjs(i,p,j,r)*K(q,p)*Girjs(i,q,j,r))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!!$  !fourth term first Fock
  DO m=1,noccfull
     DO q=1,ncabsAO
        DO j=1,nocc
           DO p=1,ncabsAO
              DO i=1,nocc
                 Bijij(i,j) = Bijij(i,j) - (Girjm(i,p,j,m)*Frr(q,p)*Grimj(q,i,m,j) &
                      & + Girjm(j,p,i,m)*Frr(q,p)*Grimj(q,j,m,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!!$  !Fpp
  DO a=1,nvirt
     DO q=1,nbasis
        DO j=1,nocc
           DO p=1,nbasis
              DO i=1,nocc
                 Bijij(i,j) = Bijij(i,j) - (Gipja(i,p,j,a)*Fpp(q,p)*Gpiaj(q,i,a,j) &
                      & + Gipja(j,p,i,a)*Fpp(q,p)*Gpiaj(q,j,a,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !Fmm
  DO n=1,noccfull
     DO m=1,noccfull
        DO j=1,nocc
           DO a=1,ncabs
              DO i=1,nocc
                 Bijij(i,j) = Bijij(i,j) + (Gicjm(i,a,j,m)*Fmm(n,m)*Gcimj(a,i,n,j) &
                      & + Gicjm(j,a,i,m)*Fmm(n,m)*Gcimj(a,j,n,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!!$  !Frm
  DO p=1,ncabsAO
     DO m=1,noccfull
        DO j=1,nocc
           DO a=1,ncabs
              DO i=1,nocc
                 Bijij(i,j) = Bijij(i,j) - D2*(Gicjm(i,a,j,m)*Frm(p,m)*Gcirj(a,i,p,j) &
                      &  + Gicjm(j,a,i,m)*Frm(p,m)*Gcirj(a,j,p,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !Fcp
  DO p=1,nbasis
     DO b=1,nvirt
        DO j=1,nocc
           DO a=1,ncabs
              DO i=1,nocc
                 Bijij(i,j) = Bijij(i,j) - D2*(Gipja(i,p,j,b)*Fcp(a,p)*Gciaj(a,i,b,j) &
                      & + Gipja(j,p,i,b)*Fcp(a,p)*Gciaj(a,j,b,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  end subroutine mp2f12_Bijij

  !maybe i and j is reverted
  !> \brief Viija contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viija(Viija,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Viija(nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Ripaq(nocc,nbasis,nvirt,nbasis)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Fijka(nocc,nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Rimac(nocc,noccfull,nvirt,ncabs)
  Real(realk),intent(IN)  :: Ramic(nvirt,noccfull,nocc,ncabs)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  !
  Integer :: i,j,a

  call ccsdf12_Viija0(Viija,Fijka,nocc,nvirt)
  call ccsdf12_Viija1(Viija,Ripaq,Gipjq,nocc,noccfull,nbasis,ncabs,nvirt)
  call ccsdf12_Viija2(Viija,Rimac,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  call ccsdf12_Viija3(Viija,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

  end subroutine ccsdf12_Viija

  !> \brief Viija contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viija0(Viija,Fijka,nocc,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Viija(nocc,nocc,nvirt)
  Real(realk),intent(IN)    :: Fijka(nocc,nocc,nocc,nvirt)
  Integer,intent(IN)        :: nocc,nvirt
!
  Integer :: i,j,a
  !This ordering seems (and is) inefficient but due to a 
  !compiler bug we had to write it this way
  !my theory is that it tries to vectorize the i loop
  !if this loop is the inner most loop - but it should not due to the
  !non unit stride in Fijka(i,i,j,a)
  DO i=1,nocc
     DO a=1,nvirt
        DO j=1,nocc
           Viija(i,j,a) = Fijka(i,i,j,a)
        ENDDO
     ENDDO
  ENDDO
  end subroutine ccsdf12_Viija0

  !> \brief Viija contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viija1(Viija,Ripaq,Gipjq,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Viija(nocc,nocc,nvirt)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  Real(realk),intent(IN)  :: Ripaq(nocc,nbasis,nvirt,nbasis)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
!
  Integer :: i,j,p,q,a,b,m,c
  real(realk) :: tmp(nocc)
  DO a=1,nvirt
   DO j=1,nocc
    DO i=1,nocc
       tmp(i) = 0.0E0_realk
    ENDDO
    DO q=1,nbasis
     DO p=1,nbasis
      DO i=1,nocc
       tmp(i) = tmp(i) + Ripaq(i,p,a,q)*Gipjq(i,p,j,q)
      ENDDO
     ENDDO
    ENDDO
    DO i=1,nocc
     Viija(i,j,a) = Viija(i,j,a) - tmp(i)
    ENDDO
   ENDDO
  ENDDO
  end subroutine ccsdf12_Viija1

  !> \brief Viija contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viija2(Viija,Rimac,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Viija(nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Rimac(nocc,noccfull,nvirt,ncabs)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  !
  Integer :: i,j,p,q,a,b,m,c
  real(realk) :: tmp(nocc)
  DO a=1,nvirt
   DO j=1,nocc
    DO i=1,nocc
       tmp(i) = 0.0E0_realk
    ENDDO
    DO c=1,ncabs
     DO m=1,noccfull
      DO i=1,nocc
         tmp(i) = tmp(i) + Rimac(i,m,a,c)*Gimjc(i,m,j,c) 
      ENDDO
     ENDDO
    ENDDO
    DO i=1,nocc
       Viija(i,j,a) = Viija(i,j,a) - tmp(i)
    ENDDO
   ENDDO
  ENDDO
  end subroutine ccsdf12_Viija2

  !> \brief Viija contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viija3(Viija,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Viija(nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Ramic(nvirt,noccfull,nocc,ncabs)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  !
  Integer :: i,j,p,q,a,b,m,c
  real(realk) :: tmp
  DO a=1,nvirt
   DO j=1,nocc
    DO i=1,nocc
     tmp = 0.0E0_realk
     DO c=1,ncabs
      DO m=1,noccfull
       tmp = tmp + Ramic(a,m,i,c)*Gimjc(j,m,i,c)
      ENDDO
     ENDDO
     Viija(i,j,a) = Viija(i,j,a) - tmp
    ENDDO
   ENDDO
  ENDDO
end subroutine ccsdf12_Viija3

  !> \brief Viajj contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viajj(Viajj,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Viajj(nocc,nvirt,nocc)
  Real(realk),intent(IN)  :: Ripaq(nocc,nbasis,nvirt,nbasis)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Fijka(nocc,nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Rimac(nocc,noccfull,nvirt,ncabs)
  Real(realk),intent(IN)  :: Ramic(nvirt,noccfull,nocc,ncabs)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  !
  Integer :: i,j,p,q,a,b,m,c

  DO a=1,nvirt
   DO j=1,nocc
    DO i=1,nocc
     Viajj(i,a,j) = Fijka(j,j,i,a)
    ENDDO
   ENDDO
  ENDDO
  DO q=1,nbasis
   DO a=1,nvirt
    DO j=1,nocc
     DO p=1,nbasis
      DO i=1,nocc
       Viajj(i,a,j) = Viajj(i,a,j) - Ripaq(j,q,a,p)*Gipjq(i,p,j,q)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  DO c=1,ncabs
   DO a=1,nvirt
    DO j=1,nocc
     DO m=1,noccfull
      DO i=1,nocc
       Viajj(i,a,j) = Viajj(i,a,j) - Ramic(a,m,j,c)*Gimjc(i,m,j,c)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  DO c=1,ncabs
   DO a=1,nvirt
    DO j=1,nocc
     DO m=1,noccfull
      DO i=1,nocc
       Viajj(i,a,j) = Viajj(i,a,j) - Rimac(j,m,a,c)*Gimjc(j,m,i,c)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  end subroutine ccsdf12_Viajj

  !> \brief Viaji contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viaji(Viaji,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Viaji(nocc,nvirt,nocc)
  Real(realk),intent(IN)  :: Ripaq(nocc,nbasis,nvirt,nbasis)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Fijka(nocc,nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Rimac(nocc,noccfull,nvirt,ncabs)
  Real(realk),intent(IN)  :: Ramic(nvirt,noccfull,nocc,ncabs)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  !
  Integer :: i,j,p,q,a,b,c,m

  DO a=1,nvirt
   DO j=1,nocc
    DO i=1,nocc
     Viaji(i,a,j) = Fijka(j,i,i,a)
    ENDDO
   ENDDO
  ENDDO
  DO q=1,nbasis
   DO a=1,nvirt
    DO j=1,nocc
     DO p=1,nbasis
      DO i=1,nocc
       Viaji(i,a,j) = Viaji(i,a,j) - Ripaq(i,q,a,p)*Gipjq(i,p,j,q)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  DO c=1,ncabs
   DO a=1,nvirt
    DO j=1,nocc
     DO m=1,noccfull
      DO i=1,nocc
       Viaji(i,a,j) = Viaji(i,a,j) - Ramic(a,m,i,c)*Gimjc(i,m,j,c)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  DO c=1,ncabs
   DO a=1,nvirt
    DO j=1,nocc
     DO m=1,noccfull
      DO i=1,nocc
       Viaji(i,a,j) = Viaji(i,a,j) - Rimac(i,m,a,c)*Gimjc(j,m,i,c)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  end subroutine ccsdf12_Viaji

  !> \brief Viija contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Vijja(Vijja,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(INOUT) :: Vijja(nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Ripaq(nocc,nbasis,nvirt,nbasis)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Fijka(nocc,nocc,nocc,nvirt)
  Real(realk),intent(IN)  :: Rimac(nocc,noccfull,nvirt,ncabs)
  Real(realk),intent(IN)  :: Ramic(nvirt,noccfull,nocc,ncabs)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  !
  Integer :: i,j,p,q,a,b,c,m,k
  real(realk) :: tmp(nocc),tmp2
  DO a=1,nvirt
     DO j=1,nocc
        DO i=1,nocc
           Vijja(i,j,a) = Fijka(i,j,j,a)
        ENDDO
     ENDDO
  ENDDO
  DO a=1,nvirt
     DO j=1,nocc
        DO i=1,nocc
           tmp(i) = 0.0E0_realk
        ENDDO
        DO q=1,nbasis
           DO p=1,nbasis
              tmp2 = Ripaq(j,p,a,q)
              DO i=1,nocc
                 tmp(i) = tmp(i) + tmp2*Gipjq(i,p,j,q)
              ENDDO
           ENDDO
        ENDDO
        DO i=1,nocc
           Vijja(i,j,a) = Vijja(i,j,a) - tmp(i)
        ENDDO
     ENDDO
  ENDDO
  DO a=1,nvirt
     DO j=1,nocc
        DO i=1,nocc
           tmp(i) = 0.0E0_realk
        ENDDO
        DO c=1,ncabs
           DO m=1,noccfull
              tmp2 = Rimac(j,m,a,c)
              DO i=1,nocc
                 tmp(i) = tmp(i) + tmp2*Gimjc(i,m,j,c)
              ENDDO
           ENDDO
        ENDDO
        DO i=1,nocc
           Vijja(i,j,a) = Vijja(i,j,a) - tmp(i)
        ENDDO
     ENDDO
  ENDDO
  DO a=1,nvirt
     DO i=1,nocc
        DO j=1,nocc
           tmp(j) = 0.0E0_realk
        ENDDO
        DO c=1,ncabs
           DO m=1,noccfull
              DO j=1,nocc
                 tmp(j) = tmp(j) + Ramic(a,m,j,c)*Gimjc(j,m,i,c)
              ENDDO
           ENDDO
        ENDDO
        DO j=1,nocc
           Vijja(i,j,a) = Vijja(i,j,a) - tmp(j)
        ENDDO
     ENDDO
  ENDDO
  end subroutine ccsdf12_Vijja

  !> \brief Viajb contribution for CCSD-f12
  !> \author Simen Reine
  !> \date May 2012
  subroutine ccsdf12_Viajb(Viajb,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
  implicit none
  Real(realk),intent(OUT) :: Viajb(nocc,nvirt,nocc,nvirt)
  Real(realk),intent(IN)  :: Rapbq(nvirt,nbasis,nvirt,nbasis)
  Real(realk),intent(IN)  :: Gipjq(nocc,nbasis,nocc,nbasis)
  Real(realk),intent(IN)  :: Fiajb(nocc,nvirt,nocc,nvirt)
  Real(realk),intent(IN)  :: Rambc(nvirt,noccfull,nvirt,ncabs)
  Real(realk),intent(IN)  :: Gimjc(nocc,noccfull,nocc,ncabs)
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabs,nvirt
  !
  Integer :: i,j,p,q,a,b,m,c

  DO b=1,nvirt
   DO j=1,nocc
    DO a=1,nvirt
     DO i=1,nocc
      Viajb(i,a,j,b) = Fiajb(i,a,j,b)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  DO q=1,nbasis
   DO b=1,nvirt
    DO j=1,nocc
     DO p=1,nbasis
      DO a=1,nvirt
       DO i=1,nocc
        Viajb(i,a,j,b) = Viajb(i,a,j,b) - Rapbq(a,p,b,q)*Gipjq(i,p,j,q)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  DO c=1,ncabs
   DO b=1,nvirt
    DO j=1,nocc
     DO m=1,noccfull
      DO a=1,nvirt
       DO i=1,nocc
        Viajb(i,a,j,b) = Viajb(i,a,j,b) - Rambc(a,m,b,c)*Gimjc(i,m,j,c)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  DO c=1,ncabs
   DO b=1,nvirt
    DO j=1,nocc
     DO m=1,noccfull
      DO a=1,nvirt
       DO i=1,nocc
        Viajb(i,a,j,b) = Viajb(i,a,j,b) - Rambc(b,m,a,c)*Gimjc(j,m,i,c)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  end subroutine ccsdf12_Viajb

  subroutine mp2f12_Cjaib(Cjaib,Givic,Fvc,nocc,nvirt,ncabs)
  implicit none
  real(realk), intent(INOUT) :: Cjaib(nocc,nvirt,nocc,nvirt)
  real(realk), intent(IN)    :: Givic(nocc,nvirt,nocc,ncabs)
  real(realk), intent(IN)    :: Fvc(nvirt,ncabs)
  integer, intent(IN)        :: nocc,nvirt,ncabs
  !
  integer :: i,j,a,b,c

  Cjaib = 0E0_realk
  DO c=1,ncabs
   DO b=1,nvirt
    DO j=1,nocc
     DO a=1,nvirt
      DO i=1,nocc
       Cjaib(i,a,j,b) = Cjaib(i,a,j,b) + Givic(j,a,i,c)*Fvc(b,c) &
    &                                  + Givic(i,b,j,c)*Fvc(a,c)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  end subroutine mp2f12_Cjaib

  subroutine mp2f12_Ciajb(Ciajb,Givic,Fvc,nocc,nvirt,ncabs)
  implicit none
  real(realk), intent(INOUT) :: Ciajb(nocc,nvirt,nocc,nvirt)
  real(realk), intent(IN)    :: Givic(nocc,nvirt,nocc,ncabs)
  real(realk), intent(IN)    :: Fvc(nvirt,ncabs)
  integer, intent(IN)        :: nocc,nvirt,ncabs
  !
  integer :: i,j,a,b,c

  Ciajb = 0E0_realk
  DO c=1,ncabs
   DO b=1,nvirt
    DO j=1,nocc
     DO a=1,nvirt
      DO i=1,nocc
       Ciajb(i,a,j,b) = Ciajb(i,a,j,b) + Givic(i,a,j,c)*Fvc(b,c) &
    &                                  + Givic(j,b,i,c)*Fvc(a,c)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  end subroutine mp2f12_Ciajb


  !> \brief Bjiij contribution for MP2-f12
  !> \author S Reine and T Kjaegaard and Edvard Valeev
  !> \date May 2012
  subroutine mp2f12_Bjiij(Bjiij,Dijkl,Tirjk,Tijkr,Girjs,hJ,K,Frr,Fpp,Fmm,Frm,Fcp,&
       & Girjm,Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,&
       & nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs)
  implicit none
  Integer,intent(IN)      :: nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs
  Real(realk),intent(OUT) :: Bjiij(nocc,nocc)
  Real(realk),intent(IN)  :: Dijkl(nocc,nocc,nocc,nocc) !The double commutator [[T,g],g]
  Real(realk),intent(IN)  :: Tirjk(nocc,ncabsAO,nocc,nocc) !The Gaussian geminal operator squared g^2
  Real(realk),intent(IN)  :: Tijkr(nocc,nocc,nocc,ncabsAO) !The Gaussian geminal operator squared g^2
  Real(realk),intent(IN)  :: Girjs(nocc,ncabsAO,nocc,ncabsAO) !The Gaussian geminal operator g
  Real(realk),intent(IN)  :: Girjm(nocc,ncabsAO,nocc,noccfull)
  Real(realk),intent(IN)  :: Grimj(ncabsAO,nocc,noccfull,nocc)
  Real(realk),intent(IN)  :: Gipja(nocc,nbasis,nocc,nvirt)
  Real(realk),intent(IN)  :: Gpiaj(nbasis,nocc,nvirt,nocc)

  Real(realk),intent(IN)  :: Gicjm(nocc,ncabs,nocc,noccfull)
  Real(realk),intent(IN)  :: Gcimj(ncabs,nocc,noccfull,nocc)
  Real(realk),intent(IN)  :: Gcirj(ncabs,nocc,ncabsAO,nocc)
  Real(realk),intent(IN)  :: Gciaj(ncabs,nocc,nvirt,nocc)

  Real(realk),intent(IN)  :: hJ(nocc,ncabsAO) !(h+J)
  Real(realk),intent(IN)  :: K(ncabsAO,ncabsAO) !exchange matrix
  Real(realk),intent(IN)  :: Frr(ncabsAO,ncabsAO) !Fock matrix in RI MO basis
  Real(realk),intent(IN)  :: Fpp(nbasis,nbasis) !Fock matrix in full MO basis
  Real(realk),intent(IN)  :: Fmm(noccfull,noccfull) !Fock matrix in full MO occ basis
  Real(realk),intent(IN)  :: Frm(ncabsAO,noccfull)
  Real(realk),intent(IN)  :: Fcp(ncabs,nbasis)
  real(realk),parameter :: D05=0.5E0_realk,D2=2.0E0_realk
  !
  Integer :: i,j,p,q,r,m,a,n,b

  DO j=1,nocc
    DO i=1,nocc
      Bjiij(i,j) = Dijkl(i,j,i,j)
    ENDDO
  ENDDO
  DO j=1,nocc
     DO p=1,ncabsAO
        DO i=1,nocc
           Bjiij(i,j) = Bjiij(i,j) + Tirjk(j,p,i,j)*hJ(i,p)
        ENDDO
     ENDDO
  ENDDO
  DO j=1,nocc
     DO p=1,ncabsAO
        DO i=1,nocc
           Bjiij(i,j) = Bjiij(i,j) + Tijkr(i,j,i,p)*hJ(j,p)
        ENDDO
     ENDDO
  ENDDO
  DO r=1,ncabsAO
     DO q=1,ncabsAO
        DO j=1,nocc
           DO p=1,ncabsAO
              DO i=1,nocc
                 Bjiij(i,j) = Bjiij(i,j) + (Girjs(j,r,i,p)*K(q,p)*Girjs(i,r,j,q) &
                      & + Girjs(j,p,i,r)*K(q,p)*Girjs(i,q,j,r))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!!$  !fourth term first Fock
  DO m=1,noccfull
     DO q=1,ncabsAO
        DO j=1,nocc
           DO p=1,ncabsAO
              DO i=1,nocc
                 Bjiij(i,j) = Bjiij(i,j) - (Girjm(j,p,i,m)*Frr(q,p)*Grimj(q,i,m,j) &
                      & + Girjm(i,p,j,m)*Frr(q,p)*Grimj(q,j,m,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!!$  !Fpp
  DO a=1,nvirt
     DO q=1,nbasis
        DO j=1,nocc
           DO p=1,nbasis
              DO i=1,nocc
                 Bjiij(i,j) = Bjiij(i,j) - (Gipja(j,p,i,a)*Fpp(q,p)*Gpiaj(q,i,a,j) &
                      & + Gipja(i,p,j,a)*Fpp(q,p)*Gpiaj(q,j,a,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !Fmm
  DO n=1,noccfull
     DO m=1,noccfull
        DO j=1,nocc
           DO a=1,ncabs
              DO i=1,nocc
                 Bjiij(i,j) = Bjiij(i,j) + (Gicjm(j,a,i,m)*Fmm(n,m)*Gcimj(a,i,n,j) &
                      & + Gicjm(i,a,j,m)*Fmm(n,m)*Gcimj(a,j,n,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!!$  !Frm
  DO p=1,ncabsAO
     DO m=1,noccfull
        DO j=1,nocc
           DO a=1,ncabs
              DO i=1,nocc
                 Bjiij(i,j) = Bjiij(i,j) - D2*(Gicjm(j,a,i,m)*Frm(p,m)*Gcirj(a,i,p,j) &
                      & + Gicjm(i,a,j,m)*Frm(p,m)*Gcirj(a,j,p,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !Fcp
  DO p=1,nbasis
     DO b=1,nvirt
        DO j=1,nocc
           DO a=1,ncabs
              DO i=1,nocc
                 Bjiij(i,j) = Bjiij(i,j) - D2*(Gipja(j,p,i,b)*Fcp(a,p)*Gciaj(a,i,b,j) &
                      & + Gipja(i,p,j,b)*Fcp(a,p)*Gciaj(a,j,b,i))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  end subroutine mp2f12_Bjiij
#else
  subroutine full_f12contractions_dummy_sub12()
    implicit none
  end subroutine full_f12contractions_dummy_sub12
#endif
end module full_f12contractions

