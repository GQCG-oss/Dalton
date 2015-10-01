!> Module with files used for grid generation. Responsible people, please document these subroutines!
!> \author Variuos people, moved here by Kasper Kristensen

Module grid_utilities_module
Use precision
use TYPEDEFTYPE
use files
use basis_typetype
contains 


  !> \brief Determine molecule-specific gridbox to use for calculation of 
  !> orbital, density, electrostatic potential etc. at given points in space.
  !> \author Kasper Kristensen
  !> \date 2013
  subroutine DETERMINE_GRIDBOX(natoms,ATOMXYZ,MyPlt)
    implicit none
    !> Number of atoms in molecule
    integer, intent(in)          :: natoms
    !> (X,Y,Z) coordinates for all atoms in molecule
    real(4), intent(in)      :: ATOMXYZ(3,natoms)
    !> PLT info where gridbox is set
    type(pltinfo),intent(inout) :: MyPlt
    real(4)  :: deltax,deltay,deltaz,mybuffer
    real(4)                  :: Xn, Yn, Zn,distX,distY,distZ
    integer                      :: I


    ! *******************
    ! *    GRID BOX     *
    ! *******************
    ! * The first point in the grid box is (X1,Y1,Z1).
    ! * The remaining grid points are then defined by going out in the X,Y, and Z directions
    !   with step sizes deltax,deltay, and deltaz, until there are nX,nY, and nZ points
    !   in the X,Y, and Z directions (giving a total number of gridpoints: nGRIDPOINTS=nX*nY*nZ).
    !
    ! We define the grid box parameters such that (i) all atoms in the molecule are contained
    ! within the grid box, and (ii) there is a buffer zone (mybuffer) around the outermost atoms.

    ! Set distances between gridpoints from input
    deltax = MyPlt%deltax
    deltay = MyPlt%deltay
    deltaz = MyPlt%deltaz

    ! Set buffer zone from input
    mybuffer = MyPlt%buffer

    ! Minimum and maximum values in the gridbox. 
    Xn = -HUGE(1_4); MyPlt%X1 = HUGE(1_4)
    Yn = -HUGE(1_4); MyPlt%Y1 = HUGE(1_4)
    Zn = -HUGE(1_4); MyPlt%Z1 = HUGE(1_4)
    do I = 1, natoms
       MyPlt%X1=min(ATOMXYZ(1,I),MyPlt%X1); Xn=max(ATOMXYZ(1,I),Xn)
       MyPlt%Y1=min(ATOMXYZ(2,I),MyPlt%Y1); Yn=max(ATOMXYZ(2,I),Yn)
       MyPlt%Z1=min(ATOMXYZ(3,I),MyPlt%Z1); Zn=max(ATOMXYZ(3,I),Zn)
    enddo

    ! Subtract/add buffer from minimum/maximum X,Y,Z values
    MyPlt%X1 = MyPlt%X1 - mybuffer; MyPlt%Y1 = MyPlt%Y1 - mybuffer; MyPlt%Z1 = MyPlt%Z1 - mybuffer
    Xn = Xn + mybuffer; Yn = Yn + mybuffer; Zn = Zn + mybuffer

    ! Distances from min to max values
    distX = Xn-MyPlt%X1
    distY = Yn-MyPlt%Y1
    distZ = Zn-MyPlt%Z1

    ! Number of points
    MyPlt%nX = nint(distX/deltax) + 1
    MyPlt%nY = nint(distY/deltay) + 1
    MyPlt%nZ = nint(distZ/deltaz) + 1

    ! Total number of gridpoints
    MyPlt%ngridpoints = MyPlt%nX*MyPlt%nY*MyPlt%nZ

  end subroutine DETERMINE_GRIDBOX



  subroutine determine_orbitals(I,GAO,orbnr, nORBITALS, X, Y, Z, ls)
    ! Determine value of all atomix orbitals on atom I at point (X,Y,Z)
    ! measured relative to atom center.
    ! The values are stored in the GAO vector.
    implicit none
    type(lsitem)                   :: ls
    integer, intent(in)         :: I, nORBITALS
    integer, intent(inout)       :: orbnr
    real(4), intent(in)     :: X,Y,Z
    real(4), intent(inout)  :: GAO(nORBITALS)
    ! Local
    integer                     :: basindex,type,set,nAngmom,J,K,L,M,nsegments,P2
    integer                     :: nEXPONENTS,nORBJ,P,Q,R,icharge,itype,angmom,PN,spsize,MMM
    real(4)                 :: coeff,exponent,R2,contracted_gaussians,GAX,GAY,GAZ,GAXX,GAYY,CINT,STMP
    integer, allocatable        :: AVALUE(:), BVALUE(:),CVALUE(:)
    real(4), allocatable           :: fxyz(:)
    real(4), allocatable           :: sphmat2(:,:)
#if 1
    real(4), parameter          :: SPHMAT(9,15) = reshape((/ &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.36596254E-01, &
& 0.00000000E+00,-0.54554474E-01, 0.00000000E+00, 0.72168790E-01, 0.28867516E+00, &
& 0.00000000E+00,-0.10910895E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,-0.23145503E+00, 0.00000000E+00, &
& 0.20412414E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.73192507E-01, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
&-0.43301272E+00, 0.00000000E+00, 0.61237240E+00, 0.00000000E+00,-0.23145503E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,-0.29277003E+00, &
& 0.00000000E+00, 0.32732683E+00, 0.00000000E+00, 0.00000000E+00,-0.28867516E+00, &
& 0.00000000E+00,-0.10910895E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,-0.23145503E+00, 0.00000000E+00, &
&-0.61237240E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.65465367E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.30860671E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.36596254E-01, &
& 0.00000000E+00, 0.54554474E-01, 0.00000000E+00, 0.72168790E-01, 0.00000000E+00, &
&-0.20412414E+00, 0.00000000E+00,-0.23145503E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00,-0.29277003E+00, 0.00000000E+00,-0.32732683E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.30860671E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
& 0.97590007E-01, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00 /),(/9,15/))
#endif

if(ls%setting%IntegralTransformGC) then
stop 'determine_orbitals: MOs and AOs are not expressed in the same basis!'
end if
    ! total distance from gridpoint to atom I, squared
    R2 = X**2+Y**2+Z**2
    ! basindex = 1 or 2 (2 for auxillary basis sets)
    ! set: integer for describing which basis set is used for atom I.
    ! type: integer for describing the atomtype of atom I.
    IF(LS%setting%basis(1)%p%BINFO(RegBasParam)%labelindex .EQ. 0)THEN
       icharge = INT(ls%setting%MOLECULE(1)%p%ATOM(i)%charge) 
       itype = ls%setting%basis(1)%p%BINFO(RegBasParam)%chargeindex(icharge)
    ELSE
       itype = ls%setting%MOLECULE(1)%p%ATOM(i)%IDtype(1)
    ENDIF
    nAngmom = LS%setting%BASIS(1)%p%BINFO(RegBasParam)%ATOMTYPE(itype)%nAngmom
    DO J =1,nAngmom ! Loop over angular moments
       nsegments=LS%setting%BASIS(1)%p%BINFO(RegBasParam)%ATOMTYPE(itype)%SHELL(J)%nsegments
       DO K = 1,nsegments ! Loop over segments in basis set
          nEXPONENTS = LS%setting%BASIS(1)%p%BINFO(RegBasParam)%ATOMTYPE(itype)%SHELL(J)%segment(K)%nrow
          nORBJ = LS%setting%BASIS(1)%p%BINFO(RegBasParam)%ATOMTYPE(itype)%SHELL(J)%segment(K)%ncol
          DO L = 1,nORBJ 
             DO M = 1,nEXPONENTS 
                coeff = LS%setting%basis(1)%p%BINFO(RegBasParam)%ATOMTYPE(itype)%SHELL(J)&
                     &%segment(K)%elms(M+(L-1)*nEXPONENTS)
                exponent = LS%setting%basis(1)%p%BINFO(RegBasParam)%ATOMTYPE(itype)%SHELL(J)&
                     &%segment(K)%exponents(M)
                GAO(orbnr)=GAO(orbnr)+coeff*EXP(-exponent*R2)
             ENDDO
             IF (J.EQ. 1) THEN !1 S-orbital
                ! Orbital is already constructed and stored in GAO(orbnr).
                orbnr=orbnr+1
             ELSEIF (J.EQ. 2) THEN ! 3 P-orbitals
                ! Multiply orbital by x, y, or z for px, py, and pz orbitals.
                contracted_gaussians = GAO(orbnr)
                GAO(orbnr)=X*contracted_gaussians
                orbnr=orbnr+1
                GAO(orbnr)=Y*contracted_gaussians
                orbnr=orbnr+1
                GAO(orbnr)=Z*contracted_gaussians
                orbnr=orbnr+1
             ELSEIF (J.EQ. 3) THEN ! 5 D spherical orbitals
                GAX  = X*GAO(orbnr)
                GAY  = Y*GAO(orbnr)
                GAZ  = Z*GAO(orbnr)

                GAXX = X*GAX
                GAYY = Y*GAY

                GAO(orbnr) = Y*GAX
                orbnr=orbnr+1

                GAO(orbnr) = Y*GAZ
                orbnr=orbnr+1

                GAO(orbnr) = -0.288675134594813*(GAXX +GAYY)&
                     &+ 0.577350269189626*Z*GAZ
                orbnr=orbnr+1

                GAO(orbnr) = X*GAZ
                orbnr=orbnr+1

                GAO(orbnr) = 0.5*GAXX - 0.5*GAYY
                orbnr=orbnr+1

             ELSEIF (J.EQ. 4) THEN ! 7 F-orbitals
                ! Building the cartesian f-orbitals, fxyz
                ALLOCATE(AVALUE(7*7), BVALUE(7*7), CVALUE(7*7))
                ALLOCATE(fxyz(10))
                ! Constructing the xyz-power of the 10 cartesian f-orbitals
                ! X**AVALUE * Y**BVALUE * Z*CVALUE
                P=0
                DO R = 1,4
                   DO Q = 1,R
                      P=P+1
                      AVALUE(P)=4-R
                      BVALUE(P)=R-Q
                      CVALUE(P)=Q-1
                   ENDDO
                ENDDO
                DO Q = 1, 10
                   fxyz(Q) =&
                        & (X**AVALUE(Q))*(Y**BVALUE(Q))*(Z**CVALUE(Q))*GAO(orbnr)
                ENDDO
                ! The 7 spheric f-orbitals are found as
                ! linear combinations of the "xyz"-f orbitals:
                ! (F1, f2, ..., f7) = (fxyz1, fxyz2, ..., fxyz10) * TM (dimension 10x7)
                ! The transformation matrix TM contains many zeros.
                ! Therefore, the transformations for the 7 spheric 
                ! f-orbitals are written out explicitly.
                GAO(orbnr) = 0.612372435695794*fxyz(2) - 0.204124145231932*fxyz(7)
                orbnr=orbnr+1
                GAO(orbnr) = fxyz(5)
                orbnr=orbnr+1
                GAO(orbnr) =  -0.158113883008419*(fxyz(2)+fxyz(7)) + 0.632455532033676*fxyz(9)
                orbnr=orbnr+1
                GAO(orbnr) =  -0.387298334620742*(fxyz(3)+fxyz(8)) + 0.258198889747161*fxyz(10)
                orbnr=orbnr+1
                GAO(orbnr) = -0.158113883008419*(fxyz(1)+fxyz(4)) + 0.632455532033676*fxyz(6)
                orbnr=orbnr+1
                GAO(orbnr) = 0.5*fxyz(3) - 0.5*fxyz(8)
                orbnr=orbnr+1
                GAO(orbnr) = 0.204124145231932*fxyz(1) - 0.612372435695794*fxyz(4)
                orbnr=orbnr+1
                DEALLOCATE(AVALUE, BVALUE, CVALUE, fxyz)
             ELSEIF (J.EQ. 5) THEN
                !GENERAL CASE. Thomas Kjærgaard
                angmom= J-1
                P = J*(J+1)/2
                PN = 2*angmom+1
                spSIZE=PN*P
!                allocate(SPHMAT2(PN,P))
!                SPHMAT2 = 0_4
!                CALL BUILD_SPHMAT(angmom,spSIZE,SPHMAT2)
                ! Constructing the xyz-power of the P cartesian f-orbitals
                ! X**AVALUE * Y**BVALUE * Z*CVALUE
                MMM = J*(J+1)/2*J*(J+1)/2
                ALLOCATE(AVALUE(MMM), BVALUE(MMM), CVALUE(MMM))
                P2=0
                Rloop23: DO R = 1,P
                   DO Q = 1,R
                      P2=P2+1
                      AVALUE(P2)=J-R
                      BVALUE(P2)=R-Q
                      CVALUE(P2)=Q-1
                      IF(P2.EQ.P)EXIT Rloop23
                   ENDDO
                ENDDO Rloop23
                ALLOCATE(fxyz(P))
                DO Q = 1, P
                   fxyz(Q) = (X**AVALUE(Q))*(Y**BVALUE(Q))*(Z**CVALUE(Q))*GAO(orbnr)
                ENDDO
                DO R = 1, PN                                      
                   Q = 1
                   GAO(orbnr+R-1) = SPHMAT(R,Q)*fxyz(Q)
                   DO Q = 2, P
                      GAO(orbnr+R-1) = GAO(orbnr+R-1) + SPHMAT(R,Q)*fxyz(Q)
                   ENDDO
                END DO
                DEALLOCATE(AVALUE, BVALUE, CVALUE, fxyz)
                orbnr = orbnr+PN
             ELSE
                !GENERAL CASE Thomas Kjærgaard may be slow (but only needed 
                !for g functions. TK
                angmom= J-1
                P = J*(J+1)/2
                PN = 2*angmom+1
                spSIZE=PN*P
                allocate(SPHMAT2(PN,P))
                SPHMAT2 = 0_4
                CALL BUILD_SPHMAT(angmom,spSIZE,SPHMAT2)
                ! Constructing the xyz-power of the P cartesian f-orbitals
                ! X**AVALUE * Y**BVALUE * Z*CVALUE
                MMM = J*(J+1)/2*J*(J+1)/2
                ALLOCATE(AVALUE(MMM), BVALUE(MMM), CVALUE(MMM))
                P2=0
                Rloop24: DO R = 1,P
                   DO Q = 1,R
                      P2=P2+1
                      AVALUE(P2)=J-R
                      BVALUE(P2)=R-Q
                      CVALUE(P2)=Q-1
                      IF(P2.EQ.P)EXIT Rloop24
                   ENDDO
                ENDDO Rloop24
                ALLOCATE(fxyz(P))
                DO Q = 1, P
                   fxyz(Q) = (X**AVALUE(Q))*(Y**BVALUE(Q))*(Z**CVALUE(Q))*GAO(orbnr)
                ENDDO
                DO R = 1, PN                                      
                   Q = 1
                   GAO(orbnr+R-1) = SPHMAT2(R,Q)*fxyz(Q)
                   DO Q = 2, P
                      GAO(orbnr+R-1) = GAO(orbnr+R-1) + SPHMAT2(R,Q)*fxyz(Q)
                   ENDDO
                END DO
                DEALLOCATE(AVALUE, BVALUE, CVALUE, fxyz)
                DEALLOCATE(SPHMAT2)
                orbnr = orbnr+PN
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  end subroutine determine_orbitals



  !> \brief build spherical transformation matrices
  !> \author T. Kjaergaard
  !> \date 2010
  SUBROUTINE Build_SPHMAT(L,SIZE,SPHMAT) 
    use math_fun
    IMPLICIT NONE
    INTEGER          :: nAngmom,SIZE
    REAL(4)      :: SPHMAT(SIZE)
    !
    INTEGER          :: L,I
    Real(4), parameter :: DM1 = -1.0_4, DO = 0.0_4, D1 = 1.0_4, D2 = 2.0_4
    INTEGER  :: M1,MADR,MABS,V0, NDER, IOFF
    INTEGER  :: EE,FF,GG,BX,BY,BZ,II,JJ,KK,PX,PY,PZ
    REAL(4)  :: FACNRM,FAC1,FAC2, FAC3,FAC4, FACTOR
    INTEGER  :: T,U,V,A,B,C,P,Q,R,X,Y,Z,TOT,IADR,AX0,AY0,AZ0,NSIZE
    INTEGER  :: M,Ncol,Nrow,INDEX,NDIM,STARTINDEX
    
    IF(L .LE. 1)CALL LSQUIT('ERROR IN Build_PRECALCULATED_SPHMAT',-1)
    NSIZE=0
    NRow = 2*L+1
    NCol = (L+1)*(L+2)/2
    DO M1 = 0, 2*L 
       M = M1 - L
       IF (L.EQ. 1) THEN
          IF (M .EQ. -1) MADR =  0  
          IF (M .EQ.  0) MADR =  1 
          IF (M .EQ.  1) MADR = -1 
       ELSE
          MADR = M
       END IF
       MABS = ABS(M)
       V0 = 0
       IF (M .LT. 0) V0 = 1 
       FACNRM = D1
       IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(6,L+MABS)*&
            &FACULT(6,L-MABS))/(FACULT(6,L)*(D2**MABS))
       FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
       FACNRM = FACNRM/SQRT(FACUL2(6,2*L-1))
       DO T = 0, L - MABS, 2
          DO U = 0, T, 2
             DO V = V0, MABS, 2
                ! almost 6.4.48 in the book
                FAC3 = FACNRM*BINOM(6,L,T/2)*BINOM(6,L-T/2,MABS+T/2)&
                     &                    *BINOM(6,T/2,U/2)*BINOM(6,MABS,V)
                DO A = 0, MIN(0,T+MABS-U-V) 
                   DO B = 0, MIN(0,U+V)
                      DO C = 0, MIN(0,L-T-MABS)
                         ! 6.4.47 in the book
                         DO P = 0, - A, 2
                            DO Q = 0, - B, 2
                               DO R = 0, - C, 2
                                  FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
                                       &   D2**(-A-B-C-P-Q-R-T)*FAC3
                                  X = T+MABS-U-V-2*A-P
                                  Y = U+V-2*B-Q
                                  Z = L-T-MABS-2*C-R
                                  TOT = X + Y + Z
                                  IADR = 1 + (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR+NSIZE
                                  SPHMAT(IADR) = SPHMAT(IADR) + FACTOR 
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    NSIZE= NSIZE+Nrow*Ncol 
  END SUBROUTINE BUILD_SPHMAT

 


end Module grid_utilities_module
