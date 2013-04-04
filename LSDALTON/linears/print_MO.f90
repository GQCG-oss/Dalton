module print_moorb_grid_mod
  use fundamental
  use BUILDAOBATCH
  use integralinterfaceMod
  use matrix_module
  use matrix_operations
  use typedef
  use typedefTYPE
  use grid_utilities_module
contains
  subroutine print_moorb(filename,frmt,CMO,ls,iorb,jorb)
    implicit none
    type(matrix)                   :: CMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    Character*(*)                  :: frmt
    integer                        :: iorb,jorb
!
    integer                        :: I, J, nATOMS
    real(4), allocatable       :: ATOMXYZ(:,:)

    nATOMS = ls%setting%MOLECULE(1)%p%natoms
    allocate(ATOMXYZ(3,nATOMS))
    ! Loop over all atoms (natoms) and store atomic coordinates in ATOMXYZ
    DO I=1,nATOMS !Atom coordinates        
       DO J=1,3                                                         
          ATOMXYZ(J,I) = real(ls%setting%molecule(1)%p%atom(I)%CENTER(J),4)
       ENDDO
    ENDDO

    SELECT CASE (frmt)
!#if 0
!    CASE('CUBE')
!     call calculate_gcube(filename,CMO,ls,natoms,ATOMXYZ)
!#endif
    CASE('DPLT')
     call calculate_pdaltongrid(filename,iorb,CMO,ls,natoms,ATOMXYZ)
    CASE('PLT')
     call calculate_pplt(filename,iorb,CMO,ls,natoms,ATOMXYZ)
    CASE('CHPLT')
     call calculate_charge(filename,iorb,jorb,CMO,ls,natoms,ATOMXYZ)
    CASE('DEPLT')
     call calculate_density2(filename,CMO,ls,natoms,ATOMXYZ)
    CASE('EP')
     call calculate_ep3(filename,CMO,ls,natoms,ATOMXYZ)
    CASE DEFAULT
     write(ls%lupri,*) 'print_moorb() unknown format: ', frmt
    END SELECT

    deallocate(ATOMXYZ)
  end subroutine print_moorb

  subroutine print_xyz(filename,ls)
    implicit none
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer                        :: nATOMS
    integer                        :: I,J, IUNIT, ANUM
    real(4)                        :: ATOMXYZ(3)
    real(4)                        :: br
    br=real(bohr_to_angstrom,4)
      IUNIT=111
      nATOMS = ls%setting%MOLECULE(1)%p%natoms

      OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='FORMATTED')

      WRITE(IUNIT,'(I5)') nATOMS

      WRITE(IUNIT,*)

      DO I=1,nATOMS !Atom coordinates        
          ANUM=NINT(ls%setting%molecule(1)%p%atom(I)%Charge)
       DO J=1,3                                                     
          ATOMXYZ(J) = ls%setting%molecule(1)%p%atom(I)%CENTER(J)
       ENDDO

       WRITE(IUNIT,'(I3,F10.5,F10.5,F10.5)') ANUM, (ATOMXYZ(J)*br,J=1,3)
      
      ENDDO

      CLOSE(IUNIT,STATUS='KEEP')

  end subroutine print_xyz

  subroutine calculate_daltongrid(filename,nORBITALS,dCMO,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: GAO(:,:),CMO(:)
    real(4), allocatable           :: moorb(:,:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn
    integer                        :: II,TD, ijk,irec, nLeftover
    integer                        :: reclen
    real(4)                        :: br
    integer, parameter             :: nGRS = 2048
    br=real(bohr_to_angstrom,4)

    nBASIS    = dCMO%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)


    !Open daltongrid file and write header
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111

    reclen=nGRS

#ifdef GFORTRAN
    reclen=reclen*4
#endif

    OPEN (IUNIT,FILE=trim(filename),STATUS='UNKNOWN',FORM='UNFORMATTED',&
        & ACCESS='DIRECT',RECL=reclen)

    irec=1
    WRITE(IUNIT,REC=irec) nZ,nY,nX,fZ1,fZn,fY1,fYn,fX1,fXn,nORBITALS,nGRIDPOINTS,nGRS

    !allocate
    allocate(moorb(nGRS,nORBITALS), GAO(nBASIS,nGRS),CMO(nORBITALS*nBASIS)) !MO orbitals

    moorb = 0e0; ijk = 0
    CMO = real(dCMO%elms,4)

    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 
             

             
             orbnr=1
             ijk = ijk + 1

             GAO(:,ijk) = 0.0_4

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO(:,ijk),orbnr, nBASIS, X,Y,Z,ls)
             ENDDO



             !Compute and write nGRD grid points             
             IF (mod(ijk,nGRS).eq. 0) THEN

                 call SGEMM('t','n',nGRS,nORBITALS,nBASIS,1.0,GAO,nBASIS,CMO,nBASIS,0.0,moorb,nGRS)

                 DO iorb=1,nORBITALS
                  irec = irec + 1
                  WRITE(IUNIT,REC=irec) moorb(:,iorb)
                 ENDDO

                 ijk = 0
             ENDIF

          ENDDO
       ENDDO
    ENDDO


    !Write leftover gridpoints
    nLeftover = mod(nGRIDPOINTS,nGRS)
    IF (nLeftover.ne. 0) THEN

        call SGEMM('t','n',nGRS,nORBITALS,nBASIS,1.0,GAO,nBASIS,CMO,nBASIS,0.0,moorb,nGRS)

        DO iorb=1,nORBITALS
           irec = irec + 1
           WRITE(IUNIT,REC=irec) moorb(:,iorb)
        ENDDO
    ENDIF


    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb,GAO,CMO)


  end subroutine calculate_daltongrid

  subroutine calculate_pdaltongrid(filename,nORBITALS,dCMO,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,iunit,ig
    real(4)                        :: X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: GAO(:,:),CMO(:)
    real(4), allocatable           :: moorb(:,:),XgridB(:),YgridB(:),ZgridB(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,X,Y,Z
    integer                        :: II,TD, ijk,klm,irec, nLeftover
    integer                        :: reclen
    real(4)                        :: br
    integer, parameter             :: nGRS = 2048
    br=real(bohr_to_angstrom,4)

    nBASIS    = dCMO%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)


    !Open daltongrid file and write header
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111

    reclen=nGRS

#ifdef GFORTRAN
    reclen=reclen*4
#endif

    OPEN (IUNIT,FILE=trim(filename),STATUS='UNKNOWN',FORM='UNFORMATTED',&
        & ACCESS='DIRECT',RECL=reclen)

    irec=1
    WRITE(IUNIT,REC=irec) nZ,nY,nX,fZ1,fZn,fY1,fYn,fX1,fXn,nORBITALS,nGRIDPOINTS,nGRS

    !allocate
    allocate(moorb(nGRS,nORBITALS), GAO(nBASIS,nGRS),CMO(nORBITALS*nBASIS),XgridB(nGRS),YgridB(nGRS),ZgridB(nGRS)) !MO orbitals

    moorb = 0e0; ijk = 0
    CMO = real(dCMO%elms,4)

    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 
             
             
             ijk = ijk + 1

             XgridB(ijk)=Xgrid
             YgridB(ijk)=Ygrid
             ZgridB(ijk)=Zgrid
             

             !Compute and write nGRD grid points             
             IF (mod(ijk,nGRS).eq. 0) THEN


             GAO = 0.0_4

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm,I,orbnr,X,Y,Z)
                 DO klm=1,nGRS 
                  orbnr=1

                  DO I=1,natoms
                    ! X-, Y-, and Z-distances from gridpoint to atom I
                    X = XgridB(klm) - ATOMXYZ(1,I)
                    Y = YgridB(klm) - ATOMXYZ(2,I)
                    Z = ZgridB(klm) - ATOMXYZ(3,I)
                    call determine_orbitals(I,GAO(:,klm),orbnr, nBASIS, X,Y,Z,ls)
                  ENDDO
                 ENDDO 
!$OMP END PARALLEL DO


                 call SGEMM('t','n',nGRS,nORBITALS,nBASIS,1.0,GAO,nBASIS,CMO,nBASIS,0.0,moorb,nGRS)

                 DO iorb=1,nORBITALS
                  irec = irec + 1
                  WRITE(IUNIT,REC=irec) moorb(:,iorb)
                 ENDDO

                 ijk = 0

             ENDIF

          ENDDO
       ENDDO
    ENDDO


    !Write leftover gridpoints
    nLeftover = mod(nGRIDPOINTS,nGRS)
    IF (nLeftover.ne. 0) THEN

    GAO = 0.0_4
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm,I,orbnr,X,Y,Z)
                DO klm=1,nLeftover
                 orbnr=1

                 DO I=1,natoms
                   ! X-, Y-, and Z-distances from gridpoint to atom I
                    X = XgridB(klm) - ATOMXYZ(1,I)
                    Y = YgridB(klm) - ATOMXYZ(2,I)
                    Z = ZgridB(klm) - ATOMXYZ(3,I)
                    call determine_orbitals(I,GAO(:,klm),orbnr, nBASIS, X,Y,Z,ls)
                 ENDDO
                ENDDO 
!$OMP END PARALLEL DO


        call SGEMM('t','n',nGRS,nORBITALS,nBASIS,1.0,GAO,nBASIS,CMO,nBASIS,0.0,moorb,nGRS)

        DO iorb=1,nORBITALS
           irec = irec + 1
           WRITE(IUNIT,REC=irec) moorb(:,iorb)
        ENDDO
    ENDIF


    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb,GAO,CMO,XgridB,YgridB,ZgridB)


  end subroutine calculate_pdaltongrid


  subroutine calculate_charge(filename,iorb,jorb,dCMO,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,jorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:),CMO(:,:),GAO(:),tmp(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn, SDOT
    integer                        :: II,TD, ijk
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = dCMO%ncol
    nBASIS    = dCMO%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

    allocate(moorb(nGRIDPOINTS), GAO(nBASIS),CMO(nBASIS,nORBITALS),tmp(nBASIS*nORBITALS)) !MO orbitals
    moorb = 0.0; ijk = 0

    tmp=real(dCMO%elms,4); call SCOPY(nBASIS*nORBITALS,tmp,1,CMO,1)
    deallocate(tmp)

    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 
             

             
             GAO = 0.0_4
             orbnr=1

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO,orbnr, nBASIS, X,Y,Z,ls)
             ENDDO


             ijk = ijk + 1
             
             moorb(ijk) = SDOT(nBASIS,CMO(:,iorb),1,GAO,1) * SDOT(nBASIS,CMO(:,jorb),1,GAO,1)
             

          ENDDO
       ENDDO
    ENDDO

    deallocate(GAO,CMO)


    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
#ifdef GFORTRAN
    reclen=reclen*4
#endif

      OPEN (IUNIT,FILE=trim(filename),STATUS='NEW',FORM='UNFORMATTED',&
          & ACCESS='DIRECT',RECL=reclen)

      WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

      CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_charge


  subroutine calculate_pplt(filename,iorb,dCMO,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4),     pointer           :: GAO(:)
    real(4), allocatable           :: moorb(:),CMO(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,SDOT
    integer                        :: II,TD, ijk
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = dCMO%ncol
    nBASIS    = dCMO%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

    allocate(moorb(nGRIDPOINTS),CMO(nBASIS)) !MO orbitals
    moorb = 0; ijk = 0

    CMO=real((dCMO%elms(1+(iorb-1)*nBASIS:iorb*nBASIS)),4)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Zg,Zgrid,Yg,Ygrid,Xg,Xgrid,orbnr,I,X,Y,Z,ijk,GAO)

    allocate(GAO(nBASIS))

!$OMP DO SCHEDULE(STATIC)
    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 
             

             
             GAO = 0.0_4
             orbnr=1

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO,orbnr, nBASIS, X,Y,Z,ls)
             ENDDO

             ijk = Xg + (Yg-1)*nX + (Zg-1)*nX*nY
             
             moorb(ijk) = SDOT(nBASIS,CMO,1,GAO,1)
             

          ENDDO
       ENDDO
    ENDDO
!$OMP END DO

    deallocate(GAO)
!$OMP END PARALLEL



    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
    reclen=reclen*4

      OPEN (IUNIT,FILE=trim(filename),STATUS='NEW',FORM='UNFORMATTED',&
          & ACCESS='DIRECT',RECL=reclen)

      WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

      CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb,CMO)


  end subroutine calculate_pplt

#if 0
! subroutine calculate_ep now obsolete
  subroutine calculate_ep(filename,dD,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dD
    type(matrix)                   :: dEPint
    type(lsitem)                   :: ls
    TYPE(LSSETTING)                :: psetting
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    real(realk)                    :: R(3)
    integer                        :: I,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:),D(:),EPint(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,epnuc,SDOT
    integer                        :: II,TD, ijk, icount
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = dD%ncol
    nBASIS    = dD%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

    allocate(moorb(nGRIDPOINTS),D(nORBITALS*nBASIS)) !MO orbitals
    moorb = 0_4; ijk = 0

    D=real(dD%elms,4)


    call mem_TurnONThread_Memory()
    ls%setting%scheme%noOMP = .TRUE. 


    write(6,*) 'Computing ', nGRIDPOINTS, ' gridpoints'


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(psetting,Zg,Zgrid,Yg,Ygrid,Xg,Xgrid,epnuc,I,X,Y,Z,Q,R,dEPint,EPint,ijk)

    call init_threadmemvar()

    call mat_init(dEPint,nBASIS,nBASIS)
    allocate(EPint(nBASIS*nBASIS))

!$OMP CRITICAL
    call copy_setting(psetting,ls%setting,ls%lupri)
!$OMP END CRITICAL


!$OMP DO SCHEDULE(STATIC)
    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 
             

             epnuc = 0_4
             
             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                Q = real(ls%setting%molecule(1)%p%atom(I)%Charge,4)
                epnuc = epnuc + (Q/SQRT(X*X + Y*Y + Z*Z))
             ENDDO

             R = real((/Xgrid, Ygrid, Zgrid/),realk)
             call II_get_ep_integrals(ls%lupri,ls%luerr,psetting,dEPint,R) 
             EPint=real(dEPint%elms,4)

             ijk = Xg + (Yg-1)*nX + (Zg-1)*nX*nY
             
             moorb(ijk) = epnuc - 2.0*SDOT(nBASIS*nBASIS,D,1,EPint,1)

          ENDDO
       ENDDO
    ENDDO
!$OMP END DO

!$OMP CRITICAL
    call typedef_free_setting(psetting)
!$OMP END CRITICAL

    deallocate(EPint)
    call mat_free(dEPint)

    call collect_thread_memory()

!$OMP END PARALLEL

    call mem_TurnOffThread_Memory()
    ls%setting%scheme%noOMP = .FALSE. 


    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
#ifdef GFORTRAN
    reclen=reclen*4
#endif

      OPEN (IUNIT,FILE=trim(filename),STATUS='NEW',FORM='UNFORMATTED',&
          & ACCESS='DIRECT',RECL=reclen)

      WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

      CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb,D)


  end subroutine calculate_ep
#endif

  subroutine calculate_ep3(filename,D,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: D
    type(matrix)                   :: dEPint
    type(lsitem)                   :: ls
    TYPE(LSSETTING)                :: psetting
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    real(realk), allocatable       :: R(:,:), emoorb(:)
    integer                        :: I,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,epnuc
    integer                        :: II,TD, ijk, icount
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = D%ncol
    nBASIS    = D%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

    allocate(moorb(nGRIDPOINTS),emoorb(nGRIDPOINTS),R(3,nGRIDPOINTS)) !MO orbitals
    moorb = 0_4; emoorb = 0_realk;

    write(6,*) 'Computing ', nGRIDPOINTS, ' gridpoints'

!   Nuclear contribution

!$OMP PARALLEL DO SCHEDULE(STATIC)  DEFAULT(SHARED) PRIVATE(Zg,Zgrid,Yg,Ygrid,Xg,Xgrid,epnuc,I,X,Y,Z,Q,ijk)
    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 
             

             ijk = Xg + (Yg-1)*nX + (Zg-1)*nX*nY

             epnuc = 0_4
             
             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                Q = real(ls%setting%molecule(1)%p%atom(I)%Charge,4)
                epnuc = epnuc + (Q/SQRT(X*X + Y*Y + Z*Z))
             ENDDO

             R(:,ijk) = real((/Xgrid, Ygrid, Zgrid/),realk)
             
             moorb(ijk) = epnuc 

          ENDDO
       ENDDO
    ENDDO
!$OMP END PARALLEL DO

!   Electronic contribution

    call II_get_ep_integrals3(ls%lupri,ls%luerr,ls%setting,R,nGRIDPOINTS,D,emoorb)

!   Add the two

    deallocate(R)

    moorb = moorb - real(emoorb,4)

    deallocate(emoorb)


    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
#ifdef GFORTRAN
    reclen=reclen*4
#endif

      OPEN (IUNIT,FILE=trim(filename),STATUS='NEW',FORM='UNFORMATTED',&
          & ACCESS='DIRECT',RECL=reclen)


      WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

      CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_ep3



#if 0
!!$  subroutine calculate_gcube(filename, CMO,ls,natoms,ATOMXYZ)
!!$    !*********************************************************
!!$    ! Driver routine for calculating the molecular orbitals
!!$    ! in specified points r in space 
!!$    ! \phi(r)_{m} = \sum_{n} C_{mn} X(r)_{n}
!!$    ! where X(r)_[n] are atomic orbitals.
!!$    !*********************************************************
!!$    implicit none
!!$    type(matrix)                   :: CMO
!!$    type(lsitem)                   :: ls
!!$    Character*(*)                  :: filename
!!$    integer,intent(in)             :: natoms
!!$    real(realk), intent(in)        :: ATOMXYZ(3,nATOMS)
!!$    !Local
!!$    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
!!$    integer                        :: orbnr, nORBITALS,nBASIS,iorb,iunit,ig
!!$    real(realk)                    :: GRIDCOOR(3),charge
!!$    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,deltax, deltay, deltaz 
!!$    real(4), allocatable           :: GAO(:),moorb(:)
!!$    integer                        :: nX,nY,nZ
!!$
!!$    nORBITALS = CMO%ncol
!!$    nBASIS    = CMO%nrow
!!$    ! Determine grid-box. It extends from
!!$    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
!!$    ! and the number of gridpoints in these directions are nX/nY/nZ.
!!$    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
!!$    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
!!$         &ATOMXYZ,natoms,nGRIDPOINTS)
!!$
!!$    ! Open Cube file
!!$    IUNIT = -1
!!$    CALL LSOPEN(IUNIT,trim(filename),'UNKNOWN','FORMATTED')
!!$    ! Write cube file header
!!$    WRITE(IUNIT,*)'Dalton Cube File  '
!!$    WRITE(IUNIT,*)'Molecular Orbitals'
!!$    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')-natoms,X1,Y1,Z1
!!$    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')nX,deltax,0E0_realk,0E0_realk
!!$    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')nY,0E0_realk,deltay,0E0_realk
!!$    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')nZ,0E0_realk,0E0_realk,deltaz
!!$    do i=1,natoms
!!$       charge = ls%setting%molecule(1)%p%ATOM(i)%charge
!!$       X = ls%setting%molecule(1)%p%ATOM(i)%center(1)
!!$       Y = ls%setting%molecule(1)%p%ATOM(i)%center(2)
!!$       Z = ls%setting%molecule(1)%p%ATOM(i)%center(3)
!!$       WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6,F12.6)')INT(charge),charge,X,Y,Z
!!$    enddo
!!$    WRITE(IUNIT,'(10I5)')nORBITALS,(I,I =1,nORBITALS)   
!!$
!!$
!!$
!!$    allocate(moorb(nORBITALS), GAO(nBASIS)) !MO orbitals
!!$    moorb = 0
!!$
!!$    ! gridnr: Counter for grid points, orbnr: Counter for orbitals
!!$    !Loop over xyz-grid points to determine 3D box of grid points
!!$    DO Xg = 1,nX ! Xg: grid-coordinate number in x direction
!!$       Xgrid = X1 + deltax*(Xg-1) !x grid-coordinate
!!$       DO Yg = 1,nY ! Yg: grid-coordinate number in y direction
!!$          Ygrid = Y1 + deltay*(Yg-1) !y grid-coordinate
!!$          DO Zg = 1,nZ ! Zg: grid-coordinate number in z direction
!!$             GRIDCOOR(1) = Xgrid
!!$             GRIDCOOR(2) = Ygrid
!!$             GRIDCOOR(3) = Z1 + deltaz*(Zg-1) !z grid-coordinate
!!$             GAO = 0_4
!!$             orbnr=1
!!$
!!$             DO I=1,natoms
!!$                ! X-, Y-, and Z-distances from gridpoint to atom I
!!$                X = GRIDCOOR(1) - ATOMXYZ(1,I)
!!$                Y = GRIDCOOR(2) - ATOMXYZ(2,I)
!!$                Z = GRIDCOOR(3) - ATOMXYZ(3,I)
!!$                call determine_orbitals(I,GAO,orbnr, nBASIS, X,Y,Z,ls)
!!$             ENDDO
!!$
!!$             do iorb = 1,nORBITALS
!!$             moorb(Iorb) = 0E0_4
!!$                DO P= 1, nBASIS
!!$                   moorb(Iorb) = moorb(iorb) + CMO%elms(p+(iorb-1)*CMO%ncol)*GAO(p)
!!$                ENDDO
!!$             enddo
!!$             
!!$             write(iunit,'(6E13.5)') (moorb(Iorb),iorb=1,nORBITALS)
!!$
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$
!!$    CALL LSCLOSE(IUNIT,'KEEP')
!!$    deallocate(moorb,GAO)
!!$
!!$  end subroutine calculate_gcube
#endif

#if 0
! calculate_density() now obsolete, kept for debugging purposes
!!$  subroutine calculate_density(filename,dD,ls,natoms,ATOMXYZ)
!!$    implicit none
!!$    type(matrix)                   :: dD
!!$    type(lsitem)                   :: ls
!!$    Character*(*)                  :: filename
!!$    integer,intent(in)             :: natoms
!!$    integer                        :: nX,nY,nZ
!!$    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
!!$    !Local
!!$    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
!!$    integer                        :: orbnr,nBASIS,iorb,iunit,ig
!!$    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
!!$    real(4), allocatable           :: GAO(:),D(:),tmp(:)
!!$    real(4), allocatable           :: moorb(:)
!!$    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,SDOT
!!$    integer                        :: II,TD, ijk
!!$    integer                        :: reclen
!!$    real(4)                       :: br
!!$    br=real(bohr_to_angstrom,4)
!!$
!!$    nBASIS    = dD%nrow
!!$    ! Determine grid-box. It extends from
!!$    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
!!$    ! and the number of gridpoints in these directions are nX/nY/nZ.
!!$    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
!!$    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
!!$         &ATOMXYZ,natoms,nGRIDPOINTS)
!!$
!!$
!!$    !allocate
!!$    allocate(moorb(nGRIDPOINTS), D(nBASIS*nBASIS))
!!$
!!$    moorb = 0e0; ijk = 0
!!$    D = real(dD%elms,4)
!!$
!!$!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Zg,Zgrid,Yg,Ygrid,Xg,Xgrid,orbnr,I,X,Y,Z,ijk,GAO,tmp)
!!$
!!$   allocate(GAO(nBASIS),tmp(nBASIS))
!!$
!!$!$OMP DO SCHEDULE(STATIC)
!!$    DO Zg = 1,nZ 
!!$       Zgrid = Z1 + deltaz*(Zg-1) 
!!$       DO Yg = 1,nY 
!!$          Ygrid = Y1 + deltay*(Yg-1) 
!!$          DO Xg = 1,nX
!!$             Xgrid = X1 + deltax*(Xg-1) 
!!$             
!!$
!!$             
!!$             GAO = 0E0_4
!!$             orbnr=1
!!$
!!$             DO I=1,natoms
!!$                ! X-, Y-, and Z-distances from gridpoint to atom I
!!$                X = Xgrid - ATOMXYZ(1,I)
!!$                Y = Ygrid - ATOMXYZ(2,I)
!!$                Z = Zgrid - ATOMXYZ(3,I)
!!$                call determine_orbitals(I,GAO,orbnr, nBASIS, X,Y,Z,ls)
!!$             ENDDO
!!$
!!$             ijk = Xg + (Yg-1)*nX + (Zg-1)*nX*nY
!!$
!!$             CALL SSYMV('L',nBASIS,1.0,D,nBASIS,GAO,1,0.0,tmp,1)
!!$             moorb(ijk) = 2.0*SDOT(nBASIS,GAO,1,tmp,1)
!!$
!!$
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$!$OMP END DO
!!$
!!$    deallocate(GAO,tmp)
!!$!$OMP END PARALLEL
!!$
!!$    !Write plt file
!!$
!!$    II=3; TD=200
!!$    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
!!$    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
!!$    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);
!!$
!!$    IUNIT=111
!!$    reclen=(nGRIDPOINTS + 11)
!!$#ifdef GFORTRAN
!!$    reclen=reclen*4
!!$#endif
!!$
!!$      OPEN (IUNIT,FILE=trim(filename),STATUS='NEW',FORM='UNFORMATTED',&
!!$          & ACCESS='DIRECT',RECL=reclen)
!!$
!!$      WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb
!!$
!!$      CLOSE(IUNIT,STATUS='KEEP')
!!$
!!$    deallocate(moorb,D)
!!$
!!$
!!$  end subroutine calculate_density
#endif

  subroutine calculate_density2(filename,dD,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dD
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: GAO(:,:),D(:),tmp(:,:)
    real(4), allocatable           :: moorb(:),XgridB(:),YgridB(:),ZgridB(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,SDOT
    integer                        :: II,TD, ijk, klm, nLeftover
    integer                        :: reclen
    real(4)                        :: br
    br=real(bohr_to_angstrom,4)

    nBASIS    = dD%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

    write(*,*)  'calculate_density2 computing',nGRIDPOINTS, ' gridpoints' 

    !allocate
    allocate(moorb(nGRIDPOINTS), D(nBASIS*nBASIS))
    allocate(XgridB(nBASIS),YgridB(nBASIS),ZgridB(nBASIS))
    allocate(GAO(nBASIS,nBASIS),tmp(nBASIS,nBASIS))

    moorb = 0e0; ijk = 0; klm = 0;
    D = real(dD%elms,4)


    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 
          
   
             ijk = ijk + 1
             klm = klm + 1

             XgridB(klm)=Xgrid
             YgridB(klm)=Ygrid
             ZgridB(klm)=Zgrid
             

          IF (mod(ijk,nBASIS).eq. 0) THEN

             
             GAO = 0E0_4

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm,I,orbnr,X,Y,Z)
            DO klm=1,nBASIS
             orbnr=1

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = XgridB(klm) - ATOMXYZ(1,I)
                Y = YgridB(klm) - ATOMXYZ(2,I)
                Z = ZgridB(klm) - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO(:,klm),orbnr, nBASIS, X,Y,Z,ls)
             ENDDO
            ENDDO
!$OMP END PARALLEL DO


             CALL SSYMM('L','L',nBASIS,nBASIS,1.0,D,nBASIS,GAO,nBASIS,0.0,tmp,nBASIS)


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm)
            DO klm=1,nBASIS
             moorb(ijk-nBASIS+klm) = 2.0*SDOT(nBASIS,GAO(:,klm),1,tmp(:,klm),1)
            ENDDO
!$OMP END PARALLEL DO

          klm = 0

          ENDIF


          ENDDO
       ENDDO
    ENDDO

    !Compute leftover gridpoints
    nLeftover = mod(nGRIDPOINTS,nBASIS)
    IF (nLeftover.ne. 0) THEN

    GAO=0E0_4

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm,I,orbnr,X,Y,Z)
            DO klm=1,nLeftover
             orbnr=1

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = XgridB(klm) - ATOMXYZ(1,I)
                Y = YgridB(klm) - ATOMXYZ(2,I)
                Z = ZgridB(klm) - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO(:,klm),orbnr, nBASIS, X,Y,Z,ls)
             ENDDO
            ENDDO 
!$OMP END PARALLEL DO

             CALL SSYMM('L','L',nBASIS,nLeftover,1.0,D,nBASIS,GAO,nBASIS,0.0,tmp,nBASIS)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm)
            DO klm=1,nLeftover
             moorb(ijk-nLeftover+klm) = 2.0*SDOT(nBASIS,GAO(:,klm),1,tmp(:,klm),1)
            ENDDO
!$OMP END PARALLEL DO

    ENDIF


    deallocate(GAO,tmp,XgridB,YgridB,ZgridB,D)

    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
!removed the #ifdef VAR_GFORTRAN as this is also needed for ifort. 
    reclen=reclen*4

      OPEN (IUNIT,FILE=trim(filename),STATUS='NEW',FORM='UNFORMATTED',&
          & ACCESS='DIRECT',RECL=reclen)

      WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

      CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_density2



end module print_moorb_grid_mod


