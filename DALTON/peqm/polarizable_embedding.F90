!
!   Polarizable Embedding (PE) library
!   Copyright (C) 2013 Jógvan Magnus Haugaard Olsen
!
!   This file is part of the PE library.
!
!   The PE library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as
!   published by the Free Software Foundation, either version 3 of the
!   License, or (at your option) any later version.
!
!   The PE library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
!
!   Contact information:
!
!   Jógvan Magnus Haugaard Olsen
!   E-mail: foeroyingur@gmail.com
!
module polarizable_embedding

    use pe_precision
    use pe_blas_wrappers
    use pe_lapack_wrappers
    use pe_variables

#if defined(VAR_MPI)
    use mpi
#endif

    implicit none

    private

    ! public subroutines/functions
    public :: pe_init, pe_master, pe_dalton_input
#if defined(VAR_MPI)
    public :: pe_mpi
#endif

contains

!------------------------------------------------------------------------------

subroutine pe_init(lupri, coords, charges)

    ! Initialization routine for the PE library.
    integer :: lupri
    real(dp), dimension(:), intent(in), optional :: charges
    real(dp), dimension(:,:), intent(in), optional :: coords

    integer :: i, j, k, l
    integer :: idx, jdx, kdx, nidx
    integer, dimension(:), allocatable :: idxs
    logical, dimension(:), allocatable :: redists
    logical :: lexist
    real(dp) :: rclose, redist

    if (allocated(Rm) .and. allocated(Zm)) then
        Rm(:,:) = coords
        synced = .false.
        scfcycle = 0
        return
    end if

    luout = lupri

    if (present(coords) .and. present(charges)) then
        qmnucs = size(charges)
        allocate(Rm(3,qmnucs), Zm(1,qmnucs))
        Rm(:,:) = coords
        Zm(1,:) = charges
    else if (present(coords) .and. .not. present(charges)) then
        stop 'ERROR in pe_init: coords present but charges missing'
    else if (.not. present(coords) .and. present(charges)) then
        stop 'ERROR in pe_init: charges present but coords missing'
    end if

    ! setting up grid for MEP and CUBE calculation
    if (pe_cube) then
        origin(1) = minval(Rm(1,:)) - xsize
        origin(2) = minval(Rm(2,:)) - ysize
        origin(3) = minval(Rm(3,:)) - zsize
        step(1) = 1.0 / xgrid
        step(2) = 1.0 / ygrid
        step(3) = 1.0 / zgrid
        xsteps = int((maxval(Rm(1,:)) + xsize - origin(1)) / step(1))
        ysteps = int((maxval(Rm(2,:)) + ysize - origin(2)) / step(2))
        zsteps = int((maxval(Rm(3,:)) + zsize - origin(3)) / step(3))
        npoints = 0
        npoints = xsteps * ysteps * zsteps
        allocate(Rp(3,npoints))
        l = 1
        do i = 1, xsteps
            do j = 1, ysteps
                do k = 1, zsteps
                    Rp(1,l) = origin(1) + (i - 1) * step(1)
                    Rp(2,l) = origin(2) + (j - 1) * step(2)
                    Rp(3,l) = origin(3) + (k - 1) * step(3)
                    l = l + 1
                end do
            end do
        end do
    end if

    call read_potential(trim(potfile))

    write(luout,'(//2x,a)') 'Polarizable embedding information'
    write(luout,'(2x,a)')   '---------------------------------'
    if (nsites > 0) then
        write(luout,'(/4x,a,i6)') 'Number of classical sites: ', nsites
    end if
    if (mulorder == 5) then
        write(luout,'(/4x,a)') 'Multipole moments upto 5th order.'
    else if (mulorder == 4) then
        write(luout,'(/4x,a)') 'Multipole moments upto 4th order.'
    else if (mulorder == 3) then
        write(luout,'(/4x,a)') 'Multipole moments upto 3rd order.'
    else if (mulorder == 2) then
        write(luout,'(/4x,a)') 'Multipole moments upto 2nd order.'
    else if (mulorder == 1) then
        write(luout,'(/4x,a)') 'Multipole moments upto 1st order.'
    else if (mulorder == 0) then
        write(luout,'(/4x,a)') 'Multipole moments upto 0th order.'
    end if
    if (pe_polar) then
        if (lpol(1)) write(luout,'(/4x,a)') 'Dipole-dipole polarizabilities.'
        if (pe_gspol) then
            write(luout,'(/4x,a)') 'Dynamic response from environment will be&
                                   & neglected during response calculation.'
        end if
        if (pe_nomb) then
            write(luout,'(/4x,a)') 'Many-body interactions will be neglected.'
        end if
        if (pe_iter) then
            write(luout,'(/4x,a)') 'Iterative solver for induced moments will&
                                   & be used'
            write(luout,'(4x,a,es7.1)') 'with convergence threshold: ', thriter
        else
            write(luout,'(/4x,a)') 'Direct solver for induced moments will be&
                                  & used.'
        end if
        if (pe_damp) then
            write(luout,'(/4x,a)') 'Interactions between inducible moments&
                                   & will be damped'
            write(luout,'(4x,a,f8.4)') 'using damping coefficient:', damp
        end if
    end if
    if (pe_restart) then
         write(luout,'(/4x,a)') 'Existing files will be used to restart if&
                                & possible.'
    end if
    if (pe_cube) then
        write(luout,'(/4x,a)') 'Cube files containing the potential and&
                               & electric field from the embedding potential&
                               & will be written.'
    end if

    ! handling sites near border
    ! -----------------------------------------------
    if (pe_border) then
        ! first locate all sites within given threshold of QM nuclei
        allocate(idxs(nsites))
        idxs = 0; nidx = 0
        do i = 1, qmnucs
            do j = 1, nsites
                lexist = .false.
                do k = 1, nidx
                    if (j == idxs(k)) then
                        lexist = .true.
                        exit
                    end if
                end do
                if (lexist) cycle
                if (nrm2(Rm(:,i) - Rs(:,j)) <= Rmin) then
                    nidx = nidx + 1
                    idxs(nidx) = j
                end if
            end do
        end do

        if (border_type == 'REMOVE') then
            write(luout,*) ''
            do i = 1, nidx
                write(luout,'(4x,a,i6,2x,a)') 'Removing parameters on site:',&
                                              & idxs(i) !, elems(idxs(i))
                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0d0
                endif
                if (lpol(1)) then
                    P1s(:,idxs(i)) = 0.0d0
                end if
            end do
        else if ((border_type == 'REDIST') .or. (border_type == 'REDISA') .or.&
                &(border_type == 'REDISC')) then
            if (border_type == 'REDIST') then
                border_type = 'REDISA'
                nredist = 1
            end if
            redist = real(nredist, dp)
            allocate(redists(nsites))
            redists = .false.
            do i = 1, nidx
                rclose = 1.0d10
                do j = 1, nsites
                    lexist = .false.
                    do k = 1, nidx
                        if (j == idxs(k)) then
                            lexist = .true.
                            exit
                        end if
                    end do
                    if (lexist) cycle
                    if (nrm2(Rs(:,idxs(i)) - Rs(:,j)) <= rclose) then
                        rclose = nrm2(Rs(:,idxs(i)) - Rs(:,j))
                        idx = j
                    end if
                end do
                if (lmul(0)) then
                    M0s(:,idx) = M0s(:,idx) + M0s(:,idxs(i)) / redist
                endif
                if (border_type == 'REDISA') then
                    if (lmul(1)) then
                        M1s(:,idx) = M1s(:,idx) + M1s(:,idxs(i)) / redist
                    endif
                    if (lmul(2)) then
                        M2s(:,idx) = M2s(:,idx) + M2s(:,idxs(i)) / redist
                    endif
                    if (lmul(3)) then
                        M3s(:,idx) = M3s(:,idx) + M3s(:,idxs(i)) / redist
                    endif
                    if (lmul(4)) then
                        M4s(:,idx) = M4s(:,idx) + M4s(:,idxs(i)) / redist
                    endif
                    if (lmul(5)) then
                        M5s(:,idx) = M5s(:,idx) + M5s(:,idxs(i)) / redist
                    endif
                    if (lpol(1)) then
                        P1s(:,idx) = P1s(:,idx) + P1s(:,idxs(i)) / redist
                    end if
                end if
                if (nredist > 1) then
                    rclose = 1.0d10
                    do j = 1, nsites
                        if (j == idx) cycle
                        lexist = .false.
                        do k = 1, nidx
                            if (j == idxs(k)) then
                                lexist = .true.
                                exit
                            end if
                        end do
                        if (lexist) cycle
                        if (nrm2(Rs(:,idxs(i)) - Rs(:,j)) <= rclose) then
                            rclose = nrm2(Rs(:,idxs(i)) - Rs(:,j))
                            jdx = j
                        end if
                    end do
                    if (lmul(0)) then
                        M0s(:,jdx) = M0s(:,jdx) + M0s(:,idxs(i)) / redist
                    endif
                    if (border_type == 'REDISA') then
                        if (lmul(1)) then
                            M1s(:,jdx) = M1s(:,jdx) + M1s(:,idxs(i)) / redist
                        endif
                        if (lmul(2)) then
                            M2s(:,jdx) = M2s(:,jdx) + M2s(:,idxs(i)) / redist
                        endif
                        if (lmul(3)) then
                            M3s(:,jdx) = M3s(:,jdx) + M3s(:,idxs(i)) / redist
                        endif
                        if (lmul(4)) then
                            M4s(:,jdx) = M4s(:,jdx) + M4s(:,idxs(i)) / redist
                        endif
                        if (lmul(5)) then
                            M5s(:,jdx) = M5s(:,jdx) + M5s(:,idxs(i)) / redist
                        endif
                        if (lpol(1)) then
                            P1s(:,jdx) = P1s(:,jdx) + P1s(:,idxs(i)) / redist
                        end if
                    end if
                end if

                if (nredist > 2) then
                    rclose = 1.0d10
                    do j = 1, nsites
                        if (j == idx .or. j == jdx) cycle
                        lexist = .false.
                        do k = 1, nidx
                            if (j == idxs(k)) then
                                lexist = .true.
                                exit
                            end if
                        end do
                        if (lexist) cycle
                        if (nrm2(Rs(:,idxs(i)) - Rs(:,j)) <= rclose) then
                            rclose = nrm2(Rs(:,idxs(i)) - Rs(:,j))
                            kdx = j
                        end if
                    end do
                    if (lmul(0)) then
                        M0s(:,kdx) = M0s(:,kdx) + M0s(:,idxs(i)) / redist
                    endif
                    if (border_type == 'REDISA') then
                        if (lmul(1)) then
                            M1s(:,kdx) = M1s(:,kdx) + M1s(:,idxs(i)) / redist
                        endif
                        if (lmul(2)) then
                            M2s(:,kdx) = M2s(:,kdx) + M2s(:,idxs(i)) / redist
                        endif
                        if (lmul(3)) then
                            M3s(:,kdx) = M3s(:,kdx) + M3s(:,idxs(i)) / redist
                        endif
                        if (lmul(4)) then
                            M4s(:,kdx) = M4s(:,kdx) + M4s(:,idxs(i)) / redist
                        endif
                        if (lmul(5)) then
                            M5s(:,kdx) = M5s(:,kdx) + M5s(:,idxs(i)) / redist
                        endif
                        if (lpol(1)) then
                            P1s(:,kdx) = P1s(:,kdx) + P1s(:,idxs(i)) / redist
                        end if
                    end if
                end if

                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0d0
                endif
                if (lpol(1)) then
                    P1s(:,idxs(i)) = 0.0d0
                end if
                if (border_type == 'REDISC') then
                    write(luout,'(/4x,a,i6)') 'Redistributing charges from site:',&
                                              & idxs(i)
                    if (nredist == 3) then
                        write(luout,'(4x,a,3i6)') 'to neighbouring sites:',&
                                                  & idx, jdx, kdx
                    else if (nredist == 2) then
                        write(luout,'(4x,a,3i6)') 'to neighbouring sites:',&
                                                  & idx, jdx
                    else if (nredist == 1) then
                        write(luout,'(4x,a,3i6)') 'to neighbouring sites:',&
                                                  & idx
                    end if
                    write(luout,'(4x,a)') 'and removing all other parameters.'
                else if (border_type == 'REDISA') then
                    write(luout,'(/4x,a,i6)') 'Redistributing parameters from site:',&
                                              & idxs(i)
                    if (nredist == 3) then
                        write(luout,'(4x,a,3i6)') 'to neighbouring sites:',&
                                                  & idx, jdx, kdx
                    else if (nredist == 2) then
                        write(luout,'(4x,a,3i6)') 'to neighbouring sites:',&
                                                  & idx, jdx
                    else if (nredist == 1) then
                        write(luout,'(4x,a,3i6)') 'to neighbouring sites:',&
                                                  & idx
                    end if
                end if
                redists(idx) = .true.
                if (nredist > 1) redists(jdx) = .true.
                if (nredist > 2) redists(kdx) = .true.
            end do
            if (lmul(0)) then
                write(luout,'(/6x,a)') ' Resulting monopoles: '
                write(luout,'(6x,a)') '----------------------'
                do i = 1, nsites
                    if (redists(i)) then
                        write(luout,'(7x,a,1x,i6,2x,f9.4)') elems(i), i,&
                                                            & M0s(:,i)
                    end if
                end do
            end if
            if (border_type == 'REDISA') then
                if (lmul(1)) then
                    write(luout,'(/6x,a)') ' Resulting dipoles: '
                    write(luout,'(6x,a)') '--------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            write(luout,'(7x,a,1x,i6,2x,3f9.4)') elems(i), i,&
                                                                 & M1s(:,i)
                        end if
                    end do
                end if
                if (lmul(2)) then
                    write(luout,'(/6x,a)') ' Resulting quadrupoles: '
                    write(luout,'(6x,a)') '------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            write(luout,'(7x,a,1x,i6,2x,6f9.4)') elems(i), i,&
                                                                 & M2s(:,i)
                        end if
                    end do
                end if
                if (lmul(3)) then
                    write(luout,'(/6x,a)') ' Resulting octopoles: '
                    write(luout,'(6x,a)') '----------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            write(luout,'(7x,a,1x,i6,2x,10f9.4)') elems(i), i,&
                                                                  & M3s(:,i)
                        end if
                    end do
                end if
                if (lmul(4)) then
                    write(luout,'(/6x,a)') ' Resulting hexadecapoles: '
                    write(luout,'(6x,a)') '--------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            write(luout,'(7x,a,1x,i6,2x,15f9.4)') elems(i), i,&
                                                                  & M4s(:,i)
                        end if
                    end do
                end if
                if (lmul(5)) then
                    write(luout,'(/6x,a)') ' Resulting ditriacontapoles: '
                    write(luout,'(6x,a)') '-----------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            write(luout,'(7x,a,1x,i6,2x,21f9.4)') elems(i), i,&
                                                                  & M5s(:,i)
                        end if
                    end do
                end if
                if (lpol(1)) then
                    write(luout,'(/6x,a)') ' Resulting polarizabilities: '
                    write(luout,'(6x,a)') '-----------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            write(luout,'(7x,a,1x,i6,2x,6f9.4)') elems(i), i,&
                                                                 & P1s(:,i)
                        end if
                    end do
                end if
            end if
            deallocate(redists)
        end if
        deallocate(idxs)
    end if

    ! number of polarizabilities different from zero
    if (lpol(1)) then
        allocate(zeroalphas(nsites))
        do i = 1, nsites
            if (abs(maxval(P1s(:,i))) <= zero) then
                zeroalphas(i) = .true.
            else
                zeroalphas(i) = .false.
                npols = npols + 1
            end if
        end do
    end if

    nullify(work)

end subroutine pe_init

!------------------------------------------------------------------------------

subroutine pe_dalton_input(word, luinp, lupri)

    character(len=7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

    integer :: i, j
    character(len=7) :: option
    character(len=2) :: auoraa

    luout = lupri

    do
        read(luinp,'(a7)') option
        call chcase(option)

        ! do a Polarizable Embedding calculation
        if (trim(option(2:)) == 'PEQM') then
            peqm = .true.
        else if (trim(option(2:)) == 'POTENT') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) potfile
            end if
        ! direct solver for induced moments
        else if (trim(option(2:)) == 'DIRECT') then
            pe_iter = .false.
        ! iterative solver for induced moments (default)
        else if (trim(option(2:)) == 'ITERAT') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) thriter
            end if
            pe_iter = .true.
        ! handling sites near quantum-classical border
         else if (trim(option(2:)) == 'BORDER') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,'(a)',advance='no') border_type
                backspace(luinp)
                call chcase(border_type)
                if ((trim(border_type) /= 'REMOVE') .and.&
                   & (trim(border_type) /= 'REDIST') .and.&
                   & (trim(border_type) /= 'REDISC') .and.&
                   & (trim(border_type) /= 'REDISA')) then
                    stop 'ERROR: unknown handling of border sites!'
                else if (trim(border_type) == 'REMOVE') then
                    read(luinp,*) border_type, Rmin, auoraa
                else if (trim(border_type) == 'REDIST') then
                    read(luinp,*) border_type, Rmin, auoraa
                else if ((trim(border_type) == 'REDISC') .or.&
                        & (trim(border_type) == 'REDISA')) then
                    read(luinp,*) border_type, Rmin, auoraa, nredist
                    if ((nredist > 3) .or. (nredist < 1)) then
                        stop 'ERROR: parameters cannot only be distributed to&
                             & minimum one site and maximum three sites.'
                    end if
                end if
                call chcase(auoraa)
                if (trim(auoraa) == 'AA') Rmin = Rmin * aa2au
            end if
            pe_border = .true.
        ! damping interactions between inducible moments
        else if (trim(option(2:)) == 'DAMP') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) damp
            end if
            pe_damp = .true.
        ! neglect dynamic response from environment
        else if (trim(option(2:)) == 'GSPOL') then
            pe_gspol = .true.
        ! neglect many-body interactions
        else if (trim(option(2:)) == 'NOMB') then
            pe_nomb = .true.
        ! Use existing files for restart
        else if (trim(option(2:)) == 'RESTAR') then
            pe_restart = .true.
        ! Write cube files
        else if (trim(option(2:)) == 'CUBE') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                do
                    read(luinp,*) option
                    call chcase(option)
                    if (trim(option(1:)) == 'COARSE') then
                        xgrid = 3
                        ygrid = 3
                        zgrid = 3
                    else if (trim(option(1:)) == 'MEDIUM') then
                        xgrid = 6
                        ygrid = 6
                        zgrid = 6
                    else if (trim(option(1:)) == 'FINE') then
                        xgrid = 12
                        ygrid = 12
                        zgrid = 12
                    else if (trim(option(1:)) == 'GRID') then
                        read(luinp,*) xsize, xgrid, ysize, ygrid, zsize, zgrid
                    else if (trim(option(1:)) == 'FIELD') then
                        cube_field = .true.
                    else if (option(1:1) == '.' .or. option(1:1) == '*') then
                        backspace(luinp)
                        exit
                    else if (option(1:1) == '!' .or. option(1:1) == '#') then
                        cycle
                    else
                        stop 'ERROR: unknown option present in .CUBE section.'
                    end if
                end do
            end if
            pe_cube = .true.
        ! verbose output
        else if (trim(option(2:)) == 'VERBOS') then
            pe_verbose = .true.
        ! debug output
        else if (trim(option(2:)) == 'DEBUG') then
            pe_debug = .true.
            pe_verbose = .true.
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        else
            write(luout,*) 'Unknown option:', option
        end if
    end do

end subroutine pe_dalton_input

!------------------------------------------------------------------------------

subroutine read_potential(filename)

    character(len=*) :: filename

    integer :: i, j, s
    integer :: nlines
    integer :: lupot
    real(dp) :: trace
    real(dp), dimension(21) :: temp
    character(len=2) :: auoraa
    character(len=80) :: word
    logical :: lexist

    inquire(file=filename, exist=lexist)
    if (lexist) then
        call openfile(filename, lupot, 'old', 'formatted')
    else
        write(luout,*) 'ERROR: potential input file not found'
        stop 'ERROR: potential input file not found'
    end if

    do
        read(lupot,*, end=100) word

        if (trim(word) == 'coordinates') then
            read(lupot,*) nsites
            read(lupot,*) auoraa
            allocate(elems(nsites), Zs(1,nsites), Rs(3,nsites))
            do i = 1, nsites
                read(lupot,*) elems(i), (Rs(j,i), j = 1, 3)
                Zs(1,i) = elem2charge(elems(i))
            end do
            call chcase(auoraa)
            if (auoraa == 'AA') then
                Rs = Rs * aa2au
            end if
        else if (trim(word) == 'monopoles') then
            lmul(0) = .true.
            if (mulorder < 0) mulorder = 0
            allocate(M0s(1,nsites))
            M0s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                M0s(1,s) = temp(1)
            end do
        else if (trim(word) == 'dipoles') then
            lmul(1) = .true.
            if (mulorder < 1) mulorder = 1
            allocate(M1s(3,nsites))
            M1s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 3)
                M1s(:,s) = temp(1:3)
            end do
        else if (trim(word) == 'quadrupoles') then
            lmul(2) = .true.
            if (mulorder < 2) mulorder = 2
            allocate(M2s(6,nsites))
            M2s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                ! remove trace
                trace = (temp(1) + temp(4) + temp(6)) / 3.0d0
                temp(1) = temp(1) - trace
                temp(4) = temp(4) - trace
                temp(6) = temp(6) - trace
                M2s(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'octopoles') then
            lmul(3) = .true.
            if (mulorder < 3) mulorder = 3
            allocate(M3s(10,nsites))
            M3s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 10)
                trace = (temp(1) + temp(4) + temp(6)) / 5.0d0
                temp(1) = temp(1) - 3.0d0 * trace
                temp(4) = temp(4) - trace
                temp(6) = temp(6) - trace
                trace = (temp(2) + temp(7) + temp(9)) / 5.0d0
                temp(2) = temp(2) - trace
                temp(7) = temp(7) - 3.0d0 * trace
                temp(9) = temp(9) - trace
                trace = (temp(3) + temp(8) + temp(10)) / 5.0d0
                temp(3) = temp(3) - trace
                temp(8) = temp(8) - trace
                temp(10) = temp(10) - 3.0d0 * trace
                M3s(:,s) = temp(1:10)
            end do
        else if (trim(word) == 'hexadecapoles') then
!            stop 'Hexadecapoles not supported currently'
            write(luout,*) 'WARNING: results will be wrong if non-traceless&
                           & hexadecapoles (16-poles) are used'
            lmul(4) = .true.
            if (mulorder < 4) mulorder = 4
            allocate(M4s(15,nsites))
            M4s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 15)
                M4s(:,s) = temp(1:15)
            end do
        else if (trim(word) == 'ditriacontapoles') then
!            stop 'Ditriacontapoles not supported currently'
            write(luout,*) 'WARNING: results will be wrong if non-traceless&
                           & ditriacontapoles (32-poles) are used'
            lmul(5) = .true.
            if (mulorder < 5) mulorder = 5
            allocate(M5s(21,nsites))
            M5s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 21)
                M5s(:,s) = temp(1:21)
            end do
        else if (trim(word) == 'isoalphas') then
            lpol(1) = .true.
            polorder = 1
            pe_polar = .true.
            if (.not. allocated(P1s)) then
                allocate(P1s(6,nsites))
                P1s = 0.0d0
            end if
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                P1s(1,s) = temp(1)
                P1s(4,s) = temp(1)
                P1s(6,s) = temp(1)
            end do
        else if (trim(word) == 'alphas') then
            lpol(1) = .true.
            polorder = 1
            pe_polar = .true.
            if (.not. allocated(P1s)) then
                allocate(P1s(6,nsites))
                P1s = 0.0d0
            end if
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                P1s(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'exclists' .or. trim(word) == 'exlists') then
            read(lupot,*) lexlst
            allocate(exclists(lexlst,nsites))
            do i = 1, nsites
                read(lupot,*) (exclists(j,i), j = 1, lexlst)
            end do
        else if (word(1:1) == '!' .or. word(1:1) == '#') then
            cycle
        end if
    end do

100 continue

    close(lupot)

   ! default exclusion list (everything polarizes everything)
    if (.not. allocated(exclists)) then
        lexlst = 1
        allocate(exclists(lexlst,nsites))
        do i = 1, nsites
            exclists(1,i) = i
        end do
    end if

end subroutine read_potential

!------------------------------------------------------------------------------

subroutine pe_master(runtype, denmats, fckmats, nmats, energies, dalwrk)

    character(*), intent(in) :: runtype
    integer, intent(in) :: nmats
    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: energies
    real(dp), dimension(:), target, intent(inout) :: dalwrk

#if defined(VAR_MPI)
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)
#else
    myid = 0
    nprocs = 1
#endif

    work => dalwrk

    ! determine what to calculate and do consistency check
    if (runtype == 'fock') then
        fock = .true.
        energy = .false.
        response = .false.
        scfcycle = scfcycle + 1
        if (.not. present(fckmats)) then
            stop 'Output matrices are missing from input'
        else if (.not. present(energies)) then
            stop 'The energy variable is missing from input'
        end if
    else if (runtype == 'energy') then
        fock = .false.
        energy = .true.
        response = .false.
    else if (runtype == 'response') then
        if (pe_gspol) return
        fock = .false.
        energy = .false.
        response = .true.
        if (.not. present(fckmats)) then
            stop 'Output matrices are missing from input'
        end if
    else
        stop 'Could not determine calculation type.'
    end if

    ndens = nmats
    nnbas = size(denmats) / ndens
    nbas = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * nnbas) - 1.0d0))

    if (.not. allocated(Epe)) allocate(Epe(ndens))
    Epe = 0.0d0
    if (.not. allocated(Ees)) allocate(Ees(0:5,ndens))
    Ees = 0.0d0
    if (.not. allocated(Epol)) allocate(Epol(3,ndens))
    Epol = 0.0d0

    if (nprocs == 1) then
        site_start = 1
        site_finish = nsites
        cube_start = 1
        cube_finish = npoints
    end if

#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        call mpi_bcast(44, 1, impi, 0, comm, ierr)
        if (fock) then
            call mpi_bcast(1, 1, impi, 0, comm, ierr)
        else if (energy) then
            call mpi_bcast(2, 1, impi, 0, comm, ierr)
        else if (response) then
            call mpi_bcast(3, 1, impi, 0, comm, ierr)
        end if

        call mpi_bcast(nbas, 1, impi, 0, comm, ierr)
        call mpi_bcast(ndens, 1, impi, 0, comm, ierr)
        call mpi_bcast(nnbas, 1, impi, 0, comm, ierr)
        call mpi_bcast(denmats, nnbas * ndens, rmpi, 0, comm, ierr)

        call mpi_bcast(scfcycle, 1, impi, 0, comm, ierr)
        call mpi_bcast(synced, 1, lmpi, 0, comm, ierr)
        if (.not. synced) then
            call pe_sync()
        end if
    end if
#endif

    if (fock) then
        call pe_fock(denmats, fckmats, energies)
    else if (energy) then
        call pe_fock(denmats)
        call pe_print_info()
        if (pe_cube) then
            if (ndens > 1) stop 'ERROR: CUBE not implemented for more than 1 density matrix'
            call compute_cube(denmats)
        end if
    else if (response) then
        call pe_fock(denmats, fckmats)
    end if

    nullify(work)

end subroutine pe_master

!------------------------------------------------------------------------------

#if defined(VAR_MPI)
subroutine pe_mpi(dalwrk, runtype)

    real(dp), dimension(:), target, intent(inout) :: dalwrk
    integer :: runtype

    integer :: i, nwrk

    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)

    work => dalwrk

    nwrk = size(work)

    if (runtype == 1) then
        fock = .true.
        energy = .false.
        response = .false.
    else if (runtype == 2) then
        fock = .false.
        energy = .true.
        response = .false.
    else if (runtype == 3) then
        fock = .false.
        energy = .false.
        response = .true.
    end if

    call mpi_bcast(nbas, 1, impi, 0, comm, ierr)
    call mpi_bcast(ndens, 1, impi, 0, comm, ierr)
    call mpi_bcast(nnbas, 1, impi, 0, comm, ierr)
    call mpi_bcast(work(1), nnbas * ndens, rmpi, 0, comm, ierr)

    if (.not. allocated(Epe)) allocate(Epe(ndens))
    Epe = 0.0d0
    if (.not. allocated(Ees)) allocate(Ees(0:5,ndens))
    Ees = 0.0d0
    if (.not. allocated(Epol)) allocate(Epol(3,ndens))
    Epol = 0.0d0

    call mpi_bcast(scfcycle, 1, impi, 0, comm, ierr)
    call mpi_bcast(synced, 1, lmpi, 0, comm, ierr)
    if (.not. synced) then
        call pe_sync()
    end if

    if (fock) then
        call pe_fock(work(1:ndens*nnbas), work(ndens*nnbas+1:2*ndens*nnbas),&
                    & work(2*ndens*nnbas+1:2*ndens*nnbas+ndens))
    else if (energy) then
        call pe_fock(work(1:ndens*nnbas))
        if (pe_cube) then
            call compute_cube(work(1:ndens*nnbas))
        end if
    else if (response) then
        call pe_fock(work(1:ndens*nnbas), work(ndens*nnbas+1:2*ndens*nnbas))
    end if

    nullify(work)

end subroutine pe_mpi

!------------------------------------------------------------------------------

subroutine pe_sync()

    integer :: i, j, k
    integer :: quotient, remainder

    if (.not. allocated(siteloops)) allocate(siteloops(0:nprocs))
    if (.not. allocated(sitedists)) allocate(sitedists(0:nprocs-1))
    if (myid == 0) then
        quotient = nsites / nprocs
        sitedists = quotient
        if (nprocs * quotient < nsites) then
            remainder = nsites - nprocs * quotient
            do i = 1, remainder
                sitedists(i-1) = sitedists(i-1) + 1
            end do
        end if
        siteloops(0) = 0
        do i = 1, nprocs
            siteloops(i) = sum(sitedists(0:i-1))
        end do
        if (.not. allocated(displs)) allocate(displs(0:nprocs))
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + 3 * sitedists(i-1)
        end do
    end if

    call mpi_bcast(nsites, 1, impi, 0, comm, ierr)
    call mpi_bcast(sitedists, nprocs, impi, 0, comm, ierr)
    call mpi_bcast(siteloops, nprocs + 1, impi, 0, comm, ierr)

    site_start = siteloops(myid) + 1
    site_finish = siteloops(myid+1)

    call mpi_bcast(qmnucs, 1, impi, 0, comm, ierr)

    if (myid /= 0 .and. .not. allocated(Zm)) allocate(Zm(1,qmnucs))
    call mpi_bcast(Zm, qmnucs, rmpi, 0, comm, ierr)

    if (myid /= 0 .and. .not. allocated(Rm)) allocate(Rm(3,qmnucs))
    call mpi_bcast(Rm, 3 * qmnucs, rmpi, 0, comm, ierr)

    if (myid /= 0 .and. .not. allocated(Rs)) allocate(Rs(3,nsites))
    call mpi_bcast(Rs, 3 * nsites, rmpi, 0, comm, ierr)

    call mpi_bcast(pe_polar, 1, lmpi, 0, comm, ierr)
    call mpi_bcast(lpol, 1, lmpi, 0, comm, ierr)

    if (pe_polar) then
        if (lpol(1)) then
            if (.not. allocated(poldists)) allocate(poldists(0:nprocs-1))
            if (myid == 0) then
                poldists = 0
                do i = 1, nprocs
                    do j = siteloops(i-1) + 1, siteloops(i)
                        if (zeroalphas(j)) then
                            continue
                        else
                            poldists(i-1) = poldists(i-1) + 1
                        end if
                    end do
                end do
            end if
            call mpi_bcast(poldists, nprocs, impi, 0, comm, ierr)
            call mpi_bcast(npols, 1, impi, 0, comm, ierr)
            call mpi_bcast(lexlst, 1, impi, 0, comm, ierr)
            if (myid /= 0 .and. .not. allocated(exclists)) then
                allocate(exclists(lexlst,nsites))
            end if
            call mpi_bcast(exclists, lexlst * nsites, impi, 0, comm, ierr)
            if (myid /= 0 .and. .not. allocated(zeroalphas)) then
                allocate(zeroalphas(nsites))
            end if
            call mpi_bcast(zeroalphas, nsites, lmpi, 0, comm, ierr)
        end if
        call mpi_bcast(pe_iter, 1, lmpi, 0, comm, ierr)
        if (pe_iter) then
            if (myid /= 0 .and. .not. allocated(P1s)) allocate(P1s(6,nsites))
            call mpi_bcast(P1s, 6 * nsites, rmpi, 0, comm, ierr)
        end if
        call mpi_bcast(pe_nomb, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(pe_damp, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(damp, 1, rmpi, 0, comm, ierr)
    end if

    call mpi_bcast(pe_restart, 1, lmpi, 0, comm, ierr)

    call mpi_bcast(mulorder, 1, impi, 0, comm, ierr)
    call mpi_bcast(lmul, 6, lmpi, 0, comm, ierr)

    if (lmul(0)) then
        if (myid /= 0 .and. .not. allocated(M0s)) allocate(M0s(1,nsites))
        call mpi_bcast(M0s, nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(1)) then
        if (myid /= 0 .and. .not. allocated(M1s)) allocate(M1s(3,nsites))
        call mpi_bcast(M1s, 3 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(2)) then
        if (myid /= 0 .and. .not. allocated(M2s)) allocate(M2s(6,nsites))
        call mpi_bcast(M2s, 6 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(3)) then
        if (myid /= 0 .and. .not. allocated(M3s)) allocate(M3s(10,nsites))
        call mpi_bcast(M3s, 10 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(4)) then
        if (myid /= 0 .and. .not. allocated(M4s)) allocate(M4s(15,nsites))
        call mpi_bcast(M4s, 15 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(5)) then
        if (myid /= 0 .and. .not. allocated(M5s)) allocate(M5s(21,nsites))
        call mpi_bcast(M5s, 21 * nsites, rmpi, 0, comm, ierr)
    end if

    call mpi_bcast(pe_cube, 1, lmpi, 0, comm, ierr)

    if (pe_cube) then
        call mpi_bcast(cube_field, 1, lmpi, 0, comm, ierr)
        if (.not. allocated(cubeloops)) allocate(cubeloops(0:nprocs))
        if (.not. allocated(cubedists)) allocate(cubedists(0:nprocs-1))
        if (myid == 0) then
            quotient = npoints / nprocs
            cubedists = quotient
            if (nprocs * quotient < npoints) then
                remainder = npoints - nprocs * quotient
                do i = 1, remainder
                    cubedists(i-1) = cubedists(i-1) + 1
                end do
            end if
            cubeloops(0) = 0
            do i = 1, nprocs
                cubeloops(i) = sum(cubedists(0:i-1))
            end do
        end if
        call mpi_bcast(npoints, 1, impi, 0, comm, ierr)
        call mpi_bcast(cubedists, nprocs, impi, 0, comm, ierr)
        call mpi_bcast(cubeloops, nprocs + 1, impi, 0, comm, ierr)
        cube_start = cubeloops(myid) + 1
        cube_finish = cubeloops(myid+1)
        if (myid /= 0 .and. .not. allocated(Rp)) allocate(Rp(3,npoints))
        call mpi_bcast(Rp, 3 * npoints, rmpi, 0, comm, ierr)
!        if (myid == 0) then
!            displs(0) = 0
!            do i = 1, nprocs
!             displs(i) = displs(i-1) + 3 * cubedists(i-1)
!            end do
!            call mpi_scatterv(Rp, 3 * cubedists, displs, rmpi,&
!                             & mpi_in_place, 0, rmpi, 0, comm, ierr)
!        else if (myid /= 0) then
!            allocate(Rp(3,cubedists(myid)))
!            call mpi_scatterv(0, 0, 0, rmpi, Rp, 3 * cubedists(myid),&
!                             & rmpi, 0, comm, ierr)
!        end if
    end if

    synced = .true.

end subroutine pe_sync
#endif

!------------------------------------------------------------------------------

subroutine pe_print_info()

    integer :: i, j, k
    integer :: lu
    real(dp), dimension(3) :: indtot
    real(dp), dimension(3*npols,ndens) :: Mkinds
    logical :: lexist

    inquire(file='pe_induced_dipoles.bin', exist=lexist)

    if (lexist) then
        call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
        rewind(lu)
        read(lu) Mkinds
        close(lu)
    end if

    write(luout,'(/6x,a)') 'Polarizable embedding energy contributions:'
    write(luout,'(5x,a)') '---------------------------------------------'
    if (mulorder >= 0) write(luout,'(/7x,a)') 'Electrostatic contributions:'
    if (lmul(0)) write(luout,'(9x,a16,5x,f20.12)') 'Monopoles       ',&
                                                   & Ees(0,1)
    if (lmul(1)) write(luout,'(9x,a16,5x,f20.12)') 'Dipoles         ',&
                                                   & Ees(1,1)
    if (lmul(2)) write(luout,'(9x,a16,5x,f20.12)') 'Quadrupoles     ',&
                                                   & Ees(2,1)
    if (lmul(3)) write(luout,'(9x,a16,5x,f20.12)') 'Octopoles       ',&
                                                   & Ees(3,1)
    if (lmul(4)) write(luout,'(9x,a16,5x,f20.12)') 'Hexadecapoles   ',&
                                                   & Ees(4,1)
    if (lmul(5)) write(luout,'(9x,a16,5x,f20.12)') 'Ditriacontapoles',&
                                                   & Ees(5,1)
    if (lpol(1)) then
        write(luout,'(/7x,a)') 'Polarization contributions:'
        write(luout,'(9x,a16,5x,f20.12)') 'Electronic      ', Epol(1,1)
        write(luout,'(9x,a16,5x,f20.12)') 'Nuclear         ', Epol(2,1)
        if (mulorder >= 0) then
            write(luout,'(9x,a16,5x,f20.12)') 'Multipole       ', Epol(3,1)
        end if
    end if
    write(luout,'(/6x,a18,6x,f20.12)') 'Total PE energy: ', Epe(1)
    if (lpol(1) .and. pe_verbose) then
        write(luout,'(/6x,a)') 'Polarizable embedding information'
        write(luout,'(5x,a/)') '-----------------------------------'
        write(luout,'(25x,a)') 'Induced dipole moments'
        write(luout,'(9x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
        k = 1
        do j = 1, npols
            write(luout,'(7x,i6,3f15.8)') j, Mkinds(k:k+2,1)
            k = k + 3
        end do
        indtot = 0.0d0
        k = 1
        do j = 1, npols
            indtot(1) = indtot(1) + Mkinds(k,1)
            indtot(2) = indtot(2) + Mkinds(k+1,1)
            indtot(3) = indtot(3) + Mkinds(k+2,1)
        end do
        write(luout,'(/24x,a)') 'Total induced dipole moment'
        write(luout,'(23x,a,14x,a,14x,a)') 'X', 'Y', 'Z'
        write(luout,'(13x,3f15.8/)') indtot
    end if

end subroutine pe_print_info

!------------------------------------------------------------------------------

subroutine pe_fock(denmats, fckmats, energies)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: energies

    integer :: i
    logical :: es = .false.
    logical :: pol = .false.

    if ((mulorder >= 0) .and. .not. response) es = .true.
    if (pe_polar) pol = .true.

    if (fock .or. response) fckmats = 0.0d0

    if (fock) then
        if (es) call pe_electrostatic(denmats, fckmats)
        if (pol) call pe_polarization(denmats, fckmats)
    else if (energy) then
        if (es) call pe_electrostatic(denmats)
        if (pol) call pe_polarization(denmats)
    else if (response) then
        if (pol) call pe_polarization(denmats, fckmats)
    end if

    if (fock .or. energy) then
        if (myid == 0) then
            do i = 1, ndens
                Epe(i) = sum(Ees(:,i)) + sum(Epol(:,i))
            end do
            if (fock) energies = Epe
        end if
    end if

#if defined(VAR_MPI)
    if (fock .or. response) then
        if (myid == 0 .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, fckmats, ndens * nnbas, rmpi,&
                           & mpi_sum, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_reduce(fckmats, 0, ndens * nnbas, rmpi, mpi_sum, 0,&
                           & comm, ierr)
        end if
    end if
#endif

end subroutine pe_fock

!------------------------------------------------------------------------------

subroutine pe_electrostatic(denmats, fckmats)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    logical :: lexist
    integer :: lu
    integer :: i, j, k
    real(dp) :: Enuc, Etmp
    real(dp), dimension(ndens) :: Eel
    real(dp), dimension(:), allocatable :: tmpfcks

    if (myid == 0) then
        inquire(file='pe_electrostatics.bin', exist=lexist)
    end if
#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif
    if (lexist .and. fock .and. ((scfcycle > 1) .or. pe_restart)) then
        if (myid == 0) then
            call openfile('pe_electrostatics.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Etmp, fckmats
            close(lu)
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(0,i) = Ees(0,i) + dot(denmats(j:k), fckmats(j:k)) + Etmp
            end do
        end if
    else
        Etmp = 0.0d0
        if (lmul(0)) then
            if (fock) then
                call es_multipoles(M0s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M0s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(0,i) = Ees(0,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(1)) then
            if (fock) then
                call es_multipoles(M1s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M1s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(1,i) = Ees(1,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(2)) then
            if (fock) then
                call es_multipoles(M2s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M2s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(2,i) = Ees(2,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(3)) then
            if (fock) then
                call es_multipoles(M3s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M3s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(3,i) = Ees(3,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(4)) then
            if (fock) then
                call es_multipoles(M4s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M4s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(4,i) = Ees(4,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(5)) then
            if (fock) then
                call es_multipoles(M5s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M5s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(5,i) = Ees(5,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (fock) then
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                allocate(tmpfcks(ndens * nnbas))
                tmpfcks = fckmats
                call mpi_reduce(mpi_in_place, fckmats, ndens * nnbas, rmpi,&
                               & mpi_sum, 0, comm, ierr)
                call mpi_reduce(mpi_in_place, Etmp, 1, rmpi, mpi_sum, 0,&
                               & comm, ierr)
            else if (myid /= 0) then
                call mpi_reduce(fckmats, 0, ndens * nnbas, rmpi, mpi_sum, 0,&
                               & comm, ierr)
                call mpi_reduce(Etmp, 0, 1, rmpi, mpi_sum, 0, comm, ierr)
            end if
#endif
            if (myid == 0) then
                call openfile('pe_electrostatics.bin', lu, 'unknown', 'unformatted')
                rewind(lu)
                write(lu) Etmp, fckmats
                close(lu)
            end if
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                fckmats = tmpfcks
                deallocate(tmpfcks)
            end if
#endif
        end if
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Ees, 6 * ndens, rmpi, mpi_sum, 0,&
                           & comm, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Ees, 0, 6 * ndens, rmpi, mpi_sum, 0, comm, ierr)
        end if
#endif
    end if

end subroutine pe_electrostatic

!------------------------------------------------------------------------------

subroutine es_multipoles(Mks, denmats, Eel, Enuc, fckmats)

    real(dp), dimension(:,:), intent(in) :: Mks
    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc

    integer :: site, ncomps
    integer :: i, j, k, l, m, n
    real(dp) :: taylor
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(:), allocatable :: Tsm, symfacs
    real(dp), dimension(:,:), allocatable :: Mk_ints

    ncomps = size(Mks,1)

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * ncomps) - 1.0d0)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k)
    end if

    allocate(Tsm(ncomps))
    allocate(symfacs(ncomps))
    allocate(Mk_ints(nnbas,ncomps))

    Eel = 0.0d0; Enuc = 0.0d0

    i = 1
    do site = site_start, site_finish
        if (abs(maxval(Mks(:,site))) < zero) then
            i = i + 1
            cycle
        end if

        ! nuclei - multipole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,site)
            call Tk_tensor(Tsm, Rsm)
            call symmetry_factors(symfacs)
            do l = 1, ncomps
                Enuc = Enuc + taylor * symfacs(l) * Zm(1,j) * Mks(l,site) *&
                     & Tsm(l)
            end do
        end do

        ! electron - multipole interaction energy
        call Mk_integrals(Mk_ints, Rs(:,site), Mks(:,site))
        do l = 1, ndens
            m = (l - 1) * nnbas + 1
            n = l * nnbas
            Eel(l) = Eel(l) + dot(denmats(m:n), sum(Mk_ints, 2))
            if (fock) fckmats(m:n) = fckmats(m:n) + sum(Mk_ints, 2)
        end do
        i = i + 1
    end do

    deallocate(Tsm, symfacs, Mk_ints)

end subroutine es_multipoles

!------------------------------------------------------------------------------

subroutine pe_polarization(denmats, fckmats)

    external :: Tk_integrals

    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    integer :: site, ndist, nrest
    integer :: i, j, k, l, m, n
    logical :: skip
    real(dp), dimension(3) :: indtot
    real(dp), dimension(:,:), allocatable :: Fels
    real(dp), dimension(:), allocatable :: Fnucs, Fmuls
    real(dp), dimension(:,:), allocatable :: Mkinds, Fktots
    real(dp), dimension(:,:), allocatable :: Fel_ints, Vel_ints

    allocate(Mkinds(3*npols,ndens), Fktots(3*npols,ndens))
    Mkinds = 0.0d0; Fktots = 0.0d0
    if (lpol(1)) then
        allocate(Fels(3*npols,ndens))
        allocate(Fnucs(3*npols), Fmuls(3*npols))
    end if

    if (response) then
        if (lpol(1)) then
            call electron_fields(Fels, denmats)
            do i = 1, ndens
                Fktots(:3*npols,i) = Fels(:,i)
            end do
        end if
        call induced_moments(Mkinds, Fktots)
    else
        if (lpol(1)) then
            call electron_fields(Fels, denmats)
            call nuclear_fields(Fnucs)
            call multipole_fields(Fmuls)
        end if
        if (myid == 0) then
            do i = 1, ndens
                if (lpol(1)) then
                    Fktots(:3*npols,i) = Fels(:,i) + Fnucs + Fmuls
                end if
            end do
        end if
        call induced_moments(Mkinds, Fktots)
        if (myid == 0) then
            do i = 1, ndens
                if (lpol(1)) then
                    Epol(1,i) = - 0.5d0 * dot(Mkinds(:3*npols,i), Fels(:,i))
                    Epol(2,i) = - 0.5d0 * dot(Mkinds(:3*npols,i), Fnucs)
                    Epol(3,i) = - 0.5d0 * dot(Mkinds(:3*npols,i), Fmuls)
                end if
            end do
        end if
    end if
    if (myid == 0 .and. pe_debug .and. .not. energy) then
        write(luout,'(/6x,a)') 'Polarizable embedding information'
        write(luout,'(5x,a/)') '-----------------------------------'
        do i = 1, ndens
            write(luout,'(7x,a,i3)') 'Input density no.: ', i
            if (lpol(1)) then
                write(luout,'(/15x,a)') 'Total electric field at polarizable sites:'
                write(luout,'(9x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
                k = 1
                do site = 1, nsites
                    if (zeroalphas(site)) cycle
                    write(luout,'(7x,i6,3f15.8)') site, Fktots(k:k+2,i)
                    k = k + 3
                end do
                write(luout,'(/25x,a)') 'Induced dipole moments'
                write(luout,'(9x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
                k = 1
                do site = 1, nsites
                    if (zeroalphas(site)) cycle
                    write(luout,'(7x,i6,3f15.8)') site, Mkinds(k:k+2,i)
                    k = k + 3
                end do
            end if
            if (lpol(1)) then
                indtot = 0.0d0
                k = 1
                do j = 1, npols
                    indtot(1) = indtot(1) + Mkinds(k,i)
                    indtot(2) = indtot(2) + Mkinds(k+1,i)
                    indtot(3) = indtot(3) + Mkinds(k+2,i)
                    k = k + 3
                end do
                write(luout,'(/24x,a)') 'Total induced dipole moment'
                write(luout,'(23x,a,14x,a,14x,a)') 'X', 'Y', 'Z'
                write(luout,'(13x,3f15.8/)') indtot
            end if
        end do
    end if
#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        do i = 1, ndens
            if (lpol(1)) then
                displs(0) = 0
                do j = 1, nprocs
                    displs(j) = displs(j-1) + 3 * poldists(j-1)
                end do
                call mpi_scatterv(Mkinds(1:3*npols,i), 3 * poldists, displs,&
                                 & rmpi, mpi_in_place, 0, rmpi, 0, comm, ierr)
            end if
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            if (lpol(1)) then
                call mpi_scatterv(0, 0, 0, rmpi, Mkinds(1:3*npols,i),&
                                 & 3 * poldists(myid), rmpi, 0, comm, ierr)
            end if
        end do
    end if
#endif

    if (fock .or. response) then
        if (lpol(1)) then
            allocate(Fel_ints(nnbas,3))
            i = 0
            do site = site_start, site_finish
                if (zeroalphas(site)) cycle
                call Tk_integrals('es', Fel_ints, nnbas, 3, Rs(:,site))
                do j = 1, 3
                    do k = 1, ndens
                        l = (k - 1) * nnbas + 1
                        m = k * nnbas
                        fckmats(l:m) = fckmats(l:m) - Mkinds(i+j,k) *&
                                     & Fel_ints(:,j)
                    end do
                end do
                i = i + 3
            end do
            deallocate(Fel_ints)
        end if
    end if

    deallocate(Mkinds, Fktots)
    if (lpol(1)) then
        deallocate(Fels)
        deallocate(Fnucs, Fmuls)
    end if

end subroutine pe_polarization

!------------------------------------------------------------------------------

subroutine induced_moments(Mkinds, Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: i, j, k
    real(dp) :: M1_size, M1ind_size

    if (pe_iter) then
        call iterative_solver(Mkinds, Fs)
    else
        if (myid == 0) then
            call direct_solver(Mkinds, Fs)
        end if
    end if

    ! check induced dipoles
    if (myid == 0 .and. .not. response) then
        do i = 1, ndens
            k = 1
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                M1ind_size = nrm2(Mkinds(k:k+2,i))
                if (mulorder >= 1) M1_size = nrm2(M1s(:,j))
                if ((mulorder >= 1) .and. (M1_size > zero)) then
                    if (M1ind_size > M1_size) then
                        write(luout,10) 'WARNING: induced dipole is larger&
                                        & than the permanent dipole at&
                                        & site: ', j
                        write(luout,11) 'induced dipole magnitude:', M1ind_size
                        write(luout,11) 'permanent dipole magnitude:', M1_size
                    end if
                else if (M1ind_size > 1.0d0) then
                    write(luout,10) 'WARNING: large induced dipole at site:', j
                    write(luout,11) 'induced dipole magnitude:', M1ind_size
                end if
                k = k + 3
            end do
        end do
    end if

 10 format(/a,i6)
 11 format(9x,a,f10.4)

end subroutine induced_moments

!------------------------------------------------------------------------------

subroutine iterative_solver(Mkinds, Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: lu, iter
    integer :: i, j, k, l, m, n, o, p, q
    logical :: exclude, lexist
    logical :: converged = .false.
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: R, R3, R5, Rd, ai, aj, norm
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp), dimension(:), allocatable :: T, Rij, Ftmp, M1tmp

    if (myid == 0) then
        inquire(file='pe_induced_dipoles.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist .and. (fock .or. energy) .and. .not. pe_nomb) then
        if (myid == 0) then
            call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if
    end if

    allocate(T(6), Rij(3), Ftmp(3), M1tmp(3))
    do n = 1, ndens
        if (.not. lexist .or. response .or. pe_nomb) then
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * poldists(i-1)
                end do
                call mpi_scatterv(Fs(:,n), 3 * poldists, displs, rmpi,&
                                 &mpi_in_place, 0, rmpi, 0, comm, ierr)
            else if (myid /= 0) then
                call mpi_scatterv(0, 0, 0, rmpi, Fs(:,n), 3 * poldists(myid),&
                                 & rmpi, 0, comm, ierr)
            end if
#endif

            l = 1
            do i = site_start, site_finish
                if (zeroalphas(i)) cycle
                call spmv(P1s(:,i), Fs(l:l+2,n), Mkinds(l:l+2,n), 'L')
                l = l + 3
            end do

#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * poldists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, 0, rmpi, Mkinds(:,n),&
                                & 3 * poldists, displs, rmpi, 0, comm, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Mkinds(:,n), 3 * poldists(myid), rmpi, 0, 0,&
                                & 0, rmpi, 0, comm, ierr)
            end if
#endif
        end if

        if (pe_nomb) cycle

#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + 3 * poldists(i-1)
            end do
            call mpi_scatterv(Mkinds(:,n), 3 * poldists, displs, rmpi,&
                             & mpi_in_place, 0, rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_scatterv(0, 0, 0, rmpi, Mkinds(:,n), 3 * poldists(myid),&
                             & rmpi, 0, comm, ierr)
        end if
#endif

        iter = 1
        do
            norm = 0.0d0
            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                if (pe_damp) then
                    ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
                end if
                m = 1
                Ftmp = 0.0d0
                do j = site_start, site_finish
                    if (zeroalphas(j)) cycle
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (i == j .or. exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    R = nrm2(Rij)
                    R3 = R**3
                    R5 = R**5
                    ! damping parameters
                    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                    ! a = 2.1304 = damp
                    ! u = R / (alpha_i * alpha_j)**(1/6)
                    ! fe = 1-(a²u²/2+au+1)*exp(-au)
                    ! ft = 1-(a³u³/6+a²u²/2+au+1)*exp(-au)
                    if (pe_damp) then
                        aj = (P1s(1,j) + P1s(4,j) + P1s(6,j)) * d3i
                        Rd = damp * R / (ai * aj)**d6i
                        fe = 1.0d0 - (0.5d0 * Rd**2 + Rd + 1.0d0) * exp(-Rd)
                        ft = fe - d6i * Rd**3 * exp(-Rd)
                    end if
                    q = 1
                    do o = 1, 3
                        do p = o, 3
                            T(q) = 3.0d0 * Rij(o) * Rij(p) * ft / R5
                            if (o == p) then
                                T(q) = T(q) - fe / R3
                            end if
                            q = q + 1
                        end do
                    end do
                    call spmv(T, Mkinds(m:m+2,n), Ftmp, 'L', 1.0d0, 1.0d0)
                    m = m + 3
                end do

#if defined(VAR_MPI)
                if (myid == 0 .and. nprocs > 1) then
                    call mpi_reduce(mpi_in_place, Ftmp, 3, rmpi, mpi_sum, 0,&
                                   & comm, ierr)
                else if (myid /= 0) then
                    call mpi_reduce(Ftmp, 0, 3, rmpi, mpi_sum, 0, comm, ierr)
                end if
#endif

                if (myid == 0) then
                    M1tmp = Mkinds(l:l+2,n)
                    Ftmp = Ftmp + Fs(l:l+2,n)
                    call spmv(P1s(:,i), Ftmp, Mkinds(l:l+2,n), 'L')
                    norm = norm + nrm2(Mkinds(l:l+2,n) - M1tmp)
                end if

#if defined(VAR_MPI)
                if (myid == 0 .and. nprocs > 1) then
                    displs(0) = 0
                    do j = 1, nprocs
                        displs(j) = displs(j-1) + 3 * poldists(j-1)
                    end do
                    call mpi_scatterv(Mkinds(:,n), 3*poldists, displs, rmpi,&
                                     & mpi_in_place, 0, rmpi, 0, comm, ierr)
                else if (myid /= 0) then
                    call mpi_scatterv(0, 0, 0, rmpi, Mkinds(:,n),&
                                     & 3*poldists(myid), rmpi, 0, comm, ierr)
                end if
#endif
                l = l + 3
            end do

            if (myid == 0) then
                if (norm < thriter) then
                    if (pe_verbose .and. .not. energy) then
                        write (luout,'(6x,a,i2,a)') 'Induced dipole moments&
                                                     & converged in ', iter,&
                                                     & ' iterations.'
                    end if
                    converged = .true.
                else if (iter > 50) then
                    write(luout,*) 'ERROR: could not converge induced dipole&
                                   & moments.'
                    stop 'ERROR: could not converge induced dipole moments.'
                else
                    converged = .false.
                    iter = iter + 1
                end if
            end if

#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_bcast(converged, 1, lmpi, 0, comm, ierr)
            end if
#endif
            if (converged) exit
        end do
    end do

    if (fock) then
        if (myid == 0) then
            call openfile('pe_induced_dipoles.bin', lu, 'unknown',&
                         & 'unformatted')
            rewind(lu)
            write(lu) Mkinds
            close(lu)
        end if
    end if

    deallocate(T, Rij, Ftmp, M1tmp)

end subroutine iterative_solver

!------------------------------------------------------------------------------

subroutine direct_solver(Mkinds, Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    logical :: lexist
    integer :: lu, info, i
    integer, dimension(:), allocatable :: ipiv
    real(dp) :: anorm, rcond
    real(dp), dimension(:), allocatable :: B

    allocate(B((3*npols)*(3*npols+1)/2))

    inquire(file='pe_response_matrix.bin', exist=lexist)
    if (lexist .and. ((scfcycle > 1) .or. pe_restart)) then
        call openfile('pe_response_matrix.bin', lu, 'old', 'unformatted')
        rewind(lu)
        if (chol) then
            read(lu) B
        else
            allocate(ipiv(3*npols))
            read(lu) B, ipiv
        end if
        close(lu)
    else
        call response_matrix(B)
        if (pe_debug) then
            anorm = lansp('1', B, 'L')
            write(luout,'(/4x,a,f15.8)') '1-norm of response matrix B: ', anorm
        end if

        if (chol) then
            call pptrf(B, 'L', info)
            if (info /= 0) then
                print *, 'Cholesky factorization failed. Trying regular...'
                allocate(ipiv(3*npols))
                call sptrf(B, 'L', ipiv, info)
                if (info /= 0) then
                    stop 'ERROR: cannot create classical response matrix.'
                else
                    chol = .false.
                end if
            end if
        else
            allocate(ipiv(3*npols))
            call sptrf(B, 'L', ipiv, info)
            if (info /= 0) then
                stop 'ERROR: cannot create classical response matrix.'
            end if
        end if

        if (pe_debug) then
            if (chol) then
                call ppcon(B, anorm, rcond, 'L')
                write(luout,'(4x,a,f15.8/)') 'Condition number of response&
                                             & matrix B: ', 1.0d0 / rcond
            else
                call spcon(B, ipiv, anorm, rcond, 'L')
                write(luout,'(4x,a,f15.8/)') 'Condition number of response&
                                             & matrix B: ', 1.0d0 / rcond
            end if
        end if

        call openfile('pe_response_matrix.bin', lu, 'unknown', 'unformatted')
        rewind(lu)
        if (chol) then
            write(lu) B
        else
            write(lu) B, ipiv
        end if
        close(lu)

    end if

    Mkinds = Fs
    if (chol) then
        call pptrs(B, Mkinds, 'L', info)
        if ((info /= 0) .and. (lexist .and. ((scfcycle > 1) .or. pe_restart))) then
            call openfile('pe_response_matrix.bin', lu, 'old', 'unformatted')
            rewind(lu)
            allocate(ipiv(3*npols))
            read(lu) B, ipiv
            close(lu)
            call sptrs(B, Mkinds, ipiv, 'L', info)
            if (info /= 0) then
                stop 'ERROR: cannot create classical response matrix.'
            end if
            deallocate(ipiv)
            chol = .false.
        else if ((info /= 0) .and. .not. (lexist .and. ((scfcycle > 1) .or. pe_restart))) then
            stop 'ERROR: cannot create classical response matrix.'
        end if
        deallocate(B)
    else
        call sptrs(B, Mkinds, ipiv, 'L', info)
        if (info /= 0) then
            stop 'ERROR: cannot create classical response matrix.'
        end if
        deallocate(B, ipiv)
    end if

    if (fock) then
        call openfile('pe_induced_dipoles.bin', lu, 'unknown', 'unformatted')
        rewind(lu)
        write(lu) Mkinds
        close(lu)
    end if

end subroutine direct_solver

!------------------------------------------------------------------------------

subroutine electron_fields(Fels, denmats)

    external :: Tk_integrals

    real(dp), dimension(:,:), intent(out) :: Fels
    real(dp), dimension(:), intent(in) :: denmats

    logical :: skip
    integer :: site
    integer :: i, j, k, l, m
    real(dp), dimension(:,:), allocatable :: Fel_ints

    allocate(Fel_ints(nnbas,3))

    Fels = 0.0d0

    i = 0
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        call Tk_integrals('es', Fel_ints, nnbas, 3, Rs(:,site))
        do j = 1, 3
            do k = 1, ndens
                l = (k - 1) * nnbas + 1
                m = k * nnbas
                Fels(i+j,k) = dot(denmats(l:m), Fel_ints(:,j))
            end do
        end do
        i = i + 3
    end do

#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + 3 * poldists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(mpi_in_place, 0, rmpi, Fels(:,i), 3 * poldists,&
                            & displs, rmpi, 0, comm, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_gatherv(Fels(:,i), 3 * poldists(myid), rmpi, 0, 0, 0,&
                            & rmpi, 0, comm, ierr)
        end do
    end if
#endif

    deallocate(Fel_ints)

end subroutine electron_fields

!------------------------------------------------------------------------------

subroutine nuclear_fields(Fnucs)

    real(dp), dimension(:), intent(out) :: Fnucs

    logical :: lexist, skip
    integer :: lu, site
    integer :: i, j, k
    real(dp), dimension(3) :: Rms, Tms

    if (myid == 0) then
        inquire(file='pe_nuclear_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist .and. ((scfcycle > 1) .or. pe_restart)) then
        if (myid == 0) then
            call openfile('pe_nuclear_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fnucs
            close(lu)
        end if
    else
        Fnucs = 0.0d0
        i = 0
        do site = site_start, site_finish
            if (zeroalphas(site)) cycle
            do j = 1, qmnucs
                Rms = Rs(:,site) - Rm(:,j)
                call Tk_tensor(Tms, Rms)
                do k = 1, 3
                    Fnucs(i+k) = Fnucs(i+k) - Zm(1,j) * Tms(k)
                end do
            end do
            i = i + 3
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + 3 * poldists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, 0, rmpi, Fnucs, 3*poldists, displs,&
                            & rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Fnucs, 3*poldists(myid), rmpi, 0, 0, 0, rmpi, 0,&
                            & comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('pe_nuclear_field.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Fnucs
            close(lu)
        end if
    end if

end subroutine nuclear_fields

!------------------------------------------------------------------------------

subroutine multipole_fields(Fmuls)

    real(dp), dimension(:), intent(out) :: Fmuls

    logical :: exclude, lexist
    integer :: lu
    integer :: i, j, k, l
    real(dp), dimension(3) :: Rji

    if (myid == 0) then
        inquire(file='pe_multipole_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif
    if (lexist .and. ((scfcycle > 1) .or. pe_restart)) then
        if (myid == 0) then
            call openfile('pe_multipole_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fmuls
            close(lu)
        end if
    else
        Fmuls = 0.0d0
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            do j = site_start, site_finish
                if (i == j) then
                    cycle
                end if
                exclude = .false.
                do k = 1, lexlst
                    if (exclists(k,i) == exclists(1,j)) then
                        exclude = .true.
                        exit
                    end if
                end do
                if (exclude) cycle
                Rji = Rs(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (abs(maxval(M0s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M0s(:,j))
                    end if
                end if
                if (lmul(1)) then
                    if (abs(maxval(M1s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M1s(:,j))
                    end if
                end if
                if (lmul(2)) then
                    if (abs(maxval(M2s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M2s(:,j))
                    end if
                end if
                if (lmul(3)) then
                    if (abs(maxval(M3s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M3s(:,j))
                    end if
                end if
                if (lmul(4)) then
                    if (abs(maxval(M4s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M4s(:,j))
                    end if
                end if
                if (lmul(5)) then
                    if (abs(maxval(M5s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M5s(:,j))
                    end if
                end if
            end do
            l = l + 3
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Fmuls, 3*npols, rmpi, mpi_sum, 0,&
                           & comm, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Fmuls, 0, 3*npols, rmpi, mpi_sum, 0, comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('pe_multipole_field.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Fmuls
            close(lu)
        end if
     end if

end subroutine multipole_fields

!------------------------------------------------------------------------------

subroutine multipole_potential(Vi, Rji, Mkj)

    real(dp), intent(inout) :: Vi
    real(dp), dimension(3), intent(in) :: Rji
    real(dp), dimension(:), intent(in) :: Mkj

    integer :: k
    integer :: a, x, y, z
    real(dp) :: taylor

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Mkj)) - 1.0d0)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k)
    end if

    a = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                Vi = Vi + taylor * symfac(x,y,z) * T(Rji,x,y,z) * Mkj(a)
                a = a + 1
            end do
        end do
     end do

end subroutine multipole_potential

!------------------------------------------------------------------------------

subroutine multipole_field(Fi, Rji, Mkj)

    real(dp), dimension(3), intent(inout) :: Fi
    real(dp), dimension(3), intent(in) :: Rji
    real(dp), dimension(:), intent(in) :: Mkj

    integer :: k
    integer :: a, b, c, x, y, z
    real(dp) :: taylor

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Mkj)) - 1.0d0))

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k-1)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k-1)
    end if

    a = 1; b = 1; c = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                if (x /= 0) then
                    Fi(1) = Fi(1) + taylor * symfac(x-1,y,z) * T(Rji,x,y,z) *&
                          & Mkj(a)
                    a = a + 1
                end if
                if (y /= 0) then
                    Fi(2) = Fi(2) + taylor * symfac(x,y-1,z) * T(Rji,x,y,z) *&
                          & Mkj(b)
                    b = b + 1
                end if
                if (z /= 0) then
                    Fi(3) = Fi(3) + taylor * symfac(x,y,z-1) * T(Rji,x,y,z) *&
                          & Mkj(c)
                    c = c + 1
                end if
            end do
        end do
     end do

end subroutine multipole_field

!------------------------------------------------------------------------------

subroutine response_matrix(B)

    real(dp), dimension(:), intent(out) :: B

    logical :: exclude
    integer :: info
    integer :: i, j, k, l, m, n, o
    integer, dimension(3) :: ipiv
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: Rd, ai, aj
    real(dp) :: R, R3, R5, T
    real(dp), dimension(3) :: Rij
    real(dp), dimension(6) :: P1inv

    B = 0.0d0

    m = 0
    do i = 1, nsites
        if (zeroalphas(i)) cycle
        P1inv = P1s(:,i)
        call pptrf(P1inv, 'L', info)
        if (info /= 0) then
            P1inv = P1s(:,i)
            call sptrf(P1inv, 'L', ipiv, info)
            if (info /= 0) then
                stop 'ERROR: could not factorize polarizability.'
            else if (chol) then
                chol = .false.
            end if
            call sptri(P1inv, ipiv, 'L')
        else
            call pptri(P1inv, 'L')
        end if
        if (pe_damp) then
            ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
        end if
        do l = 3, 1, -1
            do j = i, nsites
                if (zeroalphas(j)) cycle
                if (j == i) then
                    if (l == 3) then
                        do k = 1, l
                            B(m+k) = P1inv(k)
                        end do
                    else if (l == 2) then
                        do k = 1, l
                            B(m+k) = P1inv(3+k)
                        end do
                    else if (l == 1) then
                        do k = 1, l
                            B(m+k) = P1inv(5+k)
                        end do
                    end if
                    m = m + l
                else
                    if (pe_nomb) then
                        m = m + 3
                        cycle
                    end if
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    R = nrm2(Rij)
                    R3 = R**3
                    R5 = R**5
                    ! damping parameters
                    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                    ! a = 2.1304 = damp
                    ! u = R / (alpha_i * alpha_j)**(1/6)
                    ! fe = 1-(a²u²/2+au+1)*exp(-au)
                    ! ft = 1-(a³u³/6+a²u²/2+au+1)*exp(-au)
                    if (pe_damp) then
                        aj = (P1s(1,j) + P1s(4,j) + P1s(6,j)) * d3i
                        Rd = damp * R / (ai * aj)**d6i
                        fe = 1.0d0 - (0.5d0 * Rd**2 + Rd + 1.0d0) * exp(-Rd)
                        ft = fe - d6i * Rd**3 * exp(-Rd)
                    end if
                    if (l == 3) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(1) * Rij(k) * ft / R5
                            if (k == 1) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    else if (l == 2) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(2) * Rij(k) * ft / R5
                            if (k == 2) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    else if (l == 1) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(3) * Rij(k) * ft / R5
                            if (k == 3) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    end if
                    m = m + 3
                end if ! i /= j
            end do  ! do j = i, nsites
        end do  ! do l = 3, 1, -1
    end do  ! do i = 1, nsites

    if (pe_debug) then
        do i = 1, (3 * npols) * (3 * npols + 1) / 2
            write (luout,*) 'Response matrix(i)',i, B(i)
        end do
    end if

end subroutine response_matrix

!------------------------------------------------------------------------------

subroutine Tk_coefficients

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    integer :: i, j, k, l, m, n

! TODO
!    i = max(mulorder, polorder)
!    allocate(Cnij(2*i+3,0:i+1,0:i+1))
    allocate(Cnij(2*5+3,0:5+1,0:5+1))

    Cnij = 0
    Cnij(:,0,0) = 1
    do n = 1, 2*5+3
        if (mod(n,2) == 0) cycle
        do i = 1, 5+1
            if (mod(i,2) /= 0) then
                k = i - 1
            else if (mod(i,2) == 0) then
                k = i
            end if
            do j = 0, i
                if (mod(i+j,2) /= 0) cycle
                if (j == 0) then
                    Cnij(n,i,j) = Cnij(n,i-1,j+1)
                else if (j /= i) then
                    Cnij(n,i,j) = (j + 1) * Cnij(n,i-1,j+1)
                    Cnij(n,i,j) = Cnij(n,i,j) - (n + k) * Cnij(n,i-1,j-1)
                    k = k + 2
                else if (j == i) then
                    Cnij(n,i,j) = - (n + k) * Cnij(n,i-1,j-1)
                end if
            end do
        end do
    end do

end subroutine Tk_coefficients

!------------------------------------------------------------------------------

function T(Rij, x, y, z)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    integer, intent(in) :: x, y, z
    real(dp), dimension(3), intent(in) :: Rij

    integer :: l, m, n
    real(dp) :: T
    real(dp) :: R, Cx, Cy, Cz

    if (.not. allocated(Cnij)) call Tk_coefficients

    R = nrm2(Rij)

    T = 0.0d0

    do l = 0, x
        Cx = Cnij(1,x,l)*(Rij(1)/R)**l
        do m = 0, y
            Cy = Cx * Cnij(l+x+1,y,m)*(Rij(2)/R)**m
            do n = 0, z
                Cz = Cy * Cnij(l+x+m+y+1,z,n)*(Rij(3)/R)**n
                T = T + Cz
            end do
        end do
    end do
    T = T / R**(x+y+z+1)

end function T

!------------------------------------------------------------------------------

subroutine Tk_tensor(Tk, Rij)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    real(dp), dimension(:), intent(out) :: Tk
    real(dp), dimension(3), intent(in) :: Rij

    integer :: k, i
    integer :: x, y, z

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Tk)) - 1.0d0)) - 1

    i = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                Tk(i) = T(Rij, x, y, z)
                i = i + 1
            end do
        end do
    end do

end subroutine Tk_tensor

!------------------------------------------------------------------------------

subroutine Mk_integrals(Mk_ints, Rij, Mk)

    external :: Tk_integrals

    real(dp), dimension(:,:), intent(out) :: Mk_ints
    real(dp), dimension(:), intent(in) :: Mk
    real(dp), dimension(3), intent(in) :: Rij

    integer :: i, k
    integer :: ncomps
    real(dp) :: taylor
    real(dp), dimension(:), allocatable :: factors

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Mk)) - 1.0d0)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k)
    end if

    ncomps = size(Mk_ints, 2)

    call Tk_integrals('es', Mk_ints, nnbas, ncomps, Rij)

    ! get symmetry factors
    allocate(factors(ncomps))
    call symmetry_factors(factors)

    ! dot T^(k) integrals with multipole to get M^(k) integrals
    do i = 1, ncomps
        Mk_ints(:,i) = taylor * factors(i) * Mk(i) * Mk_ints(:,i)
    end do

    deallocate(factors)

end subroutine Mk_integrals

!------------------------------------------------------------------------------

subroutine symmetry_factors(factors)

    real(dp), dimension(:), intent(out) :: factors

    integer :: idx, x, y, z, k

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(factors)) - 1.0d0)) - 1

    idx = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                factors(idx) = symfac(x, y, z)
                idx = idx + 1
            end do
        end do
     end do

end subroutine symmetry_factors

!------------------------------------------------------------------------------

function symfac(i, j, k)

    ! trinomial coefficient

    integer, intent(in) :: i, j, k

    integer :: symfac

    symfac = factorial(i+j+k) / (factorial(i) * factorial(j) * factorial(k))

end function symfac

!------------------------------------------------------------------------------

recursive function factorial(n) result(nfact)

    ! Clive Page, http://www.star.le.ac.uk/~cgp/f90course/f90.html

    integer, intent(in) :: n

    integer :: nfact

    if (n > 0) then
        nfact = n * factorial(n-1)
    else
        nfact = 1
    end if

end function factorial

!------------------------------------------------------------------------------

subroutine openfile(filename, lunit, stat, frmt)

    character(*), intent(in) :: filename, stat, frmt
    integer, intent(out) :: lunit
    integer :: i
    logical :: lexist, lopen

    if (stat == 'old') then
        inquire(file=filename, exist=lexist)

        if (.not. lexist) then
            print *, filename, ' not found!'
            stop
        end if
    end if

    do i = 21, 99
        inquire(unit=i, opened=lopen)
        if (lopen) then
            cycle
        else
            lunit = i
            open(unit=lunit, file=filename, status=stat, form=frmt)
            exit
        end if
    end do

    return

end subroutine openfile

!------------------------------------------------------------------------------

function elem2charge(elem) result(charge)

    character(len=*), intent(in) :: elem

    integer :: i
    real(dp) :: charge
    character(len=2), dimension(112) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',&
                & 'Rg', 'Cn' /)

    if (elem == 'X') then
        charge = 0.0d0
        return
    end if

    do i = 1, 112
        if (elem == trim(elements(i))) then
            charge = real(i, dp)
            exit
        else
            charge = 0.0d0
        end if
    end do

end function elem2charge

!------------------------------------------------------------------------------

subroutine chcase(string, uplo)

    character(len=*), intent(inout) :: string
    character(len=*), intent(in), optional :: uplo

    integer :: i, gap
    character(len=1) :: a, z, o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'u'
    end if

    gap = iachar('a') - iachar('A')

    if (o_uplo == 'u' .or. o_uplo == 'U') then
        a = 'a'
        z = 'z'
    else if (o_uplo == 'l' .or. o_uplo == 'L') then
        a = 'A'
        z = 'Z'
        gap = - gap
    else
        stop 'Unknown case specified'
    end if

    do i = 1, len_trim(string)
        if (lge(string(i:i), a) .and. lle(string(i:i), z)) then
            string(i:i) = achar(iachar(string(i:i)) - gap)
        end if
    end do

end subroutine chcase

!------------------------------------------------------------------------------

subroutine compute_cube(denmats)

    real(dp), dimension(:), intent(in) :: denmats

    character(len=1) :: tcl
    character(len=99) :: cl
    integer :: i, j, k, l
    integer :: lu
    logical :: exclude, lexist
    real(dp), dimension(3) :: Rsp
    real(dp), dimension(:), allocatable :: Vpe
    real(dp), dimension(:,:), allocatable :: Fpe, M1inds

    if (myid == 0) then
        allocate(Vpe(npoints))
    else if (myid /= 0) then
        allocate(Vpe(cubedists(myid)))
    end if

    Vpe = 0.0d0

    if (mulorder >= 0) then
        k = 1
        do i = cube_start, cube_finish
            do j = 1, nsites
                Rsp = Rp(:,i) - Rs(:,j)
                if (lmul(0)) then
                    call multipole_potential(Vpe(k), Rsp, M0s(:,j))
                end if
                if (lmul(1)) then
                    call multipole_potential(Vpe(k), Rsp, M1s(:,j))
                end if
                if (lmul(2)) then
                    call multipole_potential(Vpe(k), Rsp, M2s(:,j))
                end if
                if (lmul(3)) then
                    call multipole_potential(Vpe(k), Rsp, M3s(:,j))
                end if
                if (lmul(4)) then
                    call multipole_potential(Vpe(k), Rsp, M4s(:,j))
                end if
                if (lmul(5)) then
                    call multipole_potential(Vpe(k), Rsp, M5s(:,j))
                end if
            end do
            k = k + 1
        end do
    end if

    if (lpol(1)) then
        allocate(M1inds(3*npols,1))
        if (myid == 0) then
            inquire(file='pe_induced_dipoles.bin', exist=lexist)
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
        end if
#endif
        if (lexist) then
            if (myid == 0) then
                call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
                rewind(lu)
                read(lu) M1inds
                close(lu)
            end if
        else
            write(luout,*) 'WARNING: pe_induced_dipoles.bin not found'
            write(luout,*) '         cannot create cube file.'
            return
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds, 3*npols, rmpi, 0, comm, ierr)
        end if
#endif
        k = 1
        do i = cube_start, cube_finish
            l = 1
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsp = Rp(:,i) - Rs(:,j)
                call multipole_potential(Vpe(k), Rsp, M1inds(l:l+2,1))
                l = l + 3
            end do
            k = k + 1
        end do
    end if

#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + cubedists(i-1)
        end do
        call mpi_gatherv(mpi_in_place, 0, rmpi, Vpe, cubedists, displs,&
                        & rmpi, 0, comm, ierr)
    else if (myid /= 0) then
        call mpi_gatherv(Vpe, cubedists(myid), rmpi, 0, 0, 0,&
                        & rmpi, 0, comm, ierr)
    end if
#endif
    if (myid == 0) then
        call openfile('embedding_potential.cube', lu, 'unknown', 'formatted')
        write(lu,'(a)') 'Polarizable embedding potential'
        write(lu,'(a)') 'Generated by the PE library in DALTON2013'
        write(lu,'(i5,3f12.6)') qmnucs, origin
        write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
        write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
        write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
        do j = 1, qmnucs
            write(lu,'(i5,4f12.6)') int(Zm(1,j)), Zm(1,j), Rm(:,j)
        end do
        do i = 1, xsteps * ysteps
            j = (i - 1) * zsteps + 1
            k = j - 1 + zsteps
            write(lu,'(6e13.5)') Vpe(j:k)
        end do
        close(lu)
    end if
 
    if (cube_field) then
        if (myid == 0) then
            allocate(Fpe(npoints,3))
        else if (myid /= 0) then
            allocate(Fpe(cubedists(myid),3))
        end if

        Fpe = 0.0d0

        if (mulorder >= 0) then
            k = 1
            do i = cube_start, cube_finish
                do j = 1, nsites
                    Rsp = Rp(:,i) - Rs(:,j)
                    if (lmul(0)) then
                        call multipole_field(Fpe(k,:), Rsp, M0s(:,j))
                    end if
                    if (lmul(1)) then
                        call multipole_field(Fpe(k,:), Rsp, M1s(:,j))
                    end if
                    if (lmul(2)) then
                        call multipole_field(Fpe(k,:), Rsp, M2s(:,j))
                    end if
                    if (lmul(3)) then
                        call multipole_field(Fpe(k,:), Rsp, M3s(:,j))
                    end if
                    if (lmul(4)) then
                        call multipole_field(Fpe(k,:), Rsp, M4s(:,j))
                    end if
                    if (lmul(5)) then
                        call multipole_field(Fpe(k,:), Rsp, M5s(:,j))
                    end if
                end do
                k = k + 1
            end do
        end if

        if (lpol(1)) then
            k = 1
            do i = cube_start, cube_finish
                l = 1
                do j = 1, nsites
                    if (zeroalphas(j)) cycle
                    Rsp = Rp(:,i) - Rs(:,j)
                    call multipole_field(Fpe(k,:), Rsp, M1inds(l:l+2,1))
                    l = l + 3
                end do
                k = k + 1
            end do
            deallocate(M1inds)
        end if

#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + 3 * cubedists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, 0, rmpi, Fpe, 3 * cubedists,&
                            & displs, rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Fpe, 3 * cubedists(myid), rmpi, 0, 0, 0, rmpi,&
                            & 0, comm, ierr)
        end if
#endif
        if (myid == 0) then
            do l = 1, 3
                write(cl,*) l
                tcl = trim(adjustl(cl))
                call openfile('embedding_field_'//tcl//'.cube', lu, 'unknown', 'formatted')
                write(lu,'(a)') 'Polarizable embedding electric field component '//tcl
                write(lu,'(a)') 'Generated by the PE library in DALTON2013'
                write(lu,'(i5,3f12.6)') qmnucs, origin
                write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
                write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
                write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
                do j = 1, qmnucs
                    write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
                end do
                do i = 1, xsteps * ysteps
                    j = (i - 1) * zsteps + 1
                    k = j - 1 + zsteps
                    write(lu,'(6e13.5)') Fpe(j:k,l)
                end do
                close(lu)
            end do
        end if
    end if

end subroutine compute_cube

end module polarizable_embedding
