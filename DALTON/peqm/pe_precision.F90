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
module pe_precision

!    integer, parameter :: sp = selected_real_kind(6, 37)
!    integer, parameter :: dp = selected_real_kind(15, 307)
!    integer, parameter :: qp = selected_real_kind(33, 4931)

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
    integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))

!    FORTRAN 2008
!    use, intrinsic :: iso_fortran_env
!    integer, parameter :: sp = REAL32
!    integer, parameter :: dp = REAL64
!    integer, parameter :: qp = REAL128

end module pe_precision
