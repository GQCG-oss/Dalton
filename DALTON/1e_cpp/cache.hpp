/*
    Copyright 2013 Jaime Axel Rosal Sandberg

    This file is part of the EFS library.

    The EFS library is free software:  you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The EFS library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the EFS library.  If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef __CACHE__
#define __CACHE__

#define CACHE_LINE_SIZE       64
#define AVXD_PER_CACHE_LINE    2
#define MM128_PER_CACHE_LINE   4
#define DOUBLES_PER_CACHE_LINE 8
#define FLOATS_PER_CACHE_LINE 16

class cacheline {
  public:
    char d2[CACHE_LINE_SIZE];
}  __attribute__((aligned(CACHE_LINE_SIZE)));

#endif

