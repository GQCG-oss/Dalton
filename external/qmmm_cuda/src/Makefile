# The GPU Optimized QMMM Library, version 0.1, a Portable GPU library for evaluation of Coulomb interaction
# integrals between electrons in QM region of the system and point charges in MM
# region of the system. This library is designed to accelerate hybrid QM/MM computations
# in quantum chemistry programs. It can be integrated into any quantum chemistry program,
# which supports QM/MM computations, via standard C interface.
 
# Copyright (C) 2013 Mikael Engbom
 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
# For further information about The GPU Optimized Library, see http://www.scalalife.eu/content/portable-gpu-libary-qm-mm-calculations

#makefile for qmmm_library_v0.1

CC = g++
OPT = -O3 -Wall -fopenmp -c
AR = ar
AR_OPT = cvq

NVCC = nvcc
NVCC_OPT = -arch=sm_20 -Xcompiler -fopenmp -c
NVCC_INC = -I/usr/local/cuda/include/ -I.
NVCC_LIB = -L/usr/local/cuda/lib64 -lcudart

qmmm_lib.ar: integrals_1el_potential_lib.o qmmm_utility.o integrals_1el_potential_cuda_lib.o
	$(AR) $(AR_OPT) -o $@ $+

integrals_1el_potential_lib.o: integrals_1el_potential_lib.cc 
	$(CC) $(OPT) -o $@ $+

integrals_1el_potential_cuda_lib.o: integrals_1el_potential_cuda_lib.cu
	$(NVCC) $(NVCC_OPT) $(NVCC_INC) $(NVCC_LIB) -o $@ $+
	
qmmm_utility.o: qmmm_utility.cc 
	$(CC) $(OPT) -o $@ $+

clean:
	rm *.o qmmm_lib.ar
