/* The GPU Optimized QMMM Library, version 0.1, a Portable GPU library for evaluation of Coulomb interaction
 * integrals between electrons in QM region of the system and point charges in MM
 * region of the system. This library is designed to accelerate hybrid QM/MM computations
 * in quantum chemistry programs. It can be integrated into any quantum chemistry program,
 * which supports QM/MM computations, via standard C interface.
 *
 * Copyright (C) 2013 Mikael Engbom
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information about The GPU Optimized Library, see http://www.scalalife.eu/content/portable-gpu-libary-qm-mm-calculations
 */

/* The GPU optimized QMMM library, version 0.1, is based on the Ergo 3.1 source code */

/* Ergo, version 3.1, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2011 Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 *
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

#ifndef INTEGRALS_1EL_POTENTIAL_CUDA_LIB
#define INTEGRALS_1EL_POTENTIAL_CUDA_LIB

#include "qmmm_utility.h"

#define DX 0
#define DY 1
#define DZ 2

#define DX_SHIFT 0
#define DY_SHIFT 1
#define DZ_SHIFT 2

#define IS_DX(byte) (byte & 0x04 ? 1 : 0)
#define IS_DY(byte) (byte & 0x02 ? 1 : 0)
#define IS_DZ(byte) (byte & 0x01 ? 1 : 0)

#define GET_DX(byte, dx) (byte & 0x04 ? dx : 1)
#define GET_DY(byte, dy) (byte & 0x02 ? dy : 1)
#define GET_DZ(byte, dz) (byte & 0x01 ? dz : 1)

#define NEW_NMAX 3333
#define ADD 4444
#define NEW_CENTER 5555
#define NEW_ALPHA 6666

#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

int compute_V_matrix_full_gpu_lib(const qmmm_basis_info_struct& basis_info,
		const int num_qm_atoms,
		const qmmm_atom* qm_atom_list,
		const int num_mm_atoms,
		const qmmm_atom* mm_atom_list,
		const qmmm_real threshold,
		qmmm_real* result);


#endif /* INTEGRALS_1EL_POTENTIAL_CUDA_LIB */
