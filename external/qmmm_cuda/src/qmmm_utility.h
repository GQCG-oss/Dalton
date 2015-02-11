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

#ifndef QMMM_UTILITY
#define QMMM_UTILITY

#include <vector>

typedef double qmmm_real;

#define QMMM_BASIS_FUNC_POLY_MAX_DEGREE 5
#define QMMM_BOYS_X_MAX 160.0
#define QMMM_BOYS_NO_OF_INTERVALS 200
#define QMMM_BOYS_TAB_DEGREE 12
#define QMMM_BOYS_N_MAX (QMMM_BASIS_FUNC_POLY_MAX_DEGREE*4+1)
#define QMMM_MONOMIAL_N_MAX (QMMM_BASIS_FUNC_POLY_MAX_DEGREE*4)
#define QMMM_HERMITE_CONVERSION_MAX_N (QMMM_BASIS_FUNC_POLY_MAX_DEGREE*2)

#define QMMM_BOYS_TEMP QMMM_BOYS_X_MAX / QMMM_BOYS_NO_OF_INTERVALS
#define QMMM_BOYS_INTERVAL_WIDTH_INV BOYS_NO_OF_INTERVALS / QMMM_BOYS_X_MAX;

struct qmmm_hermite_conversion_contrib_struct {
	int destIndex;
	int sourceIndex;
	int a_power;
	qmmm_real coeff;
};

struct qmmm_hermite_conversion_info_struct {
	qmmm_hermite_conversion_contrib_struct* list_right[QMMM_HERMITE_CONVERSION_MAX_N+1][QMMM_HERMITE_CONVERSION_MAX_N+1];
	int counters_right[QMMM_HERMITE_CONVERSION_MAX_N+1][QMMM_HERMITE_CONVERSION_MAX_N+1];

	qmmm_hermite_conversion_info_struct();
	~qmmm_hermite_conversion_info_struct();
};

struct qmmm_distribution_spec_struct {
	qmmm_real coeff;           /**< Coefficient A */
	qmmm_real exponent;        /**< exponent alfa */
	qmmm_real center_coords[3]; /**< x0, y0, z0    */
	char monomialInts[4];  /**< nx, ny, nz    */
};


struct qmmm_basis_info_struct {
	int no_of_basis_funcs;
	int* start_index; //no_of_basis_funcs + 1 (last element = last adress)
	qmmm_distribution_spec_struct* simple_primitive_list;
};

/* Simple atom representation by its charge and cartesian coordinates. */
struct qmmm_atom {
	qmmm_real charge;
	qmmm_real coords[3];
};

#ifndef pi
#define pi 3.14159265358979323846
#endif
#ifndef sqrtpi
#define sqrtpi 1.77245385090551602730
#endif

struct work_struct{
	int x;
	int y;
	int z;
	qmmm_real factor;
};

struct psi {
	qmmm_real coord_x;
	qmmm_real coord_y;
	qmmm_real coord_z;
	qmmm_real exponent;
	std::vector<int> monomial_l;
	std::vector<int> monomial_m;
	std::vector<int> monomial_k;
	std::vector<int> coeff;
};

/* integrals_general.h */
#define K_MAX_DIM 44
struct polydeg1struct {
	qmmm_real a0;
	qmmm_real a1;
};

struct qmmm_monomial_struct {
  int ix;
  int iy;
  int iz;
};

struct qmmm_monomial_info_struct {
  qmmm_monomial_struct* monomial_list;
  int no_of_monomials_list[QMMM_MONOMIAL_N_MAX+1];
  int monomial_index_list[QMMM_MONOMIAL_N_MAX+1][QMMM_MONOMIAL_N_MAX+1][QMMM_MONOMIAL_N_MAX+1];

  qmmm_monomial_info_struct();
  ~qmmm_monomial_info_struct();
};

int get_product_simple_prims(const qmmm_distribution_spec_struct& primA_in,
		const qmmm_distribution_spec_struct& primB_in,
		std::vector<qmmm_distribution_spec_struct*> &result_list,
		qmmm_real threshold);

struct boys_func_interval_struct {
	qmmm_real midx;
	qmmm_real a[QMMM_BOYS_TAB_DEGREE];
};

int qmmm_boysfunction_init(boys_func_interval_struct boys_list[QMMM_BOYS_N_MAX][QMMM_BOYS_NO_OF_INTERVALS]);
qmmm_real qmmm_boys_pre_factor(int n);

#endif
