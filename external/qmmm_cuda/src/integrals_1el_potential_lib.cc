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
#include <stdlib.h>
#include <string.h> //for memcpy
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "integrals_1el_potential_lib.h"

static qmmm_real boys_function_pretabulated_opt(const int n_max, const qmmm_real arg,
		const boys_func_interval_struct boys_list[QMMM_BOYS_N_MAX][QMMM_BOYS_NO_OF_INTERVALS]) {

	if(arg < 0) {
		std::cout << "error in BoysFunction_pretabulated: (x < 0)\n" << std::endl;
		exit(0);
		return 0;
	}

	if(n_max >= QMMM_BOYS_X_MAX) {
		std::cout << "error in BoysFunction_pretabulated: (n >= BOYS_N_MAX)\n" << std::endl;
		exit(0);
		return 0;
	}

	if(arg >= QMMM_BOYS_X_MAX) {
		/* use "large x formula" */
		//return (semiFactorial(2*n_max-1) / std::pow((ergo_real)2, n_max+1)) * std::sqrt(pi / std::pow(arg, 2*n_max+1));

		//For small n_max, it is faster to loop
		qmmm_real arg_exp = 1;
		for(int k = 0; k < 2*n_max + 1; k++) {
			arg_exp *= arg;
		}

		return qmmm_boys_pre_factor(n_max) / sqrt(arg_exp);
		//return boys_pre_factor(n_max) / sqrt(std::pow(arg, 2*n_max +1));
	}

	/* choose which interval to use */
	const int intervalIndex = arg * (qmmm_real) (QMMM_BOYS_NO_OF_INTERVALS / QMMM_BOYS_X_MAX);
	if((intervalIndex < 0) || (intervalIndex >= QMMM_BOYS_NO_OF_INTERVALS)) 	{
		std::cout << "error in BoysFunction_pretabulated: bad intervalIndex" << std::endl;
		std::cout << "intervalIndex = " << intervalIndex << std::endl;
		std::cout << "x = " << (double)arg << std::endl;
		std::cout << "count = " << (qmmm_real)arg / ((qmmm_real)QMMM_BOYS_X_MAX / QMMM_BOYS_NO_OF_INTERVALS) << std::endl;
		exit(0);
		return 0;
	}

	const boys_func_interval_struct* interval = &boys_list[n_max][intervalIndex];
	const qmmm_real deltax = arg - interval->midx;
	qmmm_real deltaxtopowk = 1;
	qmmm_real sum = 0;
	for(int k = 0; k < QMMM_BOYS_TAB_DEGREE; k++) {
		sum += interval->a[k] * deltaxtopowk;
		deltaxtopowk *= deltax;
	}
	return sum;
}

static qmmm_real R_recursive(const int n, const int n_max, const int ix, const int iy, const int iz, const qmmm_real* R,
		const qmmm_real dx0, const qmmm_real dx1, const qmmm_real dx2, const qmmm_real factor, const qmmm_real alpha0)
{
	if(ix == 0 && iy == 0 && iz == 0) {
		return (factor*R[n]);
	}

	qmmm_real Rval = 0;
	qmmm_real fac = factor*-2*alpha0;
	if(ix > 0) {
			Rval = dx0 * R_recursive(n+1,n_max,ix-1,iy,iz,R,dx0,dx1,dx2,fac,alpha0);
		if(ix > 1)
			Rval += (ix - 1) * R_recursive(n+1,n_max,ix-2,iy,iz,R,dx0,dx1,dx2,fac,alpha0);
	} else if(iy > 0) {
		Rval = dx1 * R_recursive(n+1,n_max,ix,iy-1,iz,R,dx0,dx1,dx2,fac,alpha0);
		if(iy > 1)
			Rval += (iy - 1) * R_recursive(n+1,n_max,ix,iy-2,iz,R,dx0,dx1,dx2,fac,alpha0);
	} else if(iz > 0) {
		Rval = dx2 * R_recursive(n+1,n_max,ix,iy,iz-1,R,dx0,dx1,dx2,fac,alpha0);
		if(iz > 1)
			Rval += (iz - 1) * R_recursive(n+1,n_max,ix,iy,iz-2,R,dx0,dx1,dx2,fac,alpha0);
	}
	return Rval;
}

int compute_V_matrix_full_opt_lib(const qmmm_basis_info_struct& basis_info,
		const int num_qm_atoms,
		const qmmm_atom* qm_atom_list,
		const int num_mm_atoms,
		const qmmm_atom* mm_atom_list,
		const qmmm_real threshold,
		qmmm_real* result) {

	const int nbast = basis_info.no_of_basis_funcs;

	std::cout << "Enter compute_V_matrix_full_opt_lib" << std::endl;
	std::cout << "Number of basis functions = " << nbast << ",  number of qm atoms = " << num_qm_atoms << ", number of mm atoms " <<  num_mm_atoms << std::endl;

	const qmmm_monomial_info_struct* monomial_info = new qmmm_monomial_info_struct();
	const qmmm_hermite_conversion_info_struct* hermite_conversion_info =  new qmmm_hermite_conversion_info_struct();

	boys_func_interval_struct boys_list[QMMM_BOYS_N_MAX][QMMM_BOYS_NO_OF_INTERVALS];
	qmmm_boysfunction_init(boys_list);

	int cntr = 0;


#ifdef _OPENMP
	double time_start, time_stop;
	time_start = omp_get_wtime();
#endif

	const int iterations = nbast*(nbast+1)/2;

#ifdef _OPENMP
#pragma omp parallel reduction(+:cntr)
	{
#endif
		std::vector<qmmm_distribution_spec_struct*> psi_list;
#ifdef _OPENMP
#pragma omp master
		{
			std::cout << "Number of threads used by openMP = " << omp_get_num_threads() << std::endl;
		}
#endif

#pragma omp for schedule(dynamic, 200)
		for(int i = 0; i < iterations; i++) {
			const int mu = sqrt(2*i + 0.25) - 0.5;
			const int nu = i - (mu * (mu + 1)) / 2;
			//const int result_index = mu*nbast + nu;

			//Get first contracted gaussian
			int start_prim_mu = basis_info.start_index[mu];
			int stop_prim_mu = basis_info.start_index[mu+1];

			//Get second contracted gaussian
			int start_prim_nu = basis_info.start_index[nu];
			int stop_prim_nu = basis_info.start_index[nu+1];

			qmmm_real sum = 0;
			/* compute matrix element [mu,nu] */
			psi_list.clear();
			for(int j = start_prim_mu; j < stop_prim_mu; j++) {
				for(int k = start_prim_nu; k < stop_prim_nu; k++) {
					get_product_simple_prims(basis_info.simple_primitive_list[j], basis_info.simple_primitive_list[k], psi_list, threshold);
				}
			}

			for(unsigned int p = 0; p < psi_list.size(); p++) {
				// Center_coords and alpha are not same for all items in psi_list any more (longer list)

				const qmmm_distribution_spec_struct* psi = psi_list[p];

				// Center coords and alpha are the same for many iterations
				const qmmm_real* center_coords = psi->center_coords;
				const qmmm_real alpha = psi->exponent;
				const qmmm_real inv_alpha = 1 / alpha;
				const qmmm_real resultPreFactor = 2 * pi * inv_alpha;

				const qmmm_real coeff = psi->coeff;

				const int n1x = psi->monomialInts[0];
				const int n1y = psi->monomialInts[1];
				const int n1z = psi->monomialInts[2];

				const int n1_max = n1x + n1y + n1z;
				const int n2_max = 0;
				const int Nmax = n1_max + n2_max;

				const int monomialIndex = monomial_info->monomial_index_list[n1x][n1y][n1z];

				const int noOfMonomials_1 = monomial_info->no_of_monomials_list[n1_max];
				const int noOfMonomials_2 = monomial_info->no_of_monomials_list[n2_max];

				int no_of_contribs = hermite_conversion_info->counters_right[n1_max][n2_max];
				qmmm_hermite_conversion_contrib_struct* hermite_conversion_list = hermite_conversion_info->list_right[n1_max][n2_max];

				int minus_1_to_pow_list[Nmax+1];
				int ifactor = 1;
				for(int n = 0; n <= Nmax; n++) 	{
					minus_1_to_pow_list[n] = ifactor;
					ifactor *= -1;
				}

				int Ntot = n1_max + n2_max;
				qmmm_real inv_alpha_pow_list[Ntot+1];
				inv_alpha_pow_list[0] = 1;
				for(int m = 1; m <= Ntot; m++) {
					inv_alpha_pow_list[m] = inv_alpha_pow_list[m-1] * inv_alpha;
				}

				std::vector<work_struct> work;
				for(int i=0; i < no_of_contribs; i++) {
					if(hermite_conversion_list[i].destIndex == monomialIndex) {
						if(noOfMonomials_1 > 0 && noOfMonomials_2 > 0) {

							const int i1 = hermite_conversion_list[i].sourceIndex / noOfMonomials_2;
							const int i2 = hermite_conversion_list[i].sourceIndex % noOfMonomials_2;

							int x = monomial_info->monomial_list[i1].ix;
							int y = monomial_info->monomial_list[i1].iy;
							int z = monomial_info->monomial_list[i1].iz;
							const int n = x+y+z;

							const qmmm_real prefactor = minus_1_to_pow_list[n] * resultPreFactor;

							x += monomial_info->monomial_list[i2].ix;
							y += monomial_info->monomial_list[i2].iy;
							z += monomial_info->monomial_list[i2].iz;

							const qmmm_real tmp = prefactor * hermite_conversion_list[i].coeff * inv_alpha_pow_list[-hermite_conversion_list[i].a_power];
							//work_struct tmp_struct = {x, y, z, tmp};//(x, y, z, tmp);
							work_struct tmp_struct;
							tmp_struct.x = x;
							tmp_struct.y = y;
							tmp_struct.z = z;
							tmp_struct.factor = tmp;

							work.push_back(tmp_struct);
							//work.push_back(new work_struct());
						}
					}
				}

				//TODO  Create special case if my_vector.size == 0, if sum x, y, z == 0, if noOfContribs == 0
				for(int l = 0; l < num_mm_atoms; l++) {
					const qmmm_atom* atom = &mm_atom_list[l];
					const qmmm_real point_charge = atom->charge;
					const qmmm_real* point_charge_coords = atom->coords;

					const qmmm_real dx0 = point_charge_coords[0] - center_coords[0];
					const qmmm_real dx1 = point_charge_coords[1] - center_coords[1];
					const qmmm_real dx2 = point_charge_coords[2] - center_coords[2];
					const qmmm_real R_12_squared = dx0*dx0 + dx1*dx1 + dx2*dx2;

					/* Compute all Boys function values needed */
					/* Use downward recursion to get Boys function values */
					const qmmm_real arg = alpha * R_12_squared;

					qmmm_real BoysList[Nmax+1];
					BoysList[Nmax] = boys_function_pretabulated_opt(Nmax, arg, &boys_list[0]);

					if(Nmax > 0) {
						const qmmm_real exp_minus_arg = std::exp(-arg);
						for(int n = Nmax-1; n >= 0; n--) {
							BoysList[n] = (2*arg*BoysList[n+1] + exp_minus_arg) / (2*n+1);
						}
					}

					qmmm_real res = 0;
					for(unsigned int r = 0; r < work.size(); r++) {
						int x = work[r].x;
						int y = work[r].y;
						int z = work[r].z;
						qmmm_real factor  = work[r].factor;
						res +=  R_recursive(0,Nmax,x,y,z,BoysList,dx0,dx1,dx2,1,alpha) * factor;
					}

					sum += res * point_charge * coeff;
				} /* END MM ATOM */

			} /* END PRIMS */
			result[mu*nbast+nu] = -1 * sum;
		} /* END FOR iterations */
#ifdef _OPENMP
	} // end parallel region
#endif

	std::cout << "cntr = " << cntr << std::endl;

	// copy values to the other triangle
	for(int mu = 0; mu < nbast; mu++)
		for(int nu = mu+1; nu < nbast; nu++)
			result[mu*nbast+nu] = result[nu*nbast+mu];

#ifdef _OPENMP
	time_stop = omp_get_wtime();
	std::cout << "Execution time of compute_V_matrix_full_opt_lib = " << time_stop - time_start << std::endl;
#endif

	std::cout <<"Exit compute_V_matrix_full_opt_lib" << std::endl;

	return 0;
}
