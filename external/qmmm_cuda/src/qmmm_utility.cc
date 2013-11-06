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

#include "integrals_1el_potential_lib.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>


qmmm_monomial_info_struct::qmmm_monomial_info_struct() {
	// first get count
	int count = 0;
	for(int n1 = 0; n1 <= QMMM_MONOMIAL_N_MAX; n1++)
	{
		for(int n1x = 0; n1x <= n1; n1x++)
			for(int n1y = 0; n1y <= n1; n1y++)
				for(int n1z = 0; n1z <= n1; n1z++)
				{
					if(n1x+n1y+n1z != n1)
						continue;
					count++;
				} /* END FOR n1x n1y n1z */
	} /* END FOR n1 */

	int noOfMonomialsTot = count;
	monomial_list = new qmmm_monomial_struct[noOfMonomialsTot];

	count = 0;
	for(int n1 = 0; n1 <= QMMM_MONOMIAL_N_MAX; n1++)
	{
		for(int n1x = 0; n1x <= n1; n1x++)
			for(int n1y = 0; n1y <= n1; n1y++)
				for(int n1z = 0; n1z <= n1; n1z++)
				{
					if(n1x+n1y+n1z != n1)
						continue;
					assert(count < noOfMonomialsTot);
					monomial_list[count].ix = n1x;
					monomial_list[count].iy = n1y;
					monomial_list[count].iz = n1z;
					monomial_index_list[n1x][n1y][n1z] = count;
					count++;
				} /* END FOR n1x n1y n1z */
		no_of_monomials_list[n1] = count;
	} /* END FOR n1 */
	assert(count == noOfMonomialsTot);
}


qmmm_monomial_info_struct::~qmmm_monomial_info_struct() {
  delete []monomial_list;
}


static int multiply_polynomials(qmmm_real* result, polydeg1struct* polydeg1, int dim, qmmm_real* a) {
	qmmm_real p1[K_MAX_DIM + 1];
	qmmm_real p2[K_MAX_DIM + 1];
	if(dim >= (K_MAX_DIM-1))
		return -1;
	for(int i = 0; i <= dim; i++)
		p1[i] = a[i]*polydeg1->a0;
	p1[dim+1] = 0;
	p2[0] = 0;
	for(int i = 0; i <= dim; i++)
		p2[i+1] = a[i]*polydeg1->a1;
	for(int i = 0; i <= (dim+1); i++)
		result[i] = p1[i] + p2[i];
	return 0;
} /* END multiply_polynomials */


/*
get_product_simple_prims
This function calculates the product of two simple primitives.
The result is a list of simple primitives.
*/
int get_product_simple_prims(const qmmm_distribution_spec_struct& primA_in, const qmmm_distribution_spec_struct& primB_in,
		std::vector<qmmm_distribution_spec_struct*> &result_list, qmmm_real threshold) {

	// Use a coordinate system with primA at the origin.
	// This solves the problem with extreme positions of the primitives.

	qmmm_distribution_spec_struct primA = primA_in;
	qmmm_distribution_spec_struct primB = primB_in;

	for(int kk = 0; kk < 3; kk++) {
		primA.center_coords[kk] -= primA_in.center_coords[kk];
		primB.center_coords[kk] -= primA_in.center_coords[kk];
	}

	qmmm_real CxCyCz, AiAj, alphaNew;
	qmmm_real newCenter[3];

	qmmm_real poly0[K_MAX_DIM];
	qmmm_real poly1[K_MAX_DIM];
	qmmm_real poly2[K_MAX_DIM];

	qmmm_real tempPoly[K_MAX_DIM];
	qmmm_real tempPoly2[K_MAX_DIM];
	qmmm_real tempPoly3[K_MAX_DIM];

	int tempPolyDegree, tempPoly2Degree;

	int poly0degree, poly1degree, poly2degree;

	polydeg1struct polyDeg1;
	qmmm_real* poly;
	int* degreePtr;

	/* use the Gaussian product rule */
	qmmm_real sum = 0;
	for(int k = 0; k < 3; k++) {
		qmmm_real temp = primA.center_coords[k] - primB.center_coords[k];
		sum += temp * temp;
	} /* END FOR k */
	CxCyCz = std::exp(-primA.exponent * primB.exponent * sum / (primA.exponent + primB.exponent));

	// FIXME: do this screening properly!
	if(std::fabs(CxCyCz) < threshold)
		return 0;

	AiAj = primA.coeff * primB.coeff;
	alphaNew = primA.exponent + primB.exponent;

	newCenter[0] = (primA.exponent * primA.center_coords[0] + primB.exponent * primB.center_coords[0]) / (alphaNew);
	newCenter[1] = (primA.exponent * primA.center_coords[1] + primB.exponent * primB.center_coords[1]) / (alphaNew);
	newCenter[2] = (primA.exponent * primA.center_coords[2] + primB.exponent * primB.center_coords[2]) / (alphaNew);

	/* do product of polynomials */
	/* one coordinate at a time */
	for(int k = 0; k < 3; k++) {
		switch(k)
		{
		case 0: poly = poly0; degreePtr = &poly0degree; break;
		case 1: poly = poly1; degreePtr = &poly1degree; break;
		case 2: poly = poly2; degreePtr = &poly2degree; break;
		default: return -1;
		} /* END SWITCH k */
		tempPoly[0] = 1;
		tempPolyDegree = 0;
		for(int m = 0; m < primA.monomialInts[k]; m++) {
			polyDeg1.a0 = -primA.center_coords[k];
			polyDeg1.a1 = 1;
			if(multiply_polynomials(tempPoly2, &polyDeg1, tempPolyDegree, tempPoly) != 0)
				return -1;
			tempPolyDegree++;
			memcpy(tempPoly, tempPoly2, (tempPolyDegree+1)*sizeof(qmmm_real));
		} /* END FOR m */
		for(int m = 0; m < primB.monomialInts[k]; m++) {
			polyDeg1.a0 = -primB.center_coords[k];
			polyDeg1.a1 = 1;
			if(multiply_polynomials(tempPoly2, &polyDeg1, tempPolyDegree, tempPoly) != 0)
				return -1;
			tempPolyDegree++;
			memcpy(tempPoly,tempPoly2,(tempPolyDegree+1)*sizeof(qmmm_real));
		} /* END FOR m */

		/* now do variable change */
		for(int m = 0; m < K_MAX_DIM; m++)
			poly[m] = 0;

		tempPoly2Degree = 0;
		for(int m = 0; m <= tempPolyDegree; m++) {
			tempPoly2[0] = tempPoly[m];
			tempPoly2Degree = 0;
			for(int l = 0; l < m; l++) {
				polyDeg1.a0 = newCenter[k];
				polyDeg1.a1 = 1;
				if(multiply_polynomials(tempPoly3, &polyDeg1, tempPoly2Degree, tempPoly2) != 0)
					return -1;
				tempPoly2Degree++;
				memcpy(tempPoly2, tempPoly3, (tempPoly2Degree+1)*sizeof(qmmm_real));
			} /* END FOR l */
			for(int l = 0; l <= tempPoly2Degree; l++) {
				poly[l] += tempPoly2[l];
			} /* END FOR l */
		} /* END FOR m */
		*degreePtr = tempPoly2Degree;
	} /* END FOR k */

	for(int k = 0; k <= poly0degree; k++) {
		for(int l = 0; l <= poly1degree; l++) {
			for(int m = 0; m <= poly2degree; m++) {
				qmmm_real newCoeff = AiAj * CxCyCz * poly0[k] * poly1[l] * poly2[m];

				qmmm_real sqrtValue = std::sqrt(pi / alphaNew);
				qmmm_real absvalue = newCoeff * sqrtValue * sqrtValue * sqrtValue;
				if(absvalue < 0) absvalue *= -1;

				/* add one function to final list */
				result_list.push_back(new qmmm_distribution_spec_struct());

				// Translate this term of result back to original coordinate system
				result_list.back()->center_coords[0] = newCenter[0] + primA_in.center_coords[0];
				result_list.back()->center_coords[1] = newCenter[1] + primA_in.center_coords[1];
				result_list.back()->center_coords[2] = newCenter[2] + primA_in.center_coords[2];

				result_list.back()->coeff = newCoeff;
				result_list.back()->exponent = alphaNew;

				result_list.back()->monomialInts[0]= k;
				result_list.back()->monomialInts[1]= l;
				result_list.back()->monomialInts[2]= m;

			} /* END FOR m */
		} /* END FOR l */
	} /* END FOR k */

	return 0;
}


double semi_factorial(int n) {
	switch(n)
	{
	case -1:
		return 1;
	case 0:
		return 1;
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		return 3;
	case 4:
		return 4*2;
	case 5:
		return 5*3;
	case 6:
		return 6*4*2;
	case 7:
		return 7*5*3;
	case 8:
		return 8*6*4*2;
	case 9:
		return 9*7*5*3;
	case 10:
		return 10*8*6*4*2;
	case 11:
		return 11*9*7*5*3;
	case 12:
		return 12*10*8*6*4*2;
	case 13:
		return 13*11*9*7*5*3;
	case 14:
		return 14*12*10*8*6*4*2;
	case 15:
		return 15*13*11*9*7*5*3;
	case 16:
		return 16*14*12*10*8*6*4*2;
	case 17:
		return 17*15*13*11*9*7*5*3;
	case 18:
		return 18*16*14*12*10*8*6*4*2;
	case 19:
		return 19*(double)17*15*13*11*9*7*5*3;
	case 20:
		return 20*(double)18*16*14*12*10*8*6*4*2;
	case 21:
		return 21*19*(double)17*15*13*11*9*7*5*3;
	case 22:
		return 22*20*(double)18*16*14*12*10*8*6*4*2;
	case 23:
		return 23*21*19*(double)17*15*13*11*9*7*5*3;
	case 24:
		return 24*22*20*(double)18*16*14*12*10*8*6*4*2;
	default:
		std::cout << "error: semiFactorial not implemented for n > 24" << std::endl;
		std::cout << "n = "<< n << std::endl;
		exit(0);
		return 0;
	}
}

qmmm_real qmmm_boys_pre_factor(int n) {
	switch(n)
	{
	case 0:
		return 0.5 * sqrtpi;
	case 1:
		return 0.25 * sqrtpi;
	case 2:
		return 0.375 * sqrtpi;
	case 3:
		return 0.9375 * sqrtpi;
	case 4:
		return 3.28125 * sqrtpi;
	case 5:
		return 14.765625 * sqrtpi;
	default:
		return (semi_factorial(2*n-1) / std::pow((qmmm_real)2, n+1)) * sqrtpi;
	}
}

qmmm_real boys_function_raw_simpson(int n, qmmm_real x) {
	const int N = 100;
	qmmm_real h = (qmmm_real)0.5 / N;
	qmmm_real sum = 0;
	for(int k = 0; k <= 2*N; k++)
	{
		qmmm_real tk = (qmmm_real)k / (2*N);
		// Compute f(tk) = exp(-x*tk*tk) * pow(tk, 2*n)
		qmmm_real foftk = std::exp(-x*tk*tk);
		if(n != 0) {
			if(k != 0)
				foftk *= std::pow(tk, 2*n);
			else
				foftk = 0;
		}
		// OK, foftk done, now add to sum.
		if(k == 0 || k == 2*N) {
			sum += foftk;
			continue;
		}
		if(k % 2 == 1) {
			sum += 4 * foftk;
			continue;
		}
		sum += 2 * foftk;
	}
	return (h/3) * sum;
}


#if BASIS_FUNC_POLY_MAX_DEGREE<6
const int QMMM_MAX_NO_OF_CONTRIBS = 1000000;
#else
const int QMMM_MAX_NO_OF_CONTRIBS = 10000000;
#endif

struct symb_matrix_element{
  int ia; // power of a
  qmmm_real coeff;
};

struct poly_1d_term_struct_symb {
  int ix; // power of x
  int ia; // power of a
  qmmm_real coeff;
};

#define MAX_NO_OF_1D_TERMS 888

struct poly_1d_struct_symb {
  int noOfTerms;
  poly_1d_term_struct_symb termList[MAX_NO_OF_1D_TERMS];
};


static int get_1d_hermite_poly_inv_symb(poly_1d_struct_symb* result, int n) {
	switch(n)
	{
	case 0:
		result->noOfTerms = 1;
		result->termList[0].ix = 0;
		result->termList[0].ia = 0;
		result->termList[0].coeff = 1;
		break;
	case 1:
		result->noOfTerms = 1;
		result->termList[0].ix = 1;
		result->termList[0].ia = -1;
		result->termList[0].coeff = 0.5;
		break;
	default:
	{
		// Create polys for n-1 and n-2
		poly_1d_struct_symb poly_n_m_1;
		poly_1d_struct_symb poly_n_m_2;
		get_1d_hermite_poly_inv_symb(&poly_n_m_1, n - 1);
		get_1d_hermite_poly_inv_symb(&poly_n_m_2, n - 2);
		assert(poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms < MAX_NO_OF_1D_TERMS);
		// Now the result is 0.5*(1/a)*x*poly_n_m_1 + (n-1)*(1/a)*0.5*poly_n_m_2
		for(int i = 0; i < poly_n_m_1.noOfTerms; i++)
		{
			result->termList[i] = poly_n_m_1.termList[i];
			result->termList[i].ix++;
			result->termList[i].ia--;
			result->termList[i].coeff *= 0.5;
		}
		int nn = poly_n_m_1.noOfTerms;
		for(int i = 0; i < poly_n_m_2.noOfTerms; i++)
		{
			result->termList[nn+i] = poly_n_m_2.termList[i];
			result->termList[nn+i].ia--;
			result->termList[nn+i].coeff *= (n-1) * 0.5;
		}
		result->noOfTerms = poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms;
	}
	}
	return 0;
}

static int get_1d_hermite_poly_symb(poly_1d_struct_symb* result, int n) {
	switch(n)
	{
	case 0:
		result->noOfTerms = 1;
		result->termList[0].ix = 0;
		result->termList[0].ia = 0;
		result->termList[0].coeff = 1;
		break;
	case 1:
		result->noOfTerms = 1;
		result->termList[0].ix = 1;
		result->termList[0].ia = 1;
		result->termList[0].coeff = 2;
		break;
	default:
	{
		// Create polys for n-1 and n-2
		poly_1d_struct_symb poly_n_m_1;
		poly_1d_struct_symb poly_n_m_2;
		get_1d_hermite_poly_symb(&poly_n_m_1, n - 1);
		get_1d_hermite_poly_symb(&poly_n_m_2, n - 2);
		assert(poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms < MAX_NO_OF_1D_TERMS);
		// Now the result is 2*a*x*poly_n_m_1 - (n-1)*2*a*poly_n_m_2
		for(int i = 0; i < poly_n_m_1.noOfTerms; i++) {
			result->termList[i] = poly_n_m_1.termList[i];
			result->termList[i].ix++;
			result->termList[i].ia++;
			result->termList[i].coeff *= 2;
		}
		int nn = poly_n_m_1.noOfTerms;
		for(int i = 0; i < poly_n_m_2.noOfTerms; i++) {
			result->termList[nn+i] = poly_n_m_2.termList[i];
			result->termList[nn+i].ia++;
			result->termList[nn+i].coeff *= -2 * (n-1);
		}
		result->noOfTerms = poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms;
	}
	}
	return 0;
}

struct poly_3d_term_struct_symb {
  int monomialInts[3];
  int ia; // power of a
  qmmm_real coeff;
};

#define MAX_NO_OF_3D_TERMS 888

struct poly_3d_struct_symb {
  int noOfTerms;
  poly_3d_term_struct_symb termList[MAX_NO_OF_3D_TERMS];
};

static int create_3d_poly_from_1d_poly_symb(poly_3d_struct_symb* poly_3d,
				 poly_1d_struct_symb* poly_1d,
				 int coordIndex) {
  memset(poly_3d, 0, sizeof(poly_3d_struct_symb));
  for(int i = 0; i < poly_1d->noOfTerms; i++)
    {
      poly_3d->termList[i].coeff = poly_1d->termList[i].coeff;
      poly_3d->termList[i].monomialInts[coordIndex] = poly_1d->termList[i].ix;
      poly_3d->termList[i].ia = poly_1d->termList[i].ia;
    }
  poly_3d->noOfTerms = poly_1d->noOfTerms;
  return 0;
}

static int compute_product_of_3d_polys_symb(poly_3d_struct_symb* result,
				 poly_3d_struct_symb* poly_1,
				 poly_3d_struct_symb* poly_2)
{
	int termCount = 0;
	int termidx_1, termidx_2;
	for(termidx_1 = 0; termidx_1 < poly_1->noOfTerms; termidx_1++)
		for(termidx_2 = 0; termidx_2 < poly_2->noOfTerms; termidx_2++)
		{
			poly_3d_term_struct_symb* term_1 = &poly_1->termList[termidx_1];
			poly_3d_term_struct_symb* term_2 = &poly_2->termList[termidx_2];

			// Create product term
			poly_3d_term_struct_symb newTerm;
			newTerm.coeff = term_1->coeff * term_2->coeff;
			for(int k = 0; k < 3; k++)
				newTerm.monomialInts[k] = term_1->monomialInts[k] + term_2->monomialInts[k];
			newTerm.ia = term_1->ia + term_2->ia;

			result->termList[termCount] = newTerm;
			termCount++;
			assert(termCount < MAX_NO_OF_3D_TERMS);
		} // END FOR termidx_1 termidx_2
	result->noOfTerms = termCount;
	return 0;
}


static int get_hermite_conversion_matrix_symb(const qmmm_monomial_info_struct* monomial_info,
				   int nmax,
				   int inverseFlag,
				   symb_matrix_element* result) {
	int noOfMonomials = monomial_info->no_of_monomials_list[nmax];
	memset(result, 0, noOfMonomials*noOfMonomials*sizeof(symb_matrix_element));

	int monomialIndex;
	for(monomialIndex = 0; monomialIndex < noOfMonomials; monomialIndex++) {

		// get monomialInts
		int ix = monomial_info->monomial_list[monomialIndex].ix;
		int iy = monomial_info->monomial_list[monomialIndex].iy;
		int iz = monomial_info->monomial_list[monomialIndex].iz;

		// Get x y z 1-d Hermite polynomials
		poly_1d_struct_symb hermitePoly_1d_x;
		poly_1d_struct_symb hermitePoly_1d_y;
		poly_1d_struct_symb hermitePoly_1d_z;
		if(inverseFlag == 1)
		{
			get_1d_hermite_poly_inv_symb(&hermitePoly_1d_x, ix);
			get_1d_hermite_poly_inv_symb(&hermitePoly_1d_y, iy);
			get_1d_hermite_poly_inv_symb(&hermitePoly_1d_z, iz);
		}
		else
		{
			get_1d_hermite_poly_symb(&hermitePoly_1d_x, ix);
			get_1d_hermite_poly_symb(&hermitePoly_1d_y, iy);
			get_1d_hermite_poly_symb(&hermitePoly_1d_z, iz);
		}

		// Store x y z Hermite polys as 3-d polys
		poly_3d_struct_symb hermitePoly_3d_x;
		poly_3d_struct_symb hermitePoly_3d_y;
		poly_3d_struct_symb hermitePoly_3d_z;
		create_3d_poly_from_1d_poly_symb(&hermitePoly_3d_x, &hermitePoly_1d_x, 0);
		create_3d_poly_from_1d_poly_symb(&hermitePoly_3d_y, &hermitePoly_1d_y, 1);
		create_3d_poly_from_1d_poly_symb(&hermitePoly_3d_z, &hermitePoly_1d_z, 2);

		// Compute product
		poly_3d_struct_symb hermitePoly_3d_xy;
		poly_3d_struct_symb hermitePoly_3d_xyz;
		compute_product_of_3d_polys_symb(&hermitePoly_3d_xy,
				&hermitePoly_3d_x,
				&hermitePoly_3d_y);
		compute_product_of_3d_polys_symb(&hermitePoly_3d_xyz,
				&hermitePoly_3d_xy,
				&hermitePoly_3d_z);

		// Go through result product poly, for each term get monomialIndex and add
		// coeff to final result at position given by monomialIndex.
		for(int i = 0; i < hermitePoly_3d_xyz.noOfTerms; i++)
		{
			poly_3d_term_struct_symb* currTerm = &hermitePoly_3d_xyz.termList[i];
			// Get monomialIndex2
			int ix = currTerm->monomialInts[0];
			int iy = currTerm->monomialInts[1];
			int iz = currTerm->monomialInts[2];
			int monomialIndex2 = monomial_info->monomial_index_list[ix][iy][iz];
			result[monomialIndex * noOfMonomials + monomialIndex2].coeff += currTerm->coeff;
			if(result[monomialIndex * noOfMonomials + monomialIndex2].ia != 0)
				assert(result[monomialIndex * noOfMonomials + monomialIndex2].ia == currTerm->ia);
			result[monomialIndex * noOfMonomials + monomialIndex2].ia = currTerm->ia;
		} // END FOR i

	} // END FOR monomialIndex

	return 0;
}

qmmm_hermite_conversion_info_struct::qmmm_hermite_conversion_info_struct() {

	const int nmax = QMMM_HERMITE_CONVERSION_MAX_N;

	qmmm_monomial_info_struct monomial_info;

	for(int n1 = 0; n1 <= nmax; n1++)
		for(int n2 = 0; n2 <= nmax; n2++) {
			qmmm_hermite_conversion_contrib_struct* currlist = new qmmm_hermite_conversion_contrib_struct[QMMM_MAX_NO_OF_CONTRIBS];
			int count = 0;
			int nMon1 = monomial_info.no_of_monomials_list[n1];
			int nMon2 = monomial_info.no_of_monomials_list[n2];
			symb_matrix_element* list = new symb_matrix_element[nMon1*nMon1];
			int inverseFlag = 1;
			get_hermite_conversion_matrix_symb(&monomial_info, n1, inverseFlag, list);
			for(int j = 0; j < nMon1; j++)
				for(int i = 0; i < nMon2; i++) {
					// Now take care of matrix element (i j)
					for(int k = 0; k < nMon1; k++) {
						int idx = j*nMon1+k;
						if(std::fabs(list[idx].coeff) > 1e-5) {
							assert(count < QMMM_MAX_NO_OF_CONTRIBS);
							currlist[count].destIndex = j*nMon2+i;
							currlist[count].sourceIndex = k*nMon2+i;
							currlist[count].a_power = list[idx].ia;
							currlist[count].coeff = list[idx].coeff;
							count++;
						}
					}
				} // END FOR i j
			list_right[n1][n2] = new qmmm_hermite_conversion_contrib_struct[count];
			memcpy(list_right[n1][n2], currlist, count*sizeof(qmmm_hermite_conversion_contrib_struct));
			counters_right[n1][n2] = count;
			delete []currlist;
			delete []list;
		} // END FOR n1 n2

}


qmmm_hermite_conversion_info_struct::~qmmm_hermite_conversion_info_struct() {
	const int nmax = QMMM_HERMITE_CONVERSION_MAX_N;

	for(int n1 = 0; n1 <= nmax; n1++)
		for(int n2 = 0; n2 <= nmax; n2++) {
			delete []list_right[n1][n2];
		}
}


int qmmm_boysfunction_init(boys_func_interval_struct boys_list[QMMM_BOYS_N_MAX][QMMM_BOYS_NO_OF_INTERVALS]) {
	//  Util::TimeMeter timeMeter;
	qmmm_real halfstep, kfactorial, BoysFuncRawResult, Ak, midx;
	halfstep = (qmmm_real)QMMM_BOYS_X_MAX / QMMM_BOYS_NO_OF_INTERVALS * 0.5;
	for(int n = 0; n < QMMM_BOYS_N_MAX; n++) {
		for(int j = 0; j < QMMM_BOYS_NO_OF_INTERVALS; j++) {
			midx = (qmmm_real)QMMM_BOYS_X_MAX * j / QMMM_BOYS_NO_OF_INTERVALS + halfstep;
			boys_list[n][j].midx = midx;
			kfactorial = 1;
			int minusOneToPowk = 1;
			for(int k = 0; k < QMMM_BOYS_TAB_DEGREE; k++) {
				BoysFuncRawResult = boys_function_raw_simpson(n+k, midx);
				Ak = minusOneToPowk * BoysFuncRawResult / kfactorial;
				boys_list[n][j].a[k] = Ak;
				kfactorial *= k+1;
				minusOneToPowk *= -1;
			} /* END FOR k */
		} /* END FOR j */
	} /* END FOR n */
	//  global_Boys_init_flag = 1;
	//  timeMeter.print(LOG_AREA_INTEGRALS, "boysfunction_init");
	return 0;
}

