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

#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <cmath>
#include <vector>
#include <iomanip>
#include <omp.h>

#include "qmmm_utility.h"
#include "integrals_1el_potential_cuda_lib.h"

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

using namespace std;

struct __align__(16) cuda_atom_upper {
	qmmm_real charge;
	qmmm_real coords_x;
};

struct __align__(16) cuda_atom_lower {
	qmmm_real coords_y;
	qmmm_real coords_z;
};

struct cuda_atom {
	cuda_atom_upper upper;
	cuda_atom_lower lower;
};

struct cuda_boys_func_interval_struct {
	qmmm_real midx;
	double2 a[QMMM_BOYS_TAB_DEGREE / 2];
};

__constant__ int d_number_of_atoms;
__constant__ cuda_atom* d_cuda_atoms;
__constant__ cuda_boys_func_interval_struct* d_cuda_boys_list;

static bool isPow2(unsigned int x) {
	return ((x&(x-1))==0);
}

template <unsigned int blockSize, bool nIsPow2>
__global__ void reduce_kernel(double *d_in, double *d_out, unsigned int n) {
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
	unsigned int gridSize = blockSize*2*gridDim.x;

	extern __shared__ double smem[];

	double sum = 0.0f;

	while (i < n) {
		sum += d_in[i];
		if (nIsPow2 || i + blockSize < n) {
			sum += d_in[i+blockSize];
		}
		i += gridSize;
	}

	smem[tid] = sum;
	__syncthreads();


	if (blockSize >= 512) { if (tid < 256) { smem[tid] = sum = sum + smem[tid + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { smem[tid] = sum = sum + smem[tid + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (tid <  64) { smem[tid] = sum = sum + smem[tid +  64]; } __syncthreads(); }

	if (tid < 32) {
		volatile double* smemv = smem;
		if (blockSize >=  64) { smemv[tid] = sum = sum + smemv[tid + 32]; __syncthreads(); }
		if (blockSize >=  32) { smemv[tid] = sum = sum + smemv[tid + 16]; __syncthreads(); }
		if (blockSize >=  16) { smemv[tid] = sum = sum + smemv[tid +  8]; __syncthreads(); }
		if (blockSize >=   8) { smemv[tid] = sum = sum + smemv[tid +  4]; __syncthreads(); }
		if (blockSize >=   4) { smemv[tid] = sum = sum + smemv[tid +  2]; __syncthreads(); }
		if (blockSize >=   2) { smemv[tid] = sum = sum + smemv[tid +  1]; __syncthreads(); }
	}

	if (tid == 0) {
		d_out[blockIdx.x] = smem[0];
	}
}

// Wrapper function for kernel launch
void reduce(int size, int threads, int blocks, double *d_in, double *d_out, cudaStream_t stream) {

	int smemSize = (threads <= 32) ? 2 * threads * sizeof(double) : threads * sizeof(double);

	if (isPow2(size)) {
		switch (threads) {
		case 512:
			reduce_kernel<512, true><<< blocks, threads, smemSize, stream >>>(d_in, d_out, size); break;
		case 256:
			reduce_kernel<256, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 128:
			reduce_kernel<128, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 64:
			reduce_kernel<64, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 32:
			reduce_kernel<32, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 16:
			reduce_kernel<16, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  8:
			reduce_kernel<8, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  4:
			reduce_kernel<4, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  2:
			reduce_kernel<2, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  1:
			reduce_kernel<1, true><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		}
	}
	else {
		switch (threads) {
		case 512:
			reduce_kernel<512, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 256:
			reduce_kernel<256, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 128:
			reduce_kernel<128, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 64:
			reduce_kernel<64, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 32:
			reduce_kernel<32, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case 16:
			reduce_kernel<16, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  8:
			reduce_kernel<8, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  4:
			reduce_kernel<4, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  2:
			reduce_kernel<2, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		case  1:
			reduce_kernel<1, false><<< blocks, threads, smemSize,stream >>>(d_in, d_out, size); break;
		}
	}
}


/* To reduce number of register, center_arguments are passed into variables dx, dy, dz */
__global__ void compute_V_matrix_full_gpu_kernel_0(qmmm_real dx,
		qmmm_real dy,
		qmmm_real dz,
		const qmmm_real alpha,
		const qmmm_real boys_pre_factor,
		const qmmm_real factor,
		qmmm_real* d_result) {

	/* loop over atoms */
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	qmmm_real r_start;

	if(i < d_number_of_atoms) {
		cuda_atom atom = d_cuda_atoms[i];

		dx = atom.upper.coords_x - dx;
		dy = atom.lower.coords_y - dy;
		dz = atom.lower.coords_z - dz;

		const qmmm_real dxyz_squared = dx*dx + dy*dy + dz*dz;
		const qmmm_real arg = alpha * dxyz_squared;

		if(arg  >= QMMM_BOYS_X_MAX) {
			r_start = boys_pre_factor * rsqrt(arg);
		} else {
			/* choose which interval to use */
			const int intervalIndex = (int) arg * (qmmm_real) (QMMM_BOYS_NO_OF_INTERVALS / QMMM_BOYS_X_MAX);

			const cuda_boys_func_interval_struct* interval = &d_cuda_boys_list[intervalIndex];
			const qmmm_real deltax = arg - (qmmm_real)QMMM_BOYS_TEMP * (intervalIndex + 0.5);// arg - interval->midx;

			qmmm_real deltaxtopowk = 1;
			r_start = 0;

			for(int k = 0; k < QMMM_BOYS_TAB_DEGREE / 2; k++) {
				r_start += interval->a[k].x * deltaxtopowk;
				deltaxtopowk *= deltax;

				r_start += interval->a[k].y * deltaxtopowk;
				deltaxtopowk *= deltax;
			}
		}

		d_result[i] =  atom.upper.charge * factor * r_start;
	} // END IF ATOMS < NUMBER OF ATOMS
}

__global__ void compute_V_matrix_full_gpu_kernel(qmmm_real* d_result,
		const qmmm_real* d_work_list,
		int size) {

	/*  loop over atoms */
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	qmmm_real dx, dy, dz, dxyz_squared, alpha;

	qmmm_real arg, expMinusArg, r_start;

	int n_max;
	qmmm_real res = 0;
	qmmm_real tmp_xyz = 1;

	if(i < d_number_of_atoms) {
		cuda_atom atom = d_cuda_atoms[i];
		int j = 0;
		while (j < size) {
			switch ((int)d_work_list[j++]) {
			case NEW_CENTER :
				dx = atom.upper.coords_x - d_work_list[j++];
				dy = atom.lower.coords_y - d_work_list[j++];
				dz = atom.lower.coords_z - d_work_list[j++];
				dxyz_squared = dx*dx + dy*dy + dz*dz;
				break;
			case NEW_ALPHA :
				alpha = d_work_list[j++];
				arg = alpha * dxyz_squared;
				expMinusArg = std::exp(-arg);
				break;
			case  NEW_NMAX :
			{
				n_max = (int)d_work_list[j++];

				if(arg  >= QMMM_BOYS_X_MAX) {
					const qmmm_real boys_pre_factor = d_work_list[j];
					j++;
					qmmm_real arg_exp = 1;
					for(int k = 0; k < 2*n_max + 1; k++) {
						arg_exp *= arg;
					}
					r_start = boys_pre_factor * rsqrt(arg_exp);
				} else {
					j++; //skip boys_pre_factor

					/* choose which interval to use */
					const int intervalIndex = arg * (qmmm_real) (QMMM_BOYS_NO_OF_INTERVALS / QMMM_BOYS_X_MAX);

					const cuda_boys_func_interval_struct* interval = &d_cuda_boys_list[n_max*QMMM_BOYS_NO_OF_INTERVALS + intervalIndex];
					const qmmm_real deltax = arg - (qmmm_real)QMMM_BOYS_TEMP * (intervalIndex + 0.5);// arg - interval->midx;

					qmmm_real deltaxtopowk = 1;
					r_start = 0;

					for(int k = 0; k < QMMM_BOYS_TAB_DEGREE / 2; k++) {
						r_start += interval->a[k].x * deltaxtopowk;
						deltaxtopowk *= deltax;

						r_start += interval->a[k].y * deltaxtopowk;
						deltaxtopowk *= deltax;
					}
				}
				break;
			}
			case ADD :
			{
				const qmmm_real factor = d_work_list[j++];
				const int n = (int)d_work_list[j++];

				qmmm_real r = r_start;

				// get boys factor
				if(n_max > 0) {
					for(int k = n_max-1; k >= n; k--)
						r = (2*arg*r + expMinusArg) / (2*k+1);
				}

				res += tmp_xyz * factor * r;
				tmp_xyz = 1;
				break;
			}
			case DX :
				tmp_xyz *= dx;
				break;
			case DY :
				tmp_xyz *= dy;
				break;
			case DZ :
				tmp_xyz *= dz;
				break;
			}
		}
		d_result[i] =  atom.upper.charge * res;
	} // END IF ATOMS < NUMBER OF ATOMS
}

class monomial_help_class {
public:
	qmmm_real factor;
	int length;
	int x;
	int y;
	int z;

	monomial_help_class(qmmm_real factor, int length, int x, int y, int z)
	:factor(factor), length(length), x(x), y(y), z(z)
	{};

	std::vector<monomial_help_class> operator*(std::vector<monomial_help_class> rhs);
};

std::vector<monomial_help_class> monomial_help_class::operator*(std::vector<monomial_help_class>  rhs) {
	std::vector<monomial_help_class> m;

	for(unsigned int i = 0; i < rhs.size(); i++) {
		m.push_back(monomial_help_class(factor * rhs[i].factor, length + rhs[i].length, x + rhs[i].x, y + rhs[i].y, z + rhs[i].z));
	}
	return m;
}


static unsigned int nextPow2( unsigned int x ) {
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

static vector<int> get_compatible_devices() {

	unsigned int cuda_success = 0;

	int number_of_GPUs;
	cuda_success |= (unsigned int)(cudaGetDeviceCount(&number_of_GPUs));
	cout << "Found " << number_of_GPUs << " GPUs (cuda_success = " << cuda_success << ")" << endl;
	vector<int> devices;
	for(int i = 0; i < number_of_GPUs; i++) {
		int dev_ID;
		cudaDeviceProp props;
		cuda_success |= (unsigned int) cudaGetDevice(&dev_ID);
		cuda_success |= (unsigned int) cudaGetDeviceProperties(&props, i);
		printf("Device %d: \"%s\" with Compute %d.%d capability found\n", dev_ID, props.name, props.major, props.minor);
		//if(props.major == 2 && props.minor == 0) {
		if(props.major >= 2) {
			devices.push_back(i);
		}
	}

	if(cuda_success != 0) {
		cout << "[Function get_compatible_devices] CUDA error = " << cuda_success << endl;
	}
	return devices;
}


static void get_num_blocks_and_threads(int n, int maxBlocks, int maxThreads, int &blocks, int &threads) {
	threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
	blocks = (n + (threads * 2 - 1)) / (threads * 2);
	blocks = MIN(maxBlocks, blocks);
}

static vector<monomial_help_class>  genereate_polynom(int n, const int x, const int y, const int z) {

	vector<monomial_help_class> monomial_vector;
	if(x == 0 && y == 0 && z == 0) {
		monomial_help_class m(1, n, 0, 0, 0);
		vector<monomial_help_class> m_v;
		m_v.push_back(m);
		return m_v;
	}

	if(x > 0) {
		monomial_vector = monomial_help_class(1, 0, 1, 0, 0) * genereate_polynom(n +1, x-1 ,y, z);
		if(x > 1) {
			vector<monomial_help_class> vec_tmp = monomial_help_class(x-1, 0, 0, 0, 0) * genereate_polynom(n +1, x-2, y, z);
			monomial_vector.insert(monomial_vector.end(), vec_tmp.begin(), vec_tmp.end());
		}
	} else if(y > 0) {
		monomial_vector = monomial_help_class(1, 0, 0, 1, 0) * genereate_polynom(n +1, x, y-1, z);
		if(y > 1) {
			vector<monomial_help_class> vec_tmp = monomial_help_class(y-1, 0, 0, 0, 0) * genereate_polynom(n +1, x, y-2, z);
			monomial_vector.insert(monomial_vector.end(), vec_tmp.begin(), vec_tmp.end());
		}
	} else if(z > 0) {
		monomial_vector = monomial_help_class(1, 0, 0, 0, 1) * genereate_polynom(n +1, x, y, z-1);
		if(z > 1) {
			vector<monomial_help_class> vec_tmp = monomial_help_class(z-1, 0, 0, 0, 0) * genereate_polynom(n +1, x, y, z-2);
			monomial_vector.insert(monomial_vector.end(), vec_tmp.begin(), vec_tmp.end());
		}
	}
	return monomial_vector;
}

static int clean_work_list(vector <double> &list) {
	/* Remove ADD with factor == 0 (will result in multiply with 0 => result = 0 (this is done in the collection code)
	 * Remove NEW_NMAX if new NEW_MAX is found without an ADD between
	 * Remove NEW ALPHA if new NEW_ALPHA is found without an NEW_NMAX is found between
	 * Remove NEW_CENTERCOORDS if NEW_CENTERCOORDS is found without an NEW_NMAX is found between
	 * */
	int old_alpha_index = -1;
	int old_n_max_index = -1;
	bool n_max_found = false;
	bool add_found = true;
	bool center_found = true;

	int i = 0;
	int size = list.size();

	//Clean up N_MAX
	while (i < size) {
		switch ((int)list[i]) {
		case NEW_CENTER :
			i +=  4;
			break;
		case NEW_ALPHA :
			i += 2;
			break;
		case  NEW_NMAX :
			if(!add_found) {
				list.erase(list.begin() + old_n_max_index, list.begin() + old_n_max_index + 3);
				old_n_max_index = i-3;
				size = list.size();
			} else {
				old_n_max_index = i;
				i += 3;
			}
			add_found = false;
			break;
		case ADD :
			add_found = true;
			i += 3;
			break;
		default :
			i++;
			break;
		}
	}

	i = 0;
	size = list.size();
	n_max_found = true;
	//Clean up ALPHA
	while (i < size) {
		switch ((int)list[i]) {
		case NEW_CENTER :
			center_found = true;
			i += 4;
			break;
		case NEW_ALPHA :
			if(!n_max_found && !center_found) {
				cout << "Erasing an alpha" << endl;
				list.erase(list.begin() + old_alpha_index, list.begin() + old_alpha_index + 2);
				old_alpha_index = i-2;
				size = list.size();
			} else {
				old_alpha_index = i;
				i += 2;
			}
			center_found = false;
			n_max_found = false;
			break;
		case  NEW_NMAX :
			n_max_found = true;
			i += 3;
			break;
		case ADD :
			i += 3;
			break;
		default :
			i++;
			break;
		}
	}
	return 0;
}

static bool is_case_0(vector <double> &list) {
	bool case_0 = false;
	int i = 0;
	int size = list.size();
	int number_of_new_center = 0;
	int number_of_new_alpha = 0;
	int number_of_new_n_max = 0;
	int number_of_new_add = 0;
	int number_of_dx = 0;
	int sum_n_max = 0;
	int sum_length = 0;

	//Clean up N_MAX
	while (i < size) {
		switch ((int)list[i]) {
		case NEW_CENTER :
			number_of_new_center++;
			i = i + 4;
			break;
		case NEW_ALPHA :
			number_of_new_alpha++;
			i = i +2;
			break;
		case  NEW_NMAX :
			number_of_new_n_max++;
			sum_n_max += (int)list[i+1];
			i += 3;
			break;
		case ADD :
			number_of_new_add++;
			sum_length += (int)list[i+2];
			i += 3;
			break;
		default :
			number_of_dx++;
			i++;
			break;
		}
	}
	if(number_of_new_center == 1 && number_of_new_alpha == 1 && number_of_new_n_max == 1 && number_of_new_add == 1 && number_of_dx == 0 && sum_n_max == 0 && sum_length == 0) {
		case_0 = true;
	}
	return case_0;
}

// Library version (library input, cleaned version)
int compute_V_matrix_full_gpu_lib(const qmmm_basis_info_struct& basis_info,
		const int num_qm_atoms,
		const qmmm_atom* qm_atom_list,
		const int num_mm_atoms,
		const qmmm_atom* mm_atom_list,
		const qmmm_real threshold,
		qmmm_real* result) {

	const int number_of_basis_functions = basis_info.no_of_basis_funcs;
	const int iterations = (number_of_basis_functions * (number_of_basis_functions + 1)) / 2;
	int cntr[20];

	for(int i = 0; i < 20; i++) {
		cntr[i] = 0;
	}

	std::cout << "Enter compute_V_matrix_full_gpu" << std::endl;
	std::cout << "Number of basis functions = " << number_of_basis_functions << ",  number of qm atoms = " << num_qm_atoms << ", number of mm atoms " <<  num_mm_atoms << std::endl;
//	std::cout << "Number of iterations = " << iterations << std::endl;

#ifdef _OPENMP
	double time_start, time_stop;
	time_start = omp_get_wtime();
#endif

	const qmmm_monomial_info_struct* monomial_info = new qmmm_monomial_info_struct();
	const qmmm_hermite_conversion_info_struct* hermite_conversion_info =  new qmmm_hermite_conversion_info_struct();

	boys_func_interval_struct boys_list[QMMM_BOYS_N_MAX][QMMM_BOYS_NO_OF_INTERVALS];
	qmmm_boysfunction_init(boys_list);

	// We use "cuda_atom" to be able to do 128bit copy from global memory to each thread on the GPU
	cuda_atom* cuda_atoms = (cuda_atom*) malloc(sizeof(cuda_atom) * num_mm_atoms);
	for(int i = 0; i < num_mm_atoms; i++) {
		cuda_atoms[i].upper.charge = mm_atom_list[i].charge;
		cuda_atoms[i].upper.coords_x = mm_atom_list[i].coords[0];
		cuda_atoms[i].lower.coords_y = mm_atom_list[i].coords[1];
		cuda_atoms[i].lower.coords_z = mm_atom_list[i].coords[2];
	}

	// We use "cuda_boys_func_interval_struct" to be able to do 128bit copy from global memory to each thread on the GPU
	cuda_boys_func_interval_struct * cuda_boys_list = (cuda_boys_func_interval_struct*) malloc(sizeof(cuda_boys_func_interval_struct) * QMMM_BOYS_N_MAX * QMMM_BOYS_NO_OF_INTERVALS);
	for(int i = 0; i < QMMM_BOYS_N_MAX; i ++) {
		for(int j = 0; j < QMMM_BOYS_NO_OF_INTERVALS; j++) {
			cuda_boys_list[i*QMMM_BOYS_NO_OF_INTERVALS + j].midx  = boys_list[i][j].midx;
			for(int k = 0; k < QMMM_BOYS_TAB_DEGREE / 2; k++) {
				cuda_boys_list[i*QMMM_BOYS_NO_OF_INTERVALS + j].a[k].x = boys_list[i][j].a[2*k];
				cuda_boys_list[i*QMMM_BOYS_NO_OF_INTERVALS + j].a[k].y = boys_list[i][j].a[2*k+1];
			}
		}
	}

	unsigned int cuda_success = 0;

	vector<int> devices = get_compatible_devices();

	if(devices.size() == 0) {
		cout << "Error! No devices  with compute capabilyt >= 2 found!\n" << endl;
		free(cuda_atoms);
		free(cuda_boys_list);
		return 0;
	}

	const int nr_streams = 4;

#ifdef _OPENMP

#pragma omp parallel reduction(+:cuda_success) num_threads(devices.size())
	{
		int thread_id = omp_get_thread_num();
		cudaSetDevice(devices[thread_id]);
#endif

		int dev_ID;
		cuda_success |= (unsigned int) cudaGetDevice(&dev_ID);

		cudaSetDeviceFlags(cudaDeviceScheduleSpin);
		cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

		std::vector<qmmm_distribution_spec_struct*> psi_list;

		qmmm_real* h_result;
		h_result = (qmmm_real*)malloc(sizeof(qmmm_real)*iterations);

		int estimated_max_size_work_list = 100;
		qmmm_real* d_work_list;
		qmmm_real* h_work_list;

		//		std::cout << "Allocating "<< sizeof(qmmm_real) * estimated_max_size_work_list*nr_streams / (1024*1024) <<" MB (d_work_list) on the GPU" << std::endl;
		cuda_success |= (unsigned int) cudaMalloc(&d_work_list, sizeof(qmmm_real) * estimated_max_size_work_list*nr_streams);
		cuda_success |= (unsigned int) cudaMallocHost(&h_work_list,sizeof(qmmm_real)*estimated_max_size_work_list*nr_streams);

		//		std::cout << "Allocating "<< sizeof(cuda_boys_func_interval_struct) * BOYS_N_MAX * BOYS_NO_OF_INTERVALS / (1024*1024) <<" MB (d_cuda_boys_list) on the GPU" << std::endl;
		cuda_boys_func_interval_struct* d_cuda_boys_list;
		cuda_success |= (unsigned int) cudaMalloc(&d_cuda_boys_list, sizeof(cuda_boys_func_interval_struct) * QMMM_BOYS_N_MAX * QMMM_BOYS_NO_OF_INTERVALS);
		cuda_success |= (unsigned int) cudaMemcpy(d_cuda_boys_list, cuda_boys_list, sizeof(cuda_boys_func_interval_struct) * QMMM_BOYS_N_MAX * QMMM_BOYS_NO_OF_INTERVALS, cudaMemcpyHostToDevice);

		//		std::cout << "Allocating "<< sizeof(cuda_atom) * number_of_atoms / (1024*1024) <<" MB (d_cuda_atoms) on the GPU" << std::endl;
		cuda_atom* d_cuda_atoms;
		cuda_success |= (unsigned int) cudaMalloc(&d_cuda_atoms, sizeof(cuda_atom) * num_mm_atoms);
		cuda_success |= (unsigned int) cudaMemcpy(d_cuda_atoms, cuda_atoms, sizeof(cuda_atom) * num_mm_atoms, cudaMemcpyHostToDevice);

		qmmm_real* partial_result;
		partial_result = (qmmm_real*) malloc (sizeof(qmmm_real) * num_mm_atoms);

		//		std::cout << "Allocating "<< sizeof(qmmm_real) * number_of_atoms / (1024*1024) <<" MB (d_partial_result) on the GPU" << std::endl;
		qmmm_real* d_partial_result;
		cuda_success |= (unsigned int) cudaMalloc(&d_partial_result, sizeof(qmmm_real) * num_mm_atoms);

		//		std::cout << "Allocating "<< sizeof(qmmm_real) * iterations / (1024*1024) << " MB (d_result) on the GPU" << std::endl;
		qmmm_real* d_result;
		cuda_success |= (unsigned int) cudaMalloc(&d_result, sizeof(qmmm_real) * iterations);
		cuda_success |= cudaMemset(d_result,0,iterations*sizeof(qmmm_real));

		//create event for sync.
		int current_stream = 0;
		cudaEvent_t done[nr_streams];
		cudaStream_t stream[nr_streams];
		for(int i=0; i < nr_streams; i++)
			cuda_success |= cudaEventCreate(&done[i]);

		for (int i = 0; i < nr_streams; ++i)
			cudaStreamCreate(&stream[i]);

		const int threadsPerBlock = 64;//256;
		const int blocksPerGrid = (num_mm_atoms + threadsPerBlock -1) / threadsPerBlock;

		int threadsInReduction, threadsInReduction2;
		int blocksInReduction, blocksInReduction2;
		const int maxThreads = 512;
		const int maxBlocks = 65535;

		get_num_blocks_and_threads(num_mm_atoms, maxBlocks, maxThreads, blocksInReduction, threadsInReduction);
		get_num_blocks_and_threads(blocksInReduction, maxBlocks, maxThreads, blocksInReduction2, threadsInReduction2);

		qmmm_real* partial_reduce_result = (qmmm_real*) malloc(sizeof(qmmm_real)*blocksInReduction);

		//		cout << "blocksInReduction = " << blocksInReduction << ", threadsInReduction = " << threadsInReduction << endl;
		//		cout << "blocksInReduction2 = " << blocksInReduction2 << ", threadsInReduction2 = " << threadsInReduction2 << endl;

		cuda_success |= cudaMemcpyToSymbol("d_number_of_atoms", &num_mm_atoms, sizeof(int));
		cuda_success |= cudaMemcpyToSymbol("d_cuda_atoms", &d_cuda_atoms, sizeof(cuda_atom*));
		cuda_success |= cudaMemcpyToSymbol("d_cuda_boys_list", &d_cuda_boys_list, sizeof(cuda_boys_func_interval_struct*));

		//trigger event
		for(int i=0; i < nr_streams; i++) {
			cudaEventRecord(done[i],stream[i]);
		}

		vector<int> iteration_indexes;
		/* This loop replaces mu and nu */
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 800)
#endif
		for(int i = 0; i < iterations; i++) {

			const int mu = sqrt(2*i + 0.25) - 0.5;
			const int nu = i - (mu * (mu + 1)) / 2;
			const int result_index = mu*number_of_basis_functions + nu;

			iteration_indexes.push_back(i);

			//Get first contracted gaussian
			int start_prim_mu = basis_info.start_index[mu];
			int stop_prim_mu = basis_info.start_index[mu+1];

			//Get second contracted gaussian
			int start_prim_nu = basis_info.start_index[nu];
			int stop_prim_nu = basis_info.start_index[nu+1];

			/* compute matrix element [mu,nu] */
			qmmm_real old_alpha = -1;
			qmmm_real old_coord_x = -1;
			qmmm_real old_coord_y = -1;
			qmmm_real old_coord_z = -1;

			qmmm_real* center_coords;
			qmmm_real alpha;

			vector<qmmm_real> work_list;
			int old_n_max = -1;

			for(int j = start_prim_mu; j < stop_prim_mu; j++) {
				for(int k = start_prim_nu; k < stop_prim_nu; k++) {
					psi_list.clear();
					get_product_simple_prims(basis_info.simple_primitive_list[j], basis_info.simple_primitive_list[k], psi_list, threshold);

					if(psi_list.size() > 0) {
						// Alpha and centerCoords is same for all new prims
						center_coords = psi_list[0]->center_coords;
						alpha = psi_list[0]->exponent;
						const qmmm_real inv_alpha = 1 / alpha;
						const qmmm_real resultPreFactor = 2 * pi * inv_alpha;

						//Is center_corrds and alpha from last j, k same as this?
						if(center_coords[0] != old_coord_x || center_coords[1] != old_coord_y || center_coords[2] != old_coord_z) {
							work_list.push_back(NEW_CENTER);
							work_list.push_back(center_coords[0]);
							work_list.push_back(center_coords[1]);
							work_list.push_back(center_coords[2]);
							old_coord_x = center_coords[0];
							old_coord_y = center_coords[1];
							old_coord_z = center_coords[2];
							old_n_max = -1;
							old_alpha = -1; // need to add alpha after new center coords so that arg gets updated correctly
						}

						if(alpha != old_alpha) {
							work_list.push_back((NEW_ALPHA));
							work_list.push_back(alpha);
							old_alpha = alpha;
							old_n_max = -1;
						}

						for(int m = 0; m < psi_list.size(); m++) {
							const qmmm_distribution_spec_struct* psi = psi_list[m];

							const int n1x = psi->monomialInts[0];
							const int n1y = psi->monomialInts[1];
							const int n1z = psi->monomialInts[2];

							const int n1_max = n1x + n1y + n1z;
							const int n2_max = 0;
							const int n_max = n1_max + n2_max;

							if(n_max >= QMMM_BOYS_N_MAX) {
								cout << "Error, n_max >= BOYS_N_MAX!" << endl;
								exit(0);
							}

							if(old_n_max != n_max) {
								const qmmm_real boys_pre_factor_value = qmmm_boys_pre_factor(n_max);
								work_list.push_back(NEW_NMAX);
								work_list.push_back(n_max);
								work_list.push_back(boys_pre_factor_value);
								old_n_max = n_max;
							}

							const int monomialIndex = monomial_info->monomial_index_list[n1x][n1y][n1z];
							const int no_of_monomials_n1 = monomial_info->no_of_monomials_list[n1_max];
							const int no_of_monomials_n2 = monomial_info->no_of_monomials_list[n2_max];

							const int no_of_contribs = hermite_conversion_info->counters_right[n1_max][n2_max];
							const qmmm_hermite_conversion_contrib_struct* list = hermite_conversion_info->list_right[n1_max][n2_max];

							qmmm_real pre_factor;
							int x, y, z, sum_n, i1, i2;
							for(int n=0; n < no_of_contribs; n++) {
								if(list[n].destIndex == monomialIndex) {
									if(no_of_monomials_n1 > 0 && no_of_monomials_n2 > 0) {

										i1 = list[n].sourceIndex / no_of_monomials_n2;
										i2 = list[n].sourceIndex % no_of_monomials_n2;

										x = monomial_info->monomial_list[i1].ix;
										y = monomial_info->monomial_list[i1].iy;
										z = monomial_info->monomial_list[i1].iz;
										sum_n = x+y+z;

										x += monomial_info->monomial_list[i2].ix;
										y += monomial_info->monomial_list[i2].iy;
										z += monomial_info->monomial_list[i2].iz;

										if(sum_n % 2 == 1) {
											pre_factor = resultPreFactor * list[n].coeff * std::pow(inv_alpha, -list[n].a_power)*psi->coeff;
										} else {
											pre_factor = -resultPreFactor * list[n].coeff * std::pow(inv_alpha, -list[n].a_power)*psi->coeff;
										}

										vector<monomial_help_class> monomial_vector_tmp = genereate_polynom(0, x, y, z);

										for(int z = 0; z < monomial_vector_tmp.size(); z++) {
											const qmmm_real boys_factor = std::pow(-2*alpha, monomial_vector_tmp[z].length);
											monomial_vector_tmp[z].factor *= pre_factor * boys_factor;
											//calculation is multiplied with factor, if factor == 0 => then result == 0 as well
											if(fabs(monomial_vector_tmp[z].factor) != 0) {
												for(int w = 0; w < monomial_vector_tmp[z].x; w++) {
													work_list.push_back(DX);
												}
												for(int w = 0; w < monomial_vector_tmp[z].y; w++) {
													work_list.push_back(DY);
												}
												for(int w = 0; w < monomial_vector_tmp[z].z; w++) {
													work_list.push_back(DZ);
												}
												work_list.push_back(ADD); // end of monomial, use add
												work_list.push_back(monomial_vector_tmp[z].factor);
												work_list.push_back((qmmm_real)monomial_vector_tmp[z].length);
											}
										}
									}
								}
							} /* END FOR NEW PRIMS */
						}// end if number_of_new_prims > 0
					} /* END FOR NEW PRIMS */
				} /* END FOR k */
			} /* END FOR j */


			/* Anything to calculate for this index? */
			if(work_list.size() > 0) {

				clean_work_list(work_list);

				cudaEventSynchronize(done[current_stream]);

				if(work_list.size() > estimated_max_size_work_list){
					estimated_max_size_work_list = work_list.size();
					cuda_success |= (unsigned int) cudaFree(d_work_list);
					cuda_success |= (unsigned int) cudaMalloc(&d_work_list, sizeof(qmmm_real) * estimated_max_size_work_list*nr_streams);
					cuda_success |= (unsigned int) cudaFreeHost(h_work_list);
					cuda_success |= (unsigned int) cudaMallocHost(&h_work_list, sizeof(qmmm_real)*estimated_max_size_work_list*nr_streams);
				}

				if(is_case_0(work_list)) {

					compute_V_matrix_full_gpu_kernel_0<<<blocksPerGrid, threadsPerBlock, 0, stream[current_stream]>>>(work_list[1], // center coords
							work_list[2],
							work_list[3],
							work_list[5],  // alpha
							work_list[8],  // boys_pre_factor
							work_list[10], // factor
							d_partial_result);
				} else {
					memcpy(h_work_list+current_stream*estimated_max_size_work_list, (void*)&work_list[0],sizeof(qmmm_real)*work_list.size());
					cuda_success |= (unsigned int) cudaMemcpyAsync(d_work_list+current_stream*estimated_max_size_work_list, h_work_list+current_stream*estimated_max_size_work_list, sizeof(qmmm_real) * work_list.size(), cudaMemcpyHostToDevice,stream[current_stream]);

					cudaEventRecord(done[current_stream],stream[current_stream]);

					compute_V_matrix_full_gpu_kernel<<<blocksPerGrid, threadsPerBlock, 0, stream[current_stream]>>>(d_partial_result,
							d_work_list+current_stream*estimated_max_size_work_list,
							work_list.size());
				}

				if(blocksInReduction == 1){
					reduce(num_mm_atoms, threadsInReduction, blocksInReduction, d_partial_result, d_result+i,stream[current_stream]);
				} else {
					reduce(num_mm_atoms, threadsInReduction, blocksInReduction, d_partial_result, d_partial_result,stream[current_stream]);
					reduce(blocksInReduction,threadsInReduction2,blocksInReduction2, d_partial_result,d_result+i,stream[current_stream]);
				}

				current_stream = (current_stream + 1)%nr_streams;
			}
		} // END ITERATIONS (mu & nu), END PARALLEL LOOP

		//TODO Only memcpy the number of iterations this thread has done
		cuda_success |= (unsigned int) cudaMemcpy(h_result, d_result, sizeof(qmmm_real) * iterations, cudaMemcpyDeviceToHost);

		cout << "Device " << dev_ID << " did " << iteration_indexes.size() << " iterations ("<< 100 * iteration_indexes.size() / iterations << "%)"<< endl;

		for(int j=0; j < iteration_indexes.size(); j++) {
			const int i = iteration_indexes[j];
			const int mu = sqrt(2*i + 0.25) - 0.5;
			const int nu = i - (mu * (mu + 1)) / 2;
			const int result_index = mu*number_of_basis_functions + nu;
			result[result_index] = h_result[i];
		}

		cuda_success |= (unsigned int) cudaFree(d_cuda_atoms);
		cuda_success |= (unsigned int) cudaFree(d_work_list);
		cuda_success |= (unsigned int) cudaFree(d_partial_result);
		cuda_success |= (unsigned int) cudaFree(d_cuda_boys_list);
		cuda_success |= (unsigned int) cudaFree(d_result);
		cuda_success |= (unsigned int) cudaFreeHost(h_work_list);
		for(int i=0; i < nr_streams; i++)
			cuda_success |= (unsigned int) cudaEventDestroy(done[i]);
		for (int i = 0; i < nr_streams; ++i)
			cudaStreamDestroy(stream[i]);
		cuda_success |= (unsigned int) cudaThreadExit();

		if(cuda_success != 0) {
			cout << "[Device " << dev_ID << "] CUDA error = " << cuda_success << endl;
		}
		free(partial_result);
		free(partial_reduce_result);
		free(h_result);
#ifdef _OPENMP
	} // end parallel section
#endif


	// copy values to the other triangle
	for(int mu = 0; mu < number_of_basis_functions; mu++)
		for(int nu = mu+1; nu < number_of_basis_functions; nu++)
			result[mu*number_of_basis_functions+nu] = result[nu*number_of_basis_functions+mu];

#ifdef _OPENMP
	time_stop = omp_get_wtime();
	std::cout << "Execution time of compute_V_matrix_gpu = " << time_stop - time_start << std::endl;
#endif

	free(cuda_atoms);
	free(cuda_boys_list);
	delete hermite_conversion_info;

	//std::cout <<"Exit compute_V_matrix_gpu" << std::endl;

	return 0;
}
