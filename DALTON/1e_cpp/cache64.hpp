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

#ifndef __CACHE_64__
#define __CACHE_64__

#include <iostream>
#include <pmmintrin.h>
#include <xmmintrin.h>
#include <immintrin.h>
#include "cache.hpp"


//#define __SSE2__
//#define __AVX__

//#ifdef __AVX__
// radovan: deactivated this section
//          this crashes with GNU compilers due to invalid mem access
#ifdef __AVX__UNDEFINED

#define DOUBLES_PER_AVX 4

class cacheline64 {
  public:
    __m256d d4[AVXD_PER_CACHE_LINE];

    inline __m256d & operator[](int i) {
        return d4[i];
    }

    inline const __m256d & operator[](int i) const {
        return d4[i];
    }

    inline double & operator()(int i) {
        return ((double*)d4)[i];
    }

    inline const double & operator()(int i) const {
        return ((double*)d4)[i];
    }

    inline cacheline64 & operator=(const double & rhs) {
        d4[0] = _mm256_set1_pd(rhs);
        d4[1] = _mm256_set1_pd(rhs);

        return *this;
    }

    inline cacheline64 & operator+=(const cacheline64 & rhs) {
        for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
            d4[i] = _mm256_add_pd(d4[i], rhs.d4[i]);
        return *this;
    }

    inline cacheline64 & operator-=(const cacheline64 & rhs) {
        for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
            d4[i] = _mm256_sub_pd(d4[i], rhs.d4[i]);
        return *this;
    }

    inline cacheline64 & operator*=(const cacheline64 & rhs) {
        for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
            d4[i] = _mm256_mul_pd(d4[i], rhs.d4[i]);
        return *this;
    }

    inline cacheline64 & operator*=(const __m256d & rhs) {
        for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
            d4[i] = _mm256_mul_pd(d4[i], rhs);
        return *this;
    }

    inline cacheline64 & operator*=(double rhs) {
        __m256d v = _mm256_set1_pd(rhs);
        for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
            d4[i] = _mm256_mul_pd(d4[i], v);
        return *this;
    }


    inline cacheline64 operator+(const cacheline64 & rhs) const {
        cacheline64 ret = *this;
        ret += rhs;
        return ret;
    }

    inline cacheline64 operator-(const cacheline64 & rhs) const {
        cacheline64 ret = *this;
        ret -= rhs;
        return ret;
    }

    inline cacheline64 operator*(const cacheline64 & rhs) const {
        cacheline64 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline cacheline64 operator*(const __m256d & rhs) const {
        cacheline64 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline cacheline64 operator*(double rhs) const {
        cacheline64 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline void set (double d) {
        d4[0] = _mm256_set1_pd(d);
        d4[1] = _mm256_set1_pd(d);
    }

    inline void set (const double & a0, const double & a1, const double & a2, const double & a3, const double & a4, const double & a5, const double & a6, const double & a7) {
        d4[0] = _mm256_set_pd(a3,a2,a1,a0);
        d4[1] = _mm256_set_pd(a7,a6,a5,a4);
    }

}  __attribute__((aligned(CACHE_LINE_SIZE)));


static inline void store(double * p, const cacheline64 & rhs) {
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        _mm256_store_pd(p+DOUBLES_PER_AVX*i  , rhs.d4[i]);
}

static inline void store(cacheline64 * p, const cacheline64 & rhs) {
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        _mm256_store_pd((double*)p+DOUBLES_PER_AVX*i  , rhs.d4[i]);
}

static inline cacheline64 load(const double * p) {
    cacheline64 ret;
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        ret.d4[i] = _mm256_load_pd(p+DOUBLES_PER_AVX*i);
    return ret;
}

static inline cacheline64 load(const cacheline64 * p) {
    cacheline64 ret;
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        ret.d4[i] = _mm256_load_pd((double*)p+DOUBLES_PER_AVX*i);
    return ret;
}


inline cacheline64 sqrt(const cacheline64 & rhs) {
    cacheline64 ret;
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        ret.d4[i] = _mm256_sqrt_pd(rhs.d4[i]);
    return ret;
}

inline cacheline64 inv(const cacheline64 & rhs) {
    cacheline64 ret; __m256d one  = _mm256_set1_pd(1.);
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        ret.d4[i] = _mm256_div_pd(one, rhs.d4[i]);
    return ret;
}

//absolute value using SSE intrinsics
inline cacheline64 abs(const cacheline64 & rhs) {
    static const __m256d sign_mask = _mm256_set1_pd(-0.); // -0. = 1 << 63
    cacheline64 ret;
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        ret.d4[i] = _mm256_andnot_pd(sign_mask, rhs.d4[i]);
    return ret;
}

inline cacheline64 max(const cacheline64 & lhs, const cacheline64 & rhs) {
    cacheline64 ret;
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        ret.d4[i] = _mm256_max_pd(lhs.d4[i], rhs.d4[i]);
    return ret;
}

inline cacheline64 min(const cacheline64 & lhs, const cacheline64 & rhs) {
    cacheline64 ret;
    for (int i=0; i<AVXD_PER_CACHE_LINE; ++i)
        ret.d4[i] = _mm256_min_pd(lhs.d4[i], rhs.d4[i]);
    return ret;
}

#elif defined(__SSE2__)

class cacheline64 {
  public:
    __m128d d2[MM128_PER_CACHE_LINE];


    inline __m128d & operator[](int i) {
        return d2[i];
    }

    inline const __m128d & operator[](int i) const {
        return d2[i];
    }

    inline double & operator()(int i) {
        return ((double*)d2)[i];
    }

    inline const double & operator()(int i) const {
        return ((double*)d2)[i];
    }

    inline cacheline64 & operator=(const double & rhs) {
        d2[0] = _mm_load1_pd(&rhs);
        d2[1] = _mm_load1_pd(&rhs);
        d2[2] = _mm_load1_pd(&rhs);
        d2[3] = _mm_load1_pd(&rhs);

        return *this;
    }

    inline cacheline64 & operator+=(const cacheline64 & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_add_pd(d2[i], rhs.d2[i]);
        return *this;
    }

    inline cacheline64 & operator-=(const cacheline64 & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_sub_pd(d2[i], rhs.d2[i]);
        return *this;
    }

    inline cacheline64 & operator*=(const cacheline64 & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_mul_pd(d2[i], rhs.d2[i]);
        return *this;
    }

    inline cacheline64 & operator*=(const __m128d & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_mul_pd(d2[i], rhs);
        return *this;
    }

    inline cacheline64 & operator*=(double rhs) {
        __m128d v = _mm_load1_pd(&rhs);
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_mul_pd(d2[i], v);
        return *this;
    }


    inline cacheline64 operator+(const cacheline64 & rhs) const {
        cacheline64 ret = *this;
        ret += rhs;
        return ret;
    }

    inline cacheline64 operator-(const cacheline64 & rhs) const {
        cacheline64 ret = *this;
        ret -= rhs;
        return ret;
    }

    inline cacheline64 operator*(const cacheline64 & rhs) const {
        cacheline64 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline cacheline64 operator*(const __m128d & rhs) const {
        cacheline64 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline cacheline64 operator*(double rhs) const {
        cacheline64 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline void set (double d) {
        d2[0] = _mm_set1_pd(d);
        d2[1] = _mm_set1_pd(d);
        d2[2] = _mm_set1_pd(d);
        d2[3] = _mm_set1_pd(d);
    }

    inline void set (const double & a0, const double & a1, const double & a2, const double & a3, const double & a4, const double & a5, const double & a6, const double & a7) {
        d2[0] = _mm_loadl_pd(d2[0], &a0);
        d2[0] = _mm_loadh_pd(d2[0], &a1);
        d2[1] = _mm_loadl_pd(d2[1], &a2);
        d2[1] = _mm_loadh_pd(d2[1], &a3);
        d2[2] = _mm_loadl_pd(d2[2], &a4);
        d2[2] = _mm_loadh_pd(d2[2], &a5);
        d2[3] = _mm_loadl_pd(d2[3], &a6);
        d2[3] = _mm_loadh_pd(d2[3], &a7);
    }

}  __attribute__((aligned(CACHE_LINE_SIZE)));



static inline void store(double * p, const cacheline64 & rhs) {
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        _mm_store_pd(p+2*i  , rhs.d2[i]);
}

static inline void store(cacheline64 * p, const cacheline64 & rhs) {
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        _mm_store_pd((double*)p+2*i  , rhs.d2[i]);
}

static inline cacheline64 load(const double * p) {
    cacheline64 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_load_pd(p+2*i);
    return ret;
}

static inline cacheline64 load(const cacheline64 * p) {
    cacheline64 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_load_pd((double*)p+2*i);
    return ret;
}


inline cacheline64 sqrt(const cacheline64 & rhs) {
    cacheline64 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_sqrt_pd(rhs.d2[i]);
    return ret;
}

inline cacheline64 inv(const cacheline64 & rhs) {
    cacheline64 ret; __m128d one  = _mm_set1_pd(1.);
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_div_pd(one, rhs.d2[i]);
    return ret;
}

//absolute value using SSE intrinsics
inline cacheline64 abs(const cacheline64 & rhs) {
    static const __m128d sign_mask = _mm_set1_pd(-0.); // -0. = 1 << 63
    cacheline64 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_andnot_pd(sign_mask, rhs.d2[i]);
    return ret;
}

inline cacheline64 max(const cacheline64 & lhs, const cacheline64 & rhs) {
    cacheline64 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_max_pd(lhs.d2[i], rhs.d2[i]);
    return ret;
}

inline cacheline64 min(const cacheline64 & lhs, const cacheline64 & rhs) {
    cacheline64 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_min_pd(lhs.d2[i], rhs.d2[i]);
    return ret;
}
#endif



void PackArrays     (const double * __restrict__  arrays, unsigned long int ArraySize, cacheline64 * __restrict__  array8);
void UnPackArrays   (const cacheline64 * __restrict__ array8,  unsigned long int ArraySize, double * arrays);
std::ostream & operator<<(std::ostream & os, const cacheline64 & rhs);

#endif
