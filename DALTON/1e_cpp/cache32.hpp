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


#ifndef __CACHE_32__
#define __CACHE_32__

#include <iostream>
#include <pmmintrin.h>
#include <xmmintrin.h>
#include "../low/cache.hpp"
#include "../low/cache64.hpp"


class cacheline32 {
  public:
    __m128  d2[MM128_PER_CACHE_LINE];

    inline __m128 & operator[](int i) {
        return d2[i];
    }

    inline const __m128 & operator[](int i) const {
        return d2[i];
    }

    inline float & operator()(int i) {
        return ((float*)d2)[i];
    }

    inline const float & operator()(int i) const {
        return ((float*)d2)[i];
    }

    inline cacheline32 & operator=(const float & rhs) {
        d2[0] = _mm_load1_ps(&rhs);
        d2[1] = _mm_load1_ps(&rhs);
        d2[2] = _mm_load1_ps(&rhs);
        d2[3] = _mm_load1_ps(&rhs);

        return *this;
    }

    inline void set (const cacheline64 & rhs1, const cacheline64 & rhs2) {
        float v[FLOATS_PER_CACHE_LINE]  __attribute__((aligned(CACHE_LINE_SIZE)));

        for (int i=0;i<DOUBLES_PER_CACHE_LINE; ++i) v[i] = rhs1(i);
        for (int i=0;i<DOUBLES_PER_CACHE_LINE; ++i) v[DOUBLES_PER_CACHE_LINE+i] = rhs2(i);

        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_load_ps(v+4*i);
    }

    inline cacheline32 & operator+=(const cacheline32 & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_add_ps(d2[i], rhs.d2[i]);
        return *this;
    }

    inline cacheline32 & operator-=(const cacheline32 & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_sub_ps(d2[i], rhs.d2[i]);
        return *this;
    }

    inline cacheline32 & operator*=(const cacheline32 & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_mul_ps(d2[i], rhs.d2[i]);
        return *this;
    }

    inline cacheline32 & operator*=(const __m128 & rhs) {
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_mul_ps(d2[i], rhs);
        return *this;
    }

    inline cacheline32 & operator*=(float rhs) {
        __m128 v = _mm_load1_ps(&rhs);
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_mul_ps(d2[i], v);
        return *this;
    }


    inline cacheline32 operator+(const cacheline32 & rhs) const {
        cacheline32 ret = *this;
        ret += rhs;
        return ret;
    }

    inline cacheline32 operator-(const cacheline32 & rhs) const {
        cacheline32 ret = *this;
        ret -= rhs;
        return ret;
    }

    inline cacheline32 operator*(const cacheline32 & rhs) const {
        cacheline32 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline cacheline32 operator*(const __m128 & rhs) const {
        cacheline32 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline cacheline32 operator*(float rhs) const {
        cacheline32 ret = *this;
        ret *= rhs;
        return ret;
    }

    inline void set (float d) {
        d2[0] = _mm_set1_ps(d);
        d2[1] = _mm_set1_ps(d);
        d2[2] = _mm_set1_ps(d);
        d2[3] = _mm_set1_ps(d);
    }

    inline void set (float a[16]) {
    /*
        d2[0] = _mm_loadl_pd(d2[0], &a0);
        d2[0] = _mm_loadh_pd(d2[0], &a1);
        d2[1] = _mm_loadl_pd(d2[1], &a2);
        d2[1] = _mm_loadh_pd(d2[1], &a3);
        d2[2] = _mm_loadl_pd(d2[2], &a4);
        d2[2] = _mm_loadh_pd(d2[2], &a5);
        d2[3] = _mm_loadl_pd(d2[3], &a6);
        d2[3] = _mm_loadh_pd(d2[3], &a7);
        */
        for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
            d2[i] = _mm_load_ps(a+4*i);
    }

}  __attribute__((aligned(CACHE_LINE_SIZE)));


static inline void store(float * p, const cacheline32 & rhs) {
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        _mm_store_ps(p+4*i  , rhs.d2[i]);
}

static inline void store(cacheline32 * p, const cacheline32 & rhs) {
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        _mm_store_ps((float*)p+4*i, rhs.d2[i]);
}

static inline cacheline32 load(const float * p) {
    cacheline32 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_load_ps(p+4*i);
    return ret;
}

static inline cacheline32 load(const cacheline32 * p) {
    cacheline32 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_load_ps((float*)p+4*i);
    return ret;
}


inline cacheline32 sqrt(const cacheline32 & rhs) {
    cacheline32 ret;
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_sqrt_ps(rhs.d2[i]);
    return ret;
}

inline cacheline32 inv(const cacheline32 & rhs) {
    cacheline32 ret; __m128 one  = _mm_set1_ps(1.);
    for (int i=0; i<MM128_PER_CACHE_LINE; ++i)
        ret.d2[i] = _mm_div_ps(one, rhs.d2[i]);
    return ret;
}


void PackArrays     (const float * __restrict__  arrays, UI32 ArraySize, cacheline32 * __restrict__  array16);
void UnPackArrays   (const cacheline32 * __restrict__ array16,  UI32 ArraySize, float * arrays);

std::ostream & operator<<(std::ostream & os, const cacheline32 & rhs);

#endif
