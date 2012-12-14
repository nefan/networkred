/*
 * This file is part of networkred.
 *
 * Copyright (C) 2012, Technical University of Denmark
 * https://github.com/nefan/networkred
 *
 * networkred is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * networkred is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with networkred.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "kronred_factor.hpp"

#ifndef USE_AVX_SIMD

#define CS_LONG
#define CS_COMPLEX
#include "cs.h"
#include <complex.h>

void kronred_factor 
(
    double Mx[], 
    double Mz[], 
    int ops[], 
    int Nops
)
{
    for (int k=0; k<Nops; k++) {
        int K = -ops[4*k+0];

        std::complex<double> nodev(Mx[ops[4*k+1]],Mz[ops[4*k+1]]);

        int i = k+1;
        for (int j=0; j<K; j++, i++) {

            std::complex<double> node = 
                std::complex<double>(Mx[ops[4*i+1]],Mz[ops[4*i+1]])*
                std::complex<double>(Mx[ops[4*i+2]],Mz[ops[4*i+2]])/
                nodev;


            Mx[ops[4*i+0]] -= real(node);
            Mz[ops[4*i+0]] -= imag(node);
        }

        k += K;
    }
}

#else

#include <smmintrin.h>

#define SCALAR_IS_LTZERO(x)     ((x) < 0.)
/* scalar absolute value macro. If x is NaN, the result is NaN: */
#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO (x)) ? -(x) : (x))

typedef float v4sf __attribute__ ((vector_size (32))); // 4 vector
typedef double v4sd __attribute__ ((vector_size (32))); // 4 vector
typedef int v4si __attribute__ ((vector_size (32))); // 4 vector
typedef union 
{
    v4sf simd;
    float v[4];
} v4f;
typedef union 
{
    v4sd simd;
    double v[4];
} v4d;
typedef union 
{
    v4si simd;
    int v[4];
} v4i;

void kronred_factor 
(
    double Mx[], 
    double Mz[], 
    int ops[], 
    int Nops
)
{

    for (int k=0; k<Nops; k++) {

        int K = -ops[4*k+0];

        int i3 = ops[4*k+1];
        double dr = Mx[i3];
        double di = Mz[i3];

        int branch = SCALAR_ABS (dr) >= SCALAR_ABS (di);
        double r, den;
        if (branch) {
            r = di / dr;
            den = dr + r * di;
        } else {
            r = dr / di;
            den = r * dr + di;
        }
        v4d v4r = { r, r, r, r };
        v4d v4den = { den, den, den, den };

        int i = k+1;
        for (int j=0; j<K; j +=4, i +=4) {
            int ip0 = i;
            int ip1 = ip0+1;
            int ip2 = ip1+1;
            int ip3 = ip2+1;

            v4d v4cr = { Mx[ops[4*ip0+2]], Mx[ops[4*ip1+2]], Mx[ops[4*ip2+2]], Mx[ops[4*ip3+2]] };
            v4d v4ci = { Mz[ops[4*ip0+2]], Mz[ops[4*ip1+2]], Mz[ops[4*ip2+2]], Mz[ops[4*ip3+2]] };
            v4d v4vr, v4vi;

            if (branch) {
                v4vr.simd = (v4cr.simd + v4ci.simd * v4r.simd) / v4den.simd;
                v4vi.simd = (v4ci.simd - v4cr.simd * v4r.simd) / v4den.simd;
            } else {
                v4vr.simd = (v4cr.simd * v4r.simd + v4ci.simd) / v4den.simd;
                v4vi.simd = (v4ci.simd * v4r.simd - v4cr.simd) / v4den.simd;
            }

            int i0p0 = ops[4*ip0+0];
            int i0p1 = ops[4*ip1+0];
            int i0p2 = ops[4*ip2+0];
            int i0p3 = ops[4*ip3+0];
            v4d v4ar = { Mx[i0p0], Mx[i0p1], Mx[i0p2], Mx[i0p3] };
            v4d v4ai = { Mz[i0p0], Mz[i0p1], Mz[i0p2], Mz[i0p3] };
            v4d v4br = { Mx[ops[4*ip0+1]], Mx[ops[4*ip1+1]], Mx[ops[4*ip2+1]], Mx[ops[4*ip3+1]] };
            v4d v4bi = { Mz[ops[4*ip0+1]], Mz[ops[4*ip1+1]], Mz[ops[4*ip2+1]], Mz[ops[4*ip3+1]] };

            v4d v4resx, v4resz;
            v4resx.simd = v4ar.simd - v4vr.simd * v4br.simd + v4vi.simd * v4bi.simd;
            v4resz.simd = v4ai.simd - v4vi.simd * v4br.simd - v4vr.simd * v4bi.simd;
            Mx[i0p0] = v4resx.v[0]; Mx[i0p1] = v4resx.v[1]; Mx[i0p2] = v4resx.v[2]; Mx[i0p3] = v4resx.v[3];
            Mz[i0p0] = v4resz.v[0]; Mz[i0p1] = v4resz.v[1]; Mz[i0p2] = v4resz.v[2]; Mz[i0p3] = v4resz.v[3];
        }
        
        k += K;
    }
}

#endif
