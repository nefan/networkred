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

#include "kronred.hpp"
#include "kronred_common.hpp"
#include "kronred.cljmex.hpp"


#define COMPLEX
#include "klu_internal.h"

#include "kronred_factor.hpp"
#include <vector>
#include <cstring>
#include <omp.h>

void kronred(double *Mx, double *Mz, EncodedOps& ops) {
    kronred_factor(Mx, Mz, ops[0].data(), ops.size()) ;
}


#ifdef MEX

cljmex_start()

    const int n = M.rows;
    const int Nops = ops_in.cols;
    const int Nkr = index_in.cols;

    // Input
    EncodedOps ops;
    ops.reserve(Nops);
    {
        double* p = ops_in.x;
        for (int i=0; i<Nops; i++) {
            EncodedOp eop;
            eop[0] = *p++; eop[1] = *p++;
            eop[2] = *p++; eop[3] = *p++;

            ops.push_back(eop);
        }
    }
 
 
    double *Mx = new double[M.nz+1]; // temp arrays
    double *Mz = new double[M.nz+1]; // +1 for 'trash' location

    // Timing
    double time = -1;
    if (NITER > 1) {
        std::memcpy(Mx,M.x,M.nz*sizeof(double));
        std::memcpy(Mz,M.z,M.nz*sizeof(double));

        double t0factor = omp_get_wtime();

        for (int i=0; i<NITER; i++)
            kronred_factor(Mx, Mz, ops[0].data(), Nops) ;

        double t1factor = omp_get_wtime();
        time = (t1factor-t0factor)/NITER;
        // mexPrintf("kronred: factorization took %f ms.\n",time/1e-3);

    }

    // Factorize
    std::memcpy(Mx,M.x,M.nz*sizeof(double));
    std::memcpy(Mz,M.z,M.nz*sizeof(double));
    kronred(Mx, Mz, ops) ;
 
    // Output
 
    // Mkr
    std::vector<int> index;
    index.reserve(index_in.cols);
    for (int i=0; i<index_in.cols; i++)
        index.push_back((int)index_in.x[i]-1); // -1 for matlab indexing
    kronred_extract(M,Mx,Mz,index,Mkr);
 
     // Cleanup
     delete Mx;
     delete Mz;

cljmex_end()

#endif // MEX
