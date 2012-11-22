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

#include "networkred_common.h"
#include "thevenin_reach.cljmex.hpp"

#define CS_LONG
#define CS_COMPLEX
#include "cs.h"
#define Long SuiteSparse_long
#include "ssparse.h"

#include <omp.h>

cljmex_start()

    // Setup data structures

    CLJ_CS_TO_CS_CL(Lcs,L);

    Long *P = (Long*)mxMalloc(n * sizeof(Long));
    for (int i=0; i<n; i++)
        P[i] = Pone.x[i]-1;
    Long *Pinv = cs_pinv(P,n);
    mxAssert(Pinv != NULL, "pinv allocation failed");

    cs_complex_t *B = (cs_complex_t*)mxMalloc (Link.nz * sizeof (cs_complex_t)) ;
    Long *Bi = (Long*)mxMalloc (Link.nz * sizeof (Long)) ;
    for (int c=0; c<M; c++) {
        for (int k=Link.p[c]; k<Link.p[c+1]; k++) {
            B[k] = cs_complex_t(Link.x[k],Link.z[k]);
            Bi[k] = Pinv[Link.i[k]];
        }
    }
    CLJ_CS_TO_CS_CL(Bcs,Link);
    Bcs.i = Bi;
    Bcs.x = B;

    // Global results, workspace, etc.

    CS_INT *reach = (CS_INT*)mxMalloc ((1 + 2*n) * M * sizeof (CS_INT)) ;
    for (int i=0; i<(1+2*n)*M; i++)
        reach[i] = 0;

    // Compute
    
    for (int k=0; k<M; k++) {
        CS_INT top = cs_reach (&Lcs, &Bcs, k, reach+(1+2*n)*k+1, NULL);
        reach[(1+2*n)*k] = top;
    }

    // Output

    for (int k=0; k<M; k++)
        for (int j=0; j<n+1; j++)
            outReach.x[(1+n)*k+j] = (double)reach[(1+2*n)*k+j];

    cs_free(Pinv);

cljmex_end()
