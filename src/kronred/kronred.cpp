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

cljmex_start()

    const int n = M.rows;
    const int Nops = ops.cols;
    const int Nkr = index.cols;

     // Analyze
 
     double *Mmx = (double*)mxMalloc(M.nz * sizeof (double) + 1); // + 1 for 'trash' location
     double *Mmz = (double*)mxMalloc(M.nz * sizeof (double) + 1);
     int *ops2 = (int *)mxMalloc(Nops * sizeof (int) * 4);
 
     double t0symbolic = omp_get_wtime();
     kronred_analyze (M.p, M.i, M.x, M.z, M.nz, Mmx, Mmz, ops.x, Nops, ops2) ;
     double t1symbolic = omp_get_wtime();
     // mexPrintf("kronred: symbolic analysis took %f ms.\n",(t1symbolic-t0symbolic)/1e-3);
 
     // Factorize
     
     double time = -1;
     // more precise timing
     if (NITER > 1) {
         double t0factor = omp_get_wtime();

         for (int i=0; i<NITER; i++)
             kronred_factor(Mmx, Mmz, ops2, Nops) ;

         double t1factor = omp_get_wtime();
         time = (t1factor-t0factor)/NITER;
         // mexPrintf("kronred: factorization took %f ms.\n",time/1e-3);
     }

     kronred_analyze (M.p, M.i, M.x, M.z, M.nz, Mmx, Mmz, ops.x, Nops, ops2) ;
     kronred_factor(Mmx, Mmz, ops2, Nops) ;
 
    // Output
 
     /* Mkr */
     kronred_extract(Mmx,Mmz,M.p,M.i,Nkr,index.x,filter.x,Mkr.p,Mkr.i,Mkr.x,Mkr.z);
 
cljmex_end()

