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

#include "kronred_common.hpp"

#include "mex.h"

int linear_index
(
 int row,
 int col,
 long *Mp,
 long *Mi
)
{
    int k=Mp[col];
    for (; k<Mp[col+1]; k++)
        if (Mi[k] == row)
            break;
    // mexPrintf("%d %d %d %d %f %f\n",row,col,k,Mi[k],Mx[k],Mz[k]);
    mxAssert(Mi[k] == row, "indexing failure");

    return k;
}

int kronred_analyze
(
    long *Mp,
    long *Mi,
    double *Mx,
    double *Mz,
    int nz,
    double *Mmx,
    double *Mmz,
    double *ops1,
    int Nops,
    int *ops2
)
{
    /* M is complex */
    for (long k=0 ; k<nz ; k++)
    {
        if (!mxIsNaN(Mx[k])) {
            Mmx[k] = Mx [k] ;        /* real part */
            Mmz[k] = Mz [k] ;        /* imaginary part */
        } else {
            Mmx[k] = 0;
            Mmz[k] = 0;
        }
    }

    for (int k=0; k<Nops; k++) {
        int K = -ops1[2*k+0];
        mxAssert(K >= 0, "ops node indexing failure");
        mxAssert(K % 4 == 0, "ops node indexing failure (fail mod 4)");
        ops2[4*k+0] = -K;

        int node = (int)ops1[2*k+1];
        mxAssert(node >= 0, "ops indexing failure");
        ops2[4*k+1] = linear_index(node-1,node-1,Mp,Mi);

        int i = k+1;
        for (int j=0; j<K; j +=1, i +=1) {
            if (ops1[2*i+1] < 0) { // trash
                ops2[4*i+0] = nz;
                ops2[4*i+1] = nz;
                ops2[4*i+2] = nz;
            } else {
                ops2[4*i+0] = linear_index((int)ops1[2*i+0]-1,(int)ops1[2*i+1]-1,Mp,Mi);
                ops2[4*i+1] = linear_index((int)ops1[2*i+0]-1,node-1,Mp,Mi);
                ops2[4*i+2] = linear_index(node-1,(int)ops1[2*i+1]-1,Mp,Mi);
            }
        }

        k += K;
    }
}

// extract from M to Mkr
int kronred_extract
(
    double *Mmx,
    double *Mmz,
    long *Mp,
    long *Mi,
    int Nkr,
    double *index,
    double *filter,
    long *Mkrp,
    long *Mkri,
    double *Mkrx,
    double *Mkrz
)
{
    int k=0;
    int col=0;
    for (; col < Nkr; col++)
    {
        int colM = (int)index[col]-1;

        // index
        Mkrp[col] = (long)k;

        for (int kM=Mp[colM]; kM<Mp[colM+1]; kM++) {
            int rowM = Mi[kM];

            long row = (long)filter[rowM]-1; // -1 for matlab indexing
            if (row > -1) {
                Mkrx[k] = Mmx[kM]; // copy
                Mkrz[k] = Mmz[kM];
                Mkri[k] = row; 
                k++;
            }
        }
    }

    // nz
    Mkrp[col] = k;
}

