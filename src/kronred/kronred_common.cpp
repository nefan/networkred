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
#include <vector>

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
    assert(Mi[k] == row); // indexing failure

    return k;
}

// extract from M to Mkr
int kronred_extract(const cljmexComplex_sparse_matrix& M, 
        double *Mx, double *Mz,
        const std::vector<int>& index,
        cljmexComplex_sparse_matrix& Mkr) {

    // dimensions
    Mkr.rows = index.size();
    Mkr.cols = index.size();

    // filter
    std::vector<int> filter(M.cols,-1);
    for (int i=0; i<Mkr.cols; i++) {
        assert(index[i] < M.cols);
        filter[index[i]] = i;
    }

    // copy data from M to Mkr
    int k=0;
    int col=0;
    for (; col < Mkr.cols; col++)
    {
        int colM = index[col];

        // index
        Mkr.p[col] = k;

        for (int kM=M.p[colM]; kM<M.p[colM+1]; kM++) {
            int rowM = M.i[kM];

            int row = filter[rowM];
            if (row > -1) {
                Mkr.x[k] = Mx[kM]; // copy
                Mkr.z[k] = Mz[kM];
                Mkr.i[k] = row; 
                k++;
            }
        }
    }

    // nz
    Mkr.p[col] = k;
}

